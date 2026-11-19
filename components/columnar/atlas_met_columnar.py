"""Typed NumPy/Awkward transport. Individual event rows are built in native code."""
from __future__ import annotations

from bisect import bisect_right
import ctypes as C
import json
from pathlib import Path

import awkward as ak
import numpy as np

ASSOC = "METAssoc_AnalysisMETAux."
OWNERS = {
    "jet": "AnalysisJetsAuxDyn.",
    "electron": "AnalysisElectronsAuxDyn.",
    "photon": "AnalysisPhotonsAuxDyn.",
    "tau": "AnalysisTauJetsAuxDyn.",
    "muon": "AnalysisMuonsAuxDyn.",
}
VERTICES = "PrimaryVerticesAuxDyn.vertexType"
CORE = "MET_Core_AnalysisMETAuxDyn."
TRUTH = "MET_TruthAuxDyn."
XYZEST = ("px", "py", "pz", "e", "sumpt")
P4 = ("pt", "eta", "phi", "m")

RAW_BRANCHES = [
    *(stem + k for owner, stem in OWNERS.items()
      for k in (P4 if owner in ("jet", "tau") else P4[:3])),
    *(OWNERS["jet"] + "JetConstitScaleMomentum_" + k for k in P4),
    OWNERS["jet"] + "SumPtTrkPt500", OWNERS["jet"] + "JVFCorr", VERTICES,
    *(CORE + k for k in ("name", "source", "mpx", "mpy", "sumet")),
]

DTYPES = {k: np.ctypeslib.ndpointer(dtype=k, flags="C_CONTIGUOUS") for k in
          ("u1", "i1", "u4", "u8", "i8", "f4", "f8")}


class Column(C.Structure):
    _fields_ = (("name", C.c_char_p), ("dtype", C.c_char_p), ("size", C.c_size_t),
                ("data", C.c_void_p))


class Native:
    def __init__(self, library: Path):
        self.lib = C.CDLL(str(library.resolve(strict=True)))
        size = C.c_size_t
        self.lib.met_columnar_error.restype = C.c_char_p
        self.lib.met_columnar_open.argtypes = (C.c_char_p,) * 3
        self.lib.met_columnar_open.restype = C.c_void_p
        self.lib.met_columnar_close.argtypes = (C.c_void_p,)
        self.lib.met_decode_jet_links.argtypes = (
            DTYPES["u1"], size, DTYPES["i8"], size, DTYPES["i8"],
            DTYPES["u4"], DTYPES["u4"], size,
        )
        self.lib.met_columnar_raw_smoke.argtypes = (
            C.c_void_p, size, C.POINTER(Column), size, DTYPES["f4"], size,
            DTYPES["u8"], size, DTYPES["f4"], DTYPES["i1"], size,
        )

        self.lib.met_columnar_configure.argtypes = (C.c_void_p, C.c_char_p, C.c_char_p, C.c_char_p, C.c_char_p)
        self.lib.met_columnar_reset.argtypes = (C.c_void_p,)
        self.lib.met_columnar_reset.restype = None
        self.lib.met_prep_evaluate.argtypes = (C.c_void_p, C.c_char_p, size, C.POINTER(Column), size, DTYPES["f8"], size)
        self.lib.met_columnar_full.argtypes = (C.c_void_p, size, C.c_uint64, C.POINTER(C.c_char_p), size,
            C.POINTER(Column), size, DTYPES["f4"], DTYPES["u8"], DTYPES["f4"], C.c_char_p)

    def configure(self, context, soft, significance, prep_provider, preparation):
        self.preparation = {row["name"]:row for row in json.loads(preparation.read_text())["corrections"]}
        self.check(self.lib.met_columnar_configure(context, *(
            str(path.resolve(strict=True)).encode() for path in (soft, significance, prep_provider, preparation)
        )))

    def prepare(self, context, name, columns):
        definition = self.preparation[name]
        arrays = {entry["name"]: np.ascontiguousarray(columns[entry["name"]],
                  dtype="f8" if entry["type"]=="real" else "i8") for entry in definition["inputs"]}
        rows = len(next(iter(arrays.values())))
        width = definition["output"].get("shape", [1])
        if len(width)!=1:
            raise ValueError("Preparation ABI requires a scalar or one-dimensional output row")
        out = np.empty((rows,width[0]),dtype="f8")
        batch = describe(arrays)
        self.check(self.lib.met_prep_evaluate(context,name.encode(),rows,batch,len(batch),out,width[0]))
        return out

    def full(self, context, arrays, first, events, labels, trace):
        batch = describe(arrays)
        states = (C.c_char_p * len(labels))(*(state.encode() for state in labels))
        terms = np.empty((events,len(labels),8,5),dtype="f4")
        sources = np.empty((events,len(labels),8),dtype="u8")
        significance = np.empty((events,len(labels),4),dtype="f4")
        self.check(self.lib.met_columnar_full(context,events,first,states,len(states),batch,len(batch),
                                          terms,sources,significance,None if trace is None else str(trace).encode()))
        return terms,sources,significance

    def check(self, returned: int):
        if returned:
            raise RuntimeError(self.lib.met_columnar_error().decode())

    def open(self, provider: Path, maker: Path, nnjvt: Path):
        native = self.lib.met_columnar_open(
            *(str(path.resolve(strict=True)).encode() for path in (provider, maker, nnjvt))
        )
        if not native:
            raise RuntimeError(self.lib.met_columnar_error().decode())
        return native

    def raw_maker(self, context, arrays, events):
        descriptors = describe(arrays)
        n_assoc = int(arrays["assoc.o"][-1])
        n_jets = int(arrays["jet.o"][-1])
        terms = np.empty((events, 7, 3), dtype="f4")
        flags = np.empty(n_assoc, dtype="u8")
        scores = np.empty(n_jets, dtype="f4")
        jvt = np.empty(n_jets, dtype="i1")
        self.check(self.lib.met_columnar_raw_smoke(
            context, events, descriptors, len(descriptors), terms, terms.size,
            flags, n_assoc, scores, jvt, n_jets,
        ))
        return terms, flags, scores, jvt


def describe(arrays):
    return (Column * len(arrays))(*(
        Column(name.encode(), value.dtype.str[1:].encode(), value.size, value.ctypes.data)
        for name, value in arrays.items()
    ))


def flat(array, depth, dtype):
    layout = ak.to_layout(ak.to_packed(array))
    offsets = []
    for _ in range(depth):
        if not isinstance(layout, ak.contents.ListOffsetArray):
            raise TypeError(f"Unexpected Athena jagged layout {type(layout)}")
        offsets.append(np.ascontiguousarray(layout.offsets, dtype="i8"))
        layout = layout.content
    if not isinstance(layout, ak.contents.NumpyArray):
        raise TypeError(f"Unexpected Athena leaf layout {type(layout)}")
    return offsets, np.ascontiguousarray(layout.data, dtype=dtype)


def put(output, name, array, dtype, *axes):
    offsets, data = flat(array, len(axes), dtype)
    for axis, value in zip(axes, offsets, strict=True):
        if axis in output:
            if not np.array_equal(output[axis], value):
                raise ValueError(f"Unaligned event/owner offsets at {name} and {axis}")
        else:
            output[axis] = value
    output[name] = data


class MemberwiseJets:
    def __init__(self, native, branch):
        self.native = native
        self.branch = branch
        self.offsets = branch.entry_offsets
        self.cached_index = None
        self.cached = None

    def read(self, start, stop):
        start_at = bisect_right(self.offsets, start) - 1
        keys, indices, offsets = [], [], [np.array([0], dtype="i8")]
        length = 0
        for basket_index in range(start_at, len(self.offsets) - 1):
            first, last = self.offsets[basket_index:basket_index+2]
            if first >= stop:
                break
            a, b = max(start, first) - first, min(stop, last) - first
            if basket_index != self.cached_index:
                self.cached_index = basket_index
                self.cached = self.branch.basket(basket_index)
            basket = self.cached
            source = np.ascontiguousarray(basket.data).view("u1")
            byte_offsets = np.ascontiguousarray(basket.byte_offsets[a:b+1], dtype="i8")
            count = (int(byte_offsets[-1] - byte_offsets[0]) - 16*(b-a)) // 8
            if count < 0:
                raise ValueError("Unexpected memberwise basket")
            decoded = np.empty(b-a+1, dtype="i8")
            key = np.empty(count, dtype="u4")
            index = np.empty(count, dtype="u4")
            self.native.check(self.native.lib.met_decode_jet_links(
                source, source.size, byte_offsets, b-a, decoded, key, index, count,
            ))
            if int(decoded[-1]) != count:
                raise ValueError("Unexpected decoded link count")
            keys.append(key)
            indices.append(index)
            offsets.append(decoded[1:] + length)
            length += count
        combined = np.concatenate(offsets)
        if len(combined) != stop - start + 1:
            raise ValueError("Memberwise basket entry count mismatch")
        return combined, np.concatenate(keys), np.concatenate(indices)


def pack_raw(raw, assoc, memberwise):
    output = {"assoc.o": memberwise[0], "assoc.jet_key": memberwise[1], "assoc.jet_index": memberwise[2]}
    for owner, stem in OWNERS.items():
        for name in P4 if owner in ("jet", "tau") else P4[:3]:
            put(output, owner+"."+name, raw[stem+name], "f4", owner+".o")
    for name in P4:
        put(output, "jet.const."+name, raw[OWNERS["jet"]+"JetConstitScaleMomentum_"+name], "f4", "jet.o")
    put(output,"jet.sum_pt",raw[OWNERS["jet"]+"SumPtTrkPt500"],"f4","jet.o","jet.sum_pt.o")
    put(output,"jet.jvf",raw[OWNERS["jet"]+"JVFCorr"],"f4","jet.o")
    put(output,"vertices.type",raw[VERTICES],"u1","vertices.o")
    pack_met(output, raw, CORE, "core")
    put(output,"assoc.misc",assoc["isMisc"],"i1","assoc.o")
    put(output,"assoc.object_key",assoc["objectLinks"]["m_persKey"],"u4","assoc.o","assoc.slot_o")
    put(output,"assoc.object_index",assoc["objectLinks"]["m_persIndex"],"u4","assoc.o","assoc.slot_o")
    put(output,"assoc.edge_slot",assoc["overlapIndices"],"u8","assoc.o","assoc.slot_o","assoc.edge_o")
    put(output,"assoc.edge_type",assoc["overlapTypes"],"u1","assoc.o","assoc.slot_o","assoc.edge_o")
    for vec in ("cal", "trk"):
        for coord in XYZEST:
            put(output,"assoc."+vec+coord,assoc[vec+coord],"f4","assoc.o","assoc."+vec+"_o")
        put(output,"assoc."+vec+"key",assoc[vec+"key"],"u8","assoc.o","assoc."+vec+"_o")
    for coord in XYZEST:
        put(output,"assoc.jettrk"+coord,assoc["jettrk"+coord],"f4","assoc.o")
    return output


def pack_met(output, raw, source_prefix, stem):
    for field in ("mpx", "mpy", "sumet"):
        put(output, stem+"."+field, raw[source_prefix+field], "f4", stem+".o")
    put(output, stem+".source", raw[source_prefix+"source"], "u8", stem+".o")
    put(output, stem+".name", raw[source_prefix+"name"], "u1", stem+".o", stem+".name_o")


ASSOC_FIELDS = (
    "isMisc", "objectLinks", "overlapIndices", "overlapTypes",
    *(vec + coord for vec in ("cal", "trk", "jettrk") for coord in XYZEST), "calkey", "trkkey",
)


def read_associations(tree, start, stop):
    parent = tree[ASSOC]
    return {field: parent[ASSOC+field].array(entry_start=start, entry_stop=stop, library="ak")
            for field in ASSOC_FIELDS}
