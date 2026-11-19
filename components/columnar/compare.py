#!/usr/bin/env python3
"""Compare independent AnalysisBase A with original-input columnar C diagnostics."""
from __future__ import annotations

from collections import Counter
from itertools import groupby, zip_longest
import json
from pathlib import Path
import struct
import sys

import awkward as ak
import numpy as np
import typer

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from compare_framework import compare_numeric
from histogram_reference import met_gev, root_file
from summarize_observer import COLLECTIONS, FINAL_QUANTITIES, MAKER_LEPTONS, TERMS

ROLES = dict(zip(COLLECTIONS, ("jet", "electron", "photon", "muon", "tau"), strict=True))
EXACT = {"associationFlag", "jetJvt", "rebuildMET", "rebuildJetMET", "makerAlg",
         "builderAlg", "significanceAlg", "setRandomSeed", "applySystematicVariation"}
NORMALIZED = EXACT | {"gaussian", "applyCorrection", "varianceMET"}


def summarize(test):
    failed = np.flatnonzero(~np.asarray(test))
    return {"count": len(failed), "first_entries": failed[:8].tolist()}


def prepared_fields(stock, config, entries):
    maker = config["seq1_MetMakerAlg_AnaMET"]["Properties"]
    selection = {prefix: maker[species+"Selection"] for species,prefix in MAKER_LEPTONS.items()}
    known_owner_keys = {}
    prepared = {}
    for prefix, role in ROLES.items():
        fetch = stock[f"{prefix}_originalObjectLink"].__getitem__
        key = fetch(f"{prefix}_originalObjectLink.m_persKey").array(entry_stop=entries)
        index = fetch(f"{prefix}_originalObjectLink.m_persIndex").array(entry_stop=entries)
        if not ak.all(index == ak.local_index(index)):
            raise ValueError(f"A {prefix} changed original collection order")
        keys = np.unique(ak.to_numpy(ak.flatten(key)))
        if len(keys) != 1:
            raise ValueError(f"No unique raw owner for {prefix}")
        known_owner_keys[role] = int(keys[0])
        selected = ak.ones_like(index,dtype=bool)
        gates = {}
        for token in selection.get(prefix, "").split("&&"):
            if token:
                name,_,accessor = token.partition(",")
                gate = stock[f"{prefix}_{name}"].array(entry_stop=entries)
                gates[name] = gate
                selected = selected & ((gate != 0) if accessor == "as_char" else (gate == 0xffffffff))
        p4 = {name: stock[f"{prefix}_{name}"].array(entry_stop=entries) for name in ("pt","eta","phi","m")}
        prepared[role] = ak.to_list(ak.zip(dict(key=key,index=index,selected=selected,**p4,**gates)))
    return prepared, known_owner_keys


def frame_trace(handle, entries):
    def parse(line):
        entry,state,action,*fields = line.rstrip("\n").split("\t")
        return int(entry),int(state),action,fields

    for entry,rows in groupby(map(parse,handle),key=lambda row: row[0]):
        if entry >= entries:
            break
        yield entry, [(state,action,fields) for _,state,action,fields in rows]


def numeric_match(a, b, quantity):
    reference, given = (np.asarray([value],dtype="f4") for value in (a,b))
    return compare_numeric(reference,given,quantity)["total_failures"] == 0


def getter_rows(path, entries):
    reverse_roles = {prefix: role for prefix,role in ROLES.items()}
    with path.open() as file:
        rows = (line.rstrip("\n").split("\t") for line in file)
        for entry, group in groupby(rows,key=lambda row: int(row[0])):
            if entry >= entries:
                break
            yield entry, [[reverse_roles[role],*fields] for _,role,*fields in group]


def trace_match(stock, config, entries, states, stock_file, given_file, stock_getters):
    prepared, raw_key = prepared_fields(stock,config,entries)
    # Stock's key is its shallow-container hash. C's current object key is the
    # original owner with bit 31 set. The event-local index and original owner
    # are independently proved, and each shallow stock owner must stay unique.
    stock_key = {}
    stock_links = {
        state: (
            ak.to_list(stock[f"fullMet_constituentLinks_{state}"].array(entry_stop=entries)),
            ak.to_list(stock[f"fullMet_constituentWeights_{state}"].array(entry_stop=entries)),
        ) for state in states
    }
    compared = Counter()
    failures = Counter()
    examples = {}
    worst_p4_delta = {}
    accessor_mass = {"electron":0.510998910,"muon":105.658367}
    getter_stream = getter_rows(stock_getters,entries) if stock_getters is not None else None

    def check(kind, entry, equal, reference=None, actual=None):
        compared[kind] += 1
        if not equal:
            failures[kind] += 1
            examples.setdefault(kind,{"entry":entry,"stock":reference,"columnar":actual})

    def normalize(rows):
        by_action = {}
        for state,action,values in rows:
            if action in NORMALIZED:
                by_action.setdefault(action,[]).append((state,values))
        return by_action

    with stock_file.open() as s, given_file.open() as c:
        for entry, reference, observed in zip_longest(range(entries),frame_trace(s,entries),frame_trace(c,entries)):
            if reference is None or observed is None or reference[0] != entry or observed[0] != entry:
                raise ValueError(f"Missing or reordered trace event {entry}")
            a, b = (normalize(row[1]) for row in (reference,observed))
            for action in EXACT:
                lhs, rhs = (rows.get(action,[]) for rows in (a,b))
                check(action+"_length",entry,len(lhs) == len(rhs),len(lhs),len(rhs))
                for left,right in zip_longest(lhs,rhs):
                    check(action,entry,left == right,left,right)
            for action in ("gaussian","applyCorrection","varianceMET"):
                lhs, rhs = (rows.get(action,[]) for rows in (a,b))
                check(action+"_length",entry,len(lhs) == len(rhs),len(lhs),len(rhs))
                for (left_state,left),(right_state,right) in zip(lhs,rhs):
                    if action == "gaussian":
                        same = list(map(float.fromhex,left)) == list(map(float.fromhex,right))
                    elif action == "varianceMET":
                        same = left[0] == right[0]
                    else:
                        same = left[0] == right[0] and all(
                            numeric_match(float.fromhex(x),float.fromhex(y),name)
                            for x,y,name in zip(left[1:],right[1:],["mpx","mpy","sumet"]*2,strict=True))
                    check(action,entry,left_state == right_state and same,left,right)

            observed_seeds = [values for _,action,values in observed[1] if action=="seed"]
            for state, seed in enumerate(observed_seeds):
                setter = int(dict(a.get("setRandomSeed",[]))[state][0])
                root_values = [int(values[0]) for slot,action,values in reference[1]
                               if slot == state and action=="rootSeed"]
                expected = [str(setter if setter < 2**63 else setter-2**64),
                            str(root_values[-1] % 2**32),
                            "called" if state in dict(a.get("gaussian",[])) else "skipped"]
                check("seed",entry,seed == expected,expected,seed)
            check("seed_length",entry,len(observed_seeds) == len(states),len(states),len(observed_seeds))

            stock_getter_lookup = None
            if getter_stream is not None:
                stock_entry, getter_expected = next(getter_stream, (None, []))
                check("getter_entry",entry,stock_entry==entry,entry,stock_entry)
                getter_actual = [values for _,action,values in observed[1] if action=="preparedGetter"]
                sort_key = lambda row: (row[0],int(row[1]))
                getter_expected.sort(key=sort_key)
                getter_actual.sort(key=sort_key)
                stock_getter_lookup = {(row[0],int(row[1])): row for row in getter_expected}
                check("getter_count",entry,len(getter_expected)==len(getter_actual),
                      len(getter_expected),len(getter_actual))
                for left,right in zip_longest(getter_expected,getter_actual):
                    if left is None or right is None:
                        check("getter_identity",entry,False,left,right)
                        continue
                    check("getter_identity",entry,left[:5]==right[:5],left[:5],right[:5])
                    for field,want,have in zip(("pt","eta","phi","m","e"),left[5:],right[5:],strict=True):
                        exact = struct.pack("=d",float.fromhex(want)) == struct.pack("=d",float.fromhex(have))
                        label = f"getter_{left[0]}_{field}"
                        check(label,entry,exact,want,have)
                        if left[4] == "1":
                            check(label+"_selected",entry,exact,want,have)

            objects = []
            for role, rows in prepared.items():
                objects.extend((role,row) for row in rows[entry])
            # C preserves order within each container. Cross-container order is immaterial.
            seen = {}
            for _,action,values in observed[1]:
                if action == "prepared":
                    role,key,index,selected,*outputs = values
                    p4,gates = list(map(float,outputs[:4])),list(map(int,outputs[4:]))
                    seen.setdefault(role,[]).append((int(key),int(index),bool(int(selected)),p4,gates))
            for role,data in seen.items():
                check("prepared_count",entry,len(data) == len(prepared[role][entry]),len(prepared[role][entry]),len(data))
            for role,row in objects:
                rows = seen.get(role,[])
                if row["index"] >= len(rows):
                    check("prepared_identity",entry,False,[role,row["index"]],None)
                    continue
                key,index,selected,p4,gates = rows[row["index"]]
                want = (row["key"],row["index"],row["selected"])
                have = (key,index,selected)
                check("prepared_identity",entry,want == have,want,have)
                ignored = {"key","index","selected","pt","eta","phi","m"}
                expected_gates = {name:value for name,value in row.items() if name not in ignored}
                check("prepared_gate_count",entry,len(gates)==len(expected_gates),len(expected_gates),len(gates))
                for (name,expected),given in zip(expected_gates.items(),gates):
                    check(f"prepared_{role}_{name}",entry,given==expected,expected,given)
                for name,given in zip(("pt","eta","phi","m"),p4,strict=True):
                    stock_p4 = row[name]
                    label = f"prepared_{role}_{name}"
                    # A writes the absent mutable aux mass as zero for these owners.
                    # Their physical masses are taken from the pinned release m() accessors.
                    if name == "m" and role in accessor_mass and stock_p4 == 0:
                        stock_p4 = accessor_mass[role]
                        label += "_release_accessor"
                    quantity = "mpx" if name in ("pt","m") else name
                    worst_p4_delta[label] = max(abs(stock_p4-given),worst_p4_delta.get(label,0))
                    check(label,entry,numeric_match(stock_p4,given,quantity),stock_p4,given)

            by_state = {}
            egamma_by_state = {}
            for state,action,values in observed[1]:
                if action == "contributor":
                    by_state.setdefault(state,[]).append(values)
                if action == "egammaSignificanceGetter":
                    egamma_by_state.setdefault(state,[]).append(values)
            for state,label in enumerate(states):
                links,weights = stock_links[label]
                want = []
                want_egamma = []
                for i,(term,(prefix,_)) in enumerate(TERMS.items()):
                    for j,link in enumerate(links[entry][i]):
                        role = ROLES[prefix]
                        key = link["m_persKey"]
                        check("stock_shallow_owner",entry,stock_key.setdefault(role,key)==key,stock_key[role],key)
                        weight = weights[entry][i][j] if weights[entry][i] else None
                        want.append((term,raw_key[role] | 2**31,raw_key[role],link["m_persIndex"],weight))
                        if stock_getter_lookup is not None and role in ("electron","photon"):
                            index = link["m_persIndex"]
                            values = stock_getter_lookup[(role,index)]
                            want_egamma.append([term,role,str(index),*values[-2:]])
                have = [(term,int(current),int(raw),int(index),float(weight) if weight else None)
                        for term,current,raw,index,weight in by_state.get(state,[])]
                check("ordered_constituent_count",entry,len(want)==len(have),len(want),len(have))
                for left,right in zip_longest(want,have):
                    check("ordered_constituent_link_weight",entry,left == right,left,right)
                if stock_getter_lookup is not None:
                    got_egamma = egamma_by_state.get(state,[])
                    check("significance_egamma_getter_count",entry,len(want_egamma)==len(got_egamma),
                          len(want_egamma),len(got_egamma))
                    for left,right in zip_longest(want_egamma,got_egamma):
                        if left is None or right is None:
                            check("significance_egamma_getter_identity",entry,False,left,right)
                            continue
                        check("significance_egamma_getter_identity",entry,left[:3]==right[:3],left[:3],right[:3])
                        for name,want_value,got_value in zip(("mass","energy"),left[-2:],right[-2:],strict=True):
                            exact = struct.pack("=d",float.fromhex(want_value)) == struct.pack("=d",float.fromhex(got_value))
                            check(f"significance_egamma_production_{name}",entry,exact,want_value,got_value)
    return {"compared":compared,"failures":failures,"first_failures":examples,
            "stock_getters":str(stock_getters) if stock_getters is not None else None,
            "worst_p4_delta":worst_p4_delta,
            "owner_map":{role: {"original":key,"stock_shallow":stock_key.get(role),"columnar_current":key | 2**31}
                         for role,key in raw_key.items()}}


def main(stock_dir: Path, stock_histogram: Path, candidate_observables: Path, candidate_histogram: Path,
         output: Path, stock_trace: Path | None = None, candidate_trace: Path | None = None,
         stock_getters: Path | None = None,
         expected_events: int | None = typer.Option(None,min=1)):
    with root_file(stock_dir/"hist-reference.root") as file:
        states = file["systematicsList"].axis().labels()
    with (root_file(stock_dir/"reference.root") as reference, np.load(candidate_observables) as candidate,
          np.load(stock_histogram) as stock_hist, np.load(candidate_histogram) as candidate_hist):
        stock = reference["analysis"]
        entries = len(candidate["entry"])
        if not 0 < entries <= stock.num_entries or (expected_events is not None and entries != expected_events):
            raise ValueError(f"Invalid candidate event count {entries}")
        if candidate["states"].tolist() != states or candidate["term_names"].tolist() != list(TERMS):
            raise ValueError("Incorrect state or term order")
        report = {"entries":entries,"states":states,"numerical":{}, "exact":{}, "event_bins":{}}
        for a,b in (("eventNumber","event_number"),("mcChannelNumber","mc_channel_number"),("runNumber","run_number")):
            expected = stock[a].array(entry_stop=entries,library="np")
            report["exact"][a] = summarize(expected == candidate[b])
        report["exact"]["entry"] = summarize(np.arange(entries) == candidate["entry"])
        for j,state in enumerate(states):
            for i,quantity in enumerate(candidate["term_fields"]):
                branch = f"fullMet_{quantity}_{state}"
                expected = stock[branch].array(entry_stop=entries)
                report["numerical"][branch] = compare_numeric(expected,candidate["terms"][:,j,:,i],quantity)
            branch = f"fullMet_source_{state}"
            expected = ak.to_numpy(stock[branch].array(entry_stop=entries))
            report["exact"][branch] = summarize(np.all(expected == candidate["sources"][:,j,:],axis=1))
            for i,quantity in enumerate(FINAL_QUANTITIES):
                branch = f"finalMet_{quantity}_{state}"
                expected = stock[branch].array(entry_stop=entries)
                report["numerical"][branch] = compare_numeric(expected,candidate["significance"][:,j,i],quantity)
            edges = stock_hist["bin_edges_gev"]
            expected = met_gev(stock[f"finalMet_met_{state}"].array(entry_stop=entries,library="np"))
            given = met_gev(candidate["terms"][:,j,-1,3])
            report["event_bins"][state] = summarize(np.searchsorted(edges,expected,side="right")
                                                   == np.searchsorted(edges,given,side="right"))

        report["histogram"] = {key: np.array_equal(stock_hist[key],candidate_hist[key]) for key in stock_hist.files}
        if set(stock_hist.files) != set(candidate_hist.files):
            raise ValueError("Unexpected histogram fields")
        if bool(stock_trace) != bool(candidate_trace) or (stock_getters is not None and stock_trace is None):
            raise ValueError("Provide independent stock and candidate status traces, optionally with the stock getter trace")
        if stock_trace is not None:
            config = json.loads((stock_dir/"full_config.json").read_text())
            report["trace"] = trace_match(stock,config,entries,states,stock_trace,candidate_trace,stock_getters)
    report["passed"] = (not any(row["total_failures"] for row in report["numerical"].values())
                        and not any(row["count"] for row in (*report["exact"].values(),*report["event_bins"].values()))
                        and all(report["histogram"].values())
                        and not report.get("trace",{}).get("failures"))
    with output.open("x") as file:
        json.dump(report,file,indent=2,allow_nan=False)
        file.write("\n")
    if not report["passed"]:
        raise typer.Exit(1)


if __name__ == "__main__":
    typer.run(main)
