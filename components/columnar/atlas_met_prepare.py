"""Project raw PHYSLITE fields into the independently loaded preparation nodes."""
from __future__ import annotations

import awkward as ak
import numpy as np

from atlas_met_columnar import OWNERS, P4, TRUTH, VERTICES, flat, pack_met, put

EVENT = "EventInfoAuxDyn."
JET = OWNERS["jet"]
TAU = OWNERS["tau"]
MUON = OWNERS["muon"]
ELECTRON = OWNERS["electron"]
PHOTON = OWNERS["photon"]
CLUSTER = "egammaClusters"
CONVERSION = "GSFConversionVertices"
ID_OWNERS = ("InDetTrackParticles", "InDetForwardTrackParticles", "GSFTrackParticles")
TRACK_OWNERS = (*ID_OWNERS, "CombinedMuonTrackParticles", "ExtrapolatedMuonTrackParticles",
                "MuonSpectrometerTrackParticles", "MSOnlyExtrapolatedMuonTrackParticles")
P4_OWNERS = [ELECTRON, PHOTON, MUON, TAU]
MUON_LINKS = {
    "id": "inDetTrackParticleLink",
    "me": "extrapolatedMuonSpectrometerTrackParticleLink",
    "cb": "combinedTrackParticleLink",
    "msonly": "msOnlyExtrapolatedMuonSpectrometerTrackParticleLink",
    "ms": "muonSpectrometerTrackParticleLink",
}
SUMMARY = {
    "pixel_hits": "numberOfPixelHits", "pixel_dead": "numberOfPixelDeadSensors",
    "sct_hits": "numberOfSCTHits", "sct_dead": "numberOfSCTDeadSensors",
    "pixel_holes": "numberOfPixelHoles", "sct_holes": "numberOfSCTHoles",
}
MUON_FIELDS = {
    "precision_layers": "numberOfPrecisionLayers", "precision_hole_layers": "numberOfPrecisionHoleLayers",
    "good_precision_layers": "numberOfGoodPrecisionLayers", "inner_small_hits": "innerSmallHits",
    "inner_large_hits": "innerLargeHits", "middle_small_hits": "middleSmallHits",
    "middle_large_hits": "middleLargeHits", "outer_small_hits": "outerSmallHits",
    "outer_large_hits": "outerLargeHits", "extended_small_hits": "extendedSmallHits",
    "extended_large_hits": "extendedLargeHits", "extended_small_holes": "extendedSmallHoles",
    "is_small_good_sectors": "isSmallGoodSectors", "eta_layer1_stgc_hits": "etaLayer1STGCHits",
    "eta_layer2_stgc_hits": "etaLayer2STGCHits", "mm_hits": "MMHits",
    "combined_out_of_bounds_hits": "combinedTrackOutBoundsPrecisionHits",
}
TRACK_FIELDS = ("qOverP", "theta", "phi", "z0", "d0", "vz", "definingParametersCovMatrixDiag", *SUMMARY.values())
ADDITIONAL = [
    *(EVENT + field for field in (
        "eventNumber", "runNumber", "mcChannelNumber", "RandomRunNumber", "eventTypeBitmask",
        "actualInteractionsPerCrossing", "averageInteractionsPerCrossing",
        "beamPosSigmaX", "beamPosSigmaY", "beamPosSigmaXY",
    )),
    *(TRUTH + field for field in ("name","source","mpx","mpy","sumet")),
    "PrimaryVerticesAuxDyn.z",
    *(JET + field for field in ("NumTrkPt500","EMFrac","Width","PSFrac","EnergyPerSampling")),
    *(TAU + field for field in (
        "charge","ptFinalCalib","etaFinalCalib","ptTauEnergyScale","etaTauEnergyScale","PanTau_DecayMode",
        "isTauFlags","TESCompatibility","tauTrackLinks","truthParticleLink","truthJetLink",
    )),
    "TruthTausAuxDyn.pdgId", "TruthTausAuxDyn.IsHadronicTau", "TauTracksAuxDyn.flagSet",
    *(MUON + field for field in (
        "muonType","author","allAuthors","EnergyLoss","energyLossType","InnerDetectorPt","MuonSpectrometerPt",
        "CaloMuonScore","momentumBalanceSignificance","scatteringCurvatureSignificance",
        "scatteringNeighbourSignificance","ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000","topoetcone20",
        *MUON_LINKS.values(), *MUON_FIELDS.values(),
    )),
    *(owner+"AuxDyn."+field for owner in TRACK_OWNERS for field in TRACK_FIELDS),
    *(CLUSTER+"AuxDyn."+field for field in (
        "calEta","calPhi","ETACALOFRAME","PHICALOFRAME","etaCalo","samplingPattern","e_sampl","eta_sampl",
    )),
    *(CONVERSION+"AuxDyn."+field for field in ("x","y","pt1","pt2","trackParticleLinks")),
    *(owner+field for owner in (ELECTRON,PHOTON) for field in (
        "author","OQ","caloClusterLinks","Eadded_Lr2","Eadded_Lr3",
        "topoetcone20","topoetcone30","topoetcone40",
        "topoetcone20ptCorrection","topoetcone30ptCorrection","topoetcone40ptCorrection",
    )),
    *(ELECTRON+field for field in (
        "trackParticleLinks","DFCommonElectronsLHLooseBL",
        "ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000",
    )),
    *(PHOTON+field for field in ("vertexLinks","DFCommonPhotonsIsEMTight","DFCommonPhotonsCleaning","ptcone20")),
]


def dense(array, dtype, default=0):
    return np.ascontiguousarray(ak.to_numpy(ak.fill_none(array, default)), dtype=dtype)


def roles_for_events(offs):
    return np.repeat(np.arange(len(offs)-1, dtype="i8"), np.diff(offs))


class Link:
    def __init__(self, array):
        if len(array)==0:
            self.key = np.empty(0,"u4")
            self.index = np.empty(0,"u4")
            return
        names = ak.fields(array)
        key = next((name for name in names if name.endswith("m_persKey")), None)
        index = next((name for name in names if name.endswith("m_persIndex")), None)
        if key is None or index is None:
            raise ValueError(f"Unrecognized original ElementLink fields {names}")
        self.key = dense(array[key], "u4")
        self.index = dense(array[index], "u4", 0xffffffff)


class Resolved:
    def __init__(self, data, link, event_rows, candidates):
        self.data = data
        self.size = len(event_rows)
        self.locations = []
        self.valid = np.zeros(self.size, dtype=bool)
        if link is None:
            return
        if len(link.key) != self.size or len(link.index) != self.size:
            raise ValueError("ElementLink and owning-row shape mismatch")
        for target in candidates:
            key = data.keys.get(target)
            anchor = data.anchors.get(target)
            if key is None or anchor is None or anchor not in data.raw:
                continue
            local, = np.nonzero((link.key==key) & (link.index!=0xffffffff))
            offsets = data.offset(anchor)
            events = event_rows[local]
            inside = link.index[local] < (offsets[events+1]-offsets[events])
            source = local[inside]
            rows = offsets[events[inside]] + link.index[source]
            self.locations.append((target, source, rows))
            self.valid[source] = True

    def column(self, field, dtype="f4", default=0, slot=None, lengths=False):
        output = np.full(self.size, default, dtype=dtype)
        present = np.zeros(self.size, dtype=bool)
        for target, local, rows in self.locations:
            name = target+"AuxDyn."+field
            if name not in self.data.raw:
                continue
            source = ak.flatten(self.data.raw[name], axis=1)[rows]
            if lengths:
                source = ak.num(source, axis=1)
            elif slot is not None:
                source = ak.pad_none(source, slot+1, axis=1)[:,slot]
            missing = np.asarray(ak.is_none(source, axis=0))
            output[local] = dense(source, dtype, default)
            present[local] = ~missing
        return output, present

    def nested_link(self, field, slot):
        keys, indices = self.link_fields(field, slot)
        obj = Link.__new__(Link)
        obj.key, obj.index = keys, indices
        return obj

    def link_fields(self, field, slot):
        keys = np.zeros(self.size, dtype="u4")
        indices = np.full(self.size, 0xffffffff, dtype="u4")
        for target, local, rows in self.locations:
            name = target+"AuxDyn."+field
            if name in self.data.raw:
                source = ak.flatten(self.data.raw[name], axis=1)[rows]
                link = Link(ak.pad_none(source, slot+1, axis=1)[:,slot])
                keys[local], indices[local] = link.key, link.index
        return keys, indices

    def sampling(self, field, sampling_pattern, slot):
        # The packed array stores one value per set sampling bit, in enum order.
        output = np.zeros(self.size, dtype="f4")
        for target, local, rows in self.locations:
            source = ak.flatten(self.data.raw[target+"AuxDyn."+field], axis=1)[rows]
            relevant = (sampling_pattern[local] & (1<<slot)) != 0
            if not relevant.any():
                continue
            chosen = local[relevant]
            masks = sampling_pattern[chosen] & ((1<<slot)-1)
            rank = np.zeros(len(chosen), dtype="i8")
            for bit in range(slot):
                rank += (masks>>bit)&1
            selected = source[relevant]
            count = dense(ak.num(selected, axis=1), "i8")
            start = np.concatenate(([0], np.cumsum(count, dtype="i8")))
            good = rank < count
            output[chosen] = -999
            content = dense(ak.flatten(selected, axis=1), "f4")
            output[chosen[good]] = content[start[:-1][good]+rank[good]]
        return output


class Data:
    anchors = {
        **{name: name+"AuxDyn.qOverP" for name in TRACK_OWNERS},
        CLUSTER: CLUSTER+"AuxDyn.calEta",
        CONVERSION: CONVERSION+"AuxDyn.x",
        "TauTracks": "TauTracksAuxDyn.flagSet",
        "TruthTaus": "TruthTausAuxDyn.pdgId",
    }
    def __init__(self, raw, output, keys, evaluate, capture_preparation):
        # Split scalar ElementLinks are grouped as child-name dictionaries by uproot.
        # Zip their same-length event/object leaves into typed link records before indexing.
        self.raw = {name: ak.zip(value) if isinstance(value,dict) else value for name,value in raw.items()}
        self.output, self.keys, self.evaluate = output, keys, evaluate
        self.capture_preparation = capture_preparation
        self.counts = {}
        self.events = {owner: roles_for_events(output[short+".o"]) for short,owner in OWNERS.items()}

    def offset(self, anchor):
        if anchor not in self.counts:
            self.counts[anchor] = flat(self.raw[anchor], 1, "f4")[0][0]
        return self.counts[anchor]

    def own(self, owner, field, dtype="f4", default=None):
        n = len(self.events[owner])
        name = owner+field
        if name not in self.raw:
            if default is None:
                raise ValueError(f"Required original leaf {name} is missing")
            return np.full(n, default, dtype=dtype)
        return dense(ak.flatten(self.raw[name], axis=1), dtype, default or 0)

    def off(self, field, owner, dtype="f4", default=0):
        return self.own(owner,field,dtype,default), int(owner+field in self.raw)

    def resolve(self, owner, field, targets, slot=None):
        arr = self.raw.get(owner+field)
        if arr is not None:
            arr = ak.flatten(arr, axis=1)
            if slot is not None:
                arr = ak.pad_none(arr, slot+1, axis=1)[:,slot]
        return Resolved(self,None if arr is None else Link(arr),self.events[owner],targets)

    def event(self, owner, field, dtype="f4", default=None):
        name = EVENT+field
        if name not in self.raw:
            if default is None:
                raise ValueError(f"Required event leaf {name} is missing")
            return np.full(len(self.events[owner]),default,dtype=dtype)
        return dense(self.raw[name],dtype)[self.events[owner]]

    def publish(self, name, p4, selected):
        for k,i in zip(P4,range(4),strict=True):
            self.output[name+".prep."+k] = np.ascontiguousarray(p4[:,i],dtype="f4")
        self.output[name+".selected"] = np.ascontiguousarray(selected,dtype="u1")

    def gate_rows(self, name, *decisions):
        if self.capture_preparation:
            self.output[name+".gate"] = np.ascontiguousarray(np.column_stack(decisions),dtype="u4").ravel()


def availability(out, name, values, valid):
    out[name+".value"] = np.ascontiguousarray(values,dtype="f4")
    out[name+".valid"] = np.ascontiguousarray(valid,dtype="u1")


def prepare_tau(data):
    track_links = data.raw[TAU+"tauTrackLinks"]
    per_tau = ak.flatten(track_links,axis=1)
    track_events = np.repeat(data.events[TAU],dense(ak.num(per_tau,axis=1),"i8"))
    tracks = Resolved(data, Link(ak.flatten(per_tau,axis=1)),track_events,("TauTracks",))
    flags = tracks.column("flagSet","i8")[0]
    charged = data.evaluate("tau_track_charged",{"flag_set":flags})[:,0].astype("i8")
    counts = dense(ak.sum(ak.unflatten(charged,dense(ak.num(per_tau,axis=1),"i8")),axis=1),"i8")
    truth = data.resolve(TAU,"truthParticleLink",("TruthTaus",))
    pdg = truth.column("pdgId","i8")[0]
    hadronic = truth.column("IsHadronicTau","i8",-1)[0]
    fields = {name:data.own(TAU,leaf,dtype) for name,leaf,dtype in (
        ("pt","pt","f8"),("eta","eta","f8"),("phi","phi","f8"),("mass","m","f8"),
        ("charge","charge","f8"),("pt_final","ptFinalCalib","f8"),("eta_final","etaFinalCalib","f8"),
        ("pt_tes","ptTauEnergyScale","f8"),("eta_tes","etaTauEnergyScale","f8"),
        ("is_tau_flags","isTauFlags","i8"),("pan_tau_mode","PanTau_DecayMode","i8"),
    )}
    n = len(counts)
    fields.update(charged_tracks=counts,has_tes_pt=np.ones(n,"i8"),
                  tes_compatibility=data.own(TAU,"TESCompatibility","i8",-1),
                  truth_pdg_id=pdg,truth_is_hadronic=hadronic,
                  truth_jet_link_available=np.full(n,int(TAU+"truthJetLink" in data.raw),"i8"))
    prepared=data.evaluate("tau_prepare",fields)
    data.publish("tau",prepared[:,:4],prepared[:,4])
    data.gate_rows("tau",prepared[:,4])
    for name, field, dtype in (
        ("mode","PanTau_DecayMode","i4"),("tes_eta","etaTauEnergyScale","f4"),
        ("final_eta","etaFinalCalib","f4"),("final_pt","ptFinalCalib","f4"),
    ):
        data.output["tau."+name]=data.own(TAU,field,dtype)

def primary_positions(data):
    types = data.output["vertices.type"]
    eligible = data.evaluate("muon_selection_vertex",{"vertex_type":types})[:,0] != 0
    counts = np.diff(data.output["vertices.o"])
    indices = ak.unflatten(np.arange(len(types),dtype="i8"),counts)
    positions = ak.firsts(indices[ak.unflatten(eligible,counts)])
    return dense(positions,"i8",-1)


def primary_z(data, positions, owner):
    flat_z = dense(ak.flatten(data.raw["PrimaryVerticesAuxDyn.z"],axis=1),"f4")
    z = np.zeros(len(positions),"f4")
    valid = positions>=0
    z[valid] = flat_z[positions[valid]]
    return z[data.events[owner]], valid[data.events[owner]]


def prepare_muon(data, positions):
    n = len(data.events[MUON])
    muon_type = data.own(MUON,"muonType","u4")
    role_names = ("id","me","cb","msonly","ms")
    roles = {role:data.resolve(MUON,MUON_LINKS[role],TRACK_OWNERS) for role in role_names}
    features = {"muon_type":muon_type}
    for role, target in roles.items():
        features[role+"_link_valid"] = target.valid
    chosen = data.evaluate("muon_selection_track_roles", features).astype("i8")
    features.pop("msonly_link_valid")
    features.pop("ms_link_valid")
    features["primary_role"] = chosen[:,0]
    features.update(pt=data.own(MUON,"pt"),eta=data.own(MUON,"eta"),phi=data.own(MUON,"phi"),
        author=data.own(MUON,"author","i8"), all_authors=data.own(MUON,"allAuthors","i8",-1))
    features["eta_bits"] = data.evaluate("muon_selection_eta",{"raw_eta":features["eta"]})[:,0]
    id_pt, id_pt_available = data.off("InnerDetectorPt",MUON)
    me_pt, me_pt_available = data.off("MuonSpectrometerPt",MUON)
    features.update(id_pt=id_pt,me_pt=me_pt,
                    calibrated_track_momenta_available=np.full(n,id_pt_available & me_pt_available,"i8"))
    cov = "definingParametersCovMatrixDiag"
    for role in ("id","me","cb"):
        target = roles[role]
        features[role+"_q_over_p"] = target.column("qOverP")[0]
        features[role+"_theta"] = target.column("theta")[0]
        features[role+"_var_q_over_p"] = target.column(cov,slot=4)[0]
        features[role+"_cov_diag_size"] = target.column(cov,"i8",-1,lengths=True)[0]
        phi = target.column("phi")[0]
        if role == "cb":
            features["cb_phi"] = phi
        for k, v in (("qop",features[role+"_q_over_p"]),("theta",features[role+"_theta"]),("phi",phi)):
            data.output[f"muon.{role}.{k}"] = v
        data.output[f"muon.{role}.valid"] = target.valid.astype("u1")
    mu_id_pt = data.evaluate("met_track_pt",{
        "qop":features["id_q_over_p"], "theta":features["id_theta"],
        "phi":data.output["muon.id.phi"],"valid":roles["id"].valid,
    })[:,0]
    availability(data.output,"muon.id_pt",mu_id_pt,roles["id"].valid)

    def choose(out, role_index, field, dtype="f4", default=0, slot=None, lengths=False):
        values = np.full(n,default,dtype)
        for i, role in enumerate(role_names,1):
            rows = chosen[:,role_index]==i
            target, _ = roles[role].column(field,dtype,default,slot,lengths)
            values[rows] = target[rows]
        features[out] = values

    for geometric in ("d0","phi","z0","vz"):
        choose("primary_"+geometric,0,geometric)
    choose("primary_var_d0",0,cov,slot=0)
    choose("primary_cov_diag_size",0,cov,"i8",-1,lengths=True)
    for name, field in SUMMARY.items():
        choose(name,1,field,"i8",-1)
    for name, field in MUON_FIELDS.items():
        features[name] = data.own(MUON,field,"i8",-1)
    features["pv_z"] = primary_z(data,positions,MUON)[0]
    for suffix in ("X","Y","XY"):
        features["beam_sigma_"+suffix.lower()] = data.event(MUON,"beamPosSigma"+suffix)
    features["calo_score"], has_calo = data.off("CaloMuonScore",MUON)
    features["calo_score_available"] = np.full(n,has_calo,"i8")
    mask = sum(1<<i for i,leaf in enumerate(("momentumBalanceSignificance",
        "scatteringCurvatureSignificance","scatteringNeighbourSignificance")) if MUON+leaf in data.raw)
    features["lowpt_parameters_mask"] = np.full(n,mask,"i8")
    features["track_iso"], track_iso = data.off("ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000",MUON)
    features["calo_iso"], calo_iso = data.off("topoetcone20",MUON)
    features["isolation_available_mask"] = np.full(n,track_iso | calo_iso<<1,"i8")
    decisions = data.evaluate("muon_selection_after_calibration",features)
    prepared = np.column_stack((features["pt"],features["eta"],features["phi"],np.full(n,105.658367)))
    data.publish("muon",prepared,decisions[:,5])
    data.gate_rows("muon",features["eta_bits"],decisions[:,0],decisions[:,1],decisions[:,2],decisions[:,4])
    data.output["muon.type"] = muon_type
    data.output["muon.loss"] = data.own(MUON,"EnergyLoss")
    data.output["muon.loss_type"] = data.own(MUON,"energyLossType","u1")

def prepare_egamma(data, positions, short):
    photon = short == "photon"
    owner = OWNERS[short]
    n = len(data.events[owner])
    mass = 0.0 if photon else 0.510998910
    value = {name:data.own(owner,name,"f8") for name in ("pt","eta","phi")}
    value["mass"] = np.full(n,mass,"f8")
    value["author"] = data.own(owner,"author","i8")
    value["oq"] = data.own(owner,"OQ","i8")
    value["id_flag"] = data.own(owner,"DFCommonPhotonsIsEMTight" if photon else "DFCommonElectronsLHLooseBL","i8",-1)
    value["clean_flag"] = data.own(owner,"DFCommonPhotonsCleaning","i8",-1) if photon else np.zeros(n,"i8")
    value["run_number"] = data.event(owner,"runNumber","i8")
    value["random_run_number"] = data.event(owner,"RandomRunNumber","i8",-1)
    value["pv_z"], value["has_pv"] = primary_z(data,positions,owner)
    for suffix in ("X","Y","XY"):
        value["beam_sigma_"+suffix.lower()] = data.event(owner,"beamPosSigma"+suffix)

    cluster = data.resolve(owner,"caloClusterLinks",(CLUSTER,),slot=0)
    value["has_cluster"] = cluster.valid
    value["cluster_eta"], cluster_eta = cluster.column("calEta")
    value["cluster_phi"] = cluster.column("calPhi")[0]
    value["eta_calo"], calo_eta = cluster.column("ETACALOFRAME")
    value["phi_calo"], calo_phi = cluster.column("PHICALOFRAME")
    value["has_calo_moments"] = calo_eta & calo_phi
    value["sampling_pattern"] = cluster.column("samplingPattern","i8")[0]
    for sample, slot in (("emb1",1),("emb2",2),("emb3",3),("eme1",5),("eme2",6),("eme3",7)):
        value[sample+"_energy"] = cluster.sampling("e_sampl",value["sampling_pattern"],slot)
        if slot not in (3,7):
            value[sample+"_eta"] = cluster.sampling("eta_sampl",value["sampling_pattern"],slot)

    iso_mask, old_mask = 0, 0
    for i, cone in enumerate(("20","30","40")):
        value["iso"+cone], avail = data.off("topoetcone"+cone,owner)
        iso_mask |= avail<<i
        value["oldleak"+cone], avail = data.off("topoetcone"+cone+"ptCorrection",owner)
        old_mask |= avail<<i
    track_iso_field = "ptcone20" if photon else "ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000"
    value["track_iso"], has_track_iso = data.off(track_iso_field,owner)
    value["iso_available_mask"] = np.full(n,iso_mask | has_track_iso<<3,"i8")
    value["oldleak_available_mask"] = np.full(n,old_mask,"i8")
    value["eadded_l2"], avail_l2 = data.off("Eadded_Lr2",owner)
    value["eadded_l3"], avail_l3 = data.off("Eadded_Lr3",owner)
    value["timing_available_mask"] = np.full(n,avail_l2 | avail_l3<<1,"i8")

    track = data.resolve(owner,"trackParticleLinks",ID_OWNERS,slot=0)
    value["has_track"] = track.valid
    for leaf in ("d0","phi","z0","vz"):
        value["track_"+leaf] = track.column(leaf)[0]
    covariance = "definingParametersCovMatrixDiag"
    value["track_var_d0"] = track.column(covariance,slot=0)[0]
    value["track_cov_diag_size"] = track.column(covariance,"i8",-1,lengths=True)[0]

    conversion = data.resolve(owner,"vertexLinks",(CONVERSION,),slot=0)
    value["conv_x"] = conversion.column("x")[0]
    value["conv_y"] = conversion.column("y")[0]
    value["conv_pt1"], pt1_valid = conversion.column("pt1")
    value["conv_pt2"], pt2_valid = conversion.column("pt2")
    value["has_conv_pts"] = pt1_valid & pt2_valid
    for slot in (0,1):
        role = "conv_trk"+str(slot+1)
        linked = conversion.nested_link("trackParticleLinks",slot)
        conv_track = Resolved(data,linked,data.events[owner],TRACK_OWNERS)
        value["has_"+role] = conv_track.valid
        for suffix, leaf in (("q_over_p","qOverP"),("theta","theta"),("phi","phi")):
            value[role+"_"+suffix] = conv_track.column(leaf)[0]
        for suffix, leaf in (("pixel","numberOfPixelHits"),("sct","numberOfSCTHits")):
            value[role+"_"+suffix] = conv_track.column(leaf,"i8",-1)[0]
    for private in ("leak_et","rz_x","rz_y","leak_key"):
        value[private] = np.zeros(n,"f8")

    prepared = data.evaluate(short+"_prepare",value)
    data.publish(short,prepared[:,:4],prepared[:,12])
    gates = (4,5,6,7,9,10,11) if photon else (4,5,7,8,9,10,11)
    data.gate_rows(short,*(prepared[:,index] for index in gates))
    availability(data.output,short+".cluster_eta",value["cluster_eta"],cluster_eta)
    availability(data.output,short+".eta_calo_moment",value["eta_calo"],calo_eta)
    eta_aux, have_aux = cluster.column("etaCalo")
    availability(data.output,short+".eta_calo_aux",eta_aux,have_aux)
    # xAOD EgammaDefs AuthorFwdElectron is exactly 0x8.
    data.output[short+".forward_author"] = (value["author"]==8).astype("u1")
    if photon:
        data.output["photon.vertex_count"] = dense(ak.num(ak.flatten(data.raw[owner+"vertexLinks"],axis=1),axis=1),"u4")


def prepare_chunk(output, raw, keys, evaluate, capture_preparation=False):
    data = Data(raw,output,keys,evaluate,capture_preparation)
    for name, field, dtype in (
        ("run","runNumber","u4"),("type","eventTypeBitmask","u4"),
        ("actual_mu","actualInteractionsPerCrossing","f4"),
        ("average_mu","averageInteractionsPerCrossing","f4"),("number","eventNumber","u8"),
        ("dsid","mcChannelNumber","u4"),
    ):
        output["event."+name] = dense(raw[EVENT+field],dtype)
    pack_met(output,raw,TRUTH,"truth")
    put(output,"jet.ntrk",raw[JET+"NumTrkPt500"],"i4","jet.o","jet.ntrk.o")
    output["jet.em"] = data.own(JET,"EMFrac")
    output["jet.width"] = data.own(JET,"Width")
    ps, valid = data.off("PSFrac",JET)
    availability(output,"jet.ps",ps,np.full(len(ps),valid,"u1"))
    if not valid:
        put(output,"jet.sample",raw[JET+"EnergyPerSampling"],"f4","jet.o","jet.sample.o")
    positions = primary_positions(data)
    prepare_tau(data)
    prepare_muon(data,positions)
    prepare_egamma(data,positions,"electron")
    prepare_egamma(data,positions,"photon")
    return output
