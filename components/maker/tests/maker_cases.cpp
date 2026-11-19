// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "met_maker/maker.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {
using namespace atlas::met::maker;

constexpr Mask tst_source = 0x14020;
constexpr Mask central_tst_source = tst_source | 0x100;
constexpr std::uint8_t pflow = 1U << 3;

void expect(bool passed, const char* message)
{
    if (!passed) {
        throw std::runtime_error(message);
    }
}

ObjectId id(std::uint32_t owner, std::uint32_t index)
{
    return {{owner}, index};
}

ObjectLink link(ObjectId key)
{
    return {LinkState::valid, key};
}

Particle particle(ParticleKind kind, ObjectId key, double pt, double eta = 0.)
{
    Particle p{};
    p.id = key;
    p.kind = kind;
    p.p4 = {pt, eta, 0., 0.};
    return p;
}

const Term& term(const Result& r, TermType type)
{
    for (const auto& t : r.terms) {
        if (t.type == type) {
            return t;
        }
    }
    throw std::runtime_error("missing term");
}

bool saw(const Result& r, Issue issue)
{
    for (const auto& d : r.diagnostics) {
        if (d.issue == issue) {
            return true;
        }
    }
    return false;
}

Config synthetic_config()
{
    return Config({{20000.F, 30000.F, 60000.F},
                   {-2.5F, 0.F, 2.5F},
                   {0.2F, 0.4F, 0.6F, 0.8F}});
}

// NNJvtBinning.cxx:66-86 and JvtSelectionToolBase.cxx:29-34.
void cut_lookup()
{
    const auto config = synthetic_config();
    expect(config.jvt_cut(15000.F, 0.F) == 0.2F, "low pt folds below first edge");
    expect(config.jvt_cut(30000.F, 0.F) == 0.2F, "interior edges are right-closed");
    expect(config.jvt_cut(std::nextafter(30000.F, 31000.F), 0.F) == 0.6F,
           "one float above an edge enters the next interval");
    expect(config.jvt_cut(30000.F, std::nextafter(0.F, 1.F)) == 0.4F,
           "cut lookup uses signed eta");
    expect(config.jvt_cut(30000.F, -2.5F) == -1.F, "minimum signed eta edge is outside");
    expect(config.jvt_cut(60000.F, 2.5F) == 0.8F, "last edges are included");
    expect(config.jvt_cut(60001.F, 0.F) == -1.F, "upper overflow is -1");
}

// METMaker.cxx:311-444 and MissingETComposition.cxx:375-401.
void overlap_and_identity()
{
    const Particle e = particle(ParticleKind::electron, id(10, 0), 18000);
    const Particle p = particle(ParticleKind::photon, id(11, 0), 25000);
    const Particle missing = particle(ParticleKind::electron, id(10, 1), 5000);
    const std::array source{e, p, missing};

    std::array electrons{
        particle(ParticleKind::electron, id(20, 9), 19000),
        particle(ParticleKind::electron, id(20, 7), 6000)
    };
    electrons[0].original_link = link(e.id);
    electrons[1].original_link = link(missing.id);

    const std::array<OverlapEdge, 1> shared{{{1, pflow}}};
    std::array<std::span<const OverlapEdge>, 3> overlap_rows{};
    overlap_rows[2] = shared;
    const std::array links{ObjectLink{}, link(e.id), link(p.id)};
    const std::array<OverlapEdge, 1> misc_shared{{{0, pflow}}};
    std::array<std::span<const OverlapEdge>, 2> misc_overlaps{};
    misc_overlaps[1] = misc_shared;
    const std::array misc_links{link(e.id), link(p.id)};
    const std::array<KeyedVector, 3> misc_tracks{{
        {1, {5, 0, 0, 5, 5}}, {2, {7, 0, 0, 7, 7}}, {3, {9, 0, 0, 9, 9}}
    }};
    std::array<Association, 2> associations{};
    associations[0].objects = links;
    associations[0].overlap_rows = overlap_rows;
    associations[1].is_misc = true;
    associations[1].objects = misc_links;
    associations[1].overlap_rows = misc_overlaps;
    associations[1].trk = misc_tracks;

    const std::array core{CoreTerm{"PVSoftTrkCore", tst_source, {2.F, 3.F, 11.F}}};
    EventInput input{};
    input.associations = associations;
    input.core_terms = core;
    input.input_particles = source;
    input.electrons.selected = electrons;
    input.electrons.view.has_shallow_store = true;
    input.photons.selected = std::span(&p, 1);
    const auto config = synthetic_config();
    const Result result = evaluate(config, input);
    expect(result.status == Status::success, "overlap fixture status");
    expect(result.flags == std::vector<Mask>{2, 1}, "per-slot flags across associations");
    const auto& eterm = term(result, TermType::electron);
    expect(eterm.momentum.mpx == -19000.F && eterm.momentum.sumet == 19000.F,
           "term uses the current electron p4");
    expect(eterm.contributors.size() == 2 && eterm.contributors.at(0).current == electrons[0].id
           && eterm.contributors.at(0).weight == 1.F
           && eterm.contributors.at(1).current == electrons[1].id
           && eterm.contributors.at(1).weight == 0.F,
           "ordered calibrated-owner links include the missing electron with weight zero");
    expect(term(result, TermType::photon).contributors.empty(),
           "PFlow edge vetoes the later photon across matching associations");
    expect(saw(result, Issue::missing_electron_association), "missing electron diagnostic");
    const auto& soft = term(result, TermType::soft_track);
    expect(soft.momentum.mpx == -5.F && soft.momentum.mpy == 3.F && soft.momentum.sumet == 18.F,
           "only the unflagged miscellaneous key supplements the core");
    expect(term(result, TermType::muon_eloss).source == 0x40008,
           "empty muon view still creates MuonEloss");
    electrons[0].original_link.reset();
    const Result shallow = evaluate(config, input);
    expect(shallow.status == Status::success && saw(shallow, Issue::shallow_missing_original),
           "shallow hard collection missing its original link skips but succeeds");
}

// METMaker.cxx:714-758, 951-1068 and Jet_v1.cxx:90-97.
void weighted_jet_and_jvt()
{
    const Particle e = particle(ParticleKind::electron, id(10, 0), 11000);
    Jet jet{};
    jet.particle = particle(ParticleKind::jet, id(14, 0), 32000);
    jet.constituent_p4 = PtEtaPhiM{40000., 0., 0., 0.};
    jet.jvt_eta = 0.F;
    jet.nnjvt = 0.8F;
    const std::array source{e, jet.particle};
    const std::array objects{link(e.id), link(jet.particle.id)};
    const std::array<std::span<const OverlapEdge>, 2> overlap_rows{};
    const std::array<KeyedVector, 2> calorimeter{{
        {1, {5000, 0, 0, 5000, 5000}}, {2, {20000, 0, 0, 20000, 20000}}
    }};
    const std::array<KeyedVector, 1> tracks{{{1, {4, 0, 0, 4, 4}}}};
    Association association{};
    association.jet = link(jet.particle.id);
    association.objects = objects;
    association.overlap_rows = overlap_rows;
    association.cal = calorimeter;
    association.trk = tracks;
    association.jet_tracks = {10, 0, 0, 10, 10};
    const std::array core{CoreTerm{"PVSoftTrkCore", central_tst_source, {2, 0, 8}}};
    EventInput input{};
    input.associations = std::span(&association, 1);
    input.core_terms = core;
    input.input_particles = source;
    input.electrons.selected = std::span(&e, 1);
    input.jets = std::span(&jet, 1);
    const auto config = synthetic_config();

    Result result = evaluate(config, input);
    expect(result.status == Status::success, "clean weighted-jet event");
    const auto& hard = term(result, TermType::jet);
    expect(hard.contributors.size() == 1 && hard.contributors.front().weight == 0.75F,
           "alternate vector scales the selected calorimeter overlap");
    expect(hard.momentum.mpx == -24000.F && hard.momentum.sumet == 24000.F,
           "hard term uses calibrated weighted momentum");

    jet.nnjvt = 0.6F;  // Exactly the synthetic cut at (32 GeV, 0).
    result = evaluate(config, input);
    const auto& soft = term(result, TermType::soft_track);
    expect(result.jet_decisions.at(0).jvt_pass == false, "NNJVT cut is strict");
    expect(soft.momentum.mpx == -4.F && soft.momentum.sumet == 14.F,
           "JVT reject still adds residual tracks to central TST");
    expect(soft.contributors.size() == 1 && !soft.contributors[0].weight,
           "TST soft-jet link has no stored weight");

    jet.particle.p4.eta = 2.5;
    jet.jvt_eta = 2.5F;
    result = evaluate(config, input);
    expect(term(result, TermType::soft_track).contributors.empty(),
           "eta 2.5 is outside the central soft term");

    jet.particle.p4.eta = 4.6;
    jet.jvt_eta = 4.6F;
    jet.nnjvt.reset();
    result = evaluate(config, input);
    expect(result.status == Status::success
           && result.jet_decisions.at(0).jvt_pass == true
           && result.jet_decisions.at(0).eta_rejected_before_jvt && !result.jet_decisions.at(0).rejected,
           "out-of-range JVT accept overwrites the maximum eta rejection");
}

// METMaker.cxx:839-949 and 1025-1052.
void muon_fsr_and_removal()
{
    Particle muon = particle(ParticleKind::muon, id(12, 2), 32000);
    muon.muon_type = 0;
    muon.id_track_pt = 30000.F;
    muon.energy_loss = 10000.F;
    muon.energy_loss_type = 0;
    Jet jet{};
    jet.particle = particle(ParticleKind::jet, id(14, 0), 60001.);
    jet.constituent_p4 = PtEtaPhiM{40000., 0., 0., 0.};
    const std::array ntrk{10, 2};
    const std::array sumpt{40000.F, 32000.F};
    jet.num_trk_pt500 = ntrk;
    jet.sum_pt_trk_pt500 = sumpt;
    jet.ps_frac = 3000.F;
    jet.em_frac = 0.95F;
    jet.width = 0.02F;
    const std::array objects{link(muon.id)};
    const std::array<std::span<const OverlapEdge>, 1> overlap_rows{};
    Association association{};
    association.jet = link(jet.particle.id);
    association.objects = objects;
    association.overlap_rows = overlap_rows;
    association.jet_tracks = {10, 0, 0, 10, 10};
    const std::array core{CoreTerm{"PVSoftTrkCore", tst_source, {}}};
    const std::array<std::uint8_t, 2> vertices{0, 1};
    EventInput input{};
    input.associations = std::span(&association, 1);
    input.core_terms = core;
    input.input_particles = std::span(&muon, 1);
    input.muons.selected = std::span(&muon, 1);
    input.jets = std::span(&jet, 1);
    input.vertex_types = vertices;
    const auto config = synthetic_config();

    Result result = evaluate(config, input);
    expect(result.status == Status::success, "muon-jet FSR succeeds");
    expect(result.jet_decisions.at(0).disposition == JetDisposition::muon_eloss,
           "muon FSR routes to MuonEloss");
    const auto& ell = term(result, TermType::muon_eloss);
    expect(ell.momentum.mpx == -30000.F && ell.momentum.sumet == 30000.F,
           "EM scale jet subtracts the muon energy loss");
    expect(ell.contributors.empty() && term(result, TermType::jet).contributors.empty()
           && term(result, TermType::soft_track).contributors.empty(),
           "FSR records only the vector");

    jet.ps_frac.reset();
    const std::array sampling{2000.F, 0.F, 0.F, 0.F, 1000.F};
    jet.energy_per_sampling = sampling;
    muon.energy_loss = 40000.F;
    result = evaluate(config, input);
    expect(result.status == Status::success
           && result.jet_decisions.at(0).disposition == JetDisposition::muon_eloss
           && term(result, TermType::muon_eloss).momentum.mpx == 0.F,
           "sampling fallback and energy-loss saturation");

    jet.energy_per_sampling.reset();
    result = evaluate(config, input);
    expect(result.status == Status::failure && saw(result, Issue::missing_presampler_input),
           "active FSR evaluation needs PSFrac or EnergyPerSampling");

    jet.particle.p4.pt = 35000.;
    jet.jvt_eta = 0.F;
    jet.nnjvt = 0.8F;
    result = evaluate(config, input);
    expect(result.status == Status::success && result.jet_decisions.at(0).rejected
           && term(result, TermType::soft_track).momentum.mpx == -10.F,
           "muon removal rejects before checking the inactive presampler input");
}

// METMaker.cxx:407-424 and MissingETAssociationHelper.cxx:28-49.
void calo_muon_and_signed_flag()
{
    const Particle electron = particle(ParticleKind::electron, id(10, 0), 11000);
    Particle muon = particle(ParticleKind::muon, id(12, 0), 15000);
    muon.muon_type = 3;
    const std::array source{electron, muon};
    const std::array objects{link(electron.id), link(muon.id)};
    const std::array<OverlapEdge, 1> non_pflow_edge{{{0, 2}}};
    std::array<std::span<const OverlapEdge>, 2> overlap_rows{};
    overlap_rows[1] = non_pflow_edge;
    Association assoc{};
    assoc.objects = objects;
    assoc.overlap_rows = overlap_rows;
    const std::array core{CoreTerm{"PVSoftTrkCore", tst_source, {}}};
    EventInput input{};
    input.associations = std::span(&assoc, 1);
    input.core_terms = core;
    input.input_particles = source;
    input.electrons.selected = std::span(&electron, 1);
    input.muons.selected = std::span(&muon, 1);
    const auto config = synthetic_config();

    Result result = evaluate(config, input);
    expect(result.status == Status::success && result.flags.at(0) == 3,
           "non-PFlow edge allows flagging the muon before its additional veto");
    expect(term(result, TermType::muon).contributors.empty(),
           "CaloTagged muon veto does not clear its previously set bit");

    std::array<ObjectLink, 33> wide_objects{};
    std::array<std::span<const OverlapEdge>, 33> wide_overlaps{};
    wide_objects[31] = link(electron.id);
    assoc.objects = wide_objects;
    assoc.overlap_rows = wide_overlaps;
    input.muons.selected = {};
    result = evaluate(config, input);
    expect(result.status == Status::success && result.flags.at(0) == 0xffffffff80000000ULL,
           "release signed flag setter sign-extends slot 31");

    wide_objects[31] = {};
    wide_objects[32] = link(electron.id);
    result = evaluate(config, input);
    expect(result.status == Status::stock_undefined
           && saw(result, Issue::stock_signed_shift_out_of_range),
           "wide stock signed shift is surfaced rather than emulated");
}

// METMaker.cxx:301-309, 469-474, 597-610, and map finder 72-78 and 186-192.
void map_and_missing_input()
{
    const std::array core{CoreTerm{"PVSoftTrkCore", tst_source, {2, 0, 8}}};
    EventInput input{};
    input.core_terms = core;
    const auto config = synthetic_config();
    Result result = evaluate(config, input);
    expect(result.status == Status::success && result.terms.size() == 7
           && !result.constituent_links_available
           && term(result, TermType::soft_track).momentum.mpx == 0.F,
           "empty map creates terms and skips seeding from the core");

    Jet jet{};
    jet.particle = particle(ParticleKind::jet, id(14, 1), 35000.);
    const std::array<KeyedVector, 1> last_misc_tracks{{{1, {5, 0, 0, 5, 5}}}};
    std::array<Association, 3> associations{};
    associations[0].is_misc = true;
    associations[1].jet = link(jet.particle.id);
    associations[2].is_misc = true;
    associations[2].trk = last_misc_tracks;
    input.associations = associations;
    input.jets = std::span(&jet, 1);
    result = evaluate(config, input);
    expect(result.status == Status::success
           && result.jet_decisions.at(0).disposition == JetDisposition::missing_association,
           "the jet lookup stops at the first miscellaneous association");
    expect(term(result, TermType::soft_track).momentum.mpx == -3.F,
           "the soft complement uses the last miscellaneous association");

    input.core_terms = {};
    result = evaluate(config, input);
    expect(result.status == Status::failure && saw(result, Issue::missing_soft_core)
           && result.terms.size() == 6,
           "missing core fails after creating RefJet and before PVSoftTrk");
}

}  // namespace

int main()
{
    cut_lookup();
    overlap_and_identity();
    weighted_jet_and_jvt();
    muon_fsr_and_removal();
    calo_muon_and_signed_flag();
    map_and_missing_input();
}
