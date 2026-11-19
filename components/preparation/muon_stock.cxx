// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include <correction_plugin.h>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODMuon/MuonContainer.h>
#include <xAODRootAccess/Init.h>
#include <xAODRootAccess/TEvent.h>
#include <xAODTracking/TrackParticle.h>
#include <xAODTracking/VertexContainer.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
using Int = std::int64_t;
using Values = std::vector<correction::Variable::Type>;
using Ref = correction::Correction::Ref;
using Track = xAOD::TrackParticle;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class T, class Owner>
std::pair<T, bool> get(const Owner* owner, const char* name) {
    const SG::ConstAccessor<T> access{name};
    if (owner && access.isAvailable(*owner)) return {access(*owner), true};
    return {T{}, false};
}

template<class T, class Owner>
Int integral(const Owner* owner, const char* name) {
    const auto [value, available] = get<T>(owner, name);
    return available ? value : -1;
}

struct Diag {
    double d0 = 0, qp = 0;
    Int size = -1;
};

Diag diag(const Track* track) {
    const SG::ConstAccessor<std::vector<float>> accessor{"definingParametersCovMatrixDiag"};
    if (!track || !accessor.isAvailable(*track)) return {};
    const auto& vec = accessor(*track);
    return {vec.empty() ? 0. : vec[0], vec.size() < 5 ? 0. : vec[4], static_cast<Int>(vec.size())};
}

const Track* resolve(const ElementLink<xAOD::TrackParticleContainer>& link) {
    return link.isValid() ? *link : nullptr;
}

bool exact(float actual, float expected) {
    return (std::isnan(actual) && std::isnan(expected)) ||
           std::bit_cast<std::uint32_t>(actual) == std::bit_cast<std::uint32_t>(expected);
}

Values raw_input(const xAOD::Muon& muon, const xAOD::EventInfo& ev,
                 double pv_z, Int eta_bits, const Ref& roles) {
    const std::array<const Track*, 6> tracks{
        nullptr, resolve(muon.inDetTrackParticleLink()),
        resolve(muon.extrapolatedMuonSpectrometerTrackParticleLink()),
        resolve(muon.combinedTrackParticleLink()),
        resolve(muon.msOnlyExtrapolatedMuonSpectrometerTrackParticleLink()),
        resolve(muon.muonSpectrometerTrackParticleLink()),
    };
    const auto type = static_cast<Int>(muon.muonType());
    const auto role = roles->evaluate_array({type, Int{bool(tracks[1])}, Int{bool(tracks[2])},
        Int{bool(tracks[3])}, Int{bool(tracks[4])}, Int{bool(tracks[5])}});
    const auto* primary = tracks.at(static_cast<std::size_t>(role.at(0)));
    const auto* hits = tracks.at(static_cast<std::size_t>(role.at(1)));
    const Track* id = tracks[1];
    const Track* me = tracks[2];
    const Track* cb = tracks[3];
    const auto primary_cov = diag(primary), id_cov = diag(id), me_cov = diag(me), cb_cov = diag(cb);
    const auto [id_pt, id_pt_ok] = get<float>(&muon, "InnerDetectorPt");
    const auto [me_pt, me_pt_ok] = get<float>(&muon, "MuonSpectrometerPt");
    const auto [calo_score, calo_score_ok] = get<float>(&muon, "CaloMuonScore");
    const auto [iso_track, iso_track_ok] = get<float>(&muon, "ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000");
    const auto [iso_calo, iso_calo_ok] = get<float>(&muon, "topoetcone20");
    auto hit = [&](const char* name) { return integral<std::uint8_t>(hits, name); };
    auto muon_hit = [&](const char* name) {
        const auto value = integral<std::uint8_t>(&muon, name);
        return value >= 0 ? value : integral<std::uint8_t>(primary, name);
    };
    const Int lowpt_available = get<float>(&muon, "momentumBalanceSignificance").second |
        (get<float>(&muon, "scatteringCurvatureSignificance").second << 1) |
        (get<float>(&muon, "scatteringNeighbourSignificance").second << 2);
    return {eta_bits, muon.pt(), muon.eta(), muon.phi(), double(id_pt), double(me_pt), Int{id_pt_ok && me_pt_ok},
        type, static_cast<Int>(muon.author()), integral<std::uint16_t>(&muon, "allAuthors"),
        static_cast<Int>(role[0]), Int{bool(id)}, Int{bool(me)}, Int{bool(cb)},
        primary ? double(primary->d0()) : 0., primary ? double(primary->phi0()) : 0.,
        primary ? double(primary->z0()) : 0., primary ? double(primary->vz()) : 0.,
        primary_cov.d0, primary_cov.size, pv_z,
        double(ev.beamPosSigmaX()), double(ev.beamPosSigmaY()), double(ev.beamPosSigmaXY()),
        hit("numberOfPixelHits"), hit("numberOfPixelDeadSensors"), hit("numberOfSCTHits"),
        hit("numberOfSCTDeadSensors"), hit("numberOfPixelHoles"), hit("numberOfSCTHoles"),
        id ? double(id->qOverP()) : 0., me ? double(me->qOverP()) : 0., cb ? double(cb->qOverP()) : 0.,
        id ? double(id->theta()) : 0., me ? double(me->theta()) : 0., cb ? double(cb->theta()) : 0.,
        cb ? double(cb->phi0()) : 0., id_cov.qp, me_cov.qp, cb_cov.qp, id_cov.size, me_cov.size, cb_cov.size,
        muon_hit("numberOfPrecisionLayers"), muon_hit("numberOfPrecisionHoleLayers"), muon_hit("numberOfGoodPrecisionLayers"),
        integral<std::uint8_t>(&muon, "innerSmallHits"), integral<std::uint8_t>(&muon, "innerLargeHits"),
        integral<std::uint8_t>(&muon, "middleSmallHits"), integral<std::uint8_t>(&muon, "middleLargeHits"),
        integral<std::uint8_t>(&muon, "outerSmallHits"), integral<std::uint8_t>(&muon, "outerLargeHits"),
        integral<std::uint8_t>(&muon, "extendedSmallHits"), integral<std::uint8_t>(&muon, "extendedLargeHits"),
        integral<std::uint8_t>(&muon, "extendedSmallHoles"), integral<std::uint8_t>(&muon, "isSmallGoodSectors"),
        integral<std::uint8_t>(&muon, "etaLayer1STGCHits"), integral<std::uint8_t>(&muon, "etaLayer2STGCHits"),
        integral<std::uint8_t>(&muon, "MMHits"), integral<std::uint8_t>(&muon, "combinedTrackOutBoundsPrecisionHits"),
        double(calo_score), Int{calo_score_ok}, lowpt_available,
        double(iso_track), double(iso_calo), Int{iso_track_ok | (iso_calo_ok << 1)}};
}
} // namespace

int main(int argc, char** argv) try {
    require(argc == 5, "Usage: met_muon_prepare_stock_check PLUGIN MUON_JSON ORIGINAL OBSERVER_ROOT");
    require(xAOD::Init("met_muon_prepare_stock_check").isSuccess(), "xAOD initialize");
    std::unique_ptr<TFile> original(TFile::Open(argv[3], "READ"));
    std::unique_ptr<TFile> observer(TFile::Open(argv[4], "READ"));
    require(original && !original->IsZombie() && observer && !observer->IsZombie(), "Open input and observer");
    xAOD::TEvent reader(xAOD::TEvent::kClassAccess);
    require(reader.readFrom(original.get()).isSuccess(), "Read original PHYSLITE");
    TTreeReader reference("analysis", observer.get());
    // The observer writes an unavailable mu_m auxiliary branch as zero.
    // Muon_v1::m() returns the release constant in its methods and p4.
    TTreeReaderArray<float> pt(reference, "mu_pt"), eta(reference, "mu_eta"), phi(reference, "mu_phi");
    TTreeReaderArray<std::uint32_t> eta_bits(reference, "mu_selectEta"), pt_bits(reference, "mu_selectPt"),
        track_bits(reference, "mu_trackSelection_medium");
    TTreeReaderArray<char> quality(reference, "mu_good_muon_medium"), iso(reference, "mu_isolated_muon_medium");

    correction::PluginRegistry registry;
    registry.load_library(argv[1]);
    const auto set = correction::CorrectionSet::from_file(argv[2], registry);
    const auto eta_select = set->at("muon_selection_eta");
    const auto vertex_select = set->at("muon_selection_vertex");
    const auto roles = set->at("muon_selection_track_roles");
    const auto select = set->at("muon_selection_after_calibration");
    const auto isolation = set->at("muon_selection_isolation");
    require(isolation->evaluate({20000., -3500., 5900., Int{3}}) == 1. &&
            isolation->evaluate({20000., 3000., 6000., Int{3}}) == 1. &&
            isolation->evaluate({20000., 2900., 6001., Int{3}}) == 0. &&
            isolation->evaluate({20000., 0., 0., Int{2}}) == 0., "Loose_VarRad edge and availability fixtures");

    constexpr Int pass = 0xffffffff;
    std::array<std::size_t, 6> counts{};
    Long64_t entry = 0;
    std::size_t rows = 0;
    while (reference.Next()) {
        require(reader.getEntry(entry) >= 0, "Read original event");
        const xAOD::MuonContainer* muons = nullptr;
        const xAOD::VertexContainer* vertices = nullptr;
        const xAOD::EventInfo* event = nullptr;
        require(reader.retrieve(muons, "AnalysisMuons").isSuccess() && muons, "Retrieve original muons");
        require(reader.retrieve(vertices, "PrimaryVertices").isSuccess() && vertices, "Retrieve original vertices");
        require(reader.retrieve(event, "EventInfo").isSuccess() && event, "Retrieve original event info");
        double pv_z = 0.;
        for (const auto* pv : *vertices) {
            if (vertex_select->evaluate({static_cast<Int>(pv->vertexType())})) { pv_z = pv->z(); break; }
        }
        require(pt.GetSize() == muons->size() && eta_bits.GetSize() == muons->size(),
                "Muon row count mismatch");
        for (std::size_t i = 0; i < muons->size(); ++i) {
            const auto& muon = *muons->at(i);
            const auto raw_eta_bits = static_cast<Int>(eta_select->evaluate({muon.eta()}));
            const auto candidate = select->evaluate_array(raw_input(muon, *event, pv_z, raw_eta_bits, roles));
            const bool selected = eta_bits[i] == pass && pt_bits[i] == pass && track_bits[i] == pass && quality[i] && iso[i];
            if (raw_eta_bits != eta_bits[i] || candidate[0] != pt_bits[i] || candidate[1] != track_bits[i] ||
                candidate[2] != (quality[i] != 0) || candidate[4] != (iso[i] != 0) || candidate[5] != selected ||
                !exact(muon.pt(), pt[i]) || !exact(muon.eta(), eta[i]) ||
                !exact(muon.phi(), phi[i])) {
                std::cerr << "entry=" << entry << " muon_index=" << i << " eta=" << raw_eta_bits << '/' << eta_bits[i]
                          << " pt_mask=" << static_cast<Int>(candidate[0]) << '/' << pt_bits[i]
                          << " track=" << static_cast<Int>(candidate[1]) << '/' << track_bits[i]
                          << " pt=" << muon.pt() << '/' << pt[i] << " eta_value=" << muon.eta() << '/' << eta[i]
                          << " phi=" << muon.phi() << '/' << phi[i]
                          << " medium=" << candidate[2] << '/' << int(quality[i])
                          << " iso=" << candidate[4] << '/' << int(iso[i]) << '\n';
                return 1;
            }
            counts[0] += raw_eta_bits == pass;
            counts[1] += candidate[0] == pass;
            counts[2] += candidate[1] == pass;
            counts[3] += quality[i] != 0;
            counts[4] += iso[i] != 0;
            counts[5] += selected;
            ++rows;
        }
        ++entry;
    }
    require(rows != 0 && entry == reader.getEntries(), "Expected one observer entry for every original event");
    std::cout << "entries=" << entry << " rows=" << rows << " eta=" << counts[0]
              << " pt=" << counts[1] << " track=" << counts[2] << " medium=" << counts[3]
              << " isolation=" << counts[4] << " selected=" << counts[5] << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
