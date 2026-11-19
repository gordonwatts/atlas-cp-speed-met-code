// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include <correction_plugin.h>

#include <IsolationCorrections/IsolationCorrectionTool.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <egammaUtils/egPhotonWrtPoint.h>
#include <xAODCaloEvent/CaloCluster.h>
#include <xAODEgamma/ElectronContainer.h>
#include <xAODEgamma/PhotonContainer.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODRootAccess/Init.h>
#include <xAODRootAccess/TEvent.h>
#include <xAODTracking/TrackParticle.h>
#include <xAODTracking/VertexContainer.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {
using Int = std::int64_t;
using Values = std::vector<correction::Variable::Type>;
using Ref = correction::Correction::Ref;
using Track = xAOD::TrackParticle;

constexpr auto cones = std::array{xAOD::Iso::topoetcone20, xAOD::Iso::topoetcone30, xAOD::Iso::topoetcone40};
constexpr auto em = std::array{CaloSampling::EMB1, CaloSampling::EMB2, CaloSampling::EMB3,
                               CaloSampling::EME1, CaloSampling::EME2, CaloSampling::EME3};
constexpr auto em_eta = std::array{CaloSampling::EMB1, CaloSampling::EMB2, CaloSampling::EME1, CaloSampling::EME2};

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<typename T, class Owner>
std::pair<T, bool> get(const Owner* owner, const char* name) {
    const SG::ConstAccessor<T> access{name};
    if (owner && access.isAvailable(*owner)) return {access(*owner), true};
    return {T{}, false};
}

template<typename T, class Owner>
Int integral(const Owner* owner, const char* name) {
    const auto [value, available] = get<T>(owner, name);
    return available ? static_cast<Int>(value) : -1;
}

bool exact(float a, float b) {
    return (std::isnan(a) && std::isnan(b)) || std::bit_cast<std::uint32_t>(a) == std::bit_cast<std::uint32_t>(b);
}

std::size_t compare_tables(const correction::CorrectionSet& set, const char* path) {
    std::ifstream stock(path);
    std::string header;
    std::getline(stock, header);
    require(header == "correction key x y result_float_bits", "Expected an exported ROOT probe file");
    std::string name;
    Int key;
    double x, y;
    std::uint32_t expected;
    std::size_t rows = 0;
    while (stock >> name >> key >> x >> y >> expected) {
        const double now = name == "egamma_pointing_depth" ? set.at(name)->evaluate({x, y}) : set.at(name)->evaluate({x, key});
        if (!exact(static_cast<float>(now), std::bit_cast<float>(expected))) {
            std::cerr << "probe=" << rows << " correction=" << name << " key=" << key << " x=" << x << " y=" << y
                      << " got_bits=" << std::bit_cast<std::uint32_t>(static_cast<float>(now))
                      << " expected_bits=" << expected << '\n';
            throw std::runtime_error("ROOT leakage/profile numerical probe");
        }
        ++rows;
    }
    require(stock.eof() && rows != 0, "Expected complete ROOT probe rows");
    return rows;
}

template<class Object>
Values raw_input(const Object& particle, const xAOD::EventInfo& event, const xAOD::Vertex* pv) {
    constexpr bool photon = std::is_same_v<Object, xAOD::Photon>;
    const auto* cluster = particle.caloCluster();
    const Track* track = nullptr;
    const xAOD::Vertex* conversion = nullptr;
    if constexpr (photon) conversion = particle.vertex();
    else track = particle.trackParticle();
    const Track* ct1 = conversion && conversion->nTrackParticles() > 0 ? conversion->trackParticle(0) : nullptr;
    const Track* ct2 = conversion && conversion->nTrackParticles() > 1 ? conversion->trackParticle(1) : nullptr;
    double eta_calo = 0., phi_calo = 0.;
    const bool moments = cluster && cluster->retrieveMoment(xAOD::CaloCluster::ETACALOFRAME, eta_calo)
        && cluster->retrieveMoment(xAOD::CaloCluster::PHICALOFRAME, phi_calo);
    const auto [eadded2, have_timing2] = get<float>(&particle, "Eadded_Lr2");
    const auto [eadded3, have_timing3] = get<float>(&particle, "Eadded_Lr3");
    const auto [iso_track, have_iso_track] = get<float>(&particle, photon ? "ptcone20" :
        "ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000");
    const auto [pt1, have_pt1] = get<float>(conversion, "pt1");
    const auto [pt2, have_pt2] = get<float>(conversion, "pt2");
    const auto [cov, have_cov] = get<std::vector<float>>(track, "definingParametersCovMatrixDiag");
    Values row{
        particle.pt(), particle.eta(), particle.phi(), particle.m(), cluster ? cluster->eta() : 0.,
        cluster ? cluster->phi() : 0., eta_calo, phi_calo, pv ? double(pv->z()) : 0.,
    };
    for (auto layer : em) row.emplace_back(cluster ? double(cluster->eSample(layer)) : 0.);
    for (auto layer : em_eta) row.emplace_back(cluster ? double(cluster->etaSample(layer)) : 0.);
    Int iso_available = Int{have_iso_track} << 3, leak_available = 0;
    for (std::size_t i = 0; i < cones.size(); ++i) {
        float sum = 0.;
        if (particle.isolationValue(sum, cones[i])) iso_available |= Int{1} << i;
        row.emplace_back(double(sum));
    }
    for (std::size_t i = 0; i < cones.size(); ++i) {
        float leak = 0.;
        if (particle.isolationCaloCorrection(leak, cones[i], xAOD::Iso::ptCorrection)) leak_available |= Int{1} << i;
        row.emplace_back(double(leak));
    }
    for (auto value : std::array{
        double(eadded2), double(eadded3), double(iso_track), track ? double(track->d0()) : 0.,
        track ? double(track->phi0()) : 0., track ? double(track->z0()) : 0., track ? double(track->vz()) : 0.,
        cov.empty() ? 0. : double(cov[0]), double(event.beamPosSigmaX()), double(event.beamPosSigmaY()),
        double(event.beamPosSigmaXY()), conversion ? double(conversion->x()) : 0.,
        conversion ? double(conversion->y()) : 0., double(pt1), double(pt2)
    }) row.emplace_back(value);
    for (auto ct : {ct1, ct2}) {
        // The fixed PHYSLITE schema has no first-measurement parameter arrays.
        require(!ct || ct->numberOfParameters() == 0 || (have_pt1 && have_pt2),
                "Unexpected conversion measurement fallback schema");
        row.emplace_back(ct ? double(ct->qOverP()) : 0.);
        row.emplace_back(ct ? double(ct->theta()) : 0.);
        row.emplace_back(ct ? double(ct->phi0()) : 0.);
    }
    row.emplace_back(0.);
    row.emplace_back(0.);
    row.emplace_back(0.);
    for (auto value : std::array{
        Int{bool(cluster)}, Int{moments}, Int{cluster ? cluster->samplingPattern() : 0u}, Int{bool(pv)},
        Int{event.runNumber()}, integral<unsigned>(&event, "RandomRunNumber"), Int{particle.OQ()},
        Int{particle.author()}, integral<char>(&particle, photon ? "DFCommonPhotonsIsEMTight" : "DFCommonElectronsLHLooseBL"),
        photon ? integral<char>(&particle, "DFCommonPhotonsCleaning") : Int{0}, iso_available, leak_available,
        Int{have_timing2 | (have_timing3 << 1)}, Int{bool(track)}, have_cov ? static_cast<Int>(cov.size()) : Int{-1},
        Int{bool(ct1)}, Int{bool(ct2)}, integral<std::uint8_t>(ct1, "numberOfPixelHits"),
        integral<std::uint8_t>(ct1, "numberOfSCTHits"), integral<std::uint8_t>(ct2, "numberOfPixelHits"),
        integral<std::uint8_t>(ct2, "numberOfSCTHits"), Int{have_pt1 && have_pt2}, Int{0}
    }) row.emplace_back(value);
    return row;
}

constexpr Int pass = 0xffffffff;

struct Species {
    TTreeReaderArray<float> pt, eta, phi;
    TTreeReaderArray<std::uint32_t> gate_eta, oq, pt_gate, extra_pt, track_or_clean;
    TTreeReaderArray<char> id, isolation;

    Species(TTreeReader& tree, bool photon)
        : pt(tree, photon ? "pho_pt" : "ele_pt"), eta(tree, photon ? "pho_eta" : "ele_eta"),
          phi(tree, photon ? "pho_phi" : "ele_phi"),
          gate_eta(tree, photon ? "pho_selectEta" : "ele_selectEta"), oq(tree, photon ? "pho_goodOQ" : "ele_goodOQ"),
          pt_gate(tree, photon ? "pho_selectPt" : "ele_selectPt"), extra_pt(tree, photon ? "pho_selectPtEta" : "ele_selectPtEta"),
          track_or_clean(tree, photon ? "pho_isClean" : "ele_trackSelection_loose"),
          id(tree, photon ? "pho_selectEM_tight" : "ele_selectLikelihood_loose"),
          isolation(tree, photon ? "pho_isolated_tight" : "ele_isolated_loose") {}
};

template<class Object>
void compare_row(Long64_t entry, std::size_t i, const Object& raw, const xAOD::EventInfo& ev,
                 const xAOD::Vertex* pv, const Ref& provider, Species& stock,
                 CP::IsolationCorrectionTool& original, std::array<std::size_t, 4>& counts) {
    constexpr bool photon = std::is_same_v<Object, xAOD::Photon>;
    const auto candidate = provider->evaluate_array(raw_input(raw, ev, pv));
    const bool selected = stock.gate_eta[i] == pass && stock.oq[i] == pass && stock.track_or_clean[i] == pass &&
        stock.pt_gate[i] == pass && stock.extra_pt[i] == pass && stock.id[i] && stock.isolation[i];

    const auto p4_ok = exact(candidate[0], stock.pt[i]) && exact(candidate[1], stock.eta[i]) &&
        // Electron and photon m() are release accessor values. The observer's
        // unavailable m auxiliary branches are zero-filled by its writer.
        exact(candidate[2], stock.phi[i]) && exact(candidate[3], raw.m());
    if (!p4_ok || candidate[4] != stock.gate_eta[i] || candidate[5] != stock.oq[i] ||
        candidate[photon ? 6 : 8] != stock.track_or_clean[i] || candidate[7] != stock.pt_gate[i] ||
        candidate[9] != (stock.id[i] != 0) || candidate[10] != (stock.isolation[i] != 0) ||
        candidate[11] != stock.extra_pt[i] || candidate[12] != selected) {
        std::cerr << "entry=" << entry << " " << (photon ? "photon" : "electron") << "_index=" << i
                  << " p4_exact=" << p4_ok << " raw_pt=" << raw.pt() << " pt=" << candidate[0] << '/' << stock.pt[i]
                  << " eta=" << candidate[4] << '/' << stock.gate_eta[i] << " oq=" << candidate[5] << '/' << stock.oq[i]
                  << " track_or_clean=" << candidate[photon ? 6 : 8] << '/' << stock.track_or_clean[i]
                  << " pt_gate=" << candidate[7] << '/' << stock.pt_gate[i]
                  << " id=" << candidate[9] << '/' << int(stock.id[i])
                  << " iso=" << candidate[10] << '/' << int(stock.isolation[i])
                  << " extra_pt=" << candidate[11] << '/' << stock.extra_pt[i] << '\n';
        throw std::runtime_error("Stock nominal egamma gate or momentum differs");
    }

    Object oracle;
    oracle.makePrivateStore(raw);
    if constexpr (photon) {
        const bool reached_origin = stock.gate_eta[i] == pass && stock.oq[i] == pass && stock.track_or_clean[i] == pass && pv;
        require(candidate[20] == reached_origin, "Photon origin application status");
        if (reached_origin) {
            photonWrtPoint::correctForZ(oracle, pv->z());
            require(exact(oracle.pt(), candidate[0]) && exact(oracle.eta(), candidate[1]), "Photon origin numeric");
        }
    }
    unsigned corrected = 0;
    if (stock.pt_gate[i] == pass) {
        for (std::size_t c = 0; c < cones.size(); ++c) {
            float iso = 0., leakage = 0.;
            if (raw.isolationValue(iso, cones[c]) && raw.isolationCaloCorrection(leakage, cones[c], xAOD::Iso::ptCorrection))
                corrected |= 1u << c;
        }
        require(original.applyCorrection(oracle) == CP::CorrectionCode::Ok, "Stock egamma isolation status");
    }
    require(candidate[19] == corrected, "Egamma correction mask");
    for (std::size_t c = 0; c < cones.size(); ++c) {
        if (!(corrected & (1u << c))) continue;
        float iso = 0., leakage = 0.;
        require(oracle.isolationValue(iso, cones[c]) &&
                oracle.isolationCaloCorrection(leakage, cones[c], xAOD::Iso::ptCorrection), "Read Stock iso/decorated leakage");
        if (!exact(iso, candidate[13+c]) || !exact(leakage, candidate[16+c])) {
            std::cerr << "entry=" << entry << " " << (photon ? "photon" : "electron") << "_index=" << i
                      << " cone=" << c << " iso_bits=" << std::bit_cast<unsigned>(iso) << '/'
                      << std::bit_cast<unsigned>(float(candidate[13+c]))
                      << " leak_bits=" << std::bit_cast<unsigned>(leakage) << '/'
                      << std::bit_cast<unsigned>(float(candidate[16+c])) << '\n';
            throw std::runtime_error("Stock egamma leakage or timing correction differs");
        }
    }
    ++counts[0];
    counts[1] += selected;
    counts[2] += std::popcount(corrected);
    counts[3] += candidate[20] != 0;
}

void configure(CP::IsolationCorrectionTool& tool, bool photon) {
    auto prop = [&](const char* key, const auto& value) {
        require(tool.setProperty(key, value).isSuccess(), key);
    };
    prop("ToolVer", std::string(photon ? "REL21" : "REL22"));
    prop("CorrFile", std::string(photon ? "IsolationCorrections/v5/isolation_ptcorrections_rel20_2.root" :
                                      "IsolationCorrections/v6/isolation_ptcorrections_rel22_mc20.root"));
    prop("IsMC", true);
    prop("Trouble_categories", true);
    prop("FixTimingIssueInCore", true);
    prop("Apply_ddshifts", false);
    prop("Apply_SC_leakcorr", false);
    require(tool.initialize().isSuccess(), "Initialize stock egamma isolation tool");
    require(tool.applySystematicVariation(CP::SystematicSet{}).isSuccess(), "Initialize NOSYS isolation");
}
} // namespace

std::size_t egamma_boundaries(const correction::CorrectionSet&);

int main(int argc, char** argv) try {
    require(argc == 6, "Usage: met_egamma_prepare_stock_check PLUGIN EGAMMA_JSON EGAMMA_PROBES INPUT OBSERVER");
    correction::PluginRegistry registry;
    registry.load_library(argv[1]);
    const auto set = correction::CorrectionSet::from_file(argv[2], registry);
    const auto table_probes = compare_tables(*set, argv[3]);
    const auto electron = set->at("electron_prepare");
    const auto photon = set->at("photon_prepare");
    std::cout << "egamma_root_probe_rows=" << table_probes
              << " egamma_composed_boundary_probes=" << egamma_boundaries(*set) << '\n';

    require(xAOD::Init("met_egamma_prepare_stock_check").isSuccess(), "xAOD init");
    std::unique_ptr<TFile> original(TFile::Open(argv[4], "READ"));
    std::unique_ptr<TFile> observer(TFile::Open(argv[5], "READ"));
    require(original && !original->IsZombie() && observer && !observer->IsZombie(), "Open original and reference");
    xAOD::TEvent reader(xAOD::TEvent::kClassAccess);
    require(reader.readFrom(original.get()).isSuccess(), "Read original PHYSLITE");
    TTreeReader reference("analysis", observer.get());
    Species e_ref(reference, false), p_ref(reference, true);
    CP::IsolationCorrectionTool e_iso("met_stock_electron_isolation"), p_iso("met_stock_photon_isolation");
    configure(e_iso, false);
    configure(p_iso, true);
    std::array<std::size_t, 4> e_counts{}, p_counts{};

    Long64_t entry = 0;
    while (reference.Next()) {
        require(reader.getEntry(entry) >= 0, "Read original event");
        const xAOD::ElectronContainer* electrons = nullptr;
        const xAOD::PhotonContainer* photons = nullptr;
        const xAOD::VertexContainer* vertices = nullptr;
        const xAOD::EventInfo* event = nullptr;
        require(reader.retrieve(electrons, "AnalysisElectrons").isSuccess(), "Original electrons");
        require(reader.retrieve(photons, "AnalysisPhotons").isSuccess(), "Original photons");
        require(reader.retrieve(vertices, "PrimaryVertices").isSuccess(), "Original PVs");
        require(reader.retrieve(event, "EventInfo").isSuccess() && event->eventType(xAOD::EventInfo::IS_SIMULATION),
                "Original full simulation EventInfo");
        const xAOD::Vertex* pv = nullptr;
        for (const auto* raw : *vertices) {
            if (raw->vertexType() == xAOD::VxType::PriVtx) { pv = raw; break; }
        }
        require(e_ref.pt.GetSize() == electrons->size() && p_ref.pt.GetSize() == photons->size(),
                "Expected every stock egamma row");
        for (std::size_t i = 0; i < electrons->size(); ++i)
            compare_row(entry, i, *electrons->at(i), *event, pv, electron, e_ref, e_iso, e_counts);
        for (std::size_t i = 0; i < photons->size(); ++i)
            compare_row(entry, i, *photons->at(i), *event, pv, photon, p_ref, p_iso, p_counts);
        ++entry;
    }
    require(entry == reader.getEntries() && e_counts[0] && p_counts[0], "Expected all original events and egamma rows");
    std::cout << "entries=" << entry << " electron_rows=" << e_counts[0] << " electron_selected=" << e_counts[1]
              << " electron_iso_cones=" << e_counts[2] << " photon_rows=" << p_counts[0]
              << " photon_selected=" << p_counts[1] << " photon_iso_cones=" << p_counts[2]
              << " photon_reached_origin=" << p_counts[3] << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
