/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 *
 * Focused reference calls into the runtime owner's unmodified AnalysisBase.
 * Comparisons here supplement the event-level release observer.
 */
#include "significance.hpp"

#include <correction_plugin.h>
#include <ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h>
#include <JetCalibTools/JetCalibrationTool.h>
#include <METUtilities/METSignificance.h>
#include <MuonMomentumCorrections/MuonCalibTool.h>
#include <MuonMomentumCorrections/MuonCalibIntScaleSmearTool.h>
#include <TauAnalysisTools/TauResolutionTool.h>
#include <TFile.h>
#include <xAODCaloEvent/CaloCluster.h>
#include <xAODEgamma/ElectronContainer.h>
#include <xAODEgamma/PhotonContainer.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODJet/JetContainer.h>
#include <xAODMissingET/MissingETAuxContainer.h>
#include <xAODMissingET/MissingETContainer.h>
#include <xAODMuon/MuonContainer.h>
#include <xAODRootAccess/Init.h>
#include <xAODRootAccess/TEvent.h>
#include <xAODTau/TauJetContainer.h>
#include <xAODTracking/TrackParticle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>

namespace m = atlas::met_significance;
using ILink = ElementLink<xAOD::IParticleContainer>;

namespace {

struct Counts {
    std::size_t tests = 0;
    std::size_t failures = 0;
    double max_absolute = 0;
};

struct Check {
    std::map<std::string, Counts> counts;

    void same(std::string_view stage, double stock, double candidate,
              double absolute = 1.e-14, double relative = 5.e-13) {
        auto& c = counts[std::string(stage)];
        ++c.tests;
        if (stock == candidate || (std::isnan(stock) && std::isnan(candidate))) return;
        if (std::isfinite(stock) && std::isfinite(candidate)) {
            const double difference = std::abs(stock - candidate);
            c.max_absolute = std::max(c.max_absolute, difference);
            if (difference <= absolute + relative * std::abs(stock)) return;
        }
        if (++c.failures <= 4) {
            std::cerr << stage << " case " << c.tests << " stock=" << stock
                      << " provider=" << candidate << '\n';
        }
    }

    int report() const {
        std::size_t failed = 0;
        for (const auto& [name, c] : counts) {
            std::cout << name << '\t' << c.tests << '\t' << c.failures << '\t' << c.max_absolute << '\n';
            failed += c.failures;
        }
        return failed ? 1 : 0;
    }
};

void require(bool match, std::string_view message) {
    if (!match) throw std::runtime_error(std::string(message));
}

template<class Tool, class Value>
void property(Tool& tool, const char* name, const Value& value) {
    require(tool.setProperty(name, value).isSuccess(), name);
}

constexpr auto jes_file = "JES_data2017_2016_2015_Recommendation_PFlow_Aug2018_rel21.config";
constexpr auto jes_sequence = "JetArea_Residual_EtaJES_GSC_Smear";
constexpr auto tau_file = "tauRecTools/R22_preprod/CombinedTES_R22_Round2.5_v2.root";

struct StockTools {
    JetCalibrationTool jet{"metsig_jet_witness"};
    CP::EgammaCalibrationAndSmearingTool egamma{"metsig_egamma_witness"};
    CP::MuonCalibTool muon{"metsig_muon_witness"};
    TauAnalysisTools::TauResolutionTool tau{"metsig_tau_witness"};
    met::METSignificance significance{"metsig_event_witness"};

    StockTools() {
        property(jet, "ConfigFile", std::string(jes_file));
        property(jet, "CalibSequence", std::string(jes_sequence));
        property(jet, "CalibArea", std::string("00-04-81"));
        property(jet, "JetCollection", std::string("AntiKt4EMPFlow"));
        property(jet, "IsData", false);
        property(egamma, "ESModel", std::string("es2024_Run3_v0"));
        property(egamma, "decorrelationModel", std::string("1NP_v1"));
        property(egamma, "useFastSim", 0);
        property(egamma, "useMVACalibration", 0);
        property(muon, "calibMode", 0);
        property(tau, "WeightFileName", std::string(tau_file));
        property(significance, "JetCalibConfig", std::string(jes_file));
        property(significance, "JetCalibSequence", std::string(jes_sequence));
        property(significance, "JetCalibArea", std::string("00-04-81"));
        property(significance, "TreatPUJets", false);
        property(significance, "MuonCalibMode", 0);
        property(significance, "EgammaESModel", std::string("es2024_Run3_v0"));
        property(significance, "EgammaDecorrelationModel", std::string("1NP_v1"));
        property(significance, "TauTESConfig", std::string(tau_file));
        require(jet.initialize().isSuccess(), "jet initialize");
        require(egamma.initialize().isSuccess(), "egamma initialize");
        require(muon.initialize().isSuccess(), "muon initialize");
        require(tau.initialize().isSuccess(), "tau initialize");
        require(significance.initialize().isSuccess(), "significance initialize");
    }
};

std::string load_text(const std::filesystem::path& path) {
    std::ifstream file(path);
    require(file.is_open(), path.string());
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

template<class Container>
const Container& get(xAOD::TEvent& reader, const char* key) {
    const Container* result = nullptr;
    require(reader.retrieve(result, key).isSuccess() && result, key);
    return *result;
}

std::array<float, 3> neighbors(float center) {
    constexpr auto inf = std::numeric_limits<float>::infinity();
    return {std::nextafter(center, -inf), center, std::nextafter(center, inf)};
}

double check_jet(Check& check, StockTools& stock, const correction::Correction& provider,
                 const xAOD::Jet& jet) {
    double a = 0;
    double b = 0;
    require(stock.jet.getNominalResolutionData(jet, a).isSuccess(), "stock JER data query");
    require(stock.jet.getNominalResolutionMC(jet, b).isSuccess(), "stock JER MC query");
    const auto out = provider.evaluate_array({jet.pt(), jet.eta()});
    check.same("jet.data", a, out.at(0));
    check.same("jet.mc", b, out.at(1));
    check.same("jet.relative_float", static_cast<float>(std::max(a, b)), out.at(2), 0., 0.);
    return out.at(2);
}

void jet_edges(Check& check, StockTools& stock, const correction::Correction& provider,
               const m::JetData& tables) {
    xAOD::Jet jet;
    jet.makePrivateStore();
    const auto compare = [&](float pt, float eta) {
        jet.setJetP4(xAOD::JetFourMom_t(pt, eta, 0.5, 0.));
        check_jet(check, stock, provider, jet);
    };
    for (const auto* table : {&tables.data, &tables.mc}) {
        for (const double edge : table->abs_eta.edges) {
            for (const float eta : neighbors(static_cast<float>(edge))) {
                compare(40000, eta);
                compare(40000, -eta);
            }
        }
        for (std::size_t y = 1; y < table->projections.size(); ++y) {
            const float eta = std::midpoint(table->abs_eta.edges[y - 1], table->abs_eta.edges[y]);
            const auto& proj = table->projections[y];
            for (const double edge : proj.histogram.axis.edges) {
                for (const float pt : neighbors(static_cast<float>(edge * 1000))) compare(pt, eta);
            }
            for (const double center : proj.centers) compare(static_cast<float>(center * 1000), eta);
        }
    }
    compare(1.e7, 15.);
}

using Integer = std::int64_t;

void egamma_edges(Check& check, StockTools& stock, const correction::Correction& provider,
                  const m::EgammaData& data) {
    constexpr std::array types{PATCore::ParticleType::Electron, PATCore::ParticleType::UnconvertedPhoton,
                               PATCore::ParticleType::ConvertedPhoton};
    for (std::size_t i = 0; i < types.size(); ++i) {
        const bool with_ct = i != 0;
        for (const double edge : data.eta_edges) {
            for (const double eta : neighbors(static_cast<float>(edge))) {
                constexpr float calo = -1.5;
                for (const double et : {4999.9, 5000., 5000.1, 49999.9, 50000., 50000.1}) {
                    const double energy = et * std::cosh(eta);
                    check.same("egamma.edge", stock.egamma.resolution(energy, eta, calo, types[i], with_ct),
                               provider.evaluate({static_cast<Integer>(i), energy, eta, double(calo),
                                                  Integer{1}, static_cast<Integer>(with_ct)}));
                }
            }
        }
    }
    for (const double edge : data.constant_term.axis.edges) {
        for (const double calo : neighbors(static_cast<float>(edge))) {
            check.same("egamma.ct.edge", stock.egamma.resolution(60000., 0.42, calo, types[1], true),
                       provider.evaluate({Integer{1}, 60000., 0.42, calo, Integer{1}, Integer{1}}));
        }
    }
}

void tau_edges(Check& check, StockTools& stock, const correction::Correction& provider,
               const m::TauData& data) {
    xAOD::TauJet tau;
    tau.makePrivateStore();
    constexpr std::array etas{0.1, 0.5, 1.0, 1.5, 2.0};
    using namespace xAOD::TauJetParameters;
    const auto compare = [&](int mode, double tes, double final, double pt) {
        tau.setPanTauDetail(PanTau_DecayMode, mode);
        require(tau.setP4(TauEnergyScale, pt + 1000., tes, 0, 0), "set tau TES p4");
        require(tau.setP4(FinalCalib, pt, final, 0, 0), "set tau final p4");
        check.same("tau.edge", stock.tau.getMvaEnergyResolution(tau),
                   provider.evaluate({static_cast<Integer>(mode), Integer{1},
                                      tau.etaTauEnergyScale(), tau.etaFinalCalib(), tau.ptFinalCalib()}));
    };
    for (int mode = 0; mode < 5; ++mode) {
        for (std::size_t eta = 0; eta < etas.size(); ++eta) {
            const auto& graph = data.graphs.at(mode).at(eta);
            for (const double x : graph.x) {
                for (const float pt : neighbors(static_cast<float>(x))) compare(mode, 0.2, etas[eta], pt);
            }
            compare(mode, 0.2, etas[eta], 0.);
            compare(mode, 0.2, etas[eta], graph.max_x + 1000000.);
        }
    }
    compare(5, 0.2, 0.2, 50000.);
    compare(0, 3., 0.2, 50000.);
}

void muon_edges(Check& check, const correction::Correction& provider, const m::MuonData& data) {
    CP::MuonCalibIntScaleSmearTool stock("metsig_expected_detector_witness");
    property(stock, "calibMode", 0);
    require(stock.initialize().isSuccess(), "detector expected resolution initialize");
    const auto compare = [&](int detector, double pt, double eta, double phi) {
        check.same("muon.detector.edge." + std::to_string(detector),
                   stock.getExpectedResolution(detector, pt, eta, phi, MCP::DataYear::Data22, true),
                   provider.evaluate({static_cast<Integer>(detector), Integer{410000}, Integer{1},
                                      eta, phi, 40., pt, pt}));
    };
    constexpr double infinity = std::numeric_limits<double>::infinity();
    const auto edges = [&](int detector, const auto& histograms) {
        for (const auto& hist : histograms) {
            for (const double edge : hist.x.edges) {
                for (const double eta : {std::nextafter(edge, -infinity), edge, std::nextafter(edge, infinity)})
                    compare(detector, 43., eta, -2.0);
            }
            for (const double edge : hist.y.edges) {
                for (const double phi : {std::nextafter(edge, -infinity), edge, std::nextafter(edge, infinity)})
                    compare(detector, 43., 1.1, phi);
            }
        }
        compare(detector, 0., 0.5, 2.);
    };
    edges(1, data.ms_expected);
    edges(1, data.ms_nominal);
    edges(2, data.id_expected);
    edges(2, data.id_nominal);
}

std::optional<m::TrackParameters> track(const xAOD::Muon& muon, xAOD::Muon::TrackParticleType which) {
    const auto* input = muon.trackParticle(which);
    if (!input) return std::nullopt;
    return m::TrackParameters{input->qOverP(), input->theta(), static_cast<float>(input->phi())};
}

double check_muon(Check& check, StockTools& stock, const correction::Correction& provider,
                  const xAOD::Muon& muon, const xAOD::EventInfo& event) {
    const auto cb = m::track_pt_gev(track(muon, xAOD::Muon::CombinedTrackParticle));
    const auto id = m::track_pt_gev(track(muon, xAOD::Muon::InnerDetectorTrackParticle));
    const auto ms = m::track_pt_gev(track(muon, xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle));
    const auto detector = m::muon_detector(muon.muonType());
    const auto is_sim = static_cast<Integer>(event.eventType(xAOD::EventInfo::IS_SIMULATION));
    double relative = 0;
    for (const Integer det : {1, 2, 3}) {
        const double candidate = provider.evaluate({det, static_cast<Integer>(event.runNumber()), is_sim,
                                                   muon.eta(), muon.phi(), cb, id, ms});
        check.same("muon.detector." + std::to_string(det),
                   stock.muon.expectedResolution(static_cast<int>(det), muon, true), candidate);
        if (det == static_cast<Integer>(detector)) relative = candidate;
    }
    return relative;
}

double check_tau(Check& check, StockTools& stock, const correction::Correction& provider,
                 const xAOD::TauJet& tau) {
    int mode = 99;
    const bool found = tau.panTauDetail(xAOD::TauJetParameters::PanTau_DecayMode, mode);
    const double candidate = provider.evaluate({static_cast<Integer>(mode), static_cast<Integer>(found),
                                                tau.etaTauEnergyScale(), tau.etaFinalCalib(), tau.ptFinalCalib()});
    check.same("tau.real", stock.tau.getMvaEnergyResolution(tau), candidate);
    return candidate;
}

double check_egamma(Check& check, StockTools& stock, const correction::Correction& provider,
                    const xAOD::Egamma& object, bool is_photon) {
    const auto* cluster = object.caloCluster();
    require(cluster, "Egamma cluster link unavailable in witness input");
    const bool forward = object.author() == xAOD::EgammaParameters::AuthorFwdElectron;
    double moment = 0;
    static const SG::ConstAccessor<float> aux("etaCalo");
    const auto moment_value = cluster->retrieveMoment(xAOD::CaloCluster::ETACALOFRAME, moment)
        ? std::optional(moment) : std::nullopt;
    const auto aux_value = aux.isAvailable(*cluster) ? std::optional(aux(*cluster)) : std::nullopt;
    const float calo = (is_photon || forward || moment_value || aux_value)
        ? m::calo_frame_eta(forward, cluster->eta(), moment_value, aux_value) : 0.;
    const auto kind = is_photon ? m::photon_kind(static_cast<const xAOD::Photon&>(object).nVertices(), forward)
                                : m::EgammaKind::electron;
    const double oracle = is_photon
        ? stock.egamma.getResolution(object)
        : stock.egamma.resolution(object.e(), cluster->eta(), calo, PATCore::ParticleType::Electron, false);
    const double candidate = provider.evaluate({static_cast<Integer>(kind), object.e(), cluster->eta(),
                                                double(calo), Integer{1}, static_cast<Integer>(is_photon)});
    check.same(is_photon ? "photon.real" : "electron.real", oracle, candidate);
    return candidate;
}

struct Linked {
    const xAOD::IParticle* object;
    m::ObjectKind kind;
    float resolution;
};

void compare_event(Check& check, StockTools& stock, const xAOD::EventInfo& info,
                   const std::vector<Linked>& objects, bool with_soft, m::State& state) {
    xAOD::MissingETAuxContainer aux;
    xAOD::MissingETContainer stock_terms;
    stock_terms.setStore(&aux);
    const auto add = [&](const m::Term& term) -> xAOD::MissingET& {
        auto* p = new xAOD::MissingET();
        stock_terms.push_back(p);
        p->setName(std::string(term.name));
        p->setSource(term.source);
        p->setMpx(term.momentum.mpx_mev);
        p->setMpy(term.momentum.mpy_mev);
        p->setSumet(term.momentum.sumet_mev);
        return *p;
    };
    static const SG::Accessor<std::vector<ILink>> stock_links("ConstitObjectLinks");
    std::vector<m::Term> terms;
    std::vector<m::Constituent> constituents;
    const auto append = [&](std::string_view name, unsigned bit, m::ObjectKind kind) {
        terms.push_back({name, 0x10000u | bit, {0, 0, 0}});
        auto& links = stock_links(add(terms.back()));
        const auto begin = constituents.size();
        for (const auto& entry : objects) {
            if (entry.kind != kind) continue;
            const auto* container = static_cast<const xAOD::IParticleContainer*>(entry.object->container());
            links.emplace_back(*container, entry.object->index());
            const auto p4 = entry.object->p4();
            constituents.push_back({kind, entry.object->pt(), p4.Px(), p4.Py(), entry.resolution});
        }
        terms.back().links = std::span<const m::Constituent>(constituents).subspan(begin);
    };
    constituents.reserve(objects.size()); // Keep the span for each earlier term valid.
    append("RefEle", 0x1, m::ObjectKind::electron);
    append("RefGamma", 0x2, m::ObjectKind::photon);
    append("RefTau", 0x4, m::ObjectKind::tau);
    append("Muons", 0x8, m::ObjectKind::muon);
    append("RefJet", 0x10, m::ObjectKind::jet);
    if (with_soft) {
        terms.push_back({"PVSoftTrk", 0x14020, {-3200, 13900, 30000}});
        add(terms.back());
    }
    terms.push_back({"Final", 0x10040, {23100, -32800, 210000}});
    add(terms.back());
    require(stock.significance.varianceMET(&stock_terms, info.averageInteractionsPerCrossing(),
                                           "RefJet", "PVSoftTrk", "Final").isSuccess(), "stock significance");
    const auto got = m::reduce(terms, state);
    require(got.status == m::Status::success, "candidate significance status");
    check.same("event.var_l", stock.significance.GetVarL(), got.covariance.longitudinal);
    check.same("event.var_t", stock.significance.GetVarT(), got.covariance.transverse);
    check.same("event.rho", stock.significance.GetRho(), state.rho);
    const auto decoration = [&](std::string_view name, double a, float b) {
        check.same(name, static_cast<float>(a), b, 1.e-6, 5.e-7);
    };
    decoration("event.significance", stock.significance.GetSignificance(), got.decorations.significance);
    decoration("event.directional", stock.significance.GetSigDirectional(), got.decorations.directional);
    decoration("event.met_sqrt_sumet", stock.significance.GetMETOverSqrtSumET(), got.decorations.met_over_sqrt_sumet);
    decoration("event.met_sqrt_ht", stock.significance.GetMETOverSqrtHT(), got.decorations.met_over_sqrt_ht);
}

} // namespace

int main(int argc, char** argv) try {
    if (argc != 5) throw std::invalid_argument("Usage: met_significance_stock_compare PROVIDER BUNDLE INPUT EVENTS");
    const auto requested = std::stoll(argv[4]);
    require(requested > 0, "EVENTS must be positive");
    std::cout << std::setprecision(17);
    std::cerr << std::setprecision(17);
    require(xAOD::Init("met_significance_stock_compare").isSuccess(), "xAOD initialization");
    std::unique_ptr<TFile> source(TFile::Open(argv[3], "READ"));
    require(source && !source->IsZombie(), "input ROOT file");
    xAOD::TEvent reader(xAOD::TEvent::kClassAccess);
    require(reader.readFrom(source.get()).isSuccess(), "read input");

    correction::PluginRegistry registry;
    registry.load_library(argv[1]);
    const std::filesystem::path bundle{argv[2]};
    auto correction_set = correction::CorrectionSet::from_file((bundle / "corrections.json").string(), registry);
    const auto jet = correction_set->at("jet_resolution");
    const auto egamma = correction_set->at("egamma_resolution");
    const auto muon = correction_set->at("muon_resolution");
    const auto tau = correction_set->at("tau_resolution");

    StockTools stock;
    Check check;
    jet_edges(check, stock, *jet, m::read_jet_data(load_text(bundle / "jet.table")));
    egamma_edges(check, stock, *egamma, m::read_egamma_data(load_text(bundle / "egamma.table")));
    tau_edges(check, stock, *tau, m::read_tau_data(load_text(bundle / "tau.table")));
    muon_edges(check, *muon, m::read_muon_data(load_text(bundle / "muon.table")));

    m::State state;
    for (Long64_t entry = 0; entry < std::min<Long64_t>(requested, reader.getEntries()); ++entry) {
        require(reader.getEntry(entry) >= 0, "read event entry");
        const auto& info = get<xAOD::EventInfo>(reader, "EventInfo");
        std::vector<Linked> objects;
        const auto add = [&](const xAOD::IParticle* object, m::ObjectKind kind, double res) {
            objects.push_back({object, kind, static_cast<float>(res)});
        };
        for (auto* obj : get<xAOD::JetContainer>(reader, "AnalysisJets"))
            add(obj, m::ObjectKind::jet, check_jet(check, stock, *jet, *obj));
        for (auto* obj : get<xAOD::ElectronContainer>(reader, "AnalysisElectrons"))
            add(obj, m::ObjectKind::electron, check_egamma(check, stock, *egamma, *obj, false));
        for (auto* obj : get<xAOD::PhotonContainer>(reader, "AnalysisPhotons"))
            add(obj, m::ObjectKind::photon, check_egamma(check, stock, *egamma, *obj, true));
        for (auto* obj : get<xAOD::MuonContainer>(reader, "AnalysisMuons"))
            add(obj, m::ObjectKind::muon, check_muon(check, stock, *muon, *obj, info));
        for (auto* obj : get<xAOD::TauJetContainer>(reader, "AnalysisTauJets"))
            add(obj, m::ObjectKind::tau, check_tau(check, stock, *tau, *obj));
        compare_event(check, stock, info, objects, true, state);
        if (entry == 0) compare_event(check, stock, info, {}, false, state);
    }
    for (const auto* stage : {"electron.real", "photon.real", "muon.detector.3", "tau.real", "event.significance"}) {
        require(check.counts.contains(stage), std::string("No witness rows for ") + stage);
    }
    return check.report();
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
