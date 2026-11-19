// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include <correction_plugin.h>

#include <TauAnalysisTools/TauSelectionTool.h>
#include <TauAnalysisTools/TauSmearingTool.h>
#include <TFile.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODJet/JetContainer.h>
#include <xAODRootAccess/Init.h>
#include <xAODRootAccess/TEvent.h>
#include <xAODTau/TauJetContainer.h>
#include <xAODTau/TauTrack.h>
#include <xAODTau/TauxAODHelpers.h>
#include <xAODTruth/TruthParticle.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Int = std::int64_t;
using Ref = correction::Correction::Ref;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class Tool, class T>
void setting(Tool& tool, const char* name, const T& value) {
    require(tool.setProperty(name, value).isSuccess(), name);
}

bool exact(float actual, float expected) {
    return (std::isnan(actual) && std::isnan(expected))
        || std::bit_cast<std::uint32_t>(actual) == std::bit_cast<std::uint32_t>(expected);
}

std::vector<correction::Variable::Type> input(const xAOD::TauJet& tau, const Ref& track) {
    static const SG::ConstAccessor<char> compatibility{"TESCompatibility"};
    static const SG::ConstAccessor<float> tes_pt{"ptTauEnergyScale"};
    static const SG::ConstAccessor<char> hadronic{"IsHadronicTau"};
    static const SG::ConstAccessor<std::uint32_t> flags{"isTauFlags"};
    static const SG::ConstAccessor<ElementLink<xAOD::JetContainer>> truthjet{"truthJetLink"};

    Int count = 0;
    for (const auto& link : tau.allTauTrackLinks()) {
        if (!link.isValid()) continue;
        count += static_cast<Int>(track->evaluate({static_cast<Int>((*link)->flagSet())}));
    }
    require(count == static_cast<Int>(tau.nTracksCharged()), "Charged tau track count");
    int mode = -1;
    tau.panTauDetail(xAOD::TauJetParameters::PanTau_DecayMode, mode);
    const auto* truth = xAOD::TauHelpers::getTruthParticle(&tau);
    const Int truth_pdg = truth ? truth->pdgId() : 0;
    const Int truth_hadronic = truth && truth->isTau() && hadronic.isAvailable(*truth) ? hadronic(*truth) : -1;
    const bool has_tes = tes_pt.isAvailable(tau);
    return {tau.pt(), tau.eta(), tau.phi(), tau.m(), double(tau.charge()), tau.ptFinalCalib(), tau.etaFinalCalib(),
            has_tes ? tau.ptTauEnergyScale() : 0., has_tes ? tau.etaTauEnergyScale() : 0.,
            count, flags.isAvailable(tau) ? static_cast<Int>(flags(tau)) : Int{-1},
            static_cast<Int>(mode), static_cast<Int>(has_tes),
            static_cast<Int>(compatibility.isAvailable(tau) ? compatibility(tau) : -1), truth_pdg,
            truth_hadronic, static_cast<Int>(truthjet.isAvailable(tau))};
}

} // namespace

int main(int argc, char** argv) try {
    if (argc != 4) throw std::invalid_argument("Usage: met_tau_prepare_stock_check PLUGIN TAU_JSON INPUT_PHYSLITE");
    require(xAOD::Init("met_tau_prepare_stock_check").isSuccess(), "xAOD initialization");
    std::unique_ptr<TFile> source(TFile::Open(argv[3], "READ"));
    require(source && !source->IsZombie(), "Cannot open stock input file");
    xAOD::TEvent reader(xAOD::TEvent::kClassAccess);
    require(reader.readFrom(source.get()).isSuccess(), "Read original PHYSLITE");

    correction::PluginRegistry registry;
    registry.load_library(argv[1]);
    const auto cset = correction::CorrectionSet::from_file(argv[2], registry);
    const auto prepare = cset->at("tau_prepare");
    const auto track = cset->at("tau_track_charged");

    TauAnalysisTools::TauSelectionTool selection("tau_prepare_stock_selection");
    setting(selection, "ConfigPath", std::string("TauAnalysisAlgorithms/tau_selection_tight_noeleid.conf"));
    setting(selection, "TauContainerName", std::string("AnalysisTauJets"));
    require(selection.initialize().isSuccess(), "Initialize tau Tight oracle");

    TauAnalysisTools::TauSmearingTool smearing("tau_prepare_stock_smearing");
    setting(smearing, "Campaign", std::string("mc23"));
    setting(smearing, "RecommendationTag", std::string("2025-prerec"));
    setting(smearing, "useFastSim", false);
    setting(smearing, "useGNTau", false);
    setting(smearing, "ApplyInsituCorrection", true);
    setting(smearing, "ApplyFading", true);
    setting(smearing, "SkipTruthMatchCheck", false);
    setting(smearing, "MVATESQualityCheck", true);
    setting(smearing, "InputFilePath", std::string("TauAnalysisTools/00-04-00/Smearing/TES_TrueHadTau_RNN_mc23_v3.root"));
    require(smearing.initialize().isSuccess(), "Initialize tau nominal oracle");
    require(smearing.applySystematicVariation(CP::SystematicSet()).isSuccess(), "Set stock tau NOSYS");

    constexpr std::array cut_names{"Pt", "AbsEta", "NTrack", "AbsCharge", "JetIDWP", "EleIDWP"};
    static const SG::ConstAccessor<char> compatibility{"TESCompatibility"};
    std::array<std::size_t, 3> smear_counts{};
    std::array<std::size_t, cut_names.size()> gate_counts{};
    std::size_t selected = 0, rows = 0;
    for (Long64_t entry = 0; entry < reader.getEntries(); ++entry) {
        require(reader.getEntry(entry) >= 0, "Read original event");
        const xAOD::EventInfo* info = nullptr;
        require(reader.retrieve(info, "EventInfo").isSuccess() && info
                && info->eventType(xAOD::EventInfo::IS_SIMULATION), "Expected simulation EventInfo");
        const xAOD::TauJetContainer* taus = nullptr;
        require(reader.retrieve(taus, "AnalysisTauJets").isSuccess() && taus, "Retrieve original AnalysisTauJets");
        for (const auto* raw : *taus) {
            const auto predicted = prepare->evaluate_array(input(*raw, track));
            xAOD::TauJet varied;
            varied.makePrivateStore(*raw);
            require(smearing.applyCorrection(varied) == CP::CorrectionCode::Ok, "Stock tau correction result");
            bool momentum_ok = true;
            const std::array stock_p4{varied.pt(), varied.eta(), varied.phi(), varied.m()};
            for (std::size_t i = 0; i < stock_p4.size(); ++i) {
                momentum_ok = momentum_ok && exact(static_cast<float>(predicted.at(i)), static_cast<float>(stock_p4[i]));
            }
            const auto accept = selection.accept(varied);
            unsigned mask = 0;
            for (std::size_t i = 0; i < cut_names.size(); ++i) {
                if (accept.getCutResult(cut_names[i])) {
                    ++gate_counts[i];
                    mask |= 1u << i;
                }
            }
            const bool selected_tau = static_cast<bool>(accept);
            if (!momentum_ok || predicted.at(4) != selected_tau || predicted.at(5) != mask
                || !compatibility.isAvailable(varied) || predicted.at(6) != (compatibility(varied) != 0)) {
                std::cerr << "entry=" << entry << " tau_index=" << raw->index() << " p4_agrees=" << momentum_ok
                          << " original_bits=" << mask << " provider_bits=" << predicted.at(5)
                          << " stock_selected=" << selected_tau << " provider_selected=" << predicted.at(4) << '\n';
                return 1;
            }
            ++smear_counts.at(static_cast<std::size_t>(predicted.at(7)));
            selected += selected_tau;
            ++rows;
        }
    }
    require(rows != 0, "No stock tau rows");
    std::cout << "rows=" << rows << " selected=" << selected << " smear_branch_counts="
              << smear_counts[0] << ',' << smear_counts[1] << ',' << smear_counts[2] << '\n';
    for (std::size_t i = 0; i < cut_names.size(); ++i) std::cout << cut_names[i] << '=' << gate_counts[i] << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
