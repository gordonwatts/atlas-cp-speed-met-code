// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "bridge.hpp"
#include "maker_xaod.hpp"
#include "met_adapter/significance_correction.hpp"

#include <AsgDataHandles/ReadHandle.h>
#include <MetAnalysisAlgorithms/MetSignificanceAlg.h>
#include <xAODCaloEvent/CaloCluster.h>
#include <xAODEgamma/Egamma.h>
#include <xAODEgamma/Photon.h>
#include <xAODMuon/Muon.h>
#include <xAODTau/TauJet.h>
#include <xAODTracking/TrackParticle.h>

namespace f = atlas::met::framework;
namespace a = atlas::met::adapter;
namespace s = atlas::met_significance;

namespace {
thread_local std::map<const CP::MetSignificanceAlg*, s::State> state;
const SG::ConstAccessor<std::vector<f::ParticleLink>> constituents{"ConstitObjectLinks"};
const SG::ConstAccessor<float> eta_calo{"etaCalo"};
const SG::ConstAccessor<float> tes_eta{"etaTauEnergyScale"};
const SG::ConstAccessor<float> final_eta{"etaFinalCalib"};
const SG::ConstAccessor<float> final_pt{"ptFinalCalib"};

const correction::Correction& evaluator()
{
    struct Loaded {
        correction::PluginRegistry registry;
        std::unique_ptr<correction::CorrectionSet> set;
        Loaded()
        {
            registry.load_library(f::absolute_file("ATLAS_MET_PROVIDER"));
            set = correction::CorrectionSet::from_file(f::absolute_file("ATLAS_MET_SIGNIFICANCE_PAYLOAD"), registry);
        }
    };
    static const Loaded loaded;
    return *loaded.set->at(a::significance_correction_name);
}

std::optional<s::TrackParameters> track(const xAOD::TrackParticle* row)
{
    if (!row) return std::nullopt;
    return s::TrackParameters{row->qOverP(), row->theta(), row->phi()};
}

template <class T>
std::optional<T> optional_aux(const SG::ConstAccessor<T>& accessor, const SG::AuxElement& object)
{
    if (!accessor.isAvailable(object)) return std::nullopt;
    return accessor(object);
}

a::EgammaInput egamma(const xAOD::Egamma& object)
{
    a::EgammaInput row{};
    row.energy_mev = object.e();
    row.forward_author = object.author() == xAOD::EgammaParameters::AuthorFwdElectron;
    const auto* cluster = object.caloCluster();
    if (cluster) {
        row.cluster_eta = cluster->eta();
        double moment{};
        if (cluster->retrieveMoment(xAOD::CaloCluster::ETACALOFRAME, moment)) row.eta_calo_moment = moment;
        row.eta_calo_aux = optional_aux(eta_calo, *cluster);
    }
    if (object.type() == xAOD::Type::Photon) {
        row.vertex_link_count = static_cast<const xAOD::Photon&>(object).nVertices();
    }
    return row;
}

a::SignificanceLink link(const f::ParticleLink& source)
{
    a::SignificanceLink row{};
    if (!source.isValid()) return row;
    const auto* object = *source;
    if (!object) return row;
    row.valid = true;
    row.pt_mev = object->pt();
    row.eta = object->eta();
    row.phi = object->phi();
    const auto p4 = object->p4();
    row.p4_px_mev = p4.Px();
    row.p4_py_mev = p4.Py();
    switch (object->type()) {
    case xAOD::Type::Jet:
        row.kind = s::ObjectKind::jet;
        break;
    case xAOD::Type::Electron:
        row.kind = s::ObjectKind::electron;
        row.egamma = egamma(static_cast<const xAOD::Egamma&>(*object));
        break;
    case xAOD::Type::Photon:
        row.kind = s::ObjectKind::photon;
        row.egamma = egamma(static_cast<const xAOD::Egamma&>(*object));
        break;
    case xAOD::Type::Muon: {
        row.kind = s::ObjectKind::muon;
        const auto& muon = static_cast<const xAOD::Muon&>(*object);
        row.muon = a::MuonInput{static_cast<std::uint16_t>(muon.muonType()),
            track(muon.trackParticle(xAOD::Muon::CombinedTrackParticle)),
            track(muon.trackParticle(xAOD::Muon::InnerDetectorTrackParticle)),
            track(muon.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle))};
        break;
    }
    case xAOD::Type::Tau: {
        row.kind = s::ObjectKind::tau;
        const auto& tau = static_cast<const xAOD::TauJet&>(*object);
        int mode{};
        std::optional<int> decay;
        if (tau.panTauDetail(xAOD::TauJetParameters::PanTauDetails::PanTau_DecayMode, mode)) decay = mode;
        row.tau = a::TauInput{decay, optional_aux(tes_eta, tau),
                             optional_aux(final_eta, tau), optional_aux(final_pt, tau)};
        break;
    }
    default:
        break;
    }
    return row;
}

a::SignificanceResult run(const xAOD::MissingETContainer& met, const xAOD::EventInfo& info,
                           s::State previous)
{
    std::vector<std::vector<a::SignificanceLink>> link_storage;
    std::vector<a::SignificanceTerm> terms;
    link_storage.reserve(met.size());
    terms.reserve(met.size());
    for (const auto* term : met) {
        if (!term) throw std::runtime_error("Null xAOD MET term");
        const bool available = constituents.isAvailable(*term);
        auto& links = link_storage.emplace_back();
        if (available) {
            for (const auto& source : constituents(*term)) links.push_back(link(source));
        }
        terms.push_back({term->name(), term->source(), {term->mpx(), term->mpy(), term->sumet()}, links, available});
    }
    return a::invoke_significance(evaluator(),
        {terms, previous, info.runNumber(), info.eventType(xAOD::EventInfo::IS_SIMULATION),
         info.averageInteractionsPerCrossing()});
}

void verify_significance(const asg::AsgTool& tool)
{
    for (const char* item : {"DoPhiReso", "TreatPUJets", "ApplyBias", "DoJerForEMu", "IsDataJet",
                             "IsDataMuon", "EgammaUseFastsim"}) {
        f::expect_property(tool, item, false);
    }
    f::expect_property(tool, "SoftTermParam", 0);
    f::expect_property(tool, "SoftTermReso", 8.5);
    f::expect_property(tool, "JetPtThr", -1.F);
    f::expect_property(tool, "JetEtaThr", -1.F);
    f::expect_property(tool, "JetResoAux", std::string{});
    f::expect_property(tool, "EMuResoAux", std::string{});
    f::expect_property(tool, "JetCollection", std::string("AntiKt4EMPFlow"));
    f::expect_property(tool, "JetCalibArea", std::string("00-04-81"));
    f::expect_property(tool, "EgammaESModel", std::string("es2024_Run3_v0"));
    f::expect_property(tool, "TauTESConfig", std::string("tauRecTools/R22_preprod/CombinedTES_R22_Round2.5_v2.root"));
}
} // namespace

StatusCode CP::MetSignificanceAlg::initialize()
{
    using Original = StatusCode (*)(CP::MetSignificanceAlg*);
    static const auto original = f::original_entry<Original>("_ZN2CP18MetSignificanceAlg10initializeEv");
    const auto result = original(this);
    if (result.isSuccess() && name() == f::significance_name) state[this] = s::State{};
    return result;
}

StatusCode CP::MetSignificanceAlg::execute(const EventContext& ctx)
{
    if (name() != f::significance_name || f::stage() == f::Stage::maker_only) {
        return f::original_execute(this, ctx, "_ZN2CP18MetSignificanceAlg7executeERK12EventContext");
    }
    try {
        if (m_totalMETName.value() != "Final" || m_jetTermName.value() != "RefJet"
            || m_softTermName.value() != "PVSoftTrk") {
            throw std::runtime_error("Unsupported initialized MET significance term names");
        }
        verify_significance(f::asg_tool(&*m_significanceTool));
        SG::ReadHandle<xAOD::EventInfo> event_info(m_eventInfoKey, ctx);
        ANA_CHECK(event_info.isValid());
        // One ordered event loop in the pinned CPRun. Distinct tool instances have distinct state.
        auto& previous = state.at(this);
        for (const auto& sys : m_systematicsList.systematicsVector()) {
            xAOD::MissingETContainer* met{};
            ANA_CHECK(m_metHandle.getCopy(met, sys, ctx));
            const auto* total = (*met)[m_totalMETName.value()];
            if (!total) {
                ANA_MSG_ERROR("Missing Final term in MET significance bridge");
                return StatusCode::FAILURE;
            }
            const auto result = run(*met, *event_info.cptr(), previous);
            f::record_dispatch("significance");
            previous = result.next;
            if (result.reduction.status != s::Status::success) {
                ANA_MSG_ERROR("MET significance status " << static_cast<int>(result.reduction.status));
                return StatusCode::FAILURE;
            }
            const auto& values = result.reduction.decorations;
            m_significanceDecorHandle.set(*total, values.significance, sys);
            if (!m_sigDirectionalDecorHandle.empty()) m_sigDirectionalDecorHandle.set(*total, values.directional, sys);
            if (!m_metOverSqrtSumETDecorHandle.empty()) m_metOverSqrtSumETDecorHandle.set(*total, values.met_over_sqrt_sumet, sys);
            if (!m_metOverSqrtHTDecorHandle.empty()) m_metOverSqrtHTDecorHandle.set(*total, values.met_over_sqrt_ht, sys);
        }
    } catch (const std::exception& error) {
        ANA_MSG_ERROR("MET significance bridge: " << error.what());
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}
