// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "bridge.hpp"
#include "maker_xaod.hpp"
#include "maker_trace.hpp"
#include "met_adapter/maker_correction.hpp"

#include <MetAnalysisAlgorithms/MetMakerAlg.h>
#include <AsgDataHandles/ReadHandle.h>
#include <xAODMissingET/MissingETAuxContainer.h>
#include <xAODMissingET/MissingETAssociationHelper.h>

#if defined(ATLAS_MET_FULL_BRIDGE)
#include "soft_xaod.hpp"
#endif

namespace f = atlas::met::framework;
namespace m = atlas::met::maker;

namespace {
void verify_maker(const asg::AsgTool& tool)
{
    for (const char* enabled : {"DoPFlow", "DoRemoveElecTrks", "DoRemoveMuonJets",
                                "DoSetMuonJetEMScale", "ORCaloTaggedMuons"}) {
        f::expect_property(tool, enabled, true);
    }
    for (const char* disabled : {"CorrectJetPhi", "DoJetTruth", "DoMuonEloss", "DoRemoveElecTrksEM",
                                 "DoSoftTruth", "GreedyPhotons", "UseGhostMuons", "VeryGreedyPhotons",
                                 "skipSystematicJetSelection"}) {
        f::expect_property(tool, disabled, false);
    }
    f::expect_property(tool, "JetConstitScaleMom", std::string("JetConstitScaleMomentum"));
    f::expect_property(tool, "JetSelection", std::string("Tight"));
    f::expect_property(tool, "JetRejectionDec", std::string{});
    f::expect_property(tool, "JetEtaForw", m::Config::forward_jet_eta);
    f::expect_property(tool, "JetEtaMax", m::Config::max_jet_eta);
    f::expect_property(tool, "JetMinEFrac", m::Config::min_unique_energy_fraction);
    f::expect_property(tool, "JetMinWeightedPt", m::Config::min_weighted_jet_pt);
}

void copy_helper_flags(const m::Result& result, const xAOD::MissingETAssociationMap& map,
                       xAOD::MissingETAssociationHelper& helper)
{
    for (std::size_t row = 0; row < result.flags.size(); ++row) {
        for (std::size_t slot = 0; slot < 32; ++slot) {
            if (result.flags[row] & (m::Mask{1} << slot)) helper.setObjSelectionFlag(map.at(row), slot, true);
        }
    }
}
} // namespace

// The stock initialize method still owns all handle and tool setup.
StatusCode CP::MetMakerAlg::execute(const EventContext& ctx)
{
    if (name() != f::maker_name) {
        return f::original_execute(this, ctx, "_ZN2CP11MetMakerAlg7executeERK12EventContext");
    }
    try {
        const auto selected_stage = f::stage();
        if (m_invisHandles.size() || m_doTrackMet || m_evaluateNNMET || m_switchTauMuonOrder
            || !m_doJetJVT || !m_electronsHandle || !m_photonsHandle || !m_tausHandle || !m_muonsHandle
            || m_electronsKey.value() != "RefEle" || m_photonsKey.value() != "RefGamma"
            || m_tausKey.value() != "RefTau" || m_muonsKey.value() != "Muons"
            || m_jetsKey.value() != "RefJet" || m_softTermKey.value() != "PVSoftTrk") {
            throw std::runtime_error("Unsupported initialized MetMakerAlg settings");
        }
        const auto& tool = f::asg_tool(&*m_makerTool);
        verify_maker(tool);
        const auto& evaluator = f::maker_evaluator();
        SG::ReadHandle<xAOD::MissingETContainer> core(m_metCoreKey, ctx);
        SG::ReadHandle<xAOD::MissingETAssociationMap> association(m_metAssociationKey, ctx);
        SG::ReadHandle<xAOD::VertexContainer> pv("PrimaryVertices");
        ANA_CHECK(core.isValid());
        ANA_CHECK(association.isValid());
        const bool trace = f::maker_trace_enabled();
        const auto trace_entry = trace ? f::next_maker_trace_entry() : 0;
        std::size_t state_index{};
        for (const auto& sys : m_systematicsList.systematicsVector()) {
            f::MakerInputs input;
            auto select = [&](auto& handle, auto& selection, f::MakerInputs::Hard slot) -> StatusCode {
                const xAOD::IParticleContainer* particles{};
                ANA_CHECK(handle.retrieve(particles, sys, ctx));
                for (const auto* particle : *particles) {
                    if (selection.getBool(*particle, sys)) input.selected(slot, *particle);
                }
                return StatusCode::SUCCESS;
            };
            ANA_CHECK(select(m_electronsHandle, m_electronsSelection, f::MakerInputs::electron));
            ANA_CHECK(select(m_photonsHandle, m_photonsSelection, f::MakerInputs::photon));
            ANA_CHECK(select(m_tausHandle, m_tausSelection, f::MakerInputs::tau));
            ANA_CHECK(select(m_muonsHandle, m_muonsSelection, f::MakerInputs::muon));
            const xAOD::JetContainer* jets{};
            ANA_CHECK(m_jetsHandle.retrieve(jets, sys, ctx));
            input.jets(*jets, tool.getKey(jets) != 0);
            input.association_map(*association.cptr());
            input.core(*core.cptr());
            input.vertices(pv.isValid() ? pv.cptr() : nullptr);
            auto view = input.view();
            view.capture_flag_writes = trace;
            const auto result = atlas::met::adapter::invoke_maker(evaluator, view);
            if (trace) f::trace_maker(trace_entry, state_index, view, result, input);
            ++state_index;
            f::record_dispatch("maker");
            if (result.status != m::Status::success) {
                for (const auto& issue : result.diagnostics) {
                    ANA_MSG_ERROR("MET maker issue " << static_cast<int>(issue.issue) << " " << issue.detail);
                }
                return StatusCode::FAILURE;
            }
            auto met = std::make_unique<xAOD::MissingETContainer>();
            auto aux = std::make_unique<xAOD::MissingETAuxContainer>();
            met->setStore(aux.get());
            f::materialize(result, input, *met);

            if (!m_systematicsTool.empty()) {
                ANA_CHECK(m_systematicsTool->applySystematicVariation(sys));
                if (selected_stage == f::Stage::maker_only) {
                    xAOD::MissingETAssociationHelper helper(association.cptr());
                    copy_helper_flags(result, *association.cptr(), helper);
                    auto* soft = (*met)[m_softTermKey.value()];
                    if (!soft) return StatusCode::FAILURE;
                    m_systematicsTool->setRandomSeed(static_cast<int>(1e6 * soft->phi()));
                    ANA_CHECK(m_systematicsTool->applyCorrection(*soft, helper));
                    f::record_dispatch("stock_soft_diagnostic");
                }
#if defined(ATLAS_MET_FULL_BRIDGE)
                else {
                    ANA_CHECK(f::apply_soft(*met, sys.name(), ctx));
                    f::record_dispatch("soft");
                }
#endif
            }
            ANA_CHECK(m_metHandle.record(std::move(met), std::move(aux), sys, ctx));
        }
    } catch (const std::exception& error) {
        ANA_MSG_ERROR("MET maker bridge: " << error.what());
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}
