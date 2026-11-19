// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "bridge.hpp"
#include "soft_xaod.hpp"

#include <MetAnalysisAlgorithms/MetBuilderAlg.h>

namespace f = atlas::met::framework;

StatusCode CP::MetBuilderAlg::execute(const EventContext& ctx)
{
    if (name() != f::builder_name || f::stage() == f::Stage::maker_only) {
        return f::original_execute(this, ctx, "_ZN2CP13MetBuilderAlg7executeERK12EventContext");
    }
    try {
        if (m_finalKey.value() != "Final" || m_softTerm.value() != "PVSoftTrk") {
            throw std::runtime_error("Unsupported initialized MetBuilderAlg term names");
        }
        for (const auto& sys : m_systematicsList.systematicsVector()) {
            xAOD::MissingETContainer* met{};
            ANA_CHECK(m_metHandle.getCopy(met, sys, ctx));
            if (!(*met)[m_softTerm.value()]) {
                ANA_MSG_ERROR("Missing PVSoftTrk in MET builder bridge");
                return StatusCode::FAILURE;
            }
            ANA_CHECK(f::build_final(*met));
            f::record_dispatch("final");
        }
    } catch (const std::exception& error) {
        ANA_MSG_ERROR("MET builder bridge: " << error.what());
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}
