// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "preparation.hpp"

namespace atlas::met::preparation {

void register_tau(correction::PluginRegistry&);
void register_muon_nodes(correction::PluginRegistry&);

void register_preparation_nodes(correction::PluginRegistry& registry) {
    register_tau(registry);
    register_muon_nodes(registry);
}

} // namespace atlas::met::preparation
