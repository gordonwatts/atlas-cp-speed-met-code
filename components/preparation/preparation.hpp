// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_PREPARATION_HPP
#define ATLAS_MET_PREPARATION_HPP

namespace correction { class PluginRegistry; }

namespace atlas::met::preparation {

// Register before loading the generated correction sets. Factories and loaded
// tables are immutable. The shared provider invokes native CorrectionLib even
// when the calling adapter has already resolved the event-local element links.
void register_preparation_nodes(correction::PluginRegistry&);

} // namespace atlas::met::preparation
#endif
