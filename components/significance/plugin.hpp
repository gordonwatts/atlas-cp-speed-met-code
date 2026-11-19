/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
#ifndef ATLAS_MET_SIGNIFICANCE_PLUGIN_HPP
#define ATLAS_MET_SIGNIFICANCE_PLUGIN_HPP

namespace correction { class PluginRegistry; }

namespace atlas::met_significance {
void register_significance_nodes(correction::PluginRegistry&);
}

#endif
