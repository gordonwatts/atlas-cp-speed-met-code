/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
#include "plugin.hpp"

#include <correction_plugin.h>

extern "C" int correctionlib_plugin_api_version() { return correction::plugin_api_version; }
extern "C" const char* correctionlib_plugin_build_id() { return correction::plugin_build_id().c_str(); }
extern "C" void correctionlib_register_nodes_v4(correction::PluginRegistry& registry) {
    atlas::met_significance::register_significance_nodes(registry);
}
