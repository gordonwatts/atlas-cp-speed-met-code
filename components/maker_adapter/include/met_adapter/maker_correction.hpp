// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_MAKER_CORRECTION_HPP
#define ATLAS_MET_MAKER_CORRECTION_HPP

#include "met_maker/maker.hpp"

namespace correction {
class Correction;
class PluginRegistry;
}

namespace atlas::met::adapter {

inline constexpr char maker_correction_name[] = "build_met_terms";
void register_maker_nodes(correction::PluginRegistry& registry);
[[nodiscard]] maker::Result invoke_maker(const correction::Correction& evaluator,
                                         const maker::EventInput& event);

} // namespace atlas::met::adapter
#endif
