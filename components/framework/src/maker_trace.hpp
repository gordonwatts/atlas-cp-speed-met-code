// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_FRAMEWORK_MAKER_TRACE_HPP
#define ATLAS_MET_FRAMEWORK_MAKER_TRACE_HPP

#include "maker_xaod.hpp"

namespace atlas::met::framework {
bool maker_trace_enabled();
std::size_t next_maker_trace_entry();
void trace_maker(std::size_t entry, std::size_t state, const maker::EventInput& input,
                 const maker::Result& result, const MakerInputs& owning);
}
#endif
