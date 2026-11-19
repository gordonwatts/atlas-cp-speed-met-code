// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_SOFT_XAOD_HPP
#define ATLAS_MET_SOFT_XAOD_HPP
#include <AsgMessaging/StatusCode.h>
#include <xAODMissingET/MissingETContainer.h>
#include <string_view>
class EventContext;

namespace atlas::met::framework {
StatusCode apply_soft(xAOD::MissingETContainer&, std::string_view state, const EventContext&);
StatusCode build_final(xAOD::MissingETContainer&);
}
#endif
