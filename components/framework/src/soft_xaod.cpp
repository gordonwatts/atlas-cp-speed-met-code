// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "soft_xaod.hpp"
#include "bridge.hpp"
#include "maker_xaod.hpp"
#include "soft.hpp"

#include <AsgDataHandles/VarHandleKey.h>
#include <AsgDataHandles/ReadHandle.h>
#include <utility>
#include <AthContainers/Decorator.h>
#include <xAODEventInfo/EventInfo.h>

namespace atlas::met::framework {
namespace {
const SG::ConstAccessor<std::vector<ParticleLink>> links{"ConstitObjectLinks"};
const SG::Decorator<float> met_decoration{"met"};
const SG::Decorator<float> phi_decoration{"phi"};

const soft::SoftFile& payload()
{
    static const soft::SoftFile file(absolute_file("ATLAS_MET_SOFT_PAYLOAD"));
    return file;
}

soft::Momentum momentum(const xAOD::MissingET& term)
{
    return {term.mpx(), term.mpy(), term.sumet()};
}

std::vector<soft::Term> rows(const xAOD::MissingETContainer& container)
{
    std::vector<soft::Term> output;
    output.reserve(container.size());
    for (const auto* term : container) {
        if (!term) throw std::runtime_error("Null xAOD MET term");
        std::optional<std::size_t> count;
        if (links.isAvailable(*term)) count = links(*term).size();
        output.push_back({term->name(), term->source(), momentum(*term), count});
    }
    return output;
}

void write(xAOD::MissingET& term, soft::Momentum value)
{
    term.setMpx(value.mpx);
    term.setMpy(value.mpy);
    term.setSumet(value.sumet);
}
} // namespace

StatusCode apply_soft(xAOD::MissingETContainer& met, std::string_view state, const EventContext&)
{
    SG::ReadHandle<xAOD::EventInfo> event("EventInfo");
    SG::ReadHandle<xAOD::MissingETContainer> truth("MET_Truth");
    const bool have_event = event.isValid();
    const bool have_truth = truth.isValid();
    const auto terms = rows(met);
    const auto truth_terms = have_truth ? rows(*truth.cptr()) : std::vector<soft::Term>{};
    const soft::EventView input{terms, truth_terms,
        have_event ? event->actualInteractionsPerCrossing() : 0.F, have_truth, have_event,
        have_event && event->eventType(xAOD::EventInfo::IS_SIMULATION), true};
    const auto result = payload().apply(input, state.empty() ? "NOSYS" : state);
    if (result.status != soft::Status::success || !result.term_index) {
        throw std::runtime_error("MET soft kernel status " + std::to_string(static_cast<int>(result.status)));
    }
    write(*met.at(*result.term_index), result.momentum);
    return StatusCode::SUCCESS;
}

StatusCode build_final(xAOD::MissingETContainer& met)
{
    const auto raw = rows(met);
    const auto result = payload().build_final(raw);
    if (result.status != soft::Status::success) {
        throw std::runtime_error("MET final kernel status " + std::to_string(static_cast<int>(result.status)));
    }
    auto* total = new xAOD::MissingET;
    met.push_back(total);
    total->setName(std::string(soft::final_name));
    total->setSource(soft::final_source);
    write(*total, result.momentum);
    for (const auto* term : std::as_const(met)) {
        const auto decorations = term == total ? result.decorations : payload().decorate(momentum(*term));
        met_decoration(*term) = decorations.met;
        phi_decoration(*term) = decorations.phi;
    }
    return StatusCode::SUCCESS;
}
} // namespace atlas::met::framework
