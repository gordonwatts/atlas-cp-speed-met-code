// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "soft.hpp"

#include <algorithm>
#include <cmath>

namespace atlas::met::soft {
namespace {

constexpr std::string_view soft_name = "PVSoftTrk";
constexpr std::uint64_t soft_mask = 0x10020;
constexpr std::uint64_t truth_nonint = 0x21000;
constexpr std::uint64_t invisible = 0x100000;

bool has(std::uint64_t source, std::uint64_t pattern) {
    return (source & pattern) == pattern;
}

} // namespace

MetPhi met_phi(Momentum term) {
    // MissingET_v1 computes the square and sum with float operands before
    // narrowing the magnitude and azimuth getter results to aux float.
    return {std::sqrt(term.mpx * term.mpx + term.mpy * term.mpy),
            static_cast<float>(std::atan2(static_cast<double>(term.mpy), static_cast<double>(term.mpx)))};
}

std::optional<Seed> seed_for(Momentum term) {
    const double scaled = 1.e6 * met_phi(term).phi;
    // The stock floating-to-int cast has no reproducible NaN behavior.
    if (!std::isfinite(scaled)) return std::nullopt;
    // A finite atan2 output fits int32. The unsigned-long setter receives
    // sign extension followed by modulo conversion; TRandom stores low 32 bits.
    const auto integer = static_cast<std::int32_t>(scaled);
    const auto setter_argument = static_cast<std::uint64_t>(integer);
    return Seed{integer, setter_argument, integer ? static_cast<std::uint32_t>(setter_argument) : 1u};
}

Recoil calc_recoil(const EventView& event) {
    if (!event.truth_available) return {RecoilState::missing_truth_container};
    const auto truth = std::find_if(event.truth.begin(), event.truth.end(), [](const Term& term) {
        return has(term.source, truth_nonint);
    });
    if (truth == event.truth.end()) return {RecoilState::missing_truth_nonint};

    Recoil recoil{RecoilState::regular, truth->momentum.mpx, truth->momentum.mpy, truth->momentum.sumet};
    for (const auto& term : event.terms) {
        if (has(term.source, soft_mask) || has(term.source, final_source)) continue;
        recoil.mpx -= term.momentum.mpx;
        recoil.mpy -= term.momentum.mpy;
        recoil.sumet += term.momentum.sumet;
    }
    if (recoil.sumet < 0) return {RecoilState::negative_sumet};
    // TEnv Units is absent in the initialized track soft config, hence m_units=1.
    if (recoil.sumet < 1.e-6) {
        recoil = {RecoilState::small_sumet_floor, 1.e-6, 0., 1.e-6};
    }
    recoil.pt = std::sqrt(recoil.mpx * recoil.mpx + recoil.mpy * recoil.mpy);
    return recoil;
}

FinalResult ordered_sum(std::span<const Term> terms) {
    const auto soft = std::find_if(terms.begin(), terms.end(), [](const Term& term) {
        return term.name == soft_name;
    });
    if (soft == terms.end()) return {Status::missing_soft, {}, {}};
    if (std::any_of(terms.begin(), terms.end(), [](const Term& term) {
        return term.name == final_name;
    })) return {Status::duplicate_final, {}, {}};

    Momentum sum{};
    for (const auto& term : terms) {
        if (has(term.source, final_source) || term.source == invisible) continue;
        if (soft->source && has(term.source, soft_mask) && !has(term.source, soft->source)) continue;
        // xAOD::MissingET_v1::operator+= updates three float aux stores.
        sum.mpx += term.momentum.mpx;
        sum.mpy += term.momentum.mpy;
        sum.sumet += term.momentum.sumet;
    }
    return {Status::success, sum, met_phi(sum)};
}

} // namespace atlas::met::soft
