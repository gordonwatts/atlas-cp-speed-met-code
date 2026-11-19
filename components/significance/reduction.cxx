/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 *
 * Fixed-configuration numerical path of METUtilities::METSignificance in
 * AnalysisBase 25.2.117. TMath/TVector2/3 angle cases follow ROOT 6.40.02.
 */
#include "significance.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace atlas::met_significance {
namespace {

constexpr double mev_per_gev = 1.e3;
constexpr double soft_sigma_gev = 8.5;
constexpr std::uint64_t total_mask = 0x10040;
constexpr std::uint64_t soft_mask = 0x10020;
constexpr std::uint64_t track_mask = 0x04000;
constexpr std::uint64_t invisible_source = 0x100000;

bool has_pattern(std::uint64_t source, std::uint64_t mask) {
    return (source & mask) == mask;
}

double phi(double px, double py) {
    if (px != 0) return std::atan2(py, px);
    if (py == 0) return 0;
    return py > 0 ? std::numbers::pi / 2 : -std::numbers::pi / 2;
}

double delta_phi(double left, double right) {
    double delta = left - right;
    // Both arguments come from atan2 or an axis special case above.
    if (delta >= std::numbers::pi) delta -= 2 * std::numbers::pi;
    if (delta < -std::numbers::pi) delta += 2 * std::numbers::pi;
    return delta;
}

double met_phi(const Momentum& momentum) {
    return phi(momentum.mpx_mev / mev_per_gev, momentum.mpy_mev / mev_per_gev);
}

Covariance rotate_diagonal(double var_x, double var_y, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double cc = c * c;
    const double ss = s * s;
    const double cs = c * s;
    return {
        var_x * cc + var_y * ss - cs * (0.0 + 0.0),
        var_x * ss + var_y * cc + cs * (0.0 + 0.0),
        0.0 * cc - 0.0 * ss + cs * (var_x - var_y)
    };
}

} // namespace

Covariance& Covariance::operator+=(const Covariance& rhs) {
    longitudinal += rhs.longitudinal;
    transverse += rhs.transverse;
    off_diagonal += rhs.off_diagonal;
    return *this;
}

TermAction term_action(std::string_view term, std::uint64_t source) {
    if (has_pattern(source, total_mask) || source == invisible_source) return TermAction::skip;
    if (has_pattern(source, soft_mask)) {
        return term == "PVSoftTrk" && has_pattern(source, track_mask)
            ? TermAction::track_soft : TermAction::skip;
    }
    return TermAction::constituents;
}

bool uses_constituent(std::string_view term, ObjectKind kind) {
    return kind != ObjectKind::jet || term == "RefJet";
}

Covariance object_covariance(const Momentum& final, const Constituent& object) {
    const float square = object.relative_pt * object.relative_pt;
    const double var = square * object.pt_mev * object.pt_mev / mev_per_gev / mev_per_gev;
    return rotate_diagonal(var, 0, delta_phi(met_phi(final), phi(object.p4_px_mev, object.p4_py_mev)));
}

Covariance soft_covariance(const Momentum& final, const Momentum& soft) {
    const double var = soft_sigma_gev * soft_sigma_gev;
    return rotate_diagonal(var, var, delta_phi(met_phi(final), met_phi(soft)));
}

Decorations finalize(const Covariance& cov, double met_gev, double sumet_gev,
                     double soft_sumet_gev, State& state) {
    double ht_gev = sumet_gev;
    if (cov.longitudinal != 0) {
        state.rho = cov.off_diagonal / std::sqrt(cov.longitudinal * cov.transverse);
        const double denominator = std::abs(state.rho) >= 0.9
            ? cov.longitudinal
            : cov.longitudinal * (1 - std::pow(state.rho, 2));
        state.squared_significance = std::pow(met_gev, 2) / denominator;
        ht_gev -= soft_sumet_gev;
    }
    return {
        static_cast<float>(state.squared_significance > 0 ? std::sqrt(state.squared_significance) : -1.0),
        static_cast<float>(cov.longitudinal > 0 ? met_gev / std::sqrt(cov.longitudinal) : -1.0),
        static_cast<float>(sumet_gev > 0 ? met_gev / std::sqrt(sumet_gev) : -1.0),
        static_cast<float>(ht_gev > 0 ? met_gev / std::sqrt(ht_gev) : -1.0)
    };
}

Result reduce(std::span<const Term> terms, State& state) {
    Result result{};
    result.status = Status::success;
    const auto total = std::find_if(terms.begin(), terms.end(), [](const Term& term) {
        return term.name == "Final";
    });
    if (total == terms.end()) {
        result.final_state = FinalState::missing;
        result.status = Status::missing_final;
        result.decorations = finalize({}, -1, -1, 0, state);
        return result;
    }
    if (!has_pattern(total->source, total_mask)) {
        result.final_state = FinalState::wrong_source;
        result.decorations = finalize({}, -1, -1, 0, state);
        return result;
    }
    result.final_state = FinalState::ready;

    const auto& final = total->momentum;
    double soft_sumet_gev = 0;
    bool included_soft = false;
    for (const auto& term : terms) {
        if (included_soft && !has_pattern(term.source, total_mask)
            && has_pattern(term.source, soft_mask | track_mask)) {
            ++result.later_track_soft_terms;
        }
        switch (term_action(term.name, term.source)) {
        case TermAction::skip:
            continue;
        case TermAction::track_soft: {
            if (included_soft) continue;
            const auto contribution = soft_covariance(final, term.momentum);
            result.type_covariance[2] += contribution;
            result.covariance += contribution;
            soft_sumet_gev = term.momentum.sumet_mev / mev_per_gev;
            included_soft = true;
            continue;
        }
        case TermAction::constituents:
            break;
        }

        if (!term.links_available) {
            result.status = Status::unavailable_constituents;
            return result;
        }
        for (const auto& link : term.links) {
            if (!link.valid) {
                result.status = Status::invalid_constituent;
                return result;
            }
            switch (link.kind) {
            case ObjectKind::jet:
            case ObjectKind::electron:
            case ObjectKind::muon:
            case ObjectKind::photon:
            case ObjectKind::tau:
                break;
            default:
                result.status = Status::unsupported_object;
                return result;
            }
            if (!uses_constituent(term.name, link.kind)) continue;
            const auto contribution = object_covariance(final, link);
            result.type_covariance[static_cast<std::size_t>(link.kind)] += contribution;
            result.covariance += contribution;
        }
    }
    // MissingET_v1::met returns float after squaring and summing float inputs.
    const float squared = final.mpx_mev * final.mpx_mev + final.mpy_mev * final.mpy_mev;
    const float met_mev = std::sqrt(squared);
    result.longitudinal_nonzero = result.covariance.longitudinal != 0;
    result.decorations = finalize(result.covariance, met_mev / mev_per_gev,
                                  final.sumet_mev / mev_per_gev, soft_sumet_gev, state);
    return result;
}

} // namespace atlas::met_significance
