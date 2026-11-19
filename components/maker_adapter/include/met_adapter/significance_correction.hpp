// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_SIGNIFICANCE_CORRECTION_HPP
#define ATLAS_MET_SIGNIFICANCE_CORRECTION_HPP

#include "significance.hpp"

namespace correction {
class Correction;
class PluginRegistry;
}

namespace atlas::met::adapter {

namespace significance = atlas::met_significance;
inline constexpr char significance_correction_name[] = "met_significance_event";

struct EgammaInput {
    double energy_mev{};
    std::optional<double> cluster_eta;
    bool forward_author{};
    std::optional<double> eta_calo_moment;
    std::optional<float> eta_calo_aux;
    std::size_t vertex_link_count{};
};

struct MuonInput {
    std::uint16_t type{};
    std::optional<significance::TrackParameters> combined;
    std::optional<significance::TrackParameters> inner_detector;
    std::optional<significance::TrackParameters> extrapolated_ms;
};

struct TauInput {
    std::optional<int> decay_mode;
    std::optional<float> tes_eta;
    std::optional<float> final_eta;
    std::optional<float> final_pt_mev;
};

struct SignificanceLink {
    significance::ObjectKind kind{};
    bool valid{};
    double pt_mev{};
    double eta{};
    double phi{};
    double p4_px_mev{};
    double p4_py_mev{};
    std::optional<EgammaInput> egamma;
    std::optional<MuonInput> muon;
    std::optional<TauInput> tau;
};

struct SignificanceTerm {
    std::string_view name;
    std::uint64_t source{};
    significance::Momentum momentum;
    std::span<const SignificanceLink> links;
    bool links_available{};
};

struct SignificanceEvent {
    std::span<const SignificanceTerm> terms;
    significance::State previous;
    std::uint32_t run_number{};
    bool simulation{};
    float average_mu{};
};

struct SignificanceResult {
    significance::Result reduction;
    significance::State next;
};

void register_significance_event_node(correction::PluginRegistry& registry);
[[nodiscard]] SignificanceResult invoke_significance(const correction::Correction& evaluator,
                                                      const SignificanceEvent& event);

} // namespace atlas::met::adapter
#endif
