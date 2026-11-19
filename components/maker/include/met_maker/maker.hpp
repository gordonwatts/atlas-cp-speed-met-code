// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_MAKER_MAKER_HPP
#define ATLAS_MET_MAKER_MAKER_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace atlas::met::maker {

using Mask = std::uint64_t;

struct OwnerKey {
    std::uint32_t value{};
    bool operator==(const OwnerKey&) const = default;
};

struct ObjectId {
    OwnerKey owner;
    std::uint32_t index{};
    bool operator==(const ObjectId&) const = default;
};

enum class LinkState { null, unresolved, valid };

struct ObjectLink {
    LinkState state{LinkState::null};
    ObjectId id;
};

enum class ParticleKind : std::uint8_t {
    other,
    jet,
    electron,
    photon,
    tau,
    muon
};

struct PtEtaPhiM {
    double pt{};
    double eta{};
    double phi{};
    double mass{};
};

struct Particle {
    ObjectId id;
    ParticleKind kind{};
    PtEtaPhiM p4;
    std::optional<ObjectLink> original_link;
    std::optional<std::uint16_t> muon_type;
    std::optional<float> energy_loss;
    std::optional<std::uint8_t> energy_loss_type;
    // Link resolution is adapter work. Null or unresolved ID tracks become nullopt.
    std::optional<float> id_track_pt;
};

struct ViewInfo {
    bool available{true};
    bool has_shallow_store{};
    bool owns_elements{};
    bool owning_key_available{true};
};

struct HardObjectView {
    std::span<const Particle> selected;
    ViewInfo view;
};

struct Jet {
    Particle particle;
    std::optional<PtEtaPhiM> constituent_p4;
    std::optional<float> nnjvt;
    // Supply the initialized selector's JetEtaName aux float, including for "eta".
    std::optional<float> jvt_eta;
    std::optional<std::span<const int>> num_trk_pt500;
    std::optional<std::span<const float>> sum_pt_trk_pt500;
    std::optional<float> em_frac;
    std::optional<float> width;
    std::optional<float> ps_frac;
    std::optional<std::span<const float>> energy_per_sampling;
};

struct ConstituentVector {
    float px{};
    float py{};
    float pz{};
    float energy{};
    float sum_pt{};
};

struct KeyedVector {
    Mask key{};
    ConstituentVector vector;
};

struct OverlapEdge {
    std::size_t slot{};
    std::uint8_t types{};
};

struct Association {
    ObjectLink jet;
    bool is_misc{};
    std::span<const ObjectLink> objects;
    // One row per object slot, including unresolved and null slots.
    std::span<const std::span<const OverlapEdge>> overlap_rows;
    std::span<const KeyedVector> cal;
    std::span<const KeyedVector> trk;
    ConstituentVector jet_tracks;
};

struct MetVector {
    float mpx{};
    float mpy{};
    float sumet{};
};

struct CoreTerm {
    std::string_view name;
    Mask source{};
    MetVector momentum;
};

struct EventInput {
    // All resolvable target/current/original rows required by the event. Include
    // association-linked rows outside the selected views and original-link targets.
    std::span<const Particle> input_particles;
    HardObjectView electrons;
    HardObjectView photons;
    HardObjectView taus;
    HardObjectView muons;
    std::span<const Jet> jets;
    ViewInfo jet_view;
    std::span<const Association> associations;
    std::span<const CoreTerm> core_terms;
    // Raw xAOD::VxType::VertexType per vertex in the original container order.
    std::span<const std::uint8_t> vertex_types;
    bool association_map_available{true};
    bool core_container_available{true};
};

struct NnjvtCutMap {
    std::vector<float> pt_edges_mev;
    std::vector<float> eta_edges;
    // Row-major [pt interval][eta interval].
    std::vector<float> cuts;
};

struct JvtRange {
    float min_pt_mev{20000.F};
    float max_pt_mev{60000.F};
    float min_abs_eta{-1.F};
    float max_abs_eta{2.5F};
};

class Config {
public:
    // Construct once per loaded payload. Run evaluate with an immutable Config.
    explicit Config(NnjvtCutMap cut_map, JvtRange jvt_range = {});
    [[nodiscard]] float jvt_cut(float pt, float eta) const;
    [[nodiscard]] const JvtRange& jvt_range() const noexcept { return range_; }

    static constexpr std::string_view recipe = "AnalysisBase-25.2.117-PFlow-Tight-PVSoftTrk";
    static constexpr std::string_view working_point = "FixedEffPt";
    static constexpr std::string_view core_key = "PVSoftTrkCore";
    static constexpr double central_jet_pt = 20000.;
    static constexpr double forward_jet_pt = 30000.;
    static constexpr double forward_jet_eta = 2.5;
    static constexpr double max_jet_eta = 4.5;
    static constexpr double min_unique_energy_fraction = 0.;
    static constexpr double min_weighted_jet_pt = 20000.;
    static constexpr float missing_object_warning_pt = 7000.F;
    static constexpr double min_electron_pt = 90000.;
    static constexpr double min_electron_residual_track_pt = 10000.;
    static constexpr double muon_fraction = 0.8;
    static constexpr double muon_jet_id_pt_ratio = 2.;
    static constexpr float muon_jet_track_count = 5.F;
    static constexpr float em_jet_track_count = 3.F;
    static constexpr double em_jet_fraction = 0.9;
    static constexpr double em_jet_width = 0.1;
    static constexpr double em_jet_ps_energy = 2500.;

private:
    NnjvtCutMap cuts_;
    JvtRange range_;
};

enum class TermType { electron, photon, tau, muon, muon_eloss, jet, soft_track };

struct Contributor {
    ObjectId current;
    // Stock does not write soft jet weights to the TST term.
    std::optional<float> weight;
};

struct Term {
    TermType type{};
    std::string_view name;
    Mask source{};
    MetVector momentum;
    std::vector<Contributor> contributors;
};

enum class Status { success, failure, source_access_exception, stock_undefined, unsupported_input };
enum class Issue {
    missing_core_container,
    missing_association_map,
    empty_association_map,
    missing_hard_container,
    missing_jet_container,
    missing_owner_key,
    shallow_missing_original,
    invalid_original_link,
    unbound_valid_target,
    missing_electron_association,
    missing_other_association,
    missing_jet_association,
    missing_soft_core,
    non_track_core,
    missing_muon_eloss_term,
    missing_presampler_input,
    missing_required_aux,
    short_vertex_moment,
    short_sampling_vector,
    malformed_overlap_rows,
    stock_signed_shift_out_of_range,
    stock_mask_shift_out_of_range
};

enum class Severity { debug, warning, error };

struct Diagnostic {
    Issue issue;
    Severity severity;
    std::optional<TermType> term;
    std::optional<ObjectId> object;
    // Empty for ordinary logging/StatusCode failure.
    std::string_view detail;
};

enum class JetDisposition { missing_association, hard, soft_track, outside_soft_region, muon_eloss };

struct JetDecision {
    ObjectId current;
    JetDisposition disposition;
    std::optional<bool> jvt_pass;
    bool passes_pt_cuts{};
    bool eta_rejected_before_jvt{};
    bool rejected{};
};

struct Result {
    Status status{Status::success};
    std::vector<Term> terms;
    // Aux-vector availability applies to the whole xAOD output container.
    bool constituent_links_available{};
    bool constituent_weights_available{};
    std::vector<Mask> flags;
    std::vector<JetDecision> jet_decisions;
    std::vector<Diagnostic> diagnostics;
};

[[nodiscard]] Result evaluate(const Config& config, const EventInput& event);

}  // namespace atlas::met::maker
#endif
