/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
#ifndef ATLAS_MET_SIGNIFICANCE_HPP
#define ATLAS_MET_SIGNIFICANCE_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace atlas::met_significance {

struct Axis {
    bool uniform;
    std::vector<double> edges;
    std::size_t bins() const { return edges.size() - 1; }
    std::size_t root_bin(double value) const;
};

struct Histogram1 {
    Axis axis;
    std::vector<double> values; // includes underflow and overflow
    double lookup(double value) const;
};

struct Histogram2 {
    Axis x;
    Axis y;
    std::vector<double> values; // x runs fastest, including flow on each axis
    double muon_lookup(double eta, double phi) const;
};

struct JetProjection {
    Histogram1 histogram;
    std::vector<double> centers; // all centers, including flow-bin centers
    double interpolate(double pt_gev) const;
};

struct JetTable {
    Axis abs_eta;
    std::array<double, 2> pt_limits; // source TH2 bin-low-edge values, before ProjectionX
    std::array<double, 2> eta_limits;
    std::vector<JetProjection> projections; // same indexing as ROOT y bins, excluding overflow
    double lookup(double pt_gev, double eta) const;
};

struct JetData {
    JetTable data;
    JetTable mc;
};

struct JetResolution {
    double data;
    double mc;
    float relative_pt;
};

JetResolution jet_resolution(const JetData&, double pt_mev, double eta);

enum class EgammaKind { electron, unconverted_photon, converted_photon };

struct EgammaData {
    std::vector<double> eta_edges;
    // By EgammaKind, then sampling/noise/constant. These have regular bins only.
    std::array<std::array<std::vector<double>, 3>, 3> sigma90;
    Histogram1 constant_term;
};

EgammaKind photon_kind(std::size_t vertex_link_count, bool forward_electron_author = false);

// Mirrors get_eta_calo precedence. The missing-calorimeter-frame branch in the
// release returns an uninitialized local. This path has no reproducible value.
float calo_frame_eta(bool forward_electron, double cluster_eta,
                     std::optional<double> eta_moment, std::optional<float> eta_aux);

double egamma_resolution(const EgammaData&, EgammaKind, double energy_mev,
                         double cluster_eta, std::optional<float> eta_calo, bool with_constant_term);

enum class MuonDetector : int { ms = 1, id = 2, cb = 3 };

struct MuonData {
    // r1_ID_MC, r2_ID_MC, r2_ID_TAN2_MC; r0_MS_MC, r1_MS_MC, r2_MS_MC.
    std::array<Histogram2, 3> id_expected;
    std::array<Histogram2, 3> ms_expected;
    // Nominal smear coefficients in ID r1/r2 and ME r0/r1/r2 order.
    std::array<Histogram2, 2> id_nominal;
    std::array<Histogram2, 3> ms_nominal;
};

struct TrackParameters {
    float q_over_p;
    float theta;
    float phi;
};

struct MuonTrackPt {
    double cb_gev;
    double id_gev;
    double ms_gev;
};

MuonDetector muon_detector(std::uint16_t muon_type);
double track_pt_gev(std::optional<TrackParameters>);
double muon_resolution(const MuonData&, MuonDetector, std::uint32_t run_number,
                       bool is_simulation, double eta, double phi, const MuonTrackPt&);

struct TauGraph {
    bool sorted_bit;
    double max_x;
    std::vector<double> x;
    std::vector<double> y;
    double evaluate(double pt_mev) const;
};

struct TauData {
    std::array<std::array<TauGraph, 5>, 5> graphs; // mode, eta bin
};

double tau_resolution(const TauData&, std::optional<int> decay_mode, float tes_eta,
                      float final_eta, double final_pt_mev);

// Each parser eagerly checks one exported section and its version. The shared
// provider can load just these four files without CorrectionLib or ROOT I/O.
JetData read_jet_data(std::string_view);
EgammaData read_egamma_data(std::string_view);
MuonData read_muon_data(std::string_view);
TauData read_tau_data(std::string_view);

struct Momentum {
    float mpx_mev;
    float mpy_mev;
    float sumet_mev;
};

struct Covariance {
    double longitudinal = 0;
    double transverse = 0;
    double off_diagonal = 0;
    Covariance& operator+=(const Covariance&);
};

enum class ObjectKind { jet = 1, electron = 3, muon = 4, photon = 5, tau = 6 };
enum class TermAction { skip, track_soft, constituents };

// Fixed stock keys Final, RefJet, PVSoftTrk and selection mask predicates.
TermAction term_action(std::string_view term_name, std::uint64_t source);
bool uses_constituent(std::string_view term_name, ObjectKind);

struct Constituent {
    ObjectKind kind;
    double pt_mev;
    double p4_px_mev;
    double p4_py_mev;
    float relative_pt;
    bool valid = true;
};

struct Term {
    std::string_view name;
    std::uint64_t source;
    Momentum momentum;
    std::span<const Constituent> links = {};
    bool links_available = true;
};

enum class Status {
    success,
    missing_final,
    invalid_constituent,
    unavailable_constituents,
    unsupported_object
};

enum class FinalState { ready, missing, wrong_source };

struct State {
    // Caller-owned for one tool-equivalent call sequence. Stock leaves these
    // untouched for a wrong-source/missing total or VarL == 0.
    double squared_significance = 0;
    double rho = 0;
};

struct Decorations {
    float significance;
    float directional;
    float met_over_sqrt_sumet;
    float met_over_sqrt_ht;
};

struct Result {
    Status status;
    FinalState final_state;
    Decorations decorations; // valid to publish only on algorithm success
    Covariance covariance;
    std::array<Covariance, 7> type_covariance; // stock codes; 2 is the soft term
    std::size_t later_track_soft_terms;
    bool longitudinal_nonzero;
};

Covariance object_covariance(const Momentum& final, const Constituent&);
Covariance soft_covariance(const Momentum& final, const Momentum& soft);
Decorations finalize(const Covariance&, double met_gev, double sumet_gev,
                     double soft_sumet_gev, State&);

// Input views follow the release term/link order. A graph or link failure stops
// reduction. The adapter owns link/traversal errors and decoration publication.
Result reduce(std::span<const Term>, State&);

} // namespace atlas::met_significance
#endif
