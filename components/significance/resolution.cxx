/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 *
 * The source equivalents are JetSmearingCorrection, eg_resolution,
 * egammaEnergyCorrectionTool, MuonCalibTool/CalibContainer and TauResolutionTool
 * from AnalysisBase 25.2.117. The ROOT lookup arithmetic follows TAxis,
 * TH1::Interpolate and TGraph::Eval from ROOT 6.40.02.
 */
#include "significance.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace atlas::met_significance {
namespace {

constexpr double mev_per_gev = 1.e3;

double clamp_jet(double value, const std::array<double, 2>& limits) {
    if (value >= limits[1]) return limits[1] - 1.e-6;
    if (value <= limits[0]) return limits[0] + 1.e-6;
    return value;
}

double clamp_muon(double value, const Axis& axis) {
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    return std::max(std::min(value, axis.edges.back() - epsilon), axis.edges.front() + epsilon);
}

double add_nominal_smear(double expected, double nominal) {
    return std::sqrt(std::pow(expected, 2) + std::pow(nominal, 2));
}

double id_resolution(const MuonData& data, double pt, double eta, double phi) {
    if (pt == 0) return 1.e12;
    const auto& expected = data.id_expected;
    const auto& nominal = data.id_nominal;
    const double expected_p1 = expected[0].muon_lookup(eta, phi);
    double expected_p2 = expected[1].muon_lookup(eta, phi);
    const double p2_tan = expected[2].muon_lookup(eta, phi);
    if (p2_tan) expected_p2 = p2_tan;
    const double r1 = nominal[0].muon_lookup(eta, phi);
    const double r2 = nominal[1].muon_lookup(eta, phi);
    const double p1 = add_nominal_smear(expected_p1, r1);
    double p2 = add_nominal_smear(expected_p2, r2);
    if (p2_tan) p2 = p2 * std::sinh(eta) * std::sinh(eta);
    return std::sqrt(std::pow(p1, 2) + std::pow(p2 * pt, 2));
}

double ms_resolution(const MuonData& data, double pt, double eta, double phi) {
    if (pt == 0) return 1.e12;
    const auto& expected = data.ms_expected;
    const auto& nominal = data.ms_nominal;
    const double expected_p0 = expected[0].muon_lookup(eta, phi);
    const double expected_p1 = expected[1].muon_lookup(eta, phi);
    const double expected_p2 = expected[2].muon_lookup(eta, phi);
    const double r0 = nominal[0].muon_lookup(eta, phi);
    const double r1 = nominal[1].muon_lookup(eta, phi);
    const double r2 = nominal[2].muon_lookup(eta, phi);
    const double p0 = add_nominal_smear(expected_p0, r0);
    const double p1 = add_nominal_smear(expected_p1, r1);
    const double p2 = add_nominal_smear(expected_p2, r2);
    return std::sqrt(std::pow(p0 / pt, 2) + std::pow(p1, 2) + std::pow(p2 * pt, 2));
}

std::size_t tau_eta_index(float eta) {
    const float abs_eta = std::abs(eta);
    constexpr std::array<double, 5> limits{0.3, 0.8, 1.3, 1.6, 2.6};
    for (std::size_t bin = 0; bin < limits.size(); ++bin) {
        if (abs_eta < limits[bin]) return bin;
    }
    return 99;
}

} // namespace

std::size_t Axis::root_bin(double value) const {
    const auto count = bins();
    const double lower = edges.front();
    const double upper = edges.back();
    if (value < lower) return 0;
    if (!(value < upper)) return count + 1; // ROOT assigns NaN to overflow.
    if (!uniform) {
        return std::upper_bound(edges.begin(), edges.end(), value) - edges.begin();
    }
    const double width = (upper - lower) / static_cast<double>(count);
    const auto approximate = static_cast<int>((value - lower) / width);
    return 1 + approximate - (value < lower + width * approximate)
        + (lower + width * (approximate + 1) <= value);
}

double Histogram1::lookup(double value) const {
    return values.at(axis.root_bin(value));
}

double Histogram2::muon_lookup(double eta, double phi) const {
    const auto eta_bin = x.root_bin(clamp_muon(eta, x));
    const auto phi_bin = y.root_bin(clamp_muon(phi, y));
    return values.at(eta_bin + (x.bins() + 2) * phi_bin);
}

double JetProjection::interpolate(double pt_gev) const {
    // ROOT's overflow-neighbor arithmetic propagates a NaN abscissa.
    if (std::isnan(pt_gev)) return pt_gev;
    const auto last = histogram.axis.bins();
    if (pt_gev <= centers[1]) return histogram.values[1];
    if (pt_gev >= centers[last]) return histogram.values[last];

    const std::size_t bin = histogram.axis.root_bin(pt_gev);
    const std::size_t left = pt_gev <= centers[bin] ? bin - 1 : bin;
    const double x0 = centers[left];
    const double x1 = centers[left + 1];
    const double y0 = histogram.values[left];
    const double y1 = histogram.values[left + 1];
    return y0 + (pt_gev - x0) * ((y1 - y0) / (x1 - x0));
}

double JetTable::lookup(double pt_gev, double eta) const {
    const double x = clamp_jet(pt_gev, pt_limits);
    const double y = clamp_jet(std::abs(eta), eta_limits);
    // A NaN eta selects the overflow y-bin. The stock projection vector has
    // no overflow entry and its checked lookup throws.
    return projections.at(abs_eta.root_bin(y)).interpolate(x);
}

JetResolution jet_resolution(const JetData& payload, double pt_mev, double eta) {
    const double data = payload.data.lookup(pt_mev / mev_per_gev, eta);
    const double mc = payload.mc.lookup(pt_mev / mev_per_gev, eta);
    return {data, mc, static_cast<float>(std::max(data, mc))};
}

EgammaKind photon_kind(std::size_t vertex_link_count, bool forward_electron_author) {
    if (forward_electron_author) return EgammaKind::electron;
    return vertex_link_count ? EgammaKind::converted_photon : EgammaKind::unconverted_photon;
}

float calo_frame_eta(bool forward_electron, double cluster_eta,
                     std::optional<double> eta_moment, std::optional<float> eta_aux) {
    if (forward_electron) return static_cast<float>(cluster_eta);
    if (eta_moment) return static_cast<float>(*eta_moment);
    if (eta_aux) return *eta_aux;
    throw std::invalid_argument("MET significance etaCalo missing: release fallback is undefined");
}

double egamma_resolution(const EgammaData& data, EgammaKind kind, double energy_mev,
                         double cluster_eta, std::optional<float> eta_calo, bool with_constant_term) {
    const auto& coefficients = data.sigma90.at(static_cast<std::size_t>(kind));
    const float abs_eta = std::abs(cluster_eta);
    std::size_t bin = data.eta_edges.size() - 2;
    for (std::size_t edge = 1; edge < data.eta_edges.size(); ++edge) {
        if (abs_eta < data.eta_edges[edge]) {
            bin = edge - 1;
            break;
        }
    }

    const double energy_gev = energy_mev * 1.e-3;
    const double sampling = coefficients[0][bin];
    const double noise = coefficients[1][bin];
    const double constant = coefficients[2][bin];
    const double base_variance = sampling * sampling / energy_gev
        + noise * noise / energy_gev / energy_gev + constant * constant;
    double variance = std::pow(std::sqrt(base_variance), 2);

    const double et = energy_mev / std::cosh(cluster_eta);
    const double limited_et = std::clamp(et, 5000., 50000.);
    const double pileup_noise = std::sqrt(54.) * (60. + 40. * std::log(limited_et / 10000.) / std::log(5.));
    variance += std::pow(pileup_noise / et, 2);

    if (with_constant_term) {
        if (!eta_calo) throw std::invalid_argument("MET significance photon requires etaCalo");
        variance += std::pow(std::max(0., data.constant_term.lookup(*eta_calo)), 2);
    }
    return std::sqrt(variance);
}

MuonDetector muon_detector(std::uint16_t muon_type) {
    if (muon_type == 0) return MuonDetector::cb;
    if (muon_type == 1) return MuonDetector::ms;
    return MuonDetector::id;
}

double track_pt_gev(std::optional<TrackParameters> track) {
    if (!track) return 0;
    // xAOD TrackingDetails::genvecP4 projects in float before PxPyPzM4D::Pt.
    float p = 10.e6;
    if (std::fabs(track->q_over_p) > 0.) p = 1 / std::fabs(track->q_over_p);
    const float sin_theta = std::sin(track->theta);
    const double px = p * sin_theta * std::cos(track->phi);
    const double py = p * sin_theta * std::sin(track->phi);
    return std::sqrt(px * px + py * py) * 1.e-3;
}

double muon_resolution(const MuonData& data, MuonDetector detector, std::uint32_t run_number,
                       bool is_simulation, double eta, double phi, const MuonTrackPt& pt) {
    if (!is_simulation || (run_number != 410000 && run_number != 330000)) {
        throw std::invalid_argument("MET significance muon payload is for full-simulation Data22 periods");
    }
    switch (detector) {
    case MuonDetector::ms:
        return ms_resolution(data, pt.ms_gev, eta, phi);
    case MuonDetector::id:
        return id_resolution(data, pt.id_gev, eta, phi);
    case MuonDetector::cb: {
        const double sigma_id = id_resolution(data, pt.id_gev, eta, phi) * pt.id_gev;
        const double sigma_ms = ms_resolution(data, pt.ms_gev, eta, phi) * pt.ms_gev;
        const double denominator = pt.cb_gev * std::sqrt(sigma_id * sigma_id + sigma_ms * sigma_ms);
        return denominator ? sigma_id * sigma_ms / denominator : 0.;
    }
    }
    throw std::invalid_argument("MET significance muon detector must be MS, ID or CB");
}

double TauGraph::evaluate(double pt_mev) const {
    const double value = std::min(pt_mev, max_x);
    if (x.size() == 1) return y.front();
    std::size_t low;
    std::size_t up;
    if (sorted_bit) {
        const auto found = std::lower_bound(x.begin(), x.end(), value);
        if (found != x.end() && *found == value) return y[found - x.begin()];
        const auto index = static_cast<std::size_t>(found - x.begin());
        low = index == 0 ? 0 : std::min(index - 1, x.size() - 2);
        up = low + 1;
    } else {
        const auto missing = x.size();
        low = up = missing;
        std::size_t low_neighbor = missing;
        std::size_t up_neighbor = missing;
        for (std::size_t i = 0; i < x.size(); ++i) {
            if (x[i] < value) {
                if (low == missing || x[i] > x[low]) {
                    low_neighbor = low;
                    low = i;
                } else if (low_neighbor == missing) low_neighbor = i;
            } else if (x[i] > value) {
                if (up == missing || x[i] < x[up]) {
                    up_neighbor = up;
                    up = i;
                } else if (up_neighbor == missing) up_neighbor = i;
            } else return y[i];
        }
        if (up == missing) {
            up = low;
            low = low_neighbor;
        }
        if (low == missing) {
            low = up;
            up = up_neighbor;
        }
    }
    const double x_low = x.at(low);
    const double x_up = x.at(up);
    if (x_low == x_up) return y[low];
    return y[up] + (value - x_up) * (y[low] - y[up]) / (x_low - x_up);
}

double tau_resolution(const TauData& data, std::optional<int> decay_mode, float tes_eta,
                      float final_eta, double final_pt_mev) {
    if (!decay_mode || *decay_mode < 0 || *decay_mode > 4 || tau_eta_index(tes_eta) > 4) {
        return 1.0;
    }
    return data.graphs.at(*decay_mode).at(tau_eta_index(final_eta)).evaluate(final_pt_mev);
}

} // namespace atlas::met_significance
