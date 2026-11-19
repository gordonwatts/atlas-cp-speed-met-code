// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "preparation.hpp"
#include "significance.hpp"

#include <correction_plugin.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace atlas::met::preparation {
namespace {
using correction::PluginConfig;
using correction::PluginContext;
using correction::PluginNode;
using Values = std::vector<correction::Variable::Type>;
using Type = correction::Variable::VarType;
using Int = std::int64_t;

constexpr Int pass = 0xffffffff;
constexpr Int fail_one = pass ^ 1;
constexpr Int fail_two = pass ^ 3;
constexpr Int fail_three = pass ^ 7;
constexpr std::size_t outputs = 21;

constexpr std::array real_names{
    "pt", "eta", "phi", "mass", "cluster_eta", "cluster_phi", "eta_calo", "phi_calo", "pv_z",
    "emb1_energy", "emb2_energy", "emb3_energy", "eme1_energy", "eme2_energy", "eme3_energy",
    "emb1_eta", "emb2_eta", "eme1_eta", "eme2_eta",
    "iso20", "iso30", "iso40", "oldleak20", "oldleak30", "oldleak40",
    "eadded_l2", "eadded_l3", "track_iso", "track_d0", "track_phi", "track_z0", "track_vz",
    "track_var_d0", "beam_sigma_x", "beam_sigma_y", "beam_sigma_xy", "conv_x", "conv_y",
    "conv_pt1", "conv_pt2", "conv_trk1_q_over_p", "conv_trk1_theta", "conv_trk1_phi",
    "conv_trk2_q_over_p", "conv_trk2_theta", "conv_trk2_phi", "leak_et", "rz_x", "rz_y"
};
enum Real : std::size_t {
    pt, eta, phi, mass, cluster_eta, cluster_phi, eta_calo, phi_calo, pv_z,
    emb1_energy, emb2_energy, emb3_energy, eme1_energy, eme2_energy, eme3_energy,
    emb1_eta, emb2_eta, eme1_eta, eme2_eta,
    iso20, iso30, iso40, oldleak20, oldleak30, oldleak40,
    eadded_l2, eadded_l3, track_iso, track_d0, track_phi, track_z0, track_vz,
    track_var_d0, beam_sigma_x, beam_sigma_y, beam_sigma_xy, conv_x, conv_y,
    conv_pt1, conv_pt2, conv_trk1_q_over_p, conv_trk1_theta, conv_trk1_phi,
    conv_trk2_q_over_p, conv_trk2_theta, conv_trk2_phi, leak_et, rz_x, rz_y
};
constexpr std::array int_names{
    "has_cluster", "has_calo_moments", "sampling_pattern", "has_pv", "run_number", "random_run_number",
    "oq", "author", "id_flag", "clean_flag", "iso_available_mask", "oldleak_available_mask",
    "timing_available_mask", "has_track", "track_cov_diag_size", "has_conv_trk1", "has_conv_trk2",
    "conv_trk1_pixel", "conv_trk1_sct", "conv_trk2_pixel", "conv_trk2_sct", "has_conv_pts", "leak_key"
};
enum Integer : std::size_t {
    has_cluster, has_calo_moments, sampling_pattern, has_pv, run_number, random_run_number, oq, author,
    id_flag, clean_flag, iso_available_mask, oldleak_available_mask, timing_available_mask, has_track,
    track_cov_diag_size, has_conv_trk1, has_conv_trk2, conv_trk1_pixel, conv_trk1_sct,
    conv_trk2_pixel, conv_trk2_sct, has_conv_pts, leak_key
};
static_assert(real_names.size() == rz_y + 1 && int_names.size() == leak_key + 1);

template<class T>
T read(std::istream& in) {
    T value{};
    if (!(in >> value)) throw std::invalid_argument("Malformed MET egamma calibration table");
    return value;
}

std::vector<double> numbers(std::istream& in, std::size_t size) {
    if (size < 1 || size > 500000) throw std::invalid_argument("Invalid MET egamma calibration extent");
    std::vector<double> result(size);
    for (auto& v : result) {
        v = read<double>(in);
        if (!std::isfinite(v)) throw std::invalid_argument("Nonfinite MET egamma calibration value");
    }
    return result;
}

struct InterpolationAxis {
    atlas::met_significance::Axis axis;
    std::vector<double> centers, upper, width;
    int bins() const { return static_cast<int>(width.size()); }
    double center(int bin) const { return centers.at(bin); }

    static InterpolationAxis load(std::istream& in) {
        InterpolationAxis result;
        const int uniform = read<int>(in);
        const auto count = read<std::size_t>(in);
        if ((uniform != 0 && uniform != 1) || count < 2 || count > 500000)
            throw std::invalid_argument("Invalid shower-depth axis");
        result.axis.uniform = uniform;
        result.axis.edges = numbers(in, count);
        if (std::adjacent_find(result.axis.edges.begin(), result.axis.edges.end(),
                              std::greater_equal<double>{}) != result.axis.edges.end())
            throw std::invalid_argument("Unordered shower-depth edges");
        result.centers = numbers(in, count + 1);
        result.upper = numbers(in, count - 1);
        result.width = numbers(in, count - 1);
        return result;
    }

    int limit_to_range(float& value) const {
        constexpr float epsilon = 1.e-6;
        const int bin = axis.root_bin(value);
        if (bin < 1) value = axis.edges.front() + epsilon;
        if (bin > bins()) value = axis.edges.back() - epsilon;
        return axis.root_bin(value);
    }

    int lower_center(float value, int bin) const {
        const auto i = static_cast<std::size_t>(bin - 1);
        return upper[i] - value <= width[i] / 2 ? bin : bin - 1;
    }
};

struct ShowerDepth {
    InterpolationAxis x;
    InterpolationAxis y;
    static ShowerDepth load(const std::string& encoded) {
        std::istringstream in(encoded);
        in.imbue(std::locale::classic());
        if (read<std::string>(in) != "MetCaloDepthV1") throw std::invalid_argument("Expected MetCaloDepthV1");
        auto x = InterpolationAxis::load(in);
        auto y = InterpolationAxis::load(in);
        in >> std::ws;
        if (!in.eof()) throw std::invalid_argument("Surplus shower-depth axis fields");
        return {std::move(x), std::move(y)};
    }
};

struct Row {
    const Values& values;
    const std::array<std::size_t, real_names.size()>& real;
    const std::array<std::size_t, int_names.size()>& integer;
    float f(std::size_t slot) const { return std::get<double>(values[real[slot]]); }
    Int i(std::size_t slot) const { return std::get<Int>(values[integer[slot]]); }
    bool flag(std::size_t slot) const {
        const auto value = i(slot);
        if (value != 0 && value != 1) throw std::invalid_argument("Expected 0/1 egamma availability");
        return value == 1;
    }
};

struct Layers {
    float eta1;
    float eta2;
    float energy;
};

bool sampling(Int pattern, bool endcap, unsigned layer) {
    return ((pattern >> (layer + (endcap ? 4 : 0))) & 1u) != 0;
}

float layer_eta(const Row& row, unsigned layer) {
    const auto pattern = row.i(sampling_pattern);
    const bool barrel = sampling(pattern, false, layer);
    const bool endcap = sampling(pattern, true, layer);
    const float be = row.f(emb1_energy + layer - 1), ee = row.f(eme1_energy + layer - 1);
    const float ba = row.f(emb1_eta + layer - 1), ea = row.f(eme1_eta + layer - 1);
    if (barrel && endcap) {
        const float sum = be + ee;
        if (!(sum > 100.f)) return 0.5 * (ba + ea);
        if ((be > 0 && ee > 0) || (be < 0 && ee < 0)) return (be * ba + ee * ea) / sum;
        return be > 0 ? ba : ea;
    }
    if (barrel) return ba;
    if (endcap) return ea;
    return -999.f;
}

float layer_energy(const Row& row, unsigned layer) {
    const auto pattern = row.i(sampling_pattern);
    const auto b = sampling(pattern, false, layer) ? double(row.f(emb1_energy + layer - 1)) : 0.;
    const auto e = sampling(pattern, true, layer) ? double(row.f(eme1_energy + layer - 1)) : 0.;
    return b + e;
}

Layers layers(const Row& row) {
    return {layer_eta(row, 1), layer_eta(row, 2),
            layer_energy(row, 1) + layer_energy(row, 2) + layer_energy(row, 3)};
}

bool in_window(double value, double low, double high) { return value > low && value < high; }
bool dead_hv_pass(const Row& row) {
    const Int run = row.i(random_run_number);
    if (run < 0 || !row.flag(has_calo_moments)) return false;
    const double ec = row.f(eta_calo), pc = row.f(phi_calo);
    if (run >= 296939 && run < 305291 && in_window(ec, -1.825, -1.5)) {
        if (in_window(pc, -1.104921, -0.956748)) return false;
        const bool second_sector = (run >= 298967 && run < 299144) || run >= 300279;
        if (second_sector && in_window(pc, 1.054922, 1.203097)) return false;
    }
    return !(run >= 482374 && run < 491041 && in_window(ec, -2.050, -1.750) &&
             in_window(pc, 1.914, 2.405));
}

float shower_radius(float eta, bool first) {
    const float a = std::abs(eta);
    float depth;
    if (first) {
        if (a < 0.8) depth = 1558.859292 - 4.990838 * a - 21.144279 * a * a;
        else if (a < 1.5) depth = 1522.775373 + 27.970192 * a - 21.104108 * a * a;
        else depth = 3790.671754;
    } else {
        if (a < 1.425) depth = 1698.990944 - 49.431767 * a - 24.504976 * a * a;
        else if (a < 1.5) depth = 8027.574119 - 2717.653528 * a;
        else depth = 3473.473909 + 453.941515 * a - 119.101945 * a * a;
    }
    return eta < 0. && a > 1.5 ? -depth : depth;
}

std::pair<float,float> convert_rz(float eta, float depth) {
    const float sh = std::sinh(eta);
    if (std::abs(eta) < 1.5) return {depth, depth * sh};
    return {depth / sh, depth};
}

void photon_origin(const Row& row, float eta1, double* out) {
    double corrected_eta = row.f(eta);
    if (std::abs(eta1) < 10) {
        const auto [r, z] = convert_rz(eta1, shower_radius(eta1, true));
        if (r > 0.) corrected_eta = std::asinh((double(z) - row.f(pv_z)) / r);
    }
    // Egamma_v1::e() uses a double PtEtaPhiM4D, followed by float setters.
    const double p = double(row.f(pt)) * std::cosh(double(row.f(eta)));
    const double m = row.f(mass);
    const double e2 = p * p + (m >= 0. ? m * m : -m * m);
    out[0] = static_cast<float>(std::sqrt(e2 > 0. ? e2 : 0.) / std::cosh(corrected_eta));
    out[1] = static_cast<float>(corrected_eta);
    out[20] = 1.;
}

bool char_gate(Int value) {
    if (value < 0 || value > std::numeric_limits<unsigned char>::max())
        throw std::invalid_argument("Reached missing or invalid egamma derivation ID / cleaning char");
    return value != 0;
}

Int electron_track_bits(const Row& row) {
    if (!row.flag(has_track)) return fail_three;
    Int result = fail_three | 1;
    const double sx = row.f(beam_sigma_x), sy = row.f(beam_sigma_y), sxy = row.f(beam_sigma_xy);
    if (sx*sx + sy*sy >= 2*std::abs(sxy) && row.i(track_cov_diag_size) >= 0) {
        if (!row.i(track_cov_diag_size)) throw std::out_of_range("Released egamma d0 lookup found an empty variance vector");
        const double track_angle = row.f(track_phi), s = std::sin(track_angle), c = std::cos(track_angle);
        const double beam_variance = s * (s * (sx*sx) - c*sxy) + c * (c * (sy*sy) - s*sxy);
        const double var = row.f(track_var_d0);
        if (var > 0.) {
            const float significance = row.f(track_d0) / std::sqrt(var + beam_variance);
            if (std::abs(significance) < 5.f) result |= 2;
        }
    }
    const double p = std::abs(double(row.f(pt))), angle = row.f(phi);
    const double x = p * std::cos(angle), y = p * std::sin(angle), z = p * std::sinh(double(row.f(eta)));
    const double theta = x == 0. && y == 0. && z == 0. ? 0. : std::atan2(std::sqrt(x*x + y*y), z);
    const float sum = row.f(track_z0) + row.f(track_vz);
    const double vertex = row.flag(has_pv) ? row.f(pv_z) : 0.;
    const float z0_sin_theta = (sum - vertex) * std::sin(theta);
    if (std::abs(z0_sin_theta) < 0.5f) result |= 4;
    return result;
}

float track_pt(const Row& row, std::size_t qp) {
    const float q_over_p = row.f(qp);
    const float p = std::abs(q_over_p) > 0. ? 1 / std::abs(q_over_p) : 10.e6f;
    const float pxy = p * std::sin(row.f(qp + 1)), angle = row.f(qp + 2);
    const double x = pxy * std::cos(angle), y = pxy * std::sin(angle);
    return std::sqrt(x*x + y*y);
}

bool silicon(Int pixel, Int sct) {
    return static_cast<std::uint8_t>(std::max<Int>(0, pixel) + std::max<Int>(0, sct)) != 0;
}

int photon_category(const Row& row) {
    if (!row.flag(has_conv_trk1)) return 0;
    const bool is_double_si = row.flag(has_conv_trk2) &&
        silicon(row.i(conv_trk1_pixel), row.i(conv_trk1_sct)) &&
        silicon(row.i(conv_trk2_pixel), row.i(conv_trk2_sct));
    if (!is_double_si) return 1;
    const bool stored = row.flag(has_conv_pts);
    const float p1 = stored ? row.f(conv_pt1) : track_pt(row, conv_trk1_q_over_p);
    const float p2 = stored ? row.f(conv_pt2) : track_pt(row, conv_trk2_q_over_p);
    const auto [min, max] = std::minmax(p1, p2);
    if (max == 0.f) return -1;
    const float ratio = min / max;
    const float x = row.f(conv_x), y = row.f(conv_y);
    const float radius = std::sqrt(x*x + y*y);
    if (ratio > 0.3 || (ratio < 0.3 && radius > 140.)) return 1;
    if (radius < 140.) return 2;
    return 0;
}

class PointingProfile final : public PluginNode {
public:
    PointingProfile(const PluginConfig& config, const PluginContext& ctx)
        : depth_(ShowerDepth::load(config.string("depth_axes"))),
          x_(ctx.input_index("rz_x", Type::real)), y_(ctx.input_index("rz_y", Type::real)),
          content_(ctx.child(0)) {
        config.require_keys({"depth_axes"});
        if (ctx.child_count() != 1 || ctx.output_size() != 0 || content_->output_size() != 0)
            throw std::invalid_argument("Expected scalar shower-depth profile child");
    }

    double evaluate(const Values& original) const override {
        const float eta = std::get<double>(original[x_]), phi = std::get<double>(original[y_]);
        return interpolate(original, eta, phi);
    }

private:
    double interpolate(const Values& original, float eta, float phi) const {
        Values input = original;
        const int xb = depth_.x.limit_to_range(eta), yb = depth_.y.limit_to_range(phi);
        if (xb < 1 || xb > depth_.x.bins() || yb < 1 || yb > depth_.y.bins()) return 0.;
        const int xl = depth_.x.lower_center(eta, xb), yl = depth_.y.lower_center(phi, yb);
        const double x0 = depth_.x.center(xl), x1 = depth_.x.center(xl+1);
        const double y0 = depth_.y.center(yl), y1 = depth_.y.center(yl+1);
        auto at = [&](int xi, int yi) {
            input[x_] = depth_.x.center(std::clamp(xi, 1, depth_.x.bins()));
            input[y_] = depth_.y.center(std::clamp(yi, 1, depth_.y.bins()));
            return content_->evaluate(input);
        };
        const double area = (x1-x0)*(y1-y0);
        return at(xl, yl)/area*(x1-eta)*(y1-phi) + at(xl+1, yl)/area*(eta-x0)*(y1-phi)
             + at(xl, yl+1)/area*(x1-eta)*(phi-y0) + at(xl+1, yl+1)/area*(eta-x0)*(phi-y0);
    }

    ShowerDepth depth_;
    std::size_t x_, y_;
    std::shared_ptr<const PluginNode> content_;
};

class LeakageGraph final : public PluginNode {
public:
    LeakageGraph(const PluginConfig& config, const PluginContext& ctx)
        : x_(ctx.input_index("leak_et", Type::real)) {
        config.require_keys({"points"});
        if (ctx.child_count() || ctx.output_size())
            throw std::invalid_argument("Expected scalar egamma leakage graph");
        std::istringstream in(config.string("points"));
        in.imbue(std::locale::classic());
        if (read<std::string>(in) != "MetLeakageGraphV1")
            throw std::invalid_argument("Expected MetLeakageGraphV1");
        const int sorted = read<int>(in);
        const auto size = read<std::size_t>(in);
        if (sorted < 0 || sorted > 1) throw std::invalid_argument("Invalid egamma graph sort bit");
        graph_ = {bool(sorted), std::numeric_limits<double>::infinity(), numbers(in, size), numbers(in, size)};
        in >> std::ws;
        if (!in.eof()) throw std::invalid_argument("Surplus egamma graph fields");
        if (graph_.sorted_bit && !std::is_sorted(graph_.x.begin(), graph_.x.end()))
            throw std::invalid_argument("Unordered TGraph with sorted-X bit");
    }

    double evaluate(const Values& input) const override {
        return graph_.evaluate(std::get<double>(input[x_]));
    }

private:
    std::size_t x_;
    atlas::met_significance::TauGraph graph_;
};

class EgammaNode final : public PluginNode {
public:
    EgammaNode(const PluginConfig& cfg, const PluginContext& ctx)
        : is_photon_(cfg.integer("photon") == 1) {
        cfg.require_keys({"photon", "recipe"});
        if ((cfg.integer("photon") != 0 && cfg.integer("photon") != 1) ||
            cfg.string("recipe") != "25.2.117-fixed-fullsim" || ctx.output_size() != outputs ||
            ctx.child_count() != 3) throw std::invalid_argument("Unsupported MET egamma preparation configuration");
        for (std::size_t i = 0; i < real_.size(); ++i) real_[i] = ctx.input_index(real_names[i], Type::real);
        for (std::size_t i = 0; i < integer_.size(); ++i) integer_[i] = ctx.input_index(int_names[i], Type::integer);
        for (std::size_t i = 0; i < children_.size(); ++i) {
            children_[i] = ctx.child(i);
            if (children_[i]->output_size()) throw std::invalid_argument("Expected scalar MET egamma lookup");
        }
    }

    double evaluate(const Values&) const override {
        throw std::invalid_argument("MET egamma preparation returns a fixed array");
    }

    std::size_t output_size() const override { return outputs; }

    void evaluate_array(const Values& values, double* out, std::size_t size) const override {
        if (size != outputs) throw std::invalid_argument("Invalid MET egamma output extent");
        const Row row{values, real_, integer_};
        std::fill_n(out, size, 0.);
        for (std::size_t i = pt; i <= mass; ++i) out[i] = row.f(i);
        for (std::size_t i = 0; i < 3; ++i) {
            out[13+i] = row.f(iso20+i);
            out[16+i] = row.f(oldleak20+i);
        }
        out[4] = fail_three | 1;
        out[5] = fail_two;
        out[6] = is_photon_ ? fail_one : pass;
        out[7] = fail_one;
        out[8] = is_photon_ ? pass : fail_three;
        out[11] = fail_one;

        if (!row.flag(has_cluster)) return;
        out[4] = fail_three | 3;
        const float limit = is_photon_ ? 2.37f : 2.47f;
        const auto samplings = layers(row);
        if (!(std::abs(samplings.eta2) <= limit)) return;
        out[4] = pass;

        const Int oq_mask = is_photon_ ? 34214 : 1446;
        if (row.i(oq) < 0 || row.i(oq) > pass) throw std::invalid_argument("Missing or invalid egamma OQ");
        out[5] = fail_two | ((row.i(oq) & oq_mask) == 0 ? 1 : 0) | (dead_hv_pass(row) ? 2 : 0);
        if (out[5] != pass) return;
        if (is_photon_) {
            if (!char_gate(row.i(clean_flag))) return;
            out[6] = pass;
            if (row.flag(has_pv)) photon_origin(row, samplings.eta1, out);
        }

        if (!(std::isfinite(out[0]) && out[0] >= 10000.)) return;
        out[7] = out[11] = pass;
        correct_isolation(row, samplings, out);
        if (!is_photon_) {
            out[8] = electron_track_bits(row);
            if (out[8] != pass) return;
        }
        out[9] = char_gate(row.i(id_flag));
        if (!out[9]) return;

        const float track_cut = (is_photon_ ? 0.05 : 0.06) * out[0];
        const float calo_cut = is_photon_ ? float(0.022 * out[0] + 2450.) : float(0.06 * out[0]);
        const unsigned needed = is_photon_ ? 12 : 9;
        out[10] = out[12] = (row.i(iso_available_mask) & needed) == needed &&
            row.f(track_iso) <= track_cut && out[is_photon_ ? 15 : 13] <= calo_cut;
    }

private:
    void correct_isolation(const Row& row, const Layers& samplings, double* out) const {
        const auto which = row.i(iso_available_mask) & row.i(oldleak_available_mask) & 7;
        if (!which || row.i(run_number) < 0) return;
        const float ran_phi = row.f(cluster_phi);
        constexpr std::array edges{0.f, .1f, .6f, .8f, 1.15f, 1.37f, 1.52f, 1.81f, 2.01f, 2.37f, 2.47f};
        const float a = std::abs(samplings.eta2), a1 = std::abs(samplings.eta1);
        const auto found = std::upper_bound(edges.begin(), edges.end(), a);
        const int bin = static_cast<int>(found - edges.begin()) - 1;
        int category = -1;
        if (a1 <= 2.5 && std::abs(ran_phi) <= float(std::numbers::pi) && bin >= 0 && bin < 10 &&
            (!is_photon_ || (bin != 5 && bin != 9))) {
            if (is_photon_) category = photon_category(row);
            else if (row.i(author) == 1) category = 3;
            else if (row.i(author) == 16) category = 4;
        }
        Values child_inputs = row.values;
        if (category >= 0) {
            // IsolationCorrection passes the active isMC=true flag to an isData
            // parameter. The loaded hData is used to test the pointing geometry.
            child_inputs[real_[rz_y]] = double(ran_phi);
            child_inputs[real_[rz_x]] = double(samplings.eta1);
            const float d1 = shower_radius(samplings.eta1, true) - float(children_[2]->evaluate(child_inputs));
            child_inputs[real_[rz_x]] = double(samplings.eta2);
            const float d2 = shower_radius(samplings.eta2, false) - float(children_[2]->evaluate(child_inputs));
            const auto p1 = convert_rz(samplings.eta1, d1);
            const auto p2 = convert_rz(samplings.eta2, d2);
            if (std::abs(p2.first - p1.first) < 1.e-6f) category = -1;
        }

        const double et = samplings.energy / std::cosh(double(samplings.eta2));
        child_inputs[real_[leak_et]] = et;
        const bool lower_trouble = (category == 2 || category == 4) && et < 250000.f;
        const bool timing = row.i(run_number) > 400000;
        if (timing && row.i(timing_available_mask) != 3)
            throw std::invalid_argument("Run 3 MET egamma isolation reached missing Eadded_Lr2 or Eadded_Lr3");
        const float timing_correction = timing ? (row.f(eadded_l2) + row.f(eadded_l3)) /
            std::cosh(samplings.eta2) : 0.f;
        for (int i = 0; i < 3; ++i) {
            if (!(which & (1 << i))) continue;
            float leak = 0.f;
            if (category >= 0) {
                child_inputs[integer_[leak_key]] = Int{i * 50 + category * 10 + bin};
                leak = children_[lower_trouble ? 1 : 0]->evaluate(child_inputs);
            }
            float isolated = row.f(iso20+i);
            isolated += (row.f(oldleak20+i) - leak);
            if (timing) isolated += timing_correction;
            out[13+i] = isolated;
            out[16+i] = leak;
        }
        out[19] = double(which);
    }

    bool is_photon_;
    std::array<std::size_t, real_names.size()> real_{};
    std::array<std::size_t, int_names.size()> integer_{};
    std::array<std::shared_ptr<const PluginNode>,3> children_;
};
} // namespace

void register_egamma_nodes(correction::PluginRegistry& registry) {
    static_assert(correction::plugin_api_version == 4);
    registry.register_node("atlas.met.preparation.leakage_graph", 1,
        [](const PluginConfig& cfg, const PluginContext& ctx) { return std::make_shared<LeakageGraph>(cfg, ctx); });
    registry.register_node("atlas.met.preparation.pointing_profile", 1,
        [](const PluginConfig& cfg, const PluginContext& ctx) { return std::make_shared<PointingProfile>(cfg, ctx); });
    registry.register_node("atlas.met.preparation.egamma", 1,
        [](const PluginConfig& cfg, const PluginContext& ctx) { return std::make_shared<EgammaNode>(cfg, ctx); });
}
} // namespace atlas::met::preparation
