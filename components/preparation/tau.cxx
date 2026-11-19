// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "preparation.hpp"
#include "significance.hpp"

#include <correction_plugin.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace atlas::met::preparation {
namespace {

namespace resolution = atlas::met_significance;
using correction::PluginConfig;
using correction::PluginContext;
using Values = std::vector<correction::Variable::Type>;
using VarType = correction::Variable::VarType;

constexpr std::size_t tau_output_size = 9;
constexpr auto tau_recipe = "AnalysisBase-25.2.117-mc23-TauTight-NOSYS";

enum TauReal : std::size_t { pt, eta, phi, mass, charge, pt_final, eta_final, pt_tes, eta_tes };
enum TauInt : std::size_t { tracks, flags, decay_mode, has_tes_pt, tes_compat, truth_pdg,
                            truth_hadronic, truth_jet_available };

constexpr std::array real_fields{"pt", "eta", "phi", "mass", "charge", "pt_final", "eta_final",
                                 "pt_tes", "eta_tes"};
constexpr std::array int_fields{"charged_tracks", "is_tau_flags", "pan_tau_mode", "has_tes_pt",
                                "tes_compatibility", "truth_pdg_id", "truth_is_hadronic",
                                "truth_jet_link_available"};

enum class SmearBranch : int { skipped_truth = 0, skipped_prong = 1, written = 2 };

bool one_or_three(std::int64_t count) { return count == 1 || count == 3; }

std::size_t calo_eta_bin(float eta_tes_value) {
    const float absolute = std::abs(eta_tes_value);
    constexpr std::array limits{0.3, 0.8, 1.3, 1.6, 2.6};
    for (std::size_t i = 0; i < limits.size(); ++i) {
        if (absolute < limits[i]) return i;
    }
    return limits.size();
}

resolution::Axis read_axis(std::istream& reader) {
    int uniform;
    std::size_t count;
    if (!(reader >> uniform >> count) || (uniform != 0 && uniform != 1) || count < 2 || count > 4096) {
        throw std::invalid_argument("Invalid tau scale factor axis shape");
    }
    std::vector<double> edges(count);
    for (auto& edge : edges) {
        if (!(reader >> edge) || !std::isfinite(edge)) throw std::invalid_argument("Tau axis needs finite edges");
    }
    if (std::adjacent_find(edges.begin(), edges.end(), std::greater_equal<double>{}) != edges.end()) {
        throw std::invalid_argument("Unordered tau scale factor edges");
    }
    return {uniform != 0, std::move(edges)};
}

resolution::Histogram2 read_scale(std::string encoded) {
    std::istringstream reader(std::move(encoded));
    reader.imbue(std::locale::classic());
    std::string tag;
    if (!(reader >> tag) || tag != "TauNominalScaleV1") throw std::invalid_argument("Invalid tau scale data");
    auto x = read_axis(reader);
    auto y = read_axis(reader);
    std::vector<double> bins((x.bins() + 2) * (y.bins() + 2));
    for (auto& value : bins) {
        if (!(reader >> value) || !std::isfinite(value)) throw std::invalid_argument("Invalid tau scale bin");
    }
    if (reader >> tag) throw std::invalid_argument("Trailing tau scale data");
    return {std::move(x), std::move(y), std::move(bins)};
}

double scale_factor(const resolution::Histogram2& hist, double input_pt, double input_eta) {
    const double x = input_pt / 1000.;
    const double y = std::abs(input_eta);
    const auto bx = hist.x.root_bin(std::min(std::max(x, hist.x.edges.front()), hist.x.edges.back() * .999));
    const auto by = hist.y.root_bin(std::min(std::max(y, hist.y.edges.front()), hist.y.edges.back() * .999));
    const double value = hist.values.at(bx + (hist.x.bins() + 2) * by);
    // Fixed TF1 title validated for both sf_1p and sf_3p by the exporter.
    // The unclamped GeV coordinate is used by the fading function.
    const double fade = (x < 50.) * 1.0 + ((x >= 50.) && (x < 70.)) * (1 - (x - 50.) / 20.);
    return (value - 1.) * fade + 1.;
}

class TauPrepare final : public correction::PluginNode {
public:
    TauPrepare(const PluginConfig& config, const PluginContext& ctx)
        : sf1_(read_scale(config.string("sf_1p"))), sf3_(read_scale(config.string("sf_3p"))),
          calo_(resolution::read_tau_data(config.string("calo_graphs"))),
          min_pt_gev_(config.real("min_pt_gev")),
          eta_regions_{float(config.real("eta_low")), float(config.real("eta_crack_low")), float(config.real("eta_crack_high")),
                       float(config.real("eta_high"))} {
        for (std::size_t i = 0; i < reals_.size(); ++i) reals_[i] = ctx.input_index(real_fields[i], VarType::real);
        for (std::size_t i = 0; i < ints_.size(); ++i) ints_[i] = ctx.input_index(int_fields[i], VarType::integer);
    }

    std::size_t output_size() const override { return tau_output_size; }
    double evaluate(const Values&) const override {
        throw std::invalid_argument("Tau preparation has fixed array outputs");
    }
    void evaluate_array(const Values& values, double* output, std::size_t size) const override {
        if (!output || size != tau_output_size) throw std::invalid_argument("Invalid tau output extent");
        const auto r = [&](TauReal i) { return static_cast<float>(std::get<double>(values.at(reals_[i]))); };
        const auto n = [&](TauInt i) { return std::get<std::int64_t>(values.at(ints_[i])); };
        if (n(tracks) < 0) throw std::invalid_argument("Unavailable charged tau track count");

        float current_pt = r(pt);
        double sf = 1.;
        SmearBranch branch = SmearBranch::skipped_truth;
        bool compatibility = true;
        if (n(tes_compat) != -1) {
            compatibility = n(tes_compat) != 0;
        } else if (n(has_tes_pt)) {
            const auto mode = n(decay_mode);
            const auto bin = calo_eta_bin(r(eta_tes));
            compatibility = false;
            if (mode >= 0 && mode < 5 && bin < 5) {
                const double et = r(pt_tes);
                const double sigma = et * calo_.graphs.at(mode).at(bin).evaluate(et);
                const double difference = double(r(pt_final)) - et;
                compatibility = !(std::abs(difference) > 5 * sigma);
            }
        }

        const auto pdg = n(truth_pdg);
        const bool is_tau = pdg == 15 || pdg == -15;
        if (is_tau && n(truth_hadronic) < 0) {
            throw std::invalid_argument("Missing IsHadronicTau on a resolved truth tau");
        }
        if (!is_tau && pdg != 11 && pdg != -11 && pdg != 13 && pdg != -13 && !n(truth_jet_available)) {
            throw std::invalid_argument("Missing truthJetLink on the tau fallback path");
        }
        if (is_tau && n(truth_hadronic)) {
            branch = SmearBranch::skipped_prong;
            if (one_or_three(n(tracks))) {
                branch = SmearBranch::written;
                sf = scale_factor(n(tracks) == 1 ? sf1_ : sf3_, r(pt_final), r(eta_final));
                current_pt = static_cast<float>(r(pt_final) * sf);
            }
        }

        // TauSelectionTool follows enum order and returns on the first failed
        // cut. AsgSelectionAlg writes the boolean AcceptData result as char.
        const double absolute_eta = std::abs(double(r(eta)));
        const bool eta_pass = (absolute_eta >= eta_regions_[0] && absolute_eta <= eta_regions_[1])
            || (absolute_eta >= eta_regions_[2] && absolute_eta <= eta_regions_[3]);
        constexpr auto tight_flag = 1ULL << 31;
        unsigned gate_mask = 0;
        if (current_pt / 1000. >= min_pt_gev_) {
            gate_mask |= 1;
            if (eta_pass) {
                gate_mask |= 2;
                if (one_or_three(n(tracks))) {
                    gate_mask |= 4;
                    if (std::abs(r(charge)) == 1.f) {
                        gate_mask |= 8;
                        // TauSelectionTool catches a missing isTauFlags accessor
                        // in this step and retains only the already passed cuts.
                        if (n(flags) >= 0 && (static_cast<std::uint64_t>(n(flags)) & tight_flag)) {
                            gate_mask |= 16 | 32; // JetID, then unconditional ELEIDNONE.
                        }
                    }
                }
            }
        }

        const std::array<double, tau_output_size> result{
            current_pt, r(eta), r(phi), r(mass), gate_mask == 63 ? 1. : 0., double(gate_mask),
            compatibility ? 1. : 0., double(static_cast<int>(branch)), sf
        };
        std::copy(result.begin(), result.end(), output);
    }

private:
    resolution::Histogram2 sf1_, sf3_;
    resolution::TauData calo_;
    double min_pt_gev_;
    std::array<float, 4> eta_regions_;
    std::array<std::size_t, real_fields.size()> reals_{};
    std::array<std::size_t, int_fields.size()> ints_{};
};

class TauTrackCharged final : public correction::PluginNode {
public:
    explicit TauTrackCharged(const PluginContext& ctx)
        : flags_(ctx.input_index("flag_set", VarType::integer)) {}
    double evaluate(const Values& values) const override {
        const auto flags = std::get<std::int64_t>(values.at(flags_));
        if (flags < 0 || flags > 0xffff) throw std::invalid_argument("Missing or invalid TauTrack flagSet");
        return (flags & (1 << 5)) ? 1. : 0.;
    }
private:
    std::size_t flags_;
};

} // namespace

void register_tau(correction::PluginRegistry& registry) {
    registry.register_node("atlas.met.tau_prepare", 1, [](const PluginConfig& cfg, const PluginContext& ctx) {
        cfg.require_keys({"recipe", "sf_1p", "sf_3p", "calo_graphs", "min_pt_gev", "eta_low", "eta_crack_low",
                          "eta_crack_high", "eta_high"});
        if (cfg.string("recipe") != tau_recipe || ctx.child_count() || ctx.output_size() != tau_output_size) {
            throw std::invalid_argument("Unsupported tau preparation payload");
        }
        return std::make_shared<TauPrepare>(cfg, ctx);
    });
    registry.register_node("atlas.met.tau_track_charged", 1, [](const PluginConfig& cfg, const PluginContext& ctx) {
        cfg.require_keys({});
        if (ctx.child_count() || ctx.output_size()) throw std::invalid_argument("Invalid tau track node");
        return std::make_shared<TauTrackCharged>(ctx);
    });
}

} // namespace atlas::met::preparation
