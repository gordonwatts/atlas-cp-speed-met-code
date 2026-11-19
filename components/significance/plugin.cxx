/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
#include "significance.hpp"
#include "plugin.hpp"

#include <correction_plugin.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>

namespace atlas::met_significance {
namespace {

using Values = std::vector<correction::Variable::Type>;
using VarType = correction::Variable::VarType;
using correction::PluginConfig;
using correction::PluginContext;

double real(const Values& args, std::size_t index) { return std::get<double>(args.at(index)); }
auto integer(const Values& args, std::size_t index) { return std::get<std::int64_t>(args.at(index)); }

bool flag(const Values& args, std::size_t index) {
    const auto value = integer(args, index);
    if (value != 0 && value != 1) throw std::invalid_argument("MET significance expected a 0 or 1 flag");
    return value != 0;
}

void require_shape(const PluginContext& context, std::size_t size) {
    if (context.child_count() || context.output_size() != size) {
        throw std::invalid_argument("Unexpected MET significance node shape or children");
    }
}

class Scalar final : public correction::PluginNode {
public:
    using Evaluate = std::function<double(const Values&)>;
    explicit Scalar(Evaluate evaluate) : evaluate_(std::move(evaluate)) {}
    double evaluate(const Values& values) const override { return evaluate_(values); }
private:
    Evaluate evaluate_;
};

template<std::size_t Size>
class Array final : public correction::PluginNode {
public:
    using Evaluate = std::function<std::array<double, Size>(const Values&)>;
    explicit Array(Evaluate evaluate) : evaluate_(std::move(evaluate)) {}
    std::size_t output_size() const override { return Size; }
    double evaluate(const Values&) const override {
        throw std::invalid_argument("MET significance node has an array output");
    }
    void evaluate_array(const Values& args, double* output, std::size_t count) const override {
        if (!output || count != Size) throw std::invalid_argument("Incorrect MET significance output buffer");
        const auto result = evaluate_(args);
        std::copy(result.begin(), result.end(), output);
    }
private:
    Evaluate evaluate_;
};

} // namespace

void register_significance_nodes(correction::PluginRegistry& registry) {
    registry.register_node("atlas.met.jet_resolution", 1, [](const PluginConfig& cfg, const PluginContext& ctx) {
        cfg.require_keys({"table"});
        require_shape(ctx, 3);
        const auto pt = ctx.input_index("pt_mev", VarType::real);
        const auto eta = ctx.input_index("eta", VarType::real);
        return std::make_shared<Array<3>>([=, tables = read_jet_data(cfg.string("table"))](const Values& args) {
            const auto value = jet_resolution(tables, real(args, pt), real(args, eta));
            return std::array<double, 3>{value.data, value.mc, value.relative_pt};
        });
    });
    registry.register_node("atlas.met.egamma_resolution", 1, [](const PluginConfig& cfg, const PluginContext& ctx) {
        cfg.require_keys({"table"});
        require_shape(ctx, 0);
        const auto kind = ctx.input_index("kind", VarType::integer);
        const auto energy = ctx.input_index("energy_mev", VarType::real);
        const auto eta = ctx.input_index("cluster_eta", VarType::real);
        const auto eta_calo = ctx.input_index("eta_calo", VarType::real);
        const auto has_eta_calo = ctx.input_index("has_eta_calo", VarType::integer);
        const auto with_ct = ctx.input_index("with_constant_term", VarType::integer);
        return std::make_shared<Scalar>([=, tables = read_egamma_data(cfg.string("table"))](const Values& args) {
            const auto type = integer(args, kind);
            if (type < 0 || type > 2) throw std::invalid_argument("Unsupported MET egamma category");
            const auto calo = flag(args, has_eta_calo)
                ? std::optional<float>(static_cast<float>(real(args, eta_calo))) : std::nullopt;
            return egamma_resolution(tables, static_cast<EgammaKind>(type), real(args, energy), real(args, eta), calo,
                                     flag(args, with_ct));
        });
    });
    registry.register_node("atlas.met.muon_resolution", 1, [](const PluginConfig& cfg, const PluginContext& ctx) {
        cfg.require_keys({"table"});
        require_shape(ctx, 0);
        const auto detector = ctx.input_index("detector", VarType::integer);
        const auto run = ctx.input_index("run_number", VarType::integer);
        const auto is_mc = ctx.input_index("is_simulation", VarType::integer);
        const auto eta = ctx.input_index("eta", VarType::real);
        const auto phi = ctx.input_index("phi", VarType::real);
        const auto cb_pt = ctx.input_index("cb_pt_gev", VarType::real);
        const auto id_pt = ctx.input_index("id_pt_gev", VarType::real);
        const auto ms_pt = ctx.input_index("ms_pt_gev", VarType::real);
        return std::make_shared<Scalar>([=, tables = read_muon_data(cfg.string("table"))](const Values& args) {
            const auto type = integer(args, detector);
            const auto run_number = integer(args, run);
            if (type < 1 || type > 3 || run_number < 0 || run_number > std::numeric_limits<std::uint32_t>::max()) {
                throw std::invalid_argument("Unsupported MET muon detector or run number");
            }
            return muon_resolution(tables, static_cast<MuonDetector>(type), static_cast<std::uint32_t>(run_number),
                                   flag(args, is_mc), real(args, eta), real(args, phi),
                                   {real(args, cb_pt), real(args, id_pt), real(args, ms_pt)});
        });
    });
    registry.register_node("atlas.met.tau_resolution", 1, [](const PluginConfig& cfg, const PluginContext& ctx) {
        cfg.require_keys({"table"});
        require_shape(ctx, 0);
        const auto mode = ctx.input_index("decay_mode", VarType::integer);
        const auto has_mode = ctx.input_index("has_decay_mode", VarType::integer);
        const auto tes_eta = ctx.input_index("tes_eta", VarType::real);
        const auto final_eta = ctx.input_index("final_eta", VarType::real);
        const auto pt = ctx.input_index("final_pt_mev", VarType::real);
        return std::make_shared<Scalar>([=, tables = read_tau_data(cfg.string("table"))](const Values& args) {
            const auto raw_mode = integer(args, mode);
            const auto decay = flag(args, has_mode) && raw_mode >= 0 && raw_mode <= 4
                ? std::optional<int>(static_cast<int>(raw_mode)) : std::nullopt;
            return tau_resolution(tables, decay, static_cast<float>(real(args, tes_eta)),
                                  static_cast<float>(real(args, final_eta)), real(args, pt));
        });
    });
    registry.register_node("atlas.met.significance_finish", 1, [](const PluginConfig& cfg, const PluginContext& ctx) {
        cfg.require_keys({});
        require_shape(ctx, 6);
        const auto longitudinal = ctx.input_index("varl_gev2", VarType::real);
        const auto transverse = ctx.input_index("vart_gev2", VarType::real);
        const auto off_diagonal = ctx.input_index("cov_gev2", VarType::real);
        const auto met = ctx.input_index("met_gev", VarType::real);
        const auto sumet = ctx.input_index("sumet_gev", VarType::real);
        const auto soft_sumet = ctx.input_index("soft_sumet_gev", VarType::real);
        const auto previous_squared = ctx.input_index("previous_squared", VarType::real);
        const auto previous_rho = ctx.input_index("previous_rho", VarType::real);
        return std::make_shared<Array<6>>([=](const Values& args) {
            State state{real(args, previous_squared), real(args, previous_rho)};
            const Covariance cov{real(args, longitudinal), real(args, transverse), real(args, off_diagonal)};
            const auto result = finalize(cov, real(args, met), real(args, sumet), real(args, soft_sumet), state);
            return std::array<double, 6>{result.significance, result.directional, result.met_over_sqrt_sumet,
                result.met_over_sqrt_ht, state.squared_significance, state.rho};
        });
    });
}

} // namespace atlas::met_significance
