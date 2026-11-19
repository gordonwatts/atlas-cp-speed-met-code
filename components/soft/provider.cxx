// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "soft.hpp"
#include "root_rng.hpp"

#include <correction_plugin.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace atlas::met::soft {
namespace {

using Args = std::vector<correction::Variable::Type>;
using VarType = correction::Variable::VarType;
using correction::PluginConfig;
using correction::PluginContext;
using correction::PluginNode;

constexpr std::uint64_t track_soft = 0x14020;
constexpr std::uint64_t jet_mask = 0x10010;
constexpr std::string_view soft_name = "PVSoftTrk";
constexpr std::string_view recipe = "AnalysisBase-25.2.117-PFlow-Tight-PVSoftTrk";

struct Domain {};

struct Invocation {
    const EventView* event = nullptr;
    std::optional<SoftResult>* corrected = nullptr;
    std::span<const Term> terms;
    std::optional<FinalResult>* final = nullptr;
};

class Binding {
public:
    Binding(const Domain& owner, Invocation args) : owner_(&owner), args_(args) {
        if (active_) throw std::logic_error("Nested MET soft event binding");
        if (sequence_ == std::numeric_limits<std::int64_t>::max()) {
            throw std::overflow_error("MET soft event binding ticket exhausted");
        }
        ticket_ = ++sequence_;
        active_ = this;
    }

    ~Binding() { active_ = nullptr; }
    Binding(const Binding&) = delete;
    Binding& operator=(const Binding&) = delete;
    std::int64_t ticket() const { return ticket_; }

    static const Invocation& find(const Domain& owner, std::int64_t ticket) {
        if (!active_ || active_->owner_ != &owner || active_->ticket_ != ticket) {
            throw std::logic_error("Stale or unbound MET soft event ticket");
        }
        return active_->args_;
    }

private:
    const Domain* owner_;
    Invocation args_;
    std::int64_t ticket_{};
    static thread_local std::int64_t sequence_;
    static thread_local const Binding* active_;
};

thread_local std::int64_t Binding::sequence_ = 0;
thread_local const Binding* Binding::active_ = nullptr;

enum class Variation { nominal, scale_up, scale_down, reso_para, reso_perp, reso_corr };

std::optional<Variation> variation(std::string_view name) {
    if (name == "NOSYS") return Variation::nominal;
    if (name == "MET_SoftTrk_Scale__1up") return Variation::scale_up;
    if (name == "MET_SoftTrk_Scale__1down") return Variation::scale_down;
    if (name == "MET_SoftTrk_ResoPara") return Variation::reso_para;
    if (name == "MET_SoftTrk_ResoPerp") return Variation::reso_perp;
    if (name == "MET_SoftTrk_ResoCorr") return Variation::reso_corr;
    return std::nullopt;
}

struct JetInfo {
    Status status{};
    std::size_t size{};
};

JetInfo count_jets(std::span<const Term> terms) {
    const auto jet = std::find_if(terms.begin(), terms.end(), [](const Term& term) {
        return (term.source & jet_mask) == jet_mask;
    });
    if (jet == terms.end()) {
        // The release dereferences find(source) before checking the term pointer.
        return {Status::undefined_jet_iterator};
    }
    if (!jet->constituent_count) return {Status::unavailable_jet_links};
    return {Status::success, *jet->constituent_count};
}

struct XY {
    double x;
    double y;
};

XY project(XY input, const Recoil& recoil) {
    const double ux = recoil.mpx / recoil.pt;
    const double uy = recoil.mpy / recoil.pt;
    return {ux * input.x + uy * input.y, uy * input.x - ux * input.y};
}

using Children = std::array<std::shared_ptr<const PluginNode>, 3>;

class SoftNode final : public PluginNode {
public:
    SoftNode(const PluginContext& ctx, std::shared_ptr<const Domain> domain)
        : domain_(std::move(domain)),
          ticket_(ctx.input_index("ticket", VarType::integer)),
          state_(ctx.input_index("state", VarType::string)),
          pt_(ctx.input_index("pthard", VarType::real)),
          njet_(ctx.input_index("njet", VarType::real)),
          mu_(ctx.input_index("mu", VarType::real)),
          tables_{ctx.child(0), ctx.child(1), ctx.child(2)} {
        for (const auto& child : tables_) {
            if (child->output_size()) throw std::invalid_argument("Expected scalar MET MultiBinning child");
        }
    }

    double evaluate(const Args& args) const override {
        const auto& binding = Binding::find(*domain_, std::get<std::int64_t>(args.at(ticket_)));
        if (!binding.event || !binding.corrected || binding.corrected->has_value()) {
            throw std::logic_error("Wrong or repeated MET soft bound operation");
        }
        return static_cast<double>(binding.corrected->emplace(apply(*binding.event, args)).status);
    }

private:
    SoftResult apply(const EventView& event, const Args& args) const {
        SoftResult result{};
        const auto state = variation(std::get<std::string>(args.at(state_)));
        if (!state) {
            result.status = Status::unsupported_state;
            return result;
        }
        const auto soft = std::find_if(event.terms.begin(), event.terms.end(), [](const Term& term) {
            return term.name == soft_name;
        });
        if (soft == event.terms.end()) {
            result.status = Status::missing_soft;
            return result;
        }
        result.term_index = static_cast<std::size_t>(soft - event.terms.begin());
        result.momentum = soft->momentum;
        const auto seed = seed_for(soft->momentum);
        if (!seed) {
            result.status = Status::undefined_seed;
            return result;
        }
        result.seed = *seed;
        if (!event.event_info_available) {
            result.status = Status::missing_event_info;
            return result;
        }
        if (!event.is_simulation) {
            result.status = Status::not_simulation;
            return result;
        }
        if ((soft->source & track_soft) != track_soft) {
            result.status = Status::unsupported_soft_source;
            return result;
        }
        if (!event.soft_owned) {
            result.status = Status::unowned_soft;
            return result;
        }
        if (*state == Variation::nominal) return result;

        result.recoil = calc_recoil(event);
        const auto jet = count_jets(event.terms);
        if (jet.status != Status::success) {
            result.status = jet.status;
            return result;
        }
        result.jets = jet.size;
        Args coords = args;
        // Use the multi-binning node's overflow cell for +infinity and NaN.
        // ROOT makes this same choice before its special high-pT bin clamp.
        const auto overflow_nan = [](double value) {
            return std::isnan(value) ? std::numeric_limits<double>::infinity() : value;
        };
        coords[pt_] = overflow_nan(result.recoil.pt);
        coords[njet_] = static_cast<double>(jet.size);
        coords[mu_] = overflow_nan(event.actual_mu);
        result.shift = tables_[0]->evaluate(coords);
        result.gaussian = root_rng::first_gaussian(result.seed.effective);
        result.drew = true;

        switch (*state) {
        case Variation::scale_up:
            result.varied_para = result.shift;
            break;
        case Variation::scale_down:
            result.varied_para = -1. * result.shift;
            break;
        case Variation::reso_para:
        case Variation::reso_corr:
            result.reso_para = tables_[1]->evaluate(coords);
            result.varied_para = result.reso_para * result.gaussian;
            if (*state != Variation::reso_corr) break;
            [[fallthrough]];
        case Variation::reso_perp:
            result.reso_perp = tables_[2]->evaluate(coords);
            result.varied_perp = result.reso_perp * result.gaussian;
            break;
        case Variation::nominal:
            break;
        }

        auto rotated = project({soft->momentum.mpx, soft->momentum.mpy}, result.recoil);
        rotated.x += result.varied_para;
        rotated.y += result.varied_perp;
        const auto shifted = project(rotated, result.recoil);
        result.momentum.mpx = static_cast<float>(shifted.x);
        result.momentum.mpy = static_cast<float>(shifted.y);
        return result;
    }

    std::shared_ptr<const Domain> domain_;
    std::size_t ticket_, state_, pt_, njet_, mu_;
    Children tables_;
};

class FinalNode final : public PluginNode {
public:
    FinalNode(const PluginContext& ctx, std::shared_ptr<const Domain> domain)
        : domain_(std::move(domain)), ticket_(ctx.input_index("ticket", VarType::integer)) {}

    double evaluate(const Args& args) const override {
        const auto& binding = Binding::find(*domain_, std::get<std::int64_t>(args.at(ticket_)));
        if (!binding.final || binding.final->has_value()) {
            throw std::logic_error("Wrong or repeated MET Final bound operation");
        }
        return static_cast<double>(binding.final->emplace(ordered_sum(binding.terms)).status);
    }

private:
    std::shared_ptr<const Domain> domain_;
    std::size_t ticket_;
};

class DecorationNode final : public PluginNode {
public:
    explicit DecorationNode(const PluginContext& ctx)
        : x_(ctx.input_index("mpx", VarType::real)), y_(ctx.input_index("mpy", VarType::real)) {}
    std::size_t output_size() const override { return 2; }
    double evaluate(const Args&) const override {
        throw std::invalid_argument("MET decorations require array output");
    }
    void evaluate_array(const Args& args, double* output, std::size_t count) const override {
        if (!output || count != 2) throw std::invalid_argument("Incorrect MET decoration output buffer");
        const auto result = met_phi({
            static_cast<float>(std::get<double>(args.at(x_))),
            static_cast<float>(std::get<double>(args.at(y_))),
            0.f
        });
        output[0] = result.met;
        output[1] = result.phi;
    }
private:
    std::size_t x_, y_;
};

void require_shape(const PluginConfig& cfg, const PluginContext& ctx, std::size_t children,
                   std::size_t outputs) {
    cfg.require_keys({"recipe"});
    if (cfg.string("recipe") != recipe || ctx.child_count() != children || ctx.output_size() != outputs) {
        throw std::invalid_argument("Unsupported MET soft recipe, node shape or child count");
    }
}

void register_nodes(correction::PluginRegistry& registry, const std::shared_ptr<const Domain>& domain) {
    registry.register_node("atlas.met.soft_track", 1, [domain](const PluginConfig& cfg, const PluginContext& ctx) {
        require_shape(cfg, ctx, 3, 0);
        return std::make_shared<SoftNode>(ctx, domain);
    });
    registry.register_node("atlas.met.final_sum", 1, [domain](const PluginConfig& cfg, const PluginContext& ctx) {
        require_shape(cfg, ctx, 0, 0);
        return std::make_shared<FinalNode>(ctx, domain);
    });
    registry.register_node("atlas.met.met_phi", 1, [](const PluginConfig& cfg, const PluginContext& ctx) {
        require_shape(cfg, ctx, 0, 2);
        return std::make_shared<DecorationNode>(ctx);
    });
}

std::unique_ptr<correction::CorrectionSet> load(const std::string& path, const std::shared_ptr<const Domain>& domain) {
    correction::PluginRegistry registry;
    register_nodes(registry, domain);
    return correction::CorrectionSet::from_file(path, registry);
}

template<class Result>
Result confirm_result(double evaluated, std::optional<Result> result) {
    if (!result) throw std::logic_error("MET bound node did not execute");
    if (evaluated != static_cast<double>(result->status)) {
        throw std::logic_error("MET bound node disagrees with its written status");
    }
    return std::move(*result);
}

} // namespace

struct SoftFile::Impl {
    explicit Impl(const std::string& path) : domain(std::make_shared<const Domain>()), set(load(path, domain)),
        correct(set->at("correct_pvsofttrk")), final(set->at("build_final")), decorate(set->at("met_phi")) {}

    std::shared_ptr<const Domain> domain;
    std::unique_ptr<correction::CorrectionSet> set;
    correction::Correction::Ref correct;
    correction::Correction::Ref final;
    correction::Correction::Ref decorate;
};

SoftFile::SoftFile(const std::string& path) : implementation_(std::make_unique<Impl>(path)) {}
SoftFile::~SoftFile() = default;
SoftFile::SoftFile(SoftFile&&) noexcept = default;
SoftFile& SoftFile::operator=(SoftFile&&) noexcept = default;

SoftResult SoftFile::apply(const EventView& input, std::string_view state) const {
    if (!implementation_) throw std::logic_error("Moved-from MET soft provider");
    std::optional<SoftResult> result;
    const Binding binding{*implementation_->domain, {&input, &result, {}, nullptr}};
    const double status = implementation_->correct->evaluate({binding.ticket(), std::string(state), 0., 0., 0.});
    return confirm_result(status, std::move(result));
}

FinalResult SoftFile::build_final(std::span<const Term> terms) const {
    if (!implementation_) throw std::logic_error("Moved-from MET soft provider");
    std::optional<FinalResult> result;
    const Binding binding{*implementation_->domain, {nullptr, nullptr, terms, &result}};
    const double status = implementation_->final->evaluate({binding.ticket()});
    return confirm_result(status, std::move(result));
}

MetPhi SoftFile::decorate(Momentum term) const {
    if (!implementation_) throw std::logic_error("Moved-from MET soft provider");
    std::array<double, 2> output{};
    implementation_->decorate->evaluate_array({double(term.mpx), double(term.mpy)}, output.data(), output.size());
    return {static_cast<float>(output[0]), static_cast<float>(output[1])};
}

} // namespace atlas::met::soft
