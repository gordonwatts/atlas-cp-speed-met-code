// SPDX-License-Identifier: Apache-2.0
#include "met_adapter/significance_correction.hpp"
#include "invocation.hpp"

#include <algorithm>
#include <cmath>
#include <correction_plugin.h>

namespace atlas::met::adapter {
namespace {

using Bound = detail::Invocation<SignificanceEvent, SignificanceResult>;
using K = significance::ObjectKind;

struct ResolutionTables {
    significance::JetData jets;
    significance::EgammaData egamma;
    significance::MuonData muons;
    significance::TauData taus;
};

float relative_pt(const ResolutionTables& tables, const SignificanceEvent& event,
                  const SignificanceLink& link)
{
    switch (link.kind) {
    case K::jet:
        return significance::jet_resolution(tables.jets, link.pt_mev, link.eta).relative_pt;
    case K::electron:
    case K::photon: {
        if (!link.egamma || !link.egamma->cluster_eta) {
            throw std::runtime_error("Selected MET egamma has no first calorimeter cluster");
        }
        const auto& row = *link.egamma;
        const bool photon = link.kind == K::photon;
        // The electron evaluates get_eta_calo in the original interface too.
        const auto calo = significance::calo_frame_eta(
            row.forward_author, *row.cluster_eta, row.eta_calo_moment, row.eta_calo_aux);
        const auto kind = photon
            ? significance::photon_kind(row.vertex_link_count, row.forward_author)
            : significance::EgammaKind::electron;
        return static_cast<float>(significance::egamma_resolution(
            tables.egamma, kind, row.energy_mev, *row.cluster_eta, calo, photon));
    }
    case K::muon: {
        if (!link.muon) throw std::runtime_error("Missing MET muon resolution inputs");
        const auto& row = *link.muon;
        return static_cast<float>(significance::muon_resolution(
            tables.muons, significance::muon_detector(row.type), event.run_number, event.simulation,
            link.eta, link.phi, {significance::track_pt_gev(row.combined),
                                significance::track_pt_gev(row.inner_detector),
                                significance::track_pt_gev(row.extrapolated_ms)}));
    }
    case K::tau: {
        if (!link.tau) throw std::runtime_error("Missing MET tau resolution inputs");
        const auto& row = *link.tau;
        if (!row.decay_mode || *row.decay_mode < 0 || *row.decay_mode > 4) return 1.F;
        if (!row.tes_eta) throw std::runtime_error("Missing etaTauEnergyScale for a selected MET tau");
        if (!(std::abs(*row.tes_eta) < 2.6)) return 1.F;
        if (!row.final_eta || !row.final_pt_mev) {
            throw std::runtime_error("Missing final calibration coordinates for a selected MET tau");
        }
        return static_cast<float>(significance::tau_resolution(
            tables.taus, row.decay_mode, *row.tes_eta, *row.final_eta, *row.final_pt_mev));
    }
    default:
        return 0.F;
    }
}

bool known(K kind)
{
    return kind == K::jet || kind == K::electron || kind == K::photon
        || kind == K::muon || kind == K::tau;
}

SignificanceResult reduce(const ResolutionTables& tables, const SignificanceEvent& event)
{
    const auto total = std::ranges::find(event.terms, std::string_view("Final"), &SignificanceTerm::name);
    const bool ready = total != event.terms.end() && (total->source & 0x10040U) == 0x10040U;
    std::vector<std::vector<significance::Constituent>> rows(event.terms.size());
    std::vector<significance::Term> terms;
    terms.reserve(event.terms.size());
    bool failed_link = false;
    for (std::size_t i = 0; i < event.terms.size(); ++i) {
        const auto& term = event.terms[i];
        auto& links = rows[i];
        if (ready && !failed_link
            && significance::term_action(term.name, term.source) == significance::TermAction::constituents) {
            if (!term.links_available) {
                failed_link = true;
            } else {
                links.reserve(term.links.size());
                for (const auto& link : term.links) {
                    if (!link.valid || !known(link.kind)) {
                        links.push_back({link.kind, 0., 0., 0., 0.F, link.valid});
                        failed_link = true;
                        break;
                    }
                    const float reso = significance::uses_constituent(term.name, link.kind)
                        ? relative_pt(tables, event, link) : 0.F;
                    links.push_back({link.kind, link.pt_mev, link.p4_px_mev, link.p4_py_mev, reso, true});
                }
            }
        }
        terms.push_back({term.name, term.source, term.momentum, links, term.links_available});
    }
    auto state = event.previous;
    const auto result = significance::reduce(terms, state);
    return {result, state};
}

class Node final : public correction::PluginNode {
public:
    Node(std::size_t index, ResolutionTables tables)
        : index_(index), tables_(std::move(tables)) {}

    double evaluate(const std::vector<correction::Variable::Type>& args) const override
    {
        const auto token = std::get<std::int64_t>(args.at(index_));
        if (token == 0) return Bound::probe(identity_);
        auto& invocation = Bound::consume(token, identity_);
        invocation.put(reduce(tables_, invocation.input()));
        return 1.;
    }

private:
    const std::int64_t identity_{detail::next_token()};
    std::size_t index_;
    ResolutionTables tables_;
};

} // namespace

void register_significance_event_node(correction::PluginRegistry& registry)
{
    registry.register_node("atlas.met.significance_event", 1,
        [](const correction::PluginConfig& config, const correction::PluginContext& context) {
            config.require_keys({"jet", "egamma", "muon", "tau"});
            if (context.child_count() || context.output_size() != 0) {
                throw std::invalid_argument("Invalid MET significance event correction definition");
            }
            return std::make_shared<Node>(
                context.input_index("invocation", correction::Variable::VarType::integer),
                ResolutionTables{significance::read_jet_data(config.string("jet")),
                                 significance::read_egamma_data(config.string("egamma")),
                                 significance::read_muon_data(config.string("muon")),
                                 significance::read_tau_data(config.string("tau"))});
        });
}

SignificanceResult invoke_significance(const correction::Correction& evaluator, const SignificanceEvent& event)
{
    if (evaluator.name() != significance_correction_name) {
        throw std::invalid_argument("Expected met_significance_event correction");
    }
    Bound invocation(evaluator, event);
    return invocation.evaluate();
}

} // namespace atlas::met::adapter
