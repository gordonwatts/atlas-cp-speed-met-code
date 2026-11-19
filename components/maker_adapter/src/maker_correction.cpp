// SPDX-License-Identifier: Apache-2.0
#include "met_adapter/maker_correction.hpp"
#include "invocation.hpp"

#include <correction_plugin.h>
#include <locale>
#include <sstream>

namespace atlas::met::adapter {
namespace {
using Bound = detail::Invocation<maker::EventInput, maker::Result>;

class MakerNode final : public correction::PluginNode {
public:
    MakerNode(std::size_t index, maker::Config config)
        : index_(index), config_(std::move(config)) {}

    double evaluate(const std::vector<correction::Variable::Type>& args) const override
    {
        const auto token = std::get<std::int64_t>(args.at(index_));
        if (token == 0) return Bound::probe(identity_);
        auto& invocation = Bound::consume(token, identity_);
        invocation.put(maker::evaluate(config_, invocation.input()));
        return 1.;
    }

private:
    const std::int64_t identity_{detail::next_token()};
    std::size_t index_;
    maker::Config config_;
};

maker::NnjvtCutMap read_cuts(const std::string& table)
{
    std::istringstream input(table);
    input.imbue(std::locale::classic());
    std::size_t pt{}, eta{};
    if (!(input >> pt >> eta) || pt < 2 || eta < 2 || pt > 1024 || eta > 1024) {
        throw std::invalid_argument("Invalid MET maker FixedEffPt table dimensions");
    }
    maker::NnjvtCutMap cuts{std::vector<float>(pt), std::vector<float>(eta),
                           std::vector<float>((pt - 1) * (eta - 1))};
    for (auto* block : {&cuts.pt_edges_mev, &cuts.eta_edges, &cuts.cuts}) {
        for (auto& value : *block) {
            if (!(input >> value)) {
                throw std::invalid_argument("Missing value in MET maker FixedEffPt table");
            }
        }
    }
    std::string extra;
    if (input >> extra) {
        throw std::invalid_argument("Extra values in MET maker FixedEffPt table");
    }
    return cuts;
}
} // namespace

void register_maker_nodes(correction::PluginRegistry& registry)
{
    registry.register_node("atlas.met.build_met_terms", 1,
        [](const correction::PluginConfig& config, const correction::PluginContext& context) {
            config.require_keys({"recipe", "source_sha256", "cuts",
                "min_pt_mev", "max_pt_mev", "min_abs_eta", "max_abs_eta"});
            if (context.child_count() || context.output_size() != 0
                || config.string("recipe") != maker::Config::recipe
                || config.string("source_sha256").size() != 64) {
                throw std::invalid_argument("Invalid MET maker correction definition");
            }
            const auto token = context.input_index("invocation", correction::Variable::VarType::integer);
            const maker::JvtRange jvt{
                static_cast<float>(config.real("min_pt_mev")),
                static_cast<float>(config.real("max_pt_mev")),
                static_cast<float>(config.real("min_abs_eta")),
                static_cast<float>(config.real("max_abs_eta"))};
            return std::make_shared<MakerNode>(
                token, maker::Config(read_cuts(config.string("cuts")), jvt));
        });
}

maker::Result invoke_maker(const correction::Correction& evaluator, const maker::EventInput& event)
{
    if (evaluator.name() != maker_correction_name) {
        throw std::invalid_argument("Expected the build_met_terms correction");
    }
    Bound invocation(evaluator, event);
    return invocation.evaluate();
}

} // namespace atlas::met::adapter
