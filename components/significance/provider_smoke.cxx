/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
#include "significance.hpp"

#include <correction_plugin.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace m = atlas::met_significance;

void require(bool match, const char* reason) {
    if (!match) throw std::runtime_error(reason);
}

int main(int argc, char** argv) try {
    if (argc != 3) throw std::invalid_argument("Usage: met_significance_smoke PROVIDER CORRECTIONS");
    correction::PluginRegistry registry;
    registry.load_library(argv[1]);
    auto set = correction::CorrectionSet::from_file(argv[2], registry);
    const auto jet = set->at("jet_resolution")->evaluate_array({40000., 0.4});
    require(jet.size() == 3 && jet[2] == static_cast<float>(std::max(jet[0], jet[1])),
            "wrong jet component layout");
    const auto egamma = set->at("egamma_resolution");
    require(std::isfinite(egamma->evaluate({std::int64_t{0}, 43000., 0.4, 0., std::int64_t{0}, std::int64_t{0}})),
            "electron call did not produce a finite resolution");
    try {
        egamma->evaluate({std::int64_t{1}, 43000., 0.4, 0., std::int64_t{0}, std::int64_t{1}});
        throw std::logic_error("photon accepted without etaCalo");
    } catch (const std::invalid_argument&) {}
    const auto muon = set->at("muon_resolution");
    require(muon->evaluate({std::int64_t{2}, std::int64_t{410000}, std::int64_t{1}, 0.4, 2.0, 0., 0., 0.}) == 1.e12,
            "muon missing-link sentinel differs");
    const auto tau = set->at("tau_resolution");
    require(tau->evaluate({std::int64_t{0}, std::int64_t{0}, 0.4, 0.6, 30000.}) == 1.,
            "tau unavailable mode sentinel differs");
    require(std::isfinite(tau->evaluate({std::int64_t{0}, std::int64_t{1}, 0.4, 0.6, 30000.})),
            "tau with valid mode failed");
    const auto finish = set->at("significance_finish");
    const auto result = finish->evaluate_array({0., 0., 0., 10., 100., 30., 9., 0.4});
    require(result == std::vector<double>{3., -1., 1., 1., 9., 0.4},
            "zero-variance state or HT behavior differs");

    m::State persistent{};
    std::array terms{m::Term{"PVSoftTrk", 0x14020, {2000, 10000, 18000}},
                     m::Term{"Final", 0x10040, {30000, -5000, 60000}}};
    const auto good = m::reduce(terms, persistent);
    require(good.status == m::Status::success && good.decorations.significance > 0,
            "soft-only event failed");
    terms[0].name = "SoftClus";
    const auto no_soft = m::reduce(terms, persistent);
    require(!no_soft.longitudinal_nonzero && no_soft.decorations.significance == good.decorations.significance
            && no_soft.decorations.directional == -1, "zero-covariance event lost stale significance");
    terms[1].source = 0;
    const auto wrong_total = m::reduce(terms, persistent);
    require(wrong_total.final_state == m::FinalState::wrong_source
            && wrong_total.decorations.significance == no_soft.decorations.significance,
            "wrong-source total did not preserve state");
    require(m::reduce({}, persistent).status == m::Status::missing_final, "missing Final did not fail");
    std::cout << "MET significance provider smoke passed\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
