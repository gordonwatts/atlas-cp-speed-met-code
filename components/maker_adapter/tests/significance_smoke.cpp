// SPDX-License-Identifier: Apache-2.0
#include "met_adapter/significance_correction.hpp"

#include <correction_plugin.h>
#include <array>
#include <iostream>
#include <stdexcept>

namespace a = atlas::met::adapter;
namespace s = atlas::met_significance;

int main(int argc, char** argv) try
{
    if (argc != 4) throw std::invalid_argument("Usage: met_significance_event_smoke PROVIDER EVENT_JSON PRIMITIVES_JSON");
    correction::PluginRegistry registry;
    registry.load_library(argv[1]);
    const auto event_set = correction::CorrectionSet::from_file(argv[2], registry);
    const auto primitives = correction::CorrectionSet::from_file(argv[3], registry);
    const auto operation = event_set->at(a::significance_correction_name);
    const auto egamma = primitives->at("egamma_resolution");

    a::SignificanceLink photon{};
    photon.kind = s::ObjectKind::photon;
    photon.valid = true;
    photon.pt_mev = 45000.;
    photon.p4_px_mev = photon.pt_mev;
    photon.egamma = a::EgammaInput{60000., .4, true, .02, .03F, 1};
    const std::array links{photon};
    const s::Momentum final{30000, 40000, 65000};
    const std::array terms{a::SignificanceTerm{"RefGamma", 0x10002, {}, links, true},
                           a::SignificanceTerm{"Final", 0x10040, final, {}, true}};
    const auto resolved = a::invoke_significance(*operation, {terms, {}, 410000, true, 40.F});

    const auto kind = [](s::EgammaKind k) { return static_cast<std::int64_t>(k); };
    const auto scalar = [&](s::EgammaKind k) {
        return static_cast<float>(egamma->evaluate(
            {kind(k), 60000., .4, static_cast<double>(.4F), std::int64_t{1}, std::int64_t{1}}));
    };
    const float electron = scalar(s::EgammaKind::electron);
    if (electron == scalar(s::EgammaKind::converted_photon)) {
        throw std::runtime_error("Forward-author witness does not distinguish the resolution categories");
    }
    const auto expected = s::object_covariance(final,
        {s::ObjectKind::photon, photon.pt_mev, photon.p4_px_mev, 0., electron});
    const auto got = resolved.reduction.covariance;
    if (resolved.reduction.status != s::Status::success
        || got.longitudinal != expected.longitudinal || got.transverse != expected.transverse
        || got.off_diagonal != expected.off_diagonal) {
        throw std::runtime_error("Forward-author photon did not use electron coefficients with CT");
    }
    std::cout << "MET significance forward-author event smoke passed\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
