// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "soft.hpp"

#include <correction.h>
#include <array>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace s = atlas::met::soft;

void require(bool yes, const char* message) {
    if (!yes) throw std::runtime_error(message);
}

double decode_number(const std::string& value) {
    char* end = nullptr;
    const double result = std::strtod(value.c_str(), &end);
    require(end != value.c_str() && !*end, "malformed table probe number");
    return result;
}

void check_tables(const char* payload, const char* probes) {
    const auto file = correction::CorrectionSet::from_file(payload);
    std::ifstream stream(payload);
    require(stream.is_open(), "cannot read table payload");
    const std::string json{std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
    const auto text = correction::CorrectionSet::from_string(json.c_str());
    std::ifstream witness(probes);
    require(witness.is_open(), "cannot read the stock table probes");

    std::string name;
    std::array<std::string, 4> fields;
    std::size_t count = 0;
    while (witness >> name >> fields[0] >> fields[1] >> fields[2] >> fields[3]) {
        const std::vector<correction::Variable::Type> coords{
            decode_number(fields[0]), decode_number(fields[1]), decode_number(fields[2])
        };
        const double expected = decode_number(fields[3]);
        for (const auto* evaluator : {file.get(), text.get()}) {
            const double actual = evaluator->at(name)->evaluate(coords);
            if (std::bit_cast<std::uint64_t>(actual) != std::bit_cast<std::uint64_t>(expected)) {
                std::cerr << "Table " << name << " case " << count << " expected "
                          << std::hexfloat << expected << " got " << actual << '\n';
                throw std::runtime_error("ROOT table does not match the native MultiBinning");
            }
        }
        ++count;
    }
    require(witness.eof() && count == 774, "stock table probe coverage");
    std::cout << "stock_table_rows_per_loader=" << count << '\n';
}

int main(int argc, char** argv) try {
    if (argc != 2 && argc != 4) {
        throw std::invalid_argument("Usage: met_soft_smoke CORRECTIONS [TABLES TABLE_PROBES_TSV]");
    }
    if (argc == 4) check_tables(argv[2], argv[3]);
    const s::SoftFile provider(argv[1]);
    const s::SoftFile second(argv[1]);
    constexpr std::array truth{s::Term{"NonInt", 0x21000, {-12000.f, 30000.f, 50000.f}, {}}};
    std::array terms{
        s::Term{"RefJet", 0x10010, {2200.f, -130.f, 15000.f}, 4},
        s::Term{"PVSoftTrk", 0x14020, {210.f, -900.f, 25000.f}, {}},
        s::Term{"Invisible", 0x100000, {100.f, 2.f, 3.f}, {}},
        s::Term{"SoftClus", 0x11020, {5010.f, 0.f, 5010.f}, {}}
    };
    s::EventView event{terms, truth, 46.5f};
    require(provider.apply(event, "JES").status == s::Status::unsupported_state, "unsupported variation");
    const auto nominal = provider.apply(event, "NOSYS");
    require(nominal.status == s::Status::success && !nominal.drew &&
            nominal.seed.truncated < 0 && nominal.momentum.mpx == terms[1].momentum.mpx,
            "nominal input or seed");
    const auto scale = second.apply(event, "MET_SoftTrk_Scale__1up");
    const auto corr = provider.apply(event, "MET_SoftTrk_ResoCorr");
    require(scale.status == s::Status::success && scale.drew && corr.drew &&
            scale.gaussian == corr.gaussian && nominal.seed.effective == corr.seed.effective,
            "scale and correlated state must share the draw");
    require(corr.recoil.mpx == -14300. && corr.recoil.mpy == 30128. && corr.recoil.sumet == 65003.,
            "hard-recoil term selection");
    require(corr.varied_para == corr.reso_para * corr.gaussian &&
            corr.varied_perp == corr.reso_perp * corr.gaussian, "correlated smear");

    event.truth_available = false;
    const auto no_truth = provider.apply(event, "MET_SoftTrk_ResoPara");
    require(no_truth.status == s::Status::success && no_truth.drew &&
            no_truth.recoil.state == s::RecoilState::missing_truth_container &&
            std::isnan(no_truth.momentum.mpx), "missing truth source fallback");
    event.truth_available = true;
    terms[0].constituent_count.reset();
    require(provider.apply(event, "MET_SoftTrk_Scale__1up").status == s::Status::unavailable_jet_links,
            "missing jet aux must fail");
    terms[0].constituent_count = 4;
    const auto total = provider.build_final(terms);
    require(total.status == s::Status::success && total.momentum.mpx == 2410.f
            && total.momentum.mpy == -1030.f && total.momentum.sumet == 40000.f, "ordered Final sum");
    require(total.decorations.phi == provider.decorate(total.momentum).phi, "Final decoration");
    require(s::seed_for({})->effective == 1, "zero seed must map to one");
    require(provider.build_final({}).status == s::Status::missing_soft, "missing soft term");
    std::cout << "MET soft provider smoke passed\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
