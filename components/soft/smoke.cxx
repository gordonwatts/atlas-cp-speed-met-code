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

void check_wrong_nodes(const char* fixture) {
    const s::SoftFile wrong(fixture);
    const auto rejects_missing_write = [](auto operation) {
        try {
            operation();
        } catch (const std::logic_error& error) {
            return std::string_view(error.what()) == "MET bound node did not execute";
        }
        return false;
    };
    require(rejects_missing_write([&] { wrong.apply({}, "NOSYS"); }),
            "constant-zero soft substitute passed the bound-node guard");
    require(rejects_missing_write([&] { wrong.build_final({}); }),
            "constant-zero Final substitute passed the bound-node guard");
}

void same_soft_result(const s::SoftResult& original, const s::SoftResult& cached) {
    require(original.status == cached.status && original.term_index == cached.term_index &&
            original.drew == cached.drew && original.jets == cached.jets &&
            original.recoil.state == cached.recoil.state, "cached soft status or path");
    require(original.seed.truncated == cached.seed.truncated &&
            original.seed.setter_argument == cached.seed.setter_argument &&
            original.seed.effective == cached.seed.effective, "cached soft seed");
    const auto float_bits = [](s::Momentum q) {
        return std::array{std::bit_cast<std::uint32_t>(q.mpx), std::bit_cast<std::uint32_t>(q.mpy),
                          std::bit_cast<std::uint32_t>(q.sumet)};
    };
    require(float_bits(original.momentum) == float_bits(cached.momentum), "cached soft float momentum");
    const auto double_bits = [](const s::SoftResult& value) {
        std::array output{value.recoil.mpx, value.recoil.mpy, value.recoil.sumet, value.recoil.pt,
                          value.gaussian, value.shift, value.reso_para, value.reso_perp,
                          value.varied_para, value.varied_perp};
        return std::bit_cast<std::array<std::uint64_t,10>>(output);
    };
    require(double_bits(original) == double_bits(cached), "cached soft double stages");
}

void check_cached_first_gaussian(const s::SoftFile& provider, s::EventView input, s::Term& soft,
                                 s::Term& jet) {
    using Use = s::FirstGaussianUse;
    const auto check = [&](s::FirstGaussianCache& cache, std::string_view state, Use expected) {
        auto use = Use::reused; // Verify that each call resets the disposition.
        const auto cached = provider.apply_cached_first_gaussian(input, state, cache, use);
        const auto original = provider.apply(input, state);
        same_soft_result(original,cached);
        require(use == expected && cached.drew == (expected != Use::skipped), "cached soft draw disposition");
        return cached;
    };
    constexpr std::array<std::string_view,6> ordered{"NOSYS", "MET_SoftTrk_Scale__1up",
        "MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPerp", "MET_SoftTrk_ResoCorr", "MET_SoftTrk_Scale__1down"};
    constexpr std::array<std::string_view,6> reordered{"MET_SoftTrk_ResoCorr", "MET_SoftTrk_ResoPerp",
        "NOSYS", "MET_SoftTrk_Scale__1down", "MET_SoftTrk_ResoPara", "MET_SoftTrk_Scale__1up"};
    for (const auto& states : {ordered,reordered}) {
        s::FirstGaussianCache cache;
        bool first = true;
        for (const auto state : states) {
            if (state == "NOSYS") check(cache,state,Use::skipped);
            else {
                check(cache,state,first ? Use::evaluated : Use::reused);
                first = false;
            }
        }
    }

    s::FirstGaussianCache cache;
    const auto before = soft.momentum;
    constexpr auto active = "MET_SoftTrk_ResoPara";
    const auto first = check(cache,active,Use::evaluated);
    soft.momentum.mpy = -soft.momentum.mpy;
    const auto changed = check(cache,active,Use::evaluated);
    require(first.seed.effective != changed.seed.effective, "seed-mismatch fixture must change key");
    check(cache,active,Use::reused);
    soft.momentum = before;
    check(cache,active,Use::evaluated); // The intervening seed replaced the single slot.

    soft.momentum = {1.e6f,0.f,before.sumet};
    const auto zero = check(cache,active,Use::evaluated);
    soft.momentum.mpy = 1.5f;
    const auto one = check(cache,active,Use::reused);
    require(zero.seed.truncated == 0 && one.seed.truncated == 1 &&
            zero.seed.effective == 1 && one.seed.effective == 1,
            "zero fallback and signed one share only the effective seed");
    soft.momentum = before;

    s::FirstGaussianCache next_event;
    input.event_info_available = false;
    require(check(next_event,active,Use::skipped).status == s::Status::missing_event_info,
            "missing EventInfo must skip");
    input.event_info_available = true;
    const auto jet_count = jet.constituent_count;
    jet.constituent_count.reset();
    require(check(next_event,active,Use::skipped).status == s::Status::unavailable_jet_links,
            "missing jet links must skip");
    jet.constituent_count = jet_count;
    check(next_event,active,Use::evaluated);
    require(check(next_event,"JES",Use::skipped).status == s::Status::unsupported_state,
            "unsupported state must skip");
    input.truth_available = false;
    require(check(next_event,active,Use::reused).recoil.state == s::RecoilState::missing_truth_container,
            "missing truth still has a logical draw");
}

int main(int argc, char** argv) try {
    if (argc != 2 && argc != 4 && argc != 5) {
        throw std::invalid_argument("Usage: met_soft_smoke CORRECTIONS [TABLES TABLE_PROBES_TSV [WRONG_NODE]]");
    }
    if (argc >= 4) check_tables(argv[2], argv[3]);
    if (argc == 5) check_wrong_nodes(argv[4]);
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
    check_cached_first_gaussian(provider,event,terms[1],terms[0]);
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
    const auto zero = s::seed_for({});
    require(zero && zero->truncated == 0 && zero->setter_argument == 0 && zero->effective == 1,
            "zero seed must map to one");
    require(nominal.seed.setter_argument == static_cast<std::uint64_t>(nominal.seed.truncated) &&
            nominal.seed.effective == static_cast<std::uint32_t>(nominal.seed.setter_argument),
            "negative seed conversion");
    require(provider.build_final({}).status == s::Status::missing_soft, "missing soft term");
    std::cout << "MET soft provider smoke passed\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
