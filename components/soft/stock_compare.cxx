// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "soft.hpp"
#include "root_rng.hpp"

#include <METUtilities/METHelpers.h>
#include <METUtilities/METSystematicsTool.h>
#include <TRandom3.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODMissingET/MissingETAuxContainer.h>
#include <xAODRootAccess/Init.h>
#include <xAODRootAccess/TEvent.h>
#include <xAODRootAccess/TStore.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace s = atlas::met::soft;

namespace {

constexpr std::string_view witness_truth = "METSoftWitnessTruth";
constexpr std::string_view witness_event = "METSoftWitnessEventInfo";
constexpr std::uint64_t soft = 0x14020;
constexpr std::uint64_t jet = 0x10010;

using Links = std::vector<ElementLink<xAOD::IParticleContainer>>;

void require(bool yes, std::string_view message) {
    if (!yes) throw std::runtime_error(std::string(message));
}

template<class T>
bool matches(T actual, T expected) {
    using Bits = std::conditional_t<sizeof(T) == 4, std::uint32_t, std::uint64_t>;
    return (std::isnan(actual) && std::isnan(expected))
        || std::bit_cast<Bits>(actual) == std::bit_cast<Bits>(expected);
}

void same(s::Momentum candidate, const xAOD::MissingET& original, std::string_view label) {
    if (!matches(candidate.mpx, original.mpx()) || !matches(candidate.mpy, original.mpy())
        || !matches(candidate.sumet, original.sumet())) {
        std::cerr << label << " stock=" << std::hexfloat << original.mpx() << ','
                  << original.mpy() << ',' << original.sumet() << " provider="
                  << candidate.mpx << ',' << candidate.mpy << ',' << candidate.sumet() << '\n';
        throw std::runtime_error("Stored soft or Final components differ");
    }
}

xAOD::MissingET& append(xAOD::MissingETContainer& container, const s::Term& term) {
    auto* output = new xAOD::MissingET();
    container.push_back(output);
    output->setName(std::string(term.name));
    output->setSource(term.source);
    output->setMpx(term.momentum.mpx);
    output->setMpy(term.momentum.mpy);
    output->setSumet(term.momentum.sumet);
    if (term.constituent_count) {
        static const SG::Accessor<Links> links("ConstitObjectLinks");
        links(*output).resize(*term.constituent_count);
    }
    return *output;
}

void record_info(xAOD::TStore& store, float mu, bool simulation) {
    auto info = std::make_unique<xAOD::EventInfo>();
    info->makePrivateStore();
    info->setEventTypeBitmask(simulation ? xAOD::EventInfo::IS_SIMULATION : 0u);
    info->setActualInteractionsPerCrossing(mu);
    require(store.record(std::move(info), std::string(witness_event)).isSuccess(), "Record EventInfo");
}

void record_truth(xAOD::TStore& store, std::span<const s::Term> terms) {
    auto aux = std::make_unique<xAOD::MissingETAuxContainer>();
    auto truth = std::make_unique<xAOD::MissingETContainer>();
    truth->setStore(aux.get());
    for (const auto& term : terms) append(*truth, term);
    const std::string key(witness_truth);
    require(store.record(std::move(aux), key + "Aux.").isSuccess(), "Record truth aux");
    require(store.record(std::move(truth), key).isSuccess(), "Record truth");
}

void seed_and_decorations(const s::SoftFile& provider, const xAOD::MissingET& term,
                          const s::SoftResult& got) {
    const auto decorations = provider.decorate({term.mpx(), term.mpy(), term.sumet()});
    require(matches(decorations.met, term.met()) && matches(decorations.phi, term.phi()), "stock term met/phi");
    const int seed = static_cast<int>(1.e6 * term.phi());
    require(got.seed.truncated == seed, "stock signed seed");
    require(got.seed.setter_argument == static_cast<unsigned long>(seed), "stock unsigned-long setter argument");
    const auto effective = static_cast<unsigned long>(seed == 0 ? 1 : seed);
    require(got.seed.effective == static_cast<std::uint32_t>(effective), "stock effective seed");
    if (got.drew) {
        TRandom3 random;
        random.SetSeed(effective);
        require(matches(got.gaussian, random.Gaus(0., 1.)), "stock first Gaussian");
    }
}

std::size_t stage_probes(const s::SoftFile& provider) {
    xAOD::TEvent event(xAOD::TEvent::kClassAccess);
    event.setActive();
    xAOD::TStore store;
    store.setActive();
    met::METSystematicsTool stock("METSoftStageWitness");
    require(stock.setProperty("EventInfo", std::string(witness_event)).isSuccess(), "set EventInfo");
    require(stock.setProperty("TruthContainer", std::string(witness_truth)).isSuccess(), "set truth");
    require(stock.setProperty("ConfigPrefix", std::string("METUtilities/R22_PreRecs/")).isSuccess(), "set prefix");
    require(stock.setProperty("ConfigSoftTrkFile", std::string("TrackSoftTerms-pflow_Dec24.config")).isSuccess(),
            "set TST file");
    require(stock.setProperty("ConfigSoftCaloFile", std::string()).isSuccess(), "disable CST");
    require(stock.setProperty("ConfigJetTrkFile", std::string()).isSuccess(), "disable jet-track systematics");
    require(stock.initialize().isSuccess(), "stock tool initialize");
    constexpr std::array labels{"NOSYS", "MET_SoftTrk_Scale__1up", "MET_SoftTrk_Scale__1down",
        "MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPerp", "MET_SoftTrk_ResoCorr"};
    constexpr std::array recoil_checks{
        s::RecoilState::regular, s::RecoilState::regular, s::RecoilState::small_sumet_floor,
        s::RecoilState::negative_sumet, s::RecoilState::missing_truth_container,
        s::RecoilState::missing_truth_nonint
    };

    std::size_t checks = 0;
    for (std::size_t scenario = 0; scenario < recoil_checks.size(); ++scenario) {
        for (std::string_view label : labels) {
            store.clear();
            const bool nominal = label == "NOSYS";
            const float mu = scenario == 0 ? 46.5f : (scenario == 1 ? 700.f : 0.f);
            record_info(store, mu, true);
            std::vector<s::Term> terms{
                {"RefEle", 0x10001, {9000.f, -7000.f, 15000.f}, {}},
                {"RefJet", jet, {-140.f, 5300.f, 30000.f}, 4},
                {"Invisible", 0x100000, {-512.f, 3502.f, 4000.f}, {}},
                {"PVSoftTrk", soft, {4000.f, -7853.f, 10034.f}, {}},
                {"SoftClus", 0x11020, {4999.f, 0.f, 5000.f}, {}}
            };
            std::array truth{s::Term{"NonInt", 0x21000, {-30000.f, 89000.f, 10000.f}, {}}};
            if (scenario == 1) {
                truth[0].momentum = {8348.f, 1802.f, 10000.f}; // Hard recoil exactly zero, positive sumET.
                terms[3].momentum = {};
            }
            if (scenario == 2) truth[0].momentum = {0.f, 0.f, -49000.f}; // The small-activity floor.
            if (scenario == 3) truth[0].momentum.sumet = -50000.f;
            if (scenario == 5) truth[0].source = 0x22000;
            if (scenario != 4) record_truth(store, truth);

            xAOD::MissingETAuxContainer aux;
            xAOD::MissingETContainer stock_terms;
            stock_terms.setStore(&aux);
            for (const auto& term : terms) append(stock_terms, term);
            auto& stock_soft = *stock_terms["PVSoftTrk"];
            CP::SystematicSet systematic;
            if (!nominal) systematic.insert(CP::SystematicVariation(std::string(label)));
            require(stock.applySystematicVariation(systematic).isSuccess(), label);
            stock.setRandomSeed(static_cast<int>(1.e6 * stock_soft.phi()));

            const auto got = provider.apply({terms, truth, mu, scenario != 4}, label);
            require(got.status == s::Status::success && got.term_index == 3, "provider soft status or index");
            require(got.drew == !nominal, "provider Gaussian reachability");
            require(got.recoil.state == (nominal ? s::RecoilState::unused : recoil_checks[scenario]),
                    "provider recoil state");
            seed_and_decorations(provider, stock_soft, got);
            const xAOD::MissingETAssociationHelper helper(nullptr);
            require(stock.applyCorrection(stock_soft, helper) == CP::CorrectionCode::Ok, "stock correction");
            same(got.momentum, stock_soft, label);
            terms[*got.term_index].momentum = got.momentum;
            const auto final = provider.build_final(terms);
            require(final.status == s::Status::success, "provider Final status");
            require(met::buildMETSum("Final", &stock_terms, stock_soft.source()).isSuccess(), "stock Final sum");
            const auto& stock_final = *stock_terms["Final"];
            require(stock_final.source() == s::final_source, "stock Final source");
            same(final.momentum, stock_final, "Final");
            require(matches(final.decorations.phi, stock_final.phi()) &&
                    matches(final.decorations.met, stock_final.met()), "stock Final decorations");
            ++checks;
        }
    }
    store.clear();
    return checks;
}

std::vector<std::string_view> split(const std::string& line) {
    std::vector<std::string_view> result;
    std::string_view input(line);
    for (;;) {
        const auto i = input.find('\t');
        result.push_back(input.substr(0, i));
        if (i == std::string_view::npos) return result;
        input.remove_prefix(i + 1);
    }
}

double trace_number(std::string_view x) {
    const std::string text(x);
    char* end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    require(end != text.c_str() && !*end, "Invalid trace number");
    return value;
}

struct SeedWitness {
    std::optional<std::uint64_t> seed;
    std::optional<double> gaussian;
    bool state_known = false;
    bool nominal = false;
};

void compare_release_trace(const char* path, std::ostream* witness) {
    std::ifstream source(path);
    require(source.is_open(), "Cannot open stock trace");
    using Key = std::pair<long long, long long>;
    std::map<Key, SeedWitness> pending;
    std::string line;
    std::size_t seed_rows = 0, draw_rows = 0, first = 0;
    while (std::getline(source, line)) {
        const auto parts = split(line);
        if (parts.size() < 4) continue;
        const auto action = parts[2];
        if (action != "setRandomSeed" && action != "gaussian" && action != "applyCorrection"
            && action != "applySystematicVariation") continue;
        const Key key(std::stoll(std::string(parts[0])), std::stoll(std::string(parts[1])));
        auto& row = pending[key];
        if (action == "applySystematicVariation") {
            row.state_known = true;
            row.nominal = parts[3].empty() || parts[3] == "NOSYS";
        } else if (action == "setRandomSeed") {
            require(!row.seed, "Repeated seed in stock trace");
            row.seed = std::stoull(std::string(parts[3]));
        } else if (action == "gaussian") {
            require(!row.gaussian && parts.size() == 6, "Repeated or malformed stock Gaussian trace");
            require(trace_number(parts[3]) == 0 && trace_number(parts[4]) == 1, "Nonstandard stock draw");
            row.gaussian = trace_number(parts[5]);
        } else {
            require(parts.size() == 10 && parts[3] == "2" && row.seed && row.state_known,
                    "Incomplete or unsuccessful stock correction trace");
            const s::Momentum term{static_cast<float>(trace_number(parts[4])),
                static_cast<float>(trace_number(parts[5])), static_cast<float>(trace_number(parts[6]))};
            const auto seed = s::seed_for(term);
            require(seed && seed->setter_argument == *row.seed, "Stock trace seed mismatch");
            require(row.gaussian.has_value() == !row.nominal, "Stock trace unexpected draw reachability");
            ++seed_rows;
            if (row.gaussian) {
                require(matches(s::root_rng::first_gaussian(seed->effective), *row.gaussian), "Stock trace Gaussian mismatch");
                ++draw_rows;
            }
            if (witness) {
                *witness << key.first << '\t' << key.second << '\t' << seed->truncated << '\t'
                         << *row.seed << '\t' << seed->effective << '\t';
                if (row.gaussian) *witness << std::hexfloat << *row.gaussian << std::defaultfloat;
                *witness << '\n';
            }
            if (key.first == 1 && row.nominal) {
                require(seed->truncated == 1158532, "Expected pinned float-phi anchor at entry 1");
                ++first;
            }
            pending.erase(key);
        }
    }
    require(seed_rows && draw_rows && first && pending.empty(), "Stock trace missing required observations");
    std::cout << "stock_seed_rows=" << seed_rows << " stock_draw_rows=" << draw_rows << '\n';
}

} // namespace

int main(int argc, char** argv) try {
    if (argc != 3 && argc != 4) {
        throw std::invalid_argument("Usage: met_soft_stock_compare CORRECTIONS STOCK_STATUS_TSV [SEED_WITNESS_TSV]");
    }
    require(xAOD::Init("met_soft_stock_compare").isSuccess(), "xAOD initialization");
    const s::SoftFile provider(argv[1]);
    const auto cases = stage_probes(provider);
    std::cout << "stock_stage_cases=" << cases << '\n';
    std::ofstream seeds;
    if (argc == 4) {
        require(!std::filesystem::exists(argv[3]), "Choose a fresh seed witness path");
        seeds.open(argv[3]);
        require(seeds.is_open(), "Cannot create seed witness");
        seeds << "entry\tstate_slot\tsigned_int32\tstock_unsigned_long\teffective_uint32\tstock_gaussian\n";
    }
    compare_release_trace(argv[2], argc == 4 ? &seeds : nullptr);
    if (argc == 4) require(seeds.good(), "Seed witness write failed");
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
