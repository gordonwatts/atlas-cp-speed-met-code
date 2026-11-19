// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_SOFT_HPP
#define ATLAS_MET_SOFT_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>

namespace atlas::met::soft {

struct Momentum {
    float mpx{};
    float mpy{};
    float sumet{};
};

struct Term {
    std::string_view name;
    std::uint64_t source{};
    Momentum momentum;
    // Needed only for the first source matching the stock jet mask.
    std::optional<std::size_t> constituent_count;
};

struct EventView {
    std::span<const Term> terms;
    std::span<const Term> truth;
    float actual_mu{};
    bool truth_available{true};
    bool event_info_available{true};
    bool is_simulation{true};
    bool soft_owned{true};
};

enum class Status {
    success,
    unsupported_state,
    missing_soft,
    missing_event_info,
    not_simulation,
    unsupported_soft_source,
    unowned_soft,
    undefined_seed,
    undefined_jet_iterator,
    unavailable_jet_links,
    duplicate_final,
};

enum class RecoilState { unused, regular, small_sumet_floor, missing_truth_container,
                         missing_truth_nonint, negative_sumet };

struct Recoil {
    RecoilState state{};
    double mpx{};
    double mpy{};
    double sumet{};
    double pt{};
};

struct Seed {
    std::int32_t truncated{};
    std::uint64_t setter_argument{}; // Linux unsigned-long ABI, before zero handling.
    std::uint32_t effective{};
};

struct SoftResult {
    Status status{};
    std::optional<std::size_t> term_index;
    Momentum momentum;
    Seed seed;
    Recoil recoil;
    std::size_t jets{};
    bool drew{};
    double gaussian{};
    double shift{};
    double reso_para{};
    double reso_perp{};
    double varied_para{};
    double varied_perp{};
};

inline constexpr std::string_view final_name = "Final";
inline constexpr std::uint64_t final_source = 0x10040;

struct MetPhi {
    float met{};
    float phi{};
};

struct FinalResult {
    Status status{};
    Momentum momentum;
    MetPhi decorations;
};

// The stored components are MeV and the public decorations are floats.
MetPhi met_phi(Momentum);
std::optional<Seed> seed_for(Momentum);
Recoil calc_recoil(const EventView&);
FinalResult ordered_sum(std::span<const Term> terms);

// Immutable registry and bound event domain. These public calls execute
// CorrectionLib v4 nodes. All input views need only survive the current call.
class SoftFile {
public:
    explicit SoftFile(const std::string& corrections_path);
    ~SoftFile();
    SoftFile(SoftFile&&) noexcept;
    SoftFile& operator=(SoftFile&&) noexcept;

    SoftResult apply(const EventView&, std::string_view state) const;
    FinalResult build_final(std::span<const Term> terms) const;
    MetPhi decorate(Momentum) const;

private:
    struct Impl;
    std::unique_ptr<Impl> implementation_;
};

} // namespace atlas::met::soft
#endif
