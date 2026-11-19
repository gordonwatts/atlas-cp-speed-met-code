// SPDX-License-Identifier: Apache-2.0
#include "maker_trace.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace atlas::met::framework {
namespace {
const char* path()
{
    static const char* selected = std::getenv("ATLAS_MET_MAKER_TRACE_FILE");
    return selected;
}

std::ofstream& stream()
{
    static auto output = [] {
        const std::filesystem::path trace(path());
        if (!trace.is_absolute() || std::filesystem::exists(trace)) {
            throw std::runtime_error("MET maker trace requires an absolute new file");
        }
        std::ofstream output(trace);
        if (!output) throw std::runtime_error("Cannot create MET maker trace");
        return output;
    }();
    return output;
}

std::string_view term_label(maker::TermType term)
{
    using T = maker::TermType;
    switch (term) {
    case T::electron: return "RefEle";
    case T::photon: return "RefGamma";
    case T::tau: return "RefTau";
    case T::muon: return "Muons";
    default: return {};
    }
}
} // namespace

bool maker_trace_enabled()
{
    return path() && *path();
}

std::size_t next_maker_trace_entry()
{
    static std::size_t entry{};
    return entry++;
}

void trace_maker(std::size_t entry, std::size_t state, const maker::EventInput& input,
                 const maker::Result& result, const MakerInputs& owning)
{
    auto& output = stream();
    auto prefix = [&](std::string_view action) -> std::ostream& {
        return output << entry << '\t' << state << '\t' << action << '\t';
    };
    for (const auto& flag : result.flag_writes) {
        prefix("associationFlag") << flag.association << '\t' << flag.slot << '\t' << flag.selected << '\n';
        const auto original = owning.current_link(flag.original);
        prefix("associationActor") << flag.association << '\t' << flag.slot << '\t'
            << term_label(flag.term) << '\t' << original.key() << '\t' << original.index() << '\n';
    }
    for (const auto& jet : result.jet_decisions) {
        if (jet.jvt_pass) prefix("jetJvt") << jet.current.index << '\t' << *jet.jvt_pass << '\n';
    }
    if (result.status == maker::Status::success) {
        prefix("rebuildMET") << "RefEle\t" << input.electrons.selected.size() << "\t1\n";
        prefix("rebuildMET") << "RefGamma\t" << input.photons.selected.size() << "\t1\n";
        prefix("rebuildMET") << "RefTau\t" << input.taus.selected.size() << "\t1\n";
        prefix("rebuildMET") << "Muons\t" << input.muons.selected.size() << "\t1\n";
        prefix("rebuildJetMET") << "RefJet\tPVSoftTrk\t" << input.jets.size() << "\t1\t1\n";
    } else {
        prefix("makerFailure") << static_cast<int>(result.status) << '\n';
    }
    if (!output) throw std::runtime_error("Cannot write MET maker trace");
}

} // namespace atlas::met::framework
