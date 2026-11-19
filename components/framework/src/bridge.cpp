// SPDX-License-Identifier: Apache-2.0
#include "bridge.hpp"
#include "met_adapter/maker_correction.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <mutex>

namespace atlas::met::framework {
namespace {
struct Runtime {
    correction::PluginRegistry registry;
    std::unique_ptr<correction::CorrectionSet> maker;
    Runtime()
    {
        const auto provider = absolute_file("ATLAS_MET_PROVIDER");
        const auto payload = absolute_file("ATLAS_MET_MAKER_PAYLOAD");
        registry.load_library(provider);
        maker = correction::CorrectionSet::from_file(payload, registry);
        if (std::getenv("ATLAS_MET_DISPATCH_AUDIT")) {
            std::cerr << "MET_B_COMPONENT provider " << provider << '\n'
                      << "MET_B_COMPONENT maker " << payload << '\n';
        }
    }
};

struct Audit {
    std::mutex mutex;
    std::map<std::string_view, std::size_t> calls;
    ~Audit()
    {
        for (const auto& [operation, count] : calls) {
            std::cerr << "MET_B_DISPATCH " << operation << ' ' << count << '\n';
        }
    }
};
} // namespace

std::string absolute_file(const char* variable)
{
    const char* value = std::getenv(variable);
    if (!value || !std::filesystem::path(value).is_absolute()) {
        throw std::runtime_error(std::string("Set an absolute path for ") + variable);
    }
    return std::filesystem::canonical(value).string();
}

Stage stage()
{
    static const Stage selected = [] {
        const char* option = std::getenv("ATLAS_MET_STAGE");
        const std::string_view value = option ? option : "full";
        if (value == "maker_only") return Stage::maker_only;
#if defined(ATLAS_MET_FULL_BRIDGE)
        if (value == "full") return Stage::full;
#endif
        throw std::runtime_error("Requested MET B stage is absent in the loaded bridge");
    }();
    return selected;
}

const correction::Correction& maker_evaluator()
{
    static const Runtime runtime;
    return *runtime.maker->at(adapter::maker_correction_name);
}

void record_dispatch(std::string_view operation)
{
    static const bool enabled = std::getenv("ATLAS_MET_DISPATCH_AUDIT") != nullptr;
    if (!enabled) return;
    static Audit audit;
    std::lock_guard guard(audit.mutex);
    ++audit.calls[operation];
}

} // namespace atlas::met::framework
