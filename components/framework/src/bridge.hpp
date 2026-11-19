// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_FRAMEWORK_BRIDGE_HPP
#define ATLAS_MET_FRAMEWORK_BRIDGE_HPP

#include <AsgTools/AsgTool.h>
#include <AsgMessaging/StatusCode.h>
#include <correction_plugin.h>

#include <dlfcn.h>
#include <stdexcept>
#include <string>
#include <string_view>

class EventContext;

namespace atlas::met::framework {

inline constexpr std::string_view maker_name = "seq1_MetMakerAlg_AnaMET";
inline constexpr std::string_view builder_name = "seq1_MetBuilderAlg_AnaMET";
inline constexpr std::string_view significance_name = "seq1_MetSignificanceAlg_AnaMET";

enum class Stage { maker_only, full };
Stage stage();
std::string absolute_file(const char* variable);
const correction::Correction& maker_evaluator();
void record_dispatch(std::string_view operation);

template<class Method>
Method original_entry(const char* symbol)
{
    dlerror();
    const auto method = reinterpret_cast<Method>(dlsym(RTLD_NEXT, symbol));
    if (const char* error = dlerror()) throw std::runtime_error(error);
    if (!method) throw std::runtime_error(std::string("Unavailable stock symbol: ") + symbol);
    return method;
}

template<class Algorithm>
StatusCode original_execute(Algorithm* algorithm, const EventContext& ctx, const char* symbol)
{
    using Method = StatusCode (*)(Algorithm*, const EventContext&);
    return original_entry<Method>(symbol)(algorithm, ctx);
}

template<class T>
void expect_property(const asg::AsgTool& tool, const char* name, const T& expected)
{
    const T* value = tool.getProperty<T>(name);
    if (!value || *value != expected) {
        throw std::runtime_error(std::string("Unsupported initialized property ") + tool.name()
                                 + "." + name);
    }
}

template<class Interface>
const asg::AsgTool& asg_tool(const Interface* tool)
{
    const auto* asg_tool = dynamic_cast<const asg::AsgTool*>(tool);
    if (!asg_tool) throw std::runtime_error("Expected an initialized standalone asg::AsgTool");
    return *asg_tool;
}

} // namespace atlas::met::framework
#endif
