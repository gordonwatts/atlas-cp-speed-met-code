// SPDX-License-Identifier: Apache-2.0
#include "met_adapter/maker_correction.hpp"
#include "invocation.hpp"

#include <correction_plugin.h>
#include <array>
#include <iostream>
#include <stdexcept>

namespace m = atlas::met::maker;
namespace a = atlas::met::adapter;

template<class Function>
void expect_rejected(Function&& action, const char* reason)
{
    try {
        action();
    } catch (const std::logic_error&) {
        return;
    }
    throw std::runtime_error(reason);
}

int main(int argc, char** argv) try
{
    if (argc != 3) throw std::invalid_argument("Usage: maker_adapter_smoke PROVIDER PAYLOAD");
    correction::PluginRegistry registry;
    registry.load_library(argv[1]);
    const auto set = correction::CorrectionSet::from_file(argv[2], registry);
    const auto other = correction::CorrectionSet::from_file(argv[2], registry);
    const auto op = set->at(a::maker_correction_name);
    const auto other_op = other->at(a::maker_correction_name);

    m::EventInput event{};
    if (a::invoke_maker(*op, event).status != m::Status::failure) {
        throw std::runtime_error("Missing core did not propagate failure");
    }
    expect_rejected([&] { (void)op->evaluate({std::int64_t{0}}); },
                    "Unbound identity probe accepted");

    const std::array core{m::CoreTerm{"PVSoftTrkCore", 0x14020, {}}};
    event.core_terms = core;
    using Binding = a::detail::Invocation<m::EventInput, m::Result>;
    std::int64_t stale{};
    {
        Binding outer(*op, event);
        stale = outer.token();
        expect_rejected([&] { (void)other_op->evaluate({outer.token()}); },
                        "Active token accepted by a different CorrectionSet");
        {
            Binding inner(*other_op, event);
            expect_rejected([&] { (void)op->evaluate({outer.token()}); },
                            "Outer binding accepted during nested invocation");
            if (inner.evaluate().status != m::Status::success) throw std::runtime_error("Nested evaluate failed");
            expect_rejected([&] { (void)other_op->evaluate({inner.token()}); },
                            "Consumed nested invocation was reusable");
        }
        const auto value = outer.evaluate();
        if (value.status != m::Status::success || value.terms.size() != 7) {
            throw std::runtime_error("Valid outer event lost its own result");
        }
    }
    expect_rejected([&] { (void)op->evaluate({stale}); }, "Stale invocation accepted");
    std::cout << "MET maker shared invocation smoke passed\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
