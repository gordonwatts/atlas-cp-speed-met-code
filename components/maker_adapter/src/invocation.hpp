// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_INVOCATION_HPP
#define ATLAS_MET_INVOCATION_HPP

#include <correction.h>

#include <atomic>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace atlas::met::adapter::detail {

inline std::int64_t next_token()
{
    static std::atomic<std::uint64_t> next{1};
    const auto token = next.fetch_add(1, std::memory_order_relaxed);
    if (token > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("MET invocation token space exhausted");
    }
    return static_cast<std::int64_t>(token);
}

template <class Input, class Result>
class Invocation {
public:
    Invocation(const correction::Correction& expected, const Input& input)
        : expected_(&expected), input_(input), token_(next_token()),
          node_(identify(expected)), previous_(current_)
    {
        current_ = this;
    }

    ~Invocation() { current_ = previous_; }
    Invocation(const Invocation&) = delete;
    Invocation& operator=(const Invocation&) = delete;

    // The zero input is reserved for the library's identity probe. Direct
    // evaluators have no probe or event binding, including for the zero input.
    static double probe(std::int64_t node)
    {
        if (!probe_ || *probe_) throw std::logic_error("MET node identity probe is absent or reused");
        *probe_ = node;
        return acknowledgment;
    }

    static Invocation& consume(std::int64_t token, std::int64_t node)
    {
        if (!current_ || current_->token_ != token || current_->node_ != node || current_->consumed_) {
            throw std::logic_error("MET binding is absent, stale, reused, or belongs to another correction");
        }
        current_->consumed_ = true;
        return *current_;
    }

    [[nodiscard]] std::int64_t token() const { return token_; }
    [[nodiscard]] const Input& input() const { return input_; }
    void put(Result result) { result_.emplace(std::move(result)); }

    Result evaluate()
    {
        if (expected_->evaluate({token_}) != acknowledgment || !consumed_ || !result_) {
            throw std::logic_error("MET correction did not execute the bound node exactly once");
        }
        return std::move(*result_);
    }

private:
    static constexpr double acknowledgment = 1.;
    static std::int64_t identify(const correction::Correction& expected)
    {
        std::optional<std::int64_t> identity;
        struct ProbeScope {
            std::optional<std::int64_t>* previous;
            ~ProbeScope() { probe_ = previous; }
        } scope{std::exchange(probe_, &identity)};
        if (expected.evaluate({std::int64_t{0}}) != acknowledgment || !identity) {
            throw std::logic_error("MET correction did not identify a registered event node");
        }
        return *identity;
    }

    const correction::Correction* expected_;
    const Input& input_;
    const std::int64_t token_;
    const std::int64_t node_;
    Invocation* const previous_;
    bool consumed_{};
    std::optional<Result> result_;
    static inline thread_local Invocation* current_{};
    static inline thread_local std::optional<std::int64_t>* probe_{};
};

} // namespace atlas::met::adapter::detail
#endif
