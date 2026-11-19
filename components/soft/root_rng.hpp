/*
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Behavior adapted from ROOT v6-40-02, commit d4c0003542952bec9c4edc86c419784a32875e9d.
 * ROOT TRandom credits Rene Brun and Fons Rademakers.
 * The original MT19937 implementation credits Makoto Matsumoto and Takuji Nishimura.
 * ROOT TRandom3 credits its CLHEP adaptation to Peter Malzacher.
 * The ACR method is due to Wolfgang Hörmann and Gerhard Derflinger.
 * The ROOT source credits UNU.RAN (Wolfgang Hörmann and Josef Leydold, WU Wien).
 * Adapted from the 2026-11-17 muon provider ROOT compatibility source. See NOTICE.md.
 */
#pragma once

#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>

namespace atlas::met::soft::root_rng {

template<class NextWord>
double uniform_from_words(NextWord& next) {
    std::uint32_t word;
    do {
        word = static_cast<std::uint32_t>(next());
    } while (word == 0);
    return static_cast<double>(word) * 0x1p-32;
}

class Uniforms {
public:
    explicit Uniforms(std::uint32_t seed) : engine_(nonzero(seed)) {}

    double operator()() { return uniform_from_words(engine_); }

private:
    static std::uint32_t nonzero(std::uint32_t seed) {
        if (seed == 0) {
            throw std::domain_error("MET soft wants a nonzero effective ROOT seed");
        }
        return seed;
    }

    std::mt19937 engine_;
};

// The constants and evaluation order reproduce TRandom::Gaus(0, 1).
// A callable uniform source makes branch consumption observable in the probe.
template<class Uniform>
double gaussian(Uniform& uniform) {
    constexpr double upper_region = 0.516058551;
    constexpr double lower_region = 0.107981933;
    constexpr double central_region = 0.483941449;
    constexpr double log_pi_over_two = 0.4515827053;
    const auto finish = [](double value) { return 0.0 + 1.0 * value; };

    const double pick = uniform();
    if (pick > upper_region) return finish(4.132731354 * pick - 3.132731354);
    if (pick < lower_region) {
        const double shifted = 18.52161694 * pick - 1;
        return finish(shifted > 0 ? 1 + shifted : -1 + shifted);
    }

    if (pick < central_region) {
        const double draw = uniform();
        const double centered = draw - 1 + draw;
        const double reflected = centered > 0 ? 2 - centered : -2 - centered;
        if ((1.448242853 - pick) * (1.46754004 + std::abs(reflected)) < 3.307147487) {
            return finish(reflected);
        }
        const double square = centered * centered;
        if ((pick + 1.036467755) * (3.631288474 + square) < 5.295844968) {
            return finish(centered);
        }
        if (0.591923442 - pick < std::exp(-(reflected * reflected + log_pi_over_two) / 2)) {
            return finish(reflected);
        }
        if (pick + 0.375959516 < std::exp(-(square + log_pi_over_two) / 2)) {
            return finish(centered);
        }
    }

    constexpr double tail_height = 0.187308492;
    for (;;) {
        double x = uniform();
        double y = tail_height * uniform();
        const bool positive = 0.4571828819 - 0.7270572718 * x - y > 0;
        if (!positive) {
            x = 1 - x;
            y = tail_height - y;
        }
        const double magnitude = 2 + y / x;
        const double candidate = positive ? magnitude : -magnitude;
        if ((y - 0.8853395638 + x) * (0.2770276848 + x) + 0.2452635696 < 0) {
            return finish(candidate);
        }
        if (y < x + 0.03895759111 && candidate * candidate < 4 * (0.5029324303 - std::log(x))) {
            return finish(candidate);
        }
    }
}

inline double first_gaussian(std::uint32_t seed) {
    Uniforms uniform(seed);
    return gaussian(uniform);
}

} // namespace atlas::met::soft::root_rng
