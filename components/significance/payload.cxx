/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
#include "significance.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

namespace atlas::met_significance {
namespace {

class Reader {
public:
    explicit Reader(std::string_view input) : input_(input) {}

    void tag(std::string_view expected) {
        const auto got = next();
        if (got != expected) {
            throw std::invalid_argument("MET significance payload expected " + std::string(expected)
                                        + ", got " + std::string(got));
        }
    }

    double real() { return numeric<double>(); }

    bool flag() {
        const int bit = numeric<int>();
        if (bit != 0 && bit != 1) throw std::invalid_argument("Invalid MET calibration flag");
        return bit != 0;
    }

    std::size_t size() {
        const std::size_t n = numeric<std::size_t>();
        if (!n || n > 4096) throw std::invalid_argument("Invalid MET calibration dimension");
        return n;
    }

    std::vector<double> reals(std::size_t n) {
        std::vector<double> values;
        values.reserve(n);
        for (std::size_t i = 0; i < n; ++i) values.push_back(real());
        return values;
    }

    Axis axis() {
        tag("axis");
        const bool uniform = flag();
        return {uniform, edges()};
    }

    std::vector<double> edges() {
        const auto count = size();
        auto values = reals(count);
        if (count < 2 || !std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); })
            || std::adjacent_find(values.begin(), values.end(), std::greater_equal<double>{}) != values.end()) {
            throw std::invalid_argument("Invalid MET calibration axis edges");
        }
        return values;
    }

    Histogram1 histogram1() {
        auto x = axis();
        auto values = reals(x.bins() + 2);
        return {std::move(x), std::move(values)};
    }

    Histogram2 histogram2() {
        auto x = axis();
        auto y = axis();
        auto values = reals((x.bins() + 2) * (y.bins() + 2));
        return {std::move(x), std::move(y), std::move(values)};
    }

    std::array<double, 2> limits() {
        const std::array<double, 2> limits{real(), real()};
        if (!std::isfinite(limits[0]) || !std::isfinite(limits[1]) || limits[0] >= limits[1]) {
            throw std::invalid_argument("Invalid JER clamp limits");
        }
        return limits;
    }

    void finish() {
        if (!next().empty()) throw std::invalid_argument("Trailing fields in MET significance payload");
    }

private:
    std::string_view input_;

    std::string_view next() {
        constexpr std::string_view whitespace = " \t\r\n";
        const auto start = input_.find_first_not_of(whitespace);
        if (start == std::string_view::npos) return {};
        input_.remove_prefix(start);
        const auto end = input_.find_first_of(whitespace);
        const auto result = input_.substr(0, end);
        input_.remove_prefix(result.size());
        return result;
    }

    template<typename T>
    T numeric() {
        const auto word = next();
        if (word.empty()) throw std::invalid_argument("Truncated MET significance payload");
        T value{};
        const auto end = word.data() + word.size();
        const auto parsed = std::from_chars(word.data(), end, value);
        if (parsed.ec != std::errc{} || parsed.ptr != end) {
            throw std::invalid_argument("Invalid MET significance payload number");
        }
        return value;
    }
};

JetTable read_jet_table(Reader& reader, std::string_view name) {
    reader.tag(name);
    const auto x_limits = reader.limits();
    const auto y_limits = reader.limits();
    auto y = reader.axis();
    std::vector<JetProjection> projections;
    projections.reserve(y.bins() + 1);
    for (std::size_t bin = 0; bin <= y.bins(); ++bin) {
        reader.tag("projection");
        auto h = reader.histogram1();
        auto centers = reader.reals(h.axis.bins() + 2);
        if (!std::all_of(centers.begin(), centers.end(), [](double x) { return std::isfinite(x); })
            || std::adjacent_find(centers.begin() + 1, centers.end() - 1,
                               std::greater_equal<double>{}) != centers.end() - 1) {
            throw std::invalid_argument("Invalid JER projection bin centers");
        }
        projections.push_back({std::move(h), std::move(centers)});
    }
    return {std::move(y), x_limits, y_limits, std::move(projections)};
}

} // namespace

JetData read_jet_data(std::string_view encoded) {
    Reader reader(encoded);
    reader.tag("METSignificanceJetV1");
    auto data = read_jet_table(reader, "data");
    auto mc = read_jet_table(reader, "mc");
    reader.finish();
    return {std::move(data), std::move(mc)};
}

EgammaData read_egamma_data(std::string_view encoded) {
    Reader reader(encoded);
    reader.tag("METSignificanceEgammaV1");
    EgammaData data{};
    reader.tag("eta_edges");
    data.eta_edges = reader.edges();
    const auto bin_count = data.eta_edges.size() - 1;
    constexpr std::array<std::string_view, 3> categories{"electron", "unconverted", "converted"};
    for (std::size_t i = 0; i < categories.size(); ++i) {
        reader.tag(categories[i]);
        for (auto& term : data.sigma90[i]) term = reader.reals(bin_count);
    }
    reader.tag("photon_ct");
    data.constant_term = reader.histogram1();
    reader.finish();
    return data;
}

MuonData read_muon_data(std::string_view encoded) {
    Reader reader(encoded);
    reader.tag("METSignificanceMuonData22V1");
    MuonData data{};
    reader.tag("id_expected");
    for (auto& h : data.id_expected) h = reader.histogram2();
    reader.tag("ms_expected");
    for (auto& h : data.ms_expected) h = reader.histogram2();
    reader.tag("id_nominal");
    for (auto& h : data.id_nominal) h = reader.histogram2();
    reader.tag("ms_nominal");
    for (auto& h : data.ms_nominal) h = reader.histogram2();
    reader.finish();
    return data;
}

TauData read_tau_data(std::string_view encoded) {
    Reader reader(encoded);
    reader.tag("METSignificanceTauV1");
    TauData data{};
    for (auto& mode : data.graphs) {
        for (auto& graph : mode) {
            reader.tag("graph");
            graph.sorted_bit = reader.flag();
            graph.max_x = reader.real();
            const auto size = reader.size();
            graph.x = reader.reals(size);
            graph.y = reader.reals(size);
            if (!std::all_of(graph.x.begin(), graph.x.end(), [](double x) { return std::isfinite(x); })
                || graph.max_x != *std::max_element(graph.x.begin(), graph.x.end())
                || (graph.sorted_bit && !std::is_sorted(graph.x.begin(), graph.x.end()))) {
                throw std::invalid_argument("Invalid tau TGraph abscissa metadata");
            }
        }
    }
    reader.finish();
    return data;
}

} // namespace atlas::met_significance
