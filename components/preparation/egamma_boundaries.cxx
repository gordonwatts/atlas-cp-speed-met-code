// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include <correction_plugin.h>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
using Int = std::int64_t;
using Value = correction::Variable::Type;
using Values = std::vector<Value>;
using Updates = std::initializer_list<std::pair<std::string_view, Value>>;
using Ref = correction::Correction::Ref;
constexpr double pass = 0xffffffff;

void check(bool condition, std::string_view probe) {
    if (!condition) throw std::runtime_error("Egamma boundary probe: " + std::string(probe));
}

class Probe {
public:
    explicit Probe(Ref correction) : correction_(std::move(correction)) {
        const auto& inputs = correction_->inputs();
        inputs_.reserve(inputs.size());
        for (auto& variable : inputs) {
            slots_.emplace(variable.name(), inputs_.size());
            inputs_.emplace_back(variable.type() == correction::Variable::VarType::real ? Value{0.} : Value{Int{0}});
        }
        assign(inputs_, {{"pt", 10000.}, {"has_cluster", Int{1}}, {"has_calo_moments", Int{1}},
                {"sampling_pattern", Int{6}}, {"emb2_energy", 50000.}, {"run_number", Int{410000}},
                {"random_run_number", Int{450000}}, {"id_flag", Int{1}}, {"clean_flag", Int{1}},
                {"iso_available_mask", Int{13}}, {"timing_available_mask", Int{3}},
                {"has_track", Int{1}}, {"track_cov_diag_size", Int{5}}, {"track_var_d0", 1.}});
    }

    std::vector<double> evaluate(Updates edits) const {
        Values row = inputs_;
        assign(row, edits);
        return correction_->evaluate_array(row);
    }

private:
    void assign(Values& row, Updates edits) const {
        for (const auto& [key, value] : edits) row.at(slots_.at(std::string(key))) = value;
    }

    Ref correction_;
    Values inputs_;
    std::unordered_map<std::string, std::size_t> slots_;
};
} // namespace

std::size_t egamma_boundaries(const correction::CorrectionSet& set) {
    const Probe electron(set.at("electron_prepare")), photon(set.at("photon_prepare"));
    const float e_limit = 2.47f, p_limit = 2.37f;
    const double high = std::numeric_limits<float>::infinity();
    check(electron.evaluate({})[12] == 1., "electron inclusive 10 GeV, raw nominal");
    check(photon.evaluate({})[12] == 1. && photon.evaluate({})[20] == 0., "photon inclusive 10 GeV without PriVtx");
    check(electron.evaluate({{"emb2_eta", double(e_limit)}})[4] == pass &&
          electron.evaluate({{"emb2_eta", double(std::nextafter(e_limit, float(high)))}})[4] == 0xfffffffb,
          "electron rounded layer-2 eta boundary");
    check(photon.evaluate({{"emb2_eta", double(p_limit)}})[4] == pass &&
          photon.evaluate({{"emb2_eta", double(std::nextafter(p_limit, float(high)))}})[4] == 0xfffffffb,
          "photon rounded layer-2 eta boundary");
    check(electron.evaluate({{"has_cluster", Int{0}}})[4] == 0xfffffff9, "missing cluster retains cast bit");
    check(electron.evaluate({{"has_calo_moments", Int{0}}})[5] == 0xfffffffd &&
          electron.evaluate({{"random_run_number", Int{-1}}})[5] == 0xfffffffd,
          "missing dead-HV inputs fail only the second OQ bit");
    auto dead_hv = [&](Int run, double eta, double phi) {
        return photon.evaluate({{"random_run_number", run}, {"eta_calo", eta}, {"phi_calo", phi}})[5];
    };
    check(dead_hv(482373, -1.9, 2.) == pass && dead_hv(482374, -1.9, 2.) == 0xfffffffd &&
          dead_hv(491041, -1.9, 2.) == pass, "Run-3 dead-OTx endpoints");
    check(dead_hv(298966, -1.65, 1.1) == pass && dead_hv(298967, -1.65, 1.1) == 0xfffffffd &&
          dead_hv(299144, -1.65, 1.1) == pass, "Run-2 second dead-HV sector");
    check(photon.evaluate({{"clean_flag", Int{0}}})[6] == 0xfffffffe &&
          photon.evaluate({{"clean_flag", Int{0}}})[7] == 0xfffffffe, "photon cleaning prevents downstream gates");
    check(electron.evaluate({{"has_track", Int{0}}})[8] == 0xfffffff8 &&
          electron.evaluate({{"track_d0", 5.}})[8] == 0xfffffffd, "electron missing track and strict d0 edge");
    check(electron.evaluate({{"track_iso", -800.}, {"iso20", 600.}})[10] == 1. &&
          electron.evaluate({{"track_iso", -800.}, {"iso20", 601.}})[10] == 0.,
          "signed electron isolation and inclusive calo threshold");
    check(photon.evaluate({{"track_iso", 500.}, {"iso40", 2670.}})[10] == 1. &&
          photon.evaluate({{"track_iso", 501.}, {"iso40", 2670.}})[10] == 0.,
          "photon inclusive FixedCutTight thresholds");

    auto timing = [&](Int run, Int available, Int timed) {
        return electron.evaluate({{"author", Int{0}}, {"run_number", run}, {"iso20", 100.},
            {"oldleak20", 10.}, {"oldleak_available_mask", available}, {"eadded_l2", 50.},
            {"eadded_l3", -20.}, {"timing_available_mask", timed}});
    };
    check(timing(400000, 1, 0)[13] == 110. && timing(400001, 1, 3)[13] == 140. &&
          timing(400001, 0, 0)[13] == 100. && timing(400001, 0, 0)[19] == 0.,
          "timing transition, unrecognized author and skipped old leakage");
    bool throws = false;
    try {
        timing(400001, 1, 1);
    } catch (const std::invalid_argument&) {
        throws = true;
    }
    check(throws, "missing reached timing field");

    const auto fit = set.at("photon_leakage_fit"), low = set.at("photon_leakage_low");
    auto conv = [&](double energy, double pt2, double radius, Int first_link = 1, double second_eta = 0.) {
        return photon.evaluate({{"emb2_energy", energy}, {"emb2_eta", second_eta},
            {"has_conv_trk1", first_link}, {"has_conv_trk2", Int{1}}, {"conv_trk1_pixel", Int{1}},
            {"conv_trk2_sct", Int{1}}, {"has_conv_pts", Int{1}}, {"conv_x", radius},
            {"conv_pt1", 100.}, {"conv_pt2", pt2}, {"oldleak_available_mask", Int{1}}})[16];
    };
    check(conv(249999., 10., 100.) == float(low->evaluate({249999., Int{20}})) &&
          conv(250000., 10., 100.) == float(fit->evaluate({250000., Int{20}})),
          "photon double-Si trouble graph below 250 GeV and fit at edge");
    check(conv(50000., 31., 100.) == float(fit->evaluate({50000., Int{10}})) &&
          conv(50000., 10., 141.) == float(fit->evaluate({50000., Int{10}})) &&
          conv(50000., 10., 140.) == float(fit->evaluate({50000., Int{0}})),
          "photon conversion ratio and strict 140 mm radius comparisons");
    check(conv(50000., 10., 100., 0) == float(fit->evaluate({50000., Int{0}})) &&
          conv(50000., 10., 100., 1, 1.4) == 0., "invalid first track and crack leakage skip");

    return 17;
}
