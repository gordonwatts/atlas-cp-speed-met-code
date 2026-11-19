// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "met_maker/maker.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace atlas::met::maker {
namespace {

constexpr Mask refined = 0x10000;
constexpr Mask central = 0x00100;
constexpr Mask track_source = refined | 0x04000;
constexpr std::uint8_t pflow_overlap = 1U << 3;
constexpr std::uint8_t primary_vertex = 1;
constexpr std::uint16_t calo_tagged_muon = 3;
constexpr std::size_t no_slot = std::numeric_limits<std::size_t>::max();

struct Stop {
    Status status;
    Issue issue;
    std::string_view detail{};
    Severity severity{Severity::error};
};

bool same_target(const ObjectLink& link, ObjectId id)
{
    return link.state == LinkState::valid && link.id == id;
}

std::size_t find_slot(const Association& association, ObjectId id)
{
    const auto it = std::ranges::find_if(association.objects, [id](const ObjectLink& link) {
        return same_target(link, id);
    });
    return it == association.objects.end()
        ? no_slot : static_cast<std::size_t>(it - association.objects.begin());
}

// The release setter, per-object sum, and muon mask use a signed int shift.
Mask stock_signed_mask(std::size_t slot)
{
    if (slot >= 32) {
        throw Stop{Status::stock_undefined, Issue::stock_signed_shift_out_of_range,
                   "AnalysisBase shifts an int by this slot index"};
    }
    return static_cast<Mask>(static_cast<std::int32_t>(std::uint32_t{1} << slot));
}

Mask stock_unsigned_mask(std::size_t slot)
{
    if (slot >= 64) {
        throw Stop{Status::stock_undefined, Issue::stock_mask_shift_out_of_range,
                   "AnalysisBase shifts a uint64_t by this slot index"};
    }
    return Mask{1} << slot;
}

template <class T>
const T& require_aux(const std::optional<T>& value, std::string_view name)
{
    if (!value) {
        throw Stop{Status::source_access_exception, Issue::missing_required_aux, name};
    }
    return *value;
}

void add(ConstituentVector& target, const ConstituentVector& v)
{
    target.px += v.px;
    target.py += v.py;
    target.pz += v.pz;
    target.energy += v.energy;
    target.sum_pt += v.sum_pt;
}

void subtract(ConstituentVector& target, const ConstituentVector& v)
{
    target.px -= v.px;
    target.py -= v.py;
    target.pz -= v.pz;
    target.energy -= v.energy;
    target.sum_pt -= v.sum_pt;
}

void scale(ConstituentVector& target, float factor)
{
    target.px *= factor;
    target.py *= factor;
    target.pz *= factor;
    target.energy *= factor;
    target.sum_pt *= factor;
}

ConstituentVector reduce(std::span<const KeyedVector> rows, Mask mask, bool invert = false)
{
    ConstituentVector result{};
    for (const KeyedVector& row : rows) {
        if (((mask & row.key) != 0) != invert) {
            add(result, row.vector);
        }
    }
    return result;
}

void add_particle(Term& term, float px, float py, float pt)
{
    term.momentum.mpx -= px;
    term.momentum.mpy -= py;
    term.momentum.sumet += pt;
}

void add_met(Term& term, const MetVector& value)
{
    term.momentum.mpx += value.mpx;
    term.momentum.mpy += value.mpy;
    term.momentum.sumet += value.sumet;
}

// ROOT 6.40.02 JetFourMom_t coordinates. The Cartesian-to-polar round trip
// for alternate scale energy preserves the source momentum representation.
constexpr double eta_max = 22756.0;

double normalized_phi(double phi)
{
    constexpr double pi = std::numbers::pi_v<double>;
    if (phi <= -pi || phi > pi) {
        phi = phi - std::floor(phi / (2 * pi) + .5) * 2 * pi;
    }
    return phi;
}

double spatial_magnitude(double pt, double eta)
{
    if (pt > 0) {
        return pt * std::cosh(eta);
    }
    if (eta > eta_max) {
        return eta - eta_max;
    }
    if (eta < -eta_max) {
        return -eta - eta_max;
    }
    return 0.;
}

double jet_energy(const PtEtaPhiM& p, bool restrict_negative_mass)
{
    const double mag = spatial_magnitude(p.pt, p.eta);
    const double p2 = mag * mag;
    double mass = p.mass;
    if (restrict_negative_mass && mass < 0 && p2 - mass * mass < 0) {
        mass = -mag;
    }
    const double mass2 = mass >= 0 ? mass * mass : -mass * mass;
    const double e2 = p2 + mass2;
    return std::sqrt(e2 > 0 ? e2 : 0.);
}

double cartesian_eta(double pt, double z)
{
    if (pt > 0) {
        static const double large_z = std::pow(std::numeric_limits<double>::epsilon(), -.25);
        const double ratio = z / pt;
        if (std::fabs(ratio) < large_z) {
            return std::log(ratio + std::sqrt(ratio * ratio + 1.0));
        }
        return z > 0 ? std::log(2.0 * ratio + 0.5 / ratio) : -std::log(-2.0 * ratio);
    }
    if (z == 0) {
        return 0;
    }
    return z > 0 ? z + eta_max : z - eta_max;
}

double alternate_energy(const ConstituentVector& v)
{
    const double x = v.px;
    const double y = v.py;
    const double z = v.pz;
    const double e = v.energy;
    const double pt = std::sqrt(x * x + y * y);
    const double mass2 = e * e - x * x - y * y - z * z;
    return jet_energy({pt, cartesian_eta(pt, z), 0.,
                       mass2 >= 0 ? std::sqrt(mass2) : -std::sqrt(-mass2)}, false);
}

void add_particle(Term& term, const Particle& particle)
{
    const PtEtaPhiM& p = particle.p4;
    const bool is_jet = particle.kind == ParticleKind::jet;
    // Jet_v1::p4 uses PtEtaPhiM4D. TLorentzVector::SetPtEtaPhiM uses abs(pt).
    const double pt = is_jet ? p.pt : std::abs(p.pt);
    const double phi = is_jet ? normalized_phi(p.phi) : p.phi;
    add_particle(term, static_cast<float>(pt * std::cos(phi)),
                 static_cast<float>(pt * std::sin(phi)), static_cast<float>(p.pt));
}

struct TermDefinition {
    TermType type;
    std::string_view name;
    Mask source;
};

constexpr TermDefinition electron_term{TermType::electron, "RefEle", refined | 0x1};
constexpr TermDefinition photon_term{TermType::photon, "RefGamma", refined | 0x2};
constexpr TermDefinition tau_term{TermType::tau, "RefTau", refined | 0x4};
constexpr TermDefinition muon_term{TermType::muon, "Muons", refined | 0x8};
constexpr TermDefinition muon_eloss_term{TermType::muon_eloss, "MuonEloss", 0x40008};
constexpr TermDefinition jet_term{TermType::jet, "RefJet", refined | 0x10};

enum class MuonJet { ordinary, removed, fsr };

class Evaluation {
public:
    Evaluation(const Config& config, const EventInput& event)
        : config_(config), event_(event)
    {
        result_.terms.reserve(7);
    }

    Result run()
    {
        try {
            if (!event_.core_container_available) {
                throw Stop{Status::failure, Issue::missing_core_container};
            }
            if (!event_.association_map_available) {
                throw Stop{Status::failure, Issue::missing_association_map};
            }
            result_.flags.resize(event_.associations.size());
            build_hard_term(electron_term, event_.electrons);
            build_hard_term(photon_term, event_.photons);
            build_hard_term(tau_term, event_.taus);
            build_hard_term(muon_term, event_.muons);
            build_jet_terms();
        } catch (const Stop& stop) {
            result_.status = stop.status;
            note(stop.issue, stop.severity, stop.detail);
        }
        return std::move(result_);
    }

private:
    void note(Issue issue, Severity severity, std::string_view detail = {})
    {
        result_.diagnostics.push_back({issue, severity, stage_, object_, detail});
    }

    Term& add_term(TermDefinition definition)
    {
        return result_.terms.emplace_back(Term{definition.type, definition.name,
                                               definition.source, {}, {}});
    }

    const Particle& get_particle(ObjectId id)
    {
        const auto it = std::ranges::find(event_.input_particles, id, &Particle::id);
        if (it == event_.input_particles.end()) {
            throw Stop{Status::source_access_exception, Issue::unbound_valid_target};
        }
        return *it;
    }

    const Particle& original(const Particle& current, bool use_link)
    {
        if (!use_link) {
            return current;
        }
        if (!current.original_link || current.original_link->state != LinkState::valid) {
            throw Stop{Status::source_access_exception, Issue::invalid_original_link};
        }
        return get_particle(current.original_link->id);
    }

    void check_owner(const ViewInfo& view)
    {
        if (view.owns_elements && !view.owning_key_available) {
            throw Stop{Status::failure, Issue::missing_owner_key};
        }
    }

    bool slot_selected(std::size_t association, std::size_t slot) const
    {
        return slot < 64 && (result_.flags[association] & (Mask{1} << slot));
    }

    bool selected_anywhere(ObjectId id) const
    {
        for (std::size_t i = 0; i < event_.associations.size(); ++i) {
            if (slot_selected(i, find_slot(event_.associations[i], id))) {
                return true;
            }
        }
        return false;
    }

    std::span<const OverlapEdge> overlaps(const Association& row, std::size_t slot) const
    {
        if (slot >= row.overlap_rows.size()) {
            throw Stop{Status::source_access_exception, Issue::malformed_overlap_rows};
        }
        return row.overlap_rows[slot];
    }

    void build_hard_term(TermDefinition definition, const HardObjectView& view)
    {
        stage_ = definition.type;
        object_.reset();
        if (!view.view.available) {
            throw Stop{Status::failure, Issue::missing_hard_container};
        }
        Term& term = add_term(definition);
        const bool is_muon_term = definition.type == TermType::muon;
        if (is_muon_term) {
            eloss_term_ = result_.terms.size();
            add_term(muon_eloss_term);
        }
        if (event_.associations.empty()) {
            note(Issue::empty_association_map, Severity::warning);
            return;
        }
        result_.constituent_links_available = true;
        result_.constituent_weights_available = true;
        check_owner(view.view);
        if (view.selected.empty()) {
            return;
        }
        const bool use_link = view.selected.front().original_link.has_value();
        if (view.view.has_shallow_store && !use_link) {
            note(Issue::shallow_missing_original, Severity::warning);
            return;
        }
        for (const Particle& row : view.selected) {
            object_ = row.id;
            add_hard_particle(term, row, original(row, use_link), is_muon_term);
        }
    }

    void add_hard_particle(Term& term, const Particle& row, const Particle& orig,
                           bool is_muon_term)
    {
        using Membership = std::pair<std::size_t, std::size_t>;
        std::vector<Membership> memberships;
        for (std::size_t i = 0; i < event_.associations.size(); ++i) {
            const auto slot = find_slot(event_.associations[i], orig.id);
            if (slot != no_slot) {
                memberships.emplace_back(i, slot);
            }
        }
        if (memberships.empty()) {
            if (orig.kind == ParticleKind::electron) {
                const Severity severity = orig.p4.pt > Config::missing_object_warning_pt
                    ? Severity::warning : Severity::debug;
                note(Issue::missing_electron_association, severity);
                term.contributors.push_back({row.id, 0.F});
                return;
            }
            note(Issue::missing_other_association, Severity::error);
        }
        for (const auto& [i, slot] : memberships) {
            if (slot_selected(i, slot)) {
                return;
            }
        }
        bool has_overlap = false;
        for (const auto& [i, slot] : memberships) {
            for (const OverlapEdge& edge : overlaps(event_.associations[i], slot)) {
                if (slot_selected(i, edge.slot) && (edge.types & pflow_overlap)) {
                    has_overlap = true;
                    // Stock hasOverlaps stops at the first qualifying edge per association.
                    break;
                }
            }
        }
        if (!has_overlap) {
            for (const auto& [i, slot] : memberships) {
                result_.flags[i] |= stock_signed_mask(slot);
            }
        } else if (!is_muon_term) {
            return;
        }
        if (is_muon_term && orig.kind == ParticleKind::muon
            && require_aux(orig.muon_type, "muonType") == calo_tagged_muon) {
            bool veto = false;
            for (const auto& [i, slot] : memberships) {
                const Association& association = event_.associations[i];
                for (const OverlapEdge& edge : overlaps(association, slot)) {
                    if (edge.slot >= association.objects.size()) {
                        throw Stop{Status::stock_undefined, Issue::malformed_overlap_rows};
                    }
                    const ObjectLink& link = association.objects[edge.slot];
                    if (link.state == LinkState::valid && slot_selected(i, edge.slot)
                        && get_particle(link.id).kind == ParticleKind::electron) {
                        veto = true;
                        break;
                    }
                }
            }
            if (veto) {
                return;
            }
        }
        add_particle(term, row);
        term.contributors.push_back({row.id, 1.F});
    }

    const CoreTerm& soft_core()
    {
        const auto it = std::ranges::find(event_.core_terms, Config::core_key, &CoreTerm::name);
        if (it == event_.core_terms.end()) {
            throw Stop{Status::failure, Issue::missing_soft_core, {}, Severity::warning};
        }
        if ((it->source & track_source) != track_source) {
            throw Stop{Status::unsupported_input, Issue::non_track_core,
                       "the resolved recipe requires a PVSoftTrkCore track source"};
        }
        return *it;
    }

    std::size_t jet_association(ObjectId id) const
    {
        for (std::size_t i = 0; i < event_.associations.size(); ++i) {
            const Association& row = event_.associations[i];
            // Stock MissingETAssociationMap_v1::f_findConst stops at the first misc.
            if (row.is_misc) {
                break;
            }
            if (same_target(row.jet, id)) {
                return i;
            }
        }
        return no_slot;
    }

    std::optional<std::size_t> primary_vertex_index() const
    {
        const auto it = std::ranges::find(event_.vertex_types, primary_vertex);
        if (it == event_.vertex_types.end()) {
            return std::nullopt;
        }
        return static_cast<std::size_t>(it - event_.vertex_types.begin());
    }

    template <class T>
    float pv_moment(const std::optional<std::span<const T>>& rows,
                    std::optional<std::size_t> pv, std::string_view name) const
    {
        if (!rows || !pv) {
            return 0.F;
        }
        if (*pv >= rows->size()) {
            throw Stop{Status::stock_undefined, Issue::short_vertex_moment, name};
        }
        return static_cast<float>((*rows)[*pv]);
    }

    bool pass_jvt(const Jet& jet) const
    {
        const JvtRange& range = config_.jvt_range();
        const double pt = jet.particle.p4.pt;
        if (pt < range.min_pt_mev || pt > range.max_pt_mev) {
            return true;
        }
        const float eta = require_aux(jet.jvt_eta, "JetEtaName");
        if (!(std::abs(eta) >= range.min_abs_eta && std::abs(eta) <= range.max_abs_eta)) {
            return true;
        }
        const float score = require_aux(jet.nnjvt, "NNJvt");
        return score > config_.jvt_cut(static_cast<float>(pt), eta);
    }

    bool electron_exception(const Association& row, std::size_t i,
                            std::span<const Particle* const> electrons) const
    {
        const ConstituentVector trk = reduce(row.trk, result_.flags[i]);
        const float unique_pt = row.jet_tracks.sum_pt - trk.sum_pt;
        ConstituentVector cal_sum{};
        ConstituentVector trk_sum{};
        for (const Particle* electron : electrons) {
            const std::size_t slot = find_slot(row, electron->id);
            const Mask mask = slot == no_slot ? 0 : stock_signed_mask(slot);
            add(cal_sum, reduce(row.cal, mask));
            add(trk_sum, reduce(row.trk, mask));
        }
        const float electron_trk_pt = std::sqrt(trk_sum.px * trk_sum.px + trk_sum.py * trk_sum.py);
        return electron_trk_pt > 1e-9 && unique_pt > Config::min_electron_residual_track_pt;
    }

    MuonJet classify_muon_jet(const Jet& jet, const Particle& muon,
                             std::optional<std::size_t> pv) const
    {
        const float id_pt = muon.id_track_pt.value_or(0.F);
        float track_sum = pv_moment(jet.sum_pt_trk_pt500, pv, "SumPtTrkPt500");
        if (0.9999 * id_pt > track_sum) {
            track_sum += id_pt;
        }
        const float track_count = pv_moment(jet.num_trk_pt500, pv, "NumTrkPt500");
        if (id_pt > 1e-9 && track_sum > 1e-9
            && jet.particle.p4.pt / id_pt < Config::muon_jet_id_pt_ratio
            && id_pt / track_sum > Config::muon_fraction
            && track_count < Config::muon_jet_track_count) {
            return MuonJet::removed;
        }
        float ps_energy{};
        if (jet.ps_frac) {
            ps_energy = *jet.ps_frac;
        } else if (jet.energy_per_sampling) {
            if (jet.energy_per_sampling->size() <= 4) {
                throw Stop{Status::stock_undefined, Issue::short_sampling_vector};
            }
            ps_energy = (*jet.energy_per_sampling)[0] + (*jet.energy_per_sampling)[4];
        } else {
            throw Stop{Status::failure, Issue::missing_presampler_input};
        }
        if (track_sum > 1e-9 && track_count < Config::em_jet_track_count
            && id_pt / track_sum > Config::muon_fraction
            && require_aux(jet.em_frac, "EMFrac") > Config::em_jet_fraction
            && require_aux(jet.width, "Width") < Config::em_jet_width
            && ps_energy > Config::em_jet_ps_energy) {
            return MuonJet::fsr;
        }
        return MuonJet::ordinary;
    }

    void add_jet(Term& hard, Term& soft, const Jet& jet, std::size_t i,
                 std::optional<std::size_t> pv)
    {
        const Association& row = event_.associations[i];
        const auto& p4 = jet.particle.p4;
        const double abs_eta = std::abs(p4.eta);
        const bool passes_pt_cuts =
            (abs_eta < Config::forward_jet_eta && p4.pt > Config::central_jet_pt)
            || (abs_eta >= Config::forward_jet_eta && p4.pt > Config::forward_jet_pt);

        // The active JVT assignment overwrites the initial eta rejection.
        const bool jvt_pass = pass_jvt(jet);
        bool rejected = !jvt_pass;
        ConstituentVector cal = reduce(row.cal, result_.flags[i]);
        const bool cal_overlap = cal.energy > 0;
        const PtEtaPhiM& constituent = require_aux(jet.constituent_p4, "JetConstitScaleMomentum");
        const double constituent_energy = jet_energy(constituent, true);
        const std::size_t alternate = find_slot(row, row.jet.id);
        const double denominator = alternate == no_slot ? constituent_energy
            : alternate_energy(reduce(row.cal, stock_unsigned_mask(alternate)));
        const double const_sf = denominator > 1e-9 ? constituent_energy / denominator : 0.;
        scale(cal, static_cast<float>(const_sf));
        const double cphi = normalized_phi(constituent.phi);
        double overlap_px = constituent.pt * std::cos(cphi) - cal.px;
        double overlap_py = constituent.pt * std::sin(cphi) - cal.py;

        if (!eloss_term_) {
            throw Stop{Status::failure, Issue::missing_muon_eloss_term, {}, Severity::warning};
        }
        std::vector<const Particle*> muons;
        std::vector<const Particle*> electrons;
        for (const ObjectLink& link : row.objects) {
            if (link.state != LinkState::valid) {
                continue;
            }
            const Particle& linked = get_particle(link.id);
            if (linked.kind == ParticleKind::muon) {
                const Particle& orig = original(linked, linked.original_link.has_value());
                if (selected_anywhere(orig.id)) {
                    muons.push_back(&orig);
                }
            } else if (linked.kind == ParticleKind::electron) {
                const Particle& orig = original(linked, linked.original_link.has_value());
                if (slot_selected(i, find_slot(row, orig.id))
                    && orig.p4.pt > Config::min_electron_pt) {
                    electrons.push_back(&orig);
                }
            }
        }
        const bool electron_track_exception = electron_exception(row, i, electrons);

        bool muon_fsr = false;
        Mask muon_mask{};
        for (const Particle* muon : muons) {
            const float loss = require_aux(muon->energy_loss, "EnergyLoss");
            if (!rejected) {
                const MuonJet classification = classify_muon_jet(jet, *muon, pv);
                rejected = classification == MuonJet::removed;
                if (classification == MuonJet::fsr) {
                    const float factor = loss >= constituent_energy
                        ? 0.F : static_cast<float>(1. - loss / constituent_energy);
                    overlap_px *= factor;
                    overlap_py *= factor;
                    muon_fsr = true;
                }
            }
            // All five xAOD::Muon::EnergyLossType enumerators enter this mask.
            if (require_aux(muon->energy_loss_type, "energyLossType") <= 4) {
                muon_mask |= stock_signed_mask(find_slot(row, muon->id));
            }
        }

        const ConstituentVector muon_cal = reduce(row.cal, muon_mask);
        overlap_px += muon_cal.px;
        overlap_py += muon_cal.py;
        const double overlap_pt = std::sqrt(overlap_px * overlap_px + overlap_py * overlap_py);
        const double unique_fraction = 1. - (cal.energy - muon_cal.energy) / constituent_energy;

        JetDecision decision{jet.particle.id, JetDisposition::outside_soft_region,
                             jvt_pass, passes_pt_cuts, abs_eta > Config::max_jet_eta, rejected};
        if (muon_fsr) {
            add_particle(result_.terms[*eloss_term_], static_cast<float>(overlap_px),
                         static_cast<float>(overlap_py), static_cast<float>(overlap_pt));
            decision.disposition = JetDisposition::muon_eloss;
        } else if (passes_pt_cuts && !rejected
                   && (!cal_overlap || ((unique_fraction > Config::min_unique_energy_fraction
                                         || electron_track_exception)
                                        && overlap_pt > Config::min_weighted_jet_pt))) {
            if (!cal_overlap) {
                add_particle(hard, jet.particle);
            } else {
                const double phi = normalized_phi(p4.phi);
                // Jet_v1::px()/py() return float before multiplying by the double weight.
                const auto px = static_cast<float>(p4.pt * std::cos(phi));
                const auto py = static_cast<float>(p4.pt * std::sin(phi));
                add_particle(hard, static_cast<float>(unique_fraction * px),
                             static_cast<float>(unique_fraction * py),
                             static_cast<float>(unique_fraction * p4.pt));
            }
            hard.contributors.push_back({jet.particle.id, static_cast<float>(unique_fraction)});
            decision.disposition = JetDisposition::hard;
        } else if (abs_eta < 2.5 || !(soft.source & central)) {
            ConstituentVector residual{};
            if (row.jet_tracks.energy > 1e-9) {
                residual = row.jet_tracks;
                subtract(residual, reduce(row.trk, result_.flags[i]));
            }
            add_particle(soft, residual.px, residual.py, residual.sum_pt);
            soft.contributors.push_back({jet.particle.id, std::nullopt});
            decision.disposition = JetDisposition::soft_track;
        }
        result_.jet_decisions.push_back(decision);
    }

    void build_jet_terms()
    {
        stage_ = TermType::jet;
        object_.reset();
        if (!event_.jet_view.available) {
            throw Stop{Status::failure, Issue::missing_jet_container};
        }
        Term& hard = add_term(jet_term);
        const CoreTerm& core = soft_core();
        Term& soft = add_term({TermType::soft_track, "PVSoftTrk", core.source});
        if (event_.associations.empty()) {
            note(Issue::empty_association_map, Severity::warning);
            return;
        }
        result_.constituent_links_available = true;
        result_.constituent_weights_available = true;
        add_met(soft, core.momentum);
        check_owner(event_.jet_view);
        const bool use_link = !event_.jets.empty()
            && event_.jets.front().particle.original_link.has_value();
        const auto pv = primary_vertex_index();
        for (const Jet& jet : event_.jets) {
            object_ = jet.particle.id;
            const auto association = jet_association(original(jet.particle, use_link).id);
            if (association == no_slot) {
                note(Issue::missing_jet_association, Severity::warning);
                result_.jet_decisions.push_back({jet.particle.id, JetDisposition::missing_association,
                                                 std::nullopt, false, false, false});
                continue;
            }
            add_jet(hard, soft, jet, association, pv);
        }
        object_.reset();
        for (std::size_t i = event_.associations.size(); i > 0; --i) {
            const Association& row = event_.associations[i - 1];
            if (row.is_misc) {
                const auto vec = reduce(row.trk, result_.flags[i - 1], true);
                add_particle(soft, vec.px, vec.py, vec.sum_pt);
                break;
            }
        }
    }

    const Config& config_;
    const EventInput& event_;
    Result result_;
    std::optional<TermType> stage_;
    std::optional<ObjectId> object_;
    std::optional<std::size_t> eloss_term_;
};

bool valid_edges(std::span<const float> edges)
{
    return edges.size() >= 2
        && std::ranges::all_of(edges, [](float x) { return std::isfinite(x); })
        && std::adjacent_find(edges.begin(), edges.end(), std::greater_equal<float>{}) == edges.end();
}

}  // namespace

Config::Config(NnjvtCutMap cut_map, JvtRange jvt_range)
    : cuts_(std::move(cut_map)), range_(jvt_range)
{
    if (!valid_edges(cuts_.pt_edges_mev) || !valid_edges(cuts_.eta_edges)) {
        throw std::invalid_argument("MET NNJVT cut map requires finite increasing bin edges");
    }
    const std::size_t pts = cuts_.pt_edges_mev.size() - 1;
    const std::size_t etas = cuts_.eta_edges.size() - 1;
    if (cuts_.cuts.size() / etas != pts || cuts_.cuts.size() % etas != 0
        || !std::ranges::all_of(cuts_.cuts, [](float cut) { return std::isfinite(cut); })) {
        throw std::invalid_argument("MET NNJVT cut map does not match the declared intervals");
    }
    if (!(range_.min_pt_mev <= range_.max_pt_mev
          && range_.min_abs_eta <= range_.max_abs_eta)) {
        throw std::invalid_argument("MET NNJVT selection range is invalid");
    }
}

float Config::jvt_cut(float pt, float eta) const
{
    const auto& pts = cuts_.pt_edges_mev;
    const auto& etas = cuts_.eta_edges;
    const auto p = std::lower_bound(pts.begin(), pts.end(), pt);
    const auto e = std::lower_bound(etas.begin(), etas.end(), eta);
    if (p == pts.end() || e == etas.begin() || e == etas.end()) {
        return -1.F;
    }
    const auto pt_bin = static_cast<std::size_t>(std::max(p, pts.begin() + 1) - pts.begin() - 1);
    const auto eta_bin = static_cast<std::size_t>(e - etas.begin() - 1);
    return cuts_.cuts[pt_bin * (etas.size() - 1) + eta_bin];
}

Result evaluate(const Config& config, const EventInput& event)
{
    return Evaluation(config, event).run();
}

}  // namespace atlas::met::maker
