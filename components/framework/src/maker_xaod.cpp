// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#include "maker_xaod.hpp"

#include <AthContainers/AuxElement.h>
#include <xAODCore/ShallowAuxContainer.h>
#include <xAODMissingET/MissingETAssociation.h>
#include <xAODMuon/Muon.h>
#include <xAODTracking/TrackParticle.h>

#include <limits>
#include <stdexcept>
#include <string>

namespace atlas::met::framework {
namespace {
const SG::ConstAccessor<ParticleLink> original_link{"originalObjectLink"};
const SG::ConstAccessor<float> energy_loss{"EnergyLoss"};
const SG::ConstAccessor<std::uint8_t> energy_loss_type{"energyLossType"};
const SG::ConstAccessor<std::uint16_t> muon_type{"muonType"};
const SG::ConstAccessor<ElementLink<xAOD::TrackParticleContainer>> id_track{"inDetTrackParticleLink"};
const SG::ConstAccessor<float> jvt{"NNJvt"};
const SG::ConstAccessor<float> jvt_eta{"eta"};
const SG::ConstAccessor<std::vector<int>> trk_n{"NumTrkPt500"};
const SG::ConstAccessor<std::vector<float>> trk_sum{"SumPtTrkPt500"};
const SG::ConstAccessor<float> em{"EMFrac"};
const SG::ConstAccessor<float> width{"Width"};
const SG::ConstAccessor<float> ps{"PSFrac"};
const SG::ConstAccessor<std::vector<float>> sampling{"EnergyPerSampling"};
const SG::ConstAccessor<float> const_pt{"JetConstitScaleMomentum_pt"};
const SG::ConstAccessor<float> const_eta{"JetConstitScaleMomentum_eta"};
const SG::ConstAccessor<float> const_phi{"JetConstitScaleMomentum_phi"};
const SG::ConstAccessor<float> const_mass{"JetConstitScaleMomentum_m"};

maker::ParticleKind kind(xAOD::Type::ObjectType type)
{
    using K = maker::ParticleKind;
    switch (type) {
    case xAOD::Type::Jet: return K::jet;
    case xAOD::Type::Electron: return K::electron;
    case xAOD::Type::Photon: return K::photon;
    case xAOD::Type::Tau: return K::tau;
    case xAOD::Type::Muon: return K::muon;
    default: return K::other;
    }
}

template <class T>
std::optional<T> optional_value(const SG::ConstAccessor<T>& accessor, const SG::AuxElement& p)
{
    if (!accessor.isAvailable(p)) return std::nullopt;
    return accessor(p);
}

template <class T>
std::optional<std::span<const T>> optional_view(const SG::ConstAccessor<std::vector<T>>& accessor,
                                                const SG::AuxElement& p)
{
    if (!accessor.isAvailable(p)) return std::nullopt;
    return std::span<const T>(accessor(p));
}

maker::ConstituentVector constituent(const xAOD::MissingETAssociation::ConstVec& v)
{
    return {v.cpx(), v.cpy(), v.cpz(), v.ce(), v.sumpt()};
}
} // namespace

maker::ObjectId MakerInputs::key(const xAOD::IParticle& p)
{
    if (!p.container() || p.index() > std::numeric_limits<std::uint32_t>::max()
        || owners_.size() >= std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("Invalid xAOD MET owning collection or row");
    }
    const auto owner = owners_.try_emplace(p.container(),
        maker::OwnerKey{static_cast<std::uint32_t>(owners_.size() + 1)}).first->second;
    const auto index = static_cast<std::uint32_t>(p.index());
    const auto [it, inserted] = current_.try_emplace({owner.value, index}, &p);
    if (!inserted && it->second != &p) throw std::runtime_error("Ambiguous xAOD MET owner/index");
    return {owner, index};
}

template <class Link>
maker::ObjectLink MakerInputs::resolve(const Link& link)
{
    if (link.isValid()) {
        const auto* p = *link;
        if (p) {
            const auto id = key(*p);
            remember(*p);
            return {maker::LinkState::valid, id};
        }
    }
    const bool null = (link.key() == 0 && link.index() == 0)
        || link.index() == std::numeric_limits<std::uint32_t>::max();
    return {null ? maker::LinkState::null : maker::LinkState::unresolved, {}};
}

void MakerInputs::remember(const xAOD::IParticle& p)
{
    if (seen_.insert(&p).second) particles_.push_back(snapshot(p));
}

maker::Particle MakerInputs::snapshot(const xAOD::IParticle& p)
{
    maker::Particle value{};
    value.id = key(p);
    value.kind = kind(p.type());
    value.p4 = {p.pt(), p.eta(), p.phi(), p.m()};
    if (original_link.isAvailable(p)) value.original_link = resolve(original_link(p));
    if (value.kind == maker::ParticleKind::muon) {
        value.muon_type = optional_value(muon_type, p);
        value.energy_loss = optional_value(energy_loss, p);
        value.energy_loss_type = optional_value(energy_loss_type, p);
        if (id_track.isAvailable(p) && id_track(p).isValid()) {
            if (const auto* track = *id_track(p)) value.id_track_pt = track->pt();
        }
    }
    return value;
}

void MakerInputs::selected(Hard hard, const xAOD::IParticle& p)
{
    auto& selected = hard_.at(hard);
    if (selected.empty()) {
        hard_info_.at(hard).has_shallow_store =
            dynamic_cast<const xAOD::ShallowAuxContainer*>(p.container()->getConstStore()) != nullptr;
    }
    selected.push_back(snapshot(p));
}

void MakerInputs::jets(const xAOD::JetContainer& jets, bool owning_key_available)
{
    jet_info_.owns_elements = jets.ownPolicy() == SG::OWN_ELEMENTS;
    jet_info_.owning_key_available = owning_key_available;
    jets_.reserve(jets.size());
    for (const auto* p : jets) {
        maker::Jet jet{};
        jet.particle = snapshot(*p);
        if (const_pt.isAvailable(*p) && const_eta.isAvailable(*p)
            && const_phi.isAvailable(*p) && const_mass.isAvailable(*p)) {
            jet.constituent_p4 = maker::PtEtaPhiM{const_pt(*p), const_eta(*p), const_phi(*p), const_mass(*p)};
        }
        jet.nnjvt = optional_value(jvt, *p);
        jet.jvt_eta = optional_value(jvt_eta, *p);
        jet.num_trk_pt500 = optional_view(trk_n, *p);
        jet.sum_pt_trk_pt500 = optional_view(trk_sum, *p);
        jet.em_frac = optional_value(em, *p);
        jet.width = optional_value(width, *p);
        jet.ps_frac = optional_value(ps, *p);
        jet.energy_per_sampling = optional_view(sampling, *p);
        jets_.push_back(jet);
    }
}

void MakerInputs::association_map(const xAOD::MissingETAssociationMap& map)
{
    association_storage_.resize(map.size());
    associations_.reserve(map.size());
    for (std::size_t i = 0; i < map.size(); ++i) {
        const auto& association = *map.at(i);
        auto& data = association_storage_[i];
        for (const auto& link : association.objectLinks()) data.objects.push_back(resolve(link));
        data.edges.resize(data.objects.size());
        const auto& indices = association.overlapIndices();
        const auto& types = association.overlapTypes();
        if (indices.size() != types.size() || indices.size() != data.objects.size()) {
            throw std::runtime_error("Misaligned MET association slot arrays");
        }
        for (std::size_t slot = 0; slot < indices.size(); ++slot) {
            if (indices[slot].size() != types[slot].size()) {
                throw std::runtime_error("Misaligned MET association overlap rows");
            }
            auto& row = data.edges[slot];
            row.reserve(indices[slot].size());
            for (std::size_t edge = 0; edge < indices[slot].size(); ++edge) {
                row.push_back({indices[slot][edge], types[slot][edge]});
            }
        }
        for (const auto& row : data.edges) data.overlap_rows.emplace_back(row);
        for (std::size_t key = 0; key < association.sizeCal(); ++key) {
            data.cal.push_back({association.calkey(key), constituent(association.calVec(key))});
        }
        for (std::size_t key = 0; key < association.sizeTrk(); ++key) {
            data.trk.push_back({association.trkkey(key), constituent(association.trkVec(key))});
        }
        const auto jet_tracks = association.isMisc()
            ? maker::ConstituentVector{} : constituent(association.jetTrkVec());
        associations_.push_back({resolve(association.jetLink()), association.isMisc(),
            data.objects, data.overlap_rows, data.cal, data.trk, jet_tracks});
    }
}

void MakerInputs::core(const xAOD::MissingETContainer& metcore)
{
    core_.reserve(metcore.size());
    for (const auto* row : metcore) {
        if (row) core_.push_back({row->name(), row->source(), {row->mpx(), row->mpy(), row->sumet()}});
    }
}

void MakerInputs::vertices(const xAOD::VertexContainer* pvs)
{
    if (!pvs) return;
    vertices_.reserve(pvs->size());
    for (const auto* pv : *pvs) vertices_.push_back(static_cast<std::uint8_t>(pv->vertexType()));
}

maker::EventInput MakerInputs::view() const
{
    const auto view_for = [this](Hard hard) -> maker::HardObjectView {
        return {hard_.at(hard), hard_info_.at(hard)};
    };
    return {particles_, view_for(electron), view_for(photon), view_for(tau), view_for(muon),
            jets_, jet_info_, associations_, core_, vertices_};
}

ParticleLink MakerInputs::current_link(maker::ObjectId id) const
{
    const auto* object = current_.at({id.owner.value, id.index});
    const auto* owning = static_cast<const xAOD::IParticleContainer*>(object->container());
    return ParticleLink(*owning, object->index());
}

void materialize(const maker::Result& result, const MakerInputs& input, xAOD::MissingETContainer& output)
{
    static const SG::Accessor<std::vector<ParticleLink>> links_out{"ConstitObjectLinks"};
    static const SG::Accessor<std::vector<float>> weights_out{"ConstitObjectWeights"};
    for (const auto& row : result.terms) {
        auto* term = new xAOD::MissingET;
        output.push_back(term);
        term->setName(std::string(row.name));
        term->setSource(row.source);
        term->setMpx(row.momentum.mpx);
        term->setMpy(row.momentum.mpy);
        term->setSumet(row.momentum.sumet);
        if (result.constituent_links_available) {
            auto& links = links_out(*term);
            links.reserve(row.contributors.size());
            for (const auto& contributor : row.contributors) links.push_back(input.current_link(contributor.current));
        }
        if (result.constituent_weights_available) {
            auto& weights = weights_out(*term);
            for (const auto& contributor : row.contributors) {
                if (contributor.weight) weights.push_back(*contributor.weight);
            }
        }
    }
}

} // namespace atlas::met::framework
