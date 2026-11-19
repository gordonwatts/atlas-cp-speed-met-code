// Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
// SPDX-License-Identifier: Apache-2.0
#ifndef ATLAS_MET_MAKER_XAOD_HPP
#define ATLAS_MET_MAKER_XAOD_HPP

#include "met_maker/maker.hpp"

#include <AthLinks/ElementLink.h>
#include <xAODBase/IParticleContainer.h>
#include <xAODJet/JetContainer.h>
#include <xAODMissingET/MissingETAssociationMap.h>
#include <xAODMissingET/MissingETContainer.h>
#include <xAODTracking/VertexContainer.h>

#include <array>
#include <map>
#include <set>
#include <unordered_map>

namespace atlas::met::framework {

using ParticleLink = ElementLink<xAOD::IParticleContainer>;

class MakerInputs {
public:
    enum Hard { electron, photon, tau, muon };
    void selected(Hard kind, const xAOD::IParticle& particle);
    void jets(const xAOD::JetContainer&, bool owning_key_available);
    void association_map(const xAOD::MissingETAssociationMap&);
    void core(const xAOD::MissingETContainer&);
    void vertices(const xAOD::VertexContainer*);
    [[nodiscard]] maker::EventInput view() const;
    [[nodiscard]] ParticleLink current_link(maker::ObjectId id) const;

private:
    struct AssociationRows {
        std::vector<maker::ObjectLink> objects;
        std::vector<std::vector<maker::OverlapEdge>> edges;
        std::vector<std::span<const maker::OverlapEdge>> overlap_rows;
        std::vector<maker::KeyedVector> cal;
        std::vector<maker::KeyedVector> trk;
    };
    maker::ObjectId key(const xAOD::IParticle&);
    void remember(const xAOD::IParticle&);
    maker::Particle snapshot(const xAOD::IParticle&);
    template <class Link> maker::ObjectLink resolve(const Link&);

    std::unordered_map<const void*, maker::OwnerKey> owners_;
    std::map<std::pair<std::uint32_t, std::uint32_t>, const xAOD::IParticle*> current_;
    std::set<const xAOD::IParticle*> seen_;
    std::array<std::vector<maker::Particle>, 4> hard_;
    std::array<maker::ViewInfo, 4> hard_info_{};
    std::vector<maker::Particle> particles_;
    std::vector<maker::Jet> jets_;
    maker::ViewInfo jet_info_;
    std::vector<AssociationRows> association_storage_;
    std::vector<maker::Association> associations_;
    std::vector<maker::CoreTerm> core_;
    std::vector<std::uint8_t> vertices_;
};

void materialize(const maker::Result&, const MakerInputs&, xAOD::MissingETContainer&);

} // namespace atlas::met::framework
#endif
