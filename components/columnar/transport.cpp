// SPDX-License-Identifier: Apache-2.0
#include "transport.h"
#include "met_adapter/maker_correction.hpp"
#include "met_adapter/significance_correction.hpp"
#include "soft.hpp"
#include "small_jvt.h"

#include <correction_plugin.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <type_traits>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace m = atlas::met::maker;
namespace a = atlas::met::adapter;
namespace sf = atlas::met::soft;
namespace sg = atlas::met_significance;

namespace {
thread_local std::string error_message;
void register_transport_nodes(correction::PluginRegistry& registry);
constexpr std::uint32_t jets_key = 452794649;
constexpr std::uint32_t electrons_key = 956497600;
constexpr std::uint32_t taus_key = 518718875;
constexpr std::uint32_t muons_key = 980095599;
constexpr std::uint32_t photons_key = 902907695;
constexpr std::uint32_t invalid = std::numeric_limits<std::uint32_t>::max();

std::uint32_t u32(const std::uint8_t* p)
{
    return (std::uint32_t(p[0]) << 24) | (std::uint32_t(p[1]) << 16)
        | (std::uint32_t(p[2]) << 8) | std::uint32_t(p[3]);
}

template<class T>
constexpr std::string_view code()
{
    if constexpr (std::is_same_v<T, float>) return "f4";
    if constexpr (std::is_same_v<T, double>) return "f8";
    if constexpr (std::is_same_v<T, std::uint8_t>) return "u1";
    if constexpr (std::is_same_v<T, std::int8_t>) return "i1";
    if constexpr (std::is_same_v<T, std::int32_t>) return "i4";
    if constexpr (std::is_same_v<T, std::uint32_t>) return "u4";
    if constexpr (std::is_same_v<T, std::uint64_t>) return "u8";
    if constexpr (std::is_same_v<T, std::int64_t>) return "i8";
}

class Columns {
public:
    explicit Columns(std::span<const MetColumn> columns)
    {
        for (const auto& column : columns) {
            if (!column.name || !column.dtype || (column.size && !column.data)
                || !fields_.try_emplace(column.name, &column).second) {
                throw std::invalid_argument("Invalid or repeated column buffer descriptor");
            }
        }
    }

    template<class T>
    std::span<const T> get(const std::string& name) const
    {
        const auto it = fields_.find(name);
        if (it == fields_.end() || std::string_view(it->second->dtype) != code<T>()) {
            throw std::invalid_argument("Missing or mistyped flat buffer " + name);
        }
        const auto& row = *it->second;
        return {static_cast<const T*>(row.data), row.size};
    }
private:
    std::unordered_map<std::string_view, const MetColumn*> fields_;
};

class Offset {
public:
    Offset(const Columns& cols, const std::string& name)
        : values_(cols.get<std::int64_t>(name))
    {
        if (values_.empty() || values_.front() != 0 || !std::ranges::is_sorted(values_)) {
            throw std::invalid_argument("Invalid jagged offsets " + name);
        }
    }
    std::size_t start(std::size_t row) const { return static_cast<std::size_t>(values_[row]); }
    std::size_t end(std::size_t row) const { return start(row + 1); }
    std::size_t size() const { return values_.size() - 1; }
    std::size_t length() const { return start(size()); }

    template<class T>
    std::span<const T> row(std::span<const T> flat, std::size_t i) const
    {
        if (i >= size() || end(i) > flat.size()) throw std::out_of_range("Jagged row outside flat buffer");
        return flat.subspan(start(i), end(i) - start(i));
    }
private:
    std::span<const std::int64_t> values_;
};

template<class T>
T at(std::span<const T> values, std::size_t index)
{
    if (index >= values.size()) throw std::out_of_range("MET flat column index out of range");
    return values[index];
}

struct P4 {
    std::span<const float> pt, eta, phi, mass;
    std::optional<double> fixed_mass;
    P4(const Columns& cols, const std::string& stem, std::optional<double> fixed_mass = {})
        : pt(cols.get<float>(stem + "pt")), eta(cols.get<float>(stem + "eta")),
          phi(cols.get<float>(stem + "phi")),
          mass(fixed_mass ? std::span<const float>{} : cols.get<float>(stem + "m")), fixed_mass(fixed_mass) {}
    m::PtEtaPhiM row(std::size_t index) const
    {
        return {at(pt,index), at(eta,index), at(phi,index), fixed_mass ? *fixed_mass : at(mass,index)};
    }
};

std::optional<double> xAOD_fixed_mass(m::ParticleKind kind)
{
    switch (kind) {
    case m::ParticleKind::electron: return 0.510998910;
    case m::ParticleKind::photon: return 0.;
    case m::ParticleKind::muon: return 105.658367;
    default: return {};
    }
}

struct ObjectRows {
    Offset events;
    P4 p4;
    std::optional<P4> prepared;
    std::uint32_t key;
    m::ParticleKind kind;
    ObjectRows(const Columns& cols, const std::string& stem, std::uint32_t key, m::ParticleKind kind,
               bool has_preparation=false)
        : events(cols,stem + ".o"), p4(cols,stem + ".",xAOD_fixed_mass(kind)), key(key), kind(kind)
    {
        if (has_preparation) prepared.emplace(cols,stem+".prep.",xAOD_fixed_mass(kind));
    }

    m::PtEtaPhiM current(std::size_t index) const
    {
        return (prepared ? *prepared : p4).row(index);
    }
};

m::ObjectLink object_link(std::uint32_t key, std::uint32_t index,
                          std::span<const ObjectRows> owners, std::size_t event)
{
    if ((key == 0 && index == 0) || index == invalid) return {m::LinkState::null, {}};
    for (const auto& owner : owners) {
        if (owner.key == key) {
            if (index >= owner.events.end(event) - owner.events.start(event)) {
                throw std::out_of_range("MET link escapes its event-local owning collection");
            }
            return {m::LinkState::valid, {{key}, index}};
        }
    }
    throw std::invalid_argument("Unrecognized MET association element-link owner key " + std::to_string(key));
}

struct Vectors {
    std::array<std::span<const float>, 5> fields;
    Vectors(const Columns& cols, const std::string& stem)
        : fields{cols.get<float>(stem+"px"), cols.get<float>(stem+"py"), cols.get<float>(stem+"pz"),
                 cols.get<float>(stem+"e"), cols.get<float>(stem+"sumpt")} {}
    m::ConstituentVector row(std::size_t i) const
    {
        return {at(fields[0],i),at(fields[1],i),at(fields[2],i),at(fields[3],i),at(fields[4],i)};
    }
};
} // namespace

struct MetColumnar {
    correction::PluginRegistry registry;
    std::unique_ptr<correction::CorrectionSet> maker, jvt;
    correction::Correction::Ref rebuild, update, significance;
    std::unique_ptr<correction::CorrectionSet> significance_set, preparation;
    std::unique_ptr<sf::SoftFile> soft;
    sg::State significance_state{};
    MetColumnar(const char* provider, const char* builder, const char* nnjvt)
    {
        registry.load_library(provider);
        register_small_jvt(registry);
        register_transport_nodes(registry);
        maker = correction::CorrectionSet::from_file(builder, registry);
        jvt = correction::CorrectionSet::from_file(nnjvt, registry);
        rebuild = maker->at(a::maker_correction_name);
        update = jvt->at("small_nnjvt_update");
    }
};

extern "C" const char* met_columnar_error() { return error_message.c_str(); }

extern "C" MetColumnar* met_columnar_open(const char* provider, const char* maker, const char* nnjvt)
{
    try { return new MetColumnar(provider, maker, nnjvt); }
    catch (const std::exception& error) { error_message = error.what(); return nullptr; }
}

extern "C" void met_columnar_close(MetColumnar* handle) { delete handle; }

extern "C" int met_decode_jet_links(const std::uint8_t* raw, std::size_t raw_size,
    const std::int64_t* offsets, std::size_t events, std::int64_t* output,
    std::uint32_t* keys, std::uint32_t* indices, std::size_t capacity)
{
    try {
        std::size_t position{};
        output[0] = 0;
        for (std::size_t i = 0; i < events; ++i) {
            if (offsets[i] < 0 || offsets[i+1] < offsets[i]
                || static_cast<std::size_t>(offsets[i+1]) > raw_size) {
                throw std::out_of_range("MET jetLink basket offset outside data");
            }
            const std::size_t bytes = offsets[i+1] - offsets[i];
            if (bytes < 16 || bytes >= 0x40000000) throw std::invalid_argument("MET jetLink basket row size");
            const auto* row = raw + offsets[i];
            const std::size_t count = u32(row + 12);
            if (!(row[4] & 0x40) || u32(row) != (0x40000000U | (bytes - 4))
                || count != (bytes - 16) / 8 || bytes != 16 + 8 * count
                || count > capacity - position) {
                throw std::invalid_argument("Unsupported memberwise MET jetLink basket");
            }
            for (std::size_t j = 0; j < count; ++j) {
                keys[position+j] = u32(row + 16 + 4*j);
                indices[position+j] = u32(row + 16 + 4*(count+j));
            }
            position += count;
            output[i+1] = static_cast<std::int64_t>(position);
        }
        return 0;
    } catch (const std::exception& error) { error_message = error.what(); return 1; }
}

namespace {
struct MetRows {
    Offset events, names;
    std::span<const std::uint8_t> text;
    std::span<const std::uint64_t> source;
    std::span<const float> mpx, mpy, sumet;
    MetRows(const Columns& cols, const std::string& stem)
        : events(cols, stem+".o"), names(cols, stem+".name_o"), text(cols.get<std::uint8_t>(stem+".name")),
          source(cols.get<std::uint64_t>(stem+".source")), mpx(cols.get<float>(stem+".mpx")),
          mpy(cols.get<float>(stem+".mpy")), sumet(cols.get<float>(stem+".sumet")) {}
    std::vector<m::CoreTerm> event(std::size_t i) const
    {
        std::vector<m::CoreTerm> rows;
        for (std::size_t row = events.start(i); row < events.end(i); ++row) {
            const auto label = names.row(text,row);
            rows.push_back({{reinterpret_cast<const char*>(label.data()), label.size()},at(source,row),
                            {at(mpx,row),at(mpy,row),at(sumet,row)}});
        }
        return rows;
    }
};

struct AssociationStorage {
    std::vector<m::ObjectLink> links;
    std::vector<m::OverlapEdge> edges;
    std::vector<std::span<const m::OverlapEdge>> overlaps;
    std::vector<m::KeyedVector> cal, trk;
};

struct Associations {
    Offset events, slot, edge, calo, track;
    std::span<const std::uint32_t> jet_key, jet_index, object_key, object_index;
    std::span<const std::int8_t> misc;
    std::span<const std::uint64_t> edge_slot, cal_key, trk_key;
    std::span<const std::uint8_t> edge_type;
    Vectors cal_vec, trk_vec, jet_trk;

    explicit Associations(const Columns& cols)
        : events(cols,"assoc.o"), slot(cols,"assoc.slot_o"), edge(cols,"assoc.edge_o"),
          calo(cols,"assoc.cal_o"), track(cols,"assoc.trk_o"),
          jet_key(cols.get<std::uint32_t>("assoc.jet_key")), jet_index(cols.get<std::uint32_t>("assoc.jet_index")),
          object_key(cols.get<std::uint32_t>("assoc.object_key")),
          object_index(cols.get<std::uint32_t>("assoc.object_index")),
          misc(cols.get<std::int8_t>("assoc.misc")), edge_slot(cols.get<std::uint64_t>("assoc.edge_slot")),
          cal_key(cols.get<std::uint64_t>("assoc.calkey")), trk_key(cols.get<std::uint64_t>("assoc.trkkey")),
          edge_type(cols.get<std::uint8_t>("assoc.edge_type")), cal_vec(cols,"assoc.cal"),
          trk_vec(cols,"assoc.trk"), jet_trk(cols,"assoc.jettrk") {}

    std::vector<m::Association> event(std::size_t i, std::span<const ObjectRows> owners,
                                      std::vector<AssociationStorage>& storage) const
    {
        storage.resize(events.end(i)-events.start(i));
        std::vector<m::Association> rows;
        rows.reserve(storage.size());
        for (std::size_t j = events.start(i); j < events.end(i); ++j) {
            auto& owned = storage[j-events.start(i)];
            const auto first_slot = slot.start(j);
            const auto last_slot = slot.end(j);
            owned.edges.reserve(edge.start(last_slot) - edge.start(first_slot));
            for (std::size_t k = first_slot; k < last_slot; ++k) {
                owned.links.push_back(object_link(at(object_key,k),at(object_index,k),owners,i));
                const auto start = owned.edges.size();
                for (std::size_t l = edge.start(k); l < edge.end(k); ++l) {
                    owned.edges.push_back({static_cast<std::size_t>(at(edge_slot,l)), at(edge_type,l)});
                }
                owned.overlaps.emplace_back(std::span(owned.edges).subspan(start));
            }
            for (std::size_t k = calo.start(j); k < calo.end(j); ++k) {
                owned.cal.push_back({at(cal_key,k),cal_vec.row(k)});
            }
            for (std::size_t k = track.start(j); k < track.end(j); ++k) {
                owned.trk.push_back({at(trk_key,k),trk_vec.row(k)});
            }
            rows.push_back({object_link(at(jet_key,j),at(jet_index,j),owners,i),
                at(misc,j)!=0,owned.links,owned.overlaps,owned.cal,owned.trk,jet_trk.row(j)});
        }
        return rows;
    }
};

struct Jets {
    P4 constituent;
    Offset vertex_events, track_moments;
    std::span<const std::uint8_t> vertex_types;
    std::span<const float> sums, jvf;

    explicit Jets(const Columns& cols)
        : constituent(cols,"jet.const."), vertex_events(cols,"vertices.o"),
          track_moments(cols,"jet.sum_pt.o"),vertex_types(cols.get<std::uint8_t>("vertices.type")),
          sums(cols.get<float>("jet.sum_pt")), jvf(cols.get<float>("jet.jvf")) {}

    std::vector<m::Jet> build(std::size_t event, const ObjectRows& rows, const MetColumnar& bundle,
                              std::span<float> scores) const
    {
        const auto vertices = vertex_events.row(vertex_types,event);
        const auto primary = std::ranges::find(vertices, std::uint8_t{1});
        const bool has_primary = primary != vertices.end();
        const std::size_t pv_index = primary - vertices.begin();
        const bool dummy = vertices.size()==1 && vertices.front()==0;
        if (!has_primary && !dummy && rows.events.end(event)!=rows.events.start(event)) {
            throw std::runtime_error("NNJVT original input has neither PriVtx nor its one dummy PV");
        }
        std::vector<m::Jet> output;
        for (std::size_t row = rows.events.start(event); row < rows.events.end(event); ++row) {
            const auto raw = rows.p4.row(row);
            const auto moment = track_moments.row(sums,row);
            const double sum_pt = has_primary ? at(moment,pv_index) : 0.;
            const std::vector<correction::Variable::Type> input{
                raw.pt, raw.eta, sum_pt, static_cast<double>(at(jvf,row)),
                std::int64_t{has_primary ? 2 : 1}, std::int64_t{has_primary ? 1 : 0}};
            std::array<double,4> value{};
            bundle.update->evaluate_array(input,value.data(),value.size());
            if (value[3]!=1.) throw std::runtime_error("Pinned NNJVT update did not write a score");
            const auto idx = static_cast<std::uint32_t>(row-rows.events.start(event));
            m::Jet jet{};
            jet.particle.id = {{rows.key},idx};
            jet.particle.kind = rows.kind;
            jet.particle.p4 = raw;
            jet.constituent_p4 = constituent.row(row);
            jet.nnjvt = static_cast<float>(value[0]);
            jet.jvt_eta = static_cast<float>(raw.eta);
            jet.sum_pt_trk_pt500 = moment;
            scores[row]=*jet.nnjvt;
            output.push_back(jet);
        }
        return output;
    }
};
} // namespace

extern "C" int met_columnar_raw_smoke(const MetColumnar* bundle, std::size_t events,
    const MetColumn* columns, std::size_t column_count, float* terms, std::size_t terms_count,
    std::uint64_t* masks, std::size_t mask_count, float* scores,
    std::int8_t* decisions, std::size_t jet_count)
{
    try {
        if (!bundle) throw std::invalid_argument("Missing C operation bundle");
        const Columns cols({columns,column_count});
        const Associations association(cols);
        const Jets jet_input(cols);
        const MetRows core(cols,"core");
        const std::array<ObjectRows,5> objects{
            ObjectRows(cols,"jet",jets_key,m::ParticleKind::jet),
            ObjectRows(cols,"electron",electrons_key,m::ParticleKind::electron),
            ObjectRows(cols,"photon",photons_key,m::ParticleKind::photon),
            ObjectRows(cols,"tau",taus_key,m::ParticleKind::tau),
            ObjectRows(cols,"muon",muons_key,m::ParticleKind::muon)};
        if (association.events.size()!=events || core.events.size()!=events
            || jet_input.vertex_events.size()!=events
            || mask_count!=association.events.length()
            || jet_count!=objects[0].events.length() || terms_count/21!=events || terms_count%21) {
            throw std::invalid_argument("MET raw smoke shape mismatch");
        }
        for (const auto& owner : objects) {
            if (owner.events.size()!=events || owner.events.length()!=owner.p4.pt.size()) {
                throw std::invalid_argument("Original MET object shape mismatch");
            }
        }
        if (jet_count) std::fill_n(decisions,jet_count,std::int8_t{-1});
        for (std::size_t event=0; event<events; ++event) {
            std::vector<m::Particle> particles;
            for (const auto& owner : objects) {
                for (auto row=owner.events.start(event); row<owner.events.end(event); ++row) {
                    const auto local = static_cast<std::uint32_t>(row-owner.events.start(event));
                    particles.push_back({{{owner.key},local},owner.kind,owner.p4.row(row)});
                }
            }
            auto jets = jet_input.build(event,objects[0],*bundle,{scores,jet_count});
            std::vector<AssociationStorage> association_storage;
            const auto map = association.event(event,objects,association_storage);
            const auto core_rows = core.event(event);
            m::EventInput input{};
            input.input_particles = particles;
            input.associations = map;
            input.jets = jets;
            input.core_terms = core_rows;
            input.vertex_types = jet_input.vertex_events.row(jet_input.vertex_types,event);
            const auto result = a::invoke_maker(*bundle->rebuild,input);
            if (result.status!=m::Status::success || result.terms.size()!=7) {
                const auto issue = result.diagnostics.empty() ? -1
                    : static_cast<int>(result.diagnostics.back().issue);
                throw std::runtime_error("Shared maker failed on raw chunk event "+std::to_string(event)
                                         +" issue "+std::to_string(issue));
            }
            for (std::size_t term=0; term<7; ++term) {
                const auto& value = result.terms[term].momentum;
                const auto index = 21*event+3*term;
                terms[index]=value.mpx;
                terms[index+1]=value.mpy;
                terms[index+2]=value.sumet;
            }
            if (!result.flags.empty()) std::copy(result.flags.begin(),result.flags.end(),masks+association.events.start(event));
            for (const auto& jet : result.jet_decisions) {
                if (jet.jvt_pass) decisions[objects[0].events.start(event)+jet.current.index]=*jet.jvt_pass;
            }
        }
        return 0;
    } catch (const std::exception& error) { error_message = error.what(); return 1; }
}

extern "C" int met_columnar_configure(MetColumnar* bundle, const char* soft,
    const char* significance, const char* preparation_provider, const char* preparation)
{
    try {
        if (!bundle || !soft || !significance || !preparation_provider || !preparation) {
            throw std::invalid_argument("Full C requires independent preparation and all MET payloads");
        }
        bundle->registry.load_library(preparation_provider);
        bundle->preparation = correction::CorrectionSet::from_file(preparation,bundle->registry);
        bundle->soft = std::make_unique<sf::SoftFile>(soft);
        bundle->significance_set = correction::CorrectionSet::from_file(significance,bundle->registry);
        bundle->significance = bundle->significance_set->at(a::significance_correction_name);
        bundle->significance_state = {};
        return 0;
    } catch (const std::exception& error) { error_message=error.what(); return 1; }
}

extern "C" void met_columnar_reset(MetColumnar* bundle)
{
    bundle->significance_state = {};
}

extern "C" int met_prep_evaluate(const MetColumnar* bundle, const char* name, std::size_t rows,
    const MetColumn* columns, std::size_t count, double* output, std::size_t components)
{
    try {
        if (!bundle || !bundle->preparation) throw std::invalid_argument("Preparation is not loaded");
        const auto op = bundle->preparation->at(name);
        if (count!=op->inputs().size() || components!=std::max(op->output_size(),std::size_t{1})) {
            throw std::invalid_argument("Preparation component or input count mismatch for "+std::string(name));
        }
        const Columns cols({columns,count});
        std::vector<correction::BatchInput> input;
        for (const auto& var : op->inputs()) {
            switch (var.type()) {
            case correction::Variable::VarType::real: {
                auto col=cols.get<double>(var.name());
                if (col.size()!=rows) throw std::invalid_argument("Preparation real row count mismatch");
                input.push_back(correction::BatchInput::reals(col.data(),rows));
                break;
            }
            case correction::Variable::VarType::integer: {
                auto col=cols.get<std::int64_t>(var.name());
                if (col.size()!=rows) throw std::invalid_argument("Preparation integer row count mismatch");
                input.push_back(correction::BatchInput::integers(col.data(),rows));
                break;
            }
            default: throw std::invalid_argument("Unsupported preparation variable type");
            }
        }
        if (!rows) return 0;
        if (op->output_size()) op->evaluate_batch_array(input,rows,output,rows*components);
        else op->evaluate_batch(input,rows,output,rows);
        return 0;
    } catch (const std::exception& error) { error_message=error.what(); return 1; }
}

namespace {
std::optional<float> float_field(const Columns& cols, const std::string& prefix, std::size_t row)
{
    if (!at(cols.get<std::uint8_t>(prefix+".valid"),row)) return {};
    return at(cols.get<float>(prefix+".value"),row);
}

std::optional<sg::TrackParameters> track(const Columns& cols, const std::string& prefix, std::size_t row)
{
    if (!at(cols.get<std::uint8_t>(prefix+".valid"),row)) return {};
    return sg::TrackParameters{at(cols.get<float>(prefix+".qop"),row),
        at(cols.get<float>(prefix+".theta"),row), at(cols.get<float>(prefix+".phi"),row)};
}

std::pair<std::size_t,std::size_t> owner_row(std::span<const ObjectRows> owners, m::ObjectId id,
                                         std::size_t event)
{
    for (std::size_t i=0; i<owners.size(); ++i) {
        if (owners[i].key!=(id.owner.value & 0x7fffffffU)) continue;
        auto r=owners[i].events.start(event)+id.index;
        if (r>=owners[i].events.end(event)) throw std::out_of_range("Current MET owner row escapes event");
        return {i,r};
    }
    throw std::invalid_argument("Unknown current MET owner");
}

constexpr std::array<const char*,5> stems{"jet","electron","photon","tau","muon"};

m::Particle raw_particle(const Columns& cols,const ObjectRows& owner, std::size_t global, std::size_t local)
{
    m::Particle value{{{owner.key},static_cast<std::uint32_t>(local)},owner.kind,owner.p4.row(global)};
    if (owner.kind==m::ParticleKind::muon) {
        value.muon_type=static_cast<std::uint16_t>(at(cols.get<std::uint32_t>("muon.type"),global));
        value.energy_loss=at(cols.get<float>("muon.loss"),global);
        value.energy_loss_type=at(cols.get<std::uint8_t>("muon.loss_type"),global);
        value.id_track_pt=float_field(cols,"muon.id_pt",global);
    }
    return value;
}

std::vector<sf::Term> soft_terms(std::span<const m::Term> terms)
{
    std::vector<sf::Term> result;
    for (const auto& term : terms) {
        const auto v=term.momentum;
        result.push_back({term.name,term.source,{v.mpx,v.mpy,v.sumet},term.contributors.size()});
    }
    return result;
}

auto prepared_kinematics(const MetColumnar& bundle,const m::PtEtaPhiM& q,m::ParticleKind kind)
{
    return bundle.preparation->at("met_prepared_p4")->evaluate_array({
        q.pt,q.eta,q.phi,q.mass,static_cast<std::int64_t>(kind)});
}

a::SignificanceLink significance_link(const Columns& cols, const ObjectRows& owner,
    std::size_t group, std::size_t row, const m::PtEtaPhiM& p4, const MetColumnar& bundle)
{
    using K=sg::ObjectKind;
    a::SignificanceLink object{};
    object.valid=true;
    object.pt_mev=p4.pt;
    object.eta=p4.eta;
    object.phi=p4.phi;
    const auto kinematics=prepared_kinematics(bundle,p4,owner.kind);
    object.p4_px_mev=kinematics.at(0);
    object.p4_py_mev=kinematics.at(1);
    switch (owner.kind) {
    case m::ParticleKind::jet: object.kind=K::jet; break;
    case m::ParticleKind::electron:
    case m::ParticleKind::photon: {
        const bool photon=owner.kind==m::ParticleKind::photon;
        object.kind=photon ? K::photon : K::electron;
        const std::string prefix=stems[group];
        const auto eta=float_field(cols,prefix+".cluster_eta",row);
        const auto calo=float_field(cols,prefix+".eta_calo_moment",row);
        object.egamma = a::EgammaInput{kinematics.at(2),eta,
            at(cols.get<std::uint8_t>(prefix+".forward_author"),row)!=0,calo,
            float_field(cols,prefix+".eta_calo_aux",row),
            photon ? at(cols.get<std::uint32_t>("photon.vertex_count"),row) : 0U};
        break;
    }
    case m::ParticleKind::tau:
        object.kind=K::tau;
        object.tau=a::TauInput{at(cols.get<std::int32_t>("tau.mode"),row),
            at(cols.get<float>("tau.tes_eta"),row),at(cols.get<float>("tau.final_eta"),row),
            at(cols.get<float>("tau.final_pt"),row)};
        break;
    case m::ParticleKind::muon:
        object.kind=K::muon;
        object.muon=a::MuonInput{static_cast<std::uint16_t>(at(cols.get<std::uint32_t>("muon.type"),row)),
            track(cols,"muon.cb",row),track(cols,"muon.id",row),track(cols,"muon.me",row)};
        break;
    default: throw std::invalid_argument("Unsupported linked MET object");
    }
    return object;
}

struct EventObjects {
    std::vector<m::Particle> originals;
    std::array<std::vector<m::Particle>,4> selected;
    std::vector<m::Jet> jets;
};

EventObjects prepare_views(const Columns& cols, std::span<const ObjectRows> owners, std::size_t event,
                           const Jets& jet_input, const MetColumnar& bundle, std::span<float> scores)
{
    EventObjects packed;
    const auto ntrk=cols.get<std::int32_t>("jet.ntrk");
    const Offset ntrk_offsets(cols,"jet.ntrk.o");
    for (std::size_t i=0; i<owners.size(); ++i) {
        const auto& owner=owners[i];
        const auto start=owner.events.start(event);
        const auto selected=i ? cols.get<std::uint8_t>(std::string(stems[i])+".selected")
                              : std::span<const std::uint8_t>{};
        for (auto row=start; row<owner.events.end(event); ++row) {
            auto original=raw_particle(cols,owner,row,row-start);
            packed.originals.push_back(original);
            if (i && at(selected,row)) {
                auto current=original;
                current.original_link=m::ObjectLink{m::LinkState::valid,original.id};
                current.id.owner.value|=0x80000000U;
                current.p4=owner.current(row);
                packed.selected[i-1].push_back(current);
            }
        }
    }
    packed.jets=jet_input.build(event,owners.front(),bundle,scores);
    const auto start=owners.front().events.start(event);
    for (std::size_t i=0; i<packed.jets.size(); ++i) {
        auto& jet=packed.jets[i];
        jet.particle.original_link=m::ObjectLink{m::LinkState::valid,jet.particle.id};
        jet.particle.id.owner.value|=0x80000000U;
        jet.num_trk_pt500=ntrk_offsets.row(ntrk,start+i);
        jet.em_frac=at(cols.get<float>("jet.em"),start+i);
        jet.width=at(cols.get<float>("jet.width"),start+i);
        jet.ps_frac=float_field(cols,"jet.ps",start+i);
        if (!jet.ps_frac) {
            const Offset samples(cols,"jet.sample.o");
            jet.energy_per_sampling=samples.row(cols.get<float>("jet.sample"),start+i);
        }
    }
    return packed;
}

void trace_prepared_rows(std::ostream& output,const Columns& cols,
    std::span<const ObjectRows> owners,const MetColumnar& bundle,std::uint64_t entry,std::size_t event)
{
    constexpr std::array<std::size_t,5> gate_count{0,7,7,1,5};
    for (std::size_t i=0; i<owners.size(); ++i) {
        const auto& owner=owners[i];
        const auto start=owner.events.start(event);
        const std::string stem=stems[i];
        const auto selected=i ? cols.get<std::uint8_t>(stem+".selected") : std::span<const std::uint8_t>{};
        const auto gates=i ? cols.get<std::uint32_t>(stem+".gate") : std::span<const std::uint32_t>{};
        for (auto row=start; row<owner.events.end(event); ++row) {
            const auto q=owner.current(row);
            const auto energy=prepared_kinematics(bundle,q,owner.kind).at(2);
            const auto local=row-start;
            const int pass=i ? (at(selected,row)!=0) : -1;
            output << entry << "\t0\tpreparedGetter\t" << stem << '\t' << local
                << '\t' << owner.key << '\t' << local << '\t' << pass
                << '\t' << std::hexfloat << q.pt << '\t' << q.eta << '\t'
                << q.phi << '\t' << q.mass << '\t' << energy << std::defaultfloat << '\n';
            output << entry << "\t0\tprepared\t" << stem << '\t' << owner.key << '\t'
                << row-start << '\t' << (i ? at(selected,row)!=0 : true) << '\t'
                << q.pt << '\t' << q.eta << '\t' << q.phi << '\t' << q.mass;
            for (std::size_t column=0; column<gate_count.at(i); ++column) {
                output << '\t' << at(gates,gate_count[i]*row+column);
            }
            output << '\n';
        }
    }
}

void trace_maker_rows(std::ostream& output,std::uint64_t entry,std::size_t state,
    const m::EventInput& input,const m::Result& result)
{
    const auto prefix=[&](std::string_view action) -> std::ostream& {
        return output << entry << '\t' << state << '\t' << action << '\t';
    };
    for (const auto& flag:result.flag_writes) {
        prefix("associationFlag") << flag.association << '\t' << flag.slot << '\t' << flag.selected << '\n';
        prefix("associationActor") << flag.association << '\t' << flag.slot << '\t'
            << static_cast<int>(flag.term) << '\t' << flag.original.owner.value << '\t'
            << flag.original.index << '\n';
    }
    for (const auto& row:result.jet_decisions) {
        if (row.jvt_pass) prefix("jetJvt") << row.current.index << '\t' << *row.jvt_pass << '\n';
    }
    prefix("rebuildMET") << "RefEle\t" << input.electrons.selected.size() << "\t1\n";
    prefix("rebuildMET") << "RefGamma\t" << input.photons.selected.size() << "\t1\n";
    prefix("rebuildMET") << "RefTau\t" << input.taus.selected.size() << "\t1\n";
    prefix("rebuildMET") << "Muons\t" << input.muons.selected.size() << "\t1\n";
    prefix("rebuildJetMET") << "RefJet\tPVSoftTrk\t" << input.jets.size() << "\t1\t1\n";
    for (const auto& term:result.terms) {
        for (const auto& c:term.contributors) {
            prefix("contributor") << term.name << '\t' << c.current.owner.value << '\t'
                << (c.current.owner.value&0x7fffffffU) << '\t' << c.current.index << '\t';
            if (c.weight) output << *c.weight;
            output << '\n';
        }
    }
}

a::SignificanceResult calculate_significance(const Columns& cols,std::span<const ObjectRows> owners,
    const m::Result& result,std::span<const sf::Term> corrected,std::size_t event,MetColumnar& bundle,
    std::ostream* trace,std::uint64_t entry,std::size_t state)
{
    std::vector<std::vector<a::SignificanceLink>> links(result.terms.size());
    std::vector<a::SignificanceTerm> input;
    input.reserve(corrected.size());
    for (std::size_t i=0; i<corrected.size(); ++i) {
        const auto& term=corrected[i];
        std::span<const a::SignificanceLink> linked;
        if (i<result.terms.size() && sg::term_action(term.name,term.source)==sg::TermAction::constituents) {
            for (const auto& contributor:result.terms[i].contributors) {
                const auto [group,row]=owner_row(owners,contributor.current,event);
                const auto& owner=owners[group];
                const auto prepared=owner.current(row);
                auto link=significance_link(cols,owner,group,row,prepared,bundle);
                if (trace && link.egamma) {
                    *trace << entry << '\t' << state << "\tegammaSignificanceGetter\t" << term.name
                           << '\t' << stems[group] << '\t' << contributor.current.index << '\t'
                           << std::hexfloat << prepared.mass << '\t' << link.egamma->energy_mev
                           << std::defaultfloat << '\n';
                }
                links[i].push_back(std::move(link));
            }
            linked=links[i];
        }
        input.push_back({term.name,term.source,{term.momentum.mpx,term.momentum.mpy,term.momentum.sumet},
                         linked,result.constituent_links_available});
    }
    const auto output=a::invoke_significance(*bundle.significance,{input,bundle.significance_state,
        at(cols.get<std::uint32_t>("event.run"),event),
        (at(cols.get<std::uint32_t>("event.type"),event)&1U)!=0,
        at(cols.get<float>("event.average_mu"),event)});
    bundle.significance_state=output.next;
    return output;
}
} // namespace

extern "C" int met_columnar_full(MetColumnar* bundle, std::size_t events, std::uint64_t first_entry,
    const char* const* states, std::size_t state_count, const MetColumn* columns,std::size_t column_count,
    float* terms_out,std::uint64_t* sources_out,float* significance_out,const char* trace_path)
{
    try {
        if (!bundle || !bundle->soft || !bundle->significance || !bundle->preparation) {
            throw std::invalid_argument("Full C configuration missing");
        }
        if (!state_count || state_count>6 || std::string_view(states[0])!="NOSYS") {
            throw std::invalid_argument("Full C requires an explicit nominal-first state list");
        }
        const Columns cols({columns,column_count});
        const Associations association(cols);
        const Jets jet_input(cols);
        const MetRows core(cols,"core"), truth(cols,"truth");
        const std::array<ObjectRows,5> objects{
            ObjectRows(cols,"jet",jets_key,m::ParticleKind::jet),
            ObjectRows(cols,"electron",electrons_key,m::ParticleKind::electron,true),
            ObjectRows(cols,"photon",photons_key,m::ParticleKind::photon,true),
            ObjectRows(cols,"tau",taus_key,m::ParticleKind::tau,true),
            ObjectRows(cols,"muon",muons_key,m::ParticleKind::muon,true)};
        if (association.events.size()!=events || core.events.size()!=events || truth.events.size()!=events
            || jet_input.vertex_events.size()!=events) {
            throw std::invalid_argument("Full C event shape mismatch");
        }
        for (const auto& owner:objects) {
            if (owner.events.size()!=events || owner.events.length()!=owner.p4.pt.size()) {
                throw std::invalid_argument("Full C object shape mismatch");
            }
        }
        std::ofstream trace;
        if (trace_path) {
            trace.open(trace_path,std::ios::app);
            if (!trace) throw std::runtime_error("Cannot append C trace");
            trace.precision(std::numeric_limits<double>::max_digits10);
        }
        std::vector<float> scores(objects.front().events.length());
        for (std::size_t event=0; event<events; ++event) {
            const auto view=prepare_views(cols,objects,event,jet_input,*bundle,scores);
            if (trace_path) trace_prepared_rows(trace,cols,objects,*bundle,first_entry+event,event);
            const auto core_rows=core.event(event);
            const auto truth_rows=truth.event(event);
            std::vector<sf::Term> truth_terms;
            for (const auto& row:truth_rows) {
                truth_terms.push_back({row.name,row.source,
                    {row.momentum.mpx,row.momentum.mpy,row.momentum.sumet},std::nullopt});
            }
            std::vector<AssociationStorage> association_storage;
            const auto map=association.event(event,objects,association_storage);
            m::EventInput input{};
            input.input_particles=view.originals;
            input.jets=view.jets;
            input.associations=map;
            input.core_terms=core_rows;
            input.electrons.selected=view.selected[0];
            input.photons.selected=view.selected[1];
            input.taus.selected=view.selected[2];
            input.muons.selected=view.selected[3];
            input.vertex_types=jet_input.vertex_events.row(jet_input.vertex_types,event);
            input.capture_flag_writes=trace_path!=nullptr;

            for (std::size_t state=0; state<state_count; ++state) {
                const std::string stamp="entry "+std::to_string(first_entry+event)+" state "+states[state];
                const auto maker=a::invoke_maker(*bundle->rebuild,input);
                if (maker.status!=m::Status::success || maker.terms.size()!=7) {
                    const auto issue=maker.diagnostics.empty() ? -1:static_cast<int>(maker.diagnostics.back().issue);
                    throw std::runtime_error("C maker issue "+std::to_string(issue)+" "+stamp);
                }
                if (trace_path) trace_maker_rows(trace,first_entry+event,state,input,maker);
                auto rebuilt=soft_terms(maker.terms);
                if (!maker.constituent_links_available) {
                    for (auto& term:rebuilt) term.constituent_count.reset();
                }
                const auto type=at(cols.get<std::uint32_t>("event.type"),event);
                const sf::EventView soft_event{rebuilt,truth_terms,
                    at(cols.get<float>("event.actual_mu"),event),true,true,(type&1U)!=0,true};
                const auto corrected=bundle->soft->apply(soft_event,states[state]);
                if (corrected.status!=sf::Status::success || !corrected.term_index) {
                    throw std::runtime_error("C soft status "+std::to_string(static_cast<int>(corrected.status))
                                             +" "+stamp);
                }
                rebuilt.at(*corrected.term_index).momentum=corrected.momentum;
                const auto final=bundle->soft->build_final(rebuilt);
                if (final.status!=sf::Status::success) throw std::runtime_error("C Final failed "+stamp);
                rebuilt.push_back({sf::final_name,sf::final_source,final.momentum,std::nullopt});
                const auto significance=calculate_significance(cols,objects,maker,rebuilt,event,*bundle,
                    trace_path ? &trace : nullptr,first_entry+event,state);
                if (significance.reduction.status!=sg::Status::success) {
                    throw std::runtime_error("C significance status "+
                        std::to_string(static_cast<int>(significance.reduction.status))+" "+stamp);
                }
                const auto stride=event*state_count+state;
                for (std::size_t i=0; i<rebuilt.size(); ++i) {
                    const auto& term=rebuilt[i];
                    const auto q=term.momentum;
                    const auto mp=term.name==sf::final_name ? final.decorations : bundle->soft->decorate(q);
                    const auto index=8*stride+i;
                    sources_out[index]=term.source;
                    const std::array fields{q.mpx,q.mpy,q.sumet,mp.met,mp.phi};
                    std::copy(fields.begin(),fields.end(),terms_out+5*index);
                }
                const auto d=significance.reduction.decorations;
                const std::array values{d.significance,d.directional,d.met_over_sqrt_sumet,d.met_over_sqrt_ht};
                std::copy(values.begin(),values.end(),significance_out+4*stride);
                if (trace_path) {
                    const auto prefix=[&](std::string_view action) -> std::ostream& {
                        return trace << first_entry+event << '\t' << state << '\t' << action << '\t';
                    };
                    prefix("applySystematicVariation") << (std::string_view(states[state])=="NOSYS" ? "" : states[state]) << "\t1\n";
                    prefix("setRandomSeed") << corrected.seed.setter_argument << '\n';
                    prefix("seed") << corrected.seed.truncated << '\t' << corrected.seed.effective << '\t'
                                  << (corrected.drew ? "called" : "skipped") << '\n';
                    if (corrected.drew) {
                        prefix("gaussian") << std::hexfloat << 0. << '\t' << 1. << '\t'
                            << corrected.gaussian << std::defaultfloat << '\n';
                    }
                    const auto before=maker.terms.at(*corrected.term_index).momentum;
                    const auto after=corrected.momentum;
                    prefix("applyCorrection") << "2\t" << std::hexfloat << before.mpx << '\t'
                        << before.mpy << '\t' << before.sumet << '\t' << after.mpx << '\t'
                        << after.mpy << '\t' << after.sumet << std::defaultfloat << '\n';
                    prefix("finalStatus") << "1\n";
                    prefix("varianceMET") << "1\n";
                }
            }
            if (trace_path) {
                for (const auto action:{"makerAlg","builderAlg","significanceAlg"}) {
                    trace << first_entry+event << '\t' << state_count << '\t' << action << "\t1\n";
                }
            }
        }
        if (trace_path && !trace) throw std::runtime_error("C trace write failed");
        return 0;
    } catch (const std::exception& error) { error_message=error.what(); return 1; }
}

namespace {
using Values=std::vector<correction::Variable::Type>;

class P4Access final : public correction::PluginNode {
public:
    explicit P4Access(std::array<std::size_t,5> inputs) : inputs_(inputs) {}
    std::size_t output_size() const override { return 3; }
    double evaluate(const Values&) const override { throw std::logic_error("P4 access is fixed-array"); }
    void evaluate_array(const Values& values,double* out,std::size_t size) const override
    {
        if (size!=3 || !out) throw std::invalid_argument("P4 access output size");
        const double pt=std::get<double>(values.at(inputs_[0]));
        const double eta=std::get<double>(values.at(inputs_[1]));
        double phi=std::get<double>(values.at(inputs_[2]));
        const double mass=std::get<double>(values.at(inputs_[3]));
        const auto input_kind=std::get<std::int64_t>(values.at(inputs_[4]));
        if (input_kind<static_cast<int>(m::ParticleKind::jet) ||
            input_kind>static_cast<int>(m::ParticleKind::muon)) {
            throw std::invalid_argument("Invalid prepared P4 particle kind");
        }
        const auto kind=static_cast<m::ParticleKind>(input_kind);
        const bool jet=kind==m::ParticleKind::jet;
        // JetFourMom_t normalizes phi. TLorentzVector p4 sets abs(pt).
        if (jet && (phi<=-std::numbers::pi || phi>std::numbers::pi)) {
            phi-=std::floor(phi/(2*std::numbers::pi)+.5)*2*std::numbers::pi;
        }
        const double cartesian_pt=jet ? pt : std::abs(pt);
        out[0]=cartesian_pt*std::cos(phi);
        out[1]=cartesian_pt*std::sin(phi);
        // AnalysisBase 25.2.117 Photon_v1.cxx:46 overrides e(). Jet_v1.cxx:63,
        // Egamma_v1.cxx:86, Muon_v1.cxx:62 and TauJet_v3.cxx:87 delegate to
        // ROOT 6.40.02 PtEtaPhiM4D.h:150-189 and etaMax.h for the double e().
        if (kind==m::ParticleKind::photon) {
            out[2]=pt*std::cosh(eta);
            return;
        }
        constexpr double eta_max=22756.;
        const double momentum=pt>0 ? pt*std::cosh(eta)
            : eta>eta_max ? eta-eta_max : eta < -eta_max ? -eta-eta_max : 0.;
        const double mass2=mass*mass;
        const double square=momentum*momentum+(mass>=0 ? mass2 : -mass2);
        out[2]=std::sqrt(square>0 ? square : 0.);
    }
private:
    std::array<std::size_t,5> inputs_;
};

class TrackAccess final : public correction::PluginNode {
public:
    explicit TrackAccess(std::array<std::size_t,4> inputs) : inputs_(inputs) {}
    double evaluate(const Values& values) const override
    {
        if (!std::get<std::int64_t>(values.at(inputs_[3]))) return 0.;
        const auto val=[&](std::size_t i){return static_cast<float>(std::get<double>(values.at(inputs_[i])));};
        return sg::track_pt_gev(sg::TrackParameters{val(0),val(1),val(2)})*1000.;
    }
private:
    std::array<std::size_t,4> inputs_;
};

void register_transport_nodes(correction::PluginRegistry& registry)
{
    using V=correction::Variable::VarType;
    registry.register_node("atlas.met.prepared_p4",2,[](const correction::PluginConfig& cfg,const correction::PluginContext& ctx){
        cfg.require_keys({});
        if (ctx.output_size()!=3 || ctx.child_count()) throw std::invalid_argument("Invalid p4 transport schema");
        return std::make_shared<P4Access>(std::array{
            ctx.input_index("pt",V::real),ctx.input_index("eta",V::real),ctx.input_index("phi",V::real),
            ctx.input_index("mass",V::real),ctx.input_index("particle_kind",V::integer)});
    });
    registry.register_node("atlas.met.track_pt",1,[](const correction::PluginConfig& cfg,const correction::PluginContext& ctx){
        cfg.require_keys({});
        if (ctx.output_size() || ctx.child_count()) throw std::invalid_argument("Invalid track transport schema");
        return std::make_shared<TrackAccess>(std::array{
            ctx.input_index("qop",V::real),ctx.input_index("theta",V::real),
            ctx.input_index("phi",V::real),ctx.input_index("valid",V::integer)});
    });
}
}
