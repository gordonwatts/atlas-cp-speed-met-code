// SPDX-License-Identifier: Apache-2.0
// Included only by the untimed M27 diagnostic library.
extern "C" int met_check_jet_ntrk_offsets()
{
    try {
        const auto require=[](bool condition,const char* message) {
            if (!condition) throw std::logic_error(message);
        };
        const auto require_error=[]<class Error>(auto&& operation,std::string_view message,
                                                 std::type_identity<Error>) {
            try {
                operation();
            } catch (const Error& error) {
                if (std::string_view(error.what())!=message) throw;
                return;
            }
            throw std::logic_error("Expected track-count offset validation failure");
        };
        constexpr auto invalid=std::type_identity<std::invalid_argument>{};
        constexpr auto range=std::type_identity<std::out_of_range>{};
        constexpr auto missing_values="Missing or mistyped flat buffer jet.ntrk";
        constexpr auto missing_offsets="Missing or mistyped flat buffer jet.ntrk.o";
        constexpr auto malformed="Invalid jagged offsets jet.ntrk.o";
        constexpr auto outside="Jagged row outside flat buffer";

        const Columns absent(std::span<const MetColumn>{});
        ChunkJetTrackCounts empty_chunk;
        require(empty_chunk.constructor_calls()==0,"An empty chunk validated track offsets");
        require_error([&] { empty_chunk.bind(absent); },missing_values,invalid);
        require(empty_chunk.constructor_calls()==0,"Missing values did not win before offset validation");

        const std::array values{std::int32_t{5},std::int32_t{7},std::int32_t{11}};
        const MetColumn flat{"jet.ntrk","i4",values.size(),values.data()};
        const std::array missing_specs{flat};
        const Columns missing(missing_specs);
        ChunkJetTrackCounts no_offsets;
        require_error([&] { no_offsets.bind(missing); },missing_offsets,invalid);
        require(no_offsets.constructor_calls()==1,"Missing offsets should be checked on first bind");

        const auto offsets_column=[](auto& offsets) {
            return MetColumn{"jet.ntrk.o","i8",offsets.size(),offsets.data()};
        };
        const std::array<std::int64_t,1> no_jet_offsets{0};
        const std::array no_jet_specs{MetColumn{"jet.ntrk","i4",0,nullptr},offsets_column(no_jet_offsets)};
        const Columns no_jets(no_jet_specs);
        ChunkJetTrackCounts no_jet_chunk;
        no_jet_chunk.bind(no_jets);
        no_jet_chunk.bind(no_jets);
        require(no_jet_chunk.constructor_calls()==1,"A nonempty chunk with no jets should validate once");

        const std::array<std::int64_t,4> offsets{0,2,2,3};
        const std::array good_specs{flat,offsets_column(offsets)};
        const Columns good(good_specs);
        ChunkJetTrackCounts good_chunk;
        {
            const auto [storage,rows]=good_chunk.bind(good);
            const auto first=rows.row(storage,0);
            require(first.size()==2 && first[0]==5 && first[1]==7,"First track-count row");
            require(rows.row(storage,1).empty(),"Duplicate offsets must give an empty row");
            require(rows.row(storage,2).size()==1 && rows.row(storage,2).front()==11,"Last track-count row");
            const auto [again_storage,again_rows]=good_chunk.bind(good);
            require(again_storage.data()==storage.data() && &again_rows==&rows &&
                    good_chunk.constructor_calls()==1,"Repeated event should retain the validated view");
        }
        ChunkJetTrackCounts next_chunk;
        next_chunk.bind(good);
        require(next_chunk.constructor_calls()==1,"Next chunk should validate its own offsets");

        const std::array<std::int64_t,3> decreasing{0,2,1};
        const std::array malformed_specs{flat,offsets_column(decreasing)};
        const Columns bad(malformed_specs);
        ChunkJetTrackCounts bad_chunk;
        require_error([&] { bad_chunk.bind(bad); },malformed,invalid);
        require(bad_chunk.constructor_calls()==1,"Decreasing offsets validation count");

        const std::array empty_specs{flat,MetColumn{"jet.ntrk.o","i8",0,nullptr}};
        const Columns no_offset_rows(empty_specs);
        ChunkJetTrackCounts empty_offsets;
        require_error([&] { empty_offsets.bind(no_offset_rows); },malformed,invalid);

        const std::array<std::int64_t,3> long_offsets{0,1,4};
        const std::array long_specs{flat,offsets_column(long_offsets)};
        const Columns too_long(long_specs);
        ChunkJetTrackCounts excess;
        const auto [storage,rows]=excess.bind(too_long);
        require(rows.row(storage,0).size()==1 && rows.row(storage,0).front()==5,
                "A later out-of-range row must not move ahead of earlier valid rows");
        require_error([&] { rows.row(storage,1); },outside,range);
        require_error([&] { rows.row(storage,2); },outside,range);
        require(excess.constructor_calls()==1,"Per-row range failure must not revalidate offsets");
        return 0;
    } catch (const std::exception& error) {
        error_message=error.what();
        return 1;
    }
}
