# ROOT-compatible random source

`root_rng.hpp` adapts the nonzero-seed uniform and ACR functions from the muon campaign's `correctionlib/atlas/muon/momentum/root_compat.h` (2026-11-17 campaign source, SHA256 `52a179d532a4ab60759955da883b2bc35cc6fcd72d8ffb8fe8135b20100e9ba1`). The prior campaign verified that implementation against the ROOT 6.40.02 primitives. This component calls one Gaussian after deriving the MET-specific seed, which needs its own independent proof.

- The original ROOT sources are from [ROOT 6.40.02](https://github.com/root-project/root/tree/d4c0003542952bec9c4edc86c419784a32875e9d), under the [ROOT license](https://github.com/root-project/root/blob/d4c0003542952bec9c4edc86c419784a32875e9d/LICENSE).
- [TRandom3](https://root.cern.ch/doc/v640/TRandom3_8cxx_source.html) credits Peter Malzacher's CLHEP adaptation and the MT19937 implementation by Makoto Matsumoto and Takuji Nishimura. The helper uses C++ `std::mt19937` and applies ROOT's nonzero-word conversion.
- [TRandom::Gaus](https://root.cern.ch/doc/v640/TRandom_8cxx_source.html#l00258) credits ACR to Wolfgang Hörmann and Gerhard Derflinger and refers to UNU.RAN, credited to Wolfgang Hörmann and Josef Leydold.
