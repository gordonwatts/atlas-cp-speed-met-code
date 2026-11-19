# MET significance component

The fixed AnalysisBase 25.2.117 Run 3 full-simulation Tight EMPFlow track-soft-term recipe uses the Aug2018 data and MC JER tables, Run 3 egamma with distinct electron and photon constant-term settings, Data22 muon expected ID/MS/CB resolutions, and the CombinedTES tau FinalCalib graphs. The exporter checks the original initialized configuration, resolves ROOT inputs with PathResolver, and records their checksums and object identities. The numerical library uses the exported tables and the standard C++ library.

## Integrating the provider

- Add this directory to the parent CMake project with `add_subdirectory`. Link `met_significance_numeric` into the module that implements the event node. Link `met_significance_nodes` and call `register_significance_nodes(registry)` if the individual CorrectionLib nodes are also exposed there. The independently loadable `met_significance_provider` target owns its own v4 dynamic registration symbols.
- Export the payload for the validated stock configuration. Load the four versioned `.table` sections with `read_jet_data`, `read_egamma_data`, `read_muon_data` and `read_tau_data` during immutable event-node initialization. The shared event node calls the functions in `significance.hpp`. The standalone `corrections.json` in the same bundle inlines identical bytes from those four sections for v4 component tests.
- B and C supply the linked nominal object inputs and the ordered MET term views. Use the actual linked object pT and p4 Cartesian x/y in MeV. Retain the type, original identity and validity independently of the MET weights. Use `term_action` and `uses_constituent` to avoid resolving skipped terms. Narrow the output of each resolution to float once before passing it to the reducer. The jet helper already includes that stage.
- Supply muon track link kinematics via `track_pt_gev`, the raw run number and simulation flag, and muon object eta and phi. A missing track is distinct from a linked track with zero qOverP. An optional tau decay mode preserves the invalid/missing mode return. An egamma call needs the linked first cluster and its eta. Use the particle energy at the selected nominal momentum. For photons, pass `with_constant_term=true` and derive the class from the photon vertex link count and forward-electron author branch. Obtain etaCalo from the moment then the float aux, with the stock forward-author override. Absence of both frame fields for a non-forward photon raises because the release returns an uninitialized value on that path. The electron call passes `with_constant_term=false`.
- Keep one `State` per tool-equivalent call sequence, including across ordered MET states and events. Treat null constituent links and absent link decorations as failures. An absent Final is a framework failure, a wrong-source Final still returns success in the stock tool. A wrong-source Final or zero longitudinal variance preserves the previously stored squared significance and rho. Later track-soft terms return a diagnostic count. The four decorations are already narrowed to float. Relative-resolution float roundings and current-state angles affect covariance even with fixed nominal object inputs.
- Export rejects alternate initialized physics settings. This component targets reconstructed species, the fixed soft-term recipe, and Data22 full simulation. When integrating, enforce matching live recipe properties and reject other prescriptions rather than invoking these nodes with another configuration. Maintain an independent stock reference for every active branch and configured data dataset.

`jet_resolution` returns the data, MC, and float max as an array. `egamma_resolution`, `muon_resolution`, and `tau_resolution` return the release getters' double resolutions. `significance_finish` returns the four float decorations promoted to double, then the updated squared significance and rho. The schema-v3 payload declares the field layouts. The event traversal that consumes variable-length links stays behind the campaign's event node, which shares the same numerical library.

## Runtime-owner verification

Run from the fresh MET checkout in a release shell. Use the installed `correctionlibConfig.cmake` directory from a same-compiler build of the current fork pin `3dcfa20ec01fb40c4d5fec58383785519dec1c00`. The shell variable `MET_CORE_CONFIG` below is that directory. Keep the build and outputs in a new result directory. The commands below assume CPU 0 as in the current stock reference.

```bash
cd /home/gwatts/atlas-cp-speed/met-20261119
source scripts/release_environment.sh
export PYTHONPATH="$PWD/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
out="$PWD/results/significance-component-001"
mkdir -- "$out"
taskset -c 0 python components/significance/export_resolution.py \
  --reference-config results/reference-observer-full-001/full_config.json \
  --output-dir "$out/bundle" > "$out/export.log" 2>&1

taskset -c 0 cmake -S components/significance -B "$out/build" \
  -DCMAKE_BUILD_TYPE=Release -Dcorrectionlib_DIR="$MET_CORE_CONFIG" \
  -DMET_SIGNIFICANCE_STOCK_CHECK=ON > "$out/configure.log" 2>&1
taskset -c 0 cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
taskset -c 0 "$out/build/met_significance_smoke" "$out/build/libmet_significance_provider.so" \
  "$out/bundle/corrections.json" > "$out/smoke.log" 2>&1
read -r input < config/input.txt
taskset -c 0 "$out/build/met_significance_stock_compare" \
  "$out/build/libmet_significance_provider.so" "$out/bundle" "$input" 64 > "$out/stock.log" 2>&1
ldd "$out/build/libmet_significance_provider.so" > "$out/provider.ldd"
```

The stock check uses actual release public resolution calls, edge and knot fixtures, original input objects, and a synthetic MET container with real links and source-ordered terms. A row reports its witness count, mismatches and maximum absolute error. The smoke also checks the state sentinel and provider load. Repeat the smoke with the same CorrectionLib core and compiler-family runtime in the portable environment without the release setup. Confirm the dynamic dependency inventory there. The stock witness supplements full A/B and A/C comparisons of the prepared constituent identities and the saved decorated result.

The original source trace and discussion of edge behavior are in the Windows campaign's `met-20261119/resolution-design.md`. In particular, jet projection construction follows `JetSmearingCorrection`, the muon component follows `getExpectedResolution` rather than the momentum smearing call, and tau graph points retain ROOT order and the sorted-X flag. The historical JER smearing stream and muon calibration stream are outside these deterministic resolution getters.
