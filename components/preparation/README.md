# Standalone nominal object preparation

The exported definitions use the pinned CorrectionLib API v4. Link `atlas::met_preparation` and call `atlas::met::preparation::register_preparation_nodes(registry)` before loading the preparation payloads, or load `libmet_preparation_provider.so` into a dedicated registry. Numerical evaluation stays in CorrectionLib. The adapter reads original PHYSLITE branches and resolves event-local links, preserving the original row order and link identities.

## Tau nominal stage

Load `tau.json` and call `tau_track_charged(flag_set)` on each resolvable TauTrack in an original tau's `tauTrackLinks`. Invalid/thinned links contribute zero. The scalar correction tests `TauTracksAuxDyn.flagSet` bit 5. Sum its boolean outputs by owning tau. This is the number of classified charged tracks used by both the stock scale tool and Tight.

| `tau_prepare` inputs in order | Original values |
| --- | --- |
| `pt, eta, phi, mass, charge` | `AnalysisTauJetsAuxDyn.{pt,eta,phi,m,charge}` in MeV where dimensional |
| `pt_final, eta_final, pt_tes, eta_tes` | `ptFinalCalib, etaFinalCalib, ptTauEnergyScale, etaTauEnergyScale` |
| `charged_tracks, is_tau_flags, pan_tau_mode, has_tes_pt` | Count above, widened unsigned isTauFlags, PanTau_DecayMode or -1, and availability of TES pT |
| `tes_compatibility` | Persisted char converted to bool or -1 for absent. The frozen file does not store this decoration. |
| `truth_pdg_id, truth_is_hadronic, truth_jet_link_available` | Resolve the original `truthParticleLink` key and index. Use PDG ID 0 if unresolved. For a resolved true tau use its `IsHadronicTau`, and -1 for unavailable. The final bool is presence of `truthJetLink`, distinct from validity. |

The fixed output is `[pt, eta, phi, mass, selected, passed_cut_bits, TESCompatibility, smear_branch, scale]`. Bits 0–5 record stock cut order `Pt, AbsEta, NTrack, AbsCharge, JetIDWP, EleIDWP`. Later bits remain false after the first failure. Only the overall bool enters `selected_tau_tight,as_char`. Branch values 0, 1 and 2 identify nonhadronic truth, wrong charged-prong count, and nominal pT rewrite. Source with `IsTruthMatched` does not substitute for the link test. The current p4 rewrite uses `float(ptFinalCalib * sf)` and retains eta, phi and mass. The compatibility check uses the Calo-only graphs and does not gate the six active Tight cuts.

The exporter verifies the configured Tight file, run period, truth/insitu switches and graph source. The nominal SF uses the last map instance seen in ROOT key iteration for each prong count. It preserves `FindFixBin` after the source axis protection, and the nominal fade evaluated at unclamped final pT in GeV. The Calo graph extraction shares the significance library's ROOT-style graph evaluator.

```bash
cd /home/gwatts/atlas-cp-speed/met-20261119
source scripts/release_environment.sh
export PYTHONPATH="$PWD/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
out="$PWD/results/tau-preparation-001"
mkdir -- "$out"
taskset -c 0 python components/preparation/export_tau.py \
  --reference-config results/reference-observer-full-001/full_config.json \
  --output-dir "$out/bundle" > "$out/export.log" 2>&1
taskset -c 0 cmake -S components/preparation -B "$out/build" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=23 -DCMAKE_CXX_STANDARD_REQUIRED=ON \
  -Dcorrectionlib_DIR="$MET_CORE_CONFIG" -DMET_PREPARATION_STOCK_CHECK=ON > "$out/configure.log" 2>&1
taskset -c 0 cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
read -r input < config/input.txt
taskset -c 0 "$out/build/met_tau_prepare_stock_check" \
  "$out/build/libmet_preparation_provider.so" "$out/bundle/tau.json" "$input" > "$out/stock.log" 2>&1
ldd "$out/build/libmet_preparation_provider.so" > "$out/provider.ldd"
```

`MET_CORE_CONFIG` is the directory containing `correctionlibConfig.cmake` for the current pinned evaluator, built with the same compiler. The stock check runs separately through the release TauSmearingTool and TauSelectionTool on every original tau and compares all six cut results, Tight, TES compatibility, charged track count and binary32 p4. The stock check uses release-xAOD links as a component oracle. Standalone input decoding still needs its own comparison.
