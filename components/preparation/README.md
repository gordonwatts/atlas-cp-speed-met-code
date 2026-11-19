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

## Muon nominal stage

The installed calibration algorithm skips nominal muon calibration for this PHYSLITE recipe. Use the original pT, eta, phi and mass. The stored `InnerDetectorPt` and `MuonSpectrometerPt` enter the configured release Medium quality check.

`muon_schema.json` lists the input types, order and output components. Evaluate `muon_selection_eta(raw_eta)` and preserve its uint32 gate. Evaluate `muon_selection_vertex(vertex_type)` in raw vertex order to choose the first PriVtx; use zero z when none exists. `muon_selection_track_roles` returns the primary and ID-summary role codes `[absent=0, ID=1, extrapolated MS=2, combined=3, MS-only extrapolated=4, MS=5]`. Use the available xAOD link targets, including the forward ID collection. Collect raw track and muon hit fields by those roles. An unavailable covariance vector has size -1 and an unavailable hit summary is -1. Slot 0 is d0 variance. Slot 4 is q/p variance. The three precision-layer summaries can fall back to the primary track when the muon accessor is absent.

`muon_selection_after_calibration` returns `[pt_bits, track_bits, medium, bad_muon, isolation, selected]`. The track mask has low bits for track existence, d0 significance and projected z0. Unused bits are set. Skipped track is `0xfffffff8` and skipped pT is `0xfffffffe`. `track_iso` uses raw `ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000` and `calo_iso` uses raw `topoetcone20`. Their availability is recorded as low bits 0 and 1. Compare the signed float decorations separately to float-rounded `0.15 * pt` and `0.30 * pt` with inclusive thresholds. Both have the empty suffix. `calo_score` is `CaloMuonScore`. The `lowpt_parameters_mask` bits indicate the availability of `momentumBalanceSignificance`, `scatteringCurvatureSignificance` and `scatteringNeighbourSignificance`. The prior Medium and track implementations are traced in `muon.cxx`. Stock isolation is `IsolationConditionFormula` and checks the two Loose_VarRad conditions in order.

Use the same release environment and CMake setup as for taus with a new output directory. Compare all raw rows against the precomputed independent A observer using:

```bash
out="$PWD/results/muon-preparation-001"
mkdir -- "$out"
taskset -c 0 python components/preparation/export_muon.py \
  --reference-config results/reference-observer-full-001/full_config.json \
  --output-dir "$out/bundle" > "$out/export.log" 2>&1
taskset -c 0 cmake -S components/preparation -B "$out/build" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=23 -DCMAKE_CXX_STANDARD_REQUIRED=ON \
  -Dcorrectionlib_DIR="$MET_CORE_CONFIG" -DMET_PREPARATION_STOCK_CHECK=ON > "$out/configure.log" 2>&1
taskset -c 0 cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
read -r input < config/input.txt
taskset -c 0 "$out/build/met_muon_prepare_stock_check" \
  "$out/build/libmet_preparation_provider.so" "$out/bundle/muon.json" "$input" \
  results/reference-observer-full-001/reference.root > "$out/stock.log" 2>&1
```

The witness reads only the original file for provider inputs. It compares binary32 pT, eta, phi, raw eta/pT/track masks, Medium, Loose_VarRad and final decisions against A. The observer's `mu_m` field reports an unavailable auxiliary column as zero. The muon p4 uses the original release's fixed `Muon_v1::m()` instead. It independently probes signed isolation, equality thresholds and unavailable fields. The observer did not expose the unused bad-muon diagnostic.

## Electron and photon nominal stage

Use the same registered API4 provider and load the exported `egamma.json`. `egamma_schema.py` lists in order the real and integer fields for `electron_prepare` and `photon_prepare`, and the output components. Supply binary32 source decorations widened to reals. Missing link inputs follow the availability flags. The final three reals and final integer are private child-evaluation slots. Supply zero. The data are for the recorded MC campaign.

| Fields | Original owner |
| --- | --- |
| p4, OQ, author, ID, and photon cleaning | Original AnalysisElectrons or AnalysisPhotons. ID is the corresponding `DFCommonElectronsLHLooseBL` or `DFCommonPhotonsIsEMTight` char. Clean is `DFCommonPhotonsCleaning`. Missing char is -1 and fails when reached. |
| cluster kinematics, moments and samplings | First resolvable index-0 `caloClusterLinks` target in `egammaClusters`. The eta selection uses the computed layer-2 eta. Use actual raw eta and phi for cluster coordinates. For each of sampling slots 1/2/3 and 5/6/7, use the population count of lower pattern bits to find the packed slot in `e_sampl` or `eta_sampl`. Use -999 for a present sampling with an out-of-range packed index. Supply zero for a sampling absent from the pattern. Supply `has_calo_moments` only when both calorimeter-frame moments can be read. |
| event / PV | Actual `runNumber` and separate `RandomRunNumber` (-1 if missing) on EventInfo. Choose first PriVtx in original PrimaryVertices order as for muons. Set `has_pv=0` with no such vertex. Electron track uses beamPosSigmaX/Y/XY and PV z or zero. |
| calorimeter isolation | Per-cone raw `topoetcone20/30/40`, corresponding `topoetconeXXptCorrection`, and the `Eadded_Lr2/Lr3` per-object floats. Isolation/old-leak availability bits for cones 20/30/40 are 0/1/2. Timing availability bits are 0/1. |
| track isolation and electron track | Track isolation is `ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000` for electron or `ptcone20` for photon. Its availability is bit 3 of the isolation mask. For an electron use first index-0 `trackParticleLinks` target in `GSFTrackParticles`, including d0 variance from index 0 of the covariance diagonal and size -1 for absent. |
| photon conversion | First index-0 `vertexLinks` target in `GSFConversionVertices`, with x, y and optional conversion pt1 and pt2. Follow its first two trackParticleLinks to raw tracks. Pass separate link validity, pixel and SCT summaries (-1 for missing) and qOverP, theta, phi. The fixed schema has no first-measurement perigee arrays, so the stock fallback when pt1/pt2 are unavailable computes pT from perigee. |

The shared node reproduces the uint32 release accept masks, including partial bits of early eta, OQ, cleaning and track cuts. The final boolean follows the configured METMaker predicate. Photon origin correction rewrites float pT and eta after its eta, OQ and cleaning stages when a PriVtx exists. Electron nominal momentum remains the input momentum per the recorded calibration configuration. All reached topo cones receive the old-to-new leakage and Run 3 timing updates before track/ID/isolation selection.

The three native scalar children of the preparation node evaluate TF1 translations, the trouble-category TGraph, and the hData TProfile2D displacement. TFormula payloads preserve the active REL22 power/log and REL21 quadratic forms. Profile interpolation uses bin centres and native MultiBinning contents with the release's boundary and precision rules. `egamma_pointing_depth` and the four per-species leakage corrections also expose these lookups to standalone checks. The exporter validates initialized settings and the SHA256 of each resolved input file.

Use the same pinned release and build environment above with a fresh evidence directory:

```bash
out="$PWD/results/egamma-preparation-001"
mkdir -- "$out"
taskset -c 0 python components/preparation/export_egamma.py \
  --reference-config results/reference-observer-full-001/full_config.json \
  --output-dir "$out/bundle" > "$out/export.log" 2>&1
taskset -c 0 cmake -S components/preparation -B "$out/build" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=23 -DCMAKE_CXX_STANDARD_REQUIRED=ON \
  -Dcorrectionlib_DIR="$MET_CORE_CONFIG" -DMET_PREPARATION_STOCK_CHECK=ON > "$out/configure.log" 2>&1
taskset -c 0 cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
read -r input < config/input.txt
taskset -c 0 "$out/build/met_egamma_prepare_stock_check" \
  "$out/build/libmet_preparation_provider.so" "$out/bundle/egamma.json" \
  "$out/bundle/egamma-probes.tsv" "$input" results/reference-observer-full-001/reference.root > "$out/stock.log" 2>&1
```

The witness checks native table/graph evaluations at ROOT-derived probes and then compares the provider with the A observer on all original electron and photon rows. Electron and photon mass uses the release's `m()` accessor because the A auxiliary mass branches are zero-filled. It runs the independent release photon-origin and IsolationCorrectionTool on reached rows to compare every corrected cone, including stored leakage and timing effects. Seventeen source-derived composed fixtures cover partial accept masks, eta and pT equality, two dead-HV periods, cleaning, electron d0, signed isolation, missing or skipped timing and conversion category dispatch. The original file has electron cone 20 and photon cones 20/40, so the other exported cones are checked by their independent table probes. [ROOT 6.40 TH2 interpolation](https://root.cern.ch/doc/v640/TH2_8cxx_source.html#l01338)
