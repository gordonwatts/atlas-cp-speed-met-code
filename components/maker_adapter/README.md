# Shared MET maker and significance event nodes

The module `libatlas_met_maker_provider.so` registers the maker and significance event operations and the peer significance primitives. The same binary `libatlas_met_maker_correction.so` owns the C++ entry points used by both B and C and the scoped event context used inside the registered nodes. Link this shared library in consumers to share that context with the module. A combined provider may call the named registration functions from its API-v4 registration entry point.

`invoke_maker` and `invoke_significance` take the actual resolved `correction::Correction` and typed input view. They issue a private identity probe while binding the node instance, then pass a unique invocation ID to an ordinary installed CorrectionLib evaluation. The probe is unavailable outside its scoped call. Only that node instance on the current thread may consume the ID once. Nested evaluation restores the parent scope and failed calls clear the binding on every exit. Results and vectors remain in caller-owned evaluation storage. The numeric routines run inside the registered node.

C supplies type-correct owner/index keys, complete original link targets and immutable nominal prepared object rows independently from the original file. See the maker types in [maker.hpp](../maker/include/met_maker/maker.hpp) and the raw public significance types in [significance_correction.hpp](include/met_adapter/significance_correction.hpp). The significance maker supplies ordered constituent rows and uses the current constituent pT and p4 direction. Egamma inputs distinguish the first cluster eta, calorimeter-frame moment and auxiliary value, forward author, and photon vertex-link count. Muons require optional CB, ID and extrapolated MS track parameters. Taus require the decay mode and TES and final calibrated coordinates on their executed branches. C carries the returned significance state in tool-equivalent call order.

## Payloads and boundary smokes

Use the framework release environment for a matching core and the same compiler. With the runtime owner's single CPU and build job and a fresh output directory:

```bash
build="$PWD/results/framework-build-maker-check-001"
mkdir "$PWD/results/framework-payload-check-001"
payload="$PWD/results/framework-payload-check-001"
cuts=/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/JetPileupTag/NNJvt/2022-03-22/NNJVT.Cuts.FixedEffPt.Offline.Nonprompt_All_MaxW.json
python components/maker/export_nnjvt.py "$cuts" "$payload/nnjvt.json"
python components/maker_adapter/export_maker.py "$payload/nnjvt.json" "$payload/maker.json"
provider="$build/maker_adapter/libatlas_met_maker_provider.so"
taskset -c 0 "$build/maker_adapter/atlas_met_maker_adapter_smoke" "$provider" "$payload/maker.json"

# Point this at the peer's exported jet/egamma/muon/tau primitive definitions.
primitives=/absolute/exported/significance/corrections.json
python components/maker_adapter/export_significance.py "$primitives" "$payload/significance_event.json"
taskset -c 0 "$build/maker_adapter/atlas_met_significance_event_smoke" "$provider" \
    "$payload/significance_event.json" "$primitives"
```

The maker smoke covers failed-then-successful event calls, nested evaluation, unbound and stale tokens, cross-CorrectionSet tokens and repeated consumption. The significance smoke checks a forward-author photon with a vertex link using the electron coefficients and photon constant-term policy via the event node, and verifies that the contrasting converted-photon resolution is distinguishable. The significance event JSON refers to the peer-exported tables by absolute path. Regenerate it after relocating the frozen bundle or resolve those paths when assembling the portable deliverable. This is a source checkpoint. The runtime owner builds and executes the smoke and the independent B/C comparisons.

[CorrectionLib develop API-v4 source](https://github.com/gordonwatts/correctionlib/blob/3dcfa20ec01fb40c4d5fec58383785519dec1c00/include/correction_plugin.h) / [Stock significance call](https://gitlab.cern.ch/atlas/athena/-/blob/6ce3f7ae7c41dd9eaec170c3c349df15faef7e5c/Reconstruction/MET/METUtilities/Root/METSignificance.cxx) / [Egamma type mapping](https://gitlab.cern.ch/atlas/athena/-/blob/6ce3f7ae7c41dd9eaec170c3c349df15faef7e5c/PhysicsAnalysis/ElectronPhotonID/ElectronPhotonFourMomentumCorrection/Root/EgammaCalibrationAndSmearingTool.cxx#L936-964).
