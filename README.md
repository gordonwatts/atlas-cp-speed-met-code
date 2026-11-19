# ATLAS missing transverse momentum with CorrectionLib

This repository contains the implementation used to compare the ATLAS CP tools in AnalysisBase 25.2.117 with a CorrectionLib replacement inside AnalysisBase and a standalone Python analysis reading the same PHYSLITE file. It covers MET reconstruction, track soft-term variations and Final significance with nominal upstream objects. The delivered implementation is [51d6b22](https://github.com/gordonwatts/atlas-cp-speed-met-code/commit/51d6b22bfcfb1b097952f101266cb2ed1aa6b1e9). The [report, slides and run logs](https://github.com/gordonwatts/atlas-cp-speed-met) are in the companion repository.

## Repository map

| Path | What to look for |
| --- | --- |
| [config](config/) and [observer](observer/) | Stock AnalysisBase configuration and the independent trace used to record the release’s behavior. |
| [components/maker](components/maker/) | Reconstructs hard terms and the track soft term. Its result includes the ordered constituent links and weights used downstream. |
| [components/soft](components/soft/) | Exports the track soft-term calibration, evaluates its systematic variations and random draws, and builds Final MET. |
| [components/significance](components/significance/) | Exports object resolutions and accumulates covariance in the varied Final direction. |
| [components/maker_adapter](components/maker_adapter/) and [components/framework](components/framework/) | Register the event-level CorrectionLib nodes and call them from the AnalysisBase sequence. |
| [components/preparation](components/preparation/) and [components/nominal_jvt](components/nominal_jvt/) | Prepare nominal electron, photon, muon and tau inputs and check the jet vertex tagger against the release. |
| [components/columnar](components/columnar/) | Standalone PHYSLITE reader, typed input preparation, native event loop and tests. Start at [met_demo.py](components/columnar/met_demo.py) and follow the call into [transport.cpp](components/columnar/transport.cpp). |
| [scripts](scripts/) and [fixtures](fixtures/) | Build and packaging recipes, reference jobs, acceptance comparisons, matched benchmarks and numerical fixtures. The [core precision patch](scripts/precision-correctionset.patch) accompanies the pinned CorrectionLib build. |
| [scripts/profiling](scripts/profiling/) | Elapsed-phase measurements, ROOT read accounting and sampled native CPU profiles. |
| [historical/table-check](historical/table-check/) | Preserved scripts for the earlier 747-probe soft-term table comparison. The subsequent full workflow has its own independent release reference. |

## Optimization history

The original 56 campaign commits are ancestors of `main`. [6b19c47](https://github.com/gordonwatts/atlas-cp-speed-met-code/commit/6b19c477e8d5578bb2c8c3a5e0eccc11e3437759) is the faithful standalone implementation. Follow its [changes through the delivered build](https://github.com/gordonwatts/atlas-cp-speed-met-code/compare/6b19c47...51d6b22) to see the retained optimizations and the trials that were reversed. The [experiment register](https://github.com/gordonwatts/atlas-cp-speed-met/blob/main/met-20261119/optimization.md) links each decision to its measurement.

## Running the result

For the tested original-file analysis, use the [offline distribution and instructions](https://github.com/gordonwatts/atlas-cp-speed-met/tree/main/met-20261119/standalone-offset). It targets CPython 3.14 on x86-64 Linux and requires the specified PHYSLITE input. For a source build, the original recipes use the AnalysisBase release, a [pinned CorrectionLib checkout](https://github.com/gordonwatts/correctionlib/commit/3dcfa20ec01fb40c4d5fec58383785519dec1c00), exported calibration payloads and the [small-R jet provider](https://github.com/gordonwatts/atlas-cp-speed-jets-code/tree/main/provider). The recorded shell recipes contain campaign-local dependency and result paths. The [report’s reproduction section](https://github.com/gordonwatts/atlas-cp-speed-met/blob/main/met-20261119/report.md#code-and-how-to-run-it) gives the original setup, timed endpoints and standalone invocation.
