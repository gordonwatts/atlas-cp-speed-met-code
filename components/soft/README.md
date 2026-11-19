# MET track soft term and Final

This component implements the fixed AnalysisBase 25.2.117 Tight EMPFlow track-soft-term prescription. It exports the original TEnv and TH3D inputs, applies the configured states through CorrectionLib API v4, and builds the ordered float-stored Final term. ROOT is used to export the source calibration and run witness tests.

## Framework and standalone interface

Add this directory with `add_subdirectory`, link `atlas::met_soft`, and include `soft.hpp`. Create an immutable `atlas::met::soft::SoftFile` from the generated `corrections.json` once per payload. Both adapters call the same methods.

1. Supply `EventView` with ordered rebuilt terms, original truth terms and `EventInfo.actualInteractionsPerCrossing()`. For the first term matching the jet source mask, pass the size of its raw `ConstitObjectLinks` vector, including invalid links. An available empty vector and a missing decoration have distinct inputs. Carry original source masks for every term. Projection to these records does not change the reconstruction term links and weights.
2. Call `apply(event, state)` with the initialized systematic string or `NOSYS`. The result identifies the `PVSoftTrk` ordinal, supplies the varied storage components, and records the signed and effective seed, recoil disposition, counts, selected coefficients, and draw status. Leave the associated metadata and constituent links attached to the same row. A C++ RAII binding carries the typed input into the registered node. It checks the thread-local ticket and payload instance, rejects stale and nested bindings, and clears the binding even if evaluation throws. The node executes the three standard MultiBinning definitions as ordered scalar children.
3. Pass the current ordered terms, with the selected soft row updated, to `build_final`. The node returns the storage components and float magnitude and azimuth. Append the result as `Final` with source `0x10040`. For additional decorations on other terms use `decorate(momentum)`, the shared two-output node.

The scope is nominal and the four recommended states, plus `MET_SoftTrk_ResoCorr` for a separate diagnostic. Unknown state names, non-TST soft prescriptions and configured calibration differences fail explicitly. Final requires a selected soft term and rejects a duplicate Final name. For recoil the stock function subtracts every non-soft, non-total term, including an exact invisible source if present. Final skips total and invisible sources and keeps soft terms matching the chosen soft source mask. It accumulates three float aux fields in container order.

The raw seed is `int(1.e6 * phi32)` from the pre-correction soft term. The pinned release binary's `CP::MetMakerAlg::execute` shows double `atan2`, narrowing to float, widening for multiplication, integer truncation and sign extension into the unsigned-long setter. The provider maps signed zero to effective seed one and negative int seeds to their low 32 bits. Nominal prepares the seed without a Gaussian. Scale and resolution states each draw once, and the correlated diagnostic uses that draw for both axes. A nonfinite azimuth cannot enter a defined release integer cast, so the component reports `undefined_seed`. [ROOT-compatible RNG provenance](NOTICE.md)

Missing truth and negative hard sumET return a zero recoil to the original projection. Zero recoil can yield NaN while the correction status is successful. Only the source's tiny nonnegative sumET branch floors the recoil vector; ordinary positive activity may still have zero transverse recoil. A missing first jet source is undefined in the release's dereferenced iterator, and a missing jet link accessor is reported separately. A missing EventInfo is a framework read-handle error, even though the correction source has a later null-pointer guard. Neither case should publish a corrected term.

## Export and checked activation

The actual full initialized configuration is `results/reference-observer-full-001/full_config.json` in the fresh checkout. [Relevant settings for orientation](example-stock-properties.json)

Run from `/home/gwatts/atlas-cp-speed/met-20261119` using the fresh same-compiler CorrectionLib API v4 core. Set `MET_CORE_CONFIG` to the containing directory of the installed `correctionlibConfig.cmake`. Choose an unused result directory and retain each log.

```bash
source scripts/release_environment.sh
export PYTHONPATH="$PWD/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
export ROOT_MAX_THREADS=1
out="$PWD/results/soft-component-001"
mkdir -- "$out"
taskset -c 0 python components/soft/export_soft.py export \
  --reference-config results/reference-observer-full-001/full_config.json \
  --output-dir "$out/bundle" > "$out/export.log" 2>&1
taskset -c 0 cmake -S components/soft -B "$out/build" \
  -DCMAKE_BUILD_TYPE=Release -Dcorrectionlib_DIR="$MET_CORE_CONFIG" \
  -DMET_SOFT_STOCK_CHECK=ON > "$out/configure.log" 2>&1
taskset -c 0 cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
taskset -c 0 "$out/build/met_soft_smoke" "$out/bundle/corrections.json" \
  "$out/bundle/tables.json" "$out/bundle/table-probes.tsv" > "$out/smoke.log" 2>&1
taskset -c 0 "$out/build/met_soft_stock_compare" "$out/bundle/corrections.json" \
  results/reference-status-full-001/status.tsv > "$out/stock.log" 2>&1
ldd "$out/build/met_soft_smoke" > "$out/native.ldd"
```

The export verifies the historical source hashes against the freshly resolved files, checks exact axis agreement across all three TH3Ds, and records source objects and output hashes. It emits 774 stock table probes per load path, including the 747 historical center, adjacent-edge and flow cases plus 27 infinity and NaN cases. The native smoke checks both loaders bitwise and dispatches all runtime numerical operations through the pinned native CorrectionLib core. In an approved portable environment run the same smoke and record its loaded libraries. When the pinned Python binding is available, the optional Typer `check-tables --bundle "$out/bundle" --output "$out/python-tables.json"` tests that binding's two load paths independently.

The stock comparison runs the actual release soft tool and Final helper on controlled input terms. It compares binary32 output components and decorations, recoil branches, and known seeds. It also replays the signed seed, uncorrected soft inputs and binary64 Gaussian observations already present in the full A status trace. These checks are separate from real-input B/C acceptance for corrected terms and Final results.
