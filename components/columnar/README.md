# Direct PHYSLITE MET transport

The standalone program starts from the original file with UUID `79ee6652-311f-11f1-aa37-0445010abeef`. `met_demo.py` shows the named uproot projections, input chunk and state axes, binding calls and histogram fill. NumPy and Awkward create typed flat buffers and event/owner offsets. Native code loops over events, evaluates the installed CorrectionLib operations and returns dense term/source/significance rows. There are no Python loops over events or physics objects.

The default path requires the independently loaded tau, muon, electron and photon preparation nodes. It opens the original ROOT file, resolves same-event owner keys using the raw EventFormat metadata, regenerates NNJVT and prepares all nominal candidates. The event entry then calls the shared maker, track-soft variation, Final and significance for each state. This uses the event-major, state-minor order of the original algorithm sequence. Original and current identities are distinct within the native maker. The input probe has no path from a reference ntuple or a saved A preparation row.

## Build and runnable raw smoke

Run in a fresh build directory with the pinned GCC 15.2 compiler and the runtime owner's exact core and lwtnn archives. No ROOT tool or header enters this build.

```bash
cd /home/gwatts/atlas-cp-speed/met-20261119
core="$PWD/results/core-exact-build-001/install/correctionlib"
sdk="$PWD/vendor/correctionlib-exact-001"
lw="$PWD/results/core-exact-build-001/build/lwtnn/lib/liblwtnn-stat.a"
jvt_src=/home/gwatts/atlas-cp-speed/jets-skill-c8a2374-20261117/provider
build="$PWD/results/columnar-build-full-001"
taskset -c 0 cmake -S components/columnar -B "$build" -Dcorrectionlib_DIR="$core/cmake" \
    -DATLAS_MET_SMALL_JVT="$jvt_src" -DATLAS_MET_PINNED_SDK="$sdk" -DATLAS_MET_LWTNN="$lw" \
    -DATLAS_MET_COLUMNAR_PREPARATION=ON
taskset -c 0 cmake --build "$build" --parallel 1

read -r input < config/input.txt
"$python314" components/columnar/raw_smoke.py "$input" "$fresh/raw.npz" \
    "$build/libatlas_met_columnar.so" "$build/maker_adapter/libatlas_met_maker_provider.so" \
    "$maker_json" "$nnjvt_json" --entries 16 --chunk-size 7
```

Supply the reviewed exact JVT source and the exported NNJVT JSON. Prefer the adopted exact-core replay export and record source and checksum from the runtime owner. Disable `ATLAS_MET_COLUMNAR_PREPARATION` for the standalone raw smoke build before the peer's unified preparation source is available.

The raw smoke has separately labeled empty hard-object views and checks the original memberwise jet links, object-owner transport, pinned NNJVT and maker. It cannot serve as a production path or produce the production histogram. For each basket row it checks the byte-count and memberwise framing, then reads big-endian key block followed by index block. Uproot's ordinary jetLink interpretation decodes the wrong layout on this original file. The audited independent ROOT comparison found zero discrepancies for this memberwise reading and for direct nested objectLinks reads over all 20,000 entries. The first raw C stage was run successfully by the runtime owner with 16 entries in release-free CPython 3.14.

## Complete calculation

Prepare the same correction payloads used for B and the peer's original-file nominal preparation payloads. Use a fresh installed runtime directory with this layout.

| File | Input |
| --- | --- |
| `native/libatlas_met_columnar.so`, `native/libatlas_met_maker_correction.so`, `native/libatlas_met_maker_provider.so`, `native/libmet_preparation_provider.so` | `cmake --install "$build" --prefix "$runtime"` with preparation enabled |
| `native/libcorrectionlib.so` | Adopted exact-core installation |
| `native/libstdc++.so.6`, `native/libgcc_s.so.1` | Matching GCC runtime for the wheel assembler |
| `payloads/maker.json`, `payloads/soft.json`, `payloads/significance.json` | Verified exported maker, soft and event-significance JSON |
| `payloads/nnjvt.json` | Verified exact-core shared NNJVT update |
| `payloads/preparation.json` | `combine_preparation.py` with the complete peer-exported tau, muon and egamma JSON definitions |
| `payloads/EventFormat.json` | The original-file metadata-only witness in `results/original-metadata-001/EventFormat.json` |

The tested native CorrectionLib API exposes scalar and fixed-width batch inputs. The shared maker/significance operations require typed variable-length event rows and output links beyond that Python ABI. The visible runner therefore loads the verified native library and calls its named preparation batch entry and event-state evaluator. Their C++ implementation explicitly loads the pinned CorrectionSets and uses the immutable scoped bound nodes, with no ROOT release dependency.

In a clean CPython 3.14 runtime, first supply `--runtime "$runtime"` on the source script. The process needs the matching libstdc++ directory on its library path for that source check. After assembly the wheel bootstrap loads its own C++ libraries before importing numerical packages.

```bash
"$python314" components/columnar/combine_preparation.py "$runtime/payloads/preparation.json" \
    "$tau_payload" "$muon_payload" "$egamma_payload"
taskset -c 0 "$python314" components/columnar/met_demo.py "$input" "$fresh/c.npz" \
    --runtime "$runtime" --stop 16 --chunk-size 7 \
    --trace "$fresh/c-status.tsv" --observables "$fresh/c-observables.npz"
```

Check ordered association and JVT actions against the independent stock trace. Compare event/term identity, stored kinematics, each significance output, exact event bins and all five NPZ arrays to the stock results. Exclude optional traces from timing and compare traced and ordinary outputs first. `--correlated-diagnostic` adds the non-recommended correlated variation separately. With `-n N` each pass reopens PHYSLITE, resets per-pass significance state and histograms, and rereads with the same labels and branches. Repeated passes must yield identical histograms, and only the final histogram is saved.

## Wheel and isolation

```bash
"$python314" components/columnar/assemble.py "$runtime" "$fresh/bundle"
clean=(env -u PYTHONPATH -u PYTHONHOME -u LD_LIBRARY_PATH -u LD_PRELOAD -u ROOTSYS \
    -u ROOT_INCLUDE_PATH -u CMAKE_PREFIX_PATH -u ATLAS_LOCAL_ROOT_BASE -u VIRTUAL_ENV)
"${clean[@]}" taskset -c 0 uv run --no-project --python "$python314" \
    "$fresh/bundle/met_demo.py" "$input" "$fresh/uv-16.npz" --stop 16
"${clean[@]}" taskset -c 0 uv run --no-project --python "$python314" \
    "$fresh/bundle/met_demo.py" "$input" "$fresh/uv-16-repeat.npz" --stop 16 -n 2
"${clean[@]}" "$python314" -m venv "$fresh/pip-env"
"${clean[@]}" "$fresh/pip-env/bin/python" -m pip install "$fresh/bundle"/wheels/*.whl
"${clean[@]}" taskset -c 0 "$fresh/pip-env/bin/python" \
    "$fresh/bundle/met_demo.py" "$input" "$fresh/pip-16.npz" --stop 16
```

Record public dependency pins, wheel hash and bytes, all declared and live native dependencies and outputs in the isolated runtime receipts. Match the isolated single/repeat and pip outputs to the accepted original-stream candidate after its independent A gate. This source checkpoint does not establish those full numerical or isolation gates.
