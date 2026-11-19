#!/usr/bin/env bash
set -eu
repo=$(cd "$(dirname "$0")/.." && pwd)
out=$1
pin=$2
mkdir -- "$out"
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
export PYTHONPATH="$1/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
set -euo pipefail
repo=$1
out=$2
pin=$3
cd "$repo"
mkdir "$out/source"
directories=(components/columnar components/maker components/maker_adapter components/significance components/soft components/preparation)
git archive "$pin" "${directories[@]}" | tar -xf - -C "$out/source"
git rev-parse "$pin" > "$out/source-head.txt"
core="$repo/results/core-exact-build-001/install/correctionlib"
args=(
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=23 -DCMAKE_CXX_STANDARD_REQUIRED=ON
  -Dcorrectionlib_DIR="$core/cmake"
  -DATLAS_MET_SMALL_JVT=/home/gwatts/atlas-cp-speed/jets-skill-c8a2374-20261117/provider
  -DATLAS_MET_PINNED_SDK="$repo/vendor/correctionlib-exact-001"
  -DATLAS_MET_LWTNN="$repo/results/core-exact-build-001/build/lwtnn/lib/liblwtnn-stat.a"
  -DATLAS_MET_COLUMNAR_PREPARATION=ON
)
cmake -S "$out/source/components/columnar" -B "$out/build" "${args[@]}" > "$out/configure.log" 2>&1
cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
runtime="$out/runtime"
cmake --install "$out/build" --prefix "$runtime" > "$out/install.log" 2>&1
install -D "$core/lib/libcorrectionlib.so" "$runtime/native/libcorrectionlib.so"
install "$(g++ -print-file-name=libstdc++.so.6)" "$runtime/native/libstdc++.so.6"
install "$(g++ -print-file-name=libgcc_s.so.1)" "$runtime/native/libgcc_s.so.1"
mkdir "$runtime/payloads"
cp results/maker-component-001/maker.json "$runtime/payloads/maker.json"
cp results/soft-component-exact-001/bundle/corrections.json "$runtime/payloads/soft.json"
cp results/significance-event-payload-001/met_significance_event.json "$runtime/payloads/significance.json"
cp results/nnjvt-component-exact-001/nnjvt.json "$runtime/payloads/"
cp results/original-metadata-001/EventFormat.json "$runtime/payloads/"
python "$out/source/components/columnar/combine_preparation.py" "$runtime/payloads/preparation.json"     results/preparation-frozen-001/tau-bundle/tau.json     results/preparation-frozen-001/muon-bundle/muon.json     results/egamma-preparation-export-002/bundle/egamma.json
(cd "$runtime"; find native payloads -type f -print0 | sort -z | xargs -0 sha256sum) > "$out/runtime-sha256.txt"
env LD_LIBRARY_PATH="$runtime/native" ldd "$runtime/native/libatlas_met_columnar.so" > "$out/runtime.ldd"
' bash "$repo" "$out" "$pin" > "$out/job.log" 2>&1
