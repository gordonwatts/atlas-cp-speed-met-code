#!/usr/bin/env bash
set -eu
repo=$(cd "$(dirname "$0")/.." && pwd)
output=$1
core_build=${2:-$repo/results/core-build-001}
source_tree=${3:-$repo/vendor/correctionlib}
mkdir -- "$output"
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
export PYTHONPATH="$1/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
set -eu
cd "$1"
out=$2
prior=/home/gwatts/atlas-cp-speed/jets-skill-c8a2374-20261117
core="$3/install/correctionlib"
source=$4
git -C "$prior" diff --exit-code HEAD -- provider/small_jvt.h provider/small_jvt.cxx portable/lwtnn.h exporters/small_jvt.py
git -C "$prior" rev-parse HEAD > "$out/reused-source-head.txt"
sha256sum "$prior"/{provider/small_jvt.h,provider/small_jvt.cxx,portable/lwtnn.h,exporters/small_jvt.py} > "$out/reused-source-sha256.txt"
python scripts/export_nnjvt_diagnostic.py "$prior/exporters/small_jvt.py" /cvmfs/atlas.cern.ch/repo/sw/database/GroupData results/reference-observer-full-001/full_config.json "$out/nnjvt.json"
flags=(
  -I "$core/include" -I "$prior/provider"
  -I "$source/lwtnn/include" -I "$source/eigen" -I "$source/rapidjson/include"
)
links=(
  "$3/build/lwtnn/lib/liblwtnn-stat.a"
  -L "$core/lib" -Wl,-rpath,"$core/lib" -lcorrectionlib -Wl,-z,defs
)
g++ -std=c++20 -O2 -fPIC -shared "${flags[@]}" components/nominal_jvt/evaluate.cxx "$prior/provider/small_jvt.cxx" "${links[@]}" -o "$out/libmet_nnjvt_probe.so"
g++ -print-file-name=libstdc++.so.6 > "$out/cxx-runtime-path.txt"
sha256sum "$out/libmet_nnjvt_probe.so" "$out/nnjvt.json"
' bash "$repo" "$output" "$core_build" "$source_tree" > "$output/build.log" 2>&1
