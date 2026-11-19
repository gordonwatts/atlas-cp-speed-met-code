#!/usr/bin/env bash
set -eu
repo=$(cd "$(dirname "$0")/.." && pwd)
out=$1
core=$2
mkdir -- "$out"
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
set -eu
repo=$1
out=$2
core=$3
cd "$repo"
git rev-parse HEAD > "$out/source-head.txt"
git status --porcelain -- components/framework components/maker components/maker_adapter components/significance components/soft > "$out/source-status.txt"
args=(
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=23 -DCMAKE_CXX_STANDARD_REQUIRED=ON
  -Dcorrectionlib_DIR="$core" -DAnalysisBase_DIR="$AnalysisBase_DIR/cmake" -DATLAS_MET_FULL_BRIDGE=ON
)
cmake -S components/framework -B "$out/build" "${args[@]}" > "$out/configure.log" 2>&1
cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
ldd "$out/build/libAtlasMetFramework.so" > "$out/bridge.ldd"
' bash "$repo" "$out" "$core" > "$out/setup.log" 2>&1
