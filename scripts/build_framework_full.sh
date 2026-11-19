#!/usr/bin/env bash
set -eu
repo=$(cd "$(dirname "$0")/.." && pwd)
out=$1
core=$2
pin=${3:-}
mkdir -- "$out"
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
set -euo pipefail
repo=$1
out=$2
core=$3
pin=$4
cd "$repo"
directories=(components/framework components/maker components/maker_adapter components/significance components/soft)
source_tree=$repo
if [[ -n $pin ]]; then
  mkdir "$out/source"
  git archive "$pin" "${directories[@]}" | tar -xf - -C "$out/source"
  source_tree=$out/source
  git rev-parse "$pin" > "$out/source-head.txt"
  printf "%s\n" "Build from immutable Git archive" > "$out/source-status.txt"
else
  git rev-parse HEAD > "$out/source-head.txt"
  git status --porcelain -- "${directories[@]}" > "$out/source-status.txt"
fi
(cd "$source_tree"; find "${directories[@]}" -type f -print0 | sort -z | xargs -0 sha256sum) > "$out/source-sha256.txt"
args=(
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=23 -DCMAKE_CXX_STANDARD_REQUIRED=ON
  -Dcorrectionlib_DIR="$core" -DAnalysisBase_DIR="$AnalysisBase_DIR/cmake" -DATLAS_MET_FULL_BRIDGE=ON
)
cmake -S "$source_tree/components/framework" -B "$out/build" "${args[@]}" > "$out/configure.log" 2>&1
cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
ldd "$out/build/libAtlasMetFramework.so" > "$out/bridge.ldd"
' bash "$repo" "$out" "$core" "$pin" > "$out/setup.log" 2>&1
