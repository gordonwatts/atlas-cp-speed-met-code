#!/usr/bin/env bash
set -eu
repo=$(cd "$(dirname "$0")/.." && pwd)
out=$1
pin=$2
core=$3
mkdir -- "$out"
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
export PYTHONPATH="$1/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
set -euo pipefail
repo=$1
out=$2
pin=$3
core=$4
cd "$repo"
input=$(cat config/input.txt)
mkdir "$out/source"
git archive "$pin" components/preparation components/significance | tar -xf - -C "$out/source"
git rev-parse "$pin" > "$out/source-head.txt"
component="$out/source/components/preparation"
stock_config="$repo/results/reference-observer-full-001/full_config.json"
stock_ntuple="$repo/results/reference-observer-full-001/reference.root"
for species in tau muon egamma; do
  python "$component/export_$species.py" --reference-config "$stock_config" --output-dir "$out/$species-bundle" > "$out/$species-export.log" 2>&1 || echo "Failed $species export"
done
args=(
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=23 -DCMAKE_CXX_STANDARD_REQUIRED=ON
  -Dcorrectionlib_DIR="$core" -DMET_PREPARATION_STOCK_CHECK=ON
)
cmake -S "$component" -B "$out/build" "${args[@]}" > "$out/configure.log" 2>&1 || exit
cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1 || exit
provider="$out/build/libmet_preparation_provider.so"
ldd "$provider" > "$out/provider.ldd"
for species in tau muon egamma; do
  payload="$out/$species-bundle/$species.json"
  if [[ -f $payload ]]; then
    probe=()
    witness=()
    if [[ $species == egamma ]]; then probe=("$out/$species-bundle/egamma-probes.tsv"); fi
    if [[ $species != tau ]]; then witness=("$stock_ntuple"); fi
    status=0
    "$out/build/met_${species}_prepare_stock_check" "$provider" "$payload" "${probe[@]}" "$input" "${witness[@]}" > "$out/$species-stock.log" 2>&1 || status=$?
    printf "%s %s\n" "$species" "$status" >> "$out/status.txt"
  fi
done
' bash "$repo" "$out" "$pin" "$core" > "$out/job.log" 2>&1
