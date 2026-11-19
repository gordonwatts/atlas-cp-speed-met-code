#!/usr/bin/env bash
set -eu
repo=$(cd "$(dirname "$0")/.." && pwd)
build=$1
out=$2
events=$3
mode=$4
mkdir -- "$out"
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
set -eu
unset LD_PRELOAD LD_DEBUG LD_DEBUG_OUTPUT
repo=$1
build=$2
out=$3
events=$4
mode=$5
export ATLAS_MET_FRAMEWORK_BRIDGE="$build/libAtlasMetFramework.so"
export ATLAS_MET_PROVIDER="$build/maker_adapter/libatlas_met_maker_provider.so"
export ATLAS_MET_MAKER_PAYLOAD="$repo/results/maker-component-001/maker.json"
export ATLAS_MET_SOFT_PAYLOAD="$repo/results/soft-component-exact-001/bundle/corrections.json"
export ATLAS_MET_SIGNIFICANCE_PAYLOAD="$repo/results/significance-component-exact-001/bundle/corrections.json"
export ATLAS_MET_STAGE=full ATLAS_MET_DISPATCH_AUDIT=1
case "$mode" in
  trace) export ATLAS_MET_MAKER_TRACE_FILE="$out/maker-trace.tsv";;
  control) unset ATLAS_MET_MAKER_TRACE_FILE;;
  *) exit 2;;
esac
input=$(cat "$repo/config/input.txt")
cd "$out"
timeout --signal=TERM --kill-after=5s 90s python "$repo/components/framework/bridge_entry.py" --direct-driver -i "$input" -t "$repo/config/reference-smoke.yaml" -o candidate -e "$events" --dump-full-config
' bash "$repo" "$build" "$out" "$events" "$mode" > "$out/job.log" 2>&1
