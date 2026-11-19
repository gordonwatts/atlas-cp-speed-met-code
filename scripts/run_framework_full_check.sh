#!/usr/bin/env bash
set -eu
repo=$(cd "$(dirname "$0")/.." && pwd)
build=$1
out=$2
events=$3
mode=$4
scope=${5:-rich-full}
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
scope=$6
export ATLAS_MET_FRAMEWORK_BRIDGE="$build/libAtlasMetFramework.so"
export ATLAS_MET_PROVIDER="$build/maker_adapter/libatlas_met_maker_provider.so"
export ATLAS_MET_MAKER_PAYLOAD="$repo/results/maker-component-001/maker.json"
export ATLAS_MET_SOFT_PAYLOAD="$repo/results/soft-component-exact-001/bundle/corrections.json"
export ATLAS_MET_SIGNIFICANCE_PAYLOAD="$repo/results/significance-event-payload-001/met_significance_event.json"
export ATLAS_MET_STAGE=full ATLAS_MET_DISPATCH_AUDIT=1
export ATLAS_MET_FRAMEWORK_PYTHON="$repo/components/framework/bridge_entry.py"
if [[ -f "$build/../source/components/framework/bridge_entry.py" ]]; then
  ATLAS_MET_FRAMEWORK_PYTHON="$build/../source/components/framework/bridge_entry.py"
fi
case "$mode" in
  trace) export ATLAS_MET_MAKER_TRACE_FILE="$out/maker-trace.tsv";;
  control) unset ATLAS_MET_MAKER_TRACE_FILE;;
  *) exit 2;;
esac
input=$(cat "$repo/config/input.txt")
cd "$out"
limit="$((60 + events / 100))s"
flags=()
case "$scope" in
  rich-full) config=reference-smoke.yaml;;
  rich-nominal) config=reference-smoke.yaml; flags=(--no-systematics);;
  production-full) config=production.yaml;;
  production-nominal) config=production.yaml; flags=(--no-systematics);;
  rich-corr)
    export PYTHONPATH="$repo/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
    exec timeout --signal=TERM --kill-after=5s "$limit" python "$repo/scripts/run_framework_corr_check.py" "$input" "$repo/config/reference-smoke.yaml" "$events";;
  *) exit 2;;
esac
exec timeout --signal=TERM --kill-after=5s "$limit" python "$ATLAS_MET_FRAMEWORK_PYTHON" --direct-driver -i "$input" -t "$repo/config/$config" -o candidate -e "$events" --dump-full-config "${flags[@]}"
' bash "$repo" "$build" "$out" "$events" "$mode" "$scope" > "$out/job.log" 2>&1
