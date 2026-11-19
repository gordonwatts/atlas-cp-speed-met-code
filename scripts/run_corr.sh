#!/usr/bin/env bash
set -u
if [[ $# != 2 ]]; then
  echo "Usage: $0 NEW_RESULT_DIRECTORY EVENT_COUNT" >&2
  exit 2
fi
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$1
events=$2
if [[ $work != /* || ! $events =~ ^[1-9][0-9]*$ ]]; then
  exit 2
fi
mkdir -- "$work" || exit
read -r input < "$repo/config/input.txt"
cp -- "$repo/config/reference-smoke.yaml" "$work/reference-smoke.yaml" || exit
cat > "$work/command.txt" <<CMD
taskset -c 0 python $repo/scripts/run_corr_diagnostic.py $input $work/reference-smoke.yaml --events $events
CMD
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
export PYTHONPATH="$1/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
export ROOT_MAX_THREADS=1
cd "$2" || exit
if [[ -n ${MET_TRACE_LIBRARY:-} ]]; then
  export LD_PRELOAD="$MET_TRACE_LIBRARY"
  export MET_STOCK_TRACE="$PWD/status.tsv"
fi
exec python "$1/scripts/run_corr_diagnostic.py" "$3" "$2/reference-smoke.yaml" --events "$4"
' bash "$repo" "$work" "$input" "$events" > "$work/job.log" 2>&1
