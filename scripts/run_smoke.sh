#!/usr/bin/env bash
set -u
if [[ $# < 2 || $# > 3 ]]; then
  echo "Usage: $0 NEW_RESULT_DIRECTORY EVENT_COUNT [full|nominal]" >&2
  exit 2
fi
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$1
events=$2
systematics=()
case "${3:-full}" in
  full) ;;
  nominal) systematics=(--no-systematics);;
  *) exit 2;;
esac
if [[ $work != /* || ! $events =~ ^[1-9][0-9]*$ ]]; then
  echo 'Use an absolute result path and a positive event count' >&2
  exit 2
fi
mkdir -- "$work" || exit
read -r input < "$repo/config/input.txt"
cd "$work" || exit
cp -- "$repo/config/reference-smoke.yaml" reference-smoke.yaml || exit
cat > command.txt <<CMD
CPU 0. Stock AnalysisBase 25.2.117.
CPRun.py --direct-driver -i $input -t $work/reference-smoke.yaml -o reference -e $events --dump-full-config ${systematics[*]}
CMD
taskset -c 0 bash -c '
source "$1" || exit
export ROOT_MAX_THREADS=1
if [[ -n ${MET_TRACE_LIBRARY:-} ]]; then
  export LD_PRELOAD="$MET_TRACE_LIBRARY"
  export MET_STOCK_TRACE="$PWD/status.tsv"
fi
shift
exec python "$AnalysisBase_DIR/bin/CPRun.py" "$@"
' bash "$repo/scripts/release_environment.sh" --direct-driver -i "$input" -t "$work/reference-smoke.yaml" -o reference -e "$events" --dump-full-config "${systematics[@]}" > job.log 2>&1 &
job_pid=$!
while kill -0 "$job_pid" 2>/dev/null; do
  executable=$(readlink "/proc/$job_pid/exe" 2>/dev/null || :)
  if [[ $executable == /cvmfs/*/bin/python* ]]; then
    {
      printf 'pid=%s exe=%s ' "$job_pid" "$executable"
      awk '/^(Threads|Cpus_allowed_list):/ {printf "%s%s ", $1, $2}' "/proc/$job_pid/status"
      printf 'children='
      cat "/proc/$job_pid/task/$job_pid/children"
      echo
    } >> process-observation.txt 2>/dev/null
    if [[ -r /proc/$job_pid/maps ]]; then
      { cat loaded-met-libraries.txt 2>/dev/null || :; awk '$6 ~ /[Ll]ib(METUtil|xAODMissingET|MetAnalysis|JetCalib|TauAnalysis)/ {print $6}' "/proc/$job_pid/maps" 2>/dev/null; } | sort -u > .libs-current
      mv -- .libs-current loaded-met-libraries.txt
    fi
  fi
  sleep 0.2
done
wait "$job_pid"
