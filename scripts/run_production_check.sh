#!/usr/bin/env bash
set -u
if [[ $# != 3 ]]; then
  echo "Usage: $0 NEW_RESULT_DIRECTORY EVENT_COUNT full|nominal|stored" >&2
  exit 2
fi
repo=$(cd "$(dirname "$0")/.." && pwd)
out=$1
events=$2
if [[ $out != /* || ! $events =~ ^[1-9][0-9]*$ ]]; then
  exit 2
fi
config=production.yaml
systematics=()
case "$3" in
  full) ;;
  nominal) systematics=(--no-systematics);;
  stored) config=stored-core.yaml; systematics=(--no-systematics);;
  *) exit 2;;
esac
mkdir -- "$out" || exit
read -r input < "$repo/config/input.txt"
cp -- "$repo/config/$config" "$out/config.yaml" || exit
cd "$out" || exit
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
export ROOT_MAX_THREADS=1
shift
exec CPRun.py "$@"
' bash "$repo" --direct-driver -i "$input" -t "$out/config.yaml" -o reference -e "$events" "${systematics[@]}" > job.log 2>&1
