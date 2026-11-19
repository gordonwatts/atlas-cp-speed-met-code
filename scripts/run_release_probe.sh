#!/usr/bin/env bash
set -u
repo=$(cd "$(dirname "$0")/.." && pwd)
read -r input < "$repo/config/input.txt"
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
export PYTHONPATH="$1/.venv/release-cli313${PYTHONPATH:+:$PYTHONPATH}"
export ROOT_MAX_THREADS=1
repo=$1
shift
exec python "$repo/scripts/inspect_release.py" "$@"
' bash "$repo" "$input" "$@"
