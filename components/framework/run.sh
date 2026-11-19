#!/usr/bin/env bash
# The runtime owner provides the release environment before this wrapper.
set -euo pipefail
unset LD_PRELOAD LD_DEBUG LD_DEBUG_OUTPUT
framework=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd -- "$framework/../.." && pwd)
export ATLAS_MET_FRAMEWORK_BRIDGE
ATLAS_MET_FRAMEWORK_BRIDGE=$(realpath "${1:?FRAMEWORK_SO}")
out=${2:?NEW_OUTPUT_DIRECTORY}
mode=${3:?nominal or full_systematics}
entries=${4:-16}
case "$mode" in
    nominal) flags=(--no-systematics) ;;
    full_systematics) flags=() ;;
    *) printf '%s\n' "Unknown run scope: $mode" >&2; exit 2 ;;
esac
if [[ -n ${ATLAS_MET_DISPATCH_AUDIT:-} ]]; then
    flags+=(--dump-full-config)
fi
read -r input < "$repo/config/input.txt"
mkdir "$out"
cd "$out"
exec python "$framework/bridge_entry.py" --direct-driver -i "$input" -t "$repo/config/production.yaml" -o candidate -e "$entries" "${flags[@]}"
