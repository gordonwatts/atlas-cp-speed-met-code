#!/usr/bin/env bash
set -u
repo=$(cd "$(dirname "$0")/.." && pwd)
build=$1
mkdir -- "$build" || exit
exec taskset -c 0 bash -c '
source "$1/scripts/release_environment.sh" || exit
read -ra root_flags <<< "$(root-config --cflags)"
flags=(
  -std=c++23 -O2 -fPIC -shared -DXAOD_STANDALONE -DXAOD_ANALYSIS
  -DCOLUMNAR_MODE_SUFFIX= -DCOLUMNAR_MODE_NAME=ColumnarModeXAOD
  -I "$AnalysisBase_DIR/include" -I "$(root-config --incdir)/eigen3"
)
libraries=(
  -L "$AnalysisBase_DIR/lib" -Wl,-z,defs
  -lxAODMissingET -lxAODBase -lxAODJet -lxAODRootAccess -lxAODRootAccessInterfaces -lAthLinks -lSystematicsHandlesLib -lSelectionHelpersLib -lAsgDataHandlesLib -lPATInterfaces -lAthContainers -lCxxUtils -L "$(root-config --libdir)" -lCore -ldl
)
exec g++ "${root_flags[@]}" "${flags[@]}" "$1/observer/status_trace.cxx" "${libraries[@]}" -o "$2/libmet_stock_trace.so"
' bash "$repo" "$build" > "$build/build.log" 2>&1
