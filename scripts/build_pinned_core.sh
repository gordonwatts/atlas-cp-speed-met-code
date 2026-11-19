#!/usr/bin/env bash
set -eo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
out=$1
mkdir -- "$out"
exec taskset -c 0 bash -c '
set -o pipefail
repo=$1
out=$2
source "$repo/scripts/release_environment.sh" || exit
set -eu
vendor=$repo/vendor/correctionlib
pin=3dcfa20ec01fb40c4d5fec58383785519dec1c00
[[ "$(git -C "$vendor" rev-parse HEAD)" == "$pin" ]]
compiler=$(command -v g++)
printf "core origin %s\ncore pin %s\ncompiler %s\n"   "$(git -C "$vendor" remote get-url origin)" "$pin" "$compiler" > "$out/source.txt"
git -C "$vendor" ls-tree HEAD | awk '"'"'$2 == "commit" {print $3, $4}'"'"' > "$out/submodule-lock.txt"
while read -r revision name; do
  destination="$vendor/$name"
  if [[ ! -d $destination/.git && ! -f $destination/.git ]]; then
    origin=$(git -C "$vendor" config -f .gitmodules --get "submodule.$name.url")
    git -C "$destination" init --quiet
    git -C "$destination" remote add origin "$origin"
    found=
    for prior in         /home/gwatts/atlas-cp-speed/jets-skill-c8a2374-20261117/source/correctionlib         /home/gwatts/atlas-cp-speed/muon-skill-c8a2374-20261117/correctionlib; do
      if [[ -e $prior/$name/.git ]] && git -C "$prior/$name" cat-file -e "$revision^{commit}" 2>/dev/null; then
        found=$prior/$name
        break
      fi
    done
    fetch_source=${found:-$origin}
    printf "%s %s %s\n" "$name" "$revision" "$fetch_source" >> "$out/submodule-resolved.txt"
    GIT_TERMINAL_PROMPT=0 git -C "$destination" fetch --depth 1 "$fetch_source" "$revision"
    git -C "$destination" checkout --quiet --detach FETCH_HEAD
  fi
  [[ "$(git -C "$destination" rev-parse HEAD)" == "$revision" ]]
  [[ -z "$(git -C "$destination" status --porcelain --untracked-files=no)" ]]
done < "$out/submodule-lock.txt"
cmake -S "$vendor" -B "$out/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$out/install"   -DCMAKE_CXX_COMPILER="$compiler" -DCORRECTIONLIB_BUILD_PYTHON=OFF -DCORRECTIONLIB_BUILD_TESTS=OFF   -DBUILD_DEMO=OFF > "$out/configure.log" 2>&1
cmake --build "$out/build" --parallel 1 > "$out/build.log" 2>&1
cmake --install "$out/build" > "$out/install.log" 2>&1
find "$out/install" -name "*orrectionlib*Config.cmake" -o -name libcorrectionlib.so
sha256sum "$out/install/correctionlib/lib/libcorrectionlib.so" "$vendor/include/correction.h" "$vendor/include/correction_plugin.h"
' bash "$repo" "$out" > "$out/job.log" 2>&1
