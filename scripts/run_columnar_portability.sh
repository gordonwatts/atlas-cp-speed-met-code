#!/usr/bin/env bash
# Exercise the actual local wheel from two fresh, release-independent installations.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
out=$1
snapshot=$2
reference_python=$3
uv=$4
mkdir -- "$out"
cd "$repo"
python314=$(realpath -- "$reference_python")
source="$snapshot/source/components/columnar"
runtime="$snapshot/runtime"
input=$(cat config/input.txt)
clean=(env -i "HOME=$HOME" PATH=/usr/bin:/bin LANG=C.UTF-8
       OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
       VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 RAYON_NUM_THREADS=1)
"${clean[@]}" "$reference_python" -c 'from importlib.metadata import distributions; print("\n".join(sorted(d.metadata["Name"]+"=="+d.version for d in distributions())))' > "$out/constraints.txt"
"${clean[@]}" "$uv" --version > "$out/uv-version.txt"
"${clean[@]}" "$python314" --version > "$out/python-version.txt"
sha256sum "$uv" "$python314" > "$out/executables.sha256"
"${clean[@]}" taskset -c 0 "$reference_python" "$source/assemble.py" "$runtime" "$out/bundle" \
  "$repo/vendor/correctionlib-exact-001" "$repo/results/core-exact-build-001/core-local.patch" > "$out/assemble.log" 2>&1
wheel=("$out/bundle"/wheels/*.whl)
test "${#wheel[@]}" -eq 1
"${clean[@]}" taskset -c 0 "$python314" -m venv "$out/pip-env"
pip="$out/pip-env/bin/python"
mkdir "$out/wheelhouse"
"${clean[@]}" taskset -c 0 "$pip" -m pip download --disable-pip-version-check --only-binary=:all: \
  --constraint "$out/constraints.txt" --dest "$out/wheelhouse" "${wheel[0]}" > "$out/wheel-download.log" 2>&1
(
  cd "$out/wheelhouse"
  sha256sum ./*.whl
) > "$out/wheelhouse.sha256"
"${clean[@]}" taskset -c 0 "$pip" -m pip install --disable-pip-version-check --no-index \
  --find-links "$out/wheelhouse" "${wheel[0]}" > "$out/pip-install.log" 2>&1
"${clean[@]}" "$pip" -m pip freeze > "$out/pip-freeze.txt"
visible="$out/bundle/met_demo.py"
uv_run=("${clean[@]}" "UV_CACHE_DIR=$out/uv-cache" taskset -c 0 "$uv"
        run --offline --no-index --find-links "$out/wheelhouse" --no-project --python "$python314")
"${uv_run[@]}" "$visible" "$input" "$out/uv-16.npz" --stop 16 --runtime-record "$out/uv-runtime.json" > "$out/uv-16.log" 2>&1
"${uv_run[@]}" "$visible" "$input" "$out/uv-16-repeat.npz" --stop 16 -n 2 > "$out/uv-16-repeat.log" 2>&1
"${clean[@]}" taskset -c 0 "$pip" "$visible" "$input" "$out/pip-16.npz" --stop 16 --runtime-record "$out/pip-runtime.json" > "$out/pip-16.log" 2>&1
cmp "$out/uv-16.npz" "$out/uv-16-repeat.npz"
cmp "$out/uv-16.npz" "$out/pip-16.npz"
cmp "$out/pip-16.npz" "$repo/results/columnar-full-ledger-smoke-001/stock.npz"
for bad in repeat stop; do
  bad_rc=0
  case "$bad" in
    repeat) flags=(-n 0) ;;
    stop) flags=(--stop -1) ;;
  esac
  "${clean[@]}" taskset -c 0 "$pip" "$visible" "$input" "$out/invalid-$bad.npz" "${flags[@]}" > "$out/invalid-$bad.log" 2>&1 || bad_rc=$?
  printf '%s %d\n' "$bad" "$bad_rc" >> "$out/invalid-exit-codes.txt"
  test "$bad_rc" -eq 2
  test ! -e "$out/invalid-$bad.npz"
done
for scope in full nominal; do
  options=()
  if [ "$scope" = nominal ]; then options=(--nominal-only); fi
  "${clean[@]}" taskset -c 0 "$pip" "$visible" "$input" "$out/pip-$scope.npz" --stop 20000 "${options[@]}" > "$out/pip-$scope.log" 2>&1
  cmp "$out/pip-$scope.npz" "$repo/results/reference-production-$scope-002/histogram.npz"
done
sha256sum "$out"/*.npz > "$out/histogram-sha256.txt"
du -s --apparent-size --block-size=1 "$out/pip-env" "$out/uv-cache" "$out/wheelhouse" "$out/bundle/wheels" > "$out/footprint-bytes.tsv"
"${clean[@]}" "$pip" - "$out" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import sys

out = Path(sys.argv[1])
bundle = json.loads((out / "bundle" / "bundle.json").read_text())
packaged = {key: info for key, info in bundle["inputs"].items()
            if key.startswith(("atlas_met_portable/native/", "atlas_met_portable/payloads/"))}
report = {}
for installer in ("uv", "pip"):
    live = json.loads((out / f"{installer}-runtime.json").read_text())
    libraries = [Path(path) for path in live["mapped_shared_objects"]]
    anchor, = (path for path in libraries if path.name == "libatlas_met_columnar.so")
    package = anchor.parent.parent.parent
    for name, expected in packaged.items():
        data = (package / name).read_bytes()
        assert {"bytes": len(data), "sha256": sha256(data).hexdigest()} == expected, (installer, name)
    report[installer] = {"installation": str(package), "native_and_payload_files_verified": len(packaged),
                         "all_match_wheel_manifest": True}
(out / "installed-identities.json").write_text(json.dumps(report, indent=2) + "\n")
PY
