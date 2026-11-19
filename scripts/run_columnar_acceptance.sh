#!/usr/bin/env bash
# Compare direct-original execution against independent release evidence.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
out=$1
snapshot=$2
histpy=$3
mkdir -- "$out"
cd "$repo"
source="$snapshot/source/components/columnar"
runtime="$snapshot/runtime"
cp "$snapshot/source-head.txt" "$out/source-head.txt"
input=$(cat config/input.txt)
clean=(env -i "HOME=$HOME" PATH=/usr/bin:/bin LANG=C.UTF-8 "LD_LIBRARY_PATH=$runtime/native"
       OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
       VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 RAYON_NUM_THREADS=1)
for scope in full nominal corr; do
  dest="$out/$scope"
  mkdir "$dest"
  mode=()
  trace=()
  stock="$repo/results/reference-observer-$scope-001"
  stock_hist="$repo/results/reference-production-$scope-002/histogram.npz"
  if [ "$scope" = nominal ]; then
    mode=(--nominal-only)
  else
    mode=(--trace "$dest/c-status.tsv")
    trace=(--stock-trace "$repo/results/reference-status-$scope-001/status.tsv" --candidate-trace "$dest/c-status.tsv"
           --stock-getters "$repo/results/reference-prepared-full-001/prepared-p4.tsv")
  fi
  if [ "$scope" = corr ]; then
    mode+=(--correlated-diagnostic)
    stock_hist="$dest/stock.npz"
    "${clean[@]}" taskset -c 0 "$histpy" scripts/histogram_reference.py "$stock/reference.root" "$stock/hist-reference.root" "$stock_hist" > "$dest/compare.log" 2>&1
  fi
  "${clean[@]}" taskset -c 0 "$histpy" "$source/met_demo.py" "$input" "$dest/c.npz" --runtime "$runtime" \
    --stop 20000 --observables "$dest/c-observables.npz" --runtime-record "$dest/c-runtime.json" "${mode[@]}" > "$dest/job.log" 2>&1
  status=0
  "${clean[@]}" "PYTHONPATH=$repo/scripts" taskset -c 0 "$histpy" "$source/compare.py" "$stock" "$stock_hist" \
    "$dest/c-observables.npz" "$dest/c.npz" "$dest/compare.json" "${trace[@]}" --expected-events 20000 >> "$dest/compare.log" 2>&1 || status=$?
  printf '%s %d\n' "$scope" "$status" >> "$out/status.txt"
  if [ "$status" -ne 0 ]; then tail -35 "$dest/compare.log"; exit "$status"; fi
  "${clean[@]}" taskset -c 0 "$histpy" scripts/compare_columnar.py "$stock" "$dest" 20000 >> "$dest/compare.log" 2>&1
  cmp "$dest/c.npz" "$stock_hist"
done
