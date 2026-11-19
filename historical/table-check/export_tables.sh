#!/usr/bin/env bash
set +e
source /home/gwatts/atlas/correctionlib_tests/muon/setup-reference.sh
rc=$?
[ "$rc" = 0 ] || exit "$rc"
set -e
export BLIS_NUM_THREADS=1
cd /home/gwatts/atlas/correctionlib_tests
exec taskset -c 0 python met/scripts/export_tables.py "$@"
