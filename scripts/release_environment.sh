#!/usr/bin/env bash
# Source from a fresh process with no release or replacement-tool overlay.
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet || return
asetup AnalysisBase,25.2.117 || return
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
export CMAKE_BUILD_PARALLEL_LEVEL=1 RAYON_NUM_THREADS=1
