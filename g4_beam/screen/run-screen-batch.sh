#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
phase=$1
plan=$2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1
exec "$GPAW_PYTHON" -u "$(dirname -- "$plan")/source/run_bulk_batch.py" \
 --plan "$plan" --setup-root "$PWD/setups" >"run-$phase.log" 2>&1
