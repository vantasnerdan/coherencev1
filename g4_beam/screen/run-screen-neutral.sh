#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
phase=$1
plan=$2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1 GPAW_SETUP_PATH="$PWD/setups"
exec "$GPAW_PYTHON" -u "$(dirname -- "$plan")/source/run_neutral_cycle.py" \
 --input "$plan" --output "$(dirname -- "$plan")" >"run-$phase.log" 2>&1
