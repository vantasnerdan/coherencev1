#!/usr/bin/env bash
# Caller supplies the systemd cap and kill-safe mpi-unit-finish.sh callback.
# Usage: NIH_PHASE=UNIQUE GPAW_PYTHON=python bash run-nih.sh zero|probe RANKS
set -euo pipefail
phase=${NIH_PHASE:?Required unique NIH_PHASE}
stage=${1:?Required zero or probe stage}
ranks=${2:?Required rank count}
[[ "$phase" =~ ^[a-zA-Z0-9_.-]+$ && "$phase" != . && "$phase" != .. ]]
[[ "$stage" == zero || "$stage" == probe ]]
[[ "$ranks" =~ ^[1-9][0-9]*$ ]]
cd -- "$(dirname -- "$0")"
log="$PWD/run-$phase.log"
( set -o noclobber; : >"$log" )
finish() {
 status=$?
 trap - EXIT
 printf 'NIH_ALL_EXIT phase=%s status=%s\n' "$phase" "$status" | tee -a "$log"
 exit "$status"
}
trap finish EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
export NIH_RUN_ID="$phase"
export GPAW_SETUP_PATH="$PWD/setups" GPAW_MPI_BACKEND=cgpaw
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python=${GPAW_PYTHON:-python}
command=(mpiexec --bind-to none -n "$ranks" "$python" -u -m gpaw python -- nih_dense.py "$stage" --ranks "$ranks")
{
 printf 'UTC_START '; date -u +%FT%TZ
 printf 'AFFINITY '; taskset -pc $$
 printf 'ARGV '; printf '%q ' "${command[@]}"; printf '\n'
} | tee -a "$log"
( set -o noclobber; printf '%q ' "${command[@]}" >"$phase.argv" )
printf '\n' >>"$phase.argv"
/usr/bin/time -v -o "$phase.time" "${command[@]}" 2>&1 | tee -a "$log"
