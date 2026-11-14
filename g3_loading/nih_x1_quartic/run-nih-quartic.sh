#!/usr/bin/env bash
# One direction per bounded Native unit; caller owns cost admission and callback.
set -euo pipefail
phase=${NIH_PHASE:?Required unique NIH_PHASE}
direction=${1:?Required 100 or 111 direction}
ranks=${2:?Required rank count}
[[ "$phase" =~ ^[a-zA-Z0-9_.-]+$ && "$phase" != . && "$phase" != .. ]]
[[ "$direction" == 100 || "$direction" == 111 ]]
[[ "$ranks" =~ ^[1-9][0-9]*$ ]]
cd -- "$(dirname -- "$0")"
log="$PWD/run-$phase.log"
( set -o noclobber; : >"$log" )
finish() {
 status=$?
 trap - EXIT
 printf 'NIH_QUARTIC_ALL_EXIT phase=%s status=%s\n' "$phase" "$status" | tee -a "$log"
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
command=(mpiexec --bind-to none -n "$ranks" "$python" -u -m gpaw python -- nih_quartic.py "$direction" --ranks "$ranks")
{
 printf 'UTC_START '; date -u +%FT%TZ
 printf 'AFFINITY '; taskset -pc $$
 printf 'ARGV '; printf '%q ' "${command[@]}"; printf '\n'
} | tee -a "$log"
( set -o noclobber; printf '%q ' "${command[@]}" >"$phase.argv" )
printf '\n' >>"$phase.argv"
/usr/bin/time -v -o "$phase.time" "${command[@]}" 2>&1 | tee -a "$log"
