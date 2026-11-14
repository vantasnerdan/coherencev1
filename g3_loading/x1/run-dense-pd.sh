#!/usr/bin/env bash
# One explicitly admitted stage only; no resource admission or mesh adoption.
# Usage: bash run-dense-pd.sh probe --mesh 3 --a 4.084 [--ranks 6]
# For MPI use the installed GPAW CLI under mpiexec, never plain python replicas.
# Caller supplies Native/systemd caps and arms all-exit subscription before launch.
# To cover SIGKILL/pre-shell failure, ALSO set systemd ExecStopPost to:
# /bin/bash "$PWD/unit-finish.sh" "$DENSE_PD_PHASE"
# unit-finish.sh appends its established marker to run-$DENSE_PD_PHASE.log.
# GNU time is diagnostic only for MPI; charge Native/cgroup CPUUsageNSec instead.
set -euo pipefail
phase="${DENSE_PD_PHASE:-dense-pd-$(date -u +%Y%m%dT%H%M%S)-$$}"
log=""
finish() {
    status=$?
    trap - EXIT
    printf 'DENSE_PD_ALL_EXIT phase=%s status=%s\n' "$phase" "$status"
    if [[ -n "$log" ]]; then
        printf 'DENSE_PD_ALL_EXIT phase=%s status=%s\n' "$phase" "$status" >>"$log"
    fi
    exit "$status"
}
trap finish EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
[[ "$phase" =~ ^[a-zA-Z0-9_.-]+$ ]] || { printf 'Unsafe DENSE_PD_PHASE\n' >&2; exit 2; }
cd -- "$(dirname -- "$0")"
# Refuse a reused phase without appending to any previous evidence.
[[ ! -e "$PWD/$phase.argv" && ! -e "$PWD/$phase.time" ]] || { printf 'Reused phase\n' >&2; exit 2; }
candidate_log="$PWD/run-$phase.log"
( set -o noclobber; : >"$candidate_log" )
log="$candidate_log"
export DENSE_PD_RUN_ID="$phase"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export GPAW_SETUP_PATH="$PWD/setups"
python=${GPAW_PYTHON:-python}
ranks=1
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
    if [[ "${args[i]}" == --ranks ]]; then
        ((i+1 < ${#args[@]})) || { printf 'Missing --ranks value\n' >&2; exit 2; }
        ranks="${args[i+1]}"
    elif [[ "${args[i]}" == --ranks=* ]]; then
        ranks="${args[i]#--ranks=}"
    fi
done
[[ "$ranks" =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid --ranks\n' >&2; exit 2; }
if ((ranks == 1)); then
    export GPAW_MPI_BACKEND=serial
    command=("$python" -u dense_pd.py "${args[@]}")
else
    export GPAW_MPI_BACKEND=cgpaw
    command=(mpiexec --bind-to none -n "$ranks" "$python" -u -m gpaw python -- dense_pd.py "${args[@]}")
fi
{
    printf 'UTC_START '; date -u +%FT%TZ
    printf 'AFFINITY '; taskset -pc $$
    printf 'ARGV '; printf '%q ' "${command[@]}"; printf '\n'
    printf 'OUTPUT_LOG %s\n' "$log"
} | tee -a "$log"
printf '%q ' "${command[@]}" >"$PWD/$phase.argv"
printf '\n' >>"$PWD/$phase.argv"
/usr/bin/time -v -o "$PWD/$phase.time" "${command[@]}" 2>&1 | tee -a "$log"
