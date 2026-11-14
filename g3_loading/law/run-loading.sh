#!/usr/bin/env bash
# One admitted zero|probe stage; each probe's axis is registered in its config.
# Usage: bash run-loading.sh zero --config loading-configs/CONFIG.json --mesh 3 --ranks N
# Displaced stages additionally require --probe clean|near; zero is config-global.
# Caller owns Native/systemd caps and all-exit subscription. Set ExecStopPost to
# /bin/bash "$PWD/mpi-unit-finish.sh" "$LOADING_PD_PHASE" for kill-safe accounting.
# Charge cgroup CPUUsageNSec, not GNU time or launcher/rank-only rusage.
set -euo pipefail
phase="${LOADING_PD_PHASE:-loading-pd-$(date -u +%Y%m%dT%H%M%S)-$$}"
log=""
finish() {
    status=$?
    trap - EXIT
    printf 'LOADING_PD_ALL_EXIT phase=%s status=%s\n' "$phase" "$status"
    if [[ -n "$log" ]]; then
        printf 'LOADING_PD_ALL_EXIT phase=%s status=%s\n' "$phase" "$status" >>"$log"
    fi
    exit "$status"
}
trap finish EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
[[ "$phase" =~ ^[a-zA-Z0-9_.-]+$ && "$phase" != . && "$phase" != .. ]] || { printf 'Unsafe LOADING_PD_PHASE\n' >&2; exit 2; }
cd -- "$(dirname -- "$0")"
[[ ! -e "$PWD/$phase.argv" && ! -e "$PWD/$phase.time" ]] || { printf 'Reused phase\n' >&2; exit 2; }
candidate_log="$PWD/run-$phase.log"
( set -o noclobber; : >"$candidate_log" )
log="$candidate_log"
export LOADING_PD_RUN_ID="$phase"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export GPAW_SETUP_PATH="$PWD/setups"
python=${GPAW_PYTHON:-python}
ranks=""
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
    if [[ "${args[i]}" == --ranks ]]; then
        ((i+1 < ${#args[@]})) || { printf 'Missing --ranks value\n' >&2; exit 2; }
        ranks="${args[i+1]}"
    elif [[ "${args[i]}" == --ranks=* ]]; then
        ranks="${args[i]#--ranks=}"
    fi
done
[[ "$ranks" =~ ^[1-9][0-9]*$ ]] || { printf 'Required positive --ranks\n' >&2; exit 2; }
if ((ranks == 1)); then
    export GPAW_MPI_BACKEND=serial
    command=("$python" -u loading_pd.py "${args[@]}")
else
    export GPAW_MPI_BACKEND=cgpaw
    command=(mpiexec --bind-to none -n "$ranks" "$python" -u -m gpaw python -- loading_pd.py "${args[@]}")
fi
{
    printf 'UTC_START '; date -u +%FT%TZ
    printf 'AFFINITY '; taskset -pc $$
    printf 'ARGV '; printf '%q ' "${command[@]}"; printf '\n'
    printf 'OUTPUT_LOG %s\n' "$log"
} | tee -a "$log"
( set -o noclobber; printf '%q ' "${command[@]}" >"$PWD/$phase.argv" )
printf '\n' >>"$PWD/$phase.argv"
/usr/bin/time -v -o "$PWD/$phase.time" "${command[@]}" 2>&1 | tee -a "$log"
