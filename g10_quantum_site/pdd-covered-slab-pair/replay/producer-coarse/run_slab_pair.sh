#!/usr/bin/env bash
# A single bounded native lane disjoint from TiD2 (6-10).
set -euo pipefail
root=${G10_ROOT:?Registered G10 root required}
stop=${G10_ABSOLUTE_STOP:?UTC stop required}
python=${G10_GPAW_PYTHON:?Registered GPAW executable required}
well=${G10_WELL_SOURCE:?Pinned shared code required}
source=$(dirname -- "${BASH_SOURCE[0]}")
output="$root/slab-pair/native"
ref="$root/sites/surface-relax-002/reference.json"
cutoff=$(date -ud "$stop" +%s)

point() {
  local mesh=$1 role=$2 amp=$3 ranks=4 cores=0-3 memory=6G min_gib=12 cap=3300 rc=0 free left label base
  if ((mesh == 3)); then ranks=3; cores=0-2; memory=8G; min_gib=14; cap=4800; fi
  label="m$mesh-${role,,}-$amp"
  base="$root/receipts/slab-pair-$label"
  free=$(awk '$1=="MemAvailable:" {print $2; exit}' /proc/meminfo)
  left=$((cutoff - $(date -u +%s) - 600))
  if ((free < min_gib*1024*1024 || left < 300)); then
    printf 'G10_SLAB_PAIR_NOT_DISPATCHED %s available_kib=%s remaining_s=%s\n' "$label" "$free" "$left"
    return 2
  fi
  if ((left < cap)); then cap=$left; fi
  local cmd=(systemd-run --user --wait --unit="g10-slab-pair-$label" --description='G10 fixed covered-slab pair'
    -p Type=exec -p Restart=no -p "RuntimeMaxSec=$cap" -p TimeoutStopSec=35 -p "MemoryMax=$memory" -p MemoryAccounting=yes
    --working-directory="$root" --setenv="PYTHONPATH=$source:$well" --setenv=PYTHONDONTWRITEBYTECODE=1
    --setenv=GPAW_MPI_BACKEND=cgpaw --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 --setenv=MKL_NUM_THREADS=1
    /usr/bin/taskset -c "$cores" /usr/bin/mpiexec --bind-to none -n "$ranks" /usr/bin/taskset -c "$cores"
    "$python" -u -m gpaw python -- "$source/slab_pair.py" --reference "$ref" point "$role" "$amp"
    --mesh "$mesh" --ranks "$ranks" --output "$output/$label.json")
  (set -o noclobber; printf '%s\n' "${cmd[@]}" > "$base.argv") || return 2
  date -u +%FT%TZ > "$base.started"
  printf 'G10_SLAB_PAIR_DISPATCH %s\n' "$label"
  "${cmd[@]}" > "$base.stdout" 2> "$base.stderr" || rc=$?
  printf '%s\n' "$rc" > "$base.exit"
  date -u +%FT%TZ > "$base.ended"
  if ((rc)) || [[ ! -s "$output/$label.json" ]]; then
    cat "$base.stderr" >&2
    printf 'G10_SLAB_PAIR_POINT_EXIT %s rc=%s\n' "$label" "$rc"
    return 1
  fi
  printf 'G10_SLAB_PAIR_POINT_DONE %s\n' "$label"
}

point 2 zero 0 || exit 1
for amp in 6 12 24 36 48; do
  for role in A B AB; do point 2 "$role" "$amp" || exit 1; done
  if ((amp == 36)); then
    if point 3 zero 0; then point 3 AB 36 || :; fi
  fi
done
printf 'G10_SLAB_PAIR_LANE_SETTLED mesh_complete=%s\n' "$([[ -s "$output/m3-ab-36.json" ]] && printf yes || printf no)"
