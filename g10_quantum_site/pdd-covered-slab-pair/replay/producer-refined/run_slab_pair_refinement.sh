#!/usr/bin/env bash
# Separate native lanes; reuse the sealed m3 zero and AB-36.
set -euo pipefail
root=${G10_ROOT:?Registered G10 root required}
stop=${G10_ABSOLUTE_STOP:?UTC stop required}
python=${G10_GPAW_PYTHON:?Registered GPAW executable required}
well=${G10_WELL_SOURCE:?Pinned shared code required}
source=$(dirname -- "${BASH_SOURCE[0]}")
output="$root/slab-pair/native"
ref="$root/sites/surface-relax-002/reference.json"
cutoff=$(date -ud "$stop" +%s)
lane=${1:?Use lane curve or control}

point() {
  local mesh=$1 role=$2 amp=$3 ranks=3 cores=0-2 memory=8G min_gib=10 cap=3600 rc=0 free left label base
  if [[ "$lane" == control ]]; then cores=3-5; fi
  if ((mesh == 4)); then ranks=4; cores=3-6; memory=10G; min_gib=12; cap=4500; fi
  label="m$mesh-${role,,}-$amp"
  base="$root/receipts/slab-pair-$label"
  free=$(awk '$1=="MemAvailable:" {print $2; exit}' /proc/meminfo)
  left=$((cutoff - $(date -u +%s) - 600))
  if ((free < min_gib*1024*1024 || left < 300)); then
    printf 'G10_SLAB_REFINED_NOT_DISPATCHED %s available_kib=%s remaining_s=%s\n' "$label" "$free" "$left"
    return 2
  fi
  if ((left < cap)); then cap=$left; fi
  local cmd=(systemd-run --user --wait --unit="g10-slab-pair-$label" --description='G10 covered slab k-mesh refinement'
    -p Type=exec -p Restart=no -p "RuntimeMaxSec=$cap" -p TimeoutStopSec=35 -p "MemoryMax=$memory" -p MemorySwapMax=0 -p MemoryAccounting=yes
    --working-directory="$root" --setenv="PYTHONPATH=$source:$well" --setenv=PYTHONDONTWRITEBYTECODE=1
    --setenv=GPAW_MPI_BACKEND=cgpaw --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 --setenv=MKL_NUM_THREADS=1
    /usr/bin/taskset -c "$cores" /usr/bin/mpiexec --bind-to none -n "$ranks" /usr/bin/taskset -c "$cores"
    "$python" -u -m gpaw python -- "$source/slab_pair.py" --reference "$ref" point "$role" "$amp"
    --mesh "$mesh" --ranks "$ranks" --output "$output/$label.json")
  (set -o noclobber; printf '%s\n' "${cmd[@]}" > "$base.argv") || return 2
  date -u +%FT%TZ > "$base.started"
  printf 'G10_SLAB_REFINED_DISPATCH %s available_kib=%s\n' "$label" "$free"
  "${cmd[@]}" > "$base.stdout" 2> "$base.stderr" || rc=$?
  printf '%s\n' "$rc" > "$base.exit"
  date -u +%FT%TZ > "$base.ended"
  if ((rc)) || [[ ! -s "$output/$label.json" ]]; then
    cat "$base.stderr" >&2
    printf 'G10_SLAB_REFINED_POINT_EXIT %s rc=%s\n' "$label" "$rc"
    return 1
  fi
  printf 'G10_SLAB_REFINED_POINT_DONE %s\n' "$label"
}

case "$lane" in
  curve)
    point 3 A 36; point 3 B 36
    for amp in 6 12; do
      for role in A B AB; do point 3 "$role" "$amp"; done
    done
    point 3 AB 48
    ;;
  control)
    # Complete the remaining main-curve samples even if the single control fails.
    point 3 A 24; point 3 B 24; point 3 AB 24; point 3 A 48; point 3 B 48
    if systemctl --user --quiet is-active g10-tid2-case-110m120.service g10-tid2-case-110m240.service ||
       [[ ! -s "$root/tid2/native/110m240/point.json" ]]; then
      printf 'G10_SLAB_REFINED_NOT_DISPATCHED m4-zero-0 reason=TiD2-core-6-not-released\n'
      exit 2
    fi
    point 4 zero 0
    for role in A B AB; do point 4 "$role" 36; done
    ;;
  *) printf 'Unregistered lane: %s\n' "$lane" >&2; exit 2 ;;
esac
