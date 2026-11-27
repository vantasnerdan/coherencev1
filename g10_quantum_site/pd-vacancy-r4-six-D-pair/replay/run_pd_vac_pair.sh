#!/usr/bin/env bash
# One matched frozen Pd vacancy slot; no concurrent large-cell SCFs.
set -euo pipefail
r=${G10_ROOT:?Registered repository work root}
p=${G10_GPAW_PYTHON:?Pinned Python executable}
epoch=${G10_VAC_PAIR_EPOCH:?Run receipt epoch}
a=$(dirname -- "${BASH_SOURCE[0]}")
ref="$r/pd-multivac/r4-d6-rest-001/reference.json"
out="$r/pd-multivac/pair-native"

point() {
  local mesh=$1 role=$2 amp=$3 n=8 cores=0-7 mem=14G lower=16 wall=3h label q rc=0 free
  if ((mesh==3)); then n=7; cores=0-6; mem=19G; lower=19; wall=5h; fi
  label="m$mesh-${role,,}-$amp"
  q="$r/receipts/vac-pair-$epoch-$label"
  free=$(awk '$1=="MemAvailable:" {print $2; exit}' /proc/meminfo)
  if ((free < lower*1024*1024)); then
    printf 'G10_VAC_PAIR_STOP %s insufficient-memory-headroom\n' "$label"
    return 2
  fi
  local cmd=(systemd-run --user --wait --unit="g10-pd-vac-pair-$epoch-$label" --description='G10 frozen Pd vacancy pair'
    -p Type=exec -p Restart=no -p "RuntimeMaxSec=$wall" -p TimeoutStopSec=35 -p "MemoryMax=$mem" -p MemorySwapMax=0
    --working-directory="$r" --setenv="PYTHONPATH=$a" --setenv=PYTHONDONTWRITEBYTECODE=1 --setenv=GPAW_MPI_BACKEND=cgpaw
    --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 --setenv=MKL_NUM_THREADS=1
    /usr/bin/taskset -c "$cores" /usr/bin/mpiexec --bind-to none -n "$n" /usr/bin/taskset -c "$cores"
    "$p" -u -m gpaw python -- "$a/pd_vac_pair.py" --reference "$ref" --paw "$HOME/well2/setups"
    --articles "$r/sources/pd-multivac-literature-001" point "$role" "$amp" --mesh "$mesh" --ranks "$n" --output "$out/$label.json")
  (set -o noclobber; printf '%s\n' "${cmd[@]}" > "$q.argv")
  date -u +%FT%TZ > "$q.started"
  printf 'G10_VAC_PAIR_DISPATCH %s\n' "$label"
  "${cmd[@]}" > "$q.stdout" 2> "$q.stderr" || rc=$?
  printf '%s\n' "$rc" > "$q.exit"
  date -u +%FT%TZ > "$q.ended"
  if ((rc)) || [[ ! -s "$out/$label.json" ]]; then
    cat "$q.stderr" >&2
    printf 'G10_VAC_PAIR_STOP %s native-exit=%s\n' "$label" "$rc"
    return 1
  fi
  printf 'G10_VAC_PAIR_DONE %s\n' "$label"
}

# Coarse zero is the exact resting-state record on the identical electronic setting.
for u in 6 12 24; do
  for r0 in A B AB; do point 2 "$r0" "$u"; done
done
point 3 zero 0
for r0 in A B AB; do point 3 "$r0" 24; done
