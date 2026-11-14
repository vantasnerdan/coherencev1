#!/usr/bin/env bash
# Preserve exact Native integer counters, then always deliver the existing all-exit marker.
set +e
cd -- "$(dirname -- "$0")"
phase=$1
/usr/bin/timeout 5s systemctl --user show "coherence-$phase.service" -p InvocationID -p CPUUsageNSec -p MemoryPeak -p Result > "unit-accounting-$phase.txt" 2> "unit-accounting-$phase.stderr"
accounting_status=$?
printf '%s\n' "$accounting_status" > "unit-accounting-$phase.exit"
/bin/bash ./unit-finish.sh "$phase"
finish_status=$?
if [ "$accounting_status" -ne 0 ]; then exit "$accounting_status"; fi
exit "$finish_status"
