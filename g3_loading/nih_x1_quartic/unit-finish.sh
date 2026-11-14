#!/usr/bin/env bash
set -eu
cd -- "$(dirname -- "$0")"
phase=$1
{
 printf 'SERVICE_RESULT=%s\n' "${SERVICE_RESULT:-unknown}"
 printf 'EXIT_CODE=%s\n' "${EXIT_CODE:-unknown}"
 printf 'EXIT_STATUS=%s\n' "${EXIT_STATUS:-unknown}"
 printf 'UTC_UNIT_END='; date -u +%FT%TZ
} >"unit-result-$phase.txt"
printf 'G2_WELLS_20261113_%s_UNIT_END result=%s code=%s status=%s\n' "$phase" "${SERVICE_RESULT:-unknown}" "${EXIT_CODE:-unknown}" "${EXIT_STATUS:-unknown}" >>"run-$phase.log"
