#!/usr/bin/env bash
# Usage: count_warnings.sh LOG [EXTRA_EXCLUDE_REGEX]
#
# Print the unique compiler/clang-tidy warnings and errors in a build log,
# leaving out the third-party submodules and linker notes, then a
# "Warnings: N" line. Fails when N is not zero.

set -u -o pipefail

log=$1
exclude='Submodules/amrex|Submodules/sundials|ld: warning:|lto-wrapper: warning:'
if [ -n "${2:-}" ]; then
  exclude="${exclude}|$2"
fi

warnings=$(grep -E 'Warning:|Error:|warning:|error:' "${log}" | grep -Ev "${exclude}" | sort -u)
count=0
if [ -n "${warnings}" ]; then
  printf '%s\n' "${warnings}"
  count=$(printf '%s\n' "${warnings}" | grep -c '')
fi
echo "Warnings: ${count}"
[ "${count}" -eq 0 ]
