#!/usr/bin/env bash
# Usage: count_warnings.sh LOG [EXTRA_EXCLUDE_REGEX]
#
# Print the unique compiler/clang-tidy warnings and errors in a build log,
# leaving out the third-party submodules and linker notes, then a
# "Warnings: N" line. Fails when N is not zero.

set -u -o pipefail

log=$1
if [ ! -r "${log}" ]; then
  echo "count_warnings.sh: cannot read log ${log}" >&2
  exit 2
fi
exclude='Submodules/amrex|Submodules/sundials|ld: warning:|lto-wrapper: warning:'
if [ -n "${2:-}" ]; then
  exclude="${exclude}|$2"
fi

warnings=$(grep -E 'Warning:|Error:|warning:|error:' "${log}" | grep -Ev "${exclude}" | sort -u)
[ -n "${warnings}" ] && printf '%s\n' "${warnings}"
count=$(grep -c . <<<"${warnings}")
echo "Warnings: ${count}"
[ "${count}" -eq 0 ]
