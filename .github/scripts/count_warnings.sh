#!/usr/bin/env bash
# Usage: count_warnings.sh LOG [EXTRA_EXCLUDE_REGEX]
#
# Print the unique compiler/clang-tidy warnings and errors in a build log,
# leaving out the third-party submodules and linker notes, then a
# "Warnings: N" line. Fails when N is not zero. The gcc problem matcher is on
# only while the filtered list prints, so the annotations match the count (on
# the raw build output it would also annotate the excluded warnings).
# Same script in PeleC and PelePhysics.

set -u -o pipefail

log=${1:-}
if [ ! -r "${log}" ]; then
  echo "count_warnings.sh: cannot read log ${log}" >&2
  exit 2
fi
exclude='Submodules/amrex|Submodules/sundials|Submodules/AMReX-Hydro|ld: warning:|lto-wrapper: warning:'
if [ -n "${2:-}" ]; then
  exclude="${exclude}|$2"
fi

warnings=$(grep -E 'Warning:|Error:|warning:|error:' "${log}" | grep -Ev "${exclude}" | sort -u)
count=0
if [ -n "${warnings}" ]; then
  # Annotate only if the matcher is there; the list and count print anyway.
  matcher="$(cd "$(dirname "$0")/../problem-matchers" 2>/dev/null && pwd)/gcc.json"
  if [ -r "${matcher}" ]; then
    echo "::add-matcher::${matcher}"
    # Removed only if added here: it may be registered for another step too.
    trap 'echo "::remove-matcher owner=gcc-problem-matcher::"' EXIT
  fi
  printf '%s\n' "${warnings}"
  count=$(printf '%s\n' "${warnings}" | grep -c '')
fi
echo "Warnings: ${count}"
[ "${count}" -eq 0 ]
