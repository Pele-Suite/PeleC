#!/usr/bin/env python3
"""Print the CMake targets whose executables a build tree's CTest tests run.

Usage: ctest_targets.py BUILD_DIR [CTEST_ARGS...]

The executable targets are read from the CMake File API reply (codemodel-v2),
which the build tree's configure must have written: create the empty file
BUILD_DIR/.cmake/api/v1/query/codemodel-v2 before running cmake. Every word
of every test command from `ctest --show-only=json-v1 CTEST_ARGS` is matched
by file name against those executables, so absolute paths, relative paths and
bare names (e.g. `--pele_exec PeleC-Spray-Conv`, run from a script) all
count. The target names are printed on one line, separated by spaces, ready for
`cmake --build BUILD_DIR --target ...`. Nothing is printed when no test
command names a target, so callers can fall back to building everything.
"""

import json
import pathlib
import shlex
import subprocess
import sys


def executable_targets(build):
    """Map each executable's file name to the name of the target making it."""
    reply = build / ".cmake" / "api" / "v1" / "reply"
    index = json.loads(max(reply.glob("index-*.json")).read_text())
    codemodel = json.loads(
        (reply / index["reply"]["codemodel-v2"]["jsonFile"]).read_text()
    )
    executables = {}
    for config in codemodel["configurations"]:
        for target in config["targets"]:
            info = json.loads((reply / target["jsonFile"]).read_text())
            if info["type"] != "EXECUTABLE":
                continue
            for artifact in info.get("artifacts", []):
                executables[pathlib.PurePath(artifact["path"]).name] = info["name"]
    return executables


def test_words(build, ctest_args):
    """Yield every word of every test command (shell strings and --key=value
    options are split)."""
    tests = json.loads(
        subprocess.run(
            ["ctest", "--show-only=json-v1", *ctest_args],
            cwd=build,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )
    for test in tests.get("tests", []):
        for arg in test.get("command", []):
            try:
                words = shlex.split(arg)
            except ValueError:
                words = arg.split()
            for word in words:
                yield from word.split("=")


def main():
    build = pathlib.Path(sys.argv[1]).resolve()
    executables = executable_targets(build)
    targets = {
        executables[name]
        for name in (pathlib.PurePath(w).name for w in test_words(build, sys.argv[2:]))
        if name in executables
    }
    print(" ".join(sorted(targets)))


if __name__ == "__main__":
    main()
