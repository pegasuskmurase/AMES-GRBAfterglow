#!/usr/bin/env bash
# Build and install into the selected Python environment. No files are deleted.
# Usage: bash install.sh [additional CMake configure arguments]
set -euo pipefail

source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${AMES_BUILD_DIR:-${source_dir}/build}"
if command -v python >/dev/null 2>&1; then
    python_cmd=python
else
    python_cmd=python3
fi
python_executable="$("${python_cmd}" -c 'import sys; print(sys.executable)')"

cmake -S "${source_dir}" -B "${build_dir}" \
    "-DPython_EXECUTABLE=${python_executable}" "$@"
cmake --build "${build_dir}" --parallel
cmake --install "${build_dir}"

# Use the actual configured interpreter, including a command-line override.
configured_python="$(sed -n 's/^Python_EXECUTABLE:[^=]*=//p' "${build_dir}/CMakeCache.txt")"
"${configured_python}" -I -c \
    'import AMES; print("Installed AMES import OK:", AMES.__file__)'
