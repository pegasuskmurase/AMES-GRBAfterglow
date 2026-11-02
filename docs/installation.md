# Installation

## Requirements

- CMake 3.18 or newer.
- A C++17 compiler and an OpenMP runtime/development headers.
- Python 3 with development headers; the current local setup has been checked
  with Python 3.11.15 and pybind11 3.0.1.
- pybind11 with its CMake package available.
- NumPy and Matplotlib for the Python examples.

On Linux, OpenMP is discovered through CMake's standard compiler checks.
For macOS, see [macOS and Apple Silicon](#macos-and-apple-silicon) below.

## Python environments

Activate the intended environment before configuring, building, installing, and
running AMES. For an existing Conda environment, use `conda activate ENV_NAME`.
For a virtual environment, create and activate it with your intended Python:

```sh
python3 -m venv /path/to/ames-venv
. /path/to/ames-venv/bin/activate
```

Ensure the dependencies above are available in that environment. Confirm the
selected interpreter and process architecture before configuring:

```sh
python3 -c "import sys, platform; print(sys.executable); print(platform.machine())"
```

The build commands below pass this interpreter explicitly to CMake. Check the
printed `AMES Python` and `AMES Python install directory` values. When switching
Python environments or architectures, use a new build directory to avoid mixing
cached Python headers, libraries, and compiler settings.

## Build and install

Replace the source path below with your own absolute path. CMake uses the source
checkout (`-S .`) and creates a separate build directory
(`-B AMES-GRBAfterglow-build`) inside it.
Installation copies the module into your selected Python environment.

```sh
cd /path/to/AMES-GRBAfterglow
cmake -S . -B AMES-GRBAfterglow-build \
  -DPython_EXECUTABLE="$(python3 -c 'import sys; print(sys.executable)')"
cmake --build AMES-GRBAfterglow-build --parallel
cmake --install AMES-GRBAfterglow-build
python3 -c "import AMES; print(AMES.__file__)"
```

The import check runs from the source directory, outside
`AMES-GRBAfterglow-build/`, to avoid loading the local build directly.
You do not need to add the build directory to `PATH` or `PYTHONPATH` after
installation, and no `install.sh` is required.
To place generated files outside the checkout, replace `AMES-GRBAfterglow-build`
in all three CMake commands with the same absolute build-directory path.

The printed module path should point into that environment's `site-packages`,
rather than an older build directory. Remove any previous AMES build-directory
entries from `PYTHONPATH` if they cause Python to load an old copy.

After changing C++ code, rebuild and reinstall from the same build directory:

```sh
cd /path/to/AMES-GRBAfterglow
cmake --build AMES-GRBAfterglow-build --parallel
cmake --install AMES-GRBAfterglow-build
```

Start a new Python process after reinstalling to load the rebuilt extension.

The build produces a compiled extension such as `AMES.cpython-311-darwin.so`,
not an `AMES.py` file. Installation places it in the selected Python's
`site-packages` directory, printed as `AMES Python install directory` during
configuration. The environment must be writable. `CMAKE_INSTALL_PREFIX` does
not override this project's absolute Python installation destination.

To use the local build before installing, add its directory to `PYTHONPATH`.
The compiled module uses the repository path for its data directory; keep the
source/data tree available at its original path when running it.

## macOS and Apple Silicon

Clang on macOS requires an external OpenMP runtime and matching development
headers. With Conda, use the environment's `libomp` and headers to share the
runtime with its numerical libraries. For a Homebrew setup, install them with
`brew install libomp`.

CMake selects the OpenMP prefix in this order:

1. An explicit `-DAMES_OPENMP_ROOT=/path/to/openmp/prefix`.
2. The selected Python's prefix, if it contains `lib/libomp.dylib`.
3. The prefix reported by `brew --prefix libomp`, then the standard Homebrew
   locations `/opt/homebrew/opt/libomp` and `/usr/local/opt/libomp`.

The chosen prefix must contain both `lib/libomp.dylib` and `include/omp.h`.
CMake stops if either is missing; it does not combine a runtime from one prefix
with headers from another. Check the printed `AMES OpenMP runtime` path.

For a native Apple Silicon build, select an ARM64 Python and matching compiler
and OpenMP runtime. The interpreter check above should report `arm64`. Configure
a fresh build directory with an explicit target architecture:

```sh
cd /path/to/AMES-GRBAfterglow
cmake -S . -B AMES-GRBAfterglow-build-arm64 \
  -DPython_EXECUTABLE="$(python3 -c 'import sys; print(sys.executable)')" \
  -DCMAKE_OSX_ARCHITECTURES=arm64
```

Use `AMES-GRBAfterglow-build-arm64` in the build and install commands above.
The architecture option does not change the selected Python's architecture.
Avoid mixing ARM64 and x86_64 tools or dependencies; a universal OpenMP library
must contain every requested target architecture. CMake checks that OpenMP can
compile and link with the selected compiler, headers, and runtime during configuration.

## Optimization options

Pass these options to the configuration command as needed:

- `-DCMAKE_BUILD_TYPE=Release`: the default for single-configuration generators
  when no build type is supplied. Use `Debug` for debugging; the project does not
  force `-O3` into Debug builds. With a multi-configuration generator, select
  `--config Release` for both the build and install commands.
- `-DAMES_FAST_MATH=OFF`: disable the default `-ffast-math` setting, which may
  affect floating-point results. This setting is independent of build type;
  a Debug build does not disable it automatically.
- `-DAMES_NATIVE_ARCH=ON`: opt into `-march=native` for the build CPU. It defaults
  to `OFF`; keep it off for portable binaries, cross-compilation, and universal
  macOS builds. CMake rejects unsupported compiler flags.

## Troubleshooting

- **Wrong Python version or old module:** check `AMES Python` in CMake's output.
  Pass the intended interpreter with `-DPython_EXECUTABLE=...` and configure in
  a fresh build directory. Check `AMES.__file__` from outside the build directory
  and remove stale build-directory entries from `PYTHONPATH`.
- **pybind11 not found:** CMake tries its normal package discovery, then the
  pybind11 package installed for the selected Python. Check that package with
  `python3 -m pybind11 --cmakedir`, or provide its CMake directory explicitly
  with `-Dpybind11_DIR=/path/to/pybind11/share/cmake/pybind11`. An explicit hint
  takes precedence over the Python-package fallback.
- **Installation permission error:** use a writable Python environment. The
  absolute `site-packages` destination is not changed by `CMAKE_INSTALL_PREFIX`.
- **Import crash with Conda:** rebuild using the current CMake configuration,
  which lets pybind11 handle Python extension linking without an explicit
  `libpython` link. Reinstall, restart Python, and check `AMES.__file__` to verify
  which build is loaded.
- **OpenMP setup or architecture error:** check the printed compiler, Python
  architecture, target architectures, and OpenMP runtime. Install matching
  headers and a runtime supporting the requested architecture, then configure
  a fresh build directory. Inspect CMake's configure diagnostics if its OpenMP
  compile/link check fails.
- **OpenMP Error #15:** rebuild against the same OpenMP runtime used by your
  Python environment. On macOS, check `AMES OpenMP runtime`. Restart Python
  after rebuilding; do not hide duplicate runtimes with
  `KMP_DUPLICATE_LIB_OK=TRUE`.
- **Long calculation:** a spectrum can take substantial time after
  `Spectrum ... Begin`. Ctrl+C may be delayed inside C++; on a Unix terminal,
  Ctrl+Z followed by `jobs -l` and `kill -KILL %N` stops the suspended job
  (`N` is the job number you checked). See the
  [quick start](quickstart.md#replot-without-calculating) to reuse saved results.

