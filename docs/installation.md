# Installation

## Requirements

- CMake 3.14 or newer.
- A C++17 compiler and an OpenMP runtime/development headers.
- Python with development headers; the current local setup has been checked
  with Python 3.11.15 and pybind11 3.0.1.
- pybind11 with its CMake package available.
- NumPy and Matplotlib for the Python examples.

Activate your intended Python environment before building. On macOS with Conda,
use the environment's `libomp` and headers to share the OpenMP runtime with its
numerical libraries. CMake prefers these when present, then falls back to
`/opt/homebrew/opt/libomp`. Other installations can specify
`-DAMES_OPENMP_ROOT=/path/to/openmp/prefix` on macOS. On Linux, OpenMP is discovered
through CMake's standard compiler checks.

## Build and install

Build in a separate directory, then install the module into your active Python
environment. Replace the example paths below with your own absolute paths.
Activate your intended environment first (for example, `conda activate crpropa3`).

```sh
python --version
mkdir -p /path/to/AMES-install
cd /path/to/AMES-install

cmake /path/to/AMES-GRBAfterglow \
  -DPython_EXECUTABLE="$(python -c 'import sys; print(sys.executable)')"
cmake --build . --parallel
cmake --install .
```

With the default Unix Makefiles generator, the last two commands can also be
written as `make -j4` and `make install`. No `install.sh` is required.
`AMES-install` is a build directory; the install step copies the module into
the selected Python environment. You do not need to add the build directory to
`PATH` or `PYTHONPATH` for normal use.

Verify the installation from another directory, with the same environment active:

```sh
cd /path/to/AMES-GRBAfterglow
python -c "import AMES; print(AMES.__file__)"
```

The printed module path should point into that environment's `site-packages`,
rather than an older build directory. Remove any previous AMES build-directory
entries from `PYTHONPATH` if they cause Python to load an old copy.

After changing C++ code, rebuild and reinstall from the same build directory:

```sh
cd /path/to/AMES-install
cmake --build . --parallel
cmake --install .
```

Start a new Python process after reinstalling to load the rebuilt extension.

The build produces a compiled extension such as `AMES.cpython-311-darwin.so`,
not an `AMES.py` file. Installation places it in the selected Python's
`site-packages` directory, printed as `AMES Python install directory` during
configuration. The environment must be writable. `CMAKE_INSTALL_PREFIX` does
not override this project's absolute Python installation destination.

If CMake cannot find pybind11, install its development package in your environment
or pass `-Dpybind11_DIR=/path/to/pybind11/share/cmake/pybind11`. If pybind11 is
installed as a Python package, its location is available with
`python -m pybind11 --cmakedir`.

When switching Python environments, use a new build directory to avoid mixing
cached Python headers and libraries. To use the local build before installing,
add its directory to `PYTHONPATH`. The compiled module uses the repository path
for its data directory; keep the source/data tree available when running it.

