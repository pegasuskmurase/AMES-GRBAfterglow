# AMES-GRBAfterglow

AMES-GRBAfterglow provides a C++ gamma-ray burst afterglow model with a Python
interface built with pybind11. The forward-shock examples calculate light curves
and spectra, save the returned vectors as text tables, and reload them for
plotting without repeating the calculation.

Developed by Prof. Kohta Murase, Prof. Bing Theodore Zhang, Prof. Shiqi Yu,
Dr. Yujia Wei, Dr. Riki Matsu, Dr. Yuri Sato, et al.

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

## Examples

After installation, activate the same Python environment and start from the
repository root:

```sh
cd example/GRBAfterglow
python lightcurve.py
python spectrum.py --time 1e5

# Reuse saved tables without running the physics again.
python lightcurve.py --plot-only
python spectrum.py --plot-only
```

Both scripts use the same simple workflow: **calculate a vector, save it with
NumPy, then plot the saved table**. Reverse shock is explicitly disabled.
The light-curve and spectrum examples retain different illustrative parameter
sets; use matching parameters if you want to compare their physical predictions.

See the [example guide](example/GRBAfterglow/README.md) for parameter conventions,
units, output columns, custom directories, and headless plotting. Use each
script's `--help` for command-line options.

## Documentation

For parameter tables, function reference, and the initial C++/Python interface
audit, see the [documentation](docs/index.md). To preview it locally:

```sh
python -m pip install -r requirements-docs.txt
python -m mkdocs serve
```

The documentation can be read as Markdown on GitHub without installing MkDocs.

## Troubleshooting

- **Wrong Python version:** check `AMES Python` in CMake's output. Pass the
  intended interpreter explicitly with `-DPython_EXECUTABLE=...` and configure
  in a new build directory.
- **Import crash with Conda:** rebuild using the current CMake configuration,
  which lets pybind11 handle Python extension linking without an explicit
  `libpython` link. Check `AMES.__file__` to verify which build is loaded.
- **OpenMP Error #15:** rebuild against the same OpenMP runtime used by your
  Python environment. On macOS, check the printed `AMES OpenMP runtime` path.
  Restart Python after rebuilding; do not hide duplicate runtimes with
  `KMP_DUPLICATE_LIB_OK=TRUE`.
- **Long calculation:** a spectrum can take substantial time after
  `Spectrum ... Begin`. Ctrl+C may be delayed inside C++; on a Unix terminal,
  Ctrl+Z followed by `jobs -l` and `kill -KILL %N` stops the suspended job
  (`N` is the job number you checked).

## Repository layout

- `src/`, `include/`: C++ model, physics modules, and Python bindings.
- `data/`: model data used by the C++ code.
- `example/GRBAfterglow/`: documented Python examples and shared setup.
- `CMakeLists.txt`: build and installation configuration.

See [LICENSE](LICENSE) for licensing terms.
