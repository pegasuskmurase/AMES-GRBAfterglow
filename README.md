# AMES-GRBAfterglow

[Documentation](docs/index.md) · [Read on GitHub](https://github.com/pegasuskmurase/AMES-GRBAfterglow/blob/main/docs/index.md)

AMES-GRBAfterglow provides a C++ gamma-ray burst afterglow model with a Python
interface built with pybind11.  The code generates templates for light curves and spectra, 
accounting for synchrotron emission, synchrotron self-Compton emission, external inverse-Compton emission, 
and synchrotron self-absorption. It supports both top-hat and structured jets viewed either on-axis or off-axis.
Example python files save the returned vectors as text tables, and plot them without repeating the calculation.

Developers: Kohta Murase, Bing Theodore Zhang, Shiqi Yu, Yu-Jia Wei, Yuri Sato
[CITATION: Zhang, Murase Veres & Meszaros Astrophys.J. 920 (2021) 55, Zhang, Murase, Ioka & Zhang JHEAp 45 (2025) 392, see also Sato et al. Phys.Rev.D.Lett. 110 (2024) L061307]

## Requirements

- CMake 3.18 or newer.
- A C++17 compiler and an OpenMP runtime with development headers.
- Python 3 with development headers.
- pybind11 with its CMake package available.
- NumPy and Matplotlib for the Python examples.

Activate your intended Python environment before building and running AMES;
use a fresh build directory when switching environments. See
[Conda and virtual environments](docs/installation.md#python-environments).

On macOS, provide matching `libomp` and headers in the selected Conda environment,
or use `brew install libomp` for a Homebrew setup. See
[macOS and Apple Silicon setup](docs/installation.md#macos-and-apple-silicon).

## Build and install

From the source checkout, configure into the separate `AMES-GRBAfterglow-build/` directory.
Installation writes to the selected Python environment's `site-packages`,
which must be writable. Replace the source path below with your own.

```sh
cd /path/to/AMES-GRBAfterglow  # source directory
cmake -S . -B AMES-GRBAfterglow-build \
  -DPython_EXECUTABLE="$(python3 -c 'import sys; print(sys.executable)')"
cmake --build AMES-GRBAfterglow-build --parallel
cmake --install AMES-GRBAfterglow-build
python3 -c "import AMES; print(AMES.__file__)"
```

The import check should print a path in that environment's `site-packages`.
Keep the source/data tree available at its original path when running AMES.
See the [installation guide](docs/installation.md) for rebuilding and module
discovery, [optimization options](docs/installation.md#optimization-options),
and [troubleshooting](docs/installation.md#troubleshooting).

## Examples

With the same Python environment active, start from the repository root:

```sh
cd example/GRBAfterglow
python3 lightcurve.py
python3 spectrum.py --time 1e5
```

Calculations can take substantial time; add `--plot-only` to reuse saved tables.
See the [example guide](example/GRBAfterglow/README.md) and
[quick start](docs/quickstart.md) for parameters, units, outputs, and plotting.
Reverse-shock emission is not included in this version.

## Documentation

The [documentation index](docs/index.md) links to parameter tables, the Python
API, and developer notes, including local documentation preview instructions.

Read relative Markdown links in a Markdown preview or on GitHub. The
[documentation website](https://pegasuskmurase.github.io/AMES-GRBAfterglow/)
currently returns 404; use the Markdown guides until
[GitHub Pages deployment](docs/developer-guide.md#github-pages-deployment) succeeds.

See [LICENSE](LICENSE) for licensing terms.
