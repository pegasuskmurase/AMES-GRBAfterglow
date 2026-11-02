# AMES-GRBAfterglow

AMES-GRBAfterglow is a C++ afterglow model with a pybind11 Python interface.
This first documentation version focuses on the forward-shock examples:
calculate light curves and spectra, save returned vectors, and replot results.

1. [Install AMES](installation.md) into your Python environment.
2. [Run the examples](quickstart.md).
3. Check [parameters and units](parameters.md) and the [Python API](python-api.md).
4. See [developer notes](developer-guide.md) for the initial function audit.

## Scope

The examples enable SSC, synchrotron self-absorption, and source gamma-gamma
attenuation, with reverse shock disabled. Other code paths exist, but their
presence or Python exposure is not a claim of scientific validation.
The two examples use different illustrative model parameters.

These pages describe the current source. They are a starting point for review,
not a complete description of the underlying physics. Check numerical
convergence and physical assumptions before interpreting scientific results.
