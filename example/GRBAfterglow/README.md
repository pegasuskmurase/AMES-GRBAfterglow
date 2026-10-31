# Forward-shock afterglow examples

Keep `lightcurve.py` and `spectrum.py` separate: they answer different questions.
A light curve follows selected photon energies over time; a spectrum resolves
emission components across energy at one time. `example_utils.py` shares only
engine setup, grids, and the required parameter order. Importing either example
does not start a calculation.

These are illustrative parameter sets, not fits to the same event: the two
scripts retain their different original forward-shock parameters. Edit each
script's `PARAMETERS` dictionary to change its model. To compare a spectrum with
a light curve, use the same dictionary and grids in both. The setter requires
exactly 18 parameters; do not append reverse-shock or proton parameters.
`example_utils.PARAMETER_ORDER` defines the required order. Distances are in cm,
energies in erg for `E_ej`, times in seconds, and angles in radians. The electron
spectral index is positive (for a distribution proportional to energy^-p).
The default `k_ex=0` describes a constant-density medium. For a wind setup use
`k_ex=2` with the appropriate `A_star` normalization, and verify the parameters
for your physical model.

## Run

Build and install AMES following the [project README](../../README.md).
Activate the Python environment where AMES, NumPy, and Matplotlib are installed.
From this directory:

```sh
# Returned nested vector -> NumPy array -> text file -> plot (default).
python lightcurve.py

# Returned component vectors -> text file -> plot at a single observer time.
python spectrum.py --time 1e5
```

Both scripts explicitly disable reverse shock and enable SSC, synchrotron
self-absorption, and source gamma-gamma attenuation. EIC is included in the saved
component layout, but may be zero for this configuration. Reverse-shock fields
may still appear in C++ diagnostic messages; that does not enable reverse shock.

The light-curve calculation uses 61 times from 10 to 10^7 seconds and six energy
channels: 1 GHz, 100 GHz (converted to eV), 1 eV, 1 keV, 1 GeV, and 100 GeV.
`energy_array_min == energy_array_max` requests monochromatic nu F_nu.
Unequal limits request band-integrated energy flux. Use exactly equal limits
for monochromatic requests. The current C++ `Flux` implementation treats limits
separated by less than 1e-5 eV as monochromatic; avoid such narrow intervals
when requesting integrated flux.

Paths are relative to the working directory. Use `--output-dir PATH` and
`--figure-dir PATH` to choose destinations, or to preserve different runs.
Running again in the same output directory overwrites the saved files.
Use `--no-show` to save plots without opening windows; on a headless machine,
set `MPLBACKEND=Agg` as well. Use `--help` for all options.

## Reuse saved results without computing again

```sh
python lightcurve.py --plot-only
python spectrum.py --plot-only
```

Use the same `--output-dir` that was used to save the calculation. In spectrum
plot-only mode the time is the one recorded in the file, not `--time`.
The plot-only paths read files without constructing a Source or afterglow engine.

You can also load the tables directly in another Python program:

```python
import numpy as np

lightcurve = np.loadtxt("result/lightcurve/flux_vector.txt", ndmin=2)
time_s = lightcurve[:, 0]
flux = lightcurve[:, 1:]  # Shape: [number of times, number of bands].
bands_eV = np.loadtxt("result/lightcurve/bands.txt", ndmin=2)

spectrum = np.loadtxt("result/spectrum/spectrum_vector.txt", ndmin=2)
energy_eV = spectrum[:, 0]
synchrotron, ssc, eic, total = spectrum[:, 1:].T
```

Text headers begin with `#`, which `np.loadtxt` skips. The NumPy writer saves
17-digit scientific notation, the coordinate axis, and column descriptions.
Flux columns have units erg cm^-2 s^-1; monochromatic spectra represent nu F_nu.
The files describe results and axes, but are not complete model checkpoints:
retain the parameter dictionaries and code version to reproduce a calculation.

## Output layouts

| Method | Default file | Numeric columns |
| --- | --- | --- |
| `GRB.flux()` / `Flux` | `result/lightcurve/flux_vector.txt` | time_s, six band fluxes |
| Light-curve band definitions | `result/lightcurve/bands.txt` | energy_min_eV, energy_max_eV; one row per band |
| `GRB.calc_spectrum()` / `Spectrum` | `result/spectrum/spectrum_vector.txt` | energy_eV, synchrotron_FS, SSC_FS, EIC_FS, total_FS |

`Flux` returns **[time, band]**. `Spectrum` returns **[component, energy]**, so
transpose the first three component rows before saving them beside the energy
axis. The remaining three rows are reverse-shock components, omitted here.

Both examples use one workflow: calculate a vector, save it with NumPy, then
plot the saved table. `--plot-only` skips the calculation. The C++ `FluxOutput`
API remains available in AMES for advanced use, but is not exposed as a second
mode in these introductory scripts. To save spectra at several times, call
`Spectrum` for each time and save each result in a separate output directory.

Figures are `lightcurve_vector.pdf` and `spectrum_vector.png` / `.pdf` under
`figure/`. The spectrum plot uses energy limits [1e-7, 1e15] eV and flux limits
[1e-15, 1e-7] erg cm^-2 s^-1.

## Before using results scientifically

These scripts demonstrate the API and storage formats. Check convergence with
respect to momentum and jet grids, and validate the physical parameters for
your application. A single spectrum may take substantial time in C++; a pause
after `Spectrum ... Begin` does not by itself imply a deadlock. Existing C++
calculations may delay Python's handling of Ctrl+C.
