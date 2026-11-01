# Quick start

After [installation](installation.md), activate the same Python environment.
From the repository root:

```sh
cd example/GRBAfterglow
python lightcurve.py
python spectrum.py --time 1e5
```

Each script calculates once, converts the returned nested lists to a NumPy
array, writes a table, and plots that saved table. Spectrum calculations can take
substantial time. Use `--help` for all command-line options.

## Replot without calculating

```sh
python lightcurve.py --plot-only
python spectrum.py --plot-only
```

For noninteractive use, add `--no-show`; on a headless machine also set
`MPLBACKEND=Agg`. Change destinations with `--output-dir` and `--figure-dir`.
Existing files at the same paths are overwritten.

## Saved tables

| File, relative to the example directory | Columns |
| --- | --- |
| `result/lightcurve/flux_vector.txt` | observer time [s], one flux column per band |
| `result/lightcurve/bands.txt` | minimum and maximum photon energy [eV] per band |
| `result/spectrum/spectrum_vector.txt` | energy [eV], FS synchrotron, FS SSC, FS EIC, FS total |

Flux units are erg cm^-2 s^-1. Monochromatic fluxes represent nu F_nu.
The spectrum table header records the observer time. The saved tables are not
complete model checkpoints: retain parameter settings and the source revision.

```python
import numpy as np

data = np.loadtxt("result/lightcurve/flux_vector.txt", ndmin=2)
time_s, flux = data[:, 0], data[:, 1:]

spectrum = np.loadtxt("result/spectrum/spectrum_vector.txt", ndmin=2)
energy_eV, total_flux = spectrum[:, 0], spectrum[:, 4]
```

The light-curve example has six channels (1 GHz, 100 GHz, 1 eV, 1 keV, 1 GeV,
100 GeV) and 61 times from 10 to 10^7 s. Radio frequencies are converted to eV
using `AMES.eV2Hz`. Match the parameters and grids in both examples before
comparing their predictions. Axis limits remain configured in the scripts.
