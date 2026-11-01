# Python API

This page describes `AMES.GRBafterglow`, based on
`src/module/main_bindings.cpp` and `src/GRBAfterglow.cpp`.
Vectors returned through the current STL bindings are Python lists, not NumPy
arrays. Use `np.asarray(result)` for array operations.

## Construction and lifetime

```python
import AMES

source = AMES.Source()
source.InitSource()
afterglow = AMES.GRBafterglow(source)
afterglow.haveReverseShock(False)
```

Keep `source` alive while using `afterglow`: the engine holds a C++ reference.
Set momentum grids before initialization and construction, as demonstrated in
`example/GRBAfterglow/example_utils.py`. Do not resize them after constructing
the engine, which caches grids and allocates buffers at construction.
The snippet above only illustrates construction, not a complete calculation.

## Configuration

| Method | Meaning |
| --- | --- |
| `setGRBafterglowParam(params)` | Set the 18 ordered forward-shock parameters; see [parameter table](parameters.md) |
| `haveSSCSpec(enabled)` | Toggle SSC calculation |
| `haveAttenuSSA(enabled)` | Toggle synchrotron self-absorption |
| `haveAttenuGGSource(enabled)` | Toggle internal gamma-gamma attenuation |
| `haveReverseShock(enabled)` | Toggle the reverse-shock path; examples explicitly use `False` |
| `haveOneZone(enabled)` | Select the one-zone path; outside the introductory examples |
| `showInfo(enabled)` | Toggle some diagnostics; some C++ output remains unconditional |

These boolean flags default to `False` in the class definition. The examples
explicitly enable their desired radiation and attenuation switches.

## `Flux(time_array, energy_array_min, energy_array_max)`

Compute forward-shock total flux at every combination of time and band.
Times are observer seconds; energy limits are eV. Minimum and maximum arrays
must have equal lengths, ordered limits, and energies inside the photon grid.
Return shape: **[number of times, number of bands]**.

Equal energy limits request monochromatic nu F_nu. Unequal limits request
integrated energy flux. The current implementation treats a separation below
1e-5 eV as monochromatic. All returned fluxes have units erg cm^-2 s^-1.
The forward total is synchrotron + SSC + EIC; it does not include RS components.

## `Spectrum(T)`

Compute component spectra at observer time `T` [s]. Return shape:
**[6, number of photon grid points]**.

| Row | Component |
| --- | --- |
| 0 | Forward-shock synchrotron |
| 1 | Forward-shock SSC |
| 2 | Forward-shock EIC |
| 3 | Reverse-shock synchrotron |
| 4 | Reverse-shock SSC |
| 5 | Reverse-shock EIC |

Get the matching energy axis [eV] from `source.getPhoton().getEnergy()`.
Fluxes represent nu F_nu [erg cm^-2 s^-1]. Disabled components may be zero.
Save the first three rows transposed beside the energy grid for FS output.

## Additional exposed methods

These methods are exposed but are not required by the introductory workflow.
Exposure does not establish validation of a feature.

| Method | Contract / caveat |
| --- | --- |
| `Flux1D(time, freq_min, freq_max)` | One band per time, returning a list of length N; all three arrays must have length N. Despite names, limits are **eV**, not Hz. Uses the same 1e-5 eV threshold as `Flux`. |
| `set_params_mcmc(theta)` | Reads seven entries: log10(E_ej), log10(n_ex), p, log10(epsilon_e), log10(epsilon_B), log10(fraction_e), view_angle. No length/range checks here. Set the base model first. |
| `calculate_lightcurve(time, freq)` | Paired time/energy arrays of equal length; despite `freq`, input is **eV**. Initialize via the MCMC setter first. Returns a list of N monochromatic FS fluxes. |
| `setGRBafterglowRSParam(params)` | Five entries, documented in [parameters](parameters.md); outside the FS guide. |
| `setOutputFolder(path)` | Select the existing directory used by `FluxOutput`; does not create it. |
| `FluxOutput(times, low, high)` | C++ writer, returns `None`; writes `flux.dat` and time-blocked `spectrum.dat`, plus separate RS files when enabled. FS spectrum columns: energy, synchrotron, SSC, EIC. |

`FluxOutput` uses a relative band-equality threshold of 1e-10 times the minimum
energy, unlike `Flux`. Do not assume nearly equal band limits behave identically.
Its file output uses six digits after the decimal in scientific notation.

## Known interface limitations

Input validation is incomplete: check finite positive times, array lengths,
parameter order, and energy bounds in calling code. A wrong-size forward
parameter list prints an error and returns rather than raising a Python exception.
Several binding docstrings are stale (21 FS parameters, eight RS parameters,
and NumPy return types); this guide follows the implementation instead.
C++ public declarations are not automatically available in Python.
