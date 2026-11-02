# Parameters and units

## Forward-shock parameters

`setGRBafterglowParam(params)` requires exactly 18 values in this order.
The table records validation bounds and fallback values from the setter,
not necessarily constructor defaults or values recommended for a physical model.

| Index (zero-based) | Name | Meaning | Unit | Accepted range | Out-of-range fallback |
| --- | --- | --- | --- | --- | --- |
| 0 | `z` | Redshift | dimensionless | [0.0, 10.0] | 0.0 |
| 1 | `dl` | Luminosity distance | cm | [1e10, 1e29] | 1e26 |
| 2 | `E_ej` | Isotropic-equivalent ejecta energy | erg | [1e45, 1e57] | 1e52 |
| 3 | `Gamma0` | Initial bulk Lorentz factor | dimensionless | [1.0, 1e4] | 100.0 |
| 4 | `n_ex` | External density parameter | cm^-3 for k_ex=0 | [1e-4, 1e6] | 1.0 |
| 5 | `A_star` | Wind normalization | dimensionless | [0, 1e5] | 1.0 |
| 6 | `k_ex` | External density power-law index | dimensionless | [0.0, 3.0] | 0.0 |
| 7 | `spectral_e` | Positive electron spectral index p | dimensionless | [0.0, 5.0] | 2.2 |
| 8 | `epsilon_e` | Electron energy fraction | dimensionless | [1e-8, 1.0] | 0.1 |
| 9 | `fraction_e` | Accelerated electron fraction | dimensionless | [1e-4, 1.0] | 1.0 |
| 10 | `eta_acc_e` | Electron acceleration parameter | dimensionless | [1.0, 1e6] | 1000.0 |
| 11 | `epsilon_B` | Magnetic energy fraction | dimensionless | [1e-8, 1.0] | 1e-3 |
| 12 | `open_angle` | Jet opening angle | rad | [0.01, 1.5] | 0.1 |
| 13 | `view_angle` | Viewing angle | rad | [0.0, 3.0] | 0.0 |
| 14 | `gaussian_cone` | Gaussian angular scale | rad | [0.0, 1.5] | 0.1 |
| 15 | `jet_index` | Jet structure index | dimensionless | [0.0, 10.0] | 0.0 |
| 16 | `T_ej` | Ejecta duration | s | [0.00001, 1e6] | 10.0 |
| 17 | `L_ph` | Prompt luminosity parameter | erg s^-1 | [0.0, 1e60] | 1e50 |

An out-of-range value prints a warning and is replaced with the fallback;
it does not raise an exception. Wrong-length FS input prints an error and leaves
the model unchanged. Explicitly reject NaN/infinite values in your own input
handling: the range comparisons do not reliably reject NaN.

`k_ex=0` is the constant-density case. For a wind, use `k_ex=2` and a positive
`A_star`; the setter substitutes `n_ex = 3e35 * A_star` as the radial density
normalization. Do not interpret that substituted value as a density in cm^-3.
General density profiles and jet structures require further physical review.

The example dictionaries are explicit model choices. They are not defaults and
should not be treated as fits. The positive electron index p describes a
power-law distribution proportional to energy^-p.

## Other axes

- Observer times: seconds.
- Photon energies and energy-band limits: eV, including arguments named `freq`
  in `Flux1D` and `calculate_lightcurve`.
- Momentum grids: the native p*c convention in eV.
- Output monochromatic flux: nu F_nu in erg cm^-2 s^-1.
- Output band-integrated flux: erg cm^-2 s^-1.

## Reverse-shock setter (advanced)

The current implementation reads five values in this order:
`spectral_e_RS`, `epsilon_e_RS`, `fraction_e_RS`, `eta_acc_e_RS`, `epsilon_B_RS`.
Supply exactly five, although its helper currently accepts and ignores trailing
entries. The binding docstring saying eight entries is outdated.
Reverse shock is disabled in the introductory examples. This section records
the input layout; it is not validation of that physical branch.
