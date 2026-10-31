"""Forward-shock light curves: calculate, save returned vectors, and replot.

Run with --help for options. See README.md for column layouts and reload usage.
All output paths are relative to the current working directory unless specified.
"""
import argparse
from pathlib import Path

import AMES
import matplotlib.pyplot as plt
import numpy as np

from example_utils import build_model

PARAMETERS = {
            'z':           0.15,
            'dl':          723 * AMES.Mpc,
            'E_ej':        1e51,
            'Gamma0':      3.16e2,
            'n_ex':        100,
            'A_star':      0,
            'k_ex':        0,
            'spectral_e':  2.8,
            'epsilon_e':   0.01,
            'fraction_e':  1.0,
            'eta_acc_e':   1.0,
            'epsilon_B':   0.01,
            'open_angle':  1.,
            'view_angle':  0.,
            'gaussian_cone': 0.,
            'jet_index':   0.,
            'T_ej':        1e0,
            'L_ph':         0
        }


class GRB:
    """Keep the Source alive for as long as its GRBafterglow engine is used."""

    def __init__(self, output_folder="result/lightcurve", figure_folder="figure"):
        self.p = PARAMETERS.copy()
        self.source, self.grb = build_model(self.p)
        self.time_array = np.logspace(1, 7, 61)  # Observer time [s].
        # Convert radio frequencies [Hz] to photon energies [eV].
        self.energy_array_min = np.array([
            1e9 / AMES.eV2Hz, 1e11 / AMES.eV2Hz, 1.0, 1e3, 1e9, 1e11,
        ])
        # Equal limits request monochromatic nu*F_nu; different limits request
        # band-integrated flux. Both are in erg cm^-2 s^-1.
        self.energy_array_max = self.energy_array_min.copy()
        self.output_folder = Path(output_folder)
        self.figure_folder = Path(figure_folder)
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.figure_folder.mkdir(parents=True, exist_ok=True)

    def save_bands(self):
        """Save the energy limits separately so column meanings survive reload."""
        np.savetxt(self.output_folder / "bands.txt",
                   np.column_stack((self.energy_array_min, self.energy_array_max)),
                   header="energy_min_eV energy_max_eV", fmt="%.17e")

    def flux(self):
        """Calculate a [time, band] vector and save time plus every flux column.

        The C++ binding returns nested Python lists. Convert them to an ndarray
        for slicing, saving, or subsequent analysis. Saving and plotting reuse
        this result without calculating the light curve a second time.
        """
        self.flux_vector = np.asarray(self.grb.Flux(
            self.time_array, self.energy_array_min, self.energy_array_max))
        self.save_bands()
        path = self.output_folder / "flux_vector.txt"
        columns = " ".join(f"flux_band_{i + 1}" for i in range(len(self.energy_array_min)))
        np.savetxt(path, np.column_stack((self.time_array, self.flux_vector)),
                   header=f"time_s {columns}\nFlux units: erg cm^-2 s^-1; bands: bands.txt",
                   fmt="%.17e")
        print(f"Saved {path}", flush=True)
        return self.flux_vector


def plot_flux(path, bands_path, figure_path, show=True):
    """Plot the saved vector table without running the afterglow calculation."""
    data = np.loadtxt(path, ndmin=2)
    bands = np.loadtxt(bands_path, ndmin=2)
    if data.shape[1] != len(bands) + 1:
        raise ValueError("Flux columns do not match the saved band definitions")
    fig, ax = plt.subplots(figsize=(7, 5))
    for j, (low, high) in enumerate(bands):
        label = f"{low:.3g} eV" if low == high else f"{low:.3g}–{high:.3g} eV"
        ax.plot(data[:, 0], data[:, j + 1], "--", lw=2, label=label)
    ax.set(xscale="log", yscale="log", xlabel="Observer time [s]",
           ylabel=r"Flux [erg cm$^{-2}$ s$^{-1}$]")
    ax.legend(fontsize=8)
    fig.tight_layout()
    figure_path = Path(figure_path)
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path)
    if show:
        plt.show()
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("result/lightcurve"))
    parser.add_argument("--figure-dir", type=Path, default=Path("figure"))
    parser.add_argument("--plot-only", action="store_true", help="Reload files; skip computation")
    parser.add_argument("--no-show", action="store_true", help="Save figure without opening a window")
    args = parser.parse_args()
    if not args.plot_only:
        model = GRB(args.output_dir, args.figure_dir)
        model.flux()
    plot_flux(args.output_dir / "flux_vector.txt", args.output_dir / "bands.txt",
              args.figure_dir / "lightcurve_vector.pdf", show=not args.no_show)


if __name__ == "__main__":
    main()
