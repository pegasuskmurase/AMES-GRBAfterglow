"""Forward-shock spectrum at one observer time, saved from a returned vector.

The saved text table can be reloaded and plotted without rerunning the model.
See README.md for component ordering, units, and command examples.
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
            'E_ej':        1.5e55,
            'Gamma0':      560,
            'n_ex':        1,
            'A_star':      0,
            'k_ex':        0,
            'spectral_e':  2.2,
            'epsilon_e':   0.025,
            'fraction_e':  1.0,
            'eta_acc_e':   1.0,
            'epsilon_B':   0.0006,
            'open_angle':  0.1,
            'view_angle':  0.,
            'gaussian_cone': 0.,
            'jet_index':   0.,
            'T_ej':        1e2,
            'L_ph':         0,
        }


class GRB:
    def __init__(self, output_folder="result/spectrum"):
        self.p = PARAMETERS.copy()
        self.source, self.grb = build_model(self.p)
        self.output_folder = Path(output_folder)
        self.output_folder.mkdir(parents=True, exist_ok=True)

    def calc_spectrum(self, time):
        """Save the energy grid and forward-shock components in column order.

        Spectrum returns [component, energy], unlike Flux's [time, band].
        Components 0/1/2 are forward-shock synchrotron/SSC/EIC; 3/4/5 are
        reverse-shock components and are omitted from this forward-shock example.
        """
        print(f"Computing spectrum at observer time {time:g} s", flush=True)
        self.flux_vector = np.asarray(self.grb.Spectrum(time))
        energy = np.asarray(self.source.getPhoton().getEnergy())
        forward = self.flux_vector[:3].T
        data = np.column_stack((energy, forward, forward.sum(axis=1)))
        path = self.output_folder / "spectrum_vector.txt"
        np.savetxt(path, data, fmt="%.17e", header=(
            f"Observer time_s = {time:.17e}\n"
            "energy_eV synchrotron_FS SSC_FS EIC_FS total_FS\n"
            "Flux columns: nu*F_nu [erg cm^-2 s^-1]; reverse shock disabled"))
        print(f"Saved {path}", flush=True)
        return self.flux_vector


def plot_spectrum(path, figure_dir="figure", show=True):
    """Reload a saved table; no engine construction or calculation is needed."""
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] != 5:
        raise ValueError("Expected energy, synchrotron, SSC, EIC, and total columns")
    fig, ax = plt.subplots(figsize=(7, 5))
    for column, label in enumerate(("synchrotron (FS)", "SSC (FS)", "EIC (FS)", "total (FS)"), 1):
        if np.any(data[:, column] > 0):  # A disabled component may be all zero.
            ax.plot(data[:, 0], data[:, column], lw=2, label=label)
    ax.set(xscale="log", yscale="log", xlim=[1e-7, 1e15], ylim=[1e-15, 1e-7], xlabel="Photon energy [eV]",
           ylabel=r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
    ax.legend(fontsize=8)
    fig.tight_layout()
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(figure_dir / f"spectrum_vector.{suffix}")
    if show:
        plt.show()
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time", type=float, default=1e5, help="Observer time [s] (default: 1e5)")
    parser.add_argument("--output-dir", type=Path, default=Path("result/spectrum"))
    parser.add_argument("--figure-dir", type=Path, default=Path("figure"))
    parser.add_argument("--plot-only", action="store_true", help="Reload saved data; --time is ignored")
    parser.add_argument("--no-show", action="store_true", help="Save figures without opening a window")
    args = parser.parse_args()
    if not np.isfinite(args.time) or args.time <= 0:
        parser.error("--time must be finite and positive")
    if not args.plot_only:
        model = GRB(args.output_dir)
        model.calc_spectrum(args.time)
    plot_spectrum(args.output_dir / "spectrum_vector.txt", args.figure_dir, show=not args.no_show)


if __name__ == "__main__":
    main()
