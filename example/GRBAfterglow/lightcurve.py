import numpy as np
import matplotlib.pyplot as plt
import os
import AMES

# -------------------------------
# Initialize AMES engine
# -------------------------------

def update_s(s):
 num_p, min_p, max_p = 121, 1e7, 1e19
 s.getTarget().setMomentum(161, 1e-9, 1e7);
 s.getPhoton().setMomentum(251, 1e-9, 1e16);
 s.getElectron().setMomentum(361, 1e4, 1e16);

s = AMES.Source()
update_s(s)
s.InitSource()
grb = AMES.GRBafterglow(s)


# ===============================================================
#   GRB AFTERGLOW MODEL (ISM vs WIND)
# ===============================================================
class GRB:
    """
    medium_type : "ISM"  -> constant density (k = 0)
                  "wind" -> stellar wind (k = 2)
    """
    def __init__(self):

        # -------------------------------
        # Forward shock microphysics
        # -------------------------------
        p = {
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
            'L_ph':         0,
            'acc_p':        0.1,
            'eta_acc_p':    1,
            'spectral_p':    2,
        }

        self.p = p

        param = [
            p['z'], p['dl'], p['E_ej'], p['Gamma0'], p['n_ex'], p['A_star'], p['k_ex'],
            p['spectral_e'], p['epsilon_e'], p['fraction_e'], p['eta_acc_e'],
            p['epsilon_B'], p['open_angle'], p['view_angle'],
            p['gaussian_cone'], p['jet_index'], p['T_ej'], p['L_ph']
        ]

        grb.setGRBafterglowParam(param)

        # -------------------------------
        # Reverse shock microphysics
        # -------------------------------
        p_RS = {
            'spectral_e_RS': 2.2,
            'epsilon_e_RS': 0.1,
            'fraction_e_RS': 1.,
            'eta_acc_e_RS': 1.,
            'epsilon_B_RS': 0.1,
        }

        self.p_RS = p_RS

        param_RS = [
            p_RS['spectral_e_RS'],
            p_RS['epsilon_e_RS'],
            p_RS['fraction_e_RS'],
            p_RS['eta_acc_e_RS'],
            p_RS['epsilon_B_RS']
        ]

        grb.setGRBafterglowRSParam(param_RS)

        # -------------------------------
        # Arrays in time and energy
        # -------------------------------
        self.time_array = np.logspace(np.log10(1e1), np.log10(1e7), 61)

        self.energy_array_min = [1e9 / AMES.eV2Hz,
                                 1e11 / AMES.eV2Hz,
                                 1e0,
                                 1e3,
                                 1e9,
                                 1e11]

        self.energy_array_max = self.energy_array_min.copy()

        # Create output directories if needed
        os.makedirs("result", exist_ok=True)
        os.makedirs("figure", exist_ok=True)
        self.output_folder = "result"

    # ===============================================================
    #   CALCULATE FLUX
    # ===============================================================
    def flux(self):
        grb.haveSSCSpec(True)
        grb.haveAttenuSSA(True)
        grb.haveAttenuGGSource(True)
        grb.haveReverseShock(False)
        grb.showInfo(False)

        flux_vector = grb.Flux(self.time_array,
                               self.energy_array_min,
                               self.energy_array_max)

        np.savetxt(self.output_folder + "flux.txt", flux_vector)
        print("Flux calculation complete → result/flux.txt saved.")

        fig, ax = plt.subplots(figsize=(7, 5))

        for j, ene in enumerate(self.energy_array_min):
            dum = [flux_vector[i][j] for i in range(len(flux_vector))]
            ax.plot(self.time_array, dum, '--', lw=2, label=f"E{j}")

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e1, 1e7])
        ax.set_ylim([1e-19, 1e-4])
        ax.set_xlabel('T [s]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg\ cm^{-2}\ s^{-1}]$', fontsize=15)
        ax.legend(fontsize=7)

        plt.tight_layout()
        plt.savefig("figure/flux.pdf")
        plt.show()

    # ===============================================================
    #   CALCULATE FLUX Output to file
    # ===============================================================
    def calc_flux_output(self):
        grb.haveSSCSpec(True)
        grb.haveAttenuSSA(True)
        grb.haveAttenuGGSource(True)
        grb.haveReverseShock(True)
        grb.showInfo(True)
        grb.setOutputFolder(self.output_folder)

        grb.FluxOutput(self.time_array,
                               self.energy_array_min,
                               self.energy_array_max)

    # ===============================================================
    #   PLOT RESULTS
    # ===============================================================
    def plot_flux_output(self):
        fig, ax = plt.subplots(figsize=(7, 5))

        data = np.loadtxt(self.output_folder + '/flux.dat')
        ax.plot(data[:, 0], data[:, 1], 'C0-')
        ax.plot(data[:, 0], data[:, 2], 'C1-')
        ax.plot(data[:, 0], data[:, 3], 'C2-')
        ax.plot(data[:, 0], data[:, 4], 'C3-')
        ax.plot(data[:, 0], data[:, 5], 'C4-')

        #reverse shock will be checked
        data = np.loadtxt(self.output_folder + '/flux_RS.dat')
        ax.plot(data[:, 0], data[:, 1], 'C0--')
        ax.plot(data[:, 0], data[:, 2], 'C1--')
        ax.plot(data[:, 0], data[:, 3], 'C2--')
        ax.plot(data[:, 0], data[:, 4], 'C3--')
        ax.plot(data[:, 0], data[:, 5], 'C4--')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e1, 1e7])
        ax.set_ylim([1e-19, 1e-4])
        ax.set_xlabel('T [s]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg\ cm^{-2}\ s^{-1}]$', fontsize=15)
        ax.legend(fontsize=7)

        plt.tight_layout()
        plt.savefig("figure/flux.pdf")
        plt.show()

        num = len(s.getPhoton().getMomentum())
        fig, ax = plt.subplots(figsize=(7, 5))
        data = np.loadtxt(self.output_folder + '/spectrum.dat')
        for i, x in enumerate(self.time_array):
            idx1 = num * i
            idx2 = num * (i + 1)
            ax.plot(data[idx1:idx2, 0], data[idx1:idx2, 1], '--', lw=2, label=str(x) + ' s, EATS')

        data = np.loadtxt(self.output_folder + '/spectrum_RS.dat')
        for i, x in enumerate(self.time_array):
            idx1 = num * i
            idx2 = num * (i + 1)
            #ax.plot(data[idx1:idx2, 0], data[idx1:idx2, 1], '-.', lw=2, label=str(x) + ' s, EATS')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e-8, 1e15])
        ax.set_ylim([1e-20, 1e-9])
        ax.set_xlabel('E [eV]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg \ cm^{-2} \ s^{-1}]$', fontsize=15)
        ax.legend(fontsize=7)
        plt.savefig('figure/spectrum.pdf')
        plt.show()

# ===============================================================
#   RUN MODEL
# ===============================================================

g = GRB()
#g.flux()
g.calc_flux_output()
g.plot_flux_output()
