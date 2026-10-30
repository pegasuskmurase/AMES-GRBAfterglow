import numpy as np
import matplotlib.pyplot as plt
import os
import AMES

# Initialize AMES
s = AMES.Source()
s.InitSource()
grb = AMES.GRBafterglow(s)


# ======================================================================
#   GRB AFTERGLOW MODEL (supports ISM vs Wind)
# ======================================================================
class GRB:
    def __init__(self):
        # -----------------------------------------------------------
        # Forward shock microphysics
        # -----------------------------------------------------------
        p = {
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
            'acc_p':        0.1,
            'eta_acc_p':    1,
            'spectral_p':    2,
        }

        self.p = p

        param = [
            p['z'], p['dl'], p['E_ej'], p['Gamma0'], p['n_ex'], p['A_star'], p['k_ex'],
            p['spectral_e'], p['epsilon_e'], p['fraction_e'], p['eta_acc_e'],
            p['epsilon_B'], p['open_angle'], p['view_angle'],
            p['gaussian_cone'], p['jet_index'], p['T_ej'], p['L_ph'], p['acc_p'], p['eta_acc_p'], p['spectral_p']
        ]
        grb.setGRBafterglowParam(param)

        # -----------------------------------------------------------
        # Reverse shock microphysics
        # -----------------------------------------------------------
        p_RS = {
            'spectral_e_RS': 2.4,
            'epsilon_e_RS': 0.002,
            'fraction_e_RS': 1.,
            'eta_acc_e_RS': 1.,
            'epsilon_B_RS': 0.5,
            'acc_p_RS':     0.25,
            'eta_acc_p_RS': 1.,
            'spectral_p_RS': 2.,
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

        # -----------------------------------------------------------
        # Time and energy grids
        # -----------------------------------------------------------
        self.time_array = np.logspace(np.log10(1e1), np.log10(1e5), 6)

        self.energy_array_min = [
            1e9  / AMES.eV2Hz,
            1e11 / AMES.eV2Hz,
            1e0,
            1e3,
            1e9,
            1e11
        ]

        self.energy_array_max = self.energy_array_min.copy()

        os.makedirs("figure", exist_ok=True)
        os.makedirs("result", exist_ok=True)

    # ==================================================================
    #   Dynamical evolution
    # ==================================================================
    def calc_dynamic(self):
        print("\n=== Running dynamical evolution ===")
        grb.calcDynamic(self.time_array)

    # ==================================================================
    #   Spectrum
    # ==================================================================
    def calc_spectrum(self, T):
        print(f"\n=== Computing spectrum at T = {T:.2e} s ===")
        self.T = T

        grb.haveSSCSpec(True)
        grb.haveAttenuSSA(True)
        grb.haveAttenuGGSource(True)
        grb.haveReverseShock(True)
        grb.showInfo(True)

        # optional one-zone mode
        if False:
            grb.haveOneZone(True)
            grb.haveEdgeEffect(False)

        self.flux_vector = grb.Spectrum(self.T)

    # ==================================================================
    #   Plot spectrum comparison
    # ==================================================================
    def plot_spectrum(self):
        energy = np.array(s.getPhoton().getEnergy())

        fig, ax = plt.subplots(figsize=(7,5))

        # --- AMES components ---
        ax.plot(energy, self.flux_vector[0], '--', lw=2, label='syn (FS)')
        ax.plot(energy, self.flux_vector[1], '--', lw=2, label='ssc (FS)')
        ax.plot(energy, self.flux_vector[3], ':',  lw=2, label='syn (RS)')
        ax.plot(energy, self.flux_vector[4], ':',  lw=2, label='ssc (RS)')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e-8, 1e15])
        ax.set_ylim([1e-19, 1e-3])
        ax.set_xlabel('E [eV]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg\ cm^{-2}\ s^{-1}]$', fontsize=15)
        ax.legend(loc='best', fontsize=7)

        plt.tight_layout()
        plt.savefig('figure/spectrum.png')
        plt.savefig('figure/spectrum.pdf')
        plt.show()


# ======================================================================
#   RUN
# ======================================================================

# medium_type can be "ISM" or "wind"
g = GRB()

# Spectrum example at T = 1e5 s
g.calc_spectrum(1e5)
g.plot_spectrum()

