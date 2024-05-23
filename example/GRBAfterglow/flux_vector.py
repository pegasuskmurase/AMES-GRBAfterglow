import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import AMES

s = AMES.Source()
grb = AMES.GRBAfterglow(s)

ED = AMES.ElectronDistribution(s)
syn = AMES.Synchrotron(s)
IC = AMES.InverseCompton(s)
gg = AMES.GammaGamma(s)
ph = AMES.Photonbackground(s)

class GRB:
    def __init__(self):
        p = {
            'z':           0.15,  # redshift
            'dl':         723 * AMES.Mpc,  # Luminosity distance [cm]
            'E_ej':          1.5e55,  # Isotropic-equivalent ejecta kinetic energy [erg]
            'Gamma0':          560,  # Initial Lorentz factor
            'n_ex':          0.6,  # external medium density [cm^-3]
            # 'n_ex':          1*3e35,  # external medium density [cm^-3]
            'k_ex':          0.,  # external medium density [cm^-3]
            'spectral_e':          2.2,  # external medium density spectral index
            'epsilon_e':          0.025,  #
            'fraction_e':          1.,  #
            'eta_acc_e':          1,  #
            'epsilon_B':          0.0006,
            'open_angle':          0.1,  #
        }


        self.p = p
        param = [p['z'], p['dl'], p['E_ej'], p['Gamma0'], p['n_ex'], p['k_ex'], p['spectral_e'],
                 p['epsilon_e'], p['fraction_e'], p['eta_acc_e'], p['epsilon_B'], p['open_angle']]
        grb.setGRBAfterglowParam(param)

        t_min = np.log10(1e2)
        t_max = np.log10(1e4)
        self.time_array = np.logspace(t_min, t_max, 5)
        self.energy_array_min = [1e9 / AMES.eV2Hz, 1e3, 1e11]
        self.energy_array_max = [1e9 / AMES.eV2Hz, 1e3, 1e11]

        folder = "result"
        grb.setOutputFolder(folder)
        if not os.path.exists(folder):
            try:
                os.makedirs(folder, 0o700)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

    def calc_flux(self):
        #Multi-zone calculation
        grb.haveAttenuSSA(True)
        grb.haveEdgeEffect(False)
        grb.haveSSCSpec(True)
        grb.haveAttenuGGSource(True)
        grb.haveOneZone(False)
        self.flux_vector = AMES.VecVecdouble()
        grb.Flux(self.flux_vector, ED, syn, IC, gg, ph, self.time_array, self.energy_array_min, self.energy_array_max)

    def plot_flux(self):
        fig, ax = plt.subplots()

        dum = []
        for i in range(len(self.flux_vector)):
            dum.append(self.flux_vector[i][0])
        ax.plot(self.time_array, dum, 'C0--', lw=2)

        dum = []
        for i in range(len(self.flux_vector)):
            dum.append(self.flux_vector[i][1])
        ax.plot(self.time_array, dum, 'C1--', lw=2)

        dum = []
        for i in range(len(self.flux_vector)):
            dum.append(self.flux_vector[i][2])
        ax.plot(self.time_array, dum, 'C2--', lw=2)

        #comparison
        data = np.loadtxt('result/flux.dat')
        ax.plot(data[:, 0], data[:, 1], 'C0-')
        ax.plot(data[:, 0], data[:, 2], 'C1-')
        ax.plot(data[:, 0], data[:, 3], 'C2-')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e2, 1e4])
        ax.set_ylim([1e-19, 1e-4])
        ax.set_xlabel('T [s]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg \ cm^{-2} \ s^{-1}]$', fontsize=15)
        ax.legend(fontsize=7)
        plt.show()


g = GRB()
g.calc_flux()
g.plot_flux()
