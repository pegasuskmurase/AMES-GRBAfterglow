import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import AMES

s = AMES.Source()
grb = AMES.GRBAfterglow(s)

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
            'view_angle':          0.,  #
            'gaussian_cone':         0.,
            'jet_index':         0.,
        }

        self.p = p
        param = [p['z'], p['dl'], p['E_ej'], p['Gamma0'], p['n_ex'], p['k_ex'], p['spectral_e'],
                 p['epsilon_e'], p['fraction_e'], p['eta_acc_e'], p['epsilon_B'], p['open_angle'], p['view_angle'], p['gaussian_cone'], p['jet_index']]
        grb.setGRBAfterglowParam(param)

        t_min = np.log10(1e1)
        t_max = np.log10(1e5)
        self.time_array = np.logspace(t_min, t_max, 20)
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
        grb.haveSSCSpec(False)
        grb.haveAttenuGGSource(True)
        grb.haveOneZone(False)
        grb.Flux(self.time_array, self.energy_array_min, self.energy_array_max)

        #One-zone calculation
        grb.haveOneZone(True)
        folder = "result-onezone"
        grb.setOutputFolder(folder)
        if not os.path.exists(folder):
            try:
                os.makedirs(folder, 0o700)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        grb.Flux(self.time_array, self.energy_array_min, self.energy_array_max)

    def plot_spectrum(self):
        cm_subsection = np.linspace(0.2, 1, len(self.time_array))
        colors = [cm.YlOrBr(x) for x in cm_subsection]
        fig, ax = plt.subplots()
        data = np.loadtxt('result-onezone/spectrum.dat')
        num = len(s.getPhoton().getMomentum())
        for i, x in enumerate(self.time_array):
            idx1 = num * i
            idx2 = num * (i + 1)
            ax.plot(data[idx1:idx2, 0], data[idx1:idx2, 1], '--', c=colors[i], lw=2, label=str(x) + ' s, Onezone')
            ax.plot(data[idx1:idx2, 0], data[idx1:idx2, 2], '--', c=colors[i], lw=2, label=str(x) + ' s')

        data = np.loadtxt('result/spectrum.dat')
        num = len(s.getPhoton().getMomentum())
        for i, x in enumerate(self.time_array):
            idx1 = num * i
            idx2 = num * (i + 1)
            ax.plot(data[idx1:idx2, 0], data[idx1:idx2, 1], '-', c=colors[i], lw=2, label=str(x) + ' s, EATS')
            ax.plot(data[idx1:idx2, 0], data[idx1:idx2, 2], '-', c=colors[i], lw=2)

        for i, x in enumerate(self.time_array):
            self.spectrum_GS02(ax, self.time_array[i])
            self.spectrum_afterglowpy(ax, self.time_array[i])

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e-8, 1e15])
        ax.set_ylim([1e-18, 1e-3])
        ax.set_xlabel('E [eV]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg \ cm^{-2} \ s^{-1}]$', fontsize=15)
        # ax.legend(fontsize=7)
        plt.show()

    def plot_flux(self):
        fig, ax = plt.subplots()

        data = np.loadtxt('result/flux.dat')
        ax.plot(data[:, 0], data[:, 1], 'C0-')
        ax.plot(data[:, 0], data[:, 2], 'C1-')
        ax.plot(data[:, 0], data[:, 3], 'C2-')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1, 1e4])
        ax.set_ylim([1e-19, 1e-4])
        ax.set_xlabel('T [s]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg \ cm^{-2} \ s^{-1}]$', fontsize=15)
        ax.legend(fontsize=7)
        plt.show()


g = GRB()
#g.calc_flux()
g.plot_spectrum()
#g.plot_flux()
