import numpy as np
import matplotlib.pyplot as plt
import os
import AMES
from GS02compare import spectrum_GS02 

s = AMES.Source()
s.InitSource()
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
            'T_ej':          1e4,  #
        }

        p_RS = {
            'spectral_e_RS':          2.4,
            'epsilon_e_RS':          0.002,
            'fraction_e_RS':          0.008,  #
            'eta_acc_e_RS':          1.,  #
            'epsilon_B_RS':          0.5,
            'acc_p_RS':          0.25,  # 0.6 for s_p = 2
            'eta_acc_p_RS':          1,  #
            'spectral_p_RS':          2.,
        }

        self.p = p
        param = [p['z'], p['dl'], p['E_ej'], p['Gamma0'], p['n_ex'], p['k_ex'], p['spectral_e'],
                 p['epsilon_e'], p['fraction_e'], p['eta_acc_e'], p['epsilon_B'], p['open_angle'], p['view_angle'], p['gaussian_cone'], p['jet_index'], p['T_ej']]
        grb.setGRBAfterglowParam(param)

        self.p_RS = p_RS
        param_RS = [p_RS['spectral_e_RS'], p_RS['epsilon_e_RS'], p_RS['fraction_e_RS'], p_RS['eta_acc_e_RS'],
                    p_RS['epsilon_B_RS']]
        grb.setGRBAfterglowRSParam(param_RS)

        t_min = np.log10(1e1)
        t_max = np.log10(1e5)
        self.time_array = np.logspace(t_min, t_max, 6)
        self.energy_array_min = [1e9 / AMES.eV2Hz, 1e11 / AMES.eV2Hz, 1e0, 1e3, 1e9, 1e11]
        self.energy_array_max = [1e9 / AMES.eV2Hz, 1e11 / AMES.eV2Hz, 1e0, 1e3, 1e9, 1e11]

    def calc_dynamic(self):
        grb.calcDynamic(self.time_array)

    def calc_flux(self):
        grb.haveSSCSpec(True)
        grb.haveAttenuSSA(True)
        grb.haveAttenuGGSource(True)
        flux_vector = grb.Flux(ED, syn, IC, gg, ph, self.time_array, self.energy_array_min, self.energy_array_max)
        self.flux_vector = flux_vector

    def calc_spectrum(self):
        self.T = 1e3
        grb.haveSSCSpec(True)
        grb.haveAttenuSSA(True)
        grb.haveAttenuGGSource(True)
        grb.haveReverseShock(True)

        have_onezone=False
        if have_onezone:
           grb.haveOneZone(have_onezone)
           grb.haveEdgeEffect(False)
        flux_vector = grb.Spectrum(self.T)
        self.flux_vector = flux_vector

    def plot_spectrum(self):
        fig, ax = plt.subplots()
        energy = np.array(s.getPhoton().getEnergy())

        spectrum_GS02(ax, self.p, self.T, energy, 'r')

        time = self.time_array
        ax.plot(energy, self.flux_vector[0], '--', lw=2, label='syn')
        ax.plot(energy, self.flux_vector[1], '--', lw=2, label='ssc')
        #ax.plot(energy, self.flux_vector[2], '--', lw=2)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e-8, 1e15])
        ax.set_ylim([1e-15, 1e-3])
        ax.set_xlabel('E [eV]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg \ cm^{-2} \ s^{-1}]$', fontsize=15)
        #ax.legend(fontsize=7)
        plt.savefig('spectrum.pdf')
        plt.show()

    def plot_flux(self):
        fig, ax = plt.subplots()

        #energy = s.getPhoton().getEnergy()
        time = self.time_array
        for j, x in enumerate(self.energy_array_min):
            dum = []
            for i in range(len(self.flux_vector)):
                dum.append(self.flux_vector[i][j])
            ax.plot(time, dum, '--', lw=2)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1, 1e4])
        ax.set_ylim([1e-19, 1e-4])
        ax.set_xlabel('T [s]', fontsize=15)
        ax.set_ylabel(r'Flux $[\rm erg \ cm^{-2} \ s^{-1}]$', fontsize=15)
        ax.legend(fontsize=7)
        plt.show()

g = GRB()
#g.calc_dynamic()
#g.calc_flux()
#g.plot_flux()
g.calc_spectrum()
g.plot_spectrum()
