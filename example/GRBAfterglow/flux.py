import numpy as np
import matplotlib.pyplot as plt
import os
import AMES
import afterglowpy as grbpy

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

        t_min = np.log10(1e4)
        t_max = np.log10(1e5)
        self.time_array = np.logspace(t_min, t_max, 1)
        self.energy_array_min = [1e9 / AMES.eV2Hz, 1e3, 1e11]
        self.energy_array_max = [1e9 / AMES.eV2Hz, 1e3, 1e11]

        grb.setOutputFolder("result")

    def calc_flux(self):
        grb.haveAttenuSSA(True)
        grb.haveEdgeEffect(False)
        grb.haveSSCSpec(True)
        grb.haveAttenuGGSource(True)
        grb.haveOneZone(False)
        grb.Flux(self.time_array, self.energy_array_min, self.energy_array_max)
        grb.setOutputFolder("result-onezone")
        grb.haveOneZone(True)
        grb.Flux(self.time_array, self.energy_array_min, self.energy_array_max)

    def spectrum_GS02(self, ax, time_obs):
        from GRB_GS02 import GS02
        g = GS02(self.p['E_ej'], self.p['n_ex'], self.p['epsilon_e'],
                 self.p['epsilon_B'], self.p['spectral_e'], self.p['dl'], self.p['z'])
        tdays = time_obs / 3600 / 24.
        energy = np.array(s.getPhoton().getEnergy())
        nu = energy * AMES.eV2Hz
        spectrum = g.gen_spectrum(nu, tdays)*1e-3/AMES.erg2Jy*nu
        ax.plot(energy, spectrum, '-.', c='C8', lw=3., label='Granot & Sari, 2002')

    def spectrum_afterglowpy(self, ax, time_obs):
        # For convenience, place arguments into a dict.
        Z = {'jetType':     grbpy.jet.TopHat,     # Top-Hat jet
             'specType':    0,                  # Basic Synchrotron Emission Spectrum
             'thetaObs':    self.p['view_angle'],   # Viewing angle in radians
             'E0':          self.p['E_ej'],  # Isotropic-equivalent energy in erg
             'thetaCore':   self.p['open_angle'],    # Half-opening angle in radians
             'n0':          self.p['n_ex'],    # circumburst density in cm^{-3}
             'p':           self.p['spectral_e'],    # electron energy distribution index
             'epsilon_e':   self.p['epsilon_e'],    # epsilon_e
             'epsilon_B':   self.p['epsilon_B'],   # epsilon_B
             'xi_N':        self.p['fraction_e'],    # Fraction of electrons accelerated
             'd_L':         self.p['dl'],  # Luminosity distance in cm
             'z':           self.p['z']}   # redshift
        # Space time points geometrically, from 10^3 s to 10^7 s
        # Calculate flux in a single X-ray band (all times have same frequency)
        energy = np.array(s.getPhoton().getEnergy())
        nu = energy * AMES.eV2Hz
        Fnu = grbpy.fluxDensity(time_obs, nu, **Z)
        spectrum = Fnu * 1e-3/AMES.erg2Jy*nu
        ax.plot(energy, spectrum, '-.', c='C5', lw=2., label='Afterglowpy')

    def plot_spectrum(self):
        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'r', 'g', 'b', 'c', 'y', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5',
                  'C6', 'C7', 'C8', 'C9', 'r', 'g', 'b', 'c', 'y', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'r', 'g', 'b', 'c', 'y']
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
g.calc_flux()
g.plot_spectrum()
#g.plot_flux()
