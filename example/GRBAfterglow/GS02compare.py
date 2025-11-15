from GRB_GS02 import GS02

def spectrum_GS02(ax, p, time_obs, energy, colors):
    g = GS02(p['E_ej'], p['n_ex'], p['epsilon_e'],
             p['epsilon_B'], p['spectral_e'], p['dl'], p['z'])
    tdays = time_obs / 3600 / 24.
    nu = energy * 2.4e14
    erg2Jy = 1e23;
    spectrum = g.gen_spectrum(nu, tdays)*1e-3/erg2Jy*nu
    ax.plot(energy, spectrum, '-', c=colors, lw=2., alpha=0.5, label='Granot & Sari, 2002')
