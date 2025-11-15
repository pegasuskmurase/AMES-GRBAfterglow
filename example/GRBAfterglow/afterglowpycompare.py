import afterglowpy as grbpy
import matplotlib.pyplot as plt

def spectrum_afterglowpy(ax, p, time_obs, energy, colors):
    # For convenience, place arguments into a dict.
    Z = {'jetType':     grbpy.jet.TopHat,     # Top-Hat jet
         'specType':    0,                  # Basic Synchrotron Emission Spectrum
         'thetaObs':    p['view_angle'],   # Viewing angle in radians
         'E0':          p['E_ej'],  # Isotropic-equivalent energy in erg
         'thetaCore':   p['open_angle'],    # Half-opening angle in radians
         'n0':          p['n_ex'],    # circumburst density in cm^{-3}
         'p':           p['spectral_e'],    # electron energy distribution index
         'epsilon_e':   p['epsilon_e'],    # epsilon_e
         'epsilon_B':   p['epsilon_B'],   # epsilon_B
         'xi_N':        p['fraction_e'],    # Fraction of electrons accelerated
         'd_L':         p['dl'],  # Luminosity distance in cm
         'z':           p['z']}   # redshift
    # Space time points geometrically, from 10^3 s to 10^7 s
    # Calculate flux in a single X-ray band (all times have same frequency)
    #energy = np.array(s.getPhoton().getEnergy())
    nu = energy * AMES.eV2Hz
    Fnu = grbpy.fluxDensity(time_obs, nu, **Z)
    spectrum = Fnu * 1e-3/AMES.erg2Jy*nu
    ax.plot(energy, spectrum, '-.', c=colors, lw=2., label='Afterglowpy')
