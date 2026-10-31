"""Shared engine setup for the two forward-shock examples.

The calling object must retain both Source and GRBafterglow: the C++ engine
holds a reference to Source. No calculations run when this module is imported.
"""
import AMES

# This order must match GRBafterglow::setGRBafterglowParam (18 entries).
PARAMETER_ORDER = (
    "z", "dl", "E_ej", "Gamma0", "n_ex", "A_star", "k_ex",
    "spectral_e", "epsilon_e", "fraction_e", "eta_acc_e", "epsilon_B",
    "open_angle", "view_angle", "gaussian_cone", "jet_index", "T_ej", "L_ph",
)


def build_model(parameters):
    """Return (source, engine), configured for forward-shock emission only.

    Momentum grids use the AMES native convention, p*c in eV. Grid resolution
    affects runtime and accuracy; check convergence before scientific use.
    """
    source = AMES.Source()
    source.getTarget().setMomentum(161, 1e-9, 1e7)
    source.getPhoton().setMomentum(251, 1e-9, 1e16)
    source.getElectron().setMomentum(361, 1e4, 1e16)
    source.InitSource()
    engine = AMES.GRBafterglow(source)
    engine.setGRBafterglowParam([parameters[key] for key in PARAMETER_ORDER])
    engine.haveSSCSpec(True)
    engine.haveAttenuSSA(True)
    engine.haveAttenuGGSource(True)
    engine.haveReverseShock(False)
    engine.showInfo(True)
    return source, engine
