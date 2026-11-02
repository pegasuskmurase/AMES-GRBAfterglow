// src/module/main_bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "module/compton.hpp"
#include "module/particle.hpp"
#include "module/gammagamma.hpp"
#include "module/photonbackground.hpp"
#include "module/gammagamma.hpp"
#include "module/compton.hpp"
#include "module/source.hpp"
#include "module/synchrotron.hpp"
#include "module/utility.hpp"
#include "GRBAfterglow.hpp"

namespace py = pybind11;

// Forward declarations
class Particle;
class Param;
class Source;
class Target;
class Escape;
class Adiabatic;
// ... other forward declarations ...
// Particle base class bindings

void bind_particle(py::module_& m, const std::string& name) {
    py::class_<Particle>(m, name.c_str())
        .def(py::init<>())
        // Properly bind getEnergy() with vector conversion
        .def("getEnergy", &Particle::getEnergy)
        .def("getSpectraInj", &Particle::getSpectraInj)
	.def("getMomentum", &Particle::getMomentum)
        .def("getLosstime", &Particle::getLosstime)
        .def("getIntetime", &Particle::getIntetime)
        .def("Energy2Momentum", 
        static_cast<double (Particle::*)(double)>(&Particle::Energy2Momentum),
        py::arg("energy"))
        .def("Momentum2Energy", 
        static_cast<double (Particle::*)(double)>(&Particle::Momentum2Energy),
        py::arg("momentum"))

        .def("getSpectrum", 
            static_cast<std::vector<double> (Particle::*)()>(&Particle::getSpectrum),
            py::return_value_policy::copy)
        .def("setSpectrum", static_cast<void (Particle::*)(std::vector<double>)>(&Particle::setSpectrum),
            py::arg("spectrum"),
            "Set entire spectrum vector")
        .def("setSpectrum", static_cast<void (Particle::*)(int, double)>(&Particle::setSpectrum),
            py::arg("index"),
            py::arg("value"),
            "Set single spectrum value at index")
        .def("setSpectrumPL", &Particle::setSpectrumPL,
            py::arg("p_inj"),
            py::arg("p_cut"),
            py::arg("index"),
            py::arg("norm"),
            py::arg("norm_type") = false)
        .def("setSpectrumBPL", &Particle::setSpectrumBPL,
            py::arg("p_inj"),
            py::arg("p_break"),
            py::arg("p_cut"),
            py::arg("index_l"),
            py::arg("index_h"),
            py::arg("norm"),
            py::arg("norm_type") = false)

        // Overload 1: No arguments
        .def("setMomentum", py::overload_cast<>(&Particle::setMomentum))
        
        // Overload 2: size_t, double, double
        .def("setMomentum", py::overload_cast<size_t, double, double>(&Particle::setMomentum),
             py::arg("p_num"), py::arg("p_min"), py::arg("p_max"))
        
        // Overload 3: int index and double value
        .def("setMomentum", py::overload_cast<int, double>(&Particle::setMomentum),
             py::arg("i"), py::arg("value"));
}


void bind_param(py::module_& m) {
    py::class_<Param>(m, "Param")
        // Constructor
        .def(py::init<>())
        // Setters
        .def("setName", &Param::setName, py::arg("_name"))
        .def("setGeometry", &Param::setGeometry, py::arg("_geometry"))
        .def("setRedshift", 
            [](Param& self, double z) { self.setRedshift(z); },
            py::arg("_redshift"))

        .def("setLumDistance", &Param::setLumDistance, py::arg("_lum_distance"))
        .def("setAmbientDensity", &Param::setAmbientDensity, py::arg("_ambient_density"))
        .def("setEmissionRadius", &Param::setEmissionRadius, py::arg("_emission_radius"))
        .def("setDopplerFactor", &Param::setDopplerFactor, py::arg("_doppler_factor"))
        .def("setObsAngle", &Param::setObsAngle, py::arg("_obs_angle"))
        .def("setVelocity", &Param::setVelocity, py::arg("_velocity"))
        .def("setLorentzFactor", &Param::setLorentzFactor, py::arg("_lorentz_factor"))
        .def("setShellWidth", &Param::setShellWidth, py::arg("_shell_width"))
        .def("setMagStrength", &Param::setMagStrength, py::arg("_mag_strength"))
        // Getters (you'll need to add these to your Param class)
        .def("getName", &Param::getName)
        .def("getGeometry", &Param::getGeometry)
        .def("getDopplerFactor", &Param::getDopplerFactor)
        .def("getVelocity", &Param::getVelocity)
        .def("getRedshift", &Param::getRedshift)
        .def("getLumDistance", &Param::getLumDistance)
        .def("getEmissionRadius", &Param::getEmissionRadius)
        .def("getShellWidth", &Param::getShellWidth)
        ;
}

// Source class bindings
void bind_source(py::module_& m) {
    py::class_<Source>(m, "Source")
        .def(py::init<>())
        // Particle getters
        .def("clear", &Source::clear)
	.def("InitSource", &Source::InitSource)
        .def("getParam", &Source::getParam, pybind11::return_value_policy::reference_internal)
        .def("getTarget", &Source::getTarget, py::return_value_policy::reference_internal)
        .def("getPhoton", &Source::getPhoton, py::return_value_policy::reference_internal)
        .def("getPhotonSyn", &Source::getPhotonSyn, py::return_value_policy::reference_internal)
        .def("getPhotonIC", &Source::getPhotonIC, py::return_value_policy::reference_internal)
        .def("getPhotonPM", &Source::getPhotonPM, py::return_value_policy::reference_internal)
        .def("getPhotonPA", &Source::getPhotonPA, py::return_value_policy::reference_internal)
        .def("getElectron", &Source::getElectron, py::return_value_policy::reference_internal)
        .def("getElectronPrim", &Source::getElectronPrim, py::return_value_policy::reference_internal)
        .def("getElectronPM", &Source::getElectronPM, py::return_value_policy::reference_internal)
        .def("getElectronPA", &Source::getElectronPA, py::return_value_policy::reference_internal)
        .def("getNeutrino", &Source::getNeutrino, py::return_value_policy::reference_internal)
        .def("getNeutron", &Source::getNeutron, py::return_value_policy::reference_internal)
        .def("getProton", &Source::getProton, py::return_value_policy::reference_internal)
        .def("getMuonminusL", &Source::getMuonminusL, py::return_value_policy::reference_internal)
        .def("getMuonminusR", &Source::getMuonminusR, py::return_value_policy::reference_internal)
        .def("getMuonplusL", &Source::getMuonplusL, py::return_value_policy::reference_internal)
        .def("getMuonplusR", &Source::getMuonplusR, py::return_value_policy::reference_internal)
        .def("getPionminus", &Source::getPionminus, py::return_value_policy::reference_internal)
        .def("getPionplus", &Source::getPionplus, py::return_value_policy::reference_internal)
        .def("getPionzero", &Source::getPionzero, py::return_value_policy::reference_internal)
        .def("getNucleus", &Source::getNucleus, py::return_value_policy::reference_internal)
        .def("setNucleus", &Source::setNucleus,
          py::arg("mass_number"),
          py::arg("charge_number"),
          "Set the nucleus properties (mass number and charge number)");
}


// Physics process classes
void bind_physics_processes(py::module_& m) {
    py::class_<Synchrotron>(m, "Synchrotron")
        .def(py::init<Source&>())
        .def("Losstime", py::overload_cast<double>(&Synchrotron::Losstime),
        py::arg("mag_strength"),
        "Calculate loss time with magnetic field strength")
        .def("LosstimeProton", py::overload_cast<double>(&Synchrotron::LosstimeProton),
        py::arg("mag_strength"),
        "Calculate loss time with magnetic field strength")
        .def("Emissivity", static_cast<void (Synchrotron::*)()>(&Synchrotron::Emissivity));

    py::class_<InverseCompton>(m, "InverseCompton")
        .def(py::init<Source&>())
        .def("Intetime", static_cast<void (InverseCompton::*)()>(&InverseCompton::Intetime))
        .def("Losstime", static_cast<void (InverseCompton::*)()>(&InverseCompton::Losstime))
        .def("Emissivity", static_cast<void (InverseCompton::*)()>(&InverseCompton::Emissivity));

    py::class_<GammaGamma>(m, "GammaGamma")
        .def(py::init<Source&>())
        .def("Intetime", static_cast<void (GammaGamma::*)()>(&GammaGamma::Intetime));

    py::class_<ElectronDistribution>(m, "ElectronDistribution")
        .def(py::init<Source&>());

    // Photon background
    py::class_<Photonbackground>(m, "Photonbackground")
        .def(py::init<Source&>())
        .def("GreyBody", &Photonbackground::GreyBody,
            py::arg("T_ph"),
            py::arg("u_ph"))
        .def("Powerlaw", &Photonbackground::Powerlaw,
            py::arg("eps_min"),
            py::arg("eps_max"),
            py::arg("index"),
            py::arg("u_ph"))
        .def("BandFunction", &Photonbackground::BandFunction,
            py::arg("eps_min"),
            py::arg("eps_peak"),
            py::arg("eps_max"),
            py::arg("index_1"),
            py::arg("index_2"),
            py::arg("u_ph"))
        .def("CMB", &Photonbackground::CMB, 
	    py::arg("redshift"));

    // ... bind other physics processes similarly ...
}

// Utility class and constants
void bind_utility_and_constants(py::module_& m) {
    py::class_<Utility>(m, "Utility")
        .def(py::init<>())
        .def("Interpolate", 
            static_cast<double (Utility::*)(const std::vector<double>&, const std::vector<double>&, double, int)>(&Utility::Interpolate),
            py::arg("x"),
            py::arg("y"),
            py::arg("value"),
            py::arg("idx") = -1,
            R"pbdoc(
                Linear interpolation between data points
                [rest of docstring...]
            )pbdoc")
        .def("Interpolate2D", 
            static_cast<double (Utility::*)(const std::vector<double>&, const std::vector<double>&, const std::vector<std::vector<double>>&, double, double)>(&Utility::Interpolate2D),
            py::arg("x"),
            py::arg("y"),
            py::arg("xy"),
            py::arg("value_x"),
            py::arg("value_y"),
            R"pbdoc(
                Linear interpolation between data points
                [rest of docstring...]
            )pbdoc")
        .def("Integrate", 
            static_cast<double (Utility::*)(const std::vector<double>&, const std::vector<double>&, double, double)>(&Utility::Integrate),
            py::arg("x"),
            py::arg("y"),
            py::arg("x_begin"),
            py::arg("x_end")
            )
        .def("Integrate", 
            static_cast<double (Utility::*)(const std::vector<double>&, const std::vector<double>&)>(&Utility::Integrate),
            py::arg("x"),
            py::arg("y")
            );
    // Constants
    m.attr("PI") = PI;
    m.attr("c_cnst") = c_cnst;
    m.attr("erg2eV") = erg2eV;
    m.attr("eV2erg") = eV2erg;
    m.attr("e_charge") = e_charge;
    m.attr("eV2Hz") = eV2Hz;
    m.attr("erg2Jy") = erg2Jy;
    m.attr("Hz2eV") = Hz2eV;
    m.attr("mp") = mp;
    m.attr("me") = me;
    m.attr("sigma_T") = sigma_T;
    m.attr("proton_mass") = proton_mass;
    m.attr("electron_mass") = electron_mass;
    m.attr("k_B") = k_B;
    m.attr("M_sun") = M_sun;
    m.attr("Mpc") = Mpc;
    m.attr("Kpc") = Kpc;
    m.attr("pc") = pc;
    m.attr("yr") = yr;
    m.attr("pic_Life") = pic_Life;
    m.attr("mu_Life") = mu_Life;
    m.attr("pic_Mass") = pic_Mass;
    m.attr("mu_Mass") = mu_Mass;
    m.attr("days") = days;
}

void bind_grbafterglow(py::module_& m) {
    py::class_<GRBafterglow>(m, "GRBafterglow")
        .def(py::init<Source&>())
        // --- The Parameter Setters You Requested ---
        // We bind them with snake_case names for Pythonic style
        .def("setGRBafterglowParam", &GRBafterglow::setGRBafterglowParam, 
             "Set Forward Shock Parameters. Expects list/array of 21 floats.\n"
             "Order: z, dl, E_iso, Gamma0, n_ex, ...",
             py::arg("params"))

        .def("setGRBafterglowRSParam", &GRBafterglow::setGRBafterglowRSParam,
             "Set Reverse Shock Parameters. Expects list/array of 8 floats.\n"
             "Order: spec_e_RS, eps_e_RS, frac_e_RS, ...",
             py::arg("params"))

// --- MCMC Parameter Setter ---
        // Expects a single list of parameters (log10 or linear as defined in your C++)
        .def("set_params_mcmc", &GRBafterglow::setParamsFromMCMC,
             "Update physics parameters from MCMC walker.\n"
             "Input: List of doubles (E_iso, n_ex, p, eps_e, eps_B, ...)",
             py::arg("theta"))

        // --- Optimized Flux Calculator ---
        // Python usage: flux = grb.calculate_lightcurve(time_array, freq_array)
        .def("calculate_lightcurve", [](GRBafterglow &self, 
                                        const std::vector<double>& time, 
                                        const std::vector<double>& freq) {
            
            // 1. Create a local output vector
            std::vector<double> flux_out;
            
            // 2. Call the C++ function (which handles resizing and parallelization)
            // Note: Ensure your C++ implementation allows flux_out to be resized internally
            // or we pre-resize it here:
            flux_out.resize(time.size());
            
            self.CalculateLightCurveMCMC(flux_out, time, freq);
            
            // 3. Return vector
            // Pybind11/stl.h automatically converts std::vector<double> -> NumPy array
            return flux_out; 
        }, 
        "Calculate Flux for MCMC.\n"
        "Returns: NumPy array of Flux values corresponding to inputs.",
        py::arg("time"), py::arg("freq"))

        .def("showInfo", &GRBafterglow::showInfo) // showinfo
        .def("setOutputFolder", &GRBafterglow::setOutputFolder, py::arg("_output"))
        .def("haveSSCSpec", &GRBafterglow::haveSSCSpec)
        .def("haveOneZone", &GRBafterglow::haveOneZone)
        .def("haveAttenuSSA", &GRBafterglow::haveAttenuSSA)
        .def("haveAttenuGGSource", &GRBafterglow::haveAttenuGGSource)
        .def("haveReverseShock", &GRBafterglow::haveReverseShock)
        .def("Flux",
            [](GRBafterglow& self, 
               const std::vector<double>& time_array,
               const std::vector<double>& energy_array_min, const std::vector<double>& energy_array_max) 
            {
	        std::vector<std::vector<double>> flux_vector;
                self.Flux(flux_vector, time_array, energy_array_min, energy_array_max);
		return flux_vector;
            },
            py::arg("time_array"),
            py::arg("energy_array_min"), py::arg("energy_array_max"),
            "Calculate flux using various physical processes (modifies flux_vector in-place)"
        )

        // --- Optimized Flux Calculation for MCMC ---
        // Instead of passing an empty vector by reference (C++ style),
        // we write a Lambda to return the result as a NumPy array (Python style).
        .def("Flux1D", [](GRBafterglow &self, 
                                const std::vector<double>& time, 
                                const std::vector<double>& freq_min,
                                const std::vector<double>& freq_max
                            ) {
            
            // 1. Create output buffer
            std::vector<double> flux_out;
            
            // 2. Call the C++ logic (freq_min = freq_max for monochromatic)
            self.Flux1D(flux_out, time, freq_min, freq_max);
            
            // 3. Return vector (Pybind11 converts this to NumPy array automatically)
            return flux_out;
        }, 
        "Calculate Lightcurve. Returns Flux array.",
        py::arg("time"), py::arg("freq_min"), py::arg("freq_max"))

        // --- MCMC Helper (If you implemented the optimized setter) ---
        // .def("set_mcmc_params", &GRBafterglow::setParamsFromMCMC)

	.def("FluxOutput",
            [](GRBafterglow& self,
               const std::vector<double>& time_array,
               const std::vector<double>& energy_array_min,
               const std::vector<double>& energy_array_max)
            {
                self.FluxOutput(time_array, energy_array_min, energy_array_max);
            },
            py::arg("time_array"),
            py::arg("energy_array_min"),
            py::arg("energy_array_max"),
            R"pbdoc(
                Calculate flux for given time and energy ranges
            )pbdoc")
        .def("Spectrum",
            [](GRBafterglow& self, 
		    double T
               ) 
            {
	        std::vector<std::vector<double>> flux_vector;
                self.Spectrum(flux_vector, T);
		return flux_vector;
            },
            py::arg("T"),
            "Calculate flux using various physical processes (modifies flux_vector in-place)"
        );
}

PYBIND11_MODULE(AMES, m) {
    bind_param(m);
    bind_utility_and_constants(m);

    // Bind particle types first
    bind_particle(m, "Particle");
    // Then bind Source
    bind_source(m);

    // Bind physics processes
    bind_physics_processes(m);

    bind_grbafterglow(m);
}

