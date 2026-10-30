// AMES
// Copyright(C) 2024, licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file GRBafterglow.cpp
// \brief Calculate emissions GRB afterglow

#include "GRBAfterglow.hpp"
#include "module/rk5.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <sstream>

GRBafterglow::GRBafterglow(Source& s)
    : s(s), syn(s.GetSynchrotron()), IC(s.GetIC()), gg(s.GetGG()),
      ElectronD(s.GetElectronD()), ph(s.GetPH())
      {


    // init
    num_t = s.getTarget().getEnergy().size();
    num_p = s.getPhoton().getEnergy().size();
    num_e = s.getElectron().getEnergy().size();

    target_energy = s.getTarget().getEnergy();
    target_syn.resize(num_t, 0.0);
    target_prompt.resize(num_t, 0.0);
    target_tot.resize(num_t, 0.0);

    photon_energy = s.getPhoton().getEnergy();
    photon_syn.resize(num_p, 0.0);
    photon_ssc.resize(num_p, 0.0);
    photon_eic.resize(num_p, 0.0);

    photon_syn_temp.resize(num_p, 0.0);
    photon_syn_absorption_SSA.resize(num_p, 0.0);
    photon_syn_obs.resize(num_p, 0.0);
    photon_ssc_temp.resize(num_p, 0.0);
    photon_ssc_obs.resize(num_p, 0.0);
    photon_eic_temp.resize(num_p, 0.0);
    photon_eic_obs.resize(num_p, 0.0);
    photon_tot.resize(num_p, 0.0);
    photon_tot_diff.resize(num_p, 0.0);
    photon_syn_diff.resize(num_p, 0.0);
    photon_ssc_diff.resize(num_p, 0.0);
    photon_eic_diff.resize(num_p, 0.0);

    photon_losstime.resize(num_p, 0.0);

    photon_syn_RS.resize(num_p, 0.0);
    photon_ssc_RS.resize(num_p, 0.0);
    photon_eic_RS.resize(num_p, 0.0);
    photon_syn_RS_obs.resize(num_p, 0.0);
    photon_ssc_RS_obs.resize(num_p, 0.0);
    photon_eic_RS_obs.resize(num_p, 0.0);

    electron_energy = s.getElectron().getEnergy();
    }

void GRBafterglow::help() {
    std::cout << "----Welcome to use AMES-GRBafterglow------" << std::endl;
    std::cout << "Please see the following information for "
                 "proper use."
              << std::endl;
    std::cout << "If you have any questions please contact "
                 "Bing Theodore "
                 "Zhang<bing.zhang @yukawa.kyoto-u.ac.jp> "
                 "and Kohta "
                 "Murase<murase@psu.edu>."
              << std::endl;
}


// Internal helper to bind a variable to its constraints and name
struct ParamBinding {
    double* target;       // Pointer to the struct member
    Range range;          // Valid range
    const char* name;     // Human-readable name for logging
};

// The validation function now takes the name directly
double validate_and_log(double value, const Range& r, const char* name) {
    if (value < r.min || value > r.max) {
        std::cerr << "[Warning] Parameter '" << name << "' value " << value
                  << " is out of bounds [" << r.min << ", " << r.max 
                  << "]. Resetting to default: " << r.def << std::endl;
        return r.def;
    }
    return value;
}

// Helper to process a list of bindings against an input vector
void process_bindings(const std::vector<double>& inputs, const std::vector<ParamBinding>& bindings) {
    if (inputs.size() < bindings.size()) {
        std::cerr << "[Error] Input vector size (" << inputs.size() 
                  << ") matches insufficient parameters (expected " << bindings.size() << ")\n";
        return;
    }

    for (size_t i = 0; i < bindings.size(); ++i) {
        *bindings[i].target = validate_and_log(inputs[i], bindings[i].range, bindings[i].name);
    }
}


// ---------------------------------------------------------
// Implementation of Setters
// ---------------------------------------------------------

void GRBafterglow::setGRBafterglowParam(const std::vector<double>& v) {
// 1. Define Bindings
    // NOTE: Order must match the input vector 'v' exactly.
    std::vector<ParamBinding> bindings = {
        // ... [Previous Params] ...
        { &param.z,             {0.0, 10.0, 0.0},     "z (redshift)" },
        { &param.dl,            {1e10, 1e29, 1e26},   "dL (luminosity distance)" },
        { &param.E_ej,          {1e45, 1e57, 1e52},   "E_ej (Isotropic Energy)" },
        { &param.Gamma0,        {1.0, 1e4, 100.0},    "Gamma0 (Lorentz factor)" },
        { &param.n_ex,          {1e-4, 1e6, 1.0},     "n_ex (Density ISM)" },
        { &param.A_star,        {0, 1e5, 1.0},     "A_star (Wind param)" },
        { &param.k_ex,          {0.0, 3.0, 0.0},      "k_ex (Density profile)" },
        { &param.spectral_e,    {0.0, 5.0, 2.2},      "spectral_e (FS)" },
        { &param.epsilon_e,     {1e-8, 1.0, 0.1},     "epsilon_e (FS)" },
        { &param.fraction_e,    {1e-4, 1.0, 1.0},     "fraction_e (FS)" },
        { &param.eta_acc_e,     {1.0, 1e6, 1000.0},   "eta_acc_e (FS)" },
        { &param.epsilon_B,     {1e-8, 1.0, 1e-3},    "epsilon_B (FS)" },
        { &param.open_angle,    {0.01, 1.5, 0.1},     "open_angle" },
        { &param.view_angle,    {0.0, 3.0, 0.0},      "view_angle" },
        
        // FIX: Min value set to 1e-4 to prevent division by zero in structured jet code
        { &param.gaussian_cone, {0.0, 1.5, 0.1},     "gaussian_cone" },
        
        { &param.jet_index,     {0.0, 10.0, 0.0},     "jet_index" },
        { &param.T_ej,          {0.00001, 1e6, 10.0},     "T_ej (Duration)" },
        
        // New merged parameters
        { &param.L_ph,          {0.0, 1e60, 1e50},    "L_ph (Prompt)" },
    };

// 2. Validate Input Size
    if (v.size() != bindings.size()) {
        std::cerr << "Error: Input parameter size mismatch. Expected " 
                  << bindings.size() << ", got " << v.size() << ".\n";
        // Do not process to avoid out-of-bounds access
        return; 
    }

    // 3. Process
    process_bindings(v, bindings);

    // 4. Post-Processing Logic
    
    // Wind Environment Logic (Chevalier & Li 2000)
    // If k ~ 2 (Wind), n_ex implies the normalization constant A.
    // A = 3.0e35 * A_star [cm^-1]
    if (param.k_ex > 0.1 && param.A_star > 0.0) {
        param.n_ex = param.A_star * 3.0e35;
        // Optional: Log this override so the user isn't confused why n_ex changed
        // std::cout << "Info: Wind environment detected. n_ex set to normalization A = " << param.n_ex << "\n";
    }    process_bindings(v, bindings);

    if (param.k_ex > 0.0 && param.A_star > 0.0) {
        param.n_ex = param.A_star * 3.0e35;

        if (param.n_ex < 1e31 && param.n_ex > 1e40) {
            std::cerr << "[Warning] Computed n_ex (" << param.n_ex 
                      << ") from A_star may be unphysical for k_ex = " 
                      << param.k_ex << ". Please verify your parameters.\n";
        }   
    }

}

void GRBafterglow::setGRBafterglowRSParam(const std::vector<double>& v) {
    std::vector<ParamBinding> bindings = {
        { &param.spectral_e_RS, {0.0, 5.0, 2.2},      "spectral_e (RS)" },
        { &param.epsilon_e_RS,  {1e-8, 1.0, 0.1},     "epsilon_e (RS)" },
        { &param.fraction_e_RS, {1e-4, 1.0, 1.0},     "fraction_e (RS)" },
        { &param.eta_acc_e_RS,  {1.0, 1e6, 1000.0},   "eta_acc_e (RS)" },
        { &param.epsilon_B_RS,  {1e-8, 0.5, 1e-3},    "epsilon_B (RS)" },
    };

    process_bindings(v, bindings);
}

void GRBafterglow::setExtendedEmissionParam(std::vector<double> _param) {
    L_EE = _param[0];
    T_a_EE = _param[1];
    alpha_fl_EE = _param[2];
    alpha_st_EE = _param[3];
    energy_min_EE = _param[4];
    energy_peak_EE = _param[5];
    energy_max_EE = _param[6];
    spectra_l_EE = _param[7];
    spectra_h_EE = _param[8];
}

void GRBafterglow::setPlateauEmissionParam(std::vector<double> _param) {
    L_PL = _param[0];
    T_a_PL = _param[1];
    alpha_fl_PL = _param[2];
    alpha_st_PL = _param[3];
    energy_min_PL = _param[4];
    energy_peak_PL = _param[5];
    energy_max_PL = _param[6];
    spectra_l_PL = _param[7];
    spectra_h_PL = _param[8];
}

void GRBafterglow::setOutputFolder(std::string _output_folder) {
    output_folder = _output_folder;
}

void GRBafterglow::showInfo(bool _show_info) {
    show_info = _show_info;
} //
  //
void GRBafterglow::haveOneZone(bool _have_onezone) {
    have_onezone = _have_onezone;
} // One zone
void GRBafterglow::haveSSCSpec(bool _have_SSCSpec) {
    have_SSCSpec = _have_SSCSpec;
}
void GRBafterglow::haveEICSpec(bool _have_EICSpec) {
    have_EICSpec = _have_EICSpec;
}
void GRBafterglow::haveEdgeEffect(bool _have_edge_effect) {
    have_edge_effect = _have_edge_effect;
}
void GRBafterglow::haveAttenuSSA(bool _have_attenu_SSA) {
    have_attenu_SSA = _have_attenu_SSA;
}
void GRBafterglow::haveAttenuFF(bool _have_attenu_FF) {
    have_attenu_FF = _have_attenu_FF;
}
void GRBafterglow::haveAttenuGGSource(bool _have_attenu_GG_source) {
    have_attenu_GG_source = _have_attenu_GG_source;
}
void GRBafterglow::haveAttenuGGCosmic(bool _have_attenu_GG_cosmic) {
    have_attenu_GG_cosmic = _have_attenu_GG_cosmic;
}

void GRBafterglow::haveReverseShock(bool _have_RS) {
    have_RS = _have_RS;
}; // have RS


void GRBafterglow::Spectrum(std::vector<std::vector<double>>& flux_vector,
                            double T) {
                                 std::cout << "Parameter check before output:" << std::endl;
    std::cout << "Redshift (z): " << param.z << std::endl;
    std::cout << "Luminosity Distance (dL): " << param.dl << " cm" << std::endl;
    std::cout << "Isotropic Energy (E ej): " << param.E_ej << " erg" << std::endl;
    std::cout << "Initial Lorentz Factor (Gamma0): " << param.Gamma0 << std::endl;
    std::cout << "External Density (n ex): " << param.n_ex << " cm^-3" << std::endl;
    std::cout << "Wind Parameter (A star): " << param.A_star << std::endl;
    std::cout << "Density Profile Index (k ex): " << param.k_ex << std::endl;
    std::cout << "Electron Spectral Index (spectral e): " << param.spectral_e << std::endl;
    std::cout << "Electron Energy Fraction (epsilon e): " << param.epsilon_e << std::endl;
    std::cout << "Electron Acceleration Fraction (fraction e): " << param.fraction_e << std::endl;
    std::cout << "Electron Acceleration Efficiency (eta acc e): " << param.eta_acc_e << std::endl;
    std::cout << "Magnetic Energy Fraction (epsilon B): " << param.epsilon_B << std::endl;
    std::cout << "Jet Opening Angle (open angle): " << param.open_angle << " rad" << std::endl;
    std::cout << "Viewing Angle (view angle): " << param.view_angle << " rad" << std::endl;
    std::cout << "Gaussian Cone Width (gaussian cone): " << param.gaussian_cone << " rad" << std::endl;
    std::cout << "Jet Structure Index (jet index): " << param.jet_index << std::endl;
    std::cout << "Ejecta Duration (T ej): " << param.T_ej << " s" << std::endl;                          


    std::cout << "Electron RS Spectral Index (spectral e RS): " << param.spectral_e_RS << std::endl;
    std::cout << "Electron RS Energy Fraction (epsilon e RS): " << param.epsilon_e_RS << std::endl;
    std::cout << "Electron RS Acceleration Fraction (fraction e RS): " << param.fraction_e_RS << std::endl;
    std::cout << "Electron RS Acceleration Efficiency (eta acc e RS): " << param.eta_acc_e_RS << std::endl;
    std::cout << "Magnetic RS Energy Fraction (epsilon B RS): " << param.epsilon_B_RS << std::endl; 
    if (show_info) {
        std::cout << "Init jet " << std::endl;
    }
    if (have_onezone) {
    } else {
        InitJet(jet);
    }

    int num = 6;

    flux_vector.resize(num);
    for (size_t i = 0; i < num; i++) {
        flux_vector[i].assign(num_p, 0.0);
    }

    if (show_info) {
        std::cout << "Spectrum T " << T << " : Begin" << std::endl;
    }

    Spectrum(T);

    if (show_info) {
        std::cout << "Spectrum T " << T << " : Output" << std::endl;
    }

        for (size_t j = 0; j < num_p; j++) {
            flux_vector[0][j] = photon_syn[j];
            flux_vector[1][j] = photon_ssc[j];
            flux_vector[2][j] = photon_eic[j];
            flux_vector[3][j] = photon_syn_RS[j];
            flux_vector[4][j] = photon_ssc_RS[j];
            flux_vector[5][j] = photon_eic_RS[j];
        }

    if (show_info) {
        std::cout << "Spectrum T " << T << " : Finish" << std::endl;
    }
}

void GRBafterglow::setParamsFromMCMC(const std::vector<double>& theta) {
    // Example MCMC Vector Mapping:
    // 0: log10(E_iso)
    // 1: log10(n_ex)
    // 2: p (spectral index)
    // 3: log10(epsilon_e)
    // 4: log10(epsilon_B)
    // 5: log10(fraction_e) (xi_e)
    // 6: view_angle (rad)
    
    // 1. Update Parameters (Convert Log -> Linear if needed)
    param.E_ej       = pow(10.0, theta[0]);
    param.n_ex       = pow(10.0, theta[1]);
    param.spectral_e = theta[2];
    param.epsilon_e  = pow(10.0, theta[3]);
    param.epsilon_B  = pow(10.0, theta[4]);
    param.fraction_e = pow(10.0, theta[5]);
    param.view_angle = theta[6];

    // ... map others as needed ...

    // 2. Intelligent Jet Initialization
    // InitJet is expensive (memory allocation). Only run it if geometry changes.
    // Note: InitJet depends on Gamma0, open_angle, view_angle, and grid resolution.
    bool geometry_changed = (std::abs(param.view_angle - prev_view_angle) > 1e-6) ||
                            (std::abs(param.open_angle - prev_open_angle) > 1e-6) ||
                            (std::abs(param.Gamma0 - prev_Gamma0) > 1e-2);

    if (geometry_changed || !have_onezone) {
        InitJet(jet);
        
        // Update history
        prev_view_angle = param.view_angle;
        prev_open_angle = param.open_angle;
        prev_Gamma0     = param.Gamma0;
    }
}

void GRBafterglow::CalculateLightCurveMCMC(std::vector<double>& flux_vector,
                                           const std::vector<double>& time_array,
                                           const std::vector<double>& energy_array) {
    
    std::cout << "Redshift (z): " << param.z << std::endl;
    std::cout << "Luminosity Distance (dL): " << param.dl << " cm" << std::endl;
    std::cout << "Isotropic Energy (E ej): " << param.E_ej << " erg" << std::endl;
    std::cout << "Initial Lorentz Factor (Gamma0): " << param.Gamma0 << std::endl;
    std::cout << "External Density (n ex): " << param.n_ex << " cm^-3" << std::endl;
    std::cout << "Wind Parameter (A star): " << param.A_star << std::endl;
    std::cout << "Density Profile Index (k ex): " << param.k_ex << std::endl;
    std::cout << "Electron Spectral Index (spectral e): " << param.spectral_e << std::endl;
    std::cout << "Electron Energy Fraction (epsilon e): " << param.epsilon_e << std::endl;
    std::cout << "Electron Acceleration Fraction (fraction e): " << param.fraction_e << std::endl;
    std::cout << "Electron Acceleration Efficiency (eta acc e): " << param.eta_acc_e << std::endl;
    std::cout << "Magnetic Energy Fraction (epsilon B): " << param.epsilon_B << std::endl;
    std::cout << "Jet Opening Angle (open angle): " << param.open_angle << " rad" << std::endl;
    std::cout << "Viewing Angle (view angle): " << param.view_angle << " rad" << std::endl;
    std::cout << "Gaussian Cone Width (gaussian cone): " << param.gaussian_cone << " rad" << std::endl;
    std::cout << "Jet Structure Index (jet index): " << param.jet_index << std::endl;
    std::cout << "Ejecta Duration (T ej): " << param.T_ej << " s" << std::endl;                          


    std::cout << "Electron RS Spectral Index (spectral e RS): " << param.spectral_e_RS << std::endl;
    std::cout << "Electron RS Energy Fraction (epsilon e RS): " << param.epsilon_e_RS << std::endl;
    std::cout << "Electron RS Acceleration Fraction (fraction e RS): " << param.fraction_e_RS << std::endl;
    std::cout << "Electron RS Acceleration Efficiency (eta acc e RS): " << param.eta_acc_e_RS << std::endl;
    std::cout << "Magnetic RS Energy Fraction (epsilon B RS): " << param.epsilon_B_RS << std::endl; 

    // 1. Pre-allocate memory (Critical for speed)
    if (flux_vector.size() != time_array.size()) {
        flux_vector.resize(time_array.size());
    }

    // 2. Physics Loop
    // We assume InitJet was handled in the setter to avoid repeated checks here
    
    int idx = 0;
    for (double t : time_array) {
        
        // A. Evolve and Calculate Spectrum
        // NOTE: show_info must be FALSE during MCMC
        Spectrum(t); 

        // B. Sum Components (Vectorized operation)
        // Calculating 1.0/E once is faster than division in loop
        #pragma omp simd
        for (size_t i = 0; i < num_p; i++) {
            photon_tot[i] = photon_syn[i] + photon_ssc[i] + photon_eic[i];
            // photon_tot_diff is likely not needed if we just interpolate photon_tot
            // unless your Integrate function strictly requires it.
        }

        // C. Interpolate Flux at the specific observation frequency
        // energy_array[idx] corresponds to the observation band at time t
        double obs_energy = energy_array[idx];
        
        // Direct Interpolation (Faster than Integration for monochromatic observation)
        // Ensure utility.Interpolate handles log-log interpolation for spectra!
        double flux_val = utility.Interpolate(photon_energy, photon_tot, obs_energy);
        
        flux_vector[idx] = flux_val;
        idx++;
    }
}

void GRBafterglow::Flux1D(std::vector<double>& flux_vector,
                          const std::vector<double>& time_array,
                          const std::vector<double>& energy_array_min,
                          const std::vector<double>& energy_array_max) {

    if (have_onezone) {
    } else {
        InitJet(jet);
    }

    flux_vector.resize(time_array.size());

    int idx = 0;
    for (auto t : time_array) {
        if (show_info) {
            std::cout << "Flux: t = " << t << std::endl;
        }

        Spectrum(t);

        for (size_t i = 0; i < num_p; i++) {
            photon_tot[i] = photon_syn[i] + photon_ssc[i] + photon_eic[i];
            photon_tot_diff[i] = photon_tot[i] / photon_energy[i];
        }

        double dum;
        if (fabs(energy_array_min[idx] - energy_array_max[idx]) < 1e-5) {
            dum = utility.Interpolate(photon_energy, photon_tot,
                                      energy_array_min[idx]);
        } else {
            dum =
                utility.Integrate(photon_energy, photon_tot_diff,
                                  energy_array_min[idx], energy_array_max[idx]);
        }
        flux_vector[idx] = dum;

        idx++;
    }
}

void GRBafterglow::Flux(std::vector<std::vector<double>>& flux_vector,
                        const std::vector<double>& time_array,
                        const std::vector<double>& energy_array_min,
                        const std::vector<double>& energy_array_max) {

    if (have_onezone) {
    } else {
        InitJet(jet);
    }

    flux_vector.resize(time_array.size());
    for (size_t i = 0; i < time_array.size(); i++) {
        flux_vector[i].assign(energy_array_min.size(), 0.0);
    }

    int idx = 0;
    for (auto t : time_array) {
        if (show_info) {
            std::cout << "Flux: t = " << t << std::endl;
        }

        Spectrum(t);

        for (size_t i = 0; i < num_p; i++) {
            photon_tot[i] = photon_syn[i] + photon_ssc[i] + photon_eic[i];
            photon_tot_diff[i] = photon_tot[i] / photon_energy[i];
        }

        double dum;
        for (size_t i = 0; i < energy_array_min.size(); i++) {
            if (fabs(energy_array_min[i] - energy_array_max[i]) < 1e-5) {
                dum = utility.Interpolate(photon_energy, photon_tot,
                                          energy_array_min[i]);
            } else {
                dum =
                    utility.Integrate(photon_energy, photon_tot_diff,
                                      energy_array_min[i], energy_array_max[i]);
            }
            flux_vector[idx][i] = dum;
        }

        idx++;
    }
}

void GRBafterglow::FluxOutput(const std::vector<double>& time_array,
                              const std::vector<double>& energy_array_min,
                              const std::vector<double>& energy_array_max) {
    
                                 std::cout << "Parameter check before output:" << std::endl;
    std::cout << "Redshift (z): " << param.z << std::endl;
    std::cout << "Luminosity Distance (dL): " << param.dl << " cm" << std::endl;
    std::cout << "Isotropic Energy (E ej): " << param.E_ej << " erg" << std::endl;
    std::cout << "Initial Lorentz Factor (Gamma0): " << param.Gamma0 << std::endl;
    std::cout << "External Density (n ex): " << param.n_ex << " cm^-3" << std::endl;
    std::cout << "Wind Parameter (A star): " << param.A_star << std::endl;
    std::cout << "Density Profile Index (k ex): " << param.k_ex << std::endl;
    std::cout << "Electron Spectral Index (spectral e): " << param.spectral_e << std::endl;
    std::cout << "Electron Energy Fraction (epsilon e): " << param.epsilon_e << std::endl;
    std::cout << "Electron Acceleration Fraction (fraction e): " << param.fraction_e << std::endl;
    std::cout << "Electron Acceleration Efficiency (eta acc e): " << param.eta_acc_e << std::endl;
    std::cout << "Magnetic Energy Fraction (epsilon B): " << param.epsilon_B << std::endl;
    std::cout << "Jet Opening Angle (open angle): " << param.open_angle << " rad" << std::endl;
    std::cout << "Viewing Angle (view angle): " << param.view_angle << " rad" << std::endl;
    std::cout << "Gaussian Cone Width (gaussian cone): " << param.gaussian_cone << " rad" << std::endl;
    std::cout << "Jet Structure Index (jet index): " << param.jet_index << std::endl;
    std::cout << "Ejecta Duration (T ej): " << param.T_ej << " s" << std::endl;                          


    std::cout << "Electron RS Spectral Index (spectral e RS): " << param.spectral_e_RS << std::endl;
    std::cout << "Electron RS Energy Fraction (epsilon e RS): " << param.epsilon_e_RS << std::endl;
    std::cout << "Electron RS Acceleration Fraction (fraction e RS): " << param.fraction_e_RS << std::endl;
    std::cout << "Electron RS Acceleration Efficiency (eta acc e RS): " << param.eta_acc_e_RS << std::endl;
    std::cout << "Magnetic RS Energy Fraction (epsilon B RS): " << param.epsilon_B_RS << std::endl; 

    // 1. Safety check
    if (energy_array_min.size() != energy_array_max.size()) {
        std::cerr << "Error: Energy array sizes do not match." << std::endl;
        return;
    }

    if (!have_onezone) {
        InitJet(jet);
    }

    // 2. Open files with check
    auto open_file = [&](const std::string& name) {
        std::ofstream fs(output_folder + "/" + name);
        if (!fs.is_open()) std::cerr << "Error opening " << name << "\n";
        fs << std::scientific << std::setprecision(6); // Set generic precision
        return fs;
    };

    std::ofstream f_spec = open_file("spectrum.dat");
    std::ofstream f_flux = open_file("flux.dat");
    
    // Only open RS files if needed to avoid creating empty files
    std::ofstream f_spec_RS, f_flux_RS;
    if (have_RS) {
        f_spec_RS = open_file("spectrum_RS.dat");
        f_flux_RS = open_file("flux_RS.dat");
    }

    if (show_info) {
        f_flux << "# t (s)  Flux_Band_1  Flux_Band_2 ...\n";
        if(have_RS) f_flux_RS << "# t (s)  Flux_Band_1  Flux_Band_2 ...\n";
    }

    // ----------------------------------------------------------------
    // Helper Lambda: Calculates Flux for specific bands
    // Reduces code duplication between FS and RS
    // ----------------------------------------------------------------
    auto write_band_flux = [&](std::ofstream& stream, double t, 
                               const std::vector<double>& p_tot, 
                               const std::vector<double>& p_tot_diff) {
        stream << t; 
        for (size_t i = 0; i < energy_array_min.size(); i++) {
            double val = 0.0;
            // Check for Monochromatic (narrow band) vs Integrated (broad band)
            if (std::abs(energy_array_min[i] - energy_array_max[i]) < 1e-10 * energy_array_min[i]) {
                val = utility.Interpolate(photon_energy, p_tot, energy_array_min[i]);
            } else {
                val = utility.Integrate(photon_energy, p_tot_diff, 
                                      energy_array_min[i], energy_array_max[i]);
            }
            stream << " " << val;
        }
        stream << "\n"; // Use \n, not endl
    };

    // ----------------------------------------------------------------
    // Main Time Loop
    // ----------------------------------------------------------------
    for (auto t : time_array) {

        if (show_info) {
            std::cout << "Flux Output: t = " << t << std::endl;
        }
        
        // Compute Physics
        Spectrum(t); 

        // --- 1. Forward Shock (FS) Spectrum Output ---
        
        // Add a header/separator so you can distinguish time steps in the file
        f_spec << "# Time = " << t << "\n"; 

        for (size_t i = 0; i < num_p; i++) {
            // Recalculate totals
            photon_tot[i] = photon_syn[i] + photon_ssc[i] + photon_eic[i];
            
            // Pre-calculate diffs (avoid division in multiple places)
            double inv_E = 1.0 / photon_energy[i];
            photon_tot_diff[i] = photon_tot[i] * inv_E;
            
            // Write Spectrum components
            f_spec << photon_energy[i] << " " 
                   << photon_syn[i] << " " 
                   << photon_ssc[i] << " "
                   << photon_eic[i] << "\n"; 
        }
        f_spec << "\n\n"; // Gnuplot "Index" separator (essential for plotting)

        // --- 2. Forward Shock (FS) Flux Output ---
        write_band_flux(f_flux, t, photon_tot, photon_tot_diff);

        // --- 3. Reverse Shock (RS) Output ---
        if (have_RS) {
            f_spec_RS << "# Time = " << t << "\n";
            
            for (size_t i = 0; i < num_p; i++) {
                // Update total arrays with RS data for the helper function
                photon_tot[i] = photon_syn_RS[i] + photon_ssc_RS[i] + photon_eic_RS[i];
                
                double inv_E = 1.0 / photon_energy[i];
                photon_tot_diff[i] = photon_tot[i] * inv_E;

                f_spec_RS << photon_energy[i] << " " 
                          << photon_syn_RS[i] << " "
                          << photon_ssc_RS[i] << " " 
                          << photon_eic_RS[i] << "\n";
            }
            f_spec_RS << "\n\n"; // Gnuplot separator

            write_band_flux(f_flux_RS, t, photon_tot, photon_tot_diff);
        }
    }

    // Files close automatically when destructors are called.
    // No manual flush needed; OS handles it efficiently.
}

void GRBafterglow::SpectrumOneZone(double T) {
    double shell_thickness_factor = 4. * (3. - param.k_ex);
        double Gamma, radius;
        EvolutionThinShellAnalytic(T, Gamma, radius);
        double beta = sqrt(1. - 1. / Gamma / Gamma);
        double dynamical_timescale = Gamma * T / (1 + param.z);
        double adiabatic_timescale = dynamical_timescale;
        double Delta =
            radius / shell_thickness_factor / Gamma; // comoving frame
        double photon_escape_timescale = Delta / c_cnst;

        double mag_strength =
            sqrt(32 * PI * mp * c_cnst * c_cnst * param.epsilon_B *
                 ExternalDensity(radius) * Gamma * (Gamma - 1));
        s.getParam().setMagStrength(mag_strength);

        double electron_inj = (4. * Gamma) * param.fraction_e *
                              ExternalDensity(radius) / dynamical_timescale;

        s.getElectron().setSpectrumPL(MomentumInj(Gamma, mag_strength, radius, beta, param.epsilon_e, param.fraction_e, param.spectral_e, param.eta_acc_e),
                                      MomentumMax(mag_strength, radius, Gamma, param.eta_acc_e),
                                      param.spectral_e, electron_inj, false);
        ElectronD.IterationSolutionSteadyState(
            syn, IC, false, dynamical_timescale, photon_escape_timescale,
            mag_strength);
        target_syn = s.getTarget().getSpectrum();

        syn.Emissivity(mag_strength, photon_energy, photon_syn);

        if (have_SSCSpec) {
            IC.Emissivity(target_energy, target_syn, photon_energy, photon_ssc);
        }

        if (have_attenu_SSA) {
            syn.SSAOpticalDepth(dynamical_timescale * c_cnst, mag_strength);
            syn.Attenuation(photon_syn);
        }

        if (have_attenu_GG_source) {
            gg.Losstime(target_energy, target_syn, photon_energy,
                        photon_losstime);
            for (size_t i = 0; i < num_p; i++) {
                s.getPhoton().setOptdepth(i, photon_losstime[i] *
                                                 dynamical_timescale);
            }
            gg.Attenuation(photon_ssc);
        }

        LorentzBoost(photon_energy, photon_syn, Gamma);
        LorentzBoost(photon_energy, photon_ssc, Gamma);
        FluxNorm(photon_energy, photon_syn, radius, Delta);
        FluxNorm(photon_energy, photon_ssc, radius, Delta);

        if (have_edge_effect) {
            if (1. / Gamma > param.open_angle) {
                for (size_t k = 0; k < num_p; k++) {
                    photon_syn[k] *=
                        param.open_angle * param.open_angle * Gamma * Gamma;
                    photon_ssc[k] *=
                        param.open_angle * param.open_angle * Gamma * Gamma;
                }
            }
        }
}

void GRBafterglow::Spectrum(double T) {

    if (have_onezone) {
        SpectrumOneZone(T);
        return;
    }

    // 1. Evolve Dynamics (Global calculation, sequential or internally parallelized)
    EvolutionThickShell(jet, T);

    // 2. Reset Global Accumulators
    photon_syn.assign(num_p, 0.0);
    photon_ssc.assign(num_p, 0.0);
    photon_eic.assign(num_p, 0.0);

    // -------------------------------------------------------------------------
    // OPENMP PARALLEL REGION
    // -------------------------------------------------------------------------
    //#pragma omp parallel
    {
        // --- A. Setup Thread-Local Objects ---
        // We create a private copy of the Source object for this thread.
        // This prevents threads from overwriting each other's magnetic field/spectra.
        Source local_s = this->s; 
        
        // Get references to the local modules from the local Source
        Synchrotron& local_syn = local_s.GetSynchrotron();
        InverseCompton& local_IC = local_s.GetIC();
        ElectronDistribution& local_ElectronD = local_s.GetElectronD();
        GammaGamma& local_gg = local_s.GetGG();
        
        // Thread-local Accumulators (to avoid race conditions during summation)
        std::vector<double> local_photon_syn(num_p, 0.0);
        std::vector<double> local_photon_ssc(num_p, 0.0);
        std::vector<double> local_photon_eic(num_p, 0.0);

        // Pre-allocate scratchpads to avoid reallocation inside loops
        std::vector<double> local_electron_N(num_e);
        std::vector<double> local_target_syn(num_t);
        std::vector<double> local_photon_temp(num_p); 
        std::vector<double> local_absorption_SSA(num_p);
        std::vector<double> local_losstime(num_p);

        // ---------------------------------------------------------------------
        // INTERNAL LAMBDA (Defined inside parallel region to capture local vars)
        // ---------------------------------------------------------------------
        auto calculate_local_physics = [&](double radius, double Gamma, double mu, 
                                           double current_theta, double theta_bin, 
                                           double phi_bin, double t_dyn, double t_adi, 
                                           bool is_on_axis) {
            
            // --- A. Geometry & Kinematics ---
            double beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));
            
            double gamma_hat = (4.0 * Gamma + 1.0) / (3.0 * Gamma);
            double num = (Gamma + 1) * pow(gamma_hat * (Gamma - 1) + 1, 2);
            double den = gamma_hat * (2 - gamma_hat) * (Gamma - 1) + 2;
            double gamma_sh = sqrt(num / den);
            double beta_sh  = sqrt(1.0 - 1.0 / (gamma_sh * gamma_sh));

            double doppler_factor = 1.0 / (Gamma * (1.0 - beta * mu));
    
            double shell_thickness_factor = 4. * (3. - param.k_ex);
            double delta_s = radius / shell_thickness_factor / (Gamma * Gamma) / std::abs(mu - beta_sh);

            double solid_angle = (is_on_axis) ? (2 * PI * theta_bin * sin(current_theta)) 
                                              : (theta_bin * sin(current_theta) * phi_bin); 
            
            double EATS_factor = solid_angle * radius * radius * std::abs(mu - beta_sh) / (1.0 - beta_sh * mu) * (1.0 + param.z) / (param.dl * param.dl);

            // --- B. Magnetic Field & Electrons ---
            double n_ext = ExternalDensity(radius);
            double mag_strength = sqrt(32 * PI * mp * c_cnst * c_cnst * param.epsilon_B * n_ext * Gamma * (Gamma - 1));
            
            local_s.getParam().setMagStrength(mag_strength); // Modify LOCAL source

            double photon_escape_timescale = radius / shell_thickness_factor / Gamma / c_cnst * 2;
            double electron_inj = (4.0 * Gamma) * param.fraction_e * n_ext / t_dyn; 

            // Debug info (Only allow master thread to print to avoid garbled output)
            //if (show_info && omp_get_thread_num() == 0 && current_theta < 1e-3) {
                 // Print only once or for specific condition
            //}

            local_s.getElectron().setSpectrumPL(
                MomentumInj(Gamma, mag_strength, radius, beta, param.epsilon_e, param.fraction_e, param.spectral_e, param.eta_acc_e),
                MomentumMax(mag_strength, radius, Gamma, param.eta_acc_e),
                param.spectral_e, electron_inj, false);

            // --- C. Electron Cooling & Emission ---
            
            // 1. Prepare Inputs
            // Get Spectrum returns a COPY or Reference. 
            // Since we modified local_s above, this gets the specific PL for this patch.
            std::vector<double> local_electron_Q = local_s.getElectron().getSpectrum(); 
            
            // Prepare Target (Prompt) - Assuming target_prompt is global read-only or calculated here
            // If PromptPhotonDensity modifies 's', it must be called on 'local_s'
            std::vector<double> local_input_target = target_prompt; 
            if(have_EICSpec) {
                // If this function modifies s.target, call it on local_s/local_ph
                // PromptPhotonDensity(local_ph, T, radius, Gamma); 
                // local_input_target = local_s.getTarget().getSpectrum();
            }

            // 2. Call Stateless Solver
            local_ElectronD.IterationSolution(local_syn, local_IC, have_EICSpec, 
                                        t_dyn, t_adi, photon_escape_timescale, mag_strength,
                                        local_electron_Q,    // Input Q
                                        local_input_target,  // Input Target
                                        local_electron_N,    // Output N
                                        local_target_syn);   // Output Total Target
            // 3. Emission

            const auto& e_energy = local_s.getElectron().getEnergy();
            
            // Synchrotron
            local_syn.Emissivity(mag_strength, e_energy, local_electron_N, photon_energy, local_photon_temp);
            for (size_t k = 0; k < num_p; k++) local_photon_temp[k] *= photon_energy[k]; // Convert to Energy Flux

            // Store Synchrotron Component (Temporary)
            std::vector<double> syn_flux_patch = local_photon_temp;

            // SSC
            std::vector<double> ssc_flux_patch(num_p, 0.0);
            if (have_SSCSpec) {
                local_IC.Emissivity(target_energy, local_target_syn, e_energy, local_electron_N, photon_energy, local_photon_temp);
                for (size_t k = 0; k < num_p; k++) ssc_flux_patch[k] = local_photon_temp[k] * photon_energy[k];
            }

            // EIC
            std::vector<double> eic_flux_patch(num_p, 0.0);
            if (have_EICSpec) {
                 double cos_theta_sc = (mu - beta) / (1 - beta * mu);
                 local_IC.EmissivityAIC(cos_theta_sc, target_energy, target_prompt, photon_energy, local_photon_temp);
                 for (size_t k = 0; k < num_p; k++) eic_flux_patch[k] = local_photon_temp[k] * photon_energy[k];
            }
            // --- E. Absorption & Accumulation ---
            
            // Calculate Absorption Coefficients (thread-local buffers)
            if (have_attenu_SSA) {
                local_syn.SSAAbsorptionCoeff(mag_strength, local_absorption_SSA);
            }
            if (have_attenu_GG_source) {
                local_gg.Losstime(target_energy, target_tot, photon_energy, local_losstime);
            }

            double shift_ratio = (1.0 + param.z) / doppler_factor;
            double geom_factor = (doppler_factor * doppler_factor) / (4 * PI) * delta_s; 

            for (size_t k = 0; k < num_p; k++) {
               double shifted_E = photon_energy[k] * shift_ratio;
               double att = 1.0;
           
               if (have_attenu_SSA) {
                   double alpha = utility.Interpolate(photon_energy, local_absorption_SSA, shifted_E) / doppler_factor;
                   double tau = alpha * delta_s;
                   att *= (tau > 1e-4) ? (1. - exp(-tau)) / tau : 1.0 - tau / 2.0; 
               }
               
               if (have_attenu_GG_source) {
                   double alpha = utility.Interpolate(photon_energy, local_losstime, shifted_E) / doppler_factor / c_cnst;
                   double tau = alpha * delta_s;
                   att *= (tau > 1e-4) ? (1. - exp(-tau)) / tau : 1.0 - tau / 2.0;
               }

               double final_weight = geom_factor * att; 
               double bin_energy_factor = EATS_factor * photon_energy[k] * eV2erg;

               // Accumulate to THREAD-LOCAL vectors
               auto get_val = [&](const std::vector<double>& source) {
                   return utility.Interpolate(photon_energy, source, shifted_E) * final_weight * bin_energy_factor;
               };

               local_photon_syn[k] += get_val(syn_flux_patch);
               if (have_SSCSpec) local_photon_ssc[k] += get_val(ssc_flux_patch);
               if (have_EICSpec) local_photon_eic[k] += get_val(eic_flux_patch);
            }
        }; 
        // --- END LAMBDA ---

        // ---------------------------------------------------------------------
        // PARALLEL FOR LOOP
        // Schedule dynamic because off-axis loop is heavier than on-axis
        // ---------------------------------------------------------------------
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < jet.theta.size(); i++) {
            
            if (param.view_angle == 0) {
                // --- ON AXIS ---
                calculate_local_physics(
                    jet.radius[i][0], jet.Gamma[i][0], jet.mu[i][0], 
                    jet.theta[i], jet.theta_bin[i], 0.0, 
                    jet.dynamical_timescale[i][0], jet.adiabatic_timescale[i][0], 
                    true
                );
            } else {
                // --- OFF AXIS ---
                for (size_t j = 0; j < jet.phi.size(); j++) {
                    calculate_local_physics(
                        jet.radius[i][j], jet.Gamma[i][j], jet.mu[i][j], 
                        jet.theta[i], jet.theta_bin[i], jet.phi_bin[j], 
                        jet.dynamical_timescale[i][j], jet.adiabatic_timescale[i][j], 
                        false
                    );
                }
            }
        } // End For

        // ---------------------------------------------------------------------
        // REDUCTION (Merge thread-local results to global)
        // ---------------------------------------------------------------------
        //#pragma omp critical
        {
            for(size_t k=0; k<num_p; ++k) {
                photon_syn[k] += local_photon_syn[k];
                photon_ssc[k] += local_photon_ssc[k];
                photon_eic[k] += local_photon_eic[k];
            }
        }

    } // --- END PARALLEL REGION ---

    if (have_RS) {
        ReverseShock(jet); // Note: ReverseShock also needs parallelization similar to above!
    }
}

void GRBafterglow::ReverseShock(Jet& jet) {

    // 1. Setup Constants & Decay Indices
    double Delta0 = param.T_ej / (1.0 + param.z) * c_cnst; 
    double R_dec = DecelerationRadius();
    
    double g_med = (param.k_ex > 0.1) ? 2.0 : 0.0; 
    
    double r_index_density, r_index_B;
    bool is_thick_shell = (Delta0 > R_dec / (2.0 * param.Gamma0 * param.Gamma0));

    if (!is_thick_shell) {
        // Thin Shell case
        Delta0 = R_dec / (2.0 * param.Gamma0 * param.Gamma0);
        r_index_density = -6.0 * (3.0 + g_med) / 7.0; 
        r_index_B       = -4.0 * (3.0 + g_med) / 7.0; 
    } else {
        // Thick Shell case
        r_index_density = -13.0 / 2.0;
        r_index_B       = -13.0 / 3.0;
    }

    // 2. Reset Global Accumulators
    photon_syn_RS.assign(num_p, 0.0);
    photon_ssc_RS.assign(num_p, 0.0);
    photon_eic_RS.assign(num_p, 0.0);

    // -------------------------------------------------------------------------
    // OPENMP PARALLEL REGION
    // -------------------------------------------------------------------------
    //#pragma omp parallel
    {
        // --- A. Setup Thread-Local Objects ---
        // Clone the Source object to ensure thread safety for internal state (spectra, B-field)
        Source local_s = this->s; 
        
        // Get references to local modules
        Synchrotron& local_syn = local_s.GetSynchrotron();
        InverseCompton& local_IC = local_s.GetIC();
        ElectronDistribution& local_ElectronD = local_s.GetElectronD();
        GammaGamma& local_gg = local_s.GetGG();
        
        // Thread-local Accumulators
        std::vector<double> local_photon_syn_RS(num_p, 0.0);
        std::vector<double> local_photon_ssc_RS(num_p, 0.0);
        std::vector<double> local_photon_eic_RS(num_p, 0.0);
        
        // Reusable Scratchpads
        std::vector<double> local_electron_N(num_e);
        std::vector<double> local_target_syn(num_t);
        std::vector<double> local_photon_temp(num_p);
        std::vector<double> local_absorption_SSA(num_p);
        std::vector<double> local_losstime(num_p);

        // ---------------------------------------------------------------------
        // INTERNAL LAMBDA: Calculates Radiation (Captures LOCAL objects)
        // ---------------------------------------------------------------------
        auto calculate_patch_radiation = [&](double mag_B, double delta_s_local, 
                                             double doppler_local, double EATS_fac,
                                             std::vector<double>& electron_dist_N,   
                                             std::vector<double>& target_photon_density 
                                            ) {
            
            // 1. Synchrotron
            // Use local_s properties (x_e) and local modules
            const auto& x_e = local_s.getElectron().getEnergy();
            
            local_syn.Emissivity(mag_B, x_e, electron_dist_N, photon_energy, local_photon_temp);
            for (size_t k = 0; k < num_p; k++) local_photon_temp[k] *= photon_energy[k]; // Convert to Flux

            // Store pure synchrotron flux for summation
            std::vector<double> syn_flux_pure = local_photon_temp;

            // 2. SSC
            std::vector<double> ssc_flux_pure(num_p, 0.0);
            if (have_SSCSpec) {
                local_IC.Emissivity(target_energy, target_photon_density, x_e, electron_dist_N, photon_energy, local_photon_temp);
                for (size_t k = 0; k < num_p; k++) ssc_flux_pure[k] = local_photon_temp[k] * photon_energy[k];
            }

            // D. Absorption & Accumulation
            double shift_ratio = (1.0 + param.z) / doppler_local;
            
            // Calculate Absorption Coefficients locally
            if (have_attenu_SSA) local_syn.SSAAbsorptionCoeff(mag_B, local_absorption_SSA);
            if (have_attenu_GG_source) local_gg.Losstime(target_energy, target_photon_density, photon_energy, local_losstime);

            for (size_t k = 0; k < num_p; k++) {
                double shifted_E = photon_energy[k] * shift_ratio;
                double att = 1.0;

                // SSA
                if (have_attenu_SSA) {
                    double alpha = utility.Interpolate(photon_energy, local_absorption_SSA, shifted_E) / doppler_local;
                    double tau = alpha * delta_s_local;
                    att *= (tau > 1e-4) ? (1. - exp(-tau))/tau : 1.0 - tau/2.0;
                }
                // GG
                if (have_attenu_GG_source) {
                     double alpha = utility.Interpolate(photon_energy, local_losstime, shifted_E) / doppler_local / c_cnst;
                     double tau = alpha * delta_s_local;
                     att *= (tau > 1e-4) ? (1. - exp(-tau))/tau : 1.0 - tau/2.0;
                }

                // Interpolation Helper
                auto get_obs = [&](const std::vector<double>& src) {
                    return utility.Interpolate(photon_energy, src, shifted_E) * (doppler_local * doppler_local) / (4 * PI) * att * delta_s_local;
                };

                // Accumulate into THREAD-LOCAL vectors
                double common = EATS_fac * photon_energy[k] * eV2erg;
                local_photon_syn_RS[k] += common * get_obs(syn_flux_pure);
                if(have_SSCSpec) local_photon_ssc_RS[k] += common * get_obs(ssc_flux_pure);
            }
        };

        // ---------------------------------------------------------------------
        // PARALLEL LOOP
        // ---------------------------------------------------------------------
        //#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < jet.theta.size(); i++) {
            
            size_t phi_limit = (param.view_angle < 1e-5) ? 1 : jet.phi.size();

            for (size_t j = 0; j < phi_limit; j++) {
                
                // --- 1. Geometry & Kinematics ---
                double Gamma = jet.Gamma_RS[i][j];
                double R_curr = jet.radius_RS[i][j];
                double mu = jet.mu[i][j];

                if (Gamma <= 1.0) continue; 

                double beta = sqrt(1.0 - 1.0/(Gamma*Gamma));
                
                double gamma_hat = (4.0 * Gamma + 1.0) / (3.0 * Gamma);
                double gamma_sh = sqrt(((Gamma+1)*(gamma_hat*(Gamma-1)+1)*(gamma_hat*(Gamma-1)+1)) / 
                                       (gamma_hat*(2-gamma_hat)*(Gamma-1)+2));
                double beta_sh = sqrt(1.0 - 1.0/(gamma_sh*gamma_sh));

                double doppler_factor = 1.0 / (Gamma * (1.0 - beta * mu));
                
                double solid_angle = (param.view_angle < 1e-5) 
                                   ? 2 * PI * jet.theta_bin[i] * sin(jet.theta[i])
                                   : jet.theta_bin[i] * sin(jet.theta[i]) * jet.phi_bin[j];

                double EATS_factor = solid_angle * R_curr * R_curr * std::abs(mu - beta_sh) / (1.0 - beta_sh * mu) * (1.0 + param.z) / (param.dl * param.dl);

                // --- 2. Determine State ---
                double mag_strength = 0.0;
                double delta_s = 0.0;
                bool active = false;
                double electron_inj = 0.0;
                double photon_escape = 0.0;

                double N0_iso = jet.E_ej[i] / (jet.Gamma0[i] * mp * c_cnst * c_cnst);

                // Prepare Inputs for Solver
                // Use GLOBAL prompt target as base (Read-Only)
                std::vector<double> local_input_target = target_prompt; 

                if (jet.mass_ejecta[i][j] <= jet.M_ej[i]) {
                    // PHASE 1: SHELL CROSSING
                    double Gamma34 = 0.5 * (Gamma / jet.Gamma0[i] + jet.Gamma0[i] / Gamma);
                    
                    if (Gamma34 > 1.0001) {
                        active = true;
                        
                        double n4 = N0_iso / (4 * PI * jet.Gamma0[i] * Delta0 * R_curr * R_curr);
                        double n3 = (4.0 * Gamma34 + 3.0) * n4; 
                        double e3 = (Gamma34 - 1.0) * n3 * mp * c_cnst * c_cnst; 
                        
                        double width_comoving = N0_iso * (jet.mass_ejecta[i][j]/jet.M_ej[i]) 
                                                / (4 * PI * R_curr * R_curr * n3);
                        delta_s = (width_comoving / Gamma) / std::abs(mu - beta_sh);

                        mag_strength = sqrt(8 * PI * param.epsilon_B_RS * e3);
                        local_s.getParam().setMagStrength(mag_strength); // Set LOCAL B-field

                        photon_escape = width_comoving / c_cnst * 2.0;
                        electron_inj = param.fraction_e_RS * n3 / jet.dynamical_timescale[i][j];

                        // Set LOCAL Electron Spectrum
                        local_s.getElectron().setSpectrumPL(
                            MomentumInj(Gamma34, mag_strength, R_curr, beta, param.epsilon_e_RS, param.fraction_e_RS, param.spectral_e_RS, param.eta_acc_e_RS),
                            MomentumMax(mag_strength, R_curr, Gamma34, param.eta_acc_e_RS),
                            param.spectral_e_RS, electron_inj, false);
                    }

                } else {    
                    // PHASE 2: POST-CROSSING
                    active = true;
                    double R_cross = jet.radius_cross[i][j];
                    double ratio = R_curr / R_cross; 
                    
                    if(R_cross > 0) {
                        double Gamma34_cross = 0.5 * (jet.Gamma_cross[i][j]/jet.Gamma0[i] + jet.Gamma0[i]/jet.Gamma_cross[i][j]);
                        double n4_cross = N0_iso / (4*PI*jet.Gamma0[i]*Delta0*R_cross*R_cross);
                        double n3_cross = (4*Gamma34_cross+3)*n4_cross;
                        double e3_cross = (Gamma34_cross-1)*n3_cross*mp*c_cnst*c_cnst;
                        double B_cross = sqrt(8*PI*param.epsilon_B_RS*e3_cross);
                        
                        mag_strength = B_cross * pow(ratio, r_index_B);
                        local_s.getParam().setMagStrength(mag_strength); // Set LOCAL B-field

                        double width_cross = N0_iso / (4*PI*R_cross*R_cross*n3_cross);
                        double width_comoving = width_cross * pow(ratio, 1.0); 
                        delta_s = (width_comoving / Gamma) / std::abs(mu - beta_sh);
                        
                        double beta = sqrt(1.0 - 1.0/(jet.Gamma_cross[i][j]*jet.Gamma_cross[i][j]));
                        double Gamma_min_new = MomentumInj(Gamma34_cross, B_cross, R_cross, beta, param.epsilon_e_RS, param.fraction_e_RS, param.spectral_e_RS, param.eta_acc_e_RS  ) * pow(ratio, -1.0);
                        double n_norm_new = (param.fraction_e_RS * n3_cross) * pow(ratio, r_index_density);
                        
                        photon_escape = width_comoving / c_cnst * 2.0;
                        
                        // Set LOCAL Electron Spectrum
                        local_s.getElectron().setSpectrumPL(
                            Gamma_min_new,
                            MomentumMax(mag_strength, R_curr, Gamma34_cross, param.eta_acc_e_RS), 
                            param.spectral_e_RS, 
                            n_norm_new / jet.dynamical_timescale[i][j], 
                            false
                        );
                    } else {
                        active = false;
                    }
                }

                // --- 3. Run Calculations if Active ---
                if (active) {
                    // Get the Injection Source Function from the LOCAL object
                    // This reflects the setSpectrumPL call above
                    std::vector<double> local_electron_Q = local_s.getElectron().getSpectrum(); 

                    // Solve Evolution
                    // IMPORTANT: Ensure IterationSolution writes to local_electron_N and local_target_syn
                    local_ElectronD.IterationSolution(
                        local_syn, local_IC, have_EICSpec, 
                        jet.dynamical_timescale[i][j], jet.adiabatic_timescale[i][j], photon_escape, mag_strength,
                        local_electron_Q,    // Input Q
                        local_input_target,  // Input Target (Prompt)
                        local_electron_N,    // Output N (Buffer)
                        local_target_syn     // Output Target (Buffer)
                    );

                    // Compute Emission and Accumulate
                    calculate_patch_radiation(mag_strength, delta_s, doppler_factor, EATS_factor,
                                            local_electron_N, local_target_syn);
                }

            } // End Phi Loop
        } // End Theta Loop

        // ---------------------------------------------------------------------
        // REDUCTION (Merge thread-local results to global)
        // ---------------------------------------------------------------------
        //#pragma omp critical
        {
            for(size_t k=0; k<num_p; ++k) {
                photon_syn_RS[k] += local_photon_syn_RS[k];
                photon_ssc_RS[k] += local_photon_ssc_RS[k];
                if(have_EICSpec) photon_eic_RS[k] += local_photon_eic_RS[k];
            }
        }

    } // --- END PARALLEL REGION ---
}

double GRBafterglow::MomentumInj(double Gamma, double mag_strength,
                                 double radius, double beta, double epsilon_e, double fraction_e, double spectral_e, double eta_acc_e) {
    // 1. Safety check for coasting/non-relativistic phase
    if (Gamma <= 1.00001) {
        return s.getElectron().Energy2Momentum(e_mass); // Return rest mass momentum (v~0)
    }

    double ave_kinetic_energy = proton_mass * (epsilon_e / fraction_e) * (Gamma - 1.0);
    double kinetic_energy = 0.0;
    double p = spectral_e;
    double gamma_min = 1.0;

    // 3. Determine Minimum Lorentz Factor (gamma_min) based on spectral index p
    if (p > 2.001) {
        // Standard case: Energy dominated by low-energy electrons
        // gamma_min = (p-2)/(p-1) * <gamma>
        kinetic_energy = (p - 2.0) / (p - 1.0) * ave_kinetic_energy;
        //std::cout << "Standard spectrum case: p=" << p << ", gamma_min approx = " << kinetic_energy / e_mass << std::endl;

    } else if (p < 1.999) {
        // Hard spectrum case: Energy dominated by high-energy electrons (gamma_max)
        // This is an approximation. Ideally, one solves for gamma_min iteratively 
        // using gamma_max, but this analytical approx is common in fast codes.
        kinetic_energy = (2.0 - p) / (p - 1.0) * ave_kinetic_energy;
        
        // Note: For very hard spectra (p < 2), gamma_min might fall below 1 
        // if not careful, though the pre-factor (mp/me) usually keeps it high.

    } else {
        // Critical case: p ~= 2.0
        // Energy is distributed equally per decade (Logarithmic divergence)
        // We need gamma_max to normalize.
        
        // Recalculate derived vars needed for Max
        double beta = std::sqrt(1.0 - 1.0 / (Gamma * Gamma));

        double gamma_max = MomentumMax(mag_strength, radius, Gamma, eta_acc_e) / e_mass; // Pass Gamma!
        double gamma_ave = proton_mass / e_mass * param.epsilon_e /
        param.fraction_e * (Gamma - 1.0);
	double log_r = std::log(gamma_max / gamma_ave);
        kinetic_energy =
        gamma_ave / (log_r + std::log(log_r)) * e_mass;
    }

    return s.getElectron().Energy2Momentum(kinetic_energy + e_mass);
}

double GRBafterglow::MomentumMax(double mag_strength, double radius, double Gamma, double eta_acc_e) {
    // Determine velocity from Gamma for transport calculations
    double beta = 0.0;
    if (Gamma > 1.0) {
        beta = std::sqrt(1.0 - 1.0 / (Gamma * Gamma));
    } else {
        return s.getElectron().Energy2Momentum(e_mass); // No acceleration
    }

    // 1. Synchrotron Limit (t_acc == t_syn)
    // Limits the energy when radiative losses equal energy gain.
    // gamma_max ~ 1/sqrt(B)
    // t_acc = eta_acc_e * gamma * me * c / (e * B)
    // t_syn = 6 * PI * me * c / (sigmaT * B^2 * gamma)
    // Setting t_acc = t_syn and solving for gamma gives:
    double gamma_max_syn = std::sqrt(6.0 * PI * e_charge / (eta_acc_e * sigmaT * mag_strength));

    // 2. Confinement / Adiabatic Limit (t_acc == t_dyn)
    //t_dyn = radius / (Gamma * beta * c_cnst);
    // Setting t_acc = t_dyn and solving for gamma gives:
    
    double gamma_max_conf = (e_charge * mag_strength * radius) / 
                            (eta_acc_e * me * c_cnst * c_cnst * Gamma * beta);

    // 3. Return the strictest limit
    double gamma_max = std::min(gamma_max_syn, gamma_max_conf);

    // Convert to momentum
    // (gamma * mc^2) -> momentum

    return s.getElectron().Energy2Momentum(gamma_max * e_mass);
}

void GRBafterglow::LorentzBoost(const std::vector<double>& photon_energy,
                                std::vector<double>& photon_spectrum,
                                double Gamma) {
    int num_p = photon_energy.size();
    std::vector<double> photon_temp(num_p, 0.0);
    for (size_t i = 0; i < num_p; i++) {
        photon_temp[i] = photon_energy[i] * photon_spectrum[i];
    }
    for (size_t i = 0; i < num_p; i++) {
        photon_spectrum[i] =
            (1 + param.z) * Gamma *
            utility.Interpolate(photon_energy, photon_temp,
                                (1 + param.z) * photon_energy[i] / Gamma);
    }
}

double GRBafterglow::ExternalDensity(double radius) {
    if (param.k_ex > 0) {
        return param.n_ex * pow(radius, -param.k_ex);
    } else {
        return param.n_ex;
    }
}

double GRBafterglow::DecelerationRadius() {
    double M_ej = param.E_ej / (param.Gamma0 * c_cnst * c_cnst);
    return pow((3 - param.k_ex) * M_ej /
                   (4 * PI * param.n_ex * mp * param.Gamma0),
               1 / (3 - param.k_ex));
}

double GRBafterglow::DecelerationRadius(double E_ej, double Gamma) {
    double M_ej = E_ej / (Gamma * c_cnst * c_cnst);
    return pow((3 - param.k_ex) * M_ej / (4 * PI * param.n_ex * mp * Gamma),
               1 / (3 - param.k_ex));
}

void GRBafterglow::EvolutionThinShellAnalytic(double T, double& Gamma,
                                              double& radius) {
    Gamma = 43.6 * exp(0.125 * log(param.E_ej / 1e53 / param.n_ex)) *
            exp(-0.375 * log(T / 1e3)) *
            exp(0.375 * log(1 + param.z)); // Sari98

    if (Gamma < 1.00001) {
        Gamma = 1.00001;
    }
    if (Gamma > param.Gamma0) {
        Gamma = param.Gamma0;
    }
    radius = 4. * Gamma * Gamma * c_cnst * T / (1 + param.z); // Sari98
}



void GRBafterglow::FluxNorm(const std::vector<double>& photon_energy,
                            std::vector<double>& photon_spectrum, double radius,
                            double Delta) {
    double flux_norm = eV2erg * (4.0 * PI * radius * radius * Delta) /
                       (4 * PI * param.dl * param.dl);
    int num_p = photon_energy.size();
    for (size_t i = 0; i < num_p; i++) {
        photon_spectrum[i] *= flux_norm * photon_energy[i];
    }
}

double GRBafterglow::MomentumMaxDyn(double mag_strength, double radius,
                                    double Gamma, double beta, double Z,
                                    double eta_acc) {
    // maximum energy, t_dyn = t_acc
    // Ref: Eq.
    return 1. / eta_acc * (radius / Gamma / beta) * Z * e_charge *
           mag_strength * erg2eV;
}

void GRBafterglow::PromptPhotonDensity(Photonbackground& ph, double T,
                                       double radius, double Gamma) {
    double u_EE = ExtendedEmissionLuminosity(T) /
                  (4 * PI * radius * radius * Gamma * Gamma * c_cnst) * erg2eV;
    ph.BrokenPowerlaw(energy_min_EE / Gamma, energy_peak_EE / Gamma,
                      energy_max_EE / Gamma, spectra_l_EE, spectra_h_EE, u_EE);
    std::vector<double> target_EE = s.getTarget().getSpectrum();

    double u_PL = PlateauEmissionLuminosity(T) /
                  (4 * PI * radius * radius * Gamma * Gamma * c_cnst) * erg2eV;
    ph.BrokenPowerlaw(energy_min_PL / Gamma, energy_peak_PL / Gamma,
                      energy_max_PL / Gamma, spectra_l_PL, spectra_h_PL, u_PL);
    std::vector<double> target_PL = s.getTarget().getSpectrum();

    for (size_t k = 0; k < num_t; k++) {
        s.getTarget().setSpectrum(k, target_EE[k] + target_PL[k]);
        target_prompt[k] = target_EE[k] + target_PL[k];
    }
}

double GRBafterglow::ExtendedEmissionLuminosity(double T) {
    double dum = 0.0;
    if (T < T_a_EE) {
        dum = L_EE * pow(T / T_a_EE, -alpha_fl_EE);
    } else {
        dum = L_EE * pow(T / T_a_EE, -alpha_st_EE);
    }
    return dum;
}

double GRBafterglow::PlateauEmissionLuminosity(double T) {
    double dum = 0.0;
    if (T < T_a_PL) {
        dum = L_PL * pow(T / T_a_PL, -alpha_fl_PL);
    } else {
        dum = L_PL * pow(T / T_a_PL, -alpha_st_PL);
    }
    return dum;
}

void GRBafterglow::InitJet(Jet& jet) {
    double sin_v = sin(param.view_angle);
    double cos_v = cos(param.view_angle);
    double sin_theta;
    double theta_min = 0.001 / param.Gamma0;
    double theta_max = param.open_angle;

    jet.theta.resize(jet.angular_num);
    jet.theta_bin.resize(jet.angular_num);

    /*
    for (size_t i = 0; i < jet.angular_num; i++) {
      jet.theta[i] = 0. + param.open_angle * i / (jet.angular_num - 1);
      jet.theta_bin[i] = param.open_angle / jet.angular_num;
    }
    */
    if (param.view_angle > 0) {
        for (size_t i = 0; i < jet.angular_num; i++) {
            jet.theta[i] = 0. + param.open_angle * i / (jet.angular_num - 1);
            jet.theta_bin[i] = param.open_angle / jet.angular_num;
            // jet.theta[i] = theta_min + (1 - theta_min) * i /
            // (jet.angular_num
            // - 1); jet.theta_bin[i] = (1 - theta_min) / jet.angular_num;
        }
    } else {

        // Calculate the logarithmic spacing factor
        double log_factor = (std::log10(theta_max) - std::log10(theta_min)) /
                            (jet.angular_num - 1);

        // Generate the logarithmically spaced grid

        for (int i = 0; i < jet.angular_num; ++i) {
            jet.theta[i] = std::pow(10, std::log10(theta_min) + i * log_factor);
        }
        double theta_l =
            std::pow(10, std::log10(theta_min) + (-1) * log_factor);
        double theta_r =
            std::pow(10, std::log10(theta_min) + jet.angular_num * log_factor);
        for (int i = 0; i < jet.angular_num; ++i) {
            if (i == 0) {
                jet.theta_bin[i] = 0.5 * (jet.theta[i + 1] - theta_l);
            } else if (i == jet.angular_num - 1) {
                jet.theta_bin[i] = 0.5 * (theta_r - jet.theta[i - 1]);
            } else {
                jet.theta_bin[i] = 0.5 * (jet.theta[i + 1] - jet.theta[i - 1]);
            }
        }
    }

    jet.phi.resize(jet.phi_num);
    jet.phi_bin.resize(jet.phi_num);
    for (size_t i = 0; i < jet.phi_num; i++) {
        jet.phi[i] = 0. + 2 * PI * i / (jet.phi_num - 1);
        jet.phi_bin[i] = 2 * PI / jet.phi_num;
    }

    jet.Gamma0.resize(jet.angular_num);
    jet.M_ej.resize(jet.angular_num);
    jet.E_ej.resize(jet.angular_num);
    jet.Gamma.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.radius.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.Gamma_RS.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.radius_RS.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.mass_ejecta.resize(jet.angular_num,
                           std::vector<double>(jet.phi_num, 0.));
    jet.radius_cross.resize(jet.angular_num,
                            std::vector<double>(jet.phi_num, 0.));
    jet.Gamma_cross.resize(jet.angular_num,
                           std::vector<double>(jet.phi_num, 0.));
    jet.mass_cross.resize(jet.angular_num,
                          std::vector<double>(jet.phi_num, 0.));
    jet.mu.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.dynamical_timescale.resize(jet.angular_num,
                                   std::vector<double>(jet.phi_num, 0.));
    jet.adiabatic_timescale.resize(jet.angular_num,
                                   std::vector<double>(jet.phi_num, 0.));
    jet.dynamical_timescale_cross.resize(jet.angular_num,
                                         std::vector<double>(jet.phi_num, 0.));
    jet.adiabatic_timescale_cross.resize(jet.angular_num,
                                         std::vector<double>(jet.phi_num, 0.));

    for (size_t i = 0; i < jet.theta.size(); i++) {
        for (size_t j = 0; j < jet.phi.size(); j++) {
            jet.mu[i][j] =
                sin(jet.theta[i]) * sin(param.view_angle) * cos(jet.phi[j]) +
                cos(jet.theta[i]) * cos(param.view_angle);
        }
    }
}

void GRBafterglow::EvolutionThickShell(const std::vector<double>& time_array,
                                       std::vector<std::vector<double>>& jet_dyn) {
    
    if (time_array.empty()) {
        std::cerr << "Warning: Empty time array passed to EvolutionThickShell." << std::endl;
        return;
    }

    const int n_vars = 11;
    jet_dyn.clear(); 
    jet_dyn.reserve(time_array.size()); 

    // Constants
    double epsilon = 0.0; 
    double C_BM = (17.0 - 4.0 * param.k_ex) > 0 ? 
                  (9.0 - 2.0 * param.k_ex) / (17.0 - 4.0 * param.k_ex) * (1.0 - epsilon) : 0.5;
    double C_ST = 16.0 / 9.0 * (3.0 - param.k_ex) * (5.0 - param.k_ex) /
                  (10.0 - 3.0 * param.k_ex);

    double Delta = param.T_ej / (1.0 + param.z) * c_cnst; 
    double R_dec = DecelerationRadius(); 
    double thin_shell_limit = R_dec / (2.0 * param.Gamma0 * param.Gamma0);
    double RS_index = (7.0 - 2.0 * param.k_ex) / 2.0; 

    if (Delta < thin_shell_limit) {
        Delta = thin_shell_limit;
        RS_index = 2.0; 
    }

    double solid_angle_factor = (1.0 - cos(param.open_angle)) / 2.0;
    double E_ej_true = param.E_ej * solid_angle_factor;
    double Gamma0 = param.Gamma0;
    double M_ej_true = E_ej_true / (Gamma0 * c_cnst * c_cnst);
    double beta0 = sqrt(1.0 - 1.0 / (Gamma0 * Gamma0));

    // Time conversion
    std::vector<double> time_array_source = time_array;
    double z_factor = 1.0 + param.z;
    for (auto& t : time_array_source) t /= z_factor;

    // --- OPTIMIZATION START ---
    // Parallelize the loop over angular segments
    // Note: The original code overwrites jet_dyn in every iteration. 
    // We assume this is intended for 1D cases, but to be safe with OpenMP, 
    // we calculate in parallel and write sequentially.

    #pragma omp parallel
    {
        // Thread-local variables to avoid race conditions and reallocation overhead
        std::vector<double> y(n_vars);
        std::vector<double> dydx(n_vars);
        std::vector<double> yout(n_vars);
        std::vector<std::vector<double>> local_jet_dyn;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < jet.theta.size(); i++) {
            
            // Re-initialize state for this thread/angle
            double R_start = 1e10;
            std::fill(y.begin(), y.end(), 0.0);
            
            y[0] = R_start;
            y[1] = 0.0;
            y[2] = Gamma0;
            y[3] = param.open_angle;
            y[4] = 0.0; 
            y[5] = M_ej_true / 1.0e6; // Seed mass
            y[6] = 0.0;
            y[7] = Gamma0; 
            y[8] = R_start / (c_cnst * beta0 * Gamma0); 
            y[9] = y[8]; 
            y[10] = R_start; 

            DynamicThickShell dyn(M_ej_true, C_BM, C_ST, param.n_ex, param.k_ex, epsilon,
                                  Gamma0, beta0, Delta, RS_index, param.open_angle);

            double x_start = 0.0; // Start Time
            dyn(x_start, y, dydx); // Prime derivatives

            rk5<DynamicThickShell> evol(dyn, y, dydx, yout, x_start, n_vars);
            
            // Perform heavy integration
            local_jet_dyn = evol.integrate(time_array_source);

            // Critical section for writing output
            // (Only needed if jet.theta.size() > 1 and we want the last result, 
            // or if we were aggregating. Based on original code, it just assigns).
            #pragma omp critical
            {
                jet_dyn = local_jet_dyn;
            }
        }
    }
    // --- OPTIMIZATION END ---
}

void GRBafterglow::EvolutionThickShell(Jet& jet, double T) {
    const int n_vars = 11;
    
    // Constants
    double epsilon = 0.0;
    double C_BM = (17.0 - 4.0 * param.k_ex) > 0 ? 
                  (9.0 - 2.0 * param.k_ex) / (17.0 - 4.0 * param.k_ex) : 0.5;
                  
    double C_ST = 16.0 / 9.0 * (3.0 - param.k_ex) * (5.0 - param.k_ex) /
                  (10.0 - 3.0 * param.k_ex);

    // --- OPTIMIZATION START ---
    // This is the heavy lifting function called by Spectrum(T).
    // We parallelize the loop over the jet surface (theta).
    
    #pragma omp parallel
    {
        // allocate thread-local vectors ONCE
        std::vector<double> y(n_vars);
        std::vector<double> dydx(n_vars);
        std::vector<double> yout(n_vars);
        std::vector<double> shock_dynamic;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < jet.theta.size(); i++) {
            
            // 1. Calculate Deceleration Radius & Shell Thickness
            double R_dec = DecelerationRadius(); 
            double Delta = param.T_ej / (1.0 + param.z) * c_cnst;
            double RS_index = (7.0 - 2.0 * param.k_ex) / 2.0;

            double min_Delta = R_dec / (2.0 * param.Gamma0 * param.Gamma0);
            if (Delta < min_Delta) {
                Delta = min_Delta;
                RS_index = 2.0; 
            }

            // 2. Structured Jet Profile
            double E_ej_local = param.E_ej;
            double Gamma_local = param.Gamma0;
            double theta = jet.theta[i];

            if (param.gaussian_cone > 0) {
                double ratio = theta / param.gaussian_cone;
                if (param.jet_index > 0) {
                    double factor = pow(1.0 + ratio * ratio, -param.jet_index / 2.0);
                    E_ej_local *= factor;
                    Gamma_local *= factor;
                } else {
                    E_ej_local *= exp(-0.5 * ratio * ratio);
                }
            }
            
            E_ej_local *= (2 * PI * (1 - cos(param.open_angle)) / (4 * PI));
            
            double M_ej_local = E_ej_local / (Gamma_local * c_cnst * c_cnst);
            double beta_local = sqrt(1.0 - 1.0 / (Gamma_local * Gamma_local));

            // Thread-safe write to unique index i
            jet.Gamma0[i] = Gamma_local;
            jet.M_ej[i] = M_ej_local;
            jet.E_ej[i] = E_ej_local;

            // 3. Initialize ODE System
            DynamicThickShell dyn(M_ej_local, C_BM, C_ST, param.n_ex, param.k_ex, epsilon,
                                  Gamma_local, beta_local, Delta, RS_index, param.open_angle);

            double R_start = 1e10; 
            y[0] = R_start;
            y[1] = 0.0;          
            y[2] = Gamma_local;  
            y[3] = param.open_angle;
            y[4] = 0.0;          
            y[5] = M_ej_local / 1e6; 
            y[6] = 0.0;          
            y[7] = Gamma_local;  
            y[8] = R_start / (c_cnst * beta_local * Gamma_local); 
            y[9] = y[8];         
            y[10] = R_start;     

            double x_start = 0; 
            
            dyn(x_start, y, dydx); 
            
            rk5<DynamicThickShell> evol(dyn, y, dydx, yout, x_start, n_vars);

            // Integrate to Target Time
            shock_dynamic = evol.integrateEATS(T / (1.0 + param.z), jet.mu[i][0]);

            // 5. Store Results
            // Loop over phi: Assuming dynamic variables are phi-independent for a given theta ring
            for (size_t j = 0; j < jet.phi.size(); j++) {
                jet.radius[i][j]              = shock_dynamic[0];
                jet.Gamma[i][j]               = shock_dynamic[2];
                jet.mass_ejecta[i][j]         = shock_dynamic[5];
                jet.dynamical_timescale[i][j] = shock_dynamic[8];
                jet.adiabatic_timescale[i][j] = shock_dynamic[9];
                jet.Gamma_RS[i][j]            = shock_dynamic[7];
                jet.radius_RS[i][j]           = shock_dynamic[10];

                jet.radius_cross[i][j]              = dyn.radius_cross;
                jet.Gamma_cross[i][j]               = dyn.Gamma_cross;
                jet.mass_cross[i][j]                = dyn.mass_cross;
                jet.dynamical_timescale_cross[i][j] = dyn.dynamical_timescale_cross;
                jet.adiabatic_timescale_cross[i][j] = dyn.adiabatic_timescale_cross;
            }
        }
    }
    // --- OPTIMIZATION END ---
}

// Constructor
DynamicThickShell::DynamicThickShell(double M_ej, double C_BM, double C_ST,
                                     double n_ex, double k_ex, double epsilon,
                                     double Gamma0, double beta0, double Delta,
                                     double RS_index, double theta_c)
    : M_ej_total(M_ej), C_BM(C_BM), C_ST(C_ST), n_ex(n_ex), k_ex(k_ex),
      epsilon(epsilon), Gamma0(Gamma0), beta0(beta0), Delta(Delta),
      RS_index(RS_index), theta_c(theta_c) {}

void DynamicThickShell::operator()(const double x, std::vector<double>& y,
                                   std::vector<double>& dydx) {
    
    // --- Unpack State Variables ---
    double R     = y[0];
    double m_sw  = y[1];
    double Gamma = y[2];
    double theta = y[3];
    double E_int = y[4];
    double M_inj = y[5];
    
    // --- 1. Geometry & Kinematics ---
    double beta = sqrt(1.0 - 1.0 / (Gamma * Gamma));
    
    // Blandford-McKee coefficients
    double gamma_hat = (4.0 * Gamma + 1.0) / (3.0 * Gamma); 
    double Gamma_eff = (gamma_hat * Gamma * Gamma - gamma_hat + 1.0) / Gamma;
    
    // Density Profile & Swept Mass Rate
    double rho_ex = mp * n_ex * pow(R, -k_ex);
    // A = dM_sw / dR
    double A = 2.0 * PI * (1.0 - cos(theta)) * R * R * rho_ex * (beta * (C_BM - C_ST) + C_ST);

    // --- 2. Base Derivatives ---
    dydx[0] = beta * c_cnst;    // dR/dt
    dydx[1] = A * dydx[0];      // dM_sw/dt

    // --- 3. Energy Injection (Thick Shell Logic) ---
    bool is_injecting = (M_inj < M_ej_total);

    if (is_injecting) {
        double injection_rate = M_ej_total / Delta * (beta0 - beta) / beta;
        if (injection_rate < 0) injection_rate = 0; 
        
        dydx[5] = injection_rate * dydx[0]; // dM_inj/dt
        dydx[6] = 0.0; // Assuming cold injection (Kinetic dominates)
    } else {
        dydx[5] = 0.0;
        dydx[6] = 0.0;
    }

    // --- 4. SOLVING THE CIRCULAR DEPENDENCY ---
    // We need dGamma/dt.
    // The energy equation is: denom * dGamma = Numerator_Terms - Gamma_eff * dE_ad
    // But dE_ad depends on dGamma.
    // Let's group terms: dE_ad = Term_Volume - Term_Gamma * dGamma
    
    double K_ad = (gamma_hat - 1.0) * E_int; // Adiabatic pre-factor
    double term_vol = -K_ad * (3.0 / R) * dydx[0]; // Volume expansion part
    double term_gam = (K_ad / Gamma);             // Feedback part related to dGamma

    // Terms for the Gamma Equation
    double denom_inertia = (M_inj + m_sw) * c_cnst * c_cnst + E_int * ((gamma_hat * Gamma * Gamma + gamma_hat - 1.0) / (Gamma * Gamma)); // Includes dGamma_eff approx
    
    double num_sw  = (Gamma_eff + 1.0) * (Gamma - 1.0) * c_cnst * c_cnst * dydx[1];
    double num_inj = (Gamma - Gamma0) * c_cnst * c_cnst * dydx[5];
    
    // Combine to solve for dGamma/dt
    // Original: D * dG = -(N_sw + N_inj) - G_eff * (Term_Vol + Term_Gam * dG)
    // Rearranged: dG * (D + G_eff * Term_Gam) = -(N_sw + N_inj) - G_eff * Term_Vol
    
    double lhs_factor = denom_inertia + Gamma_eff * term_gam;
    double rhs_force  = -(num_sw + num_inj) - Gamma_eff * term_vol;

    dydx[2] = rhs_force / lhs_factor; // This is the correct, stable dGamma/dt

    // --- 5. Internal Energy Evolution ---
    // Now we can calculate dE_ad correctly using the known dydx[2]
    // dE_ad/dt term
    double dE_ad_dt = term_vol + term_gam * dydx[2];
    
    dydx[4] = (1.0 - epsilon) * (Gamma - 1.0) * dydx[1] * c_cnst * c_cnst + dE_ad_dt;

    // --- 6. Jet Spreading ---
    if (beta * Gamma < 1.0 / theta_c) {
         double cs = c_cnst / sqrt(3.0); 
         dydx[3] = (cs / R) * dydx[0]; 
    } else {
         dydx[3] = 0.0;
    }
    // Override if you want to test spherical first
    dydx[3] = 0.0;

    // --- 7. Reverse Shock (RS) Logic ---
    if (is_injecting) {
        // Phase 1: Injection (RS is crossing shell)
        // Gamma_RS is roughly Gamma_FS (relative to lab) in simple models
        dydx[7] = dydx[2]; 
        
        // Note: For RK solver, y[7] will track Gamma via integration of dydx[7]
        // If you want to force y[7] = Gamma algebraically, you strictly shouldn't use 
        // a differential variable for it. But setting dydx[7]=dydx[2] keeps them coupled.
    } else {
        // Phase 2: Post-Crossing (RS Decays)
        
        // Latch crossing values
        if (!cross_occurred) {
            cross_occurred = true;
            Gamma_cross = Gamma;
            beta_cross = beta;
            radius_cross = R;
            mass_cross = M_inj;
            dynamical_timescale_cross = y[8];
            adiabatic_timescale_cross = y[9];
        }

        // Calculate Derivative of Power Law Decay
        // Gamma_RS(R) = Gamma_cross * (R/R_cross)^(-RS_index)
        // dGamma_RS/dt = Gamma_cross * R_cross^k * (-k) * R^(-k-1) * dR/dt
        //              = Gamma_RS_current * (-k / R) * dR/dt
        
        // Ideally use y[7] as current Gamma_RS
        double current_Gamma_RS = y[7];
        dydx[7] = current_Gamma_RS * (-RS_index / R) * dydx[0];
        
        dydx[5] = 0.0; 
        dydx[6] = 0.0; 
    }

    // --- 8. Timescales ---
    // Now safe to use dydx[2]
    dydx[8] = (1.0 / (c_cnst * beta * Gamma)) * dydx[0]; // t_dyn
    
    // t_adi approximation using calculated derivatives
    //double expansion_rate = (1.0/R) * dydx[0] - (1.0/(3.0*Gamma)) * dydx[2];
    double expansion_rate = (1.0/R) - (1.0/(3.0*Gamma)) * dydx[2] / dydx[0];
    if (std::abs(expansion_rate) > 1e-20) {
        // t_adi = 1 / (c * beta * Gamma * expansion_rate) -- simplified
        // We integrate a placeholder or algebraic
        dydx[9] = 0; 
        //y[9] = 1.0 / (c_cnst * beta * expansion_rate * 3.0); // Rough scaling
    } else {
        dydx[9] = 0;
        y[9] = y[8];
    }
    y[9] = 1.0 / (c_cnst * beta * Gamma * expansion_rate); // Rough scaling

    // --- 9. RS Radius ---
    // If RS is ultra-relativistic, beta_RS ~ 1
    double beta_RS = (y[7] > 1.001) ? sqrt(1.0 - 1.0/(y[7]*y[7])) : 0.0;
    dydx[10] = beta_RS * c_cnst; 
}

void GRBafterglow::calcDynamic(const std::vector<double>& time_array) {
    
    // 1. Prepare Data Container
    // Logic: [0]=R, [1]=m_sw, [2]=Gamma, [7]=Gamma_RS, [10]=R_RS (Based on DynamicThickShell)
    std::vector<std::vector<double>> shock_dynamic;
    
    // 2. Open File
    std::string filename = output_folder + "/dynamic_thick.dat"; // Renamed for clarity
    std::ofstream mystream(filename);
    
    if (!mystream.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return;
    }

    // 3. Write Header
    // Added Reverse Shock (RS) columns which are critical for Thick Shell analysis
    mystream << "# Time(s)  Gamma_FS  Radius_FS(cm)  Gamma_RS  Radius_RS(cm)  GammaBeta_FS  GammaShBetaSh_FS\n";
    mystream << std::scientific << std::setprecision(16);

    // 4. Calculate Dynamics (THICK SHELL)
    // Note: Argument order is (input_time, output_vector) based on your previous function definition
    EvolutionThickShell(time_array, shock_dynamic); 

    // 5. Helper Lambda for Forward Shock Physics (Blandford-McKee)
    auto calculate_gamma_shock = [](double Gamma) {
        if (Gamma <= 1.00000001) return 1.0;

        double gamma_hat = (4.0 * Gamma + 1.0) / (3.0 * Gamma);
        
        double num = (Gamma + 1.0) * pow(gamma_hat * (Gamma - 1.0) + 1.0, 2);
        double den = gamma_hat * (2.0 - gamma_hat) * (Gamma - 1.0) + 2.0;
        
        return sqrt(num / den);
    };

    // 6. Output Loop
    if (shock_dynamic.size() != time_array.size()) {
        std::cerr << "Error: Simulation steps (" << shock_dynamic.size() 
                  << ") mismatch input time array (" << time_array.size() << ")." << std::endl;
        return;
    }

    for (size_t i = 0; i < time_array.size(); i++) {
        
        const auto& state = shock_dynamic[i]; 

        // --- MAPPING INDICES (From DynamicThickShell::operator()) ---
        // y[0] = Radius
        // y[2] = Gamma (Forward Shock / Bulk)
        // y[7] = Gamma (Reverse Shock)
        // y[10] = Radius (Reverse Shock)
        
        double t         = time_array[i];
        double Radius_FS = state[0];
        double Gamma_FS  = state[2];
        double Gamma_RS  = state[7];
        double Radius_RS = state[10];

        // Safety clamp for output
        if (Gamma_FS < 1.0) Gamma_FS = 1.0;
        if (Gamma_RS < 1.0) Gamma_RS = 1.0;

        // Kinematics (Forward Shock 4-velocity)
        double u_FS = sqrt(Gamma_FS * Gamma_FS - 1.0); 

        // Shock Physics (Thermal efficiency of the shock front)
        double gamma_sh_FS = calculate_gamma_shock(Gamma_FS);
        double u_sh_FS = sqrt(gamma_sh_FS * gamma_sh_FS - 1.0);

        // Write to file
        mystream << t << " " 
                 << Gamma_FS << " " 
                 << Radius_FS << " " 
                 << Gamma_RS << " " 
                 << Radius_RS << " " 
                 << u_FS << " " 
                 << u_sh_FS << "\n";
    }

    std::cout << "Thick shell dynamics written to " << filename << std::endl;
}
