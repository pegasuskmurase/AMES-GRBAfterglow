//========================================================================================
// AMES
// Copyright(C) 2024, licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file GRBAfterglow.hpp
// \brief
// GRBafterglow.cpp is optimized for "on-axis" afterglow emission from an
// external forward shock of gamma-ray bursts (GRBs) Citation: B.T. Zhang, K.
// Murase, P. Veres, & P. Meszaros, Astrophys.J. 920 (2021) 55

#ifndef GRBAFTERGLOW_H
#define GRBAFTERGLOW_H

#include "module/compton.hpp"
#include "module/constants.hpp"
#include "module/gammagamma.hpp"
#include "module/photonbackground.hpp"
#include "module/source.hpp"
#include "module/synchrotron.hpp"
#include "module/electrondistribution.hpp"

#pragma once
#include <iostream>
#include <algorithm>
#include <string>


struct Jet {
    int angular_num = 50;
    int phi_num = 50;
    std::vector<double> theta;
    std::vector<double> theta_bin;
    std::vector<double> phi;
    std::vector<double> phi_bin;
    std::vector<double> Gamma0;
    std::vector<double> M_ej;
    std::vector<double> E_ej;
    std::vector<std::vector<double>> Gamma;
    std::vector<std::vector<double>> radius;
    std::vector<std::vector<double>> mu;
    std::vector<std::vector<double>> dynamical_timescale;
    std::vector<std::vector<double>> adiabatic_timescale;
    std::vector<std::vector<double>> Gamma_RS;
    std::vector<std::vector<double>> radius_RS;
    std::vector<std::vector<double>> mass_ejecta;
    std::vector<std::vector<double>> mass_cross;
    std::vector<std::vector<double>> Gamma_cross;
    std::vector<std::vector<double>> radius_cross;
    std::vector<std::vector<double>> dynamical_timescale_cross;
    std::vector<std::vector<double>> adiabatic_timescale_cross;
};

struct GRBafterglowParam {
    // --- Basic GRB ---
    double z = 1.0;
    double dl = 3.1e28;
    double E_ej = 1e53;
    double Gamma0 = 300.0;

    // --- Medium ---
    double n_ex = 1.0;
    double A_star = 1e-1;
    double k_ex = 0.0;

    // --- FS Electrons ---
    double spectral_e = 2.2;
    double epsilon_e = 0.1;
    double fraction_e = 1.0;
    double eta_acc_e = 10.0;

    // --- FS Magnetic ---
    double epsilon_B = 1e-3;

    // --- Geometry ---
    double open_angle = 0.1;
    double view_angle = 0.05;
    double gaussian_cone = 0.1;
    double jet_index = 1.0;

    // --- Duration ---
    double T_ej = 10.0;

    // --- RS Electrons ---
    double spectral_e_RS = 2.2;
    double epsilon_e_RS = 0.1;
    double eta_acc_e_RS = 10.0;
    double epsilon_B_RS = 1e-3;
    double fraction_e_RS = 1.0;

    // --- Prompt ---
    double L_ph = 1e51;
};

// Helper struct for validation limits
struct Range {
    double min, max, def;
};

class GRBafterglow {
  public:
    GRBafterglow(Source& s);

    GRBafterglowParam param;

    void help();
    void setGRBafterglowParam(const std::vector<double>& _param);
    void setGRBafterglowRSParam(const std::vector<double>& _param);

    // External photon density
    void setExtendedEmissionParam(std::vector<double> _param);
    void setPlateauEmissionParam(std::vector<double> _param);

    void setOutputFolder(std::string _output_folder);

    // set output
    void showInfo(bool _show_info);       // showinfo
    void haveReverseShock(bool _have_RS); // have RS
    void haveOneZone(bool _have_onezone); // One zone
    void haveSSCSpec(bool _have_SSCSpec); // SSC
    void haveEICSpec(bool _have_EICSpec); // EIC
    void haveEdgeEffect(
        bool
            _have_edge_effect); // jet break, gemetrical edge effect, only
                                // applicable for one-zone case; The integration
                                // over EATS automatically involved jet break
    void haveAttenuSSA(bool _have_attenu_SSA); // SSA
    void haveAttenuFF(bool _have_attenu_FF);   // Free-free
    void haveAttenuGGSource(
        bool _have_attenu_GG_source); // gamma-gamma attenuation inside sources
    void haveAttenuGGCosmic(
        bool _have_attenu_GG_cosmic); // gamma-gamma attenuation in the
                                      // intergalactic space

    void FluxOutput(const std::vector<double>& time_array,
              const std::vector<double>& energy_array_min,
              const std::vector<double>& energy_array_max);

    void Flux1D(std::vector<double>& flux_vector,
                const std::vector<double>& time_array,
                const std::vector<double>& energy_array_min,
                const std::vector<double>& energy_array_max);

    void Flux(std::vector<std::vector<double>>& flux_vector,
              const std::vector<double>& time_array,
              const std::vector<double>& energy_array_min,
              const std::vector<double>& energy_array_max);

    // New MCMC specific functions
    void setParamsFromMCMC(const std::vector<double>& theta);
    
    // Optimized Flux calculator
    void CalculateLightCurveMCMC(std::vector<double>& flux_out,
                                 const std::vector<double>& time_array,
                                 const std::vector<double>& freq_array);

    void Spectrum(std::vector<std::vector<double>>& flux_vector, double T);

    void Spectrum(double T);
    void SpectrumOneZone(double T);
    void ReverseShock(Jet& jet);

    void InitJet(Jet& jet);

    void EvolutionThickShell(const std::vector<double>& time_array,
                             std::vector<std::vector<double>>& jet_dyn);
    void EvolutionThickShell(Jet& jet, double T);
    void calcDynamic(const std::vector<double>& time_array);

    void PromptPhotonDensity(Photonbackground& ph, double T, double radius,
                             double Gamma);
    double ExtendedEmissionLuminosity(double T);
    double PlateauEmissionLuminosity(double T);

    double ExternalDensity(double radius);
    double MomentumInj(double Gamma, double mag_strength, double radius, double beta, double epsilon_e, double fraction_e, double spectral_e, double eta_acc_e);
    double MomentumMax(double mag_strength, double radius, double Gamma, double eta_acc_e);
    double DecelerationRadius();
    double DecelerationRadius(double E_ej, double Gamma);
    void EvolutionThinShellAnalytic(double T, double& Gamma, double& radius);
    void LorentzBoost(const std::vector<double>& photon_energy,
                      std::vector<double>& photon_spectrum, double Gamma);
    void FluxNorm(const std::vector<double>& photon_energy,
                  std::vector<double>& photon_spectrum, double radius,
                  double Delta);

    double MomentumInjProton(double Gamma);

    double MomentumMaxDyn(double mag_strength, double radius,
                          double Gamma, double beta, double Z,
                          double eta_acc);
  private:
    Source& s;
    Utility utility;
    Jet jet;

    Synchrotron& syn;
    InverseCompton& IC;
    GammaGamma& gg;
    Photonbackground& ph;
    ElectronDistribution& ElectronD;

    std::string output_folder = ".";

    bool have_onezone = false;

    // absorption factor
    bool have_SSCSpec = false;
    bool have_EICSpec = false;
    bool have_edge_effect = false;
    bool have_attenu_SSA = false;
    bool have_attenu_FF = false;
    bool have_attenu_GG_source = false;
    bool have_attenu_GG_cosmic = false;
    bool have_RS = false;
    bool show_info = false;

    // Store previous geometric params to check if InitJet is needed
    double prev_view_angle = -1.0;
    double prev_open_angle = -1.0;
    double prev_Gamma0 = -1.0;
    int prev_angular_num = -1;

    std::vector<double> target_energy;
    std::vector<double> target_syn;
    std::vector<double> target_prompt;
    std::vector<double> target_tot;
    std::vector<double> photon_energy;
    std::vector<double> photon_syn;
    std::vector<double> photon_ssc;
    std::vector<double> photon_eic;
    std::vector<double> photon_syn_temp;
    std::vector<double> photon_syn_absorption_SSA;
    std::vector<double> photon_syn_obs;
    std::vector<double> photon_ssc_temp;
    std::vector<double> photon_ssc_obs;
    std::vector<double> photon_eic_temp;
    std::vector<double> photon_eic_obs;
    std::vector<double> photon_losstime;
    std::vector<double> photon_tot;
    std::vector<double> photon_tot_diff;
    std::vector<double> photon_syn_diff;
    std::vector<double> photon_ssc_diff;
    std::vector<double> photon_eic_diff;
    std::vector<double> photon_syn_RS;
    std::vector<double> photon_ssc_RS;
    std::vector<double> photon_eic_RS;
    std::vector<double> photon_syn_RS_obs;
    std::vector<double> photon_ssc_RS_obs;
    std::vector<double> photon_eic_RS_obs;

    std::vector<double> electron_energy;
    std::vector<double> proton_energy;
    std::vector<double> timescale_dyn_out;
    std::vector<double> timescale_syn_out;
    std::vector<double> timescale_pa_out;
    std::vector<double> timescale_pm_out;
    std::vector<double> timescale_dyn_RS_out;
    std::vector<double> timescale_syn_RS_out;
    std::vector<double> timescale_pa_RS_out;
    std::vector<double> timescale_pm_RS_out;

    double syn_dum;
    double ssc_dum;
    int num_t;
    int num_p;
    int num_e;
    int num_nu;
    int num_proton;

    double beta, beta_sh, doppler_factor, delta_s, EATS_factor, delta_sh;
    double mag_strength, photon_escape_timescale, electron_inj;
    double alpha, tau, att;

    double L_EE, T_a_EE, alpha_fl_EE, alpha_st_EE;
    double energy_min_EE, energy_peak_EE, energy_max_EE, spectra_l_EE,
        spectra_h_EE;
    double L_PL, T_a_PL, alpha_fl_PL, alpha_st_PL;
    double energy_min_PL, energy_peak_PL, energy_max_PL, spectra_l_PL,
        spectra_h_PL;
};
#endif

#ifndef DYNAMICTHICKSHELL_H
#define DYNAMICTHICKSHELL_H
// Dynamics with Reverse Shock (RS) and Forward Shock (FS) energy injection
class DynamicThickShell {
  public:
    DynamicThickShell(double M_ej, double C_BM, double C_ST,
                      double n_ex, double k_ex,
                      double epsilon, double Gamma0,
                      double beta0, double Delta, double RS_index,
                      double theta_c);

    // The ODE system operator
    void operator()(const double x, std::vector<double>& y, std::vector<double>& dydx);

    // Public members to store the state at the moment of crossing (Shell crossing radius)
    // These are read by the main loop after integration or used during the post-crossing phase.
    bool cross_occurred = false;
    double Gamma_cross = 0.0;
    double beta_cross = 0.0;
    double radius_cross = 0.0;
    double mass_cross = 0.0;
    double dynamical_timescale_cross = 0.0;
    double adiabatic_timescale_cross = 0.0;

  private:
    // Model Constants
    const double M_ej_total;  // Total mass to be ejected
    const double C_BM;        // Blandford-McKee constant
    const double C_ST;        // Sedov-Taylor constant
    const double n_ex;        // External density normalization
    const double k_ex;        // Density profile index
    const double epsilon;     // Radiative efficiency
    const double Gamma0;      // Initial Lorentz factor of ejecta
    const double beta0;       // Initial velocity
    const double Delta;       // Shell thickness parameter
    const double RS_index;    // Power law index for RS decay
    const double theta_c;     // Jet opening angle

    // Physics Constants
    const double mp = 1.6726219e-24; // Proton mass [g]
    const double c_cnst = 2.99792458e10; // Speed of light [cm/s]
    const double PI = 3.14159265358979323846;
};

#endif
