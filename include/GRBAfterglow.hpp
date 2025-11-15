//========================================================================================
// AMES
// Copyright(C) 2024, licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file GRBAfterglow.hpp
// \brief 
// GRBafterglow.cpp is optimized for "on-axis" afterglow emission from an external forward shock of gamma-ray bursts (GRBs) 
// Citation:
// B.T. Zhang, K. Murase, P. Veres, & P. Meszaros, Astrophys.J. 920 (2021) 55


#ifndef GRBAFTERGLOW_H
#define GRBAFTERGLOW_H

#include "module/constants.hpp"
#include "module/electrondistribution.hpp"
#include "module/gammagamma.hpp"
#include "module/compton.hpp"
#include "module/photonbackground.hpp"
#include "module/source.hpp"
#include "module/synchrotron.hpp"

struct Jet {
  int angular_num = 30;
  int phi_num = 2;
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


struct GRBAfterglowParam {
  double z;          // redshift
  double dl;         // luminosity distance [cm]
  double E_ej;       // ejecta kinetic energy [erg]
  double Gamma0;     // ejecta initial Lorentz factor
  double n_ex;       // external medium density [cm^-3]
  double k_ex;       // external medium density spectral index
  double spectral_e; // non-thermal electron spectral index
  double epsilon_e;  // the energy fraction of internal energy that goes into
                     // electrons
  double fraction_e; // the number fraction of electrons that are accelerated
  double eta_acc_e;  // electron acceleration efficiency
  double epsilon_B;  // the energy fraction of internal energy that is converted
                     // into the magnetic energy
  double open_angle; // the openning angle of the outflow [radian]
  double view_angle; // the view angle [radian]
  double gaussian_cone;
  double jet_index;
  double T_ej; //
  double spectral_e_RS; //
  double epsilon_e_RS;  //
  double eta_acc_e_RS;  //
  double epsilon_B_RS;  //
  double fraction_e_RS; //
  double L_ph;
};

class GRBAfterglow {
public:
  GRBAfterglow(Source &s);

  void help();
  void setGRBAfterglowParam(std::vector<double> _param);
  void setGRBAfterglowRSParam(std::vector<double> _param);

  // External photon density
  void setExtendedEmissionParam(std::vector<double> _param);
  void setPlateauEmissionParam(std::vector<double> _param);

  void setOutputFolder(std::string _output_folder);

  // set output
  void haveReverseShock(bool _have_RS);                 // have RS
  void haveOneZone(bool _have_onezeon); // One zone
  void haveSSCSpec(bool _have_SSCSpec); // SSC
  void haveEICSpec(bool _have_EICSpec); // EIC
  void haveEdgeEffect(
      bool _have_edge_effect); // jet break, gemetrical edge effect, only applicable for one-zone
                               // case; The integration over EATS automatically involved jet break
  void haveAttenuSSA(bool _have_attenu_SSA);            // SSA
  void haveAttenuFF(bool _have_attenu_FF);              // Free-free
  void haveAttenuGGSource(bool _have_attenu_GG_source); // gamma-gamma attenuation inside sources
  void haveAttenuGGCosmic(
      bool _have_attenu_GG_cosmic); // gamma-gamma attenuation in the intergalactic space

  void Flux(const std::vector<double> &time_array,
            const std::vector<double> &energy_array_min,
            const std::vector<double> &energy_array_max);

  void Flux(std::vector<std::vector<double>> &flux_vector, 
                  const std::vector<double> &time_array,
                  const std::vector<double> &energy_array_min,
                  const std::vector<double> &energy_array_max);

  void Spectrum(std::vector<std::vector<double>>& flux_vector, double T);

  void Spectrum(double T);
  void ReverseShock(Jet& jet);

  void InitJet(Jet &jet);

  void EvolutionThinShell(Jet &jet, double T);

  void EvolutionThinShell(std::vector<std::vector<double>> &shock_dynamic, const std::vector<double> &time_array);
  void EvolutionThickShell(const std::vector<double> &time_array,
                           std::vector<std::vector<double>> &jet_dyn);
  void EvolutionThickShell(Jet &jet, double T);
  void calcDynamic(const std::vector<double> &time_array);

  void PromptPhotonDensity(Photonbackground &ph, double T, double radius, double Gamma);
  double ExtendedEmissionLuminosity(double T);
  double PlateauEmissionLuminosity(double T);

  double ExternalDensity(double radius);
  double MomentumInj(double Gamma);
  double MomentumMax(double mag_strength, double radius, double beta);
  double DecelerationRadius();
  double DecelerationRadius(double E_ej, double Gamma);
  void EvolutionThinShellAnalytic(double T, double &Gamma, double &radius);
  void LorentzBoost(const std::vector<double> &photon_energy, std::vector<double> &photon_spectrum,
                    double Gamma);
  void FluxNorm(const std::vector<double> &photon_energy, std::vector<double> &photon_spectrum,
                double radius, double Delta);

  double MomentumInjRS(double Gamma);
  double MomentumMaxRS(double mag_strength, double radius, double beta);

private:
  Source &s;
  GRBAfterglowParam param;
  Utility utility;
  Jet jet;

  Synchrotron &syn;
  InverseCompton &IC;
  GammaGamma &gg;
  Photonbackground &ph;
  ElectronDistribution &ED;

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
  std::vector<double> photon_psyn_RS;
  std::vector<double> photon_syn_RS_obs;
  std::vector<double> photon_ssc_RS_obs;
  std::vector<double> photon_eic_RS_obs;
  std::vector<double> photon_psyn_RS_obs;
  double syn_dum;
  double ssc_dum;
  int num_t;
  int num_p;
  int num_e;

  double beta, beta_sh, doppler_factor, delta_s, EATS_factor, delta_sh;
  double mag_strength, photon_escape_timescale, electron_inj;
  double alpha, tau, att;

  double L_EE, T_a_EE, alpha_fl_EE, alpha_st_EE;
  double energy_min_EE, energy_peak_EE, energy_max_EE, spectra_l_EE, spectra_h_EE;
  double L_PL, T_a_PL, alpha_fl_PL, alpha_st_PL;
  double energy_min_PL, energy_peak_PL, energy_max_PL, spectra_l_PL, spectra_h_PL;
};
#endif

#ifndef DYNAMICTHINSHELL_H
#define DYNAMICTHINSHELL_H
// dynamics with FS
class DynamicThinShell {
public:
  DynamicThinShell(const double M_ej, const double C_BM, const double C_ST, const double n_ex, const double k_ex, const double epsilon);
  void operator()(const double x, std::vector<double> &y, std::vector<double> &dydx);

private:
  const double M_ej;
  const double C_BM;
  const double C_ST;
  const double n_ex;
  const double k_ex;
  const double epsilon;
  double beta;
  double gamma_hat;
  double Gamma_eff;
  double dGamma_eff;
  double dEad;
  double A;
  double GammaBeta;
};
#endif

#ifndef DYNAMICTHICKSHELL_H
#define DYNAMICTHICKSHELL_H
// dynamics with RS and FS
class DynamicThickShell {
public:
  DynamicThickShell(const double M_ej, const double C_BM, const double C_ST, const double n_ex, const double k_ex,
                    const double epsilon, const double Gamma0, const double beta0, double Delta,
                    const double RS_index, const double theta_c);
  void operator()(const double x, std::vector<double> &y, std::vector<double> &dydx);

  double Gamma_cross = 0.0;
  double beta_cross = 0.0;
  double radius_cross = 0.0;
  double mass_cross = 0.0;
  double dynamical_timescale_cross = 0.0;
  double adiabatic_timescale_cross = 0.0;

private:
  const double M_ej;
  const double C_BM;
  const double C_ST;
  const double n_ex;
  const double k_ex;
  const double epsilon;
  const double Gamma0;
  const double beta0;
  double RS_index;
  double theta_c;
  double beta;
  double gamma_hat;
  double Gamma_eff;
  double dGamma_eff;
  double Gamma34;
  double dEad2, dEad3;
  double Eint3;
  double M_ej3;
  double Delta;
  double A, A_RS;
  double GammaBeta;
  bool have_cross = true;
};
#endif
