//========================================================================================
// AMES code
// Copyright(C) 2023, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file GRBAfterglow.hpp
// \brief Calculate GRB afterglow emission

#ifndef GRBAFTERGLOW_H
#define GRBAFTERGLOW_H

#include "module/constants.hpp"
#include "module/electrondistribution.hpp"
#include "module/photonbackground.hpp"
#include "module/gammagamma.hpp"
#include "module/inversecompton.hpp"
#include "module/source.hpp"
#include "module/synchrotron.hpp"

struct Jet {
  int angular_num = 30;
  int phi_num = 20;
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
};

class GRBAfterglow {
public:
  GRBAfterglow(Source &s);

  void help();
  void Init();
  void setGRBAfterglowParam(std::vector<double> _param);

  void setOutputFolder(std::string _output_folder);

  // set output
  void haveOneZone(bool _have_onezeon);                 // One zone
  void haveSSCSpec(bool _have_SSCSpec);                 // SSC
  void haveEdgeEffect(bool _have_edge_effect);          // jet break, edge effect
  void haveAttenuSSA(bool _have_attenu_SSA);            // SSA
  void haveAttenuFF(bool _have_attenu_FF);              // Free-free
  void haveAttenuGGSource(bool _have_attenu_GG_source); // gamma-gamma attenuation inside sources
  void haveAttenuGGCosmic(
      bool _have_attenu_GG_cosmic); // gamma-gamma attenuation in the intergalactic space

  void Flux(ElectronDistribution &ED, Synchrotron &syn, InverseCompton &IC, GammaGamma &gg, Photonbackground &ph, const std::vector<double> &time_array,
            const std::vector<double> &energy_array_min,
            const std::vector<double> &energy_array_max);

  void Flux(std::vector<std::vector<double>> &flux_vector, ElectronDistribution &ED,
                  Synchrotron &syn, InverseCompton &IC, GammaGamma &gg, Photonbackground &ph,
                  const std::vector<double> &time_array,
                  const std::vector<double> &energy_array_min,
                  const std::vector<double> &energy_array_max);


  void Spectrum(ElectronDistribution &ED, Synchrotron &syn, InverseCompton &IC, GammaGamma &gg,
                double T);

  void InitJet(Jet &jet);

  void EvolutionThinShell(Jet &jet, double T);

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

private:
  Source &s;
  GRBAfterglowParam param;
  Utility utility;
  Jet jet;

  std::string output_folder = ".";

  bool have_onezone = false;

  // absorption factor
  bool have_SSCSpec = false;
  bool have_edge_effect = false;
  bool have_attenu_SSA = false;
  bool have_attenu_FF = false;
  bool have_attenu_GG_source = false;
  bool have_attenu_GG_cosmic = false;

  std::vector<double> target_energy;
  std::vector<double> photon_energy;
  std::vector<double> photon_syn;
  std::vector<double> photon_ssc;
  std::vector<double> photon_syn_temp;
  std::vector<double> photon_syn_absorption_SSA;
  std::vector<double> photon_syn_obs;
  std::vector<double> photon_ssc_temp;
  std::vector<double> photon_ssc_obs;
  std::vector<double> photon_losstime;
  std::vector<double> target_syn;
  std::vector<double> photon_tot;
  std::vector<double> photon_tot_diff;
  double syn_dum;
  double ssc_dum;
  int num_t;
  int num_p;
  int num_e;
};
#endif

// dynamics with FS
class DynamicThinShell {
public:
  DynamicThinShell(const double M_ej, const double n_ex, const double k_ex, const double epsilon);
  void operator()(const double x, std::vector<double> &y, std::vector<double> &dydx);

private:
  const double M_ej;
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
