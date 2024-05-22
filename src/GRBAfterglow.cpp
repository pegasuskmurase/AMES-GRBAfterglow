// AMES code
// Copyright(C) 2022, Bing Theodore Zhang
// <bing.zhang@yukawa.kyoto-u.ac.jp> and Kohta Murase
// <murase@psu.edu> licensed under the GNU GENERAL PUBLIC
// LICENSE, see LICENSE file for details

// \file GRBAfterglow.cpp
// \brief Calculate emissions GRB afterglow

#include "GRBAfterglow.hpp"
#include "module/rk5.hpp"

#include <fstream>
#include <iostream>

GRBAfterglow::GRBAfterglow(Source &s) : s(s) {
  s.getTarget().setMomentum(321, 1e-9, 1e10);
  s.getPhoton().setMomentum(241, 1e-8, 1e16);
  s.getElectron().setMomentum(241, 1e4, 1e16);
}

void GRBAfterglow::help() {
  std::cout << "----Welcome to use AMES-GRBAfterglow------" << std::endl;
  std::cout << "Please see the following information for "
               "proper use."
            << std::endl;
  std::cout << "If you have any questions please contact "
               "Bing Theodore "
               "Zhang<bing.zhang @yukawa.kyoto-u.ac.jp> "
               "and Kohta "
               "Murase<murase@psu.edu>."
            << std::endl;
  std::cout << "***The simpliest way to use the code***" << std::endl;
}

void GRBAfterglow::Init() {
  num_t = s.getTarget().getEnergy().size();
  num_p = s.getPhoton().getEnergy().size();
  num_e = s.getElectron().getEnergy().size();
  target_energy = s.getTarget().getEnergy();
  photon_energy = s.getPhoton().getEnergy();
  photon_syn.resize(num_p, 0.0);
  photon_ssc.resize(num_p, 0.0);
  photon_syn_temp.resize(num_p, 0.0);
  photon_syn_absorption_SSA.resize(num_p, 0.0);
  photon_syn_obs.resize(num_p, 0.0);
  photon_ssc_temp.resize(num_p, 0.0);
  photon_ssc_obs.resize(num_p, 0.0);
  photon_losstime.resize(num_p, 0.0);
  target_syn.resize(num_t, 0.0);
  photon_tot.resize(num_p, 0.0);
  photon_tot_diff.resize(num_p, 0.0);
}

void GRBAfterglow::setGRBAfterglowParam(std::vector<double> _param) {
  if ((_param[0] >= 0) && (_param[0] <= 10)) {
    param.z = _param[0];
  } else {
    param.z = 0.0;
    std::cout << "Warning: the redshift parameter [z] should within the range (1, 10) "
              << std::endl;
  }
  param.dl = _param[1];
  param.E_ej = _param[2];
  param.Gamma0 = _param[3];
  param.n_ex = _param[4];
  param.k_ex = _param[5];
  param.spectral_e = _param[6];
  param.epsilon_e = _param[7];
  param.fraction_e = _param[8];
  param.eta_acc_e = _param[9];
  param.epsilon_B = _param[10];
  param.open_angle = _param[11];
  param.view_angle = _param[12];
  param.gaussian_cone = _param[13];
  param.jet_index = _param[14];
}

void GRBAfterglow::setOutputFolder(std::string _output_folder) { output_folder = _output_folder; }

void GRBAfterglow::haveOneZone(bool _have_onezone) { have_onezone = _have_onezone; } // One zone
void GRBAfterglow::haveSSCSpec(bool _have_SSCSpec) { have_SSCSpec = _have_SSCSpec; }
void GRBAfterglow::haveEdgeEffect(bool _have_edge_effect) { have_edge_effect = _have_edge_effect; }
void GRBAfterglow::haveAttenuSSA(bool _have_attenu_SSA) { have_attenu_SSA = _have_attenu_SSA; }
void GRBAfterglow::haveAttenuFF(bool _have_attenu_FF) { have_attenu_FF = _have_attenu_FF; }
void GRBAfterglow::haveAttenuGGSource(bool _have_attenu_GG_source) {
  have_attenu_GG_source = _have_attenu_GG_source;
}
void GRBAfterglow::haveAttenuGGCosmic(bool _have_attenu_GG_cosmic) {
  have_attenu_GG_cosmic = _have_attenu_GG_cosmic;
}

void GRBAfterglow::Flux(std::vector<std::vector<double>> &flux_vector, ElectronDistribution &ED,
                        Synchrotron &syn, InverseCompton &IC, GammaGamma &gg, Photonbackground &ph,
                        const std::vector<double> &time_array,
                        const std::vector<double> &energy_array_min,
                        const std::vector<double> &energy_array_max) {

  Init();

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
    std::cout << " t " << t << std::endl;

    Spectrum(ED, syn, IC, gg, t);

    for (size_t i = 0; i < num_p; i++) {
      photon_tot[i] = photon_syn[i] + photon_ssc[i];
      photon_tot_diff[i] = photon_tot[i] / photon_energy[i];
    }

    double dum;
    for (size_t i = 0; i < energy_array_min.size(); i++) {
      if (fabs(energy_array_min[i] - energy_array_max[i]) < 1e-5) {
        dum = utility.Interpolate(photon_energy, photon_tot, energy_array_min[i]);
      } else {
        dum = utility.Integrate(photon_energy, photon_tot_diff, energy_array_min[i],
                                energy_array_max[i]);
      }
      flux_vector[idx][i] = dum;
    }

    idx++;
  }
}

void GRBAfterglow::Flux(ElectronDistribution &ED, Synchrotron &syn, InverseCompton &IC,
                        GammaGamma &gg, Photonbackground &ph, const std::vector<double> &time_array,
                        const std::vector<double> &energy_array_min,
                        const std::vector<double> &energy_array_max) {

  Init();

  if (have_onezone) {
  } else {
    InitJet(jet);
  }

  std::ofstream mystream_spectrum;
  mystream_spectrum.open(output_folder + "/spectrum.dat");
  std::ofstream mystream_flux;
  mystream_flux.open(output_folder + "/flux.dat");
  mystream_flux << "# t " << std::endl;

  int idx = 0;
  for (auto t : time_array) {
    std::cout << " t " << t << std::endl;

    Spectrum(ED, syn, IC, gg, t);

    if (have_attenu_GG_cosmic) {
      ph.EBLTau(param.z, ph.EBL_Gilmore12());
      ph.EBLAttenuation(photon_energy, photon_syn);
      ph.EBLAttenuation(photon_energy, photon_ssc);
    }

    for (size_t i = 0; i < num_p; i++) {
      mystream_spectrum << photon_energy[i] << " ";
      mystream_spectrum << photon_syn[i] << " " << photon_ssc[i] << " ";
      mystream_spectrum << std::endl;
      photon_tot[i] = photon_syn[i] + photon_ssc[i];
      photon_tot_diff[i] = photon_tot[i] / photon_energy[i];
    }

    mystream_flux << t << " ";
    for (size_t i = 0; i < energy_array_min.size(); i++) {
      if (fabs(energy_array_min[i] - energy_array_max[i]) < 1e-5) {
        mystream_flux << " " << utility.Interpolate(photon_energy, photon_tot, energy_array_min[i]);
      } else {
        mystream_flux << " "
                      << utility.Integrate(photon_energy, photon_tot_diff, energy_array_min[i],
                                           energy_array_max[i]);
      }
    }
    mystream_flux << std::endl;

    idx++;
  }

  mystream_spectrum.close();
  mystream_flux.close();
}

void GRBAfterglow::Spectrum(ElectronDistribution &ED, Synchrotron &syn, InverseCompton &IC,
                            GammaGamma &gg, double T) {

  double shell_thickness_factor = 4. * (3. - param.k_ex);
  double syn_dum, ssc_dum;
  if (have_onezone) {
    double Gamma, radius;
    EvolutionThinShellAnalytic(T, Gamma, radius);
    double beta = sqrt(1. - 1. / Gamma / Gamma);
    double dynamical_timescale = Gamma * T / (1 + param.z);
    double adiabatic_timescale = dynamical_timescale;
    double Delta = radius / shell_thickness_factor / Gamma; // comoving frame
    double photon_escape_timescale = Delta / c_cnst;

    double mag_strength = sqrt(32 * PI * mp * c_cnst * c_cnst * param.epsilon_B *
                               ExternalDensity(radius) * Gamma * (Gamma - 1));
    s.getParam().setMagStrength(mag_strength);

    double electron_inj =
        (4. * Gamma) * param.fraction_e * ExternalDensity(radius) / dynamical_timescale;

    s.getElectron().setSpectrumPL(MomentumInj(Gamma), MomentumMax(mag_strength, radius, beta),
                                  param.spectral_e, electron_inj, false);
    ED.IterationSolutionSteadyState(syn, IC, false, dynamical_timescale, photon_escape_timescale);
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
      gg.Losstime(target_energy, target_syn, photon_energy, photon_losstime);
      for (size_t i = 0; i < num_p; i++) {
        s.getPhoton().setOptdepth(i, photon_losstime[i] * dynamical_timescale);
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
          photon_syn[k] *= param.open_angle * param.open_angle * Gamma * Gamma;
          photon_ssc[k] *= param.open_angle * param.open_angle * Gamma * Gamma;
        }
      }
    }

  } else {

    EvolutionThinShell(jet, T);

    photon_syn.assign(num_p, 0.0);
    photon_ssc.assign(num_p, 0.0);
    for (size_t i = 0; i < jet.theta.size(); i++) {
      if (param.view_angle == 0) {
        int j = 0;

        double beta = sqrt(1. - 1. / (jet.Gamma[i][j] * jet.Gamma[i][j]));
        double beta_sh = sqrt(1. - 1. / 2 / (jet.Gamma[i][j] * jet.Gamma[i][j]));
        double doppler_factor = 1. / jet.Gamma[i][j] / (1 - beta * jet.mu[i][j]);
        double delta_sh = jet.radius[i][j] / shell_thickness_factor / (1 - beta_sh * jet.mu[i][j]);
        double delta_s = jet.radius[i][j] / shell_thickness_factor / jet.Gamma[i][j] /
                         jet.Gamma[i][j] / fabs(jet.mu[i][j] - beta_sh);

        double EATS_factor = 2 * PI * jet.theta_bin[i] * sin(jet.theta[i]) * jet.radius[i][j] *
                             jet.radius[i][j] * fabs(jet.mu[i][j] - beta_sh) /
                             (1 - beta_sh * jet.mu[i][j]) * (1 + param.z) / param.dl / param.dl;

        double mag_strength =
            sqrt(32 * PI * mp * c_cnst * c_cnst * param.epsilon_B *
                 ExternalDensity(jet.radius[i][j]) * jet.Gamma[i][j] * (jet.Gamma[i][j] - 1));
        s.getParam().setMagStrength(mag_strength);

        double photon_escape_timescale =
            jet.radius[i][j] / shell_thickness_factor / jet.Gamma[i][j] / c_cnst * 2;
        double electron_inj = (4. * jet.Gamma[i][j]) * param.fraction_e *
                              ExternalDensity(jet.radius[i][j]) / jet.dynamical_timescale[i][j];
        s.getElectron().setSpectrumPL(MomentumInj(jet.Gamma[i][j]),
                                      MomentumMax(mag_strength, jet.radius[i][j], beta),
                                      param.spectral_e, electron_inj, false);

        ED.IterationSolution(syn, IC, false, jet.dynamical_timescale[i][j],
                             jet.adiabatic_timescale[i][j], photon_escape_timescale);

        target_syn = s.getTarget().getSpectrum();

        syn.Emissivity(mag_strength, photon_energy, photon_syn_temp);
        for (size_t k = 0; k < num_p; k++) {
          photon_syn_temp[k] *= photon_energy[k];
        }

        if (have_SSCSpec) {
          IC.Emissivity(target_energy, target_syn, photon_energy, photon_ssc_temp);
          for (size_t k = 0; k < num_p; k++) {
            photon_ssc_temp[k] *= photon_energy[k];
          }
        }

        if (have_attenu_SSA) {
          syn.SSAAbsorptionCoeff(mag_strength, photon_syn_absorption_SSA);
        }

        double alpha, tau, att;
        for (size_t k = 0; k < num_p; k++) {
          if (have_attenu_SSA) {
            alpha = utility.Interpolate(photon_energy, photon_syn_absorption_SSA,
                                        photon_energy[k] * (1 + param.z) / doppler_factor) /
                    doppler_factor;
            tau = alpha * delta_s;
            if (tau > 1e-5) {
              att = (1. - exp(-tau)) / tau;
            } else {
              att = 1. / (1 + tau);
            }
          } else {
            att = 1.0;
          }

          photon_syn_obs[k] =
              utility.Interpolate(photon_energy, photon_syn_temp,
                                  photon_energy[k] * (1 + param.z) / doppler_factor) *
              doppler_factor * doppler_factor / (4 * PI) * att * delta_s;
        }

        if (have_attenu_GG_source) {
          gg.Losstime(target_energy, target_syn, photon_energy, photon_losstime);
        }

        for (size_t k = 0; k < num_p; k++) {
          if (have_attenu_GG_source) {
            alpha = photon_losstime[k] / c_cnst;
            alpha = utility.Interpolate(photon_energy, photon_losstime,
                                        photon_energy[k] * (1 + param.z) / doppler_factor) /
                    doppler_factor / c_cnst;
            tau = alpha * delta_s;
            if (tau > 1e-5) {
              att = (1. - exp(-tau)) / tau;
            } else {
              att = 1. / (1 + tau);
            }
          } else {
            att = 1.0;
          }
          photon_ssc_obs[k] =
              utility.Interpolate(photon_energy, photon_ssc_temp,
                                  photon_energy[k] * (1 + param.z) / doppler_factor) *
              doppler_factor * doppler_factor / (4 * PI) * att * delta_s;
        }

        for (size_t k = 0; k < num_p; k++) {
          photon_syn[k] += EATS_factor * photon_energy[k] * eV2erg * photon_syn_obs[k];
          photon_ssc[k] += EATS_factor * photon_energy[k] * eV2erg * photon_ssc_obs[k];
        }

      } else {
        std::cout << "Warning: The current version is only for on-axis jet " << std::endl;
        for (size_t j = 0; j < jet.phi.size(); j++) {
          beta = sqrt(1. - 1. / (jet.Gamma[i][j] * jet.Gamma[i][j]));
          beta_sh = sqrt(1. - 1. / 2 / (jet.Gamma[i][j] * jet.Gamma[i][j]));
          doppler_factor = 1. / jet.Gamma[i][j] / (1 - beta * jet.mu[i][j]);
          delta_s = jet.radius[i][j] / shell_thickness_factor / jet.Gamma[i][j] / jet.Gamma[i][j] /
                    fabs(jet.mu[i][j] - beta_sh);

          EATS_factor = 2 * PI * jet.theta_bin[i] * sin(jet.theta[i]) * jet.radius[i][j] *
                        jet.radius[i][j] * fabs(jet.mu[i][j] - beta_sh) /
                        (1 - beta_sh * jet.mu[i][j]) * (1 + param.z) / param.dl / param.dl;

          mag_strength =
              sqrt(32 * PI * mp * c_cnst * c_cnst * param.epsilon_B *
                   ExternalDensity(jet.radius[i][j]) * jet.Gamma[i][j] * (jet.Gamma[i][j] - 1));
          s.getParam().setMagStrength(mag_strength);

          photon_escape_timescale =
              jet.radius[i][j] / shell_thickness_factor / jet.Gamma[i][j] / c_cnst * 2;
          electron_inj = (4. * jet.Gamma[i][j]) * param.fraction_e *
                         ExternalDensity(jet.radius[i][j]) / jet.dynamical_timescale[i][j];
          s.getElectron().setSpectrumPL(MomentumInj(jet.Gamma[i][j]),
                                        MomentumMax(mag_strength, jet.radius[i][j], beta),
                                        param.spectral_e, electron_inj, false);

          // ELectron cooling without External photons
          ED.IterationSolution(syn, IC, false, jet.dynamical_timescale[i][j],
                               jet.adiabatic_timescale[i][j], photon_escape_timescale);
          target_syn = s.getTarget().getSpectrum();

          syn.Emissivity(mag_strength, photon_energy, photon_syn_temp);
          for (size_t k = 0; k < num_p; k++) {
            photon_syn_temp[k] *= photon_energy[k];
          }

          // SSC process
          if (have_SSCSpec) {
            IC.Emissivity(target_energy, target_syn, photon_energy, photon_ssc_temp);
            for (size_t k = 0; k < num_p; k++) {
              photon_ssc_temp[k] *= photon_energy[k];
            }
          }

          if (have_attenu_SSA) {
            syn.SSAAbsorptionCoeff(mag_strength, photon_syn_absorption_SSA);
          }

          for (size_t k = 0; k < num_p; k++) {
            if (have_attenu_SSA) {
              alpha = utility.Interpolate(photon_energy, photon_syn_absorption_SSA,
                                          photon_energy[k] * (1 + param.z) / doppler_factor) /
                      doppler_factor;
              tau = alpha * delta_s;
              if (tau > 1e-5) {
                att = (1. - exp(-tau)) / tau;
              } else {
                att = 1. / (1 + tau);
              }
            } else {
              att = 1.0;
            }

            photon_syn_obs[k] =
                utility.Interpolate(photon_energy, photon_syn_temp,
                                    photon_energy[k] * (1 + param.z) / doppler_factor) *
                doppler_factor * doppler_factor / (4 * PI) * att * delta_s;
          }

          if (have_attenu_GG_source) {
            gg.Losstime(target_energy, target_tot, photon_energy, photon_losstime);
          }

          for (size_t k = 0; k < num_p; k++) {
            if (have_attenu_GG_source) {
              alpha = photon_losstime[k] / c_cnst;
              alpha = utility.Interpolate(photon_energy, photon_losstime,
                                          photon_energy[k] * (1 + param.z) / doppler_factor) /
                      doppler_factor / c_cnst;
              tau = alpha * delta_s;
              if (tau > 1e-5) {
                att = (1. - exp(-tau)) / tau;
              } else {
                att = 1. / (1 + tau);
              }
            } else {
              att = 1.0;
            }
            photon_ssc_obs[k] =
                utility.Interpolate(photon_energy, photon_ssc_temp,
                                    photon_energy[k] * (1 + param.z) / doppler_factor) *
                doppler_factor * doppler_factor / (4 * PI) * att * delta_s;
          }

          for (size_t k = 0; k < num_p; k++) {
            photon_syn[k] += EATS_factor * photon_energy[k] * eV2erg * photon_syn_obs[k];
            photon_ssc[k] += EATS_factor * photon_energy[k] * eV2erg * photon_ssc_obs[k];
          }
        }
      }
    }
  }
}
}

double GRBAfterglow::MomentumInj(double Gamma) {
  double energy_kin; // kinetic energy
  if (param.spectral_e > 2) {
    energy_kin = proton_mass * param.epsilon_e / param.fraction_e * (param.spectral_e - 2) /
                 (param.spectral_e - 1) * (Gamma - 1);
  } else {
    std::cout << "Warnning: spectral index should larger than 2 !!! \n" << std::endl;
    return 0;
  }
  return s.getElectron().Energy2Momentum(energy_kin + e_mass);
}

double GRBAfterglow::MomentumMax(double mag_strength, double radius, double beta) {
  // t_syn  == t_acc
  double gamma_max_syn = sqrt(6 * PI * e_charge / param.eta_acc_e / sigmaT / mag_strength);

  // t_ad  == t_acc
  double gamma_max_ad =
      e_charge * mag_strength * radius / (param.eta_acc_e * me * c_cnst * c_cnst * beta);

  return s.getElectron().Energy2Momentum(e_mass * std::min(gamma_max_syn, gamma_max_ad));
}

double GRBAfterglow::ExternalDensity(double radius) {
  if (param.k_ex > 0) {
    return param.n_ex * pow(radius, -param.k_ex);
  } else {
    return param.n_ex;
  }
}

void GRBAfterglow::InitJet(Jet &jet) {
  double theta_min = 0.1 / param.Gamma0;
  double theta_max = param.open_angle;
  jet.theta.resize(jet.angular_num);
  jet.theta_bin.resize(jet.angular_num);

  if (param.view_angle == 0) {
    for (size_t i = 0; i < jet.angular_num; i++) {
      jet.theta[i] = 0. + param.open_angle * i / (jet.angular_num - 1);
      jet.theta_bin[i] = param.open_angle / jet.angular_num;
    }
  } else {
    // Calculate the logarithmic spacing factor
    double log_factor = (std::log10(theta_max) - std::log10(theta_min)) / (jet.angular_num - 1);

    // Generate the logarithmically spaced grid

    for (int i = 0; i < jet.angular_num; ++i) {
      jet.theta[i] = std::pow(10, std::log10(theta_min) + i * log_factor);
    }
    double theta_l = std::pow(10, std::log10(theta_min) + (-1) * log_factor);
    double theta_r = std::pow(10, std::log10(theta_min) + jet.angular_num * log_factor);
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
  jet.mu.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
  jet.dynamical_timescale.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
  jet.adiabatic_timescale.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));

  for (size_t i = 0; i < jet.theta.size(); i++) {
    for (size_t j = 0; j < jet.phi.size(); j++) {
      jet.mu[i][j] = sin(jet.theta[i]) * sin(param.view_angle) * cos(jet.phi[j]) +
                     cos(jet.theta[i]) * cos(param.view_angle);
    }
  }
}

void GRBAfterglow::EvolutionThinShell(Jet &jet, double T) {
  const int n = 7;
  std::vector<double> y(n);
  std::vector<double> dydx(n);
  std::vector<double> yout(n);
  std::vector<double> shock_dynamic(n);

  double epsilon = 0.0;
  double C_BM76 = epsilon + (9.0 - 2.0 * param.k_ex) / (17.0 - 4.0 * param.k_ex) *
                                (1.0 - epsilon); // correction factor
  double E_ej, M_ej, Gamma, beta;

  // Initial values
  for (size_t i = 0; i < jet.theta.size(); i++) {
    // structured jet
    if ((param.gaussian_cone > 0) && (param.jet_index > 0)) {
      E_ej = param.E_ej *
             pow(1 + jet.theta[i] / param.gaussian_cone * jet.theta[i] / param.gaussian_cone,
                 -param.jet_index / 2) *
             (2 * PI * (1 - cos(param.open_angle)) / (4 * PI));
      Gamma = param.Gamma0 *
              pow(1 + jet.theta[i] / param.gaussian_cone * jet.theta[i] / param.gaussian_cone,
                  -param.jet_index / 2);
    } else if ((param.gaussian_cone > 0) && (param.jet_index <= 0)) {
      E_ej = param.E_ej *
             exp(-0.5 * jet.theta[i] * jet.theta[i] / param.gaussian_cone / param.gaussian_cone) *
             (2 * PI * (1 - cos(param.open_angle)) / (4 * PI));
      Gamma = param.Gamma0;
    } else {
      E_ej = param.E_ej * (2 * PI * (1 - cos(param.open_angle)) / (4.0 * PI));
      Gamma = param.Gamma0;
    }
    M_ej = E_ej / (Gamma * c_cnst * c_cnst);
    beta = sqrt(1 - 1. / Gamma / Gamma);

    jet.Gamma0[i] = Gamma;
    jet.M_ej[i] = M_ej;
    jet.E_ej[i] = E_ej;

    DynamicThinShell dyn(M_ej, C_BM76 * param.n_ex, param.k_ex, epsilon);

    y[0] = 1e10; // assume stellar radius
    y[1] = 4. / 3 * PI * y[0] * y[0] * y[0] * mp * C_BM76 * param.n_ex * pow(y[0], -param.k_ex);
    y[2] = Gamma;
    y[3] = param.open_angle;
    y[4] = (Gamma - 1) * y[1] * c_cnst * c_cnst;
    y[5] = y[0] / y[2] / c_cnst;
    y[6] = y[0] / y[2] / c_cnst;

    double x = 0;
    dyn(x, y, dydx);
    rk5<DynamicThinShell> evol(dyn, y, dydx, yout, x, n);

    if (param.view_angle == 0) {
      shock_dynamic = evol.integrateEATS(T / (1 + param.z), jet.mu[i][0]);
      for (size_t j = 0; j < jet.phi.size(); j++) {
        jet.radius[i][j] = shock_dynamic[0];
        jet.Gamma[i][j] = shock_dynamic[2];
        jet.dynamical_timescale[i][j] = shock_dynamic[5];
        jet.adiabatic_timescale[i][j] = shock_dynamic[6];
      }
    } else {
      for (size_t j = 0; j < jet.phi.size(); j++) {
        shock_dynamic = evol.integrateEATS(T / (1 + param.z), jet.mu[i][j]);
        jet.radius[i][j] = shock_dynamic[0];
        jet.Gamma[i][j] = shock_dynamic[2];
        jet.dynamical_timescale[i][j] = shock_dynamic[5];
        jet.adiabatic_timescale[i][j] = shock_dynamic[6];
      }
    }
  }
}

// dynamics with FS
DynamicThinShell::DynamicThinShell(const double M_ej, const double n_ex, const double k_ex,
                                   const double epsilon)
    : M_ej(M_ej), n_ex(n_ex), k_ex(k_ex), epsilon(epsilon) {}

void DynamicThinShell::operator()(const double x, std::vector<double> &y,
                                  std::vector<double> &dydx) {
  beta = sqrt(1 - 1. / y[2] / y[2]);
  gamma_hat = (4. + 1. / y[2]) / 3.;
  Gamma_eff = (gamma_hat * y[2] * y[2] - gamma_hat + 1) / y[2];
  dGamma_eff = (gamma_hat * y[2] * y[2] + gamma_hat - 1) / (y[2] * y[2]);
  A = 2. * PI * (1 - cos(y[3])) * y[0] * y[0] * mp * n_ex * pow(y[0], -k_ex);
  // double cs = sqrt(gamma_hat * (gamma_hat - 1) * (y[2] - 1) / (1 + gamma_hat * (y[2] - 1)) *
  //                  c_cnst * c_cnst);

  // dr/dt
  dydx[0] = beta * c_cnst;

  // dm/dt = dm/dr * dr/dt
  dydx[1] = A * dydx[0];

  // dEad / dr
  dEad = -(gamma_hat - 1) * (3. / y[0] - 1 / y[2] * dydx[2] / dydx[0]) * y[4];

  // dEint/dt = dEint / dr * dr / dt
  dydx[4] = ((1.0 - epsilon) * (y[2] - 1.0) * A * c_cnst * c_cnst + dEad) * dydx[0];

  // dGamma/dt = dGamma /dr * dr/dt
  dydx[2] =
      -(A * y[2] * (y[2] * y[2] - 1) * (gamma_hat * y[2] - gamma_hat + 1) * c_cnst * c_cnst -
        (gamma_hat - 1) * y[2] * (gamma_hat * y[2] * y[2] - gamma_hat + 1) * 3 * y[4] / y[0]) /
      (y[2] * y[2] * (M_ej + y[1]) * c_cnst * c_cnst +
       (gamma_hat * gamma_hat * y[2] * y[2] - gamma_hat * gamma_hat + 3 * gamma_hat - 2) * y[4]) *
      dydx[0];

  // dtheta/ dt
  // dydx[3] = cs / c_cnst / y[0] / y[2] * dydx[0];
  dydx[3] = 0.0;

  dydx[5] = 1.0 / c_cnst / beta / y[2] * dydx[0]; // comving dynamical timescale

  dydx[6] = 0.0;
  y[6] = 1.0 / (c_cnst * beta * y[2] *
                (1.0 / y[0] - 1.0 / 3.0 / y[2] * dydx[2] / dydx[0])); // comving adiabatic timescale
}

double GRBAfterglow::DecelerationRadius() {
  double M_ej = param.E_ej / (param.Gamma0 * c_cnst * c_cnst);
  return pow((3 - param.k_ex) * M_ej / (4 * PI * param.n_ex * mp * param.Gamma0),
             1 / (3 - param.k_ex));
}

double GRBAfterglow::DecelerationRadius(double E_ej, double Gamma) {
  double M_ej = E_ej / (Gamma * c_cnst * c_cnst);
  return pow((3 - param.k_ex) * M_ej / (4 * PI * param.n_ex * mp * Gamma), 1 / (3 - param.k_ex));
}

void GRBAfterglow::EvolutionThinShellAnalytic(double T, double &Gamma, double &radius) {
  Gamma = 43.6 * exp(0.125 * log(param.E_ej / 1e53 / param.n_ex)) * exp(-0.375 * log(T / 1e3)) *
          exp(0.375 * log(1 + param.z)); // Sari98

  if (Gamma < 1.00001) {
    Gamma = 1.00001;
  }
  if (Gamma > param.Gamma0) {
    Gamma = param.Gamma0;
  }
  radius = 4. * Gamma * Gamma * c_cnst * T / (1 + param.z); // Sari98
}

void GRBAfterglow::LorentzBoost(const std::vector<double> &photon_energy,
                                std::vector<double> &photon_spectrum, double Gamma) {
  int num_p = photon_energy.size();
  std::vector<double> photon_temp(num_p, 0.0);
  for (size_t i = 0; i < num_p; i++) {
    photon_temp[i] = photon_energy[i] * photon_spectrum[i];
  }
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] =
        (1 + param.z) * Gamma *
        utility.Interpolate(photon_energy, photon_temp, (1 + param.z) * photon_energy[i] / Gamma);
  }
}

void GRBAfterglow::FluxNorm(const std::vector<double> &photon_energy,
                            std::vector<double> &photon_spectrum, double radius, double Delta) {
  double flux_norm = eV2erg * (4.0 * PI * radius * radius * Delta) / (4 * PI * param.dl * param.dl);
  int num_p = photon_energy.size();
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] *= flux_norm * photon_energy[i];
  }
}
