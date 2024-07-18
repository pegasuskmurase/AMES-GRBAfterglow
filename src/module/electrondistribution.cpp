//==============================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file electrondistribution.hpp
// \brief numerically determine electron distribution

#include "module/electrondistribution.hpp"
#include "module/constants.hpp"
#include "module/utility.hpp"
#include <fstream>

ElectronDistribution::ElectronDistribution(Source &s) : s(s) {}

void ElectronDistribution::IterationSolutionSteadyState(Synchrotron &syn, InverseCompton &IC,
                                                        bool have_EIC, double t_dyn,
                                                        double t_esc_photon) {

  // treament of steady-state solution, original method
  double ERR = 1e-4;
  const int N = s.electron.getSize();
  std::vector<double> x = s.electron.momentum;
  std::vector<double> x_e = s.electron.energy;
  std::vector<double> Q_inj = s.electron.spectrum; // momentum space
  std::vector<double> e_loss_Syn(N);
  std::vector<double> e_loss_EIC(N);
  std::vector<double> e_loss_SSC(N);
  std::vector<double> e_loss_tot(N);
  syn.Losstime(s.param.getMagStrength(), x_e, e_loss_Syn);

  if (have_EIC) {
    IC.Losstime(s.target.energy, s.target.spectrum, x_e, e_loss_EIC);
  }
  for (size_t j = 0; j < N; j++) {
    e_loss_tot[j] = e_loss_EIC[j] + e_loss_Syn[j];
  }

  double momentum_inj;
  for (size_t j = 0; j < N; j++) {
    if (Q_inj[j] > 0.0) {
      momentum_inj = s.electron.momentum[j];
      break;
    }
  }

  double dum, momentum_cooling;
  int idx;
  for (size_t j = 0; j < N; j++) {
    dum = e_loss_tot[j] * t_dyn;
    if (dum > 1.) {
      momentum_cooling = s.electron.momentum[j];
      idx = j;
      break;
    }
  }

  double e_loss_inj = utility.Interpolate(x, e_loss_tot, momentum_inj);
  double momentum_eff;
  for (size_t j = 0; j < N; j++) {
    momentum_eff = momentum_inj * momentum_cooling / (momentum_inj + momentum_cooling);
    if (x[j] > momentum_eff) {
      dum = utility.Integrate(x, Q_inj, x[j]) / ((1. / t_dyn + e_loss_tot[j]) * x_e[j]);
    } else {
      dum = 0.0;
    }
    s.electron.setSpectrum(j, dum);
  }

  int num_t = s.target.energy.size();
  std::vector<double> target_syn(num_t);

  std::vector<double> spectrum_temp(N);
  int istep = 0;
  double u_err, num_tot, num_tot_inj;
  do {
    syn.Emissivity(s.param.getMagStrength());
    for (size_t i = 0; i < num_t; i++) {
      target_syn[i] = s.photon.Interpolate(s.target.energy[i]) * t_esc_photon;
    }

    IC.Losstime(s.target.energy, target_syn, x_e, e_loss_SSC);
    u_err = 0.0;
    for (size_t j = 0; j < N; j++) {
      e_loss_tot[j] = e_loss_EIC[j] + e_loss_Syn[j] + e_loss_SSC[j];
    }
    for (size_t j = 0; j < N; j++) {
      dum = e_loss_tot[j] * t_dyn;
      if (dum > 1.) {
        momentum_cooling = s.electron.momentum[j];
        break;
      }
    }
    double e_loss_inj = utility.Interpolate(x, e_loss_tot, momentum_inj);

    double dum;
    double momentum_eff;
    for (size_t j = 0; j < N; j++) {
      momentum_eff = momentum_inj * momentum_cooling / (momentum_inj + momentum_cooling);
      if (x[j] > momentum_eff) {
        dum = utility.Integrate(x, Q_inj, x[j]) / ((1. / t_dyn + e_loss_tot[j]) * x_e[j]);
      } else {
        dum = 0.0;
      }
      spectrum_temp[j] = dum;
      if (s.electron.spectrum[j] > 0) {
        u_err += abs(s.electron.spectrum[j] - spectrum_temp[j]) / s.electron.spectrum[j];
      }
    }
    for (size_t j = 0; j < N; j++) {
      s.electron.setSpectrum(j, spectrum_temp[j]);
    }
    istep++;
    if (istep > 10) {
      std::cout << "u_err " << istep << " " << u_err << " " << ERR << " " << momentum_inj << " "
                << momentum_cooling << " " << momentum_eff << std::endl;

      break;
    }
  } while (u_err > ERR);

  for (size_t i = 0; i < num_t; i++) {
    s.target.setSpectrum(i, target_syn[i]);
  }
}

void ElectronDistribution::IterationSolution(Synchrotron &syn, InverseCompton &IC, bool have_EIC,
                                             double t_dyn, double t_ad, double t_esc_photon) {

  // treament of quasi steady-state
  double ERR = 1e-5;
  const int N = s.electron.getSize();
  std::vector<double> x = s.electron.momentum;
  std::vector<double> x_e = s.electron.energy;
  std::vector<double> Q_inj = s.electron.spectrum; // momentum space
  std::vector<double> e_loss_Syn(N);
  std::vector<double> e_loss_EIC(N);
  std::vector<double> e_loss_SSC(N);
  std::vector<double> e_loss_tot(N);
  syn.Losstime(s.param.getMagStrength(), x_e, e_loss_Syn);

  if (have_EIC) {
    IC.Losstime(s.target.energy, s.target.spectrum, x_e, e_loss_EIC);
    for (size_t j = 0; j < N; j++) {
      e_loss_tot[j] = e_loss_EIC[j] + e_loss_Syn[j];
    }
  } else {

    for (size_t j = 0; j < N; j++) {
      e_loss_tot[j] = e_loss_Syn[j];
    }
  }

  double momentum_inj;
  for (size_t j = 0; j < N; j++) {
    if (Q_inj[j] > 0.0) {
      momentum_inj = s.electron.momentum[j];
      break;
    }
  }

  double dum, momentum_cooling;
  int idx;
  for (size_t j = 0; j < N; j++) {
    dum = e_loss_tot[j] * t_dyn;
    if (dum > 1.) {
      momentum_cooling = s.electron.momentum[j];
      idx = j;
      break;
    }
  }

  double t_eff = 0.5 * t_ad / t_dyn;
  if (t_dyn < t_ad) {
    t_eff = 0.5 * t_dyn / t_ad;
  }

  double e_loss_inj = utility.Interpolate(x, e_loss_tot, momentum_inj);
  double momentum_eff;
  for (size_t j = 0; j < N; j++) {
    momentum_eff = momentum_inj * momentum_cooling / (momentum_inj + momentum_cooling);
    if (x[j] > momentum_eff) {
      dum = utility.Integrate(x, Q_inj, x[j]) /
            ((1. / t_dyn + (1 - t_eff) / t_ad + e_loss_tot[j]) * x_e[j]);
    } else {
      dum = 0.0;
    }
    s.electron.setSpectrum(j, dum);
  }

  int num_t = s.target.energy.size();
  std::vector<double> target_syn(num_t);

  std::vector<double> spectrum_temp(N);
  int istep = 0;
  double u_err, num_tot, num_tot_inj;
  do {
    syn.Emissivity(s.param.getMagStrength());
    for (size_t i = 0; i < num_t; i++) {
      target_syn[i] = s.photon.Interpolate(s.target.energy[i]) * t_esc_photon;
    }

    IC.Losstime(s.target.energy, target_syn, x_e, e_loss_SSC);
    u_err = 0.0;
    for (size_t j = 0; j < N; j++) {
      e_loss_tot[j] = e_loss_EIC[j] + e_loss_Syn[j] + e_loss_SSC[j];
    }

    for (size_t j = 0; j < N; j++) {
      dum = e_loss_tot[j] * t_dyn;
      if (dum > 1.) {
        momentum_cooling = s.electron.momentum[j];
        break;
      }
    }
    double e_loss_inj = utility.Interpolate(x, e_loss_tot, momentum_inj);

    double dum;
    double momentum_eff;
    for (size_t j = 0; j < N; j++) {
      momentum_eff = momentum_inj * momentum_cooling / (momentum_inj + momentum_cooling);
      if (x[j] > momentum_eff) {
        dum = utility.Integrate(x, Q_inj, x[j]) /
              ((1. / t_dyn + (1 - t_eff) / t_ad + e_loss_tot[j]) * x_e[j]);
      } else {
        dum = 0.0;
      }
      spectrum_temp[j] = dum;
      if (s.electron.spectrum[j] > 0) {
        u_err += abs(s.electron.spectrum[j] - spectrum_temp[j]) / s.electron.spectrum[j];
      }
    }
    for (size_t j = 0; j < N; j++) {
      s.electron.setSpectrum(j, spectrum_temp[j]);
    }
    istep++;
    if (istep > 10) {
      std::cout << "u_err " << istep << " " << u_err << " " << ERR << " " << momentum_inj << " "
                << momentum_cooling << " " << momentum_eff << std::endl;

      break;
    }
  } while (u_err > ERR);

  for (size_t i = 0; i < num_t; i++) {
    s.target.setSpectrum(i, target_syn[i]);
  }
}

void ElectronDistribution::Test() {}
