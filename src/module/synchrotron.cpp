//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file synchrotron.hpp
// \brief Basic class for synchrotron process

#include "module/synchrotron.hpp"
#include "module/constants.hpp"
#include "module/utility.hpp"
#include <fstream>

Synchrotron::Synchrotron(Source &s) : s(s) {}

void Synchrotron::Losstime() {
  Losstime(s.param.getMagStrength(), s.electron.energy, s.electron.losstime);
}

void Synchrotron::Losstime(const double mag_strength) {
  Losstime(mag_strength, s.electron.energy, s.electron.losstime);
}

void Synchrotron::Losstime(const double mag_strength, const std::vector<double> &electron_energy,
                           std::vector<double> &electron_losstime) {
  double coeff = 4.0 / 3.0 * T_c_cnst; // coefficient
  double B = mag_strength;
  double U_B = B * B / (8.0 * PI) * erg2eV;
  double gamma;
  double beta2gamma2; // beta * beta * gamma * gamma
  int num_e = electron_energy.size();
  electron_losstime.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    gamma = electron_energy[i] / e_mass;
    if (gamma > 1) {
      beta2gamma2 = gamma * gamma - 1.;
      electron_losstime[i] = coeff * beta2gamma2 * U_B / (electron_energy[i] - e_mass);
    } else {
      electron_losstime[i] = 0;
    }
  }
}

void Synchrotron::Losstime(const double mag_strength, const int particle_charge,
                           const double particle_mass, const std::vector<double> &particle_energy,
                           std::vector<double> &particle_losstime) {
  double coeff = 4.0 / 3.0 * T_c_cnst * particle_charge * particle_charge * particle_charge *
                 particle_charge * (e_mass * e_mass / particle_mass / particle_mass); // coefficient
  double B = mag_strength;
  double U_B = B * B / (8.0 * PI) * erg2eV;
  double gamma;
  double beta2gamma2; // beta * beta * gamma * gamma
  int num_p = particle_energy.size();
  particle_losstime.resize(num_p);
  for (size_t i = 0; i < num_p; i++) {
    gamma = particle_energy[i] / particle_mass;
    if (gamma > 1) {
      beta2gamma2 = gamma * gamma - 1.;
      particle_losstime[i] = coeff * beta2gamma2 * U_B / (particle_energy[i] - particle_mass);
    } else {
      particle_losstime[i] = 0;
    }
  }
}

void Synchrotron::Spec(const double mag_strength, const std::vector<double> &primary_energy,
                       const std::vector<double> &secondary_energy) {
  Spec(mag_strength, primary_energy, secondary_energy, spec);
}

void Synchrotron::Spec(const double mag_strength, const std::vector<double> &primary_energy,
                       const std::vector<double> &secondary_energy,
                       std::vector<std::vector<double>> &spec) {
  // double coeff_critical_energy =
  //     1.5 * hbar_cnst * e_charge * c_cnst / (e_mass * e_cnst * e_cnst) *
  //     mag_strength;
  double coeff_critical_energy = 1.7365215856405376e-08 * mag_strength;
  // double coeff_emissivity = 1.81 * sqrt(3.0) * e_charge * e_charge * e_charge
  // / (hbar_cnst * 2.0 * PI * e_mass * e_cnst) * mag_strength;
  double coeff_emissivity = 64039.416604913444 * mag_strength;

  double critical_energy;
  double gamma, gamma2beta2, x;
  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].resize(num_p);
  }
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      gamma = primary_energy[j] / e_mass;
      gamma2beta2 = gamma * gamma - 1.0;
      if ((secondary_energy[i] < primary_energy[j]) && (gamma > 2.0)) {
        critical_energy = coeff_critical_energy * gamma2beta2;
        x = secondary_energy[i] / critical_energy;
        spec[i][j] = coeff_emissivity * exp(-x) /
                     sqrt(exp((2.0 / 3.0) * log(1.0 / x)) + (3.62 / PI) * (3.62 / PI)) /
                     secondary_energy[i];
      } else {
        spec[i][j] = 0.0;
      }
    }
  }
}

void Synchrotron::Spec(const double mag_strength, const double primary_charge,
                       const double primary_mass, const std::vector<double> &primary_energy,
                       const std::vector<double> &secondary_energy,
                       std::vector<std::vector<double>> &spec) {
  double coeff_critical_energy = 1.5 * hbar_cnst * primary_charge * e_charge * c_cnst /
                                 (primary_mass * e_cnst * e_cnst) * mag_strength;
  double coeff_emissivity = 1.81 * sqrt(3.0) * primary_charge * primary_charge * primary_charge *
                            e_charge * e_charge * e_charge /
                            (hbar_cnst * 2.0 * PI * primary_mass * e_cnst) * mag_strength;

  double critical_energy;
  double gamma, gamma2beta2, x;
  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].resize(num_p);
  }
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      gamma = primary_energy[j] / primary_mass;
      gamma2beta2 = gamma * gamma - 1.0;
      if ((secondary_energy[i] < primary_energy[j]) && (gamma > 2.0)) {
        critical_energy = coeff_critical_energy * gamma2beta2;
        x = secondary_energy[i] / critical_energy;
        spec[i][j] = coeff_emissivity * exp(-x) /
                     sqrt(exp((2.0 / 3.0) * log(1.0 / x)) + (3.62 / PI) * (3.62 / PI)) /
                     secondary_energy[i];
      } else {
        spec[i][j] = 0.0;
      }
    }
  }
}

void Synchrotron::SpecProton(const double mag_strength, const std::vector<double> &primary_energy,
                             const std::vector<double> &secondary_energy) {
  double coeff_critical_energy =
      1.5 * hbar_cnst * e_charge * c_cnst / (proton_mass * e_cnst * e_cnst) * mag_strength;
  double coeff_emissivity = 1.81 * sqrt(3.0) * e_charge * e_charge * e_charge /
                            (hbar_cnst * 2.0 * PI * proton_mass * e_cnst) * mag_strength;
  double critical_energy;
  double gamma, gamma2beta2, x;
  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].resize(num_p);
  }
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      gamma = primary_energy[j] / proton_mass;
      gamma2beta2 = gamma * gamma - 1.0;
      if ((secondary_energy[i] < primary_energy[j]) && (gamma > 2.0)) {
        critical_energy = coeff_critical_energy * gamma2beta2;
        x = secondary_energy[i] / critical_energy;
        spec[i][j] = coeff_emissivity * exp(-x) /
                     sqrt(exp((2.0 / 3.0) * log(1.0 / x)) + (3.62 / PI) * (3.62 / PI)) /
                     secondary_energy[i];
      } else {
        spec[i][j] = 0.0;
      }
    }
  }
}

void Synchrotron::Emissivity(const std::vector<std::vector<double>> &spec,
                             const std::vector<double> &primary_momentum,
                             const std::vector<double> &primary_spectrum,
                             const std::vector<double> &secondary_energy,
                             std::vector<double> &secondary_spectrum) {
  int num_s = secondary_energy.size();
  int num_p = primary_momentum.size();
  std::vector<double> primary_temp(num_p);
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      primary_temp[j] = primary_spectrum[j] * spec[i][j];
    }
    secondary_spectrum[i] = u.Integrate(primary_momentum, primary_temp);
  }
}

void Synchrotron::Emissivity() {
  Spec(s.param.getMagStrength(), s.electron.energy, s.photon.energy);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, s.photon.energy, s.photon.spectrum);
}

void Synchrotron::Emissivity(double mag_strength) {
  Spec(mag_strength, s.electron.energy, s.photon.energy);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, s.photon.energy, s.photon.spectrum);
}

void Synchrotron::Emissivity(double mag_strength,
                             std::vector<double> &photon_energy,
                             std::vector<double> &photon_spectrum) {
  Spec(mag_strength, s.electron.energy, photon_energy);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, photon_energy, photon_spectrum);
}

void Synchrotron::EmissivityProton(double mag_strength) {
  SpecProton(mag_strength, s.proton.energy, s.photon.energy);
  Emissivity(spec, s.proton.momentum, s.proton.spectrum, s.photon.energy, s.photon.spectrum);
}

void Synchrotron::Attenuation(const std::vector<double> &photon_energy,
                              std::vector<double> &photon_spectrum) {
  double tau, att;
  int num_p = photon_spectrum.size();
  for (size_t i = 0; i < num_p; i++) {
    tau = u.Interpolate(s.photon.energy, s.photon.optdepth, photon_energy[i]);
    if (tau > 1e-5) {
      att = (1. - exp(-tau)) / tau;
    } else {
      att = 1. / (1 + tau);
    }
    photon_spectrum[i] *= att;
  }
}

void Synchrotron::Attenuation(std::vector<double> &photon_spectrum) {
  double tau, att;
  int num_p = photon_spectrum.size();
  for (size_t i = 0; i < num_p; i++) {
    tau = s.photon.optdepth[i];
    if (tau > 1e-5) {
      att = (1. - exp(-tau)) / tau;
    } else {
      att = 1. / (1 + tau);
    }
    photon_spectrum[i] *= att;
  }
}

void Synchrotron::SSAAbsorptionCoeff(const double mag_strength, std::vector<double> &photon_syn_absorption_SSA) {
  int num_p = s.photon.momentum.size();
  int num_e = s.electron.momentum.size();
  for (size_t i = 0; i < num_p; i++) {
    for (size_t j = 0; j < num_e; j++) {
      if (s.electron.energy[j] > 1 * e_mass) {
        s.electron.temp[j] = s.electron.spectrum[j] *
                             SSACross(mag_strength, s.photon.energy[i], s.electron.energy[j]);
      } else {
        s.electron.temp[j] = 0.0;
      }
    }
    photon_syn_absorption_SSA[i] = u.Integrate(s.electron.momentum, s.electron.temp);
  }
}

void Synchrotron::SSAOpticalDepth(const double emission_size, const double mag_strength) {
  int num_p = s.photon.momentum.size();
  int num_e = s.electron.momentum.size();
  for (size_t i = 0; i < num_p; i++) {
    for (size_t j = 0; j < num_e; j++) {
      if (s.electron.energy[j] > 1 * e_mass) {
        s.electron.temp[j] = s.electron.spectrum[j] *
                             SSACross(mag_strength, s.photon.energy[i], s.electron.energy[j]);
      } else {
        s.electron.temp[j] = 0.0;
      }
    }
    s.photon.optdepth[i] = emission_size * u.Integrate(s.electron.momentum, s.electron.temp);
  }
}

double Synchrotron::SSACross(const double mag_strength, double x, double y) {
  // From Eq 2.14 of the reference Ghisellini and Svensson 1991,
  // the pitch angle theta = pi / 2 is assumed
  y /= e_mass;
  double cycenergy;
  if ((x < y * e_mass) && (y > 1.)) {
    cycenergy =
        (1.5 * hbar_cnst * e_charge * mag_strength * y * y * c_cnst / e_mass / e_cnst / e_cnst) *
        sqrt(1 - 1. / y / y) * sqrt(1 - 1. / y / y); // eV unit
    return (2 * PI / sqrt(3.) / 4. / PI) * 137.036 * 6.652e-25 * 4.413e13 *
           mBessel53(x / cycenergy) / mag_strength / pow(y, 5);
  } else {
    return 0;
  }

  /*
  // From Eq 2.17 of the reference Ghisellini and Svensson, 1991
  // the pitch angle theta = pi / 2 is assumed
  y /= e_mass;
  double cycenergy;
  if ((x < y * e_mass) && (y > 1.)) {
    cycenergy =
        (1.5 * hbar_cnst * e_charge * mag_strength * y * y * c_cnst / e_mass / e_cnst / e_cnst) *
        sqrt(1 - 1. / y / y) * sqrt(1 - 1. / y / y); // eV unit
    return (sqrt(3.) * PI / 10.0) / 4.0 / PI * 137.036 * 6.652e-25 * 4.413e13 * x / cycenergy *
           (boost::math::cyl_bessel_k(4.0 / 3.0, 0.5 * x / cycenergy) *
                boost::math::cyl_bessel_k(4.0 / 3.0, 0.5 * x / cycenergy) -
            boost::math::cyl_bessel_k(1.0 / 3.0, 0.5 * x / cycenergy) *
                boost::math::cyl_bessel_k(1.0 / 3.0, 0.5 * x / cycenergy)) /
           mag_strength / pow(y, 5);
  } else {
    return 0;
  }
  */
}

double Synchrotron::mBessel53(double x) {
  double dummy1 = -1.0194198041210243 * x + 0.28011396300530672 * sqrt(x) -
                  7.71058491739234908 * 0.01 * pow(x, (1. / 3.));
  double dummy2 = -15.761577796582387 * x;
  if (x >= 0) {
    return 0.5 * 0.902745 * pow(0.5 * x, -(5. / 3.)) * exp(dummy1) +
           (1. - exp(dummy2)) * sqrt(0.5 * PI / x) * exp(-x);
  } else {
    return 0;
  }
}

void Synchrotron::Test() { std::cout << "Test : Synchrotron module         " << std::endl; }
