//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file synchrotron.hpp
// \brief Class for synchrotron process

#include "module/constants.hpp"
#include "module/particle.hpp"

#include <vector>

#ifndef SYNCHROTRON_H
#define SYNCHROTRON_H

class Source;

class Synchrotron {
  /*
  @class: Synchrotron
  @brief: basic class for Synchrotron emission process
  */
public:
  Synchrotron(Source &s);

  void Transport(std::vector<ParticleTransport> &p, const double &mag_strength,
                 bool &have_SSAAttenuation);
  void Transport(std::vector<ParticleTransportChangCooper> &p,
                 const double &mag_strength, bool &have_SSAAttenuation);

  // energy loss time in units [s^-1]
  // void Losstime();
  void Losstime(const double mag_strength);
  void Losstime(const double mag_strength,
                const std::vector<double> &electron_energy,
                std::vector<double> &electron_losstime);
  void Losstime(const double mag_strength, const int particle_charge,
                const double particle_mass,
                const std::vector<double> &particle_energy,
                std::vector<double> &particle_losstime);

  void LosstimeProton(const double mag_strength);

  // radiation emissivity in units [eV^-1 s^-1 cm^-3]
  void Emissivity();
  void Emissivity(double mag_strength);
  void Emissivity(double mag_strength, std::vector<double> &photon_energy,
                  std::vector<double> &photon_spectrum);
  void Emissivity(const std::vector<std::vector<double>> &spec,
                  const std::vector<double> &primary_momentum,
                  const std::vector<double> &primary_spectrum,
                  const std::vector<double> &secondary_energy,
                  std::vector<double> &secondary_spectrum);
  void Emissivity(double mag_strength,
                  const std::vector<double> &primary_energy,
                  const std::vector<double> &primary_spectrum,
                  const std::vector<double> &secondary_energy,
                  std::vector<double> &secondary_spectrum);
  void EmissivityProton(double mag_strength);
  // spec in units [eV^-1 s^-1]
  void Spec(const double mag_strength,
            const std::vector<double> &primary_energy,
            const std::vector<double> &secondary_energy);
  void Spec(const double mag_strength,
            const std::vector<double> &primary_energy,
            const std::vector<double> &secondary_energy,
            std::vector<std::vector<double>> &spec);
  void Spec(const double mag_strength, const double primary_charge,
            const double primary_mass,
            const std::vector<double> &primary_energy,
            const std::vector<double> &secondary_energy,
            std::vector<std::vector<double>> &spec);
  void SpecImplicit(const double mag_strength, const double primary_charge,
                    const double primary_mass,
                    const std::vector<double> &primary_energy,
                    const std::vector<double> &secondary_energy,
                    const std::vector<double> &secondary_energy_bin,
                    std::vector<std::vector<double>> &spec);

  void SpecProton(const double mag_strength,
                  const std::vector<double> &primary_energy,
                  const std::vector<double> &secondary_energy);

  void SSAOpticalDepth(const double mag_strength, const double emission_size);
  void SSAAbsorptionCoeff(const double mag_strength,
                          std::vector<double> &photon_syn_absorption_SSA);
  void Attenuation(const std::vector<double> &photon_energy,
                   std::vector<double> &photon_spectrum);
  void Attenuation(std::vector<double> &photon_spectrum);

  void IntetimeSSA(const double mag_strength,
                   const std::vector<double> &electron_momentum,
                   const std::vector<double> &electron_energy,
                   const std::vector<double> &electron_spectrum,
                   const std::vector<double> &photon_energy,
                   std::vector<double> &photon_intetime);

  std::vector<std::vector<double>> spec;

  void Test();

private:
  Source &s;
  Utility u;

  double SSACross(const double mag_strength, double x, double y);
  double mBessel53(double x);
};

class myFunc {
public:
  static double f(double x, void *params) {
    std::vector<double> alpha = *(std::vector<double> *)params;
    double f = log(alpha[0] * x) / sqrt(x);
    return f;
  }

  static double f_diff(double x, void *params) {
    (void)(params);
    return pow(x, 1.5);
  }

  static double SynSpec(double x, void *params) {
    double coeff_emissivity = *(double *)params;
    return coeff_emissivity * exp(-x) /
           sqrt(exp((2.0 / 3.0) * log(1.0 / x)) + (3.62 / PI) * (3.62 / PI));
  }

  static double SynSpec2(double x, void *params) {
    std::vector<double> alpha = *(std::vector<double> *)params;
    double coeff_critical_energy = alpha[0];
    double secondary_energy = alpha[1];
    double critical_energy = coeff_critical_energy * x * x;
    double xc = secondary_energy / critical_energy;
    return x * x * exp(-xc) /
           sqrt(exp((2.0 / 3.0) * log(1.0 / xc)) + (3.62 / PI) * (3.62 / PI));
  }

  static double SynSpecInte(double x, void *params) {
    std::vector<double> alpha = *(std::vector<double> *)params;
    double mag_strength = alpha[0];
    double p_inj = alpha[1];
    double p_cut = alpha[2];
    double p_index = alpha[3];
    double norm = alpha[4];
    double secondary_energy = alpha[5];
    double coeff_critical_energy = 1.7365215856405376e-08 * mag_strength;
    double coeff_emissivity = 64039.416604913444 * mag_strength;
    double gamma = x;
    if (gamma <= 1.0)
      return 0.0;
    if (x < secondary_energy)
      return 0.0;
    if (x < p_inj)
      return 0.0;
    if (x > p_cut)
      return 0.0;
    double gamma2beta2 = gamma * gamma - 1.0;
    double critical_energy = coeff_critical_energy * gamma2beta2;
    double xc = secondary_energy / critical_energy;
    if (xc > 300)
      return 0.0;

    double f = norm * exp((-p_index) * log(x)) *
               SynSpec(xc, &coeff_emissivity) / secondary_energy;
    return f;
  }
};

#endif
