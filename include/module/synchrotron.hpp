//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file synchrotron.hpp
// \brief Class for synchrotron process

#include "module/source.hpp"

#ifndef SYNCHROTRON_H
#define SYNCHROTRON_H

class Synchrotron {
  /*
  @class: Synchrotron
  @brief: basic class for Synchrotron emission process
  */
public:
  Synchrotron(Source &s);

  // energy loss time in units [s^-1]
  void Losstime();
  void Losstime(const double mag_strength);
  void Losstime(const double mag_strength, const std::vector<double> &electron_energy,
                std::vector<double> &electron_losstime);
  void Losstime(const double mag_strength, const int particle_charge, const double particle_mass,
                const std::vector<double> &particle_energy, std::vector<double> &particle_losstime);

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
  void EmissivityProton(double mag_strength);

  // spec in units [eV^-1 s^-1]
  void Spec(const double mag_strength, const std::vector<double> &primary_energy,
            const std::vector<double> &secondary_energy);
  void Spec(const double mag_strength, const std::vector<double> &primary_energy,
            const std::vector<double> &secondary_energy, std::vector<std::vector<double>> &spec);
  void Spec(const double mag_strength, const double primary_charge, const double primary_mass,
            const std::vector<double> &primary_energy, const std::vector<double> &secondary_energy,
            std::vector<std::vector<double>> &spec);

  void SSAOpticalDepth(const double mag_strength, const double emission_size);
  void SSAAbsorptionCoeff(const double mag_strength, std::vector<double> &photon_syn_absorption_SSA);
  void Attenuation(const std::vector<double> &photon_energy, std::vector<double> &photon_spectrum);
  void Attenuation(std::vector<double> &photon_spectrum);

  std::vector<std::vector<double>> spec;

  void Test();

private:
  Source &s;
  Utility u;

  double SSACross(const double mag_strength, double x, double y);
  double mBessel53(double x);
};

#endif
