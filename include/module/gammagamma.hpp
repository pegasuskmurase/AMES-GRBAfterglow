//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file gammagamma.hpp
// \brief Basic class for gammagamma process

#include "module/source.hpp"

#ifndef GAMMAGAMMA_H
#define GAMMAGAMMA_H

class GammaGamma {
  /*
  @class: GammaGamma
  @brief: Class for gamma-gamma pair annihilation process
  */
public:
  GammaGamma(Source &s);

  // energy loss time in units [s^-1]
  void Losstime(); // Electron Inverse-Compton energy loss time
  void Losstime(const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &photon_energy, std::vector<double> &photon_losstime);

  // Interaction time in units [s^-1]
  void Intetime();
  void IntetimeMuon();
  void Intetime(const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &photon_energy, std::vector<double> &photon_intetime);
  void IntetimeMuon(const std::vector<double> &target_energy,
                    const std::vector<double> &target_spectrum,
                    const std::vector<double> &photon_energy, std::vector<double> &photon_intetime);
  void Intetime(const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &photon_energy, std::vector<double> &photon_intetime,
                double redshift);

  void Attenuation(const std::vector<double> &photon_energy, std::vector<double> &photon_spectrum);
  void Attenuation(std::vector<double> &photon_spectrum);

  // return Compton cross section (sigma * c)
  double GGCross(double xs) {
    return u.Interpolate(mandelstam, GGcrossT, xs, (int)(log10(xs) / mandelstam_bin));
  }
  // return Compton cross section (sigma * c)
  double EGGCross(double xs) {
    return u.Interpolate(mandelstam, EGGcrossT, xs, (int)(log10(xs) / mandelstam_bin));
  }

  void Test();

private:
  Utility u;
  Source &s;

  void Losstime(const std::vector<std::vector<double>> &table,
                const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &particle_energy, std::vector<double> &particle_losstime);

  void EGGCrossTable();
  void GGCrossTable();
  void ReadCross();
  const static size_t mandelstam_num = 601;
  static constexpr double mandelstam_bin = 0.025;
  std::vector<double> mandelstam;
  std::vector<double> GGcrossT;
  std::vector<double> EGGcrossT;
  std::vector<std::vector<double>> GGCross_table;
  std::vector<std::vector<double>> EGGCross_table;

  void setAngular();
  size_t angular_num = 51;
  std::vector<double> angular;
  std::vector<double> angular_bin;
};

#endif
