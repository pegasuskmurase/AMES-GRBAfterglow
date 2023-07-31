//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file inversecompton.hpp
// \brief Basic class for inversecompton process

#include "module/source.hpp"

#ifndef INVERSECOMPTON_H
#define INVERSECOMPTON_H

class InverseCompton {
  /*
  @class: InverseCompton
  @brief: Class for Inverse Compton emission process
  */
public:
  InverseCompton(Source &s);

  // energy loss time in units [s^-1]
  void Losstime(); // Electron Inverse-Compton energy loss time
  void Losstime(const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &electron_energy, std::vector<double> &electron_losstime);

  // Interaction time in units [s^-1]
  void Intetime();
  void Intetime(const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &electron_energy, std::vector<double> &electron_intetime);

  void Spec(const std::vector<double> &primary_energy, const std::vector<double> &secondary_energy,
            const std::vector<double> &target_energy, const std::vector<double> &target_spectrum);

  // radiation emissivity in units [eV^-1 s^-1 cm^-3]
  void Emissivity(const std::vector<std::vector<double>> &spec,
                  const std::vector<double> &primary_momentum,
                  const std::vector<double> &primary_spectrum,
                  const std::vector<double> &secondary_energy,
                  std::vector<double> &secondary_spectrum);
  void Emissivity(); // Inverse-Compton emissivity
  void Emissivity(const std::vector<double> &target_energy,
                  const std::vector<double> &target_spectrum);
  void Emissivity(const std::vector<double> &target_energy,
                  const std::vector<double> &target_spectrum, std::vector<double> &photon_energy,
                  std::vector<double> &photon_spectrum);

  void Table(); // Differential cross section for secondary photons from ICS
  std::vector<std::vector<double>> spec;
  std::vector<std::vector<std::vector<double>>> table;

  void Test();

private:
  Utility u;
  Source &s;

  void Losstime(const std::vector<std::vector<double>> &table,
                const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &particle_energy, std::vector<double> &particle_losstime);

  void Spec(const std::vector<double> &primary_energy, const std::vector<double> &secondary_energy,
            const std::vector<double> &target_energy, const std::vector<double> &target_spectrum,
            const std::vector<std::vector<std::vector<double>>> &table,
            std::vector<std::vector<double>> &spec);

  void EICCrossTable();
  void ICCrossTable();
  // return Compton cross section (sigma * c)
  double ICCross(double xs) {
    // return u.Interpolate(mandelstam, ICcrossT, xs, (int)(log10(xs) / mandelstam_bin));
    return u.Interpolate(mandelstam, ICcrossT, xs);
  }
  // return Compton cross section (sigma * c)
  double EICCross(double xs) {
    // return u.Interpolate(mandelstam, EICcrossT, xs, (int)(log10(xs) / mandelstam_bin));
    return u.Interpolate(mandelstam, EICcrossT, xs);
  }
  void ReadCross();
  const static size_t mandelstam_num = 601;
  static constexpr double mandelstam_bin = 0.025;
  std::vector<double> mandelstam;
  std::vector<double> ICcrossT;
  std::vector<double> EICcrossT;
  std::vector<std::vector<double>> ICCross_table;
  std::vector<std::vector<double>> EICCross_table;

  void setAngular();
  size_t angular_num = 51;
  std::vector<double> angular;
  std::vector<double> angular_bin;
};

#endif
