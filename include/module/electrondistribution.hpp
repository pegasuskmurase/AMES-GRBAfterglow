//==============================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file electrondistribution.hpp
// \brief numerically determine electron distribution

#include "utility.hpp"
#include "module/gammagamma.hpp"
#include "module/photonbackground.hpp"
#include "module/compton.hpp"
#include "module/synchrotron.hpp"

#ifndef ELECTRONDISTRIBUTION_H
#define ELECTRONDISTRIBUTION_H

class Source;

class ElectronDistribution {
public:
  ElectronDistribution(Source &s);

  void IterationSolutionSteadyState(Synchrotron &syn, InverseCompton &IC, bool have_EIC = false,
                         double t_dyn = 1e10, double t_esc_photon = 1e10, double mag_strength=0.);
  void IterationSolution(Synchrotron &syn, InverseCompton &IC, bool have_EIC = false,
                         double t_dyn = 1e10, double t_ad = 1e10, double t_esc_photon = 1e10, double mag_strength=0.);

  void setMagStrength(double _mag_strength) { mag_strength = _mag_strength; }
  double getMagStrength() { return mag_strength; }

  void Test();
private:
  Source &s;
  Utility utility;
  double mag_strength;
};

#endif
