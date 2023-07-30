//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file photonbackground.hpp
// \brief photonbackground class.

#ifndef PHOTONBACKGROUND_H
#define PHOTONBACKGROUND_H

#include "source.hpp"
#include "utility.hpp"

class Photonbackground {
public:
  Photonbackground(Source &s);

  void GreyBody(double T_ph, double u_ph);
  void BlackBody(double T_ph);
  void Powerlaw(double eps_min, double eps_max, double index, double u_ph);
  void BrokenPowerlaw(double eps_min, double eps_break, double eps_max, double index_1,
                      double index_2, double u_ph);

  void CMB(double redshift);
  void EBL(double redshift, std::string EBL_model = "EBL_Gilmore12");
  void EBLTau(double redshift, std::string EBL_model = "EBL_Gilmore12");
  void EBLAttenuation(const std::vector<double> &photon_energy,
                      std::vector<double> &photon_spectrum);
  double EBLAttenuation(double photon_energy, double photon_spectrum);

  void setOutputFile(std::string _output_file) { output_file = _output_file; }

  /*
  @redshift
  @energy: [eV]
  @density: [eV cm^-3] in proper volume
  */
  void EBL_Franceschini08(double redshift);
  void EBL_Gilmore12(double redshift);
  void EBL_Finke10(double redshift);
  void EBLTau_Franceschini08(double redshift);
  void EBLTau_Gilmore12(double redshift);
  std::string EBL_Franceschini08() { return "EBL_Franceschini08"; }
  std::string EBL_Gilmore12() { return "EBL_Gilmore12"; }
  std::string EBL_Finke10() { return "EBL_Finke10"; }

private:
  Source &s;
  Utility utility;
  std::string output_file;
};

#endif
