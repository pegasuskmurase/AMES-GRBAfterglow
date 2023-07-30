//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file photonbackground.hpp
// \brief photonbackground class.

#include "module/photonbackground.hpp"
#include "module/constants.hpp"
#include "module/gammagamma.hpp"

#include <cmath>
#include <fstream>

Photonbackground::Photonbackground(Source &s) : s(s) {}

void Photonbackground::GreyBody(double T_ph, double u_ph) {
  int target_size = s.target.getSize();
  double ueps;
  std::vector<double> energy = s.target.getEnergy();
  for (size_t i = 0; i < target_size; i++) {
    ueps = 15.0 * u_ph / (PI * PI * PI * PI) * pow(energy[i] / (k_B * T_ph), 4) /
           (exp(energy[i] / (k_B * T_ph)) - 1.0);
    s.target.setSpectrum(i, ueps / energy[i] / energy[i]);
  }
}

void Photonbackground::BlackBody(double T_ph) {
  int target_size = s.target.getSize();
  double CMBtemp = T_ph * 1.3806503e-16 / e_cnst;
  std::vector<double> energy = s.target.getEnergy();
  for (size_t i = 0; i < target_size; i++) {
    if ((1 - (exp(-energy[i] / T_ph))) < 1e-6) {
      s.target.setSpectrum(i, 1.318e13 * (energy[i] * CMBtemp));
    } else {
      s.target.setSpectrum(i, 1.318e13 * energy[i] * energy[i] / (exp(energy[i] / CMBtemp) - 1));
    }
  }
}

void Photonbackground::Powerlaw(double eps_min, double eps_max, double index, double u_ph) {
  double ueps;
  int target_size = s.target.getSize();
  std::vector<double> energy = s.target.getEnergy();
  for (size_t i = 0; i < target_size; i++) {
    if (energy[i] < eps_min) {
      s.target.setSpectrum(i, 0.0);
    } else {
      s.target.setSpectrum(i, energy[i] * pow(energy[i], -index) * exp(-energy[i] / eps_max));
    }
  }
  u_ph /= utility.Integrate(energy, s.target.getSpectrum());
  std::vector<double> spectrum = s.target.getSpectrum();
  for (size_t i = 0; i < target_size; i++) {
    if (energy[i] < eps_min) {
      s.target.setSpectrum(i, 0.0);
    } else {
      s.target.setSpectrum(i, u_ph * pow(energy[i], -index) * exp(-energy[i] / eps_max));
    }
  }
}

void Photonbackground::BrokenPowerlaw(double eps_min, double eps_break, double eps_max,
                                      double index_1, double index_2, double u_ph) {
  double ueps;
  int target_size = s.target.getSize();
  std::vector<double> energy = s.target.getEnergy();
  for (size_t i = 0; i < target_size; i++) {
    if (energy[i] < eps_min) {
      s.target.setSpectrum(i, 0.0);
    } else if ((energy[i] >= eps_min) && (energy[i] < eps_break)) {
      s.target.setSpectrum(i, energy[i] * pow(energy[i] / eps_break, -index_1));
    } else if ((energy[i] >= eps_break) && (energy[i] < eps_max)) {
      s.target.setSpectrum(i, energy[i] * pow(energy[i] / eps_break, -index_2));
    } else if (energy[i] >= eps_max) {
      s.target.setSpectrum(i, 0.0);
    }
  }
  u_ph /= utility.Integrate(energy, s.target.getSpectrum());
  std::vector<double> spectrum = s.target.getSpectrum();
  for (size_t i = 0; i < target_size; i++) {
    if (energy[i] < eps_min) {
      s.target.setSpectrum(i, 0.0);
    } else if ((energy[i] >= eps_min) && (energy[i] < eps_break)) {
      s.target.setSpectrum(i, u_ph * pow(energy[i] / eps_break, -index_1));
    } else if ((energy[i] >= eps_break) && (energy[i] < eps_max)) {
      s.target.setSpectrum(i, u_ph * pow(energy[i] / eps_break, -index_2));
    } else if (energy[i] >= eps_max) {
      s.target.setSpectrum(i, 0.0);
    }
  }
}

void Photonbackground::CMB(double redshift) {
  int target_size = s.target.getSize();
  double CMBtemp = 2.725 * 1.3806503e-16 / e_cnst;
  std::vector<double> energy = s.target.getEnergy();
  for (size_t i = 0; i < target_size; i++) {
    if ((1 - (exp(-energy[i] / T_CMB / (1 + redshift)))) < 1e-6) {
      s.target.setSpectrum(i, 1.318e13 * (energy[i] * CMBtemp * (1 + redshift)));
    } else {
      s.target.setSpectrum(i, 1.318e13 * energy[i] * energy[i] /
                                  (exp(energy[i] / CMBtemp / (1 + redshift)) - 1));
    }
  }
}

void Photonbackground::EBL(double redshift, std::string EBL_model) {
  int target_size = s.target.getSize();
  CMB(redshift);
  std::vector<double> cmb = s.target.getSpectrum();
  if (EBL_model == "EBL_Franceschini08") {
    EBL_Franceschini08(redshift);
  } else if (EBL_model == "EBL_Gilmore12") {
    EBL_Gilmore12(redshift);
  } else if (EBL_model == "EBL_Finke10") {
    EBL_Finke10(redshift);
  } else {
    EBL_Franceschini08(redshift);
  }
  std::vector<double> cib = s.target.getSpectrum();
  for (size_t i = 0; i < target_size; i++) {
    s.target.setSpectrum(i, cmb[i] + cib[i]);
  }
}

void Photonbackground::EBLTau(double redshift, std::string EBL_model) {
  if (EBL_model == "EBL_Franceschini08") {
    EBLTau_Franceschini08(redshift);
  } else if (EBL_model == "EBL_Gilmore12") {
    EBLTau_Gilmore12(redshift);
  } else {
    std::cout << "No such table!" << std::endl;
  }
}

void Photonbackground::EBLAttenuation(const std::vector<double> &photon_energy,
                                      std::vector<double> &photon_spectrum) {
  double tau;
  int num_p = photon_spectrum.size();
  for (size_t i = 0; i < num_p; i++) {
    tau = utility.Interpolate(s.photon.energy, s.photon.optdepth, photon_energy[i]);
    photon_spectrum[i] *= exp(-tau);
  }
}

double Photonbackground::EBLAttenuation(double photon_energy, double photon_spectrum) {
  double tau;
  tau = utility.Interpolate(s.photon.energy, s.photon.optdepth, photon_energy);
  photon_spectrum *= exp(-tau);
  return photon_spectrum;
}

void Photonbackground::EBL_Franceschini08(double redshift) {

  std::vector<double> raw_redshift;
  std::vector<double> raw_energy;
  std::vector<std::vector<double>> raw_density; // 2d
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Franceschini_2008/target_redshift.dat",
                     raw_redshift);
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Franceschini_2008/target_energies.dat",
                     raw_energy);
  utility.ReadFile2D(std::string(DATA_PATH) + "/data/EBL/Franceschini_2008/target_density.dat",
                     raw_density); // proper density

  std::vector<double> target = s.target.getEnergy();
  std::vector<double> raw_energy_z(raw_energy.size(), 0.0);
  for (size_t i = 0; i < raw_energy_z.size(); i++) {
    raw_energy_z[i] = raw_energy[i] * (1 + redshift);

    const int target_size = target.size();
    for (size_t i = 0; i < target_size; i++) {
      s.target.setSpectrum(
          i, utility.Interpolate2D(raw_redshift, raw_energy_z, raw_density, redshift, target[i]));
    }
  }
}

void Photonbackground::EBLTau_Franceschini08(double redshift) {
  std::vector<double> opdep_redshift;
  std::vector<double> opdep_energy;
  std::vector<std::vector<double>> raw_tau; // 2d
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Franceschini_2008/opdep_redshift.dat",
                     opdep_redshift);
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Franceschini_2008/opdep_energies.dat",
                     opdep_energy);
  utility.ReadFile2D(std::string(DATA_PATH) + "/data/EBL/Franceschini_2008/opdep.dat", raw_tau);

  std::vector<double> photon_energy = s.photon.getEnergy();
  const int num_p = photon_energy.size();
  for (size_t i = 0; i < num_p; i++) {
    s.photon.setOptdepth(i, utility.Interpolate2D(opdep_redshift, opdep_energy, raw_tau, redshift,
                                                  photon_energy[i]));
  }
}

void Photonbackground::EBL_Gilmore12(double redshift) {

  std::vector<double> raw_redshift;
  std::vector<double> raw_energy;
  std::vector<std::vector<double>> raw_density; // 2d
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/target_redshift.dat",
                     raw_redshift);
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/target_energies.dat",
                     raw_energy);
  utility.ReadFile2D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/target_density.dat",
                     raw_density); // proper density

  std::vector<double> target = s.target.getEnergy();
  const int target_size = target.size();
  for (size_t i = 0; i < target_size; i++) {
    s.target.setSpectrum(
        i, utility.Interpolate2D(raw_redshift, raw_energy, raw_density, redshift, target[i]));
  }
}

void Photonbackground::EBLTau_Gilmore12(double redshift) {
  std::vector<double> opdep_redshift;
  std::vector<double> opdep_energy;
  std::vector<std::vector<double>> raw_tau; // 2d
  /*
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/opdep_redshift.dat",
                     opdep_redshift);
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/opdep_energies.dat",
                     opdep_energy);
  utility.ReadFile2D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/opdep.dat", raw_tau);
  */
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/opdep_redshift_wCMB.dat",
                     opdep_redshift);
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/opdep_energies_wCMB.dat",
                     opdep_energy);
  utility.ReadFile2D(std::string(DATA_PATH) + "/data/EBL/Gilmore_2012/opdep_wCMB.dat", raw_tau);

  std::vector<double> photon_energy = s.photon.getEnergy();
  const int num_p = photon_energy.size();
  for (size_t i = 0; i < num_p; i++) {
    s.photon.setOptdepth(i, utility.Interpolate2D(opdep_redshift, opdep_energy, raw_tau, redshift,
                                                  photon_energy[i]));
  }
}

void Photonbackground::EBL_Finke10(double redshift) {
  std::vector<double> raw_redshift;
  std::vector<double> raw_energy;
  std::vector<std::vector<double>> raw_density; // 2d
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Finke_2010/target_redshift.dat",
                     raw_redshift);
  utility.ReadFile1D(std::string(DATA_PATH) + "/data/EBL/Finke_2010/target_energies.dat",
                     raw_energy);
  utility.ReadFile2D(std::string(DATA_PATH) + "/data/EBL/Finke_2010/target_density.dat",
                     raw_density);

  std::vector<double> target = s.target.getEnergy();
  const int target_size = target.size();
  for (size_t i = 0; i < target_size; i++) {
    s.target.setSpectrum(
        i, utility.Interpolate2D(raw_redshift, raw_energy, raw_density, redshift, target[i]));
  }
}
