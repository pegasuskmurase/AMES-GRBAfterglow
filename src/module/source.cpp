//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file source.hpp
// \brief This is the source class.

#include "module/source.hpp"
#include "module/constants.hpp"

Source::Source() {
  setParticle();
  setParam();
}

Source::Source(Source &s) {
  photon = s.photon;
  target = s.target;
  electron = s.electron;
  neutrino = s.neutrino;
  pionzero = s.pionzero;
  pionminus = s.pionminus;
  pionplus = s.pionplus;
  muonminusL = s.muonminusL;
  muonplusL = s.muonplusL;
  muonminusR = s.muonminusR;
  muonplusR = s.muonplusR;
  proton = s.proton;
  neutron = s.neutron;
  nucleus = s.nucleus;
  param = s.param;
  std::cout << "Calling the copy constructor " << std::endl;
}

void Source::setParticle() {
  int ID;
  std::string name;

  ID = 22;
  name = "target";
  target.setID(ID);
  target.setName(name);
  target.setCharge(0);
  target.setMassNumber(0);
  target.setMass(0.0);
  target.setLife(1e30);
  target.setMomentum(361, 1e-8, 1e10);

  ID = 22;
  name = "photon";
  photon.setID(ID);
  photon.setName(name);
  photon.setCharge(0);
  photon.setMassNumber(0);
  photon.setMass(0.0);
  photon.setLife(1e30);
  photon.setMomentum(561, 1e-8, 1e20);

  ID = 11;
  name = "electron";
  electron.setID(ID);
  electron.setName(name);
  electron.setCharge(1);
  electron.setMassNumber(0);
  electron.setMass(e_mass);
  electron.setLife(1e30);
  electron.setMomentum(381, 1e1, 1e20);

  ID = 12; // represents all flavor neutrinos
  name = "neutrino";
  neutrino.setID(ID);
  neutrino.setName(name);
  neutrino.setCharge(0);
  neutrino.setMassNumber(0);
  neutrino.setMass(0.0);
  neutrino.setLife(1e30);
  neutrino.setMomentum(151, 1e7, 1e22);

  ID = 10006;
  name = "pionzero";
  pionzero.setID(ID);
  pionzero.setName(name);
  pionzero.setCharge(0);
  pionzero.setMassNumber(0);
  pionzero.setMass(piz_Mass);
  pionzero.setLife(piz_Life);
  pionzero.setMomentum(151, 1e7, 1e22);

  ID = 10007;
  name = "pionminus";
  pionminus.setID(ID);
  pionminus.setName(name);
  pionminus.setCharge(1);
  pionminus.setMassNumber(0);
  pionminus.setMass(pic_Mass);
  pionminus.setLife(pic_Life);
  pionminus.setMomentum(151, 1e7, 1e22);

  ID = 10008;
  name = "pionplus";
  pionplus.setID(ID);
  pionplus.setName(name);
  pionplus.setCharge(1);
  pionplus.setMassNumber(0);
  pionplus.setMass(pic_Mass);
  pionplus.setLife(pic_Life);
  pionplus.setMomentum(151, 1e7, 1e22);

  ID = 100041;
  name = "muonminusL";
  muonminusL.setID(ID);
  muonminusL.setName(name);
  muonminusL.setCharge(-1);
  muonminusL.setMassNumber(0);
  muonminusL.setMass(mu_Mass);
  muonminusL.setLife(mu_Life);
  muonminusL.setMomentum(151, 1e7, 1e22);

  ID = 100042;
  name = "muonminusR";
  muonminusR.setID(ID);
  muonminusR.setName(name);
  muonminusR.setCharge(-1);
  muonminusR.setMassNumber(0);
  muonminusR.setMass(mu_Mass);
  muonminusR.setLife(mu_Life);
  muonminusR.setMomentum(151, 1e7, 1e22);

  ID = 100051;
  name = "muonplusL";
  muonplusL.setID(ID);
  muonplusL.setName(name);
  muonplusL.setCharge(1);
  muonplusL.setMassNumber(0);
  muonplusL.setMass(mu_Mass);
  muonplusL.setLife(mu_Life);
  muonplusL.setMomentum(151, 1e7, 1e22);

  ID = 100052;
  name = "muonplusR";
  muonplusR.setID(ID);
  muonplusR.setName(name);
  muonplusR.setCharge(1);
  muonplusR.setMassNumber(0);
  muonplusR.setMass(mu_Mass);
  muonplusR.setLife(mu_Life);
  muonplusR.setMomentum(151, 1e7, 1e22);

  ID = 13;
  name = "proton";
  proton.setID(ID);
  proton.setName(name);
  proton.setCharge(1);
  proton.setMassNumber(1);
  proton.setMass(proton_mass);
  proton.setLife(1e30);
  proton.setMomentum(151, 1e7, 1e22);

  ID = 14;
  name = "neutron";
  neutron.setID(ID);
  neutron.setName(name);
  neutron.setCharge(0);
  neutron.setMassNumber(1);
  neutron.setMass(neutron_mass);
  neutron.setLife(neutron_Life);
  neutron.setMomentum(151, 1e7, 1e22);

  ID = 808; // need to be set later
  name = "nucleus";
  nucleus.setID(ID);
  nucleus.setName(name);
  nucleus.setCharge(8);
  nucleus.setMassNumber(16);
  nucleus.setMass(proton_mass * 8.0 + neutron_mass * 8.0);
  nucleus.setLife(1e30);
  nucleus.setMomentum(151, 1e7, 1e22);
}

void Source::setParam() {
  param.setName("CenA");
  param.setGeometry("blob");
  param.setRedshift(0.00785);
  param.setLumDistance(2.7 * Mpc);
  param.setMagStrength(1e-4);         // G
  param.setEmissionRadius(1.0 * Kpc); // G
  param.setEmissionVolume(1.);
  param.setDopplerFactor(1.);
  param.setAmbientDensity(1.);
}

void Source::setNucleus(int mass_number, int charge_number) {
  nucleus.setMass(charge_number * proton_mass + (mass_number - charge_number) * neutron_mass);
  nucleus.setCharge(charge_number);
  nucleus.setMassNumber(mass_number);
}

void Source::showParam() {
  std::cout << " name: the name of the source "
            << " \n "
            << "geometry: can be 'blob' or 'shell'"
            << " \n "
            << "redshift: source redshift "
            << " \n "
            << "lum_distance: source luminosity distance [cm]"
            << " \n "
            << "mag_strength: the strength of magnetic field [G]"
            << " \n "
            << "emission_radiu: emission radius for spherical blob [cm]"
            << " \n "
            << "emission_volume: emission volume [cm^3]"
            << " \n "
            << "doppler_factor: source Doppler factor "
            << " \n "
            << "ambient_density: ambient matter density [cm^-3] "
            << " \n " << std::endl;
}

void Source::LorentzBoost(const std::vector<double> &photon_energy,
                          std::vector<double> &photon_spectrum) {

  int num_p = photon_energy.size();
  std::vector<double> photon_temp(num_p, 0.0);
  double D3 = param.getDopplerFactor() * param.getDopplerFactor() * param.getDopplerFactor();
  for (size_t i = 0; i < num_p; i++) {
    photon_temp[i] = photon_spectrum[i] * photon_energy[i];
  }
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] = (1 + param.getRedshift()) * D3 *
                         utility.Interpolate(photon_energy, photon_temp,
                                             (1 + param.getRedshift()) * photon_energy[i] /
                                                 param.getDopplerFactor());
  }
}

void Source::FluxNorm(const std::vector<double> &photon_energy,
                      std::vector<double> &photon_spectrum) {
  double flux_norm = (4.0 * PI / 3.0 * param.getEmissionRadius() * param.getEmissionRadius() *
                      param.getEmissionRadius()) /
                     (4 * PI * param.getLumDistance() * param.getLumDistance());
  int num_p = photon_energy.size();
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] *= flux_norm / photon_energy[i];
  }
}

void Source::LorentzBoost(const std::vector<double> &photon_energy,
                          std::vector<double> &photon_spectrum, double doppler_factor) {

  int num_p = photon_energy.size();
  std::vector<double> photon_temp(num_p, 0.0);
  double D3 = doppler_factor * doppler_factor * doppler_factor;
  for (size_t i = 0; i < num_p; i++) {
    photon_temp[i] = photon_spectrum[i] * photon_energy[i];
  }
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] =
        (1 + param.getRedshift()) * D3 *
        utility.Interpolate(photon_energy, photon_temp,
                            (1 + param.getRedshift()) * photon_energy[i] / doppler_factor);
  }
}

void Source::LuminosityNorm(const std::vector<double> &photon_energy,
                            std::vector<double> &photon_spectrum) {
  double norm = eV2erg * (4.0 * PI / 3.0 * param.getEmissionRadius() * param.getEmissionRadius() *
                          param.getEmissionRadius());
  int num_p = photon_energy.size();
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] *= norm * photon_energy[i];
  }
}

void Source::LorentzBoostShell(const std::vector<double> &photon_energy,
                               std::vector<double> &photon_spectrum, double doppler_factor) {

  int num_p = photon_energy.size();
  std::vector<double> photon_temp(num_p, 0.0);
  double D1 = doppler_factor;
  for (size_t i = 0; i < num_p; i++) {
    photon_temp[i] = photon_spectrum[i] * photon_energy[i];
  }
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] =
        (1 + param.getRedshift()) * D1 *
        utility.Interpolate(photon_energy, photon_temp,
                            (1 + param.getRedshift()) * photon_energy[i] / doppler_factor);
  }
}

void Source::LorentzBoostShell(const std::vector<double> &photon_energy,
                               std::vector<double> &photon_spectrum) {

  int num_p = photon_energy.size();
  std::vector<double> photon_temp(num_p, 0.0);
  double D1 = param.getDopplerFactor();
  for (size_t i = 0; i < num_p; i++) {
    photon_temp[i] = photon_spectrum[i] * photon_energy[i];
  }
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] = (1 + param.getRedshift()) * D1 *
                         utility.Interpolate(photon_energy, photon_temp,
                                             (1 + param.getRedshift()) * photon_energy[i] /
                                                 param.getDopplerFactor());
  }
}

void Source::FluxNormShell(const std::vector<double> &photon_energy,
                           std::vector<double> &photon_spectrum) {
  double flux_norm = eV2erg *
                     (4.0 * PI / 1.0 * param.getEmissionRadius() * param.getEmissionRadius() *
                      param.getEmissionRadius()) /
                     (4 * PI * param.getLumDistance() * param.getLumDistance());
  int num_p = photon_energy.size();
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] *= flux_norm * photon_energy[i];
  }
}

void Source::LuminosityNormShell(const std::vector<double> &photon_energy,
                                 std::vector<double> &photon_spectrum) {
  double flux_norm = eV2erg * (4.0 * PI / 1.0 * param.getEmissionRadius() *
                               param.getEmissionRadius() * param.getEmissionRadius());
  int num_p = photon_energy.size();
  for (size_t i = 0; i < num_p; i++) {
    photon_spectrum[i] *= flux_norm * photon_energy[i];
  }
}
