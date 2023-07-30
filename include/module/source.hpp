//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file source.hpp
// \brief Basic source class.

#include "particle.hpp"
#include "utility.hpp"
#include <iostream>

#ifndef SOURCE_H
#define SOURCE_H

class Param {
  /*
  @class: Param
  @brief: basic parameters of the sources
  */
public:
  void setName(std::string _name) { name = _name; }
  void setGeometry(std::string _geometry) { geometry = _geometry; }
  void setRedshift(double _redshift) { redshift = _redshift; }
  void setLumDistance(double _lum_distance) { lum_distance = _lum_distance; }
  void setMagStrength(double _mag_strength) { mag_strength = _mag_strength; }
  void setEmissionRadius(double _emission_radius) { emission_radius = _emission_radius; }
  void setVelocity(double _velocity) { velocity = _velocity; }
  void setEmissionVolume(double _emission_volume) { emission_volume = _emission_volume; }
  void setAmbientDensity(double _ambient_density) { ambient_density = _ambient_density; }
  void setDynamicalTime(double _dynamical_time) { dynamical_time = _dynamical_time; }
  void setDopplerFactor(double _doppler_factor) { doppler_factor = _doppler_factor; }
  void setLorentzFactor(double _lorentz_factor) { lorentz_factor = _lorentz_factor; }
  void setObsAngle(double _obs_angle) { obs_angle = _obs_angle; }
  std::string getName() { return name; }
  std::string getGeometry() { return geometry; }
  double getRedshift() { return redshift; }
  double getLumDistance() { return lum_distance; }
  double getMagStrength() { return mag_strength; }
  double getEmissionRadius() { return emission_radius; }
  double getEmissionVolume() { return emission_volume; }
  double getAmbientDensity() { return ambient_density; }
  double getDopplerFactor() { return doppler_factor; }
  double getVelocity() { return velocity; }
  double getLorentzFactor() { return lorentz_factor; }
  double getObsAngle() { return obs_angle; }

private:
  std::string name;
  std::string geometry;
  double redshift;
  double lum_distance;
  double mag_strength;
  double emission_radius;
  double emission_volume;
  double ambient_density;
  double dynamical_time;
  double doppler_factor;
  double velocity;
  double lorentz_factor;
  double obs_angle;
};

class Source {
  /*
  @class: Source
  @brief: The source class include many kinds of particles
  */

  friend class Photonbackground;
  friend class Synchrotron;
  friend class InverseCompton;
  friend class GammaGamma;
  friend class ElectronDistribution;
  friend class Photonbackground;

public:
  Source();
  Source(Source &s);

  void setParticle(); // set particles, default in constructor
  void setParam();    // set parameters, default in constructor
  void setNucleus(int mass_number,
                  int charge_number); // set nuclei with A and Z
  void showParam();                   // list source parameters that need to be inputed by the
                                      // user
  void LorentzBoost(const std::vector<double> &photon_energy, std::vector<double> &photon_spectrum);
  void LorentzBoost(const std::vector<double> &photon_energy, std::vector<double> &photon_spectrum,
                    double doppler_factor);
  void FluxNorm(const std::vector<double> &photon_energy, std::vector<double> &photon_spectrum);
  void LuminosityNorm(const std::vector<double> &photon_energy,
                      std::vector<double> &photon_spectrum);
  void LorentzBoostShell(const std::vector<double> &photon_energy,
                         std::vector<double> &photon_spectrum);
  void LorentzBoostShell(const std::vector<double> &photon_energy,
                         std::vector<double> &photon_spectrum, double doppler_factor);
  void FluxNormShell(const std::vector<double> &photon_energy,
                     std::vector<double> &photon_spectrum);
  void LuminosityNormShell(const std::vector<double> &photon_energy,
                           std::vector<double> &photon_spectrum);

  Particle &getPhoton() { return photon; }
  Particle &getTarget() { return target; }
  Particle &getElectron() { return electron; }
  Particle &getNeutrino() { return neutrino; }
  Particle &getPionminus() { return pionminus; }
  Particle &getPionplus() { return pionplus; }
  Particle &getPionzero() { return pionzero; }
  Particle &getMuonminusL() { return muonminusL; }
  Particle &getMuonminusR() { return muonminusR; }
  Particle &getMuonplusL() { return muonplusL; }
  Particle &getMuonplusR() { return muonplusR; }
  Particle &getProton() { return proton; }
  Particle &getNeutron() { return neutron; }
  Particle &getNucleus() { return nucleus; }
  Param &getParam() { return param; }

private:
  Particle photon;
  Particle target;
  Particle electron;
  Particle neutrino;
  Particle pionzero;
  Particle pionminus;
  Particle pionplus;
  Particle muonminusL;
  Particle muonminusR;
  Particle muonplusL;
  Particle muonplusR;
  Particle proton;
  Particle neutron;
  Particle nucleus;

  Param param;

  Utility utility;
};

#endif
