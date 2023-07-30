//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file particle.hpp
// \brief This is the particle class.

#include "utility.hpp"
#include <iostream>

#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {
  friend class EMProcess;
  friend class RadiativeProcess;
  friend class CRDistribution;
  friend class ElectronDistribution;
  friend class Photonbackground;
  friend class Synchrotron;
  friend class InverseCompton;
  friend class Bremsstrahlung;
  friend class GammaGamma;
  friend class Protonproton;
  friend class Photomeson;
  friend class Photopair;
  friend class ProtonSyn;
  friend class Photodisintegration;

public:
  Particle();
  Particle(size_t p_num, double p_min, double p_max) : p_num(p_num), p_min(p_min), p_max(p_max) {
    setMomentum();
  }
  Particle(Particle &p);
  ~Particle();

  int getID() { return ID; }
  std::string getName() { return name; }
  double getMass() { return mass; }
  double getChargeNumber() { return Z; }
  double getMassNumber() { return A; }
  std::vector<double> getMomentum() { return momentum; }
  std::vector<double> getMomentumbin() { return momentumbin; }
  std::vector<double> getMomentumghost() { return momentumghost; }
  std::vector<double> getEnergy() { return energy; }
  std::vector<double> getEnergybin() { return energybin; }
  std::vector<double> getSpectrum() { return spectrum; }
  double getSpectrum(double value) { return utility.Interpolate(momentum, spectrum, value); }

  std::vector<double> getLosstime() { return losstime; }
  std::vector<double> getIntetime() { return intetime; }
  std::vector<double> getOptdepth() { return optdepth; }
  std::vector<double> getTemp() { return temp; }

  double getSpectraInj() { return p_inj; }
  double getSpectraCut() { return p_cut; }
  double getSpectraIndex() { return p_index; }

  int getSize() { return momentum.size(); }

  void setID(int _ID);
  void setName(std::string _name);
  void setMass(double _mass);
  void setLife(double _life);
  void setCharge(double _Z);
  void setMassNumber(double _A);

  void setMomentum();
  void setMomentum(size_t _p_num, double _p_min, double _p_max);
  void setMomentum(int i, double value) { momentum[i] = value; }
  double getIndex(double value);

  void setEnergy(size_t _p_num, double _p_min, double _p_max);
  void setEnergy();

  void setSpectrum(std::vector<double> vec) { spectrum = vec; }
  void setSpectrum(int i, double value) { spectrum[i] = value; }
  void setSpectrumPL(double p_inj, double p_cut, double index, double norm, bool norm_type);
  void setSpectrumPLwExp(double _p_inj, double _p_cut, double _index, double norm,
                             bool norm_type);
  void setSpectrumPLNoExp(double p_inj, double p_cut, double index, double norm, bool norm_type);
  void setSpectrumBPL(double p_inj, double p_break, double p_cut, double index_l, double index_h,
                      double norm, bool norm_type);
  void setSpectrumBPLNoExp(double p_inj, double p_break, double p_cut, double index_l,
                           double index_h, double norm, bool norm_type);
  void setSpectrumDelta(double mu, double a, double norm);
  void setSpectrumDelta(double mu, double a, double norm, bool norm_type);
  void setSpectrumMono(double value);

  void setEnergySpectrumPL(double _p_inj, double _p_cut, double _index, double norm,
                           bool norm_type);
  void setEnergySpectrumBPL(double p_inj, double p_break, double p_cut, double index_l, double index_h,
                      double norm, bool norm_type);

  void clear();

  void setTemp(int i, double value) { temp[i] = value; }

  void setLosstime(int i, double value) { losstime[i] = value; }
  void setOptdepth(int i, double value) { optdepth[i] = value; }

  double Momentum2Energy(double momentum_x);
  double Energy2Momentum(double energy_x);
  void Momentum2Energy();
  void Energy2Momentum();

  double Interpolate(double value) { return utility.Interpolate(energy, spectrum, value); }

  double InterpolateOptdepth(double value) { return utility.Interpolate(energy, optdepth, value); }

  double Interpolate(const std::vector<double> x, const std::vector<double> y, double value) {
    return utility.Interpolate(x, y, value);
  }

private:
  Utility utility;

  int ID, Z, A;
  std::string name;
  double mass;
  double life;
  std::vector<double> momentum;      // particle momentum
  std::vector<double> momentumbin;   // particle momentum bin
  std::vector<double> momentumghost; // particle momentum ghost
  std::vector<double> energy;        // particle energy
  std::vector<double> energybin;     // particle energy
  std::vector<double> energyghost;   // particle energy
  std::vector<double> spectrum;      // particle spectrum
  std::vector<double> losstime;      // energy loss timescale [s^-1]
  std::vector<double> intetime;      // interaction timescale [s^-1]
  std::vector<double> optdepth;      // optical depth [s^-1]
  std::vector<double> temp;          //  temporary vector [In order to store temporary data]

  size_t p_num = 280;
  double p_min = 1e-7; // eV
  double p_max = 1e21; // eV

  double p_inj = 1e6;
  double p_cut = 1e16;
  double p_index = 2.2;
  double p_break = 1e12;
  double p_index_l = 2.2;
  double p_index_h = 2.2;

  void dMomentum2dEnergy();
  void dEnergy2dMomentum();
};

#endif
