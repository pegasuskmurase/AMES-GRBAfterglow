//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file particle.cpp
// \brief particle class.

#include "module/particle.hpp"
#include "module/constants.hpp"
#include "module/utility.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

Particle::Particle() {}

Particle::~Particle() {}

Particle::Particle(Particle &p) {
  energy = p.energy;
  spectrum = p.spectrum;
  losstime = p.losstime;
  optdepth = p.optdepth;
  temp = p.temp;
}

void Particle::setID(int _ID) { ID = _ID; }

void Particle::setName(std::string _name) { name = _name; }

void Particle::setMass(double _mass) { mass = _mass; }
void Particle::setLife(double _life) { life = _life; }

void Particle::setCharge(double _Z) { Z = _Z; }

void Particle::setMassNumber(double _A) { A = _A; }

double Particle::Momentum2Energy(double momentum_x) {
  return sqrt(momentum_x * momentum_x + mass * mass);
}

double Particle::Energy2Momentum(double energy_x) {
  return sqrt(energy_x * energy_x - mass * mass);
}

void Particle::Momentum2Energy() {
  for (size_t i = 0; i < p_num; i++) {
    energy[i] = sqrt(momentum[i] * momentum[i] + mass * mass);
  }
}

void Particle::Energy2Momentum() {
  for (size_t i = 0; i < p_num; i++) {
    if (energy[i] >= mass) {
      momentum[i] = sqrt(energy[i] * energy[i] - mass * mass);
    } else {
      momentum[i] = 0.0;
    }
  }
}

void Particle::dMomentum2dEnergy() {
  double factor;
  for (size_t i = 0; i < p_num; i++) {
    factor = momentum[i] / sqrt(momentum[i] * momentum[i] + mass * mass);
    spectrum[i] /= factor;
  }
}

void Particle::dEnergy2dMomentum() {
  double factor;
  for (size_t i = 0; i < p_num; i++) {
    factor = momentum[i] / sqrt(momentum[i] * momentum[i] + mass * mass);
    spectrum[i] *= factor;
  }
}

void Particle::setMomentum() {
  momentum.resize(p_num, 0.);
  momentumbin.resize(p_num, 0.);
  momentumghost.resize(3, 0.);
  energy.resize(p_num, 0.);
  energybin.resize(p_num, 0.);
  energyghost.resize(3, 0.);
  spectrum.resize(p_num, 0.);
  losstime.resize(p_num, 0.);
  intetime.resize(p_num, 0.);
  optdepth.resize(p_num, 0.);
  temp.resize(p_num, 0.);
  for (size_t i = 0; i < p_num; i++) {
    momentum[i] = p_min * pow(p_max / p_min, i / (p_num - 1.));
  }

  momentumghost[0] = p_min * pow(p_max / p_min, -1 / (p_num - 1.));           //-1
  momentumghost[1] = p_min * pow(p_max / p_min, p_num / (p_num - 1.));        // N
  momentumghost[2] = p_min * pow(p_max / p_min, (p_num + 1.) / (p_num - 1.)); // index, N + 1

  energyghost[0] = Momentum2Energy(momentumghost[0]);
  energyghost[1] = Momentum2Energy(momentumghost[1]);
  energyghost[2] = Momentum2Energy(momentumghost[2]);

  double momentum_l = p_min * pow(p_max / p_min, -1 / (p_num - 1.));
  double momentum_r = p_min * pow(p_max / p_min, p_num / (p_num - 1.));
  for (size_t i = 0; i < p_num; i++) {
    if (i == 0) {
      momentumbin[i] = 0.5 * (momentum[i + 1] - momentum_l);
    } else if (i == p_num - 1) {
      momentumbin[i] = 0.5 * (momentum_r - momentum[i - 1]);
    } else {
      momentumbin[i] = 0.5 * (momentum[i + 1] - momentum[i - 1]);
    }
  }
  Momentum2Energy();

  double energy_l = Momentum2Energy(momentum_l);
  double energy_r = Momentum2Energy(momentum_r);
  for (size_t i = 0; i < p_num; i++) {
    if (i == 0) {
      energybin[i] = 0.5 * (energy[i + 1] - energy_l);
    } else if (i == p_num - 1) {
      energybin[i] = 0.5 * (energy_r - energy[i - 1]);
    } else {
      energybin[i] = 0.5 * (energy[i + 1] - energy[i - 1]);
    }
  }
}

void Particle::setEnergy() {
  momentum.resize(p_num, 0.);
  momentumbin.resize(p_num, 0.);
  momentumghost.resize(3, 0.);
  energy.resize(p_num, 0.);
  energybin.resize(p_num, 0.);
  energyghost.resize(3, 0.);
  spectrum.resize(p_num, 0.);
  losstime.resize(p_num, 0.);
  intetime.resize(p_num, 0.);
  optdepth.resize(p_num, 0.);
  temp.resize(p_num, 0.);
  for (size_t i = 0; i < p_num; i++) {
    energy[i] = p_min * pow(p_max / p_min, i / (p_num - 1.));
  }

  energyghost[0] = p_min * pow(p_max / p_min, -1 / (p_num - 1.));           //-1
  energyghost[1] = p_min * pow(p_max / p_min, p_num / (p_num - 1.));        // N
  energyghost[2] = p_min * pow(p_max / p_min, (p_num + 1.) / (p_num - 1.)); // index, N + 1

  double energy_l = p_min * pow(p_max / p_min, -1 / (p_num - 1.));
  double energy_r = p_min * pow(p_max / p_min, p_num / (p_num - 1.));
  for (size_t i = 0; i < p_num; i++) {
    if (i == 0) {
      energybin[i] = 0.5 * (energy[i + 1] - energy_l);
    } else if (i == p_num - 1) {
      energybin[i] = 0.5 * (energy_r - energy[i - 1]);
    } else {
      energybin[i] = 0.5 * (energy[i + 1] - energy[i - 1]);
    }
  }
  Energy2Momentum();

  momentumghost[0] = Energy2Momentum(energyghost[0]);
  momentumghost[1] = Energy2Momentum(energyghost[1]);
  momentumghost[2] = Energy2Momentum(energyghost[2]);

  double momentum_l = Momentum2Energy(energy_l);
  double momentum_r = Momentum2Energy(energy_r);
  for (size_t i = 0; i < p_num; i++) {
    if (i == 0) {
      momentumbin[i] = 0.5 * (momentum[i + 1] - momentum_l);
    } else if (i == p_num - 1) {
      momentumbin[i] = 0.5 * (momentum_r - momentum[i - 1]);
    } else {
      momentumbin[i] = 0.5 * (momentum[i + 1] - momentum[i - 1]);
    }
  }
}

void Particle::setMomentum(size_t _p_num, double _p_min, double _p_max) {
  p_num = _p_num;
  p_min = _p_min;
  p_max = _p_max;
  setMomentum();
}

void Particle::setEnergy(size_t _p_num, double _p_min, double _p_max) {
  p_num = _p_num;
  p_min = _p_min;
  p_max = _p_max;
  setEnergy();
}

double Particle::getIndex(double value) {
  // value is energy, convert to momentum
  int idx = (p_num - 1) * log(Energy2Momentum(value) / p_min) / log(p_max / p_min);
  if (value <= p_min)
    idx = 0;
  if (value >= p_max)
    idx = p_num - 1;
  return idx;
}

void Particle::setSpectrumPLNoExp(double _p_inj, double _p_cut, double _index, double norm,
                                  bool norm_type) {
  // power-law distribution without exponential cutoff
  p_inj = _p_inj;
  p_cut = _p_cut;
  p_index = _index;
  for (size_t i = 0; i < p_num; i++) {
    if (momentum[i] < p_inj) {
      spectrum[i] = 0.;
    } else if (momentum[i] > p_cut) {
      spectrum[i] = 0.;
    } else {
      if (norm_type == true) {
        spectrum[i] = energy[i] * exp((-p_index) * log(momentum[i]));
      } else {
        spectrum[i] = exp((-p_index) * log(momentum[i]));
      }
    }
  }

  norm /= utility.Integrate(momentum, spectrum);

  for (size_t i = 0; i < p_num; i++) {
    if (momentum[i] < p_inj) {
      spectrum[i] = 0.;
    } else if (momentum[i] > p_cut) {
      spectrum[i] = 0.;
    } else {
      spectrum[i] = norm * exp((-p_index) * log(momentum[i]));
    }
  }
}

void Particle::setSpectrumPLwExp(double _p_inj, double _p_cut, double _index, double norm,
                                 bool norm_type) {
  // power-law distribution with exponential cutoff
  p_inj = _p_inj;
  p_cut = _p_cut;
  p_index = _index;
  double dum;
  for (size_t i = 0; i < p_num; i++) {
    if (norm_type == true) {
      spectrum[i] = energy[i] * exp(-p_index * log(momentum[i])) * exp(-p_inj / momentum[i]) *
                    exp(-momentum[i] / p_cut); // total energy
    } else {
      spectrum[i] =
          exp(-p_index * log(momentum[i])) * exp(-p_inj / momentum[i]) * exp(-momentum[i] / p_cut);
    }
  }

  norm /= (utility.Integrate(momentum, spectrum));

  for (size_t i = 0; i < p_num; i++) {
    spectrum[i] = norm * exp(-p_index * log(momentum[i])) * exp(-p_inj / momentum[i]) *
                  exp(-momentum[i] / p_cut);
  }
}

void Particle::setSpectrumPL(double _p_inj, double _p_cut, double _index, double norm,
                             bool norm_type) {
  // power-law distribution with exponential cutoff
  p_inj = _p_inj;
  p_cut = _p_cut;
  p_index = _index;
  int idx = 0;
  double dum;
  for (size_t i = 0; i < p_num; i++) {
    if (momentum[i] < p_inj) {
      idx++;
      spectrum[i] = 0.;
    } else {
      if (norm_type == true) {
        spectrum[i] = (energy[i] - mass) * exp(-p_index * log(momentum[i])) *
                      exp(-(momentum[i] / p_cut)); // total energy
      } else {
        spectrum[i] = exp(-p_index * log(momentum[i])) * exp(-momentum[i] / p_cut);
      }
    }
  }

  norm /= (utility.Integrate(momentum, spectrum));

  for (size_t i = 0; i < p_num; i++) {
    if (momentum[i] < p_inj) {
      spectrum[i] = 0.;
    } else {
      spectrum[i] = norm * exp(-p_index * log(momentum[i])) * exp(-(momentum[i] / p_cut));
    }
  }
}

void Particle::setEnergySpectrumPL(double _p_inj, double _p_cut, double _index, double norm,
                                   bool norm_type) {
  // power-law distribution with exponential cutoff
  p_inj = _p_inj;
  p_cut = _p_cut;
  p_index = _index;
  int idx = 0;
  double dum;
  for (size_t i = 0; i < p_num; i++) {
    if (energy[i] < p_inj) {
      idx++;
      spectrum[i] = 0.;
    } else {
      if (norm_type == true) {
        spectrum[i] =
            energy[i] * exp((-p_index) * log(energy[i])) * exp(-energy[i] / p_cut); // total energy
      } else {
        spectrum[i] = exp((-p_index) * log(energy[i])) * exp(-energy[i] / p_cut);
      }
    }
  }

  norm /= (utility.Integrate(energy, spectrum));

  for (size_t i = 0; i < p_num; i++) {
    if (energy[i] < p_inj) {
      spectrum[i] = 0.;
    } else {
      spectrum[i] = norm * exp(-p_index * log(energy[i])) * exp(-energy[i] / p_cut);
    }
  }
}

void Particle::setSpectrumBPLNoExp(double _p_inj, double _p_break, double _p_cut, double _index_l,
                                   double _index_h, double norm, bool norm_type) {
  // broken power-law distribution with exponential cutoff
  p_inj = _p_inj;
  p_break = _p_break;
  p_cut = _p_cut;
  p_index_l = _index_l;
  p_index_h = _index_h;
  for (size_t i = 0; i < p_num; i++) {
    if (momentum[i] < p_inj) {
      spectrum[i] = 0.;
    } else if ((momentum[i] >= p_inj) && (momentum[i] < p_break)) {
      if (norm_type == true) {
        spectrum[i] = energy[i] * exp(-p_index_l * log(momentum[i] / p_break));
      } else {
        spectrum[i] = exp(-p_index_l * log(momentum[i] / p_break));
      }
    } else if ((momentum[i] >= p_break) && (momentum[i] < p_cut)) {
      if (norm_type == true) {
        spectrum[i] = energy[i] * exp(-p_index_h * log(momentum[i] / p_break));
      } else {
        spectrum[i] = exp(-p_index_h * log(momentum[i] / p_break));
      }
    } else {
      spectrum[i] = 0.;
    }
  }

  norm /= Utility().Integrate(momentum, spectrum);

  for (size_t i = 0; i < p_num; i++) {
    if (momentum[i] < p_inj) {
      spectrum[i] = 0.;
    } else if ((momentum[i] >= p_inj) && (momentum[i] < p_break)) {
      spectrum[i] = norm * exp(-p_index_l * log(momentum[i] / p_break));
    } else if ((momentum[i] >= p_break) && (momentum[i] < p_cut)) {
      spectrum[i] = norm * exp(-p_index_h * log(momentum[i] / p_break));
    } else {
      spectrum[i] = 0.;
    }
  }
}

void Particle::setSpectrumBPL(double _p_inj, double _p_break, double _p_cut, double _index_l,
                              double _index_h, double norm, bool norm_type) {
  // broken power-law distribution with exponential cutoff
  p_inj = _p_inj;
  p_break = _p_break;
  p_cut = _p_cut;
  p_index_l = _index_l;
  p_index_h = _index_h;
  for (size_t i = 0; i < p_num; i++) {
    if (momentum[i] < p_inj) {
      spectrum[i] = 0.;
    } else if ((momentum[i] >= p_inj) && (momentum[i] < p_break)) {
      if (norm_type == true) {
        spectrum[i] =
            energy[i] * exp(-p_index_l * log(momentum[i] / p_break)) * exp(-momentum[i] / p_cut);
      } else {
        spectrum[i] = exp(-p_index_l * log(momentum[i] / p_break)) * exp(-momentum[i] / p_cut);
      }
    } else {
      if (norm_type == true) {
        spectrum[i] =
            energy[i] * exp(-p_index_h * log(momentum[i] / p_break)) * exp(-momentum[i] / p_cut);
      } else {
        spectrum[i] = exp(-p_index_h * log(momentum[i] / p_break)) * exp(-momentum[i] / p_cut);
      }
    }
  }

  norm /= Utility().Integrate(momentum, spectrum);

  for (size_t i = 0; i < p_num; i++) {
    if (momentum[i] < p_inj) {
      spectrum[i] = 0.;
    } else if ((momentum[i] >= p_inj) && (momentum[i] < p_break)) {
      spectrum[i] = norm * exp(-p_index_l * log(momentum[i] / p_break)) * exp(-momentum[i] / p_cut);
    } else {
      spectrum[i] = norm * exp(-p_index_h * log(momentum[i] / p_break)) * exp(-momentum[i] / p_cut);
    }
  }
}

void Particle::setEnergySpectrumBPL(double _p_inj, double _p_break, double _p_cut, double _index_l,
                                    double _index_h, double norm, bool norm_type) {
  // broken power-law distribution with exponential cutoff
  p_inj = _p_inj;
  p_break = _p_break;
  p_cut = _p_cut;
  p_index_l = _index_l;
  p_index_h = _index_h;
  for (size_t i = 0; i < p_num; i++) {
    if (energy[i] < p_inj) {
      spectrum[i] = 0.;
    } else if ((energy[i] >= p_inj) && (energy[i] < p_break)) {
      if (norm_type == true) {
        spectrum[i] =
            energy[i] * exp(-p_index_l * log(energy[i] / p_break)) * exp(-energy[i] / p_cut);
      } else {
        spectrum[i] = exp(-p_index_l * log(energy[i] / p_break)) * exp(-energy[i] / p_cut);
      }
    } else {
      if (norm_type == true) {
        spectrum[i] =
            energy[i] * exp(-p_index_h * log(energy[i] / p_break)) * exp(-energy[i] / p_cut);
      } else {
        spectrum[i] = exp(-p_index_h * log(energy[i] / p_break)) * exp(-energy[i] / p_cut);
      }
    }
  }

  norm /= Utility().Integrate(energy, spectrum);

  for (size_t i = 0; i < p_num; i++) {
    if (energy[i] < p_inj) {
      spectrum[i] = 0.;
    } else if ((energy[i] >= p_inj) && (energy[i] < p_break)) {
      spectrum[i] = norm * exp(-p_index_l * log(energy[i] / p_break)) * exp(-energy[i] / p_cut);
    } else {
      spectrum[i] = norm * exp(-p_index_h * log(energy[i] / p_break)) * exp(-energy[i] / p_cut);
    }
  }
}

void Particle::setSpectrumDelta(double mu, double a, double norm) {
  for (size_t i = 0; i < p_num; i++) {
    spectrum[i] = norm / a / PI * exp(-(momentum[i] - mu) * (momentum[i] - mu) / a / a);
  }
}

void Particle::setSpectrumDelta(double mu, double a, double norm, bool norm_type) {
    std::cout << norm << " " << a << " " << mu << std::endl;
  for (size_t i = 0; i < p_num; i++) {
    if (norm_type) {
      spectrum[i] = energy[i] / a / PI * exp(-(momentum[i] - mu) * (momentum[i] - mu) / a / a);
    } else {
      spectrum[i] = 1. / a / PI * exp(-(momentum[i] - mu) * (momentum[i] - mu) / a / a);
    }
  }

  norm /= utility.Integrate(momentum, spectrum);

  for (size_t i = 0; i < p_num; i++) {
    spectrum[i] = norm / a / PI * exp(-(momentum[i] - mu) * (momentum[i] - mu) / a / a);
  }
}

void Particle::setSpectrumMono(double value) {
  p_num = 1;
  setMomentum();
  energy[0] = value;
  spectrum[0] = 1;
}

void Particle::clear() {
  for (double &i : spectrum) {
    i = 0.0;
  }
}
