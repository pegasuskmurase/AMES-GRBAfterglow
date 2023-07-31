//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file gammagamma.hpp
// \brief Basic class for gamma-gamma pair annihilation process

#include "module/gammagamma.hpp"
#include "module/constants.hpp"
#include "module/utility.hpp"
#include <fstream>

GammaGamma::GammaGamma(Source &s) : s(s) {
  setAngular();
  ReadCross();
}

void GammaGamma::Losstime(const std::vector<std::vector<double>> &table,
                          const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum,
                          const std::vector<double> &particle_energy,
                          std::vector<double> &particle_losstime) {

  int num_t = target_energy.size();
  std::vector<double> target_temp(num_t);
  int num_p = particle_energy.size();
  particle_losstime.resize(num_p);
  for (size_t i = 0; i < num_p; i++) {
    for (size_t k = 0; k < num_t; k++) {
      target_temp[k] = target_spectrum[k] * table[i][k];
    }
    particle_losstime[i] = u.Integrate(target_energy, target_temp);
  }
  std::vector<double>().swap(target_temp);
}

void GammaGamma::Losstime() {
  if (EGGCross_table.size() <= 0) {
    std::cout << "table size = " << EGGCross_table.size() << std::endl;
    std::cout << "Generate GG table " << std::endl;
    EGGCrossTable();
  }

  if ((EGGCross_table.size() != s.photon.momentum.size()) ||
      (EGGCross_table[0].size() != s.target.momentum.size())) {
    std::cout << "table size : photon " << EGGCross_table.size() << std::endl;
    std::cout << "table size : target " << EGGCross_table[0].size() << std::endl;
    EGGCrossTable();
    std::cout << "table size : photon " << EGGCross_table.size() << std::endl;
    std::cout << "table size : target  " << EGGCross_table[0].size() << std::endl;
  }

  Losstime(EGGCross_table, s.target.energy, s.target.spectrum, s.photon.energy, s.photon.losstime);
}

void GammaGamma::Losstime(const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum,
                          const std::vector<double> &photon_energy,
                          std::vector<double> &photon_losstime) {
  if (EGGCross_table.size() <= 0) {
    std::cout << "table size = " << EGGCross_table.size() << std::endl;
    std::cout << "Generate GG table " << std::endl;
    EGGCrossTable();
  }

  if ((EGGCross_table.size() != s.photon.momentum.size()) ||
      (EGGCross_table[0].size() != s.target.momentum.size())) {
    std::cout << "table size : photon " << EGGCross_table.size() << std::endl;
    std::cout << "table size : target " << EGGCross_table[0].size() << std::endl;
    EGGCrossTable();
    std::cout << "table size : photon " << EGGCross_table.size() << std::endl;
    std::cout << "table size : target  " << EGGCross_table[0].size() << std::endl;
  }

  Losstime(EGGCross_table, target_energy, target_spectrum, photon_energy, photon_losstime);
}

void GammaGamma::Intetime() {
  if (GGCross_table.size() <= 0) {
    std::cout << "table size = " << GGCross_table.size() << std::endl;
    std::cout << "Generate GG table " << std::endl;
    GGCrossTable();
  }

  if ((GGCross_table.size() != s.photon.momentum.size()) ||
      (GGCross_table[0].size() != s.target.momentum.size())) {
    std::cout << "table size : photon " << GGCross_table.size() << std::endl;
    std::cout << "table size : target " << GGCross_table[0].size() << std::endl;
    GGCrossTable();
    std::cout << "table size : photon " << GGCross_table.size() << std::endl;
    std::cout << "table size : target  " << GGCross_table[0].size() << std::endl;
  }

  Losstime(GGCross_table, s.target.energy, s.target.spectrum, s.photon.energy, s.photon.intetime);
}

void GammaGamma::Intetime(const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum,
                          const std::vector<double> &photon_energy,
                          std::vector<double> &photon_intetime) {
  if (GGCross_table.size() <= 0) {
    std::cout << "table size = " << GGCross_table.size() << std::endl;
    std::cout << "Generate GG table " << std::endl;
    GGCrossTable();
  }

  if ((GGCross_table.size() != s.photon.momentum.size()) ||
      (GGCross_table[0].size() != s.target.momentum.size())) {
    std::cout << "table size : photon " << GGCross_table.size() << std::endl;
    std::cout << "table size : target " << GGCross_table[0].size() << std::endl;
    GGCrossTable();
    std::cout << "table size : photon " << GGCross_table.size() << std::endl;
    std::cout << "table size : target  " << GGCross_table[0].size() << std::endl;
  }

  Losstime(GGCross_table, target_energy, target_spectrum, photon_energy, photon_intetime);
}

void GammaGamma::Intetime(const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum,
                          const std::vector<double> &photon_energy,
                          std::vector<double> &photon_intetime, double redshift) {

  if (GGCross_table.size() <= 0) {
    std::cout << "table size = " << GGCross_table.size() << std::endl;
    std::cout << "Generate GG table " << std::endl;
    GGCrossTable();
  }

  if ((GGCross_table.size() != s.photon.momentum.size()) ||
      (GGCross_table[0].size() != s.target.momentum.size())) {
    std::cout << "table size : photon " << GGCross_table.size() << std::endl;
    std::cout << "table size : target " << GGCross_table[0].size() << std::endl;
    GGCrossTable();
    std::cout << "table size : photon " << GGCross_table.size() << std::endl;
    std::cout << "table size : target  " << GGCross_table[0].size() << std::endl;
  }

  int num_t = target_energy.size();
  int num_p = photon_energy.size();
  std::vector<double> target_temp(num_t);
  photon_intetime.resize(num_p);
  for (size_t i = 0; i < num_p; i++) {
    for (size_t k = 0; k < num_t; k++) {
      target_temp[k] =
          target_spectrum[k] * u.Interpolate2D(photon_energy, target_energy, GGCross_table,
                                               photon_energy[i] * (1 + redshift), target_energy[k]);
    }
    photon_intetime[i] = u.Integrate(target_energy, target_temp);
  }
  std::vector<double>().swap(target_temp);
}

void GammaGamma::Attenuation(std::vector<double> &photon_spectrum) {
  double tau, att;
  int num_p = photon_spectrum.size();
  for (size_t i = 0; i < num_p; i++) {
    tau = s.photon.optdepth[i];
    if (tau > 1e-5) {
      att = (1. - exp(-tau)) / tau;
    } else {
      att = 1. / (1 + tau);
    }
    photon_spectrum[i] *= att;
  }
}

void GammaGamma::Attenuation(const std::vector<double> &photon_energy,
                             std::vector<double> &photon_spectrum) {
  double tau, att;
  int num_p = photon_spectrum.size();
  for (size_t i = 0; i < num_p; i++) {
    tau = u.Interpolate(s.photon.energy, s.photon.optdepth, photon_energy[i]);
    if (tau > 1e-5) {
      att = (1. - exp(-tau)) / tau;
    } else {
      att = 1. / (1 + tau);
    }
    photon_spectrum[i] *= att;
  }
}

void GammaGamma::EGGCrossTable() {
  int num_p = s.photon.momentum.size();
  int num_t = s.target.momentum.size();
  EGGCross_table.resize(num_p);
  for (int i = 0; i < num_p; ++i) {
    EGGCross_table[i].resize(num_t);
  }
  double results, dum, mandels1, mandels2;
  for (size_t k = 0; k < num_p; k++) {
    for (size_t i = 0; i < num_t; i++) {
      dum = 0.;
      for (size_t j = 0; j < angular_num - 1; j++) {
        mandels1 =
            2. * s.photon.energy[k] * s.target.energy[i] * (1 - angular[j]) / (e_mass * e_mass);
        mandels2 =
            2. * s.photon.energy[k] * s.target.energy[i] * (1 - angular[j + 1]) / (e_mass * e_mass);
        dum += 0.5 * angular_bin[j] *
               ((1 - angular[j]) * EGGCross(mandels1) + (1 - angular[j + 1]) * EGGCross(mandels2));
      }
      dum *= 0.5;
      EGGCross_table[k][i] = dum;
    }
  }
}

void GammaGamma::GGCrossTable() {
  int num_p = s.photon.momentum.size();
  int num_t = s.target.momentum.size();
  GGCross_table.resize(num_p);
  for (int i = 0; i < num_p; ++i) {
    GGCross_table[i].resize(num_t);
  }
  double results, dum, mandels1, mandels2;
  for (size_t k = 0; k < num_p; k++) {
    for (size_t i = 0; i < num_t; i++) {
      dum = 0.;
      for (size_t j = 0; j < angular_num - 1; j++) {
        mandels1 =
            2. * s.photon.energy[k] * s.target.energy[i] * (1 - angular[j]) / (e_mass * e_mass);
        mandels2 =
            2. * s.photon.energy[k] * s.target.energy[i] * (1 - angular[j + 1]) / (e_mass * e_mass);
        dum += 0.5 * angular_bin[j] *
               ((1 - angular[j]) * GGCross(mandels1) + (1 - angular[j + 1]) * GGCross(mandels2));
      }
      dum *= 0.5;
      GGCross_table[k][i] = dum;
    }
  }
}

void GammaGamma::ReadCross() {
  mandelstam.resize(mandelstam_num);
  GGcrossT.resize(mandelstam_num);
  EGGcrossT.resize(mandelstam_num);
  double dum;
  std::ifstream gin1(std::string(DATA_PATH) + std::string("/data/EMProcess/mandelstam_table.dat"));
  if (gin1.good()) {
    for (int i = 0; i < mandelstam_num; i++) {
      gin1 >> mandelstam[i] >> GGcrossT[i] >> EGGcrossT[i] >> dum >> dum;
    }
  } else {
    std::cout << "PPset files does't exists!!!" << std::endl;
    exit(1);
  }
  gin1.close();
}

void GammaGamma::setAngular() {
  angular.resize(angular_num);
  angular_bin.resize(angular_num);
  for (size_t i = 0; i < angular_num; i++) {
    angular[i] = -1. + 2. * i / (angular_num - 1.);
    angular_bin[i] = 2. / (angular_num - 1);
  }
}

void GammaGamma::Test() { std::cout << "Test : GammaGamma module         " << std::endl; }
