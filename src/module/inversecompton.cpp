//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file inversecompton.hpp
// \brief Basic class for inversecompton process

#include "module/inversecompton.hpp"
#include "module/constants.hpp"
#include "module/utility.hpp"
#include <fstream>

InverseCompton::InverseCompton(Source &s) : s(s) {
  setAngular();
  ReadCross();
  ICCrossTable();
  EICCrossTable();
}

void InverseCompton::Losstime(const std::vector<std::vector<double>> &table,
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

void InverseCompton::Losstime() {
  Losstime(EICCross_table, s.target.energy, s.target.spectrum, s.electron.energy,
           s.electron.losstime);
}

void InverseCompton::Losstime(const std::vector<double> &target_energy,
                              const std::vector<double> &target_spectrum,
                              const std::vector<double> &electron_energy,
                              std::vector<double> &electron_losstime) {
  Losstime(EICCross_table, target_energy, target_spectrum, electron_energy, electron_losstime);
}

void InverseCompton::Intetime() {
  Losstime(ICCross_table, s.target.energy, s.target.spectrum, s.electron.energy,
           s.electron.intetime);
}

void InverseCompton::Intetime(const std::vector<double> &target_energy,
                              const std::vector<double> &target_spectrum,
                              const std::vector<double> &electron_energy,
                              std::vector<double> &electron_intetime) {
  Losstime(ICCross_table, target_energy, target_spectrum, electron_energy, electron_intetime);
}

void InverseCompton::Spec(const std::vector<double> &primary_energy,
                          const std::vector<double> &secondary_energy,
                          const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum,
                          const std::vector<std::vector<std::vector<double>>> &table,
                          std::vector<std::vector<double>> &spec) {
  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  int num_t = target_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].assign(num_p, 0.0);
  }
  std::vector<double> target_temp(num_t, 0.0);
  double dum;
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      for (size_t k = 0; k < num_t; k++) {
        target_temp[k] = target_spectrum[k] * table[i][j][k];
      }
      spec[i][j] = u.Integrate(target_energy, target_temp);
    }
  }
}

void InverseCompton::Spec(const std::vector<double> &primary_energy,
                          const std::vector<double> &secondary_energy,
                          const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum) {

  if (table.size() <= 0) {
    std::cout << "table size = " << table.size() << std::endl;
    std::cout << "Generate IC table " << std::endl;
    Table();
  }

  if ((table.size() != secondary_energy.size()) || (table[0].size() != primary_energy.size()) ||
      (table[0][0].size() != target_energy.size())) {
    std::cout << "table size : photon " << table.size() << std::endl;
    std::cout << "table size : electron " << table[0].size() << std::endl;
    std::cout << "table size : target " << table[0][0].size() << std::endl;
    Table();
    std::cout << "table size : photon " << table.size() << std::endl;
    std::cout << "table size : electron " << table[0].size() << std::endl;
    std::cout << "table size : target  " << table[0][0].size() << std::endl;
  }
  Spec(primary_energy, secondary_energy, target_energy, target_spectrum, table, spec);
}

void InverseCompton::SpecElectron(const std::vector<double> &primary_energy,
                                  const std::vector<double> &secondary_energy,
                                  const std::vector<double> &target_energy,
                                  const std::vector<double> &target_spectrum) {
  ICETable();
  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  int num_t = target_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].assign(num_p, 0.0);
  }
  std::vector<double> target_temp(num_t, 0.0);
  double dum;
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      for (size_t k = 0; k < num_t; k++) {
        target_temp[k] = target_spectrum[k] * table[i][j][k];
      }
      spec[i][j] = u.Integrate(target_energy, target_temp);
    }
  }
}

void InverseCompton::ICETable() {
  int num_e = s.electron.momentum.size();
  int num_t = s.target.momentum.size();

  table.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    table[i].resize(num_e);
    for (size_t j = 0; j < num_e; j++) {
      table[i][j].resize(num_t, 0.);
    }
  }

  double dum;
  double C = (3.0 / 8.0) * T_c_cnst;
  double s_CM, beta, r;
  for (size_t i = 0; i < num_e; i++) {
    for (size_t j = 0; j < num_e; j++) {
      for (size_t k = 0; k < num_t; k++) {
        s_CM = e_mass * e_mass + 2.0 * s.electron.energy[j] * s.target.energy[k];
        beta = (s_CM - e_mass * e_mass) / (s_CM + e_mass * e_mass);
        r = s.electron.energy[j] / s.electron.energy[i];
        if ((r >= 1) && (r <= (1 + beta) / (1 - beta))) {
          dum = 1 / s.electron.energy[j] * (1 + beta) / beta *
                (1. / r + r + 2 * (1 - beta) / beta * (1 - r) +
                 (1 - beta) * (1 - beta) / beta / beta * (1 - r) * (1 - r));
          table[i][j][k] = C * e_mass * e_mass / s_CM * dum;
        } else {
          table[i][j][k] = 0.0;
        }
      }
    }
  }
}

void InverseCompton::Emissivity(const std::vector<std::vector<double>> &spec,
                                const std::vector<double> &primary_momentum,
                                const std::vector<double> &primary_spectrum,
                                const std::vector<double> &secondary_energy,
                                std::vector<double> &secondary_spectrum) {
  int num_s = secondary_energy.size();
  int num_p = primary_momentum.size();
  std::vector<double> primary_temp(num_p);
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      primary_temp[j] = primary_spectrum[j] * spec[i][j];
    }
    secondary_spectrum[i] = u.Integrate(primary_momentum, primary_temp);
  }
}

void InverseCompton::Emissivity(const std::vector<double> &target_energy,
                                const std::vector<double> &target_spectrum,
                                std::vector<double> &photon_energy,
                                std::vector<double> &photon_spectrum) {
  Spec(s.electron.energy, s.photon.energy, target_energy, target_spectrum);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, photon_energy, photon_spectrum);
}

void InverseCompton::Emissivity(const std::vector<double> &target_energy,
                                const std::vector<double> &target_spectrum) {
  int num_t = target_spectrum.size();
  Spec(s.electron.energy, s.photon.energy, target_energy, target_spectrum);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, s.photon.energy, s.photon.spectrum);
}

void InverseCompton::Emissivity() {
  Spec(s.electron.energy, s.photon.energy, s.target.energy, s.target.spectrum);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, s.photon.energy, s.photon.spectrum);
}

void InverseCompton::Table() {
  int num_p = s.photon.momentum.size();
  int num_e = s.electron.momentum.size();
  int num_t = s.target.momentum.size();

  table.resize(num_p);
  for (size_t i = 0; i < num_p; i++) {
    table[i].resize(num_e);
    for (size_t j = 0; j < num_e; j++) {
      table[i][j].resize(num_t, 0.);
    }
  }

  /*
  double gamma;
  double dum, x, y, w;
  for (size_t i = 0; i < num_p; i++) {
    for (size_t j = 0; j < num_e; j++) {
      gamma = s.electron.energy[j] / e_mass;
      for (size_t k = 0; k < num_t; k++) {
        x = s.target.energy[k] / e_mass;
        y = 4. * x * gamma;
        if ((gamma - s.photon.energy[i] / e_mass) == 0) {
          w = 0;
        } else {
          w = s.photon.energy[i] / 4. / s.target.energy[k] / gamma /
              (gamma - s.photon.energy[i] / e_mass);
        }
        if (gamma > 1.) {
          if (((1. / 4. / gamma / gamma <= w)) && (w <= 1.)) {
            if (y * w < 1e-10) {
              dum = (0.75 * T_c_cnst / s.target.energy[k] / gamma / gamma) *
                    (2. * w * log(w) + (1 + 2 * w) * (1 - w));
            } else {
              dum = (0.75 * T_c_cnst / s.target.energy[k] / gamma / gamma) *
                    (2. * w * log(w) + (1 + 2 * w) * (1 - w) +
                     0.5 * y * w * y * w * (1 - w) / (1. + y * w));
            }
          } else {
            dum = 0;
          }
        } else {
          dum = 0;
        }
        table[i][j][k] = dum;
      }
    }
  }
  */

  double gamma;
  double dum, q, Gamma_e, x, FT, FKN, target_min, target_max;
  for (size_t i = 0; i < num_p; i++) {
    for (size_t j = 0; j < num_e; j++) {
      gamma = s.electron.energy[j] / e_mass;
      for (size_t k = 0; k < num_t; k++) {
        if (gamma >= 1) {
          x = s.photon.energy[i] / e_mass / gamma;
          // target_min = x / 4. / gamma / (1 - x);
          // target_max = x * gamma / (1 - x);
          // if ((s.target.energy[k] / e_mass > target_min) &&
          //     (s.target.energy[k] / e_mass < target_max)) {
          Gamma_e = 4. * gamma * s.target.energy[k] / e_mass;
          q = x / (Gamma_e * (1.0 - x));
          if ((q >= 1. / (4 * gamma * gamma)) && (q <= 1.)) {
            FT = 2. * q * log(q) + (1 + 2 * q) * (1. - q);
            FKN = 0.5 * (Gamma_e * Gamma_e * q * q) / (1 + Gamma_e * q) * (1. - q);
            dum = (0.75 * T_c_cnst / s.target.energy[k] / gamma / gamma) * (FT + FKN);
            table[i][j][k] = dum;
          } else {
            table[i][j][k] = 0.0;
          }
        } else {
          table[i][j][k] = 0.0;
        }
      }
    }
  }
}

void InverseCompton::EICCrossTable() {
  int num_e = s.electron.momentum.size();
  int num_t = s.target.momentum.size();
  EICCross_table.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    EICCross_table[i].resize(num_t);
  }

  double x, dum, gamma, beta, mandels1, mandels2;
  for (size_t k = 0; k < num_e; k++) {
    for (size_t i = 0; i < num_t; i++) {
      if (s.electron.energy[k] > e_mass) {
        gamma = s.electron.energy[k] / e_mass;
        beta = sqrt(1 - 1. / (gamma * gamma));
        dum = 0.;
        for (size_t j = 0; j < angular_num - 1; j++) {
          mandels1 = 1. + 2. * s.electron.energy[k] * s.target.energy[i] * (1 - beta * angular[j]) /
                              (e_mass * e_mass);
          mandels2 = 1. + 2. * s.electron.energy[k] * s.target.energy[i] *
                              (1 - beta * angular[j + 1]) / (e_mass * e_mass);
          dum += 0.5 * angular_bin[j] * (beta * beta / (1 - 1 / gamma)) *
                 ((1 - beta * angular[j]) * EICCross(mandels1) +
                  (1 - beta * angular[j + 1]) * EICCross(mandels2));
        }
        dum *= 0.5;
        EICCross_table[k][i] = dum;
      } else {
        EICCross_table[k][i] = 0.;
      }
    }
  }
}

void InverseCompton::ICCrossTable() {
  int num_e = s.electron.momentum.size();
  int num_t = s.target.momentum.size();
  ICCross_table.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    ICCross_table[i].resize(num_t);
  }
  double x, dum, gamma, beta, mandels1, mandels2;
  for (size_t k = 0; k < num_e; k++) {
    for (size_t i = 0; i < num_t; i++) {
      if (s.electron.energy[k] > e_mass) {
        gamma = s.electron.energy[k] / e_mass;
        beta = sqrt(1 - 1. / (gamma * gamma));
        dum = 0.;
        for (size_t j = 0; j < angular_num - 1; j++) {
          mandels1 = 1. + 2. * s.electron.energy[k] * s.target.energy[i] * (1 - beta * angular[j]) /
                              (e_mass * e_mass);
          mandels2 = 1. + 2. * s.electron.energy[k] * s.target.energy[i] *
                              (1 - beta * angular[j + 1]) / (e_mass * e_mass);
          dum += 0.5 * angular_bin[j] * (beta * beta / (1 - 1 / gamma)) *
                 ((1 - beta * angular[j]) * ICCross(mandels1) +
                  (1 - beta * angular[j + 1]) * ICCross(mandels2));
        }
        dum *= 0.5;
        ICCross_table[k][i] = dum;
      } else {
        ICCross_table[k][i] = 0.;
      }
    }
  }
}

void InverseCompton::ReadCross() {
  mandelstam.resize(mandelstam_num);
  ICcrossT.resize(mandelstam_num);
  EICcrossT.resize(mandelstam_num);
  double dum;
  std::ifstream gin1(std::string(DATA_PATH) + std::string("/data/EMProcess/mandelstam_table.dat"));
  if (gin1.good()) {
    for (int i = 0; i < mandelstam_num; i++) {
      gin1 >> mandelstam[i] >> dum >> dum >> ICcrossT[i] >> EICcrossT[i];
    }
  } else {
    std::cout << "ICset files does't exists!!!" << std::endl;
    exit(1);
  }
  gin1.close();
}

void InverseCompton::setAngular() {
  angular.resize(angular_num);
  angular_bin.resize(angular_num);
  for (size_t i = 0; i < angular_num; i++) {
    angular[i] = -1. + 2. * i / (angular_num - 1.);
    angular_bin[i] = 2. / (angular_num - 1);
  }
}

void InverseCompton::ComptonLosstime(const double mue, const double density,
                                     const std::vector<double> &photon_energy,
                                     std::vector<double> &photon_losstime) {
  // used in Murase et al. 2015
  int num_p = photon_energy.size();
  photon_losstime.resize(num_p);
  for (size_t i = 0; i < num_p; i++) {
    photon_losstime[i] = ComptonEffCrossSection(photon_energy[i]) * density * c_cnst / mue;
  }
}

double InverseCompton::ComptonEffCrossSection(double energy) {
  // see Eq.46 of Murase et al. 2015, ApJ
  double NucleusCharge = 1;
  double k = energy / e_mass;
  if (k > 1e-2) {
    return 2. * PI * NucleusCharge * electron_radius * electron_radius *
           (2. * (1 + k) * (1 + k) / k / k / (1 + 2. * k) -
            ((1 + 3 * k) / (1 + 2 * k) / (1 + 2 * k) +
             (1 + k) * (2 * k * k - 2 * k - 1) / k / k / (1 + 2 * k) / (1 + 2 * k) +
             4. * k * k / 3. / (1 + 2 * k) / (1 + 2 * k) / (1 + 2 * k) +
             ((1 + k) / k / k / k + 0.5 / k / k / k - 0.5 / k) * log(1 + 2 * k)));
  } else {
    return 2. * PI * NucleusCharge * electron_radius * electron_radius * (4. / 3) * k;
  }
}

double InverseCompton::Compattenuation(double energy) {
    //return (6.022*1e23*pMass/NucleusMass)*ComptonEffCrossSection(a); //cm^2/g screen region
    return (6.022*1e23)*ComptonEffCrossSection(energy); //cm^2/g screen region
}

void InverseCompton::Test() { std::cout << "Test : InverseCompton module         " << std::endl; }
