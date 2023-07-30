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

void GammaGamma::IntetimeMuon() {
  if (GGCross_table_muon.size() <= 0) {
    std::cout << "table size = " << GGCross_table_muon.size() << std::endl;
    std::cout << "Generate GG table " << std::endl;
    GGCrossTableMuon();
  }

  if ((GGCross_table_muon.size() != s.photon.momentum.size()) ||
      (GGCross_table_muon[0].size() != s.target.momentum.size())) {
    std::cout << "table size : photon " << GGCross_table_muon.size() << std::endl;
    std::cout << "table size : target " << GGCross_table_muon[0].size() << std::endl;
    GGCrossTableMuon();
    std::cout << "table size : photon " << GGCross_table_muon.size() << std::endl;
    std::cout << "table size : target  " << GGCross_table_muon[0].size() << std::endl;
  }

  Losstime(GGCross_table_muon, s.target.energy, s.target.spectrum, s.photon.energy,
           s.photon.intetime);
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

void GammaGamma::IntetimeMuon(const std::vector<double> &target_energy,
                              const std::vector<double> &target_spectrum,
                              const std::vector<double> &photon_energy,
                              std::vector<double> &photon_intetime) {
  if (GGCross_table_muon.size() <= 0) {
    std::cout << "table size = " << GGCross_table_muon.size() << std::endl;
    std::cout << "Generate GG table " << std::endl;
    GGCrossTableMuon();
  }

  if ((GGCross_table_muon.size() != s.photon.momentum.size()) ||
      (GGCross_table_muon[0].size() != s.target.momentum.size())) {
    std::cout << "table size : photon " << GGCross_table_muon.size() << std::endl;
    std::cout << "table size : target " << GGCross_table_muon[0].size() << std::endl;
    GGCrossTableMuon();
    std::cout << "table size : photon " << GGCross_table_muon.size() << std::endl;
    std::cout << "table size : target  " << GGCross_table_muon[0].size() << std::endl;
  }

  Losstime(GGCross_table_muon, target_energy, target_spectrum, photon_energy, photon_intetime);
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

void GammaGamma::Spec(const std::vector<double> &primary_energy,
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
    spec[i].resize(num_p);
  }
  std::vector<double> target_temp(num_t);
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      for (size_t k = 0; k < num_t; k++) {
        target_temp[k] = target_spectrum[k] * table[i][j][k];
      }
      spec[i][j] = u.Integrate(target_energy, target_temp);
    }
  }
  std::vector<double>().swap(target_temp);
}

void GammaGamma::Spec(const std::vector<double> &primary_energy,
                      const std::vector<double> &secondary_energy,
                      const std::vector<double> &target_energy,
                      const std::vector<double> &target_spectrum) {
  if (table.size() <= 0) {
    std::cout << "Generate GG table " << std::endl;
    Table();
  }
  if ((table.size() != s.electron.momentum.size()) ||
      (table[0].size() != s.photon.momentum.size()) ||
      (table[0][0].size() != s.target.momentum.size())) {
    std::cout << "table size : electron " << table.size() << std::endl;
    std::cout << "table size : photon " << table[0].size() << std::endl;
    std::cout << "table size : target " << table[0][0].size() << std::endl;
    Table();
    std::cout << "table size : electron " << table.size() << std::endl;
    std::cout << "table size : photon " << table[0].size() << std::endl;
    std::cout << "table size : target  " << table[0][0].size() << std::endl;
  }

  Spec(primary_energy, secondary_energy, target_energy, target_spectrum, table, spec);
}

void GammaGamma::SpecMuon(const std::vector<double> &primary_energy,
                          const std::vector<double> &secondary_energy,
                          const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum,
                          const std::vector<std::vector<std::vector<double>>> &table_muon,
                          std::vector<std::vector<double>> &spec) {
  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  int num_t = target_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].resize(num_p);
  }
  std::vector<double> target_temp(num_t);
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      for (size_t k = 0; k < num_t; k++) {
        target_temp[k] = target_spectrum[k] * table_muon[i][j][k];
      }
      spec[i][j] = u.Integrate(target_energy, target_temp);
    }
  }
  std::vector<double>().swap(target_temp);
}

void GammaGamma::SpecMuon(const std::vector<double> &primary_energy,
                          const std::vector<double> &secondary_energy,
                          const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum) {
  SpecMuon(primary_energy, secondary_energy, target_energy, target_spectrum, table_muon, spec);
}

void GammaGamma::SpecApprox(const std::vector<double> &primary_energy,
                            const std::vector<double> &secondary_energy,
                            const std::vector<double> &photon_intetime) {
  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].resize(num_p);
  }
  double r = 0.5;
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      if ((secondary_energy[i - 1] <= r * primary_energy[j]) &&
          (secondary_energy[i] > r * primary_energy[j])) {
        spec[i][j] = 2 * photon_intetime[j];
      } else {
        spec[i][j] = 0.0;
      }
    }
  }
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

void GammaGamma::Table() {
  //
  int num_p = s.photon.momentum.size();
  int num_e = s.electron.momentum.size();
  int num_t = s.target.momentum.size();

  table.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    table[i].resize(num_p);
    for (size_t j = 0; j < num_p; j++) {
      table[i][j].resize(num_t, 0.);
    }
  }

  double dum;
  double Eg, Ee, eps;
  for (size_t i = 0; i < num_e; i++) {
    Ee = s.electron.energy[i];
    for (size_t j = 0; j < num_p; j++) {
      Eg = s.photon.energy[j];
      if (Eg > Ee) {
        for (size_t k = 0; k < num_t; k++) {
          eps = s.target.energy[k];
          if ((Ee <= 0.5 * Eg * (1 + sqrt(1 - e_mass * e_mass / eps / Eg))) &&
              (Ee >= 0.5 * Eg * (1 - sqrt(1 - e_mass * e_mass / eps / Eg)))) {
            dum = (3.0 / 32.0) * T_c_cnst * e_mass * e_mass * e_mass * e_mass /
                  (eps * eps * Eg * Eg * Eg) *
                  (4 * Eg * Eg / (Eg - Ee) / Ee *
                       log(4 * eps * Ee * (Eg - Ee) / (e_mass * e_mass * Eg)) -
                   8.0 * eps * Eg / e_mass / e_mass +
                   (2.0 * Eg * Eg * (2 * eps * Eg - e_mass * e_mass)) /
                       ((Eg - Ee) * Ee * e_mass * e_mass) -
                   (1 - e_mass * e_mass / eps / Eg) *
                       (Eg * Eg * Eg * Eg / ((Eg - Ee) * (Eg - Ee) * Ee * Ee)));
            table[i][j][k] = 2.0 * dum; // two electrons
          } else {
            table[i][j][k] = 0.0;
          }
        }
      } else {
        for (size_t k = 0; k < num_t; k++) {
          table[i][j][k] = 0.0;
        }
      }
    }
  }
}

void GammaGamma::TableMuon() {
  //
  int num_p = s.photon.momentum.size();
  int num_m = s.muonminusL.momentum.size();
  int num_t = s.target.momentum.size();

  table_muon.resize(num_m);
  for (size_t i = 0; i < num_m; i++) {
    table_muon[i].resize(num_p);
    for (size_t j = 0; j < num_p; j++) {
      table_muon[i][j].resize(num_t, 0.);
    }
  }

  double dum;
  double Eg, Emu, eps;
  double r = (e_mass / muon_mass) * (e_mass / muon_mass);
  for (size_t i = 0; i < num_m; i++) {
    Emu = s.muonminusL.energy[i];
    for (size_t j = 0; j < num_p; j++) {
      Eg = s.photon.energy[j];
      if (Eg > Emu) {
        for (size_t k = 0; k < num_t; k++) {
          eps = s.target.energy[k];
          if ((Emu <= 0.5 * Eg * (1 + sqrt(1 - muon_mass * muon_mass / eps / Eg))) &&
              (Emu >= 0.5 * Eg * (1 - sqrt(1 - muon_mass * muon_mass / eps / Eg)))) {
            dum = (3.0 / 32.0) * T_c_cnst * muon_mass * muon_mass * muon_mass * muon_mass /
                  (eps * eps * Eg * Eg * Eg) *
                  (4 * Eg * Eg / (Eg - Emu) / Emu *
                       log(4 * eps * Emu * (Eg - Emu) / (muon_mass * muon_mass * Eg)) -
                   8.0 * eps * Eg / muon_mass / muon_mass +
                   (2.0 * Eg * Eg * (2 * eps * Eg - muon_mass * muon_mass)) /
                       ((Eg - Emu) * Emu * muon_mass * muon_mass) -
                   (1 - muon_mass * muon_mass / eps / Eg) *
                       (Eg * Eg * Eg * Eg / ((Eg - Emu) * (Eg - Emu) * Emu * Emu)));
            table_muon[i][j][k] = 2.0 * dum * r; // two muons
          } else {
            table_muon[i][j][k] = 0.0;
          }
        }
      } else {
        for (size_t k = 0; k < num_t; k++) {
          table_muon[i][j][k] = 0.0;
        }
      }
    }
  }
}

void GammaGamma::Table_Lee1998() {
  // from Lee 1998 Eq. 21
  int num_p = s.photon.momentum.size();
  int num_e = s.electron.momentum.size();
  int num_t = s.target.momentum.size();

  table.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    table[i].resize(num_p);
    for (size_t j = 0; j < num_p; j++) {
      table[i][j].resize(num_t, 0.);
    }
  }

  double dum;
  double Eg, Ee, eps;
  double beta2, s_CM;
  double C = (3.0 / 4.0) * T_c_cnst;
  for (size_t i = 0; i < num_e; i++) {
    Ee = s.electron.energy[i];
    for (size_t j = 0; j < num_p; j++) {
      Eg = s.photon.energy[j];
      if (Eg > Ee) {
        for (size_t k = 0; k < num_t; k++) {
          eps = s.target.energy[k];
          s_CM = 4.0 * Eg * eps; // TODO
          if (s_CM > 4.0 * e_mass * e_mass) {
            beta2 = (1 - 4.0 * e_mass * e_mass / s_CM);
            if ((Ee / Eg >= 0.5 * (1 - sqrt(beta2))) && (Ee / Eg <= 0.5 * (1 + sqrt(beta2)))) {
              dum = 1. / Eg *
                    (Ee / (Eg - Ee) + (Eg - Ee) / Ee + Eg * (1 - beta2) * (1 / Ee + 1 / (Eg - Ee)) -
                     0.25 * (Eg * Eg) * (1 - beta2) * (1 - beta2) * (1. / Ee + 1. / (Eg - Ee)) *
                         (1. / Ee + 1. / (Eg - Ee)));
              table[i][j][k] = C * e_mass * e_mass / s_CM * dum;
            } else {
              table[i][j][k] = 0.0;
            }
          } else {
            table[i][j][k] = 0.0;
          }
        }
      } else {
        for (size_t k = 0; k < num_t; k++) {
          table[i][j][k] = 0.0;
        }
      }
    }
  }
}

void GammaGamma::Table_Zdziarski1988() {
  // from Zdziarski, 1988, ApJ, Eq. B1 (Also from Kalosheve)
  int num_p = s.photon.momentum.size();
  int num_e = s.electron.momentum.size();
  int num_t = s.target.momentum.size();

  table.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    table[i].resize(num_p);
    for (size_t j = 0; j < num_p; j++) {
      table[i][j].resize(num_t, 0.);
    }
  }

  double dum;
  double Eg, Ee, eps;
  double E, Estar, r, Ee2, f;
  double C = (3.0 / 4.0) * T_c_cnst;
  for (size_t i = 0; i < num_e; i++) {
    Ee = s.electron.energy[i];
    for (size_t j = 0; j < num_p; j++) {
      Eg = s.photon.energy[j];
      if (Eg > Ee) {
        Ee2 = Eg - Ee;
        for (size_t k = 0; k < num_t; k++) {
          eps = s.target.energy[k];
          E = Eg * eps / (e_mass * e_mass);
          Estar = 0.25 * Eg * Eg / Ee / Ee2;
          f = Estar / E;
          if ((Estar < E) && (Estar > 1)) {
            r = 0.5 * (Ee / Ee2 + Ee2 / Ee);
            dum = r + f * (2 * f - (2 + r) - 2 * log(f));
            table[i][j][k] = 2.0 * C / E / Eg * dum; // two electrons
          } else {
            table[i][j][k] = 0.0;
          }
        }
      } else {
        for (size_t k = 0; k < num_t; k++) {
          table[i][j][k] = 0.0;
        }
      }
    }
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

void GammaGamma::GGCrossTableMuon() {
  int num_p = s.photon.momentum.size();
  int num_t = s.target.momentum.size();
  GGCross_table_muon.resize(num_p);
  for (int i = 0; i < num_p; ++i) {
    GGCross_table_muon[i].resize(num_t);
  }
  double results, dum, mandels1, mandels2;
  double r = e_mass / muon_mass;
  for (size_t k = 0; k < num_p; k++) {
    for (size_t i = 0; i < num_t; i++) {
      dum = 0.;
      for (size_t j = 0; j < angular_num - 1; j++) {
        mandels1 = 2. * s.photon.energy[k] * s.target.energy[i] * (1 - angular[j]) /
                   (muon_mass * muon_mass);
        mandels2 = 2. * s.photon.energy[k] * s.target.energy[i] * (1 - angular[j + 1]) /
                   (muon_mass * muon_mass);
        dum += 0.5 * angular_bin[j] *
               ((1 - angular[j]) * GGCross(mandels1) * r * r +
                (1 - angular[j + 1]) * GGCross(mandels2) * r * r);
      }
      dum *= 0.5;
      GGCross_table_muon[k][i] = dum;
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
