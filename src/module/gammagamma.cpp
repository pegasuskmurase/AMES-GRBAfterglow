//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file gammagamma.hpp
// \brief Basic class for gamma-gamma pair annihilation process

#include "module/gammagamma.hpp"
#include "module/constants.hpp"
#include "module/source.hpp"
#include <fstream>

GammaGamma::GammaGamma(Source &s) : s(s) {
  setAngular();
  ReadCross();
  GGCrossTable();
  EGGCrossTable();
}

void GammaGamma::Transport(std::vector<ParticleTransport> &p,
                           const std::vector<double> &target_energy,
                           const std::vector<double> &target_spectrum, bool have_muon,
                           bool have_spec_approx) {
  

  int idx_photon=-1;
  int idx_photon_syn=-1;
  int idx_photon_IC=-1;
  int idx_photon_pm=-1;
  int idx_photon_pa=-1;
  int idx_electron=-1;
  int idx_electron_prim=-1;
  int idx_electron_pm=-1;
  int idx_electron_pa=-1;

  for (auto &particle : p) {

    if (particle.ID == 22) {
      idx_photon = &particle - &p[0];
    }

    if (particle.ID == 2201) {
      idx_photon_syn = &particle - &p[0];
    }

    if (particle.ID == 2202) {
      idx_photon_IC = &particle - &p[0];
    }

    if (particle.ID == 2203) {
      idx_photon_pm = &particle - &p[0];
    }

    if (particle.ID == 2204) {
      idx_photon_pa = &particle - &p[0];
    }

    if (particle.ID == 11) {
      idx_electron = &particle - &p[0];
    }

    if (particle.ID == 1100) {
      idx_electron_prim = &particle - &p[0];
    }

    if (particle.ID == 1101) {
      idx_electron_pm = &particle - &p[0];
    }

    if (particle.ID == 1102) {
      idx_electron_pa = &particle - &p[0];
    }
  }

  int num_photon = p[idx_photon].x.size();
  int num_electron = p[idx_electron].x.size();

  std::vector<double> photon_intetime(num_photon, 0.0);
  Intetime(target_energy, target_spectrum, p[idx_photon].x, photon_intetime);

  for (size_t i = 0; i < num_photon; i++) {
    p[idx_photon].A[i] += photon_intetime[i];
    p[idx_photon].A_dt[i] += photon_intetime[i];

    if (idx_photon_syn > 0) {
    p[idx_photon_syn].A[i] += photon_intetime[i];
    }
    if (idx_photon_IC > 0) {
    p[idx_photon_IC].A[i] += photon_intetime[i];
    }
    if (idx_photon_pm > 0) {
    p[idx_photon_pm].A[i] += photon_intetime[i];
    }
    if (idx_photon_pa > 0) {
    p[idx_photon_pa].A[i] += photon_intetime[i];
    }
  }

  if (have_spec_approx) {
    SpecApprox(p[idx_photon].x_e, p[idx_electron].x_e, p[idx_electron].dx, photon_intetime);
  } else {
    Spec(p[idx_photon].x_e, p[idx_electron].x_e, target_energy, target_spectrum);
    for (size_t i = 0; i < num_electron; i++) {
      for (size_t j = 0; j < num_photon; j++) {
        spec[i][j] *= p[idx_electron].x[i] / p[idx_electron].x_e[i];
      }
    }
  }

  p[idx_electron].C_id.push_back(idx_photon);
  p[idx_electron].C.push_back(spec);

  if (idx_photon_syn > 0) {
  p[idx_electron_prim].C_id.push_back(idx_photon_syn);
  p[idx_electron_prim].C.push_back(spec);
  }

  if (idx_photon_IC > 0) {
  p[idx_electron_prim].C_id.push_back(idx_photon_IC);
  p[idx_electron_prim].C.push_back(spec);
  }

  if (idx_photon_pm > 0) {
  p[idx_electron_pm].C_id.push_back(idx_photon_pm);
  p[idx_electron_pm].C.push_back(spec);
  }

  if (idx_photon_pa > 0) {
  p[idx_electron_pa].C_id.push_back(idx_photon_pa);
  p[idx_electron_pa].C.push_back(spec);
  }
}

void GammaGamma::Transport(std::vector<ParticleTransportChangCooper> &p,
                           const std::vector<double> &target_energy,
                           const std::vector<double> &target_spectrum) {
  int idx_photon = 0;
  int idx_electron = 1;

  int num_photon = p[idx_photon].x_e.size();
  int num_electron = p[idx_electron].x_e.size();

  std::vector<double> photon_intetime(num_photon, 0.0);
  Intetime(target_energy, target_spectrum, p[idx_photon].x, photon_intetime);

  for (size_t i = 0; i < num_photon; i++) {
    p[idx_photon].T[i] += photon_intetime[i];
    p[idx_photon].A_dt[i] += photon_intetime[i];
  }

  /*
  Spec(p[idx_photon].x_e, p[idx_electron].x_e, target_energy, target_spectrum);
  p[idx_electron].C_id.push_back(idx_photon);
  p[idx_electron].C.push_back(spec);
  */

  SpecApprox(p[idx_photon].x_e, p[idx_electron].x_e, p[idx_electron].dx, photon_intetime);
  p[idx_electron].C_id.push_back(idx_photon);
  p[idx_electron].C.push_back(spec);
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
     //std::cout << "table size = " << EGGCross_table.size() << std::endl;
     //std::cout << "Generate EGG cross table " << std::endl;
    EGGCrossTable();
  }

  if ((EGGCross_table.size() != s.photon.momentum.size()) ||
      (EGGCross_table[0].size() != s.target.momentum.size())) {
    //std::cout << "table size : photon " << EGGCross_table.size() << std::endl;
    //std::cout << "table size : target " << EGGCross_table[0].size() << std::endl;
    EGGCrossTable();
    //std::cout << "table size : photon " << EGGCross_table.size() << std::endl;
    //std::cout << "table size : target  " << EGGCross_table[0].size() << std::endl;
  }

  Losstime(EGGCross_table, s.target.energy, s.target.spectrum, s.photon.energy, s.photon.losstime);
}

void GammaGamma::Losstime(const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum,
                          const std::vector<double> &photon_energy,
                          std::vector<double> &photon_losstime) {
  if (EGGCross_table.size() <= 0) {
    std::cout << "table size = " << EGGCross_table.size() << std::endl;
    std::cout << "Generate EGG cross table " << std::endl;
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
    std::cout << "Generate GG cross table " << std::endl;
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
    std::cout << "Generate GG cross table " << std::endl;
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
                            const std::vector<double> &secondary_momentum_bin,
                            const std::vector<double> &photon_intetime) {
  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].resize(num_p);
  }
  double r = 0.5;
  for (size_t i = 1; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      if ((secondary_energy[i - 1] <= r * primary_energy[j]) &&
          (secondary_energy[i] > r * primary_energy[j])) {
        spec[i][j] = 2 * photon_intetime[j] / secondary_momentum_bin[i];
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

void GammaGamma::Emissivity() {
  Spec(s.photon.energy, s.electron.energy, s.target.energy, s.target.spectrum);
  Emissivity(spec, s.photon.energy, s.photon.spectrum, s.electron.momentum, s.electron.spectrum);
}

void GammaGamma::EmissivityMuon() {
  SpecMuon(s.photon.energy, s.muon.energy, s.target.energy, s.target.spectrum);
  Emissivity(spec, s.photon.energy, s.photon.spectrum, s.muon.energy, s.muon.spectrum);
}

void GammaGamma::Emissivity(const std::vector<std::vector<double>> &spec,
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

double GammaGamma::TableMaxEe(double Eg, double eps) {
  double dum_max;
  double Emax = 0.5 * Eg * (1 + sqrt(1 - e_mass * e_mass / eps / Eg)) * 0.95;
  dum_max =
      (3.0 / 32.0) * T_c_cnst * e_mass * e_mass * e_mass * e_mass / (eps * eps * Eg * Eg * Eg) *
      (4 * Eg * Eg / (Eg - Emax) / Emax *
           log(4 * eps * Emax * (Eg - Emax) / (e_mass * e_mass * Eg)) -
       8.0 * eps * Eg / e_mass / e_mass +
       (2.0 * Eg * Eg * (2 * eps * Eg - e_mass * e_mass)) / ((Eg - Emax) * Emax * e_mass * e_mass) -
       (1 - e_mass * e_mass / eps / Eg) *
           (Eg * Eg * Eg * Eg / ((Eg - Emax) * (Eg - Emax) * Emax * Emax)));
  return 2 * dum_max;
}

double GammaGamma::Table(double Ee, double Eg, double eps) {
  double dum, dum_max;
  double Emax = 0.5 * Eg * (1 + sqrt(1 - e_mass * e_mass / eps / Eg));
  double Emin = 0.5 * Eg * (1 - sqrt(1 - e_mass * e_mass / eps / Eg));
  if ((Ee <= Emax) && (Ee >= Emin)) {
    dum = (3.0 / 32.0) * T_c_cnst * e_mass * e_mass * e_mass * e_mass / (eps * eps * Eg * Eg * Eg) *
          (4 * Eg * Eg / (Eg - Ee) / Ee * log(4 * eps * Ee * (Eg - Ee) / (e_mass * e_mass * Eg)) -
           8.0 * eps * Eg / e_mass / e_mass +
           (2.0 * Eg * Eg * (2 * eps * Eg - e_mass * e_mass)) / ((Eg - Ee) * Ee * e_mass * e_mass) -
           (1 - e_mass * e_mass / eps / Eg) *
               (Eg * Eg * Eg * Eg / ((Eg - Ee) * (Eg - Ee) * Ee * Ee)));
    return 2.0 * dum; // two electrons
  } else {
    return 0.0;
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

  /*
  double Emax;
  for (size_t j = 0; j < num_p; j++) {
    Eg = s.photon.energy[j];
    for (size_t k = 0; k < num_t; k++) {
      eps = s.target.energy[k];
      Emax = 0.5 * Eg * (1 + sqrt(1 - e_mass * e_mass / eps / Eg));
      if (Emax > 0.1 * Eg) {
        dum = TableMaxEe(Eg, eps);
        int idx = 0;
        for (size_t i = 0; i < num_e; i++) {
          Ee = s.electron.energy[i];
          if ((Ee < Emax)) {
            idx++;
          } else {
            break;
          }
        }
        if (idx + 1 < num_e) {
          // std::cout << "dum " << idx << " " << table[idx-1][j][k] << " " << table[idx][j][k] <<
          // std::endl;
          table[idx + 1][j][k] = dum;
        }
        // std::cout << table[idx + 1][j][k] << std::endl;
      }
    }
  }
  */
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
  // from Zdziarski, 1988, ApJ, Eq. B1 (Also from Kalasheve)
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

double GammaGamma::GGCrossTable(double Ep, double eps) {
  double dum, mandels1, mandels2;
  dum = 0.;
  for (size_t j = 0; j < angular_num - 1; j++) {
    mandels1 = 2. * Ep * eps * (1 - angular[j]) / (e_mass * e_mass);
    mandels2 = 2. * Ep * eps * (1 - angular[j + 1]) / (e_mass * e_mass);
    dum += 0.5 * angular_bin[j] *
           ((1 - angular[j]) * GGCross(mandels1) + (1 - angular[j + 1]) * GGCross(mandels2));
  }
  dum *= 0.5;
  return dum;
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

double GammaGamma::PPCross(double mandels) {
  /*
  if (mandels > 4) {
  double mandels2 = mandels * mandels;
  double lnmandels = log(mandels);
  double smandels = sqrt(mandels);
  return 1.5 * T_c_cnst / mandels2 *
         ((mandels - 1 + 0.5 / mandels - lnmandels / 2.0 + log(2.0)) *
              log(smandels + sqrt(mandels - 1)) +
          lnmandels * lnmandels / 8.0 -
          (log(smandels + sqrt(mandels - 1))) *
              (log(smandels + sqrt(mandels - 1))) / 2.0 +
          (log(2.0) * lnmandels) / 2.0 - sqrt(mandels2 - mandels));
  } else {
    return 0.0;
  }
  */
  // Lee1998
  if (mandels > 4) {
    double beta2 = (1 - 4.0 / mandels);
    double beta = sqrt(beta2);
    double beta4 = beta2 * beta2;
    if (beta < 1.0) {
      return T_c_cnst * 3.0 / 16.0 * (1 - beta2) *
             ((3 - beta4) * log((1 + beta) / (1 - beta)) - 2.0 * beta * (2 - beta2));
    } else {
      return 0.0;
    }
  } else {
    return 0.0;
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
