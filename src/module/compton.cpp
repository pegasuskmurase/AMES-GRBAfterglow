//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file inversecompton.hpp
// \brief Basic class for inversecompton process

#include "module/compton.hpp"
#include "module/constants.hpp"
#include "module/source.hpp"
#include "module/utility.hpp"

#include <fstream> // Add this line

InverseCompton::InverseCompton(Source &s) : s(s) {
  setAngular();
  ReadCross();
  ICCrossTable();
  EICCrossTable();
  Table();
}

void InverseCompton::Transport(std::vector<ParticleTransport> &p,
                               const std::vector<double> &target_energy,
                               const std::vector<double> &target_spectrum,
                               bool have_compton_attenu,
                               bool have_continue_loss, bool have_AIC,
                               bool have_muon) {
  // Find particle indices first
  int idx_photon = -1, idx_photon_IC = -1, idx_electron = -1;
  int idx_photon_pm = -1;
  int idx_photon_pa = -1;
  int idx_photon_pp = -1;
  int idx_electron_prim = -1;
  int idx_electron_pm = -1;
  int idx_electron_pa = -1;
  int idx_electron_pp = -1;
  int idx_muonminusL = -1;
  int idx_muonminusR = -1;
  int idx_muonplusL = -1;
  int idx_muonplusR = -1;
  for (size_t i = 0; i < p.size(); ++i) {
    if (p[i].ID == 22)
      idx_photon = static_cast<int>(i);
    if (p[i].ID == 2202)
      idx_photon_IC = static_cast<int>(i);
    if (p[i].ID == 2203)
      idx_photon_pm = static_cast<int>(i);
    if (p[i].ID == 2204)
      idx_photon_pa = static_cast<int>(i);
    if (p[i].ID == 2205)
      idx_photon_pp = static_cast<int>(i);
    if (p[i].ID == 11)
      idx_electron = static_cast<int>(i);
    if (p[i].ID == 1100)
      idx_electron_prim = static_cast<int>(i);
    if (p[i].ID == 1101)
      idx_electron_pm = static_cast<int>(i);
    if (p[i].ID == 1102)
      idx_electron_pa = static_cast<int>(i);
    if (p[i].ID == 1103)
      idx_electron_pp = static_cast<int>(i);
    if (p[i].ID == 100041)
      idx_muonminusL = static_cast<int>(i);
    if (p[i].ID == 100042)
      idx_muonminusR = static_cast<int>(i);
    if (p[i].ID == 100051)
      idx_muonplusL = static_cast<int>(i);
    if (p[i].ID == 100052)
      idx_muonplusR = static_cast<int>(i);
    // if (idx_photon != -1 && idx_electron != -1) break;  // Early exit if
    // both found
  }

  if (idx_electron == -1 && idx_photon == -1 && idx_muonminusL == -1 &&
      idx_muonminusR == -1 && idx_muonplusL == -1 && idx_muonplusR == -1) {
    std::cout << "IC Warning: Essential particles missing from simulation\n";
    return;
  } // No relevant particles

  // Pre-compute sizes and references
  const size_t num_electron =
      (idx_electron != -1) ? p[idx_electron].x.size() : 0;
  const size_t num_photon = (idx_photon != -1) ? p[idx_photon].x.size() : 0;
  // Electron processing
  if (idx_electron != -1) {
    auto &electron = p[idx_electron];
    std::vector<double> electron_losstime(num_electron, 0.0);
    std::vector<double> electron_intetime(num_electron, 0.0);

    Losstime(target_energy, target_spectrum, electron.x_e, electron_losstime);
    Intetime(target_energy, target_spectrum, electron.x_e, electron_intetime);

    if (have_continue_loss) {

      // Pre-compute frequently used values
      const auto &x = electron.x;
      const auto &dx = electron.dx;
      auto &A = electron.A;
      auto &A_dt = electron.A_dt;
      auto &B = electron.B;

#pragma omp parallel for
      for (size_t i = 0; i < num_electron; ++i) {
        const double term = electron_losstime[i] * x[i] / dx[i];
        A[i] += term;
        A_dt[i] += electron_losstime[i];
      }

// Process B matrix
#pragma omp parallel for
      for (size_t i = 0; i < num_electron - 1; ++i) {
        const size_t j = i + 1;
        B[i][j] += electron_losstime[j] * x[j] / dx[j];
      }

      if (idx_electron_prim > 0) {
        auto &electron = p[idx_electron_prim];
        // Pre-compute frequently used values
        const auto &x = electron.x;
        const auto &dx = electron.dx;
        auto &A = electron.A;
        auto &A_dt = electron.A_dt;
        auto &B = electron.B;

#pragma omp parallel for
        for (size_t i = 0; i < num_electron; ++i) {
          const double term = electron_losstime[i] * x[i] / dx[i];
          A[i] += term;
          A_dt[i] += electron_losstime[i];
        }

// Process B matrix
#pragma omp parallel for
        for (size_t i = 0; i < num_electron - 1; ++i) {
          const size_t j = i + 1;
          B[i][j] += electron_losstime[j] * x[j] / dx[j];
        }
      }

      if (idx_electron_pa > 0) {
        auto &electron = p[idx_electron_pa];
        // Pre-compute frequently used values
        const auto &x = electron.x;
        const auto &dx = electron.dx;
        auto &A = electron.A;
        auto &A_dt = electron.A_dt;
        auto &B = electron.B;

#pragma omp parallel for
        for (size_t i = 0; i < num_electron; ++i) {
          const double term = electron_losstime[i] * x[i] / dx[i];
          A[i] += term;
          A_dt[i] += electron_losstime[i];
        }

// Process B matrix
#pragma omp parallel for
        for (size_t i = 0; i < num_electron - 1; ++i) {
          const size_t j = i + 1;
          B[i][j] += electron_losstime[j] * x[j] / dx[j];
        }
      }

      if (idx_electron_pm > 0) {
        auto &electron = p[idx_electron_pm];
        // Pre-compute frequently used values
        const auto &x = electron.x;
        const auto &dx = electron.dx;
        auto &A = electron.A;
        auto &A_dt = electron.A_dt;
        auto &B = electron.B;

#pragma omp parallel for
        for (size_t i = 0; i < num_electron; ++i) {
          const double term = electron_losstime[i] * x[i] / dx[i];
          A[i] += term;
          A_dt[i] += electron_losstime[i];
        }

// Process B matrix
#pragma omp parallel for
        for (size_t i = 0; i < num_electron - 1; ++i) {
          const size_t j = i + 1;
          B[i][j] += electron_losstime[j] * x[j] / dx[j];
        }
      }

      if (idx_electron_pp > 0) {
        auto &electron = p[idx_electron_pp];
        // Pre-compute frequently used values
        const auto &x = electron.x;
        const auto &dx = electron.dx;
        auto &A = electron.A;
        auto &A_dt = electron.A_dt;
        auto &B = electron.B;

#pragma omp parallel for
        for (size_t i = 0; i < num_electron; ++i) {
          const double term = electron_losstime[i] * x[i] / dx[i];
          A[i] += term;
          A_dt[i] += electron_losstime[i];
        }

// Process B matrix
#pragma omp parallel for
        for (size_t i = 0; i < num_electron - 1; ++i) {
          const size_t j = i + 1;
          B[i][j] += electron_losstime[j] * x[j] / dx[j];
        }
      }

    } else {
      // Pre-compute frequently used values
      const auto &x = electron.x;
      const auto &dx = electron.dx;
      auto &A = electron.A;
      auto &A_dt = electron.A_dt;
      auto &B = electron.B;

#pragma omp parallel for
      for (size_t i = 0; i < num_electron; ++i) {
        if (electron_intetime[i] / electron_losstime[i] > 1) {
          const double term = electron_losstime[i] * x[i] / dx[i];
          A[i] += term; // Continue energy loss approximation
        } else {
          std::cout << "Warning , no continus loss assumption " << std::endl;
          A[i] += electron_intetime[i];
        }
      }

// Process B matrix
#pragma omp parallel for
      for (size_t i = 0; i < num_electron - 1; ++i) {
        const size_t j = i + 1;
        if (electron_intetime[j] / electron_losstime[j] > 1) {
          B[i][j] += electron_losstime[j] * x[j] /
                     dx[j]; // Continue energy loss approximation
        }
      }

      if (idx_electron_prim > 0) {
        auto &electron = p[idx_electron_prim];
        // Pre-compute frequently used values
        const auto &x = electron.x;
        const auto &dx = electron.dx;
        auto &A = electron.A;
        auto &A_dt = electron.A_dt;
        auto &B = electron.B;

#pragma omp parallel for
        for (size_t i = 0; i < num_electron; ++i) {
          if (electron_intetime[i] / electron_losstime[i] > 1) {
            const double term = electron_losstime[i] * x[i] / dx[i];
            A[i] += term; // Continue energy loss approximation
          } else {
            std::cout << "Warning , no continus loss assumption " << std::endl;
            A[i] += electron_intetime[i];
          }
        }

// Process B matrix
#pragma omp parallel for
        for (size_t i = 0; i < num_electron - 1; ++i) {
          const size_t j = i + 1;
          if (electron_intetime[j] / electron_losstime[j] > 1) {
            B[i][j] += electron_losstime[j] * x[j] /
                       dx[j]; // Continue energy loss approximation
          }
        }
      }

      if (idx_electron_pa > 0) {
        auto &electron = p[idx_electron_pa];
        // Pre-compute frequently used values
        const auto &x = electron.x;
        const auto &dx = electron.dx;
        auto &A = electron.A;
        auto &A_dt = electron.A_dt;
        auto &B = electron.B;

#pragma omp parallel for
        for (size_t i = 0; i < num_electron; ++i) {
          if (electron_intetime[i] / electron_losstime[i] > 1) {
            const double term = electron_losstime[i] * x[i] / dx[i];
            A[i] += term; // Continue energy loss approximation
          } else {
            std::cout << "Warning , no continus loss assumption " << std::endl;
            A[i] += electron_intetime[i];
          }
        }

// Process B matrix
#pragma omp parallel for
        for (size_t i = 0; i < num_electron - 1; ++i) {
          const size_t j = i + 1;
          if (electron_intetime[j] / electron_losstime[j] > 1) {
            B[i][j] += electron_losstime[j] * x[j] /
                       dx[j]; // Continue energy loss approximation
          }
        }
      }

      if (idx_electron_pm > 0) {
        auto &electron = p[idx_electron_pm];
        // Pre-compute frequently used values
        const auto &x = electron.x;
        const auto &dx = electron.dx;
        auto &A = electron.A;
        auto &A_dt = electron.A_dt;
        auto &B = electron.B;

#pragma omp parallel for
        for (size_t i = 0; i < num_electron; ++i) {
          if (electron_intetime[i] / electron_losstime[i] > 1) {
            const double term = electron_losstime[i] * x[i] / dx[i];
            A[i] += term; // Continue energy loss approximation
          } else {
            std::cout << "Warning , no continus loss assumption " << std::endl;
            A[i] += electron_intetime[i];
          }
        }

// Process B matrix
#pragma omp parallel for
        for (size_t i = 0; i < num_electron - 1; ++i) {
          const size_t j = i + 1;
          if (electron_intetime[j] / electron_losstime[j] > 1) {
            B[i][j] += electron_losstime[j] * x[j] /
                       dx[j]; // Continue energy loss approximation
          }
        }
      }

      if (idx_electron_pp > 0) {
        auto &electron = p[idx_electron_pp];
        // Pre-compute frequently used values
        const auto &x = electron.x;
        const auto &dx = electron.dx;
        auto &A = electron.A;
        auto &A_dt = electron.A_dt;
        auto &B = electron.B;

#pragma omp parallel for
        for (size_t i = 0; i < num_electron; ++i) {
          if (electron_intetime[i] / electron_losstime[i] > 1) {
            const double term = electron_losstime[i] * x[i] / dx[i];
            A[i] += term; // Continue energy loss approximation
          } else {
            std::cout << "Warning , no continus loss assumption " << std::endl;
            A[i] += electron_intetime[i];
          }
        }

// Process B matrix
#pragma omp parallel for
        for (size_t i = 0; i < num_electron - 1; ++i) {
          const size_t j = i + 1;
          if (electron_intetime[j] / electron_losstime[j] > 1) {
            B[i][j] += electron_losstime[j] * x[j] /
                       dx[j]; // Continue energy loss approximation
          }
        }
      }
    }
  }

  if (have_compton_attenu) {
    std::vector<double> photon_losstime(num_photon, 0.0);
    ComptonLosstime(1. / 0.85, s.getParam().getAmbientDensity(),
                    p[idx_photon].x, photon_losstime);

#pragma omp parallel for
    for (size_t i = 0; i < num_photon; i++) {
      p[idx_photon].A[i] += photon_losstime[i];
    }
  }

  // Spectrum calculations
  if (idx_electron != -1 && idx_photon != -1) {
    if (have_AIC) {
      double cosangle = cos(s.getParam().getObsAngle());
      std::cout << "cosangle = " << cosangle << std::endl;
      SpecAIC(cosangle, p[idx_electron].x_e, p[idx_photon].x, target_energy,
              target_spectrum);

      if ((idx_photon >= 0) && (p[idx_electron].ID == 11)) {
        p[idx_photon].C_id.push_back(idx_electron);
        p[idx_photon].C.push_back(spec_AIC);
      }

      if (idx_photon_IC > 0) {
        p[idx_photon_IC].C_id.push_back(idx_electron_prim);
        p[idx_photon_IC].C.push_back(spec_AIC);
      }

      if (idx_photon_pm > 0) {
        p[idx_photon_pm].C_id.push_back(idx_electron_pm);
        p[idx_photon_pm].C.push_back(spec_AIC);
      }

      if (idx_photon_pa > 0) {
        p[idx_photon_pa].C_id.push_back(idx_electron_pa);
        p[idx_photon_pa].C.push_back(spec_AIC);
      }

      if (idx_photon_pp > 0) {
        p[idx_photon_pp].C_id.push_back(idx_electron_pp);
        p[idx_photon_pp].C.push_back(spec_AIC);
      }
    } else {
      Spec(p[idx_electron].x_e, p[idx_photon].x, target_energy,
           target_spectrum);

      if ((idx_photon >= 0) && (p[idx_electron].ID == 11)) {
        p[idx_photon].C_id.push_back(idx_electron);
        p[idx_photon].C.push_back(spec);
      }

      if (idx_photon_IC > 0) {
        p[idx_photon_IC].C_id.push_back(idx_electron_prim);
        p[idx_photon_IC].C.push_back(spec);
      }

      if (idx_photon_pm > 0) {
        p[idx_photon_pm].C_id.push_back(idx_electron_pm);
        p[idx_photon_pm].C.push_back(spec);
      }

      if (idx_photon_pa > 0) {
        p[idx_photon_pa].C_id.push_back(idx_electron_pa);
        p[idx_photon_pa].C.push_back(spec);
      }

      if (idx_photon_pp > 0) {
        p[idx_photon_pp].C_id.push_back(idx_electron_pp);
        p[idx_photon_pp].C.push_back(spec);
      }
    }
  }

  // considering IC process for muon under continue energy loss assumption
  // without considering AIC and Compton attenuation
  if (have_muon) {

    const size_t num_muon = p[idx_muonminusL].x_e.size();

    std::vector<double> muon_losstime(num_muon, 0.0);
    // std::vector<double> muon_intetime(num_muon, 0.0);

    Losstime(muon_mass, target_energy, target_spectrum, p[idx_muonminusL].x_e,
             muon_losstime);
    // Intetime(muon_mass, target_energy, target_spectrum,
    // p[idx_muonminusL].x_e, muon_intetime);

    // set the value of A for muon with the energy loss rate calculated
    // above We use the continues energy loss assumpton A_i = 1/(E_i-1 - E_i)
    // dE/dt|E=E_i A_i = 1/(p_i-1 - p_i) dp/dt|p=p_i
    if (have_AIC == false && have_continue_loss) {
#pragma omp parallel for
      for (size_t i = 0; i < num_muon; i++) {
        double loss = muon_losstime[i];

        p[idx_muonminusL].A[i] +=
            loss * p[idx_muonminusL].x[i] / p[idx_muonminusL].dx[i];
        p[idx_muonminusL].A_dt[i] += loss;

        p[idx_muonminusR].A[i] +=
            loss * p[idx_muonminusR].x[i] / p[idx_muonminusR].dx[i];
        p[idx_muonminusR].A_dt[i] += loss;

        p[idx_muonplusL].A[i] +=
            loss * p[idx_muonplusL].x[i] / p[idx_muonplusL].dx[i];
        p[idx_muonplusL].A_dt[i] += loss;

        p[idx_muonplusR].A[i] +=
            loss * p[idx_muonplusR].x[i] / p[idx_muonplusR].dx[i];
        p[idx_muonplusR].A_dt[i] += loss;
      }

#pragma omp parallel for
      for (size_t i = 0; i < num_muon - 1; i++) {
        int j = i + 1;
        double loss = muon_losstime[j];

        p[idx_muonminusL].B[i][j] +=
            loss * p[idx_muonminusL].x[j] /
            p[idx_muonminusL].dx[j]; // Continue energy loss approximation
        p[idx_muonminusR].B[i][j] +=
            loss * p[idx_muonminusR].x[j] /
            p[idx_muonminusR].dx[j]; // Continue energy loss approximation
        p[idx_muonplusL].B[i][j] +=
            loss * p[idx_muonplusL].x[j] /
            p[idx_muonplusL].dx[j]; // Continue energy loss approximation
        p[idx_muonplusR].B[i][j] +=
            loss * p[idx_muonplusR].x[j] /
            p[idx_muonplusR].dx[j]; // Continue energy loss approximation
      }

      Spec(muon_mass, p[idx_muonminusL].x_e, p[idx_photon].x, target_energy,
           target_spectrum);

      p[idx_photon].C_id.push_back(idx_muonminusL);
      p[idx_photon].C.push_back(spec_particle);

      p[idx_photon].C_id.push_back(idx_muonminusR);
      p[idx_photon].C.push_back(spec_particle);

      p[idx_photon].C_id.push_back(idx_muonplusL);
      p[idx_photon].C.push_back(spec_particle);

      p[idx_photon].C_id.push_back(idx_muonplusR);
      p[idx_photon].C.push_back(spec_particle);
    }
  }
}

void InverseCompton::Transport(std::vector<ParticleTransportChangCooper> &p,
                               const std::vector<double> &target_energy,
                               const std::vector<double> &target_spectrum,
                               bool have_AIC, bool have_muon) {

  int idx_photon;
  int idx_electron;

  for (auto &particle : p) {

    if (particle.ID == 22) {
      idx_photon = &particle - &p[0];
    }

    if (particle.ID == 11) {
      idx_electron = &particle - &p[0];
    }
  }

  int num_photon = p[idx_photon].x.size();
  int num_electron = p[idx_electron].x.size();

  std::vector<double> electron_losstime(num_electron, 0.0);
  Losstime(target_energy, target_spectrum, p[idx_electron].x_e,
           electron_losstime);

  for (size_t i = 0; i < p[idx_electron].x.size(); i++) {

    p[idx_electron].B[i] += electron_losstime[i] * p[idx_electron].x[i];
    p[idx_electron].A_dt[i] += electron_losstime[i];
  }

  std::vector<double> electron_losstime_ghost(3, 0.0);
  Losstime(target_energy, target_spectrum, p[idx_electron].xghost_e,
           electron_losstime_ghost);

  for (size_t i = 0; i < p[idx_electron].xghost.size(); i++) {

    p[idx_electron].Bghost[i] +=
        electron_losstime_ghost[i] * p[idx_electron].xghost[i];
  }

  // set ParticleTransport.C --- Differential cross-section of particle from
  // parent particle with different type
  Spec(p[idx_electron].x_e, p[idx_photon].x, target_energy, target_spectrum);
  p[idx_photon].C_id.push_back(idx_electron);
  p[idx_photon].C.push_back(spec);
}

// syu
void InverseCompton::Losstime(const std::vector<std::vector<double>> &table,
                              const std::vector<double> &target_energy,
                              const std::vector<double> &target_spectrum,
                              const std::vector<double> &particle_energy,
                              std::vector<double> &particle_losstime) {
  // Get sizes and validate input
  const int num_p = particle_energy.size();
  const int num_t = target_energy.size();

  // Pre-allocate output
  particle_losstime.resize(num_p);
  // Get raw pointers for faster access
  const double *target_energy_ptr = target_energy.data();
  const double *target_spectrum_ptr = target_spectrum.data();

// Parallelize with OpenMP
#pragma omp parallel for
  for (int i = 0; i < num_p; ++i) {
    std::vector<double> target_temp(num_t);
    double *target_temp_ptr = target_temp.data();

    // Vectorized computation
    const auto &table_row = table[i];
#pragma omp simd
    for (int k = 0; k < num_t; ++k) {
      target_temp_ptr[k] = target_spectrum_ptr[k] * table_row[k];
    }

    // Perform integration
    particle_losstime[i] =
        u.Integrate(target_energy_ptr, target_temp_ptr, num_t);
  }

  // No need to swap with empty vector - modern compilers optimize this well
}

void InverseCompton::Losstime() {
  Losstime(EICCross_table, s.target.energy, s.target.spectrum,
           s.electron.energy, s.electron.losstime);
}
// syu
void InverseCompton::Losstime(const std::vector<double> &target_energy,
                              const std::vector<double> &target_spectrum,
                              const std::vector<double> &electron_energy,
                              std::vector<double> &electron_losstime) {
  Losstime(EICCross_table, target_energy, target_spectrum, electron_energy,
           electron_losstime);
}

void InverseCompton::Losstime(const double particle_mass,
                              const std::vector<double> &target_energy,
                              const std::vector<double> &target_spectrum,
                              const std::vector<double> &particle_energy,
                              std::vector<double> &particle_losstime) {
  int num_e = particle_energy.size();
  int num_t = target_energy.size();
  EICCross_table_particle.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    EICCross_table_particle[i].resize(num_t);
  }

  double x, dum, gamma, beta, mandels1, mandels2;
  double r = e_mass / particle_mass;
  for (size_t k = 0; k < num_e; k++) {
    for (size_t i = 0; i < num_t; i++) {
      if (particle_energy[k] > particle_mass) {
        gamma = particle_energy[k] / particle_mass;
        beta = sqrt(1 - 1. / (gamma * gamma));
        dum = 0.;
        for (size_t j = 0; j < angular_num - 1; j++) {
          mandels1 = 1. + 2. * particle_energy[k] * target_energy[i] *
                              (1 - beta * angular[j]) /
                              (particle_mass * particle_mass);
          mandels2 = 1. + 2. * particle_energy[k] * target_energy[i] *
                              (1 - beta * angular[j + 1]) /
                              (particle_mass * particle_mass);
          dum += 0.5 * angular_bin[j] * (beta * beta / (1 - 1 / gamma)) *
                 ((1 - beta * angular[j]) * EICCross(mandels1) +
                  (1 - beta * angular[j + 1]) * EICCross(mandels2));
        }
        dum *= 0.5;
        EICCross_table_particle[k][i] = dum * r * r;
      } else {
        EICCross_table_particle[k][i] = 0.;
      }
    }
  }

  Losstime(EICCross_table_particle, target_energy, target_spectrum,
           particle_energy, particle_losstime);
}

void InverseCompton::Intetime() {
  Losstime(ICCross_table, s.target.energy, s.target.spectrum, s.electron.energy,
           s.electron.intetime);
}

void InverseCompton::Intetime(const std::vector<double> &target_energy,
                              const std::vector<double> &target_spectrum,
                              const std::vector<double> &electron_energy,
                              std::vector<double> &electron_intetime) {
  Losstime(ICCross_table, target_energy, target_spectrum, electron_energy,
           electron_intetime);
}

void InverseCompton::Intetime(const double particle_mass,
                              const std::vector<double> &target_energy,
                              const std::vector<double> &target_spectrum,
                              const std::vector<double> &particle_energy,
                              std::vector<double> &particle_intetime) {

  int num_e = particle_energy.size();
  int num_t = target_energy.size();
  ICCross_table_particle.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    ICCross_table_particle[i].resize(num_t);
  }

  double x, dum, gamma, beta, mandels1, mandels2;
  double r = e_mass / particle_mass;
  for (size_t k = 0; k < num_e; k++) {
    for (size_t i = 0; i < num_t; i++) {
      if (particle_energy[k] > particle_mass) {
        gamma = particle_energy[k] / particle_mass;
        beta = sqrt(1 - 1. / (gamma * gamma));
        dum = 0.;
        for (size_t j = 0; j < angular_num - 1; j++) {
          mandels1 = 1. + 2. * particle_energy[k] * target_energy[i] *
                              (1 - beta * angular[j]) /
                              (particle_mass * particle_mass);
          mandels2 = 1. + 2. * particle_energy[k] * target_energy[i] *
                              (1 - beta * angular[j + 1]) /
                              (particle_mass * particle_mass);
          dum += 0.5 * angular_bin[j] * (beta * beta / (1 - 1 / gamma)) *
                 ((1 - beta * angular[j]) * ICCross(mandels1) +
                  (1 - beta * angular[j + 1]) * ICCross(mandels2));
        }
        dum *= 0.5;
        ICCross_table_particle[k][i] = dum * r * r;
      } else {
        ICCross_table_particle[k][i] = 0.;
      }
    }
  }

  Losstime(ICCross_table_particle, target_energy, target_spectrum,
           particle_energy, particle_intetime);
}

void InverseCompton::Spec(
    const std::vector<double> &primary_energy,
    const std::vector<double> &secondary_energy,
    const std::vector<double> &target_energy,
    const std::vector<double> &target_spectrum,
    const std::vector<std::vector<std::vector<double>>> &table,
    std::vector<std::vector<double>> &spec) {
  // Get sizes and validate input
  const int num_s = secondary_energy.size();
  const int num_p = primary_energy.size();
  const int num_t = target_energy.size();

  // Pre-allocate memory in one go
  spec.assign(num_s, std::vector<double>(num_p, 0.0));

  // Get raw pointers for faster access
  const double *target_energy_ptr = target_energy.data();
  const double *target_spectrum_ptr = target_spectrum.data();

  assert(table.size() == num_s && "mismatch table.size with particle bins");
  assert(target_energy.size() == num_t && "mismatch target energy length");
  if (table[0].size() != num_p) {
    std::cerr << "Assertion failed: table[0].size() == num_p\n"
              << "  table[0].size() = " << table[0].size() << "\n"
              << "  num_p = " << num_p << std::endl;
    std::abort(); // equivalent to assert failure
  }
  assert(table[0][0].size() == num_t && "mismatch table with num_t");

// Parallelize outer loops with OpenMP
#pragma omp parallel for
  for (int i = 0; i < num_s; ++i) {
    // Temporary buffer for integration
    std::vector<double> target_temp(num_t);
    double *target_temp_ptr = target_temp.data();
    for (int j = 0; j < num_p; ++j) {
      // Skip if energy thresholds not met (uncomment if needed)
      // if (secondary_energy[i] <= 1e1 || primary_energy[j] <= 1e1)
      // continue; Vectorized computation of target_temp
      const auto &table_row = table[i][j];
#pragma omp simd
      for (int k = 0; k < num_t; ++k) {
        target_temp_ptr[k] = target_spectrum_ptr[k] * table_row[k];
      }

      // Perform integration
      spec[i][j] = u.Integrate(target_energy_ptr, target_temp_ptr, num_t);
    }
  }
}

void InverseCompton::Spec(const std::vector<double> &primary_energy,
                          const std::vector<double> &secondary_energy,
                          const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum) {

  Spec(primary_energy, secondary_energy, target_energy, target_spectrum, table,
       spec);
}

void InverseCompton::Spec(const double particle_mass,
                          const std::vector<double> &primary_energy,
                          const std::vector<double> &secondary_energy,
                          const std::vector<double> &target_energy,
                          const std::vector<double> &target_spectrum) {

  if (table_particle.size() <= 0) {
    std::cout << "table size = " << table_particle.size() << std::endl;
    std::cout << "Generate IC table " << std::endl;
    Table(particle_mass, primary_energy, secondary_energy, target_energy,
          target_spectrum);
  }

  if ((table_particle.size() != secondary_energy.size()) ||
      (table_particle[0].size() != primary_energy.size()) ||
      (table_particle[0][0].size() != target_energy.size())) {
    std::cout << "table size : photon " << table_particle.size() << std::endl;
    std::cout << "table size : particle " << table_particle[0].size()
              << std::endl;
    std::cout << "table size : target " << table_particle[0][0].size()
              << std::endl;
    Table(particle_mass, primary_energy, secondary_energy, target_energy,
          target_spectrum);
    std::cout << "table size : photon " << table_particle.size() << std::endl;
    std::cout << "table size : particle " << table_particle[0].size()
              << std::endl;
    std::cout << "table size : target  " << table_particle[0][0].size()
              << std::endl;
  }

  Spec(primary_energy, secondary_energy, target_energy, target_spectrum,
       table_particle, spec_particle);
}

void InverseCompton::SpecElectron(const std::vector<double> &primary_energy,
                                  const std::vector<double> &secondary_energy,
                                  const std::vector<double> &target_energy,
                                  const std::vector<double> &target_spectrum) {

  if (ICE_table.size() <= 0) {
    std::cout << "Generate IC table " << std::endl;
    ICETable();
  }

  if ((ICE_table.size() != secondary_energy.size()) ||
      (ICE_table[0].size() != primary_energy.size()) ||
      (ICE_table[0][0].size() != target_energy.size())) {
    std::cout << "table size : photon " << ICE_table.size() << std::endl;
    std::cout << "table size : electron " << ICE_table[0].size() << std::endl;
    std::cout << "table size : target " << ICE_table[0][0].size() << std::endl;
    ICETable();
    std::cout << "table size : photon " << ICE_table.size() << std::endl;
    std::cout << "table size : electron " << ICE_table[0].size() << std::endl;
    std::cout << "table size : target  " << ICE_table[0][0].size() << std::endl;
  }

  int num_s = secondary_energy.size();
  int num_p = primary_energy.size();
  int num_t = target_energy.size();
  spec.resize(num_s);
  for (size_t i = 0; i < num_s; i++) {
    spec[i].assign(num_p, 0.0);
  }
  std::vector<double> target_temp(num_t, 0.0);
  for (size_t i = 0; i < num_s; i++) {
    for (size_t j = 0; j < num_p; j++) {
      for (size_t k = 0; k < num_t; k++) {
        target_temp[k] = target_spectrum[k] * ICE_table[i][j][k];
        // target_temp[k] = target_spectrum[k] * ICS(primary_energy[j],
        // target_energy[k],
        //                                           primary_energy[j] -
        //                                           secondary_energy[i]);
      }
      spec[i][j] = u.Integrate(target_energy, target_temp);
    }
  }
}

double InverseCompton::ICS(double E, double b, double nE) {
  double q, FT, FKN, Gamma_e;
  double x = nE / E;
  double gamma = E / e_mass;
  if (gamma >= 1) {
    Gamma_e = 4. * gamma * b / e_mass;
    q = x / (Gamma_e * (1.0 - x));
    if ((q >= 1. / (4 * gamma * gamma)) && (q <= 1.)) {
      FT = 2. * q * log(q) + (1 + 2 * q) * (1. - q);
      FKN = 0.5 * (Gamma_e * Gamma_e * q * q) / (1 + Gamma_e * q) * (1. - q);
      return (0.75 * T_c_cnst / b / gamma / gamma) * (FT + FKN);
    } else {
      return 0.0;
    }
  } else {
    return 0.0;
  }
}

void InverseCompton::ICETable() {
  int num_e = s.electron.momentum.size();
  int num_t = s.target.momentum.size();

  ICE_table.resize(num_e);
  for (size_t i = 0; i < num_e; i++) {
    ICE_table[i].resize(num_e);
    for (size_t j = 0; j < num_e; j++) {
      ICE_table[i][j].resize(num_t, 0.);
    }
  }

  // from Zdziarski, 1988, ApJ, Eq. A1
  double dum;
  double Eep, Ee, eps;
  double E, Estar, r, Ee2, f;
  double C = (3.0 / 4.0) * T_c_cnst;
  for (size_t i = 0; i < num_e; i++) {
    Eep = s.electron.energy[i];
    for (size_t j = 0; j < num_e; j++) {
      Ee = s.electron.energy[j];
      if (Ee > Eep) {
        for (size_t k = 0; k < num_t; k++) {
          eps = s.target.energy[k];
          E = Ee * eps / (e_mass * e_mass);
          Estar = 0.25 * (Ee / Eep - 1);
          f = Estar / E;
          if (Estar < E) {
            r = 0.5 * (Ee / Eep + Eep / Ee);
            dum = r + (2 - r) * f - 2 * f * f - 2 * f * log(1. / f);
            ICE_table[i][j][k] = C / E / Ee * dum;
          } else {
            ICE_table[i][j][k] = 0.0;
          }
        }
      } else {
        for (size_t k = 0; k < num_t; k++) {
          ICE_table[i][j][k] = 0.0;
        }
      }
    }
  }

  /*
  double dum;
  double C = (3.0 / 8.0) * T_c_cnst;
  double s_CM, beta, r;
  for (size_t i = 0; i < num_e; i++) {
    for (size_t j = 0; j < num_e; j++) {
      for (size_t k = 0; k < num_t; k++) {
        s_CM = e_mass * e_mass + 2.0 * s.electron.energy[j] *
  s.target.energy[k]; beta = (s_CM - e_mass * e_mass) / (s_CM + e_mass *
  e_mass); r = s.electron.energy[j] / s.electron.energy[i]; if ((r >= 1) && (r
  <= (1 + beta) / (1 - beta))) { dum = 1 / s.electron.energy[j] * (1 + beta) /
  beta * (1. / r + r + 2 * (1 - beta) / beta * (1 - r) + (1 - beta) * (1 -
  beta) / beta / beta * (1 - r) * (1 - r));
          //ICE_table[i][j][k] = C * e_mass * e_mass / s_CM * dum;
          ICE_table[i][j][k] = C * e_mass * e_mass / s_CM /
  s.electron.energy[j] * (1 + beta) / beta; } else { ICE_table[i][j][k] = 0.0;
        }
      }
    }
  }
  */
}

void InverseCompton::SpecAIC(const double costheta,
                             const std::vector<double> &electron_energy,
                             const std::vector<double> &photon_energy,
                             const std::vector<double> &target_energy,
                             const std::vector<double> &target_spectrum) {
  // Scattering of directed photon beams on isotropically distributed
  // relativistic electrons
  int num_e = electron_energy.size();
  int num_p = photon_energy.size();
  int num_t = target_energy.size();
  spec_AIC.resize(num_p);
  for (size_t i = 0; i < num_p; i++) {
    spec_AIC[i].assign(num_e, 0.0);
  }
  double gamma;
  double table_AIC, z, z1, z2, z3, zbtheta, target_min;
  std::vector<double> target_temp(num_t, 0.0);
  for (size_t i = 0; i < num_p; i++) {
    for (size_t j = 0; j < num_e; j++) {
      if ((electron_energy[j] > 2 * e_mass) &&
          (photon_energy[i] < electron_energy[j])) {
        gamma = electron_energy[j] / e_mass;
        z = photon_energy[i] / electron_energy[j];
        for (size_t k = 0; k < num_t; k++) {
          target_min =
              photon_energy[i] / (2 * (1 - costheta) * gamma * gamma *
                                  (1 - photon_energy[i] / e_mass / gamma));
          if (target_energy[k] > target_min) {
            zbtheta = photon_energy[i] / gamma /
                      (2 * (1 - costheta) * target_energy[k] * gamma);
            z1 = z * z / (2.0 * (1 - z));
            z2 = -2 * zbtheta / (1 - z);
            z3 = z2 * z2 / 2.0;
            table_AIC =
                (3. / (16 * PI) * T_c_cnst / target_energy[k] / gamma / gamma) *
                (1 + z1 + z2 + z3);
            target_temp[k] = target_spectrum[k] * table_AIC;
          } else {
            target_temp[k] = 0.0;
          }
        }
        spec_AIC[i][j] =
            4 * PI * u.Integrate(target_energy, target_temp); // times by 4 pi
      } else {
        spec_AIC[i][j] = 0.0;
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

  // Pre-allocate output size (important before parallel access)
  secondary_spectrum.assign(num_s, 0.0);

  // Get raw pointer for the invariant input array to help the compiler
  const double *primary_spectrum_ptr = primary_spectrum.data();

// 1. Start Parallel Region
#pragma omp parallel
  {
    // 2. Thread-Local Allocation
    // Allocate ONCE per thread, rather than once per iteration 'i'.
    // This significantly reduces heap allocation overhead.
    std::vector<double> primary_temp(num_p);
    double *primary_temp_ptr = primary_temp.data();

// 3. Distribute the outer loop
#pragma omp for
    for (size_t i = 0; i < num_s; i++) {
      const std::vector<double> &spec_row = spec[i];
      const double *spec_row_ptr = spec_row.data();

// 4. Vectorize the inner math (SIMD)
#pragma omp simd
      for (size_t j = 0; j < num_p; j++) {
        primary_temp_ptr[j] = primary_spectrum_ptr[j] * spec_row_ptr[j];
      }

      // 5. Integrate
      // We pass the thread-local vector 'primary_temp' to the integrator
      secondary_spectrum[i] = u.Integrate(primary_momentum, primary_temp);
    }
  }
}

void InverseCompton::Emissivity(const std::vector<double> &target_energy,
                                const std::vector<double> &target_spectrum,
                                std::vector<double> &photon_energy,
                                std::vector<double> &photon_spectrum) {
  Spec(s.electron.energy, s.photon.energy, target_energy, target_spectrum);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, photon_energy,
             photon_spectrum);
}

void InverseCompton::Emissivity(const std::vector<double> &target_energy,
                                const std::vector<double> &target_spectrum,
                                const std::vector<double> &primary_energy,
                                const std::vector<double> &primary_spectrum,
                                std::vector<double> &photon_energy,
                                std::vector<double> &photon_spectrum) {
  Spec(primary_energy, photon_energy, target_energy, target_spectrum);
  Emissivity(spec, primary_energy, primary_spectrum, photon_energy,
             photon_spectrum);
}

void InverseCompton::Emissivity(const std::vector<double> &target_energy,
                                const std::vector<double> &target_spectrum) {
  int num_t = target_spectrum.size();
  Spec(s.electron.energy, s.photon.energy, target_energy, target_spectrum);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, s.photon.energy,
             s.photon.spectrum);
}

void InverseCompton::Emissivity() {
  Spec(s.electron.energy, s.photon.energy, s.target.energy, s.target.spectrum);
  Emissivity(spec, s.electron.momentum, s.electron.spectrum, s.photon.energy,
             s.photon.spectrum);
}

void InverseCompton::EmissivityAIC(const double costheta,
                                   const std::vector<double> &target_energy,
                                   const std::vector<double> &target_spectrum) {
  SpecAIC(costheta, s.electron.energy, s.photon.energy, target_energy,
          target_spectrum);
  Emissivity(spec_AIC, s.electron.momentum, s.electron.spectrum,
             s.photon.energy, s.photon.spectrum);
}

void InverseCompton::EmissivityAIC(const double cos_theta,
                                   const std::vector<double> &target_energy,
                                   const std::vector<double> &target_spectrum,
                                   std::vector<double> &photon_energy,
                                   std::vector<double> &photon_spectrum) {
  SpecAIC(cos_theta, s.electron.energy, s.photon.energy, target_energy,
          target_spectrum);
  Emissivity(spec_AIC, s.electron.momentum, s.electron.spectrum, photon_energy,
             photon_spectrum);
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
            FKN = 0.5 * (Gamma_e * Gamma_e * q * q) / (1 + Gamma_e * q) *
                  (1. - q);
            dum = (0.75 * T_c_cnst / s.target.energy[k] / gamma / gamma) *
                  (FT + FKN);
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

void InverseCompton::Table(const double particle_mass,
                           const std::vector<double> &primary_energy,
                           const std::vector<double> &secondary_energy,
                           const std::vector<double> &target_energy,
                           const std::vector<double> &target_spectrum) {
  int num_p = secondary_energy.size();
  int num_e = primary_energy.size();
  int num_t = target_energy.size();

  table_particle.resize(num_p);
  for (size_t i = 0; i < num_p; i++) {
    table_particle[i].resize(num_e);
    for (size_t j = 0; j < num_e; j++) {
      table_particle[i][j].resize(num_t, 0.);
    }
  }

  double gamma;
  double dum, q, Gamma_e, x, FT, FKN, target_min, target_max;
  double r = e_mass / particle_mass;
  for (size_t i = 0; i < num_p; i++) {
    for (size_t j = 0; j < num_e; j++) {
      gamma = primary_energy[j] / particle_mass;
      for (size_t k = 0; k < num_t; k++) {
        if (gamma >= 1) {
          x = secondary_energy[i] / particle_mass / gamma;
          Gamma_e = 4. * gamma * target_energy[k] / particle_mass;
          q = x / (Gamma_e * (1.0 - x));
          if ((q >= 1. / (4 * gamma * gamma)) && (q <= 1.)) {
            FT = 2. * q * log(q) + (1 + 2 * q) * (1. - q);
            FKN = 0.5 * (Gamma_e * Gamma_e * q * q) / (1 + Gamma_e * q) *
                  (1. - q);
            dum = (0.75 * T_c_cnst / target_energy[k] / gamma / gamma) *
                  (FT + FKN);
            table_particle[i][j][k] = dum * r * r;
          } else {
            table_particle[i][j][k] = 0.0;
          }
        } else {
          table_particle[i][j][k] = 0.0;
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
          mandels1 = 1. + 2. * s.electron.energy[k] * s.target.energy[i] *
                              (1 - beta * angular[j]) / (e_mass * e_mass);
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
          mandels1 = 1. + 2. * s.electron.energy[k] * s.target.energy[i] *
                              (1 - beta * angular[j]) / (e_mass * e_mass);
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
  std::ifstream gin1(std::string(DATA_PATH) +
                     "/data/EMProcess/mandelstam_table.dat");

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
    photon_losstime[i] =
        ComptonEffCrossSection(photon_energy[i]) * density * c_cnst / mue;
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
             (1 + k) * (2 * k * k - 2 * k - 1) / k / k / (1 + 2 * k) /
                 (1 + 2 * k) +
             4. * k * k / 3. / (1 + 2 * k) / (1 + 2 * k) / (1 + 2 * k) +
             ((1 + k) / k / k / k + 0.5 / k / k / k - 0.5 / k) *
                 log(1 + 2 * k)));
  } else {
    return 2. * PI * NucleusCharge * electron_radius * electron_radius *
           (4. / 3) * k;
  }
}

double InverseCompton::Compattenuation(double energy) {
  // return (6.022*1e23*pMass/NucleusMass)*ComptonEffCrossSection(a); //cm^2/g
  // screen region
  return (6.022 * 1e23) *
         ComptonEffCrossSection(energy); // cm^2/g screen region
}

void InverseCompton::Test() {
  std::cout << "Test : InverseCompton module         " << std::endl;
}
