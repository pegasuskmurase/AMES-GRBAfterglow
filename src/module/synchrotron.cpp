//=======================================================================================//
//AMES
// code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file synchrotron.hpp
// \brief Basic class for synchrotron process

#include "module/synchrotron.hpp"
#include "module/constants.hpp"
#include "module/source.hpp"
#include "module/utility.hpp"
#include <fstream>

// #include <boost/math/special_functions/bessel.hpp>
// #include <gsl/gsl_integration.h>

Synchrotron::Synchrotron(Source& s) : s(s) {}

void Synchrotron::Transport(std::vector<ParticleTransport>& p,
                            const double& mag_strength, bool &have_SSAAttenuation) {
    int idx_photon = -1;
    int idx_electron = -1;
    int idx_photon_syn = -1;
    int idx_photon_IC = -1;
    int idx_photon_pa = -1;
    int idx_photon_pm = -1;
    int idx_photon_pp = -1;
    int idx_electron_prim = -1;
    int idx_electron_pm = -1;
    int idx_electron_pa = -1;
    int idx_electron_pp = -1;
    int idx_proton = -1;
    int idx_pionminus = -1;
    int idx_pionplus = -1;
    int idx_muonminus_l=-1;
    int idx_muonminus_r=-1;
    int idx_muonplus_l=-1;
    int idx_muonplus_r=-1;


    for (size_t n = 0; n < p.size(); n++) {
        if (p[n].ID == 22)
            idx_photon = n;
        if (p[n].ID == 11)
            idx_electron = n;
        if (p[n].ID == 2201)
            idx_photon_syn = n;
        if (p[n].ID == 2202)
            idx_photon_IC = n;
        if (p[n].ID == 2203)
            idx_photon_pm = n;
        if (p[n].ID == 2204)
            idx_photon_pa = n;
        if (p[n].ID == 2205)
            idx_photon_pp = n;
        if (p[n].ID == 1100)
            idx_electron_prim = n;
        if (p[n].ID == 1101)
            idx_electron_pa = n;
        if (p[n].ID == 1102)
            idx_electron_pm = n;
        if (p[n].ID == 1103)
            idx_electron_pp = n;
        if (p[n].ID == 13)
            idx_proton = n;
        if (p[n].ID == 10007)
            idx_pionminus = n;
        if (p[n].ID == 10008)
            idx_pionplus = n;
        if (p[n].ID == 100041)
            idx_muonminus_l = n;
        if (p[n].ID == 100042)
            idx_muonminus_r = n;
        if (p[n].ID == 100051)
            idx_muonplus_l = n;
        if (p[n].ID == 100052)
            idx_muonplus_r = n;
    }

    if (idx_photon < 0 || (idx_electron < 0 && 
        idx_proton < 0 &&
        idx_pionminus < 0 && idx_pionplus < 0 &&
        idx_muonminus_l < 0 && idx_muonminus_r < 0 &&
        idx_muonplus_l < 0 && idx_muonplus_r < 0))
        {
            std::cout << "SYN Warning: Essential particles missing from simulation\n";
            return;
        } // No relevant particles    

    std::vector<double> particle_losstime;
    std::vector<double> particle_intetime;

// Parallelize per-particle work
#pragma omp parallel for private(particle_losstime)
    for (size_t pid = 0; pid < p.size(); pid++) {
        auto& particle = p[pid];

        if (std::abs(particle.Z) > 0) {
            particle_losstime.assign(particle.x.size(), 0.0);

            Losstime(mag_strength, std::abs(particle.Z), particle.mass,
                     particle.x_e, particle_losstime);

// Update A and B matrices in parallel
#pragma omp parallel for
            for (size_t i = 0; i < particle.x.size(); i++) {
                particle.A[i] +=
                    particle_losstime[i] * particle.x[i] / particle.dx[i];
                if (i + 1 < particle.x.size()) {
                    size_t j = i + 1;
                    particle.B[i][j] +=
                        particle_losstime[j] * particle.x[j] / particle.dx[j];
                }
            }

            // Local spec to avoid race conditions on p[idx_photon]
            std::vector<std::vector<double>> local_spec;
            Spec(mag_strength, std::abs(particle.Z), particle.mass,
                 particle.x_e, p[idx_photon].x, local_spec);
            
            // std::cout << "Synchrotron spec computed for particle ID " << particle.ID << std::endl;
            // bool all_zeros = true;
            // for (size_t i = 0; i < local_spec.size() && all_zeros; ++i) {
            //     for (size_t j = 0; j < local_spec[i].size() && all_zeros; ++j) {
            //         if (local_spec[i][j] >= 1.0e-30) {
            //             all_zeros = false;
            //             // Optional: print the first non-zero element found
            //             std::cout << "First non-zero at [" << i << "][" << j << "]: " 
            //                       << local_spec[i][j] << std::endl;
            //         }
            //                 }
            // }
            // if (all_zeros) {
            //     std::cout << "ALL elements are zero" << std::endl;
            // } else {
            //     std::cout << "NOT all elements are zero" << std::endl;
            // }
            
// Append results safely
#pragma omp critical
            {
                if (idx_photon >= 0) {
                p[idx_photon].C_id.push_back(pid);
                p[idx_photon].C.push_back(local_spec);
		        }
		
                if ((idx_photon_syn > 0) && (particle.ID == 1100)) {
                    p[idx_photon_syn].C_id.push_back(pid);
                    p[idx_photon_syn].C.push_back(local_spec);
                }

                if ((idx_photon_pm > 0) && (particle.ID == 1101)) {
                    p[idx_photon_pm].C_id.push_back(pid);
                    p[idx_photon_pm].C.push_back(local_spec);
                }

                if ((idx_photon_pa > 0) && (particle.ID == 1102)) {
                    p[idx_photon_pa].C_id.push_back(pid);
                    p[idx_photon_pa].C.push_back(local_spec);
                }

                if ((idx_photon_pp > 0) && (particle.ID == 1103)) {
                    p[idx_photon_pp].C_id.push_back(pid);
                    p[idx_photon_pp].C.push_back(local_spec);
                }
            }
        }
    }

    // Synchrotron self-absorption (serial, depends on all particles)
    if (have_SSAAttenuation)
    {
    if (idx_photon < 0 || idx_electron < 0)
        {
            std::cout << "SSA Warning: Essential particles missing from simulation\n";
            return;
        } // No relevant particles

    particle_intetime.assign(p[idx_photon].x.size(), 0.0);
    std::vector<double> electron_spectrum(p[idx_electron].x.size(), 0.0);

    for (size_t i = 0; i < p[idx_electron].x.size(); i++) {
        electron_spectrum[i] = p[idx_electron].u[i] / p[idx_electron].dx_e[i];
    }

    IntetimeSSA(mag_strength, p[idx_electron].x, p[idx_electron].x_e,
                electron_spectrum, p[idx_photon].x_e, particle_intetime);

    for (size_t i = 0; i < p[idx_photon].x.size(); i++) {
        p[idx_photon].A[i] += particle_intetime[i];
        if (idx_photon_syn > 0) {
            p[idx_photon_syn].A[i] += particle_intetime[i];
        }
        if (idx_photon_pm > 0) {
            p[idx_photon_pm].A[i] += particle_intetime[i];
        }
        if (idx_photon_pa > 0) {
            p[idx_photon_pa].A[i] += particle_intetime[i];
        }
        if (idx_photon_pp > 0) {
            p[idx_photon_pp].A[i] += particle_intetime[i];
        }
    }
    }
}
/*
void Synchrotron::Transport(std::vector<ParticleTransport> &p, const double
&mag_strength) {

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

  std::vector<double> particle_losstime;
  std::vector<double> particle_intetime;
  for (auto &particle : p) {

    int p_id = &particle - &p[0];

    if (std::abs(particle.Z) > 0) {

      particle_losstime.assign(particle.x.size(), 0.0);
      Losstime(mag_strength, std::abs(particle.Z), particle.mass, particle.x_e,
particle_losstime);

      #pragma omp parallel for
      for (size_t i = 0; i < particle.x.size(); i++) {
          // First operation (always valid)
          particle.A[i] += particle_losstime[i] * particle.x[i] /
particle.dx[i];

          // Second operation (only if i+1 is in range)
          if (i + 1 < particle.x.size()) {
              size_t j = i + 1;
              particle.B[i][j] += particle_losstime[j] * particle.x[j] /
particle.dx[j];
          }
      }

      Spec(mag_strength, std::abs(particle.Z), particle.mass, particle.x_e,
p[idx_photon].x, spec); p[idx_photon].C_id.push_back(p_id);
      p[idx_photon].C.push_back(spec);
    }
  }

  // Synchrotron self-absorption
  particle_intetime.assign(p[idx_photon].x.size(), 0.0);
  std::vector<double> electron_spectrum(p[idx_electron].x.size(), 0.0);
  for (size_t i = 0; i < p[idx_electron].x.size(); i++) {
    electron_spectrum[i] = p[idx_electron].u[i] / p[idx_electron].dx_e[i];
  }
  IntetimeSSA(mag_strength, p[idx_electron].x, p[idx_electron].x_e,
electron_spectrum, p[idx_photon].x_e, particle_intetime); for (size_t i = 0; i <
p[idx_photon].x.size(); i++) { p[idx_photon].A[i] += particle_intetime[i];
  }

}
*/

void Synchrotron::Transport(std::vector<ParticleTransportChangCooper>& p,
                            const double& mag_strength, bool &have_SSAAttenuation) {

    int idx_photon;

    for (auto& particle : p) {

        if (particle.ID == 22) {
            idx_photon = &particle - &p[0];
        }
    }

    std::vector<double> particle_losstime;
    std::vector<double> particle_losstime_ghost(3, 0.);
    for (auto& particle : p) {

        int p_id = &particle - &p[0];

        particle_losstime.assign(particle.x_e.size(), 0.);
        if (std::abs(particle.Z) > 0) {
            // if (particle.ID == 11) {

            Losstime(mag_strength, particle.Z, particle.mass, particle.x_e,
                     particle_losstime);

            for (size_t i = 0; i < particle.x.size(); i++) {
                particle.B[i] += particle_losstime[i] * particle.x[i];
                // particle.A_dt[i] += particle_losstime[i];
            }

            Losstime(mag_strength, particle.Z, particle.mass, particle.x_e,
                     particle_losstime_ghost);
            for (size_t i = 0; i < particle.xghost.size(); i++) {
                particle.Bghost[i] +=
                    particle_losstime_ghost[i] * particle.xghost[i];
            }

            Spec(mag_strength, std::abs(particle.Z), particle.mass,
                 particle.x_e, p[idx_photon].x, spec);
            p[idx_photon].C_id.push_back(p_id);
            p[idx_photon].C.push_back(spec);
        }

        // Synchrotron self-absorption
    }
}

// void Synchrotron::Losstime() {
//   Losstime(s.param.getMagStrength(), s.electron.energy, s.electron.losstime);
// }

void Synchrotron::Losstime(const double mag_strength) {
    Losstime(mag_strength, s.electron.energy, s.electron.losstime);
}

void Synchrotron::LosstimeProton(const double mag_strength) {
    Losstime(mag_strength, 1, proton_mass, s.proton.energy, s.proton.losstime);
}

void Synchrotron::Losstime(const double mag_strength,
                           const std::vector<double>& electron_energy,
                           std::vector<double>& electron_losstime) {
    double coeff = 4.0 / 3.0 * T_c_cnst; // coefficient
    double B = mag_strength;
    double U_B = B * B / (8.0 * PI) * erg2eV;
    double gamma;
    double beta2gamma2; // beta * beta * gamma * gamma
    int num_e = electron_energy.size();
    electron_losstime.resize(num_e);
#pragma omp parallel for
    for (size_t i = 0; i < num_e; i++) {
        gamma = electron_energy[i] / e_mass;
        if (gamma > 1) {
            beta2gamma2 = gamma * gamma - 1.;
            electron_losstime[i] =
                coeff * beta2gamma2 * U_B / (electron_energy[i] - e_mass);
        } else {
            electron_losstime[i] = 0;
        }
    }
}

void Synchrotron::Losstime(const double mag_strength, const int particle_charge,
                           const double particle_mass,
                           const std::vector<double>& particle_energy,
                           std::vector<double>& particle_losstime) {
    double coeff =
        4.0 / 3.0 * T_c_cnst * particle_charge * particle_charge *
        particle_charge * particle_charge *
        (e_mass * e_mass / particle_mass / particle_mass); // coefficient
    double B = mag_strength;
    double U_B = B * B / (8.0 * PI) * erg2eV;
    double gamma;
    double beta2gamma2; // beta * beta * gamma * gamma
    int num_p = particle_energy.size();
    particle_losstime.resize(num_p);
#pragma omp parallel for
    for (size_t i = 0; i < num_p; i++) {
        gamma = particle_energy[i] / particle_mass;
        if (gamma > 1) {
            beta2gamma2 = gamma * gamma - 1.;
            particle_losstime[i] = coeff * beta2gamma2 * U_B /
                                   (particle_energy[i] - particle_mass);
        } else {
            particle_losstime[i] = 0;
        }
    }
}

void Synchrotron::Spec(const double mag_strength,
                       const std::vector<double>& primary_energy,
                       const std::vector<double>& secondary_energy) {
    Spec(mag_strength, primary_energy, secondary_energy, spec);
}

void Synchrotron::Spec(const double mag_strength,
                       const std::vector<double>& primary_energy,
                       const std::vector<double>& secondary_energy,
                       std::vector<std::vector<double>>& spec) {
    
    // --- 1. Precompute Constants ---
    // Calculate these once, outside the loops.
    // (Using your original constants)
    const double coeff_critical_energy = 1.7365215856405376e-08 * mag_strength;
    const double coeff_emissivity = 64039.416604913444 * mag_strength;
    const double pi_term_sq = (3.62 / M_PI) * (3.62 / M_PI); // Standard M_PI is better
    
    int num_s = secondary_energy.size();
    int num_p = primary_energy.size();

    // --- 2. Memory Allocation ---
    // Parallelizing resize is usually not worth the overhead unless num_s is huge.
    // Serial allocation is safer and often just as fast for vectors of vectors.
    spec.resize(num_s);
    for (size_t i = 0; i < num_s; i++) {
        spec[i].resize(num_p);
    }

    // --- 3. Parallel Calculation ---
    // collapse(2): Flattens the i and j loops into one giant loop for better balancing.
    // schedule(dynamic, 10): '10' reduces scheduler overhead compared to default '1'.
#pragma omp parallel for collapse(2) schedule(dynamic, 10)
    for (size_t i = 0; i < num_s; i++) {
        for (size_t j = 0; j < num_p; j++) {
            
            // Optimization: Early check for primary energy threshold
            double gamma = primary_energy[j] / e_mass;
            if (gamma <= 2.0) {
                 spec[i][j] = 0.0;
                 continue; 
            }

            // Optimization: Check energy bounds before heavy math
            if (secondary_energy[i] >= primary_energy[j]) {
                spec[i][j] = 0.0;
                continue;
            }

            double gamma2beta2 = gamma * gamma - 1.0;
            double critical_energy = coeff_critical_energy * gamma2beta2;
            double x = secondary_energy[i] / critical_energy;

            // Approximation formula for synchrotron function R(x)
            // Note: exp((2.0/3.0)*log(1.0/x)) is x^(-2/3)
            double denominator = sqrt(std::pow(x, -0.6666666667) + pi_term_sq); 
            
            spec[i][j] = coeff_emissivity * std::exp(-x) / denominator / secondary_energy[i];

            // assert is okay, but compiling with -DNDEBUG disables it for speed
            assert(spec[i][j] >= 0.0 && "Negative emissivity detected");
        }
    }
}

void Synchrotron::Spec(const double mag_strength, const double primary_charge,
                       const double primary_mass,
                       const std::vector<double>& primary_energy,
                       const std::vector<double>& secondary_energy,
                       std::vector<std::vector<double>>& spec) {
    double coeff_critical_energy = 1.5 * hbar_cnst * primary_charge * e_charge *
                                   c_cnst / (primary_mass * e_cnst * e_cnst) *
                                   mag_strength;

    double coeff_emissivity =
        1.81 * sqrt(3.0) * primary_charge * primary_charge * primary_charge *
        e_charge * e_charge * e_charge /
        (hbar_cnst * 2.0 * PI * primary_mass * e_cnst) * mag_strength;

    int num_s = secondary_energy.size();
    int num_p = primary_energy.size();

    spec.resize(num_s);
#pragma omp parallel for
    for (int i = 0; i < num_s; i++) {
        spec[i].resize(num_p, 0.0);
    }

// Parallelize over both loops
#pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < num_s; i++) {
        for (int j = 0; j < num_p; j++) {
            double gamma = primary_energy[j] / primary_mass;
            double gamma2beta2 = gamma * gamma - 1.0;

            if ((secondary_energy[i] < primary_energy[j]) && (gamma > 2.0)) {
                double critical_energy = coeff_critical_energy * gamma2beta2;
                double x = secondary_energy[i] / critical_energy;

                spec[i][j] = coeff_emissivity * exp(-x) /
                             sqrt(exp((2.0 / 3.0) * log(1.0 / x)) +
                                  (3.62 / PI) * (3.62 / PI)) /
                             secondary_energy[i];

                assert(spec[i][j] >= 0.0 &&
                       "Negative particle_intetime detected in sync");
            } else {
                spec[i][j] = 0.0;
            }
        }
    }
}

void Synchrotron::SpecImplicit(const double mag_strength,
                               const double primary_charge,
                               const double primary_mass,
                               const std::vector<double>& primary_energy,
                               const std::vector<double>& secondary_energy,
                               const std::vector<double>& secondary_energy_bin,
                               std::vector<std::vector<double>>& spec) {
    double coeff_critical_energy = 1.5 * hbar_cnst * primary_charge * e_charge *
                                   c_cnst / (primary_mass * e_cnst * e_cnst) *
                                   mag_strength;

    double coeff_emissivity =
        1.81 * sqrt(3.0) * primary_charge * primary_charge * primary_charge *
        e_charge * e_charge * e_charge /
        (hbar_cnst * 2.0 * PI * primary_mass * e_cnst) * mag_strength;

    double critical_energy;
    double gamma, gamma2beta2, x;
    int num_s = secondary_energy.size();
    int num_p = primary_energy.size();
    spec.resize(num_s);
    for (size_t i = 0; i < num_s; i++) {
        spec[i].resize(num_p);
    }

    for (size_t i = 0; i < num_s; i++) {
        for (size_t j = 0; j < num_p; j++) {
            gamma = primary_energy[j] / primary_mass;
            gamma2beta2 = gamma * gamma - 1.0;
            if ((secondary_energy[i] < primary_energy[j]) && (gamma > 2.0)) {
                critical_energy = coeff_critical_energy * gamma2beta2;
                x = secondary_energy[i] / critical_energy;
                spec[i][j] = coeff_emissivity * exp(-x) /
                             sqrt(exp((2.0 / 3.0) * log(1.0 / x)) +
                                  (3.62 / PI) * (3.62 / PI)) /
                             secondary_energy[i] * secondary_energy_bin[i];
            } else {
                spec[i][j] = 0.0;
            }
        }
    }
}

void Synchrotron::SpecProton(const double mag_strength,
                             const std::vector<double>& primary_energy,
                             const std::vector<double>& secondary_energy) {
    double coeff_critical_energy = 1.5 * hbar_cnst * e_charge * c_cnst /
                                   (proton_mass * e_cnst * e_cnst) *
                                   mag_strength;
    double coeff_emissivity =
        1.81 * sqrt(3.0) * e_charge * e_charge * e_charge /
        (hbar_cnst * 2.0 * PI * proton_mass * e_cnst) * mag_strength;
    double critical_energy;
    double gamma, gamma2beta2, x;
    int num_s = secondary_energy.size();
    int num_p = primary_energy.size();
    spec.resize(num_s);
    for (size_t i = 0; i < num_s; i++) {
        spec[i].resize(num_p);
    }
    for (size_t i = 0; i < num_s; i++) {
        for (size_t j = 0; j < num_p; j++) {
            gamma = primary_energy[j] / proton_mass;
            gamma2beta2 = gamma * gamma - 1.0;
            if ((secondary_energy[i] < primary_energy[j]) && (gamma > 2.0)) {
                critical_energy = coeff_critical_energy * gamma2beta2;
                x = secondary_energy[i] / critical_energy;
                spec[i][j] = coeff_emissivity * exp(-x) /
                             sqrt(exp((2.0 / 3.0) * log(1.0 / x)) +
                                  (3.62 / PI) * (3.62 / PI)) /
                             secondary_energy[i];
            } else {
                spec[i][j] = 0.0;
            }
        }
    }
}

void Synchrotron::Emissivity(const std::vector<std::vector<double>>& spec,
                             const std::vector<double>& primary_momentum,
                             const std::vector<double>& primary_spectrum,
                             const std::vector<double>& secondary_energy,
                             std::vector<double>& secondary_spectrum) {
    int num_s = secondary_energy.size();
    int num_p = primary_momentum.size();

    // 1. Start Parallel Region
    // We create the parallel region here so we can allocate the buffer 
    // ONCE per thread, rather than re-allocating it inside every loop iteration.
#pragma omp parallel
    {
        // 2. Thread-Private Buffer
        // Each thread gets its own vector to store intermediate integration values.
        std::vector<double> primary_temp_private(num_p);

        // 3. Distribute the Outer Loop (Secondary Energy Bins)
#pragma omp for
        for (size_t i = 0; i < num_s; i++) {
            
            // 4. Vectorize Inner Loop (SIMD)
            // Use SIMD to speed up the element-wise multiplication
#pragma omp simd
            for (size_t j = 0; j < num_p; j++) {
                primary_temp_private[j] = primary_spectrum[j] * spec[i][j];
            }

            // 5. Integrate using the private buffer
            // We assume u.Integrate is thread-safe (read-only inputs)
            secondary_spectrum[i] = u.Integrate(primary_momentum, primary_temp_private);
        }
    }
}

void Synchrotron::Emissivity(double mag_strength,
                             const std::vector<double>& primary_energy,
                             const std::vector<double>& primary_spectrum,
                             const std::vector<double>& secondary_energy,
                             std::vector<double>& secondary_spectrum) {

    Spec(mag_strength, primary_energy, secondary_energy);
    Emissivity(spec, primary_energy, primary_spectrum, secondary_energy,
               secondary_spectrum);
}

void Synchrotron::Emissivity() {
    Spec(s.param.getMagStrength(), s.electron.energy, s.photon.energy);
    Emissivity(spec, s.electron.momentum, s.electron.spectrum, s.photon.energy,
               s.photon.spectrum);
}

void Synchrotron::Emissivity(double mag_strength) {
    Spec(mag_strength, s.electron.energy, s.photon.energy);
    Emissivity(spec, s.electron.momentum, s.electron.spectrum, s.photon.energy,
               s.photon.spectrum);
}

void Synchrotron::Emissivity(double mag_strength,
                             std::vector<double>& photon_energy,
                             std::vector<double>& photon_spectrum) {
    Spec(mag_strength, s.electron.energy, photon_energy);
    Emissivity(spec, s.electron.momentum, s.electron.spectrum, photon_energy,
               photon_spectrum);
}

void Synchrotron::EmissivityProton(double mag_strength) {
    SpecProton(mag_strength, s.proton.energy, s.photon.energy);
    Emissivity(spec, s.proton.momentum, s.proton.spectrum, s.photon.energy,
               s.photon.spectrum);
}

void Synchrotron::Attenuation(const std::vector<double>& photon_energy,
                              std::vector<double>& photon_spectrum) {
    double tau, att;
    int num_p = photon_spectrum.size();
    for (size_t i = 0; i < num_p; i++) {
        tau =
            u.Interpolate(s.photon.energy, s.photon.optdepth, photon_energy[i]);
        if (tau > 1e-5) {
            att = (1. - exp(-tau)) / tau;
        } else {
            att = 1. / (1 + tau);
        }
        photon_spectrum[i] *= att;
    }
}

void Synchrotron::Attenuation(std::vector<double>& photon_spectrum) {
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

void Synchrotron::SSAAbsorptionCoeff(const double mag_strength, 
                                     std::vector<double>& photon_syn_absorption_SSA) {
    int num_p = s.photon.momentum.size();
    int num_e = s.electron.momentum.size();

    // 1. Start Parallel Region
#pragma omp parallel
    {
        // 2. Thread-Private Buffer
        // Allocate once per thread to replace the shared 's.electron.temp'
        std::vector<double> temp_private(num_e);

        // 3. Distribute the Outer Loop (Photon Bins)
#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < num_p; i++) {
            
            double photon_E = s.photon.energy[i]; // Optimization: cache local value

            // 4. Vectorize Inner Loop (Electron Bins)
            // Note: SSACross likely contains heavy math (pow/bessel), so SIMD speedup 
            // depends on the compiler's ability to inline that function.
#pragma omp simd
            for (size_t j = 0; j < num_e; j++) {
                // Optimization: Pre-check energy threshold
                if (s.electron.energy[j] > e_mass) {
                    temp_private[j] = s.electron.spectrum[j] *
                                      SSACross(mag_strength, photon_E, s.electron.energy[j]);
                } else {
                    temp_private[j] = 0.0;
                }
            }

            // 5. Integrate using the private buffer
            // We assume u.Integrate is thread-safe (read-only inputs)
            photon_syn_absorption_SSA[i] = 
                u.Integrate(s.electron.momentum, temp_private);
        }
    }
}

void Synchrotron::IntetimeSSA(const double mag_strength,
                              const std::vector<double>& electron_momentum,
                              const std::vector<double>& electron_energy,
                              const std::vector<double>& electron_spectrum,
                              const std::vector<double>& photon_energy,
                              std::vector<double>& photon_intetime) {
    int num_p = photon_energy.size();
    int num_e = electron_momentum.size();
    for (size_t i = 0; i < num_p; i++) {
        for (size_t j = 0; j < num_e; j++) {
            if (s.electron.energy[j] > 1 * e_mass) {
                s.electron.temp[j] = electron_spectrum[j] *
                                     SSACross(mag_strength, photon_energy[i],
                                              electron_energy[j]);
            } else {
                s.electron.temp[j] = 0.0;
            }
        }
        photon_intetime[i] =
            c_cnst * u.Integrate(electron_momentum, s.electron.temp);
    }
}

void Synchrotron::SSAOpticalDepth(const double emission_size,
                                  const double mag_strength) {
    int num_p = s.photon.momentum.size();
    int num_e = s.electron.momentum.size();
    for (size_t i = 0; i < num_p; i++) {
        for (size_t j = 0; j < num_e; j++) {
            if (s.electron.energy[j] > 1 * e_mass) {
                s.electron.temp[j] = s.electron.spectrum[j] *
                                     SSACross(mag_strength, s.photon.energy[i],
                                              s.electron.energy[j]);
            } else {
                s.electron.temp[j] = 0.0;
            }
        }
        s.photon.optdepth[i] =
            emission_size * u.Integrate(s.electron.momentum, s.electron.temp);
    }
}

double Synchrotron::SSACross(const double mag_strength, double x, double y) {
    // From Eq 2.14 of the reference Ghisellini and Svensson 1991,
    // the pitch angle theta = pi / 2 is assumed
    y /= e_mass;
    double cycenergy;
    if ((x < y * e_mass) && (y > 1.)) {
        cycenergy = (1.5 * hbar_cnst * e_charge * mag_strength * y * y *
                     c_cnst / e_mass / e_cnst / e_cnst) *
                    sqrt(1 - 1. / y / y) * sqrt(1 - 1. / y / y); // eV unit
        return (2 * PI / sqrt(3.) / 4. / PI) * 137.036 * 6.652e-25 * 4.413e13 *
               mBessel53(x / cycenergy) / mag_strength / pow(y, 5);
    } else {
        return 0;
    }

    /*
    // From Eq 2.17 of the reference Ghisellini and Svensson, 1991
    // the pitch angle theta = pi / 2 is assumed
    y /= e_mass;
    double cycenergy;
    if ((x < y * e_mass) && (y > 1.)) {
      cycenergy =
          (1.5 * hbar_cnst * e_charge * mag_strength * y * y * c_cnst / e_mass /
    e_cnst / e_cnst) * sqrt(1 - 1. / y / y) * sqrt(1 - 1. / y / y); // eV unit
      return (sqrt(3.) * PI / 10.0) / 4.0 / PI * 137.036 * 6.652e-25 * 4.413e13
    * x / cycenergy * (boost::math::cyl_bessel_k(4.0 / 3.0, 0.5 * x / cycenergy)
    * boost::math::cyl_bessel_k(4.0 / 3.0, 0.5 * x / cycenergy) -
              boost::math::cyl_bessel_k(1.0 / 3.0, 0.5 * x / cycenergy) *
                  boost::math::cyl_bessel_k(1.0 / 3.0, 0.5 * x / cycenergy)) /
             mag_strength / pow(y, 5);
    } else {
      return 0;
    }
    */
}

double Synchrotron::mBessel53(double x) {
    double dummy1 = -1.0194198041210243 * x + 0.28011396300530672 * sqrt(x) -
                    7.71058491739234908 * 0.01 * pow(x, (1. / 3.));
    double dummy2 = -15.761577796582387 * x;
    if (x >= 0) {
        return 0.5 * 0.902745 * pow(0.5 * x, -(5. / 3.)) * exp(dummy1) +
               (1. - exp(dummy2)) * sqrt(0.5 * PI / x) * exp(-x);
    } else {
        return 0;
    }
}

void Synchrotron::Test() {
    std::cout << "Test : Synchrotron module         " << std::endl;

    // For Spec and Kappa
    /*
    std::vector<double> photon_energy = s.photon.energy;
    int num_p = photon_energy.size();

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    */

    /*
    std::vector<double> my_params(6);
    my_params[0] = mag_strength;
    my_params[1] = gamma_inj;
    my_params[2] = gamma_cut;
    my_params[3] = spectral_e;
    my_params[4] = 335.307;

    gsl_function F;
    //F.function = &myFunc::SynSpecGhiseLLini91;
    // F.function = &myFunc::SynKappaGhiseLLini91;
    // F.function = &myFunc::SynSpecInte;
    //F.function = &myFunc::SynKappaInte1;
    F.function = &myFunc::SynKappaInte2;
    F.params = &my_params;

    std::vector<double> photon_spectrum(num_p);
    std::vector<double> photon_optdepth(num_p);
    for (size_t i = 0; i < num_p; i++) {
      my_params[5] = photon_energy[i];
      gsl_integration_qags(&F, 1e2, 1e4, 0, 1e-8, 1000, w, &result, &error);
      //photon_spectrum[i] = result;
      photon_optdepth[i] = result * emission_radius;
    }
    gsl_integration_workspace_free(w);

    std::ofstream mystream;
    mystream.open("result4.dat");
    for (size_t i = 0; i < num_p; i++) {
      mystream << " " << photon_energy[i] << " " << s.photon.spectrum[i] << " "
    << photon_spectrum[i]
               << " " << s.photon.optdepth[i] << " " << photon_optdepth[i] << "
    " << std::endl;
    }
    mystream.close();
    */
}
