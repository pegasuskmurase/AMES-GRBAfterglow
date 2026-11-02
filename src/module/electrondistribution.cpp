//==============================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file electrondistribution.hpp
// \brief numerically determine electron distribution

#include "module/electrondistribution.hpp"
#include "module/constants.hpp"
#include "module/source.hpp"
#include "module/utility.hpp"
#include "module/source.hpp"
#include <fstream>

ElectronDistribution::ElectronDistribution(Source &s) : s(s) {}

void ElectronDistribution::IterationSolutionSteadyState(Synchrotron &syn, InverseCompton &IC,
                                                        bool have_EIC, double t_dyn,
                                                        double t_esc_photon, double mag_strength)
{

  // treament of steady-state solution, original method
  double ERR = 1e-4;
  const int N = s.electron.getSize();
  std::vector<double> x = s.electron.momentum;
  std::vector<double> x_e = s.electron.energy;
  std::vector<double> Q_inj = s.electron.spectrum; // momentum space
  std::vector<double> e_loss_Syn(N);
  std::vector<double> e_loss_EIC(N);
  std::vector<double> e_loss_SSC(N);
  std::vector<double> e_loss_tot(N);
  syn.Losstime(mag_strength, x_e, e_loss_Syn);

  if (have_EIC)
  {
    IC.Losstime(s.target.energy, s.target.spectrum, x_e, e_loss_EIC);
  }
  for (size_t j = 0; j < N; j++)
  {
    e_loss_tot[j] = e_loss_EIC[j] + e_loss_Syn[j];
  }

  double momentum_inj;
  for (size_t j = 0; j < N; j++)
  {
    if (Q_inj[j] > 0.0)
    {
      momentum_inj = s.electron.momentum[j];
      break;
    }
  }

  double dum, momentum_cooling;
  int idx;
  for (size_t j = 0; j < N; j++)
  {
    dum = e_loss_tot[j] * t_dyn;
    if (dum > 1.)
    {
      momentum_cooling = s.electron.momentum[j];
      idx = j;
      break;
    }
  }

  double momentum_eff;
  for (size_t j = 0; j < N; j++)
  {
    momentum_eff = momentum_inj * momentum_cooling / (momentum_inj + momentum_cooling);
    if (x[j] > momentum_eff)
    {
      dum = utility.Integrate(x, Q_inj, x[j]) / ((1. / t_dyn + e_loss_tot[j]) * x_e[j]);
    }
    else
    {
      dum = 0.0;
    }
    s.electron.setSpectrum(j, dum);
  }

  int num_t = s.target.energy.size();
  std::vector<double> target_syn(num_t);

  std::vector<double> spectrum_temp(N);
  int istep = 0;
  double u_err, num_tot, num_tot_inj;
  do
  {
    syn.Emissivity(mag_strength);
    for (size_t i = 0; i < num_t; i++)
    {
      target_syn[i] = s.photon.Interpolate(s.target.energy[i]) * t_esc_photon;
    }

    IC.Losstime(s.target.energy, target_syn, x_e, e_loss_SSC);
    u_err = 0.0;

    for (size_t j = 0; j < N; j++)
    {
      e_loss_tot[j] = e_loss_EIC[j] + e_loss_Syn[j] + e_loss_SSC[j];
    }

    for (size_t j = 0; j < N; j++)
    {
      dum = e_loss_tot[j] * t_dyn;
      if (dum > 1.)
      {
        momentum_cooling = s.electron.momentum[j];
        break;
      }
    }

    double dum;
    double momentum_eff;
    for (size_t j = 0; j < N; j++)
    {
      momentum_eff = momentum_inj * momentum_cooling / (momentum_inj + momentum_cooling);
      if (x[j] > momentum_eff)
      {
        dum = utility.Integrate(x, Q_inj, x[j]) / ((1. / t_dyn + e_loss_tot[j]) * x_e[j]);
      }
      else
      {
        dum = 0.0;
      }
      spectrum_temp[j] = dum;
      if (s.electron.spectrum[j] > 0)
      {
        u_err += std::abs(s.electron.spectrum[j] - spectrum_temp[j]) / s.electron.spectrum[j];
      }
    }
    for (size_t j = 0; j < N; j++)
    {
      s.electron.setSpectrum(j, spectrum_temp[j]);
    }
    istep++;
  } while (u_err > ERR);

  for (size_t i = 0; i < num_t; i++)
  {
    s.target.setSpectrum(i, target_syn[i]);
  }
}

void ElectronDistribution::IterationSolution(
    Synchrotron &syn, InverseCompton &IC, bool have_EIC,
    double t_dyn, double t_ad, double t_esc_photon, double mag_strength,
    const std::vector<double>& input_electron_Q, 
    const std::vector<double>& input_target_spectrum,
    std::vector<double>& output_electron_N,
    std::vector<double>& output_target_photon)
{
    // Treatment of quasi steady-state
    double ERR = 1e-5;
    
    // --- 1. Setup Read-Only Access ---
    const std::vector<double>& x = s.electron.momentum;
    const std::vector<double>& x_e = s.electron.energy;
    const std::vector<double>& target_energy = s.target.energy;
    const std::vector<double>& photon_energy = s.photon.energy;
    const int N = x.size();
    const int num_t = target_energy.size();

    // --- 2. Initialize Output Vectors ---
    if (output_electron_N.size() != N) output_electron_N.resize(N, 0.0);
    if (output_target_photon.size() != num_t) output_target_photon.resize(num_t);
    
    // --- 3. Initialize Local Buffers ---
    std::vector<double> e_loss_Syn(N);
    std::vector<double> e_loss_EIC(N);
    std::vector<double> e_loss_SSC(N, 0.0);
    std::vector<double> e_loss_tot(N);
    std::vector<double> photon_syn(photon_energy.size());
    std::vector<double> spectrum_temp(N);
    
    // --- 4. Calculate Static Losses (Synchrotron + EIC) ---
    // Note: Assuming 'syn.Losstime' and 'IC.Losstime' are internally parallelized.
    syn.Losstime(mag_strength, x_e, e_loss_Syn);

    if (have_EIC) {
        IC.Losstime(target_energy, input_target_spectrum, x_e, e_loss_EIC);
        
        #pragma omp parallel for simd
        for (int j = 0; j < N; j++) {
            e_loss_tot[j] = e_loss_EIC[j] + e_loss_Syn[j];
        }
    } else {
        #pragma omp parallel for simd
        for (int j = 0; j < N; j++) {
            e_loss_tot[j] = e_loss_Syn[j];
        }
    }

    // --- 5. Find Injection Break (gamma_inj) ---
    // Serial search is fast enough; parallelism overhead is not worth it here.
    double momentum_inj = x[0];
    for (int j = 0; j < N; j++) {
        if (input_electron_Q[j] > 0.0) {
            momentum_inj = x[j];
            break;
        }
    }

    // --- 6. Initial Cooling Break Estimation ---
    double momentum_cooling = x[N - 1];
    for (int j = 0; j < N; j++) {
        if (e_loss_tot[j] * t_dyn > 1.0) {
            momentum_cooling = x[j];
            break;
        }
    }

    // --- 7. Initial Analytical Guess ---
    double momentum_eff = momentum_inj * momentum_cooling / (momentum_inj + momentum_cooling);
    
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < N; j++) {
        if (x[j] > momentum_eff) {
            // Note: utility.Integrate must be thread-safe (read-only inputs)
            double integral = utility.Integrate(x, input_electron_Q, x[j]);
            output_electron_N[j] = integral / ((1.0 / t_ad + e_loss_tot[j]) * x_e[j]);
        } else {
            output_electron_N[j] = 0.0;
        }
    }

    // --- 8. Iterative Refinement ---
    int istep = 0;
    double u_err;

    do {
        // A. Calculate Synchrotron Emissivity (Internally parallelized)
        syn.Emissivity(mag_strength, x_e, output_electron_N, photon_energy, photon_syn);

        // B. Update Target Photon Field 
        #pragma omp parallel for
        for (int i = 0; i < num_t; i++) {
            // Target = Synchrotron Self-Absorption part
            output_target_photon[i] = 
                utility.Interpolate(photon_energy, photon_syn, target_energy[i]) * t_esc_photon;
        }

        // C. Calculate SSC Loss (Internally parallelized)
        IC.Losstime(target_energy, output_target_photon, x_e, e_loss_SSC);

        // D. Update Total Loss Rate
        if (have_EIC) {
            #pragma omp parallel for simd
            for (int j = 0; j < N; j++) {
                e_loss_tot[j] = e_loss_Syn[j] + e_loss_EIC[j] + e_loss_SSC[j];
            }
        } else {
            #pragma omp parallel for simd
            for (int j = 0; j < N; j++) {
                e_loss_tot[j] = e_loss_Syn[j] + e_loss_SSC[j];
            }
        }

        // E. Re-evaluate Cooling Break (Serial scan is safer/faster for 'break')
        momentum_cooling = x[N - 1]; 
        for (int j = 0; j < N; j++) {
            if (e_loss_tot[j] * t_dyn >= 1.0) {
                momentum_cooling = x[j];
                break;
            }
        }

        // F. Re-calculate Electron Distribution & Error
        momentum_eff = momentum_inj * momentum_cooling / (momentum_inj + momentum_cooling);
        u_err = 0.0;

        // CRITICAL: reduction(+:u_err) ensures thread-safe summation of the error
        #pragma omp parallel for reduction(+:u_err) schedule(static)
        for (int j = 0; j < N; j++) {
            if (x[j] > momentum_eff) {
                double integral = utility.Integrate(x, input_electron_Q, x[j]);
                
                // Solve continuity equation
                spectrum_temp[j] = integral / ((1.0 / t_dyn + 1.0 / t_ad + e_loss_tot[j]) * x_e[j]);
            } else {
                spectrum_temp[j] = 0.0;
            }

            // Calculate Error for convergence check
            // Use local variables to avoid reading shared memory multiple times
            double diff = 0.0;
            double current_val = output_electron_N[j];
            
            if (current_val > 1e-100) {
                diff = std::abs(current_val - spectrum_temp[j]) / current_val;
            } else if (spectrum_temp[j] > 1e-100) {
                diff = 1.0; 
            }
            u_err += diff;
        }

        // Update solution for next iteration
        // (Vector assignment is optimized by STL, usually uses memcpy)
        output_electron_N = spectrum_temp; 

        istep++;
        if (istep > 10) break;

    } while (u_err > ERR);
}

void ElectronDistribution::Test() {}

