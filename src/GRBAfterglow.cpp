// AMES
// Copyright(C) 2024, licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file GRBAfterglow.cpp
// \brief Calculate emissions GRB afterglow

#include "GRBAfterglow.hpp"
#include "module/rk5.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

GRBAfterglow::GRBAfterglow(Source& s)
    : s(s), syn(s.GetSynchrotron()), IC(s.GetIC()), gg(s.GetGG()),
      ED(s.GetED()), ph(s.GetPH()) {
    //    s.getTarget().setMomentum(321, 1e-8, 1e10);
    //    s.getPhoton().setMomentum(241, 1e-8, 1e16); // Note the photon
    //    spectrum is also used
    // for synchrotron target photon.
    //    s.getElectron().setMomentum(241, 1e4, 1e16);
}

void GRBAfterglow::help() {
    std::cout << "----Welcome to use AMES-GRBAfterglow------" << std::endl;
    std::cout << "Please see the following information for "
                 "proper use."
              << std::endl;
    std::cout << "If you have any questions please contact "
                 "Bing Theodore "
                 "Zhang<bing.zhang @yukawa.kyoto-u.ac.jp> "
                 "and Kohta "
                 "Murase<murase@psu.edu>."
              << std::endl;
}

void GRBAfterglow::setGRBAfterglowParam(std::vector<double> _param) {
    if ((_param[0] >= 0) && (_param[0] <= 10)) {
        param.z = _param[0];
    } else {
        param.z = 0.0;
        std::cout << "Warning: the redshift parameter [z] should within the "
                     "range (1, 10) "
                  << std::endl;
    }
    param.dl = _param[1];
    param.E_ej = _param[2];
    param.Gamma0 = _param[3];
    param.n_ex = _param[4];
    param.k_ex = _param[5];
    param.spectral_e = _param[6];
    param.epsilon_e = _param[7];
    param.fraction_e = _param[8];
    param.eta_acc_e = _param[9];
    param.epsilon_B = _param[10];
    param.open_angle = _param[11];
    param.view_angle = _param[12];
    param.gaussian_cone = _param[13];
    param.jet_index = _param[14];
    param.T_ej = _param[15]; // for thick shell

    // init
    num_t = s.getTarget().getEnergy().size();
    num_p = s.getPhoton().getEnergy().size();
    num_e = s.getElectron().getEnergy().size();

    target_energy = s.getTarget().getEnergy();
    target_syn.resize(num_t, 0.0);
    target_prompt.resize(num_t, 0.0);
    target_tot.resize(num_t, 0.0);

    photon_energy = s.getPhoton().getEnergy();
    photon_syn.resize(num_p, 0.0);
    photon_ssc.resize(num_p, 0.0);
    photon_eic.resize(num_p, 0.0);

    photon_syn_temp.resize(num_p, 0.0);
    photon_syn_absorption_SSA.resize(num_p, 0.0);
    photon_syn_obs.resize(num_p, 0.0);
    photon_ssc_temp.resize(num_p, 0.0);
    photon_ssc_obs.resize(num_p, 0.0);
    photon_eic_temp.resize(num_p, 0.0);
    photon_eic_obs.resize(num_p, 0.0);
    photon_tot.resize(num_p, 0.0);
    photon_tot_diff.resize(num_p, 0.0);
    photon_syn_diff.resize(num_p, 0.0);
    photon_ssc_diff.resize(num_p, 0.0);
    photon_eic_diff.resize(num_p, 0.0);

    photon_losstime.resize(num_p, 0.0);

    photon_syn_RS.resize(num_p, 0.0);
    photon_ssc_RS.resize(num_p, 0.0);
    photon_eic_RS.resize(num_p, 0.0);
    photon_psyn_RS.resize(num_p, 0.0);
    photon_syn_RS_obs.resize(num_p, 0.0);
    photon_ssc_RS_obs.resize(num_p, 0.0);
    photon_eic_RS_obs.resize(num_p, 0.0);
    photon_psyn_RS_obs.resize(num_p, 0.0);
}

void GRBAfterglow::setGRBAfterglowRSParam(std::vector<double> _param) {
    param.spectral_e_RS = _param[0]; //
    param.epsilon_e_RS = _param[1];  //
    param.fraction_e_RS = _param[2];
    param.eta_acc_e_RS = _param[3];
    param.epsilon_B_RS = _param[4]; //
}

void GRBAfterglow::setExtendedEmissionParam(std::vector<double> _param) {
    L_EE = _param[0];
    T_a_EE = _param[1];
    alpha_fl_EE = _param[2];
    alpha_st_EE = _param[3];
    energy_min_EE = _param[4];
    energy_peak_EE = _param[5];
    energy_max_EE = _param[6];
    spectra_l_EE = _param[7];
    spectra_h_EE = _param[8];
}

void GRBAfterglow::setPlateauEmissionParam(std::vector<double> _param) {
    L_PL = _param[0];
    T_a_PL = _param[1];
    alpha_fl_PL = _param[2];
    alpha_st_PL = _param[3];
    energy_min_PL = _param[4];
    energy_peak_PL = _param[5];
    energy_max_PL = _param[6];
    spectra_l_PL = _param[7];
    spectra_h_PL = _param[8];
}

void GRBAfterglow::setOutputFolder(std::string _output_folder) {
    output_folder = _output_folder;
}

void GRBAfterglow::haveOneZone(bool _have_onezone) {
    have_onezone = _have_onezone;
} // One zone
void GRBAfterglow::haveSSCSpec(bool _have_SSCSpec) {
    have_SSCSpec = _have_SSCSpec;
}
void GRBAfterglow::haveEICSpec(bool _have_EICSpec) {
    have_EICSpec = _have_EICSpec;
}
void GRBAfterglow::haveEdgeEffect(bool _have_edge_effect) {
    have_edge_effect = _have_edge_effect;
}
void GRBAfterglow::haveAttenuSSA(bool _have_attenu_SSA) {
    have_attenu_SSA = _have_attenu_SSA;
}
void GRBAfterglow::haveAttenuFF(bool _have_attenu_FF) {
    have_attenu_FF = _have_attenu_FF;
}
void GRBAfterglow::haveAttenuGGSource(bool _have_attenu_GG_source) {
    have_attenu_GG_source = _have_attenu_GG_source;
}
void GRBAfterglow::haveAttenuGGCosmic(bool _have_attenu_GG_cosmic) {
    have_attenu_GG_cosmic = _have_attenu_GG_cosmic;
}

void GRBAfterglow::haveReverseShock(bool _have_RS) {
    have_RS = _have_RS;
}; // have RS

void GRBAfterglow::Spectrum(std::vector<std::vector<double>>& flux_vector,
                            double T) {
    if (have_onezone) {
    } else {
        InitJet(jet);
    }

    flux_vector.resize(3);
    for (size_t i = 0; i < 3; i++) {
        flux_vector[i].assign(num_p, 0.0);
    }

    Spectrum(T);

    for (size_t j = 0; j < num_p; j++) {
        flux_vector[0][j] = photon_syn[j];
        flux_vector[1][j] = photon_ssc[j];
        flux_vector[2][j] = photon_eic[j];
        flux_vector[3][j] = photon_syn_RS[j];
        flux_vector[4][j] = photon_ssc_RS[j];
        flux_vector[5][j] = photon_eic_RS[j];
    }
}

void GRBAfterglow::Flux(std::vector<std::vector<double>>& flux_vector,
                        const std::vector<double>& time_array,
                        const std::vector<double>& energy_array_min,
                        const std::vector<double>& energy_array_max) {

    if (have_onezone) {
    } else {
        InitJet(jet);
    }

    flux_vector.resize(time_array.size());
    for (size_t i = 0; i < time_array.size(); i++) {
        flux_vector[i].assign(energy_array_min.size(), 0.0);
    }

    int idx = 0;
    for (auto t : time_array) {
        std::cout << " t " << t << std::endl;

        Spectrum(t);

        for (size_t i = 0; i < num_p; i++) {
            photon_tot[i] = photon_syn[i] + photon_ssc[i] + photon_eic[i];
            photon_tot_diff[i] = photon_tot[i] / photon_energy[i];
        }

        double dum;
        for (size_t i = 0; i < energy_array_min.size(); i++) {
            if (fabs(energy_array_min[i] - energy_array_max[i]) < 1e-5) {
                dum = utility.Interpolate(photon_energy, photon_tot,
                                          energy_array_min[i]);
            } else {
                dum =
                    utility.Integrate(photon_energy, photon_tot_diff,
                                      energy_array_min[i], energy_array_max[i]);
            }
            flux_vector[idx][i] = dum;
        }

        idx++;
    }
}

void GRBAfterglow::Flux(const std::vector<double>& time_array,
                        const std::vector<double>& energy_array_min,
                        const std::vector<double>& energy_array_max) {

    if (have_onezone) {
    } else {
        InitJet(jet);
    }

    std::ofstream mystream_spectrum;
    mystream_spectrum.open(output_folder + "/spectrum.dat");
    std::ofstream mystream_flux;
    mystream_flux.open(output_folder + "/flux.dat");
    mystream_flux << "# t " << std::endl;

    int idx = 0;
    for (auto t : time_array) {

        Spectrum(t);

        for (size_t i = 0; i < num_p; i++) {
            mystream_spectrum << photon_energy[i] << " ";
            mystream_spectrum << photon_syn[i] << " " << photon_ssc[i] << " "
                              << photon_eic[i] << " ";
            mystream_spectrum << std::endl;
            photon_tot[i] = photon_syn[i] + photon_ssc[i] + photon_eic[i];
            photon_tot_diff[i] = photon_tot[i] / photon_energy[i];
            photon_syn_diff[i] = photon_syn[i] / photon_energy[i];
            photon_ssc_diff[i] = photon_ssc[i] / photon_energy[i];
            photon_eic_diff[i] = photon_eic[i] / photon_energy[i];
        }

        mystream_flux << t << " ";
        for (size_t i = 0; i < energy_array_min.size(); i++) {
            if (fabs(energy_array_min[i] - energy_array_max[i]) < 1e-10) {
                mystream_flux << " "
                              << utility.Interpolate(photon_energy, photon_tot,
                                                     energy_array_min[i]);
                mystream_flux << " "
                              << utility.Interpolate(photon_energy, photon_syn,
                                                     energy_array_min[i]);
                mystream_flux << " "
                              << utility.Interpolate(photon_energy, photon_ssc,
                                                     energy_array_min[i]);
                mystream_flux << " "
                              << utility.Interpolate(photon_energy, photon_eic,
                                                     energy_array_min[i]);
            } else {
                mystream_flux
                    << " "
                    << utility.Integrate(photon_energy, photon_tot_diff,
                                         energy_array_min[i],
                                         energy_array_max[i]);
                mystream_flux
                    << " "
                    << utility.Integrate(photon_energy, photon_syn_diff,
                                         energy_array_min[i],
                                         energy_array_max[i]);
                mystream_flux
                    << " "
                    << utility.Integrate(photon_energy, photon_ssc_diff,
                                         energy_array_min[i],
                                         energy_array_max[i]);
                mystream_flux
                    << " "
                    << utility.Integrate(photon_energy, photon_eic_diff,
                                         energy_array_min[i],
                                         energy_array_max[i]);
            }
        }
        mystream_flux << std::endl;

        idx++;
    }

    mystream_spectrum.close();
    mystream_flux.close();
}

void GRBAfterglow::Spectrum(double T) {

    double shell_thickness_factor = 4. * (3. - param.k_ex);

    if (have_onezone) {
        double Gamma, radius;
        EvolutionThinShellAnalytic(T, Gamma, radius);
        double beta = sqrt(1. - 1. / Gamma / Gamma);
        double dynamical_timescale = Gamma * T / (1 + param.z);
        double adiabatic_timescale = dynamical_timescale;
        double Delta =
            radius / shell_thickness_factor / Gamma; // comoving frame
        double photon_escape_timescale = Delta / c_cnst;

        double mag_strength =
            sqrt(32 * PI * mp * c_cnst * c_cnst * param.epsilon_B *
                 ExternalDensity(radius) * Gamma * (Gamma - 1));
        s.getParam().setMagStrength(mag_strength);

        double electron_inj = (4. * Gamma) * param.fraction_e *
                              ExternalDensity(radius) / dynamical_timescale;

        s.getElectron().setSpectrumPL(MomentumInj(Gamma),
                                      MomentumMax(mag_strength, radius, beta),
                                      param.spectral_e, electron_inj, false);
        ED.IterationSolutionSteadyState(syn, IC, false, dynamical_timescale,
                                        photon_escape_timescale, mag_strength);
        target_syn = s.getTarget().getSpectrum();

        syn.Emissivity(mag_strength, photon_energy, photon_syn);

        if (have_SSCSpec) {
            IC.Emissivity(target_energy, target_syn, photon_energy, photon_ssc);
        }

        if (have_attenu_SSA) {
            syn.SSAOpticalDepth(dynamical_timescale * c_cnst, mag_strength);
            syn.Attenuation(photon_syn);
        }

        if (have_attenu_GG_source) {
            gg.Losstime(target_energy, target_syn, photon_energy,
                        photon_losstime);
            for (size_t i = 0; i < num_p; i++) {
                s.getPhoton().setOptdepth(i, photon_losstime[i] *
                                                 dynamical_timescale);
            }
            gg.Attenuation(photon_ssc);
        }

        LorentzBoost(photon_energy, photon_syn, Gamma);
        LorentzBoost(photon_energy, photon_ssc, Gamma);
        FluxNorm(photon_energy, photon_syn, radius, Delta);
        FluxNorm(photon_energy, photon_ssc, radius, Delta);

        if (have_edge_effect) {
            if (1. / Gamma > param.open_angle) {
                for (size_t k = 0; k < num_p; k++) {
                    photon_syn[k] *=
                        param.open_angle * param.open_angle * Gamma * Gamma;
                    photon_ssc[k] *=
                        param.open_angle * param.open_angle * Gamma * Gamma;
                }
            }
        }

    } else {

        if (param.T_ej > 0) {
            EvolutionThickShell(jet, T);
        } else {
            EvolutionThinShell(jet, T);
        }

        photon_syn.assign(num_p, 0.0);
        photon_ssc.assign(num_p, 0.0);
        photon_eic.assign(num_p, 0.0);

        double gamma_hat;
        double gamma_sh;
        for (size_t i = 0; i < jet.theta.size(); i++) {
            if (param.view_angle == 0) {
                int j = 0;

                beta = sqrt(1. -
                            1. / (jet.Gamma[i][j] *
                                  jet.Gamma[i][j])); // The velocity of shocked
                                                     // external matter

                gamma_hat = (4. * jet.Gamma[i][j] + 1) / 3. / jet.Gamma[i][j];

                gamma_sh =
                    ((jet.Gamma[i][j] + 1) *
                     (gamma_hat * (jet.Gamma[i][j] - 1) + 1) *
                     (gamma_hat * (jet.Gamma[i][j] - 1) + 1)) /
                    (gamma_hat * (2 - gamma_hat) * (jet.Gamma[i][j] - 1) + 2);
                gamma_sh = sqrt(gamma_sh);
                beta_sh = sqrt(1. - 1. / (gamma_sh * gamma_sh));

                doppler_factor =
                    1. / jet.Gamma[i][j] / (1 - beta * jet.mu[i][j]);
                delta_s = jet.radius[i][j] / shell_thickness_factor /
                          jet.Gamma[i][j] / jet.Gamma[i][j] /
                          std::abs(jet.mu[i][j] - beta_sh);

                EATS_factor = 2 * PI * jet.theta_bin[i] * sin(jet.theta[i]) *
                              jet.radius[i][j] * jet.radius[i][j] *
                              // 2 * PI * jet.theta_bin[i] * jet.radius[i][j] *
                              // jet.radius[i][j] *
                              std::abs(jet.mu[i][j] - beta_sh) /
                              (1 - beta_sh * jet.mu[i][j]) * (1 + param.z) /
                              param.dl / param.dl;

                /*
                std::cout << "On axis, Time: " << T
                          << ",  theta: " << jet.theta[i] << " ";
                std::cout << jet.radius[i][j] << " " << jet.Gamma[i][j] << " "
                          << gamma_sh << " ";
                std::cout << " " << delta_s << " " << EATS_factor << " ";
                std::cout << " " << T + jet.mu[i][j] * jet.radius[i][j] / c_cnst
                          << " ";
                std::cout << std::endl;
                */

                mag_strength =
                    sqrt(32 * PI * mp * c_cnst * c_cnst * param.epsilon_B *
                         ExternalDensity(jet.radius[i][j]) * jet.Gamma[i][j] *
                         (jet.Gamma[i][j] - 1));
                s.getParam().setMagStrength(mag_strength);

                // the escape timescale of photons in the shocked region
                photon_escape_timescale = jet.radius[i][j] /
                                          shell_thickness_factor /
                                          jet.Gamma[i][j] / c_cnst * 2;

                electron_inj = (4. * jet.Gamma[i][j]) * param.fraction_e *
                               ExternalDensity(jet.radius[i][j]) /
                               jet.dynamical_timescale[i][j];
                s.getElectron().setSpectrumPL(
                    MomentumInj(jet.Gamma[i][j]),
                    MomentumMax(mag_strength, jet.radius[i][j], beta),
                    param.spectral_e, electron_inj, false);

                if (have_EICSpec) {

                    // Generate prompt target photon density at given radius and
                    // theta
                    PromptPhotonDensity(ph, T, jet.radius[i][j],
                                        jet.Gamma[i][j]);

                    // ELectron cooling with External photons
                    ED.IterationSolution(syn, IC, have_EICSpec,
                                         jet.dynamical_timescale[i][j],
                                         jet.adiabatic_timescale[i][j],
                                         photon_escape_timescale, mag_strength);
                    target_syn = s.getTarget().getSpectrum();

                } else {

                    // ELectron cooling without External photons
                    ED.IterationSolution(syn, IC, false,
                                         jet.dynamical_timescale[i][j],
                                         jet.adiabatic_timescale[i][j],
                                         photon_escape_timescale, mag_strength);
                    target_syn = s.getTarget().getSpectrum();
                }

                for (size_t k = 0; k < num_t; k++) {
                    target_tot[k] = target_syn[k] + target_prompt[k];
                }

                syn.Emissivity(mag_strength, photon_energy, photon_syn_temp);
                for (size_t k = 0; k < num_p; k++) {
                    photon_syn_temp[k] *= photon_energy[k];
                }

                // SSC process
                if (have_SSCSpec) {
                    IC.Emissivity(target_energy, target_syn, photon_energy,
                                  photon_ssc_temp);
                    for (size_t k = 0; k < num_p; k++) {
                        photon_ssc_temp[k] *= photon_energy[k];
                    }
                }

                // EIC process
                if (have_EICSpec) {
                    double cos_theta_sc =
                        (jet.mu[i][j] - beta) / (1 - beta * jet.mu[i][j]);
                    IC.EmissivityAIC(cos_theta_sc, target_energy, target_prompt,
                                     photon_energy, photon_eic_temp);
                    // IC.Emissivity(target_energy, target_prompt,
                    // photon_energy, photon_eic_temp); // isotropic case
                    for (size_t k = 0; k < num_p; k++) {
                        photon_eic_temp[k] *= photon_energy[k];
                    }
                }

                if (have_attenu_SSA) {
                    syn.SSAAbsorptionCoeff(mag_strength,
                                           photon_syn_absorption_SSA);
                }

                for (size_t k = 0; k < num_p; k++) {
                    if (have_attenu_SSA) {
                        alpha = utility.Interpolate(
                                    photon_energy, photon_syn_absorption_SSA,
                                    photon_energy[k] * (1 + param.z) /
                                        doppler_factor) /
                                doppler_factor;
                        tau = alpha * delta_s;
                        if (tau > 1e-5) {
                            att = (1. - exp(-tau)) / tau;
                        } else {
                            att = 1. / (1 + tau);
                        }
                    } else {
                        att = 1.0;
                    }

                    photon_syn_obs[k] =
                        utility.Interpolate(photon_energy, photon_syn_temp,
                                            photon_energy[k] * (1 + param.z) /
                                                doppler_factor) *
                        doppler_factor * doppler_factor / (4 * PI) * att *
                        delta_s;
                }

                if (have_attenu_GG_source) {
                    gg.Losstime(target_energy, target_tot, photon_energy,
                                photon_losstime);
                }

                for (size_t k = 0; k < num_p; k++) {
                    if (have_attenu_GG_source) {
                        alpha = photon_losstime[k] / c_cnst;
                        alpha = utility.Interpolate(
                                    photon_energy, photon_losstime,
                                    photon_energy[k] * (1 + param.z) /
                                        doppler_factor) /
                                doppler_factor / c_cnst;
                        tau = alpha * delta_s;
                        if (tau > 1e-5) {
                            att = (1. - exp(-tau)) / tau;
                        } else {
                            att = 1. / (1 + tau);
                        }
                    } else {
                        att = 1.0;
                    }
                    photon_ssc_obs[k] =
                        utility.Interpolate(photon_energy, photon_ssc_temp,
                                            photon_energy[k] * (1 + param.z) /
                                                doppler_factor) *
                        doppler_factor * doppler_factor / (4 * PI) * att *
                        delta_s;
                    photon_eic_obs[k] =
                        utility.Interpolate(photon_energy, photon_eic_temp,
                                            photon_energy[k] * (1 + param.z) /
                                                doppler_factor) *
                        doppler_factor * doppler_factor / (4 * PI) * att *
                        delta_s;
                }

                for (size_t k = 0; k < num_p; k++) {
                    photon_syn[k] += EATS_factor * photon_energy[k] * eV2erg *
                                     photon_syn_obs[k];
                    photon_ssc[k] += EATS_factor * photon_energy[k] * eV2erg *
                                     photon_ssc_obs[k];
                    photon_eic[k] += EATS_factor * photon_energy[k] * eV2erg *
                                     photon_eic_obs[k];
                }

            } else {

                for (size_t j = 0; j < jet.phi.size(); j++) {

                    beta = sqrt(1. - 1. / (jet.Gamma[i][j] * jet.Gamma[i][j]));
                    gamma_hat =
                        (4. * jet.Gamma[i][j] + 1) / 3. / jet.Gamma[i][j];
                    gamma_sh =
                        ((jet.Gamma[i][j] + 1) *
                         (gamma_hat * (jet.Gamma[i][j] - 1) + 1) *
                         (gamma_hat * (jet.Gamma[i][j] - 1) + 1)) /
                        (gamma_hat * (2 - gamma_hat) * (jet.Gamma[i][j] - 1) +
                         2);
                    gamma_sh = sqrt(gamma_sh);

                    beta_sh = sqrt(1. - 1. / (gamma_sh * gamma_sh));
                    doppler_factor =
                        1. / jet.Gamma[i][j] / (1 - beta * jet.mu[i][j]);
                    delta_s = jet.radius[i][j] / shell_thickness_factor /
                              jet.Gamma[i][j] / jet.Gamma[i][j] /
                              fabs(jet.mu[i][j] - beta_sh);

                    EATS_factor = jet.theta_bin[i] * sin(jet.theta[i]) *
                                  jet.phi_bin[j] * 2 * jet.radius[i][j] *
                                  // EATS_factor = jet.theta_bin[i] *
                                  // jet.phi_bin[j] * 2 * jet.radius[i][j] *
                                  jet.radius[i][j] *
                                  fabs(jet.mu[i][j] - beta_sh) /
                                  (1 - beta_sh * jet.mu[i][j]) * (1 + param.z) /
                                  param.dl / param.dl;

                    std::cout << "Off axis, Time: " << T
                              << ", theta: " << jet.theta[i]
                              << ", phi: " << jet.phi[j] << " ";

                    std::cout << jet.radius[i][j] << " " << jet.Gamma[i][j]
                              << " " << gamma_sh << " ";
                    std::cout << std::endl;

                    mag_strength =
                        sqrt(32 * PI * mp * c_cnst * c_cnst * param.epsilon_B *
                             ExternalDensity(jet.radius[i][j]) *
                             jet.Gamma[i][j] * (jet.Gamma[i][j] - 1));
                    s.getParam().setMagStrength(mag_strength);

                    photon_escape_timescale = jet.radius[i][j] /
                                              shell_thickness_factor /
                                              jet.Gamma[i][j] / c_cnst * 2;
                    electron_inj = (4. * jet.Gamma[i][j]) * param.fraction_e *
                                   ExternalDensity(jet.radius[i][j]) /
                                   jet.dynamical_timescale[i][j];
                    s.getElectron().setSpectrumPL(
                        MomentumInj(jet.Gamma[i][j]),
                        MomentumMax(mag_strength, jet.radius[i][j], beta),
                        param.spectral_e, electron_inj, false);

                    if (have_EICSpec) {

                        // Generate prompt target photon density at given radius
                        // and theta
                        PromptPhotonDensity(ph, T, jet.radius[i][j],
                                            jet.Gamma[i][j]);

                        // ELectron cooling with External photons
                        ED.IterationSolution(syn, IC, have_EICSpec,
                                             jet.dynamical_timescale[i][j],
                                             jet.adiabatic_timescale[i][j],
                                             photon_escape_timescale);
                        target_syn = s.getTarget().getSpectrum();

                    } else {

                        // ELectron cooling without External photons
                        ED.IterationSolution(syn, IC, false,
                                             jet.dynamical_timescale[i][j],
                                             jet.adiabatic_timescale[i][j],
                                             photon_escape_timescale);
                        target_syn = s.getTarget().getSpectrum();
                    }

                    for (size_t k = 0; k < num_t; k++) {
                        target_tot[k] = target_syn[k] + target_prompt[k];
                    }

                    syn.Emissivity(mag_strength, photon_energy,
                                   photon_syn_temp);
                    for (size_t k = 0; k < num_p; k++) {
                        photon_syn_temp[k] *= photon_energy[k];
                    }

                    // SSC process
                    if (have_SSCSpec) {
                        IC.Emissivity(target_energy, target_syn, photon_energy,
                                      photon_ssc_temp);
                        for (size_t k = 0; k < num_p; k++) {
                            photon_ssc_temp[k] *= photon_energy[k];
                        }
                    }

                    // EIC process
                    if (have_EICSpec) {
                        double cos_theta_sc =
                            (jet.mu[i][j] - beta) / (1 - beta * jet.mu[i][j]);
                        IC.EmissivityAIC(cos_theta_sc, target_energy,
                                         target_prompt, photon_energy,
                                         photon_eic_temp);
                        // IC.Emissivity(target_energy, target_prompt,
                        // photon_energy, photon_eic_temp); // isotropic case
                        for (size_t k = 0; k < num_p; k++) {
                            photon_eic_temp[k] *= photon_energy[k];
                        }
                    }

                    if (have_attenu_SSA) {
                        syn.SSAAbsorptionCoeff(mag_strength,
                                               photon_syn_absorption_SSA);
                    }

                    for (size_t k = 0; k < num_p; k++) {
                        if (have_attenu_SSA) {
                            alpha =
                                utility.Interpolate(
                                    photon_energy, photon_syn_absorption_SSA,
                                    photon_energy[k] * (1 + param.z) /
                                        doppler_factor) /
                                doppler_factor;
                            tau = alpha * delta_s;
                            if (tau > 1e-5) {
                                att = (1. - exp(-tau)) / tau;
                            } else {
                                att = 1. / (1 + tau);
                            }
                        } else {
                            att = 1.0;
                        }

                        photon_syn_obs[k] =
                            utility.Interpolate(photon_energy, photon_syn_temp,
                                                photon_energy[k] *
                                                    (1 + param.z) /
                                                    doppler_factor) *
                            doppler_factor * doppler_factor / (4 * PI) * att *
                            delta_s;
                    }

                    if (have_attenu_GG_source) {
                        gg.Losstime(target_energy, target_tot, photon_energy,
                                    photon_losstime);
                    }

                    for (size_t k = 0; k < num_p; k++) {
                        if (have_attenu_GG_source) {
                            alpha = photon_losstime[k] / c_cnst;
                            alpha = utility.Interpolate(
                                        photon_energy, photon_losstime,
                                        photon_energy[k] * (1 + param.z) /
                                            doppler_factor) /
                                    doppler_factor / c_cnst;
                            tau = alpha * delta_s;
                            if (tau > 1e-5) {
                                att = (1. - exp(-tau)) / tau;
                            } else {
                                att = 1. / (1 + tau);
                            }
                        } else {
                            att = 1.0;
                        }
                        photon_ssc_obs[k] =
                            utility.Interpolate(photon_energy, photon_ssc_temp,
                                                photon_energy[k] *
                                                    (1 + param.z) /
                                                    doppler_factor) *
                            doppler_factor * doppler_factor / (4 * PI) * att *
                            delta_s;
                        photon_eic_obs[k] =
                            utility.Interpolate(photon_energy, photon_eic_temp,
                                                photon_energy[k] *
                                                    (1 + param.z) /
                                                    doppler_factor) *
                            doppler_factor * doppler_factor / (4 * PI) * att *
                            delta_s;
                    }

                    for (size_t k = 0; k < num_p; k++) {
                        photon_syn[k] += EATS_factor * photon_energy[k] *
                                         eV2erg * photon_syn_obs[k];
                        photon_ssc[k] += EATS_factor * photon_energy[k] *
                                         eV2erg * photon_ssc_obs[k];
                        photon_eic[k] += EATS_factor * photon_energy[k] *
                                         eV2erg * photon_eic_obs[k];
                    }
                }
            }
        }

        if (have_RS) {
            std::cout << "begin RS" << std::endl;
            ReverseShock(jet);
            std::cout << "finished RS" << std::endl;
        }
    }
}

void GRBAfterglow::ReverseShock(Jet& jet) {

    // emission from reverse shock
    double EATS_factor, delta_s, doppler_factor, Gamma34, mag_strength, n_ej,
        electron_inj, photon_escape_timescale, beta_sh, beta, beta_cross,
        ejecta_width_comoving, ejecta_width, gamma_hat;

    // structured jet
    double Delta0 = param.T_ej / (1 + param.z) * c_cnst; // lab frame
    double R_dec = DecelerationRadius();
    double g = 2;
    double r_index, r_index_B;
    if (Delta0 < R_dec / 2. / param.Gamma0 / param.Gamma0) {
        Delta0 = R_dec / 2. / param.Gamma0 / param.Gamma0;
        r_index = -6.0 * (3 + g) / 7.0;
        r_index_B = -4.0 * (3 + g) / 7.0;
    } else {
        r_index = -13.0 / 2.0;
        r_index_B = -13.0 / 3.0;
    }

    std::vector<double> electron_spectrum_cross(num_e, 0.0);
    double ejecta_width_comoving_cross;
    double mag_strength_cross;
    double mag_strength_postshock;
    double radius_cross;
    std::vector<double> proton_spectrum_cross(num_p, 0.0);

    photon_syn_RS.assign(num_p, 0.0);
    photon_ssc_RS.assign(num_p, 0.0);
    photon_eic_RS.assign(num_p, 0.0);
    double gamma_sh;
    for (size_t i = 0; i < jet.theta.size(); i++) {
        if (param.view_angle == 0) {

            int j = 0;

            beta = sqrt(1. - 1. / (jet.Gamma_RS[i][j] * jet.Gamma_RS[i][j]));

            gamma_hat = (4. * jet.Gamma_RS[i][j] + 1) / 3. / jet.Gamma_RS[i][j];

            gamma_sh =
                ((jet.Gamma_RS[i][j] + 1) *
                 (gamma_hat * (jet.Gamma_RS[i][j] - 1) + 1) *
                 (gamma_hat * (jet.Gamma_RS[i][j] - 1) + 1)) /
                (gamma_hat * (2 - gamma_hat) * (jet.Gamma_RS[i][j] - 1) + 2);
            gamma_sh = sqrt(gamma_sh);

            beta_sh = sqrt(1. - 1. / (gamma_sh * gamma_sh));

            doppler_factor =
                1. / jet.Gamma_RS[i][j] / (1 - beta * jet.mu[i][j]);

            EATS_factor = 2 * PI * jet.theta_bin[i] * sin(jet.theta[i]) *
                          jet.radius_RS[i][j] * jet.radius_RS[i][j] *
                          std::abs(jet.mu[i][j] - beta_sh) /
                          (1 - beta_sh * jet.mu[i][j]) * (1 + param.z) /
                          param.dl / param.dl;

            // total number of protons (and electrons) in the shell
            double N0_4 = jet.E_ej[i] /
                          (2 * PI * (1 - cos(param.open_angle)) / (4.0 * PI)) /
                          (jet.Gamma0[i] * mp * c_cnst * c_cnst);

            double epsilon_B_RS = param.epsilon_B_RS;

            if (jet.mass_ejecta[i][j] <= jet.M_ej[i]) {
                Gamma34 = 0.5 * (jet.Gamma_RS[i][j] / jet.Gamma0[i] +
                                 jet.Gamma0[i] / jet.Gamma_RS[i][j]);

                if (Gamma34 > 1.0001) {

                    gamma_hat = (4. + 1. / Gamma34) / 3.;

                    double N0_3 = N0_4 * (jet.mass_ejecta[i][j] / jet.M_ej[i]);

                    // n_4 is the density of the ejecta in the ejecta comoving
                    // frame
                    double n_4 =
                        N0_4 / (4 * PI * jet.Gamma0[i] * Delta0 *
                                jet.radius_RS[i][j] * jet.radius_RS[i][j]);

                    double n_3 = 4 * Gamma34 * n_4;

                    double e_3 = (Gamma34 - 1) * n_3 * mp * c_cnst * c_cnst;

                    double ejecta_width_comoving =
                        N0_3 / (4 * PI * jet.radius_RS[i][j] *
                                jet.radius_RS[i][j] * n_3);

                    double ejecta_width =
                        ejecta_width_comoving / jet.Gamma_RS[i][j];
                    delta_s = ejecta_width / std::abs(jet.mu[i][j] - beta_sh);

                    mag_strength = sqrt(8 * PI * epsilon_B_RS * e_3);
                    ED.setMagStrength(mag_strength);

                    photon_escape_timescale =
                        ejecta_width_comoving / c_cnst * 2;
                    electron_inj = param.fraction_e_RS * n_3 /
                                   jet.dynamical_timescale[i][j];

                    s.getElectron().setSpectrumPL(
                        MomentumInjRS(Gamma34),
                        MomentumMaxRS(mag_strength, jet.radius_RS[i][j], beta),
                        param.spectral_e_RS, electron_inj, false);

                    ED.IterationSolution(syn, IC, false,
                                         jet.dynamical_timescale[i][j],
                                         jet.adiabatic_timescale[i][j],
                                         photon_escape_timescale, mag_strength);

                    // Synchrotron emission
                    target_syn = s.getTarget().getSpectrum();

                    syn.Emissivity(mag_strength, photon_energy,
                                   photon_syn_temp);

                    for (size_t k = 0; k < num_p; k++) {
                        photon_syn_temp[k] *= photon_energy[k];
                    }

                    // Inverse Compton emission
                    if (have_SSCSpec) {
                        IC.Emissivity(target_energy, target_syn,
                         photon_energy, photon_ssc_temp);
                        for (size_t k = 0; k < num_p; k++) {
                            photon_ssc_temp[k] *= photon_energy[k];
                        }
                    }

                    // SSA Absorption
                    if (have_attenu_SSA) {
                        syn.SSAAbsorptionCoeff(mag_strength,
                                               photon_syn_absorption_SSA);
                    }

                    for (size_t k = 0; k < num_p; k++) {
                        if (have_attenu_SSA) {
                            alpha =
                                utility.Interpolate(
                                    photon_energy, photon_syn_absorption_SSA,
                                    photon_energy[k] * (1 + param.z) /
                                        doppler_factor) /
                                doppler_factor;
                            tau = alpha * delta_s;

                            if (tau > 1e-5) {
                                att = (1. - exp(-tau)) / tau;
                            } else {
                                att = 1. / (1 + tau);
                            }
                        } else {
                            att = 1.0;
                        }

                        photon_syn_RS_obs[k] =
                            utility.Interpolate(photon_energy, photon_syn_temp,
                                                photon_energy[k] *
                                                    (1 + param.z) /
                                                    doppler_factor) *
                            doppler_factor * doppler_factor / (4 * PI) * att *
                            delta_s;
                    }

                    // Gamma Gamma attenuation
                    if (have_attenu_GG_source) {
                        gg.Losstime(target_energy, target_syn, photon_energy,
                                    photon_losstime);
                    }

                    for (size_t k = 0; k < num_p; k++) {
                        if (have_attenu_GG_source) {
                            alpha = photon_losstime[k] / c_cnst;
                            alpha = utility.Interpolate(
                                        photon_energy, photon_losstime,
                                        photon_energy[k] * (1 + param.z) /
                                            doppler_factor) /
                                    doppler_factor / c_cnst;
                            tau = alpha * delta_s;
                            if (tau > 1e-5) {
                                att = (1. - exp(-tau)) / tau;
                            } else {
                                att = 1. / (1 + tau);
                            }
                        } else {
                            att = 1.0;
                        }

                        photon_ssc_RS_obs[k] =
                            utility.Interpolate(photon_energy, photon_ssc_temp,
                                                photon_energy[k] *
                                                    (1 + param.z) /
                                                    doppler_factor) *
                            doppler_factor * doppler_factor / (4 * PI) * att *
                            delta_s;
                    }
                } else {
                    for (size_t k = 0; k < num_p; k++) {
                        photon_syn_RS_obs[k] = 0.0;
                        photon_ssc_RS_obs[k] = 0.0;
                        photon_eic_RS_obs[k] = 0.0;
                    }
                }
            } else {
                beta =
                    sqrt(1. - 1. / (jet.Gamma_RS[i][j] * jet.Gamma_RS[i][j]));

                gamma_hat =
                    (4. * jet.Gamma_RS[i][j] + 1) / 3. / jet.Gamma_RS[i][j];

                gamma_sh =
                    ((jet.Gamma_RS[i][j] + 1) *
                     (gamma_hat * (jet.Gamma_RS[i][j] - 1) + 1) *
                     (gamma_hat * (jet.Gamma_RS[i][j] - 1) + 1)) /
                    (gamma_hat * (2 - gamma_hat) * (jet.Gamma_RS[i][j] - 1) +
                     2);
                gamma_sh = sqrt(gamma_sh);

                beta_sh = sqrt(1. - 1. / (gamma_sh * gamma_sh));

                doppler_factor =
                    1. / jet.Gamma_RS[i][j] / (1 - beta * jet.mu[i][j]);

                EATS_factor = 2 * PI * jet.theta_bin[i] * sin(jet.theta[i]) *
                              jet.radius_RS[i][j] * jet.radius_RS[i][j] *
                              std::abs(jet.mu[i][j] - beta_sh) /
                              (1 - beta_sh * jet.mu[i][j]) * (1 + param.z) /
                              param.dl / param.dl;

                double beta_cross = sqrt(
                    1. - 1. / (jet.Gamma_cross[i][j] * jet.Gamma_cross[i][j]));
                double Gamma34_cross =
                    0.5 * (jet.Gamma_cross[i][j] / jet.Gamma0[i] +
                           jet.Gamma0[i] / jet.Gamma_cross[i][j]);

                // n_4 is the density of the ejecta in the ejecta comoving frame
                double n_4_cross =
                    N0_4 / (4 * PI * jet.Gamma0[i] * Delta0 *
                            jet.radius_cross[i][j] * jet.radius_cross[i][j]);

                double n_3_cross = 4 * Gamma34_cross * n_4_cross;

                double e_3_cross =
                    (Gamma34_cross - 1) * n_3_cross * mp * c_cnst * c_cnst;

                double ejecta_width_comoving_crossing =
                    N0_4 / (4 * PI * jet.radius_cross[i][j] *
                            jet.radius_cross[i][j] * n_3_cross);

                double ejecta_width_comoving =
                    ejecta_width_comoving_crossing *
                    pow(jet.radius_RS[i][j] / jet.radius_cross[i][j],
                        -r_index - 2);

                double ejecta_width =
                    ejecta_width_comoving / jet.Gamma_RS[i][j];
                delta_s = ejecta_width / std::abs(jet.mu[i][j] - beta_sh);

                mag_strength_cross = sqrt(8 * PI * epsilon_B_RS * e_3_cross);
                mag_strength_postshock =
                    mag_strength_cross *
                    pow(jet.radius_RS[i][j] / jet.radius_cross[i][j],
                        r_index_B);
                ED.setMagStrength(mag_strength_postshock);

                double electron_inj_cross = param.fraction_e_RS * n_3_cross /
                                            jet.dynamical_timescale[i][j];
                s.getElectron().setSpectrumPL(
                    MomentumInjRS(Gamma34_cross),
                    MomentumMaxRS(mag_strength_cross, jet.radius_cross[i][j],
                                  sqrt(1. - 1. / (jet.Gamma_cross[i][j] *
                                                  jet.Gamma_cross[i][j]))),
                    param.spectral_e_RS, electron_inj_cross, false);

                double photon_escape_timescale_crossing =
                    ejecta_width_comoving_crossing / c_cnst * 2;
                photon_escape_timescale = ejecta_width_comoving / c_cnst * 2;
                ED.IterationSolution(
                    syn, IC, false, jet.dynamical_timescale_cross[i][j],
                    jet.adiabatic_timescale_cross[i][j],
                    photon_escape_timescale_crossing, mag_strength);

                double dum;
                std::vector<double> electron_momentum =
                    s.getElectron().getMomentum();
                std::vector<double> electron_momentum_after =
                    s.getElectron().getMomentum();
                std::vector<double> electron_spectrum_after =
                    s.getElectron().getSpectrum();
                for (int k = 0; k < num_e; k++) {
                    electron_momentum_after[k] =
                        electron_momentum[k] *
                        pow(jet.radius_RS[i][j] / jet.radius_cross[i][j], -1);
                }

                for (int k = 0; k < num_e; k++) {
                    dum = utility.Interpolate(electron_momentum_after,
                                              electron_spectrum_after,
                                              electron_momentum[k]);
                    s.getElectron().setSpectrum(k, dum);
                }

                for (int k = 0; k < num_e; k++) {
                    dum = s.getElectron().getSpectrum()[k];
                    dum *= pow(jet.radius_RS[i][j] / jet.radius_cross[i][j],
                               r_index);
                    s.getElectron().setSpectrum(k, dum);
                }

                syn.Emissivity(mag_strength_postshock, photon_energy,
                               photon_syn_temp);

                for (int k = 0; k < num_t; k++) {
                    target_syn[k] =
                        utility.Interpolate(photon_energy, photon_syn_temp,
                                            target_energy[k]) *
                        photon_escape_timescale;
                }

                for (size_t k = 0; k < num_p; k++) {
                    photon_syn_temp[k] *= photon_energy[k];
                }

                if (have_SSCSpec) {
                     IC.Emissivity(target_energy, target_syn, photon_energy,
                     photon_ssc_temp);
                    for (size_t k = 0; k < num_p; k++) {
                        photon_ssc_temp[k] *= photon_energy[k];
                    }
                }

                // SSA Absorption
                if (have_attenu_SSA) {
                    syn.SSAAbsorptionCoeff(mag_strength_postshock,
                                           photon_syn_absorption_SSA);
                }

                for (size_t k = 0; k < num_p; k++) {
                    if (have_attenu_SSA) {
                        alpha = utility.Interpolate(
                                    photon_energy, photon_syn_absorption_SSA,
                                    photon_energy[k] * (1 + param.z) /
                                        doppler_factor) /
                                doppler_factor;
                        tau = alpha * delta_s;
                        if (tau > 1e-5) {
                            att = (1. - exp(-tau)) / tau;
                        } else {
                            att = 1. / (1 + tau);
                        }
                    } else {
                        att = 1.0;
                    }

                    photon_syn_RS_obs[k] =
                        utility.Interpolate(photon_energy, photon_syn_temp,
                                            photon_energy[k] * (1 + param.z) /
                                                doppler_factor) *
                        doppler_factor * doppler_factor / (4 * PI) * att *
                        delta_s;
                }

                // Gamma Gamma attenuation
                if (have_attenu_GG_source) {
                    gg.Losstime(target_energy, target_syn, photon_energy,
                                photon_losstime);
                }

                for (size_t k = 0; k < num_p; k++) {
                    if (have_attenu_GG_source) {
                        alpha = photon_losstime[k] / c_cnst;
                        alpha = utility.Interpolate(
                                    photon_energy, photon_losstime,
                                    photon_energy[k] * (1 + param.z) /
                                        doppler_factor) /
                                doppler_factor / c_cnst;
                        tau = alpha * delta_s;
                        if (tau > 1e-5) {
                            att = (1. - exp(-tau)) / tau;
                        } else {
                            att = 1. / (1 + tau);
                        }
                    } else {
                        att = 1.0;
                    }

                    photon_ssc_RS_obs[k] =
                        utility.Interpolate(photon_energy, photon_ssc_temp,
                                            photon_energy[k] * (1 + param.z) /
                                                doppler_factor) *
                        doppler_factor * doppler_factor / (4 * PI) * att *
                        delta_s;
                }
            }

            for (size_t k = 0; k < num_p; k++) {
                photon_syn_RS[k] += EATS_factor * photon_energy[k] *
                                    photon_syn_RS_obs[k] * eV2erg;
                photon_ssc_RS[k] += EATS_factor * photon_energy[k] *
                                    photon_ssc_RS_obs[k] * eV2erg;
                photon_eic_RS[k] += EATS_factor * photon_energy[k] *
                                    photon_eic_RS_obs[k] * eV2erg;
            }
        } else {
            // TODO
            // For off-axis
        }
        std::cout << "spectrum RS finish " << std::endl;
    }
}

double GRBAfterglow::MomentumInj(double Gamma) {
    double energy_kin; // kinetic energy
    if (param.spectral_e > 2) {
        energy_kin = proton_mass * param.epsilon_e / param.fraction_e *
                     (param.spectral_e - 2) / (param.spectral_e - 1) *
                     (Gamma - 1);
    } else if (param.spectral_e < 2) {
        energy_kin = proton_mass * param.epsilon_e / param.fraction_e *
                     (2.0 - param.spectral_e) / (param.spectral_e - 1) *
                     (Gamma - 1);
    } else {
        // TDDO
        // gamma_max
        // energy_kin = proton_mass * param.epsilon_e / param.fraction_e /
        // log(gamma_max / gamma_min)
        //   * (Gamma - 1);
        // std::cout << "Warnning: spectral index should larger than 2 !!! \n"
        // << std::endl; return 0;
    }
    return s.getElectron().Energy2Momentum(energy_kin + e_mass);
}

double GRBAfterglow::MomentumMax(double mag_strength, double radius,
                                 double beta) {
    // t_syn  == t_acc
    double gamma_max_syn =
        sqrt(6 * PI * e_charge / param.eta_acc_e / sigmaT / mag_strength);

    // t_ad  == t_acc
    double gamma_max_ad = e_charge * mag_strength * radius /
                          (param.eta_acc_e * me * c_cnst * c_cnst * beta);

    return s.getElectron().Energy2Momentum(
        e_mass * std::min(gamma_max_syn, gamma_max_ad));
}

double GRBAfterglow::ExternalDensity(double radius) {
    if (param.k_ex > 0) {
        return param.n_ex * pow(radius, -param.k_ex);
    } else {
        return param.n_ex;
    }
}

void GRBAfterglow::PromptPhotonDensity(Photonbackground& ph, double T,
                                       double radius, double Gamma) {
    double u_EE = ExtendedEmissionLuminosity(T) /
                  (4 * PI * radius * radius * Gamma * Gamma * c_cnst) * erg2eV;
    ph.BrokenPowerlaw(energy_min_EE / Gamma, energy_peak_EE / Gamma,
                      energy_max_EE / Gamma, spectra_l_EE, spectra_h_EE, u_EE);
    std::vector<double> target_EE = s.getTarget().getSpectrum();
    std::cout << " " << u_EE << " " << ExtendedEmissionLuminosity(T) << " "
              << radius << " " << Gamma << std::endl;

    double u_PL = PlateauEmissionLuminosity(T) /
                  (4 * PI * radius * radius * Gamma * Gamma * c_cnst) * erg2eV;
    ph.BrokenPowerlaw(energy_min_PL / Gamma, energy_peak_PL / Gamma,
                      energy_max_PL / Gamma, spectra_l_PL, spectra_h_PL, u_PL);
    std::vector<double> target_PL = s.getTarget().getSpectrum();

    for (size_t k = 0; k < num_t; k++) {
        s.getTarget().setSpectrum(k, target_EE[k] + target_PL[k]);
        target_prompt[k] = target_EE[k] + target_PL[k];
    }
}

double GRBAfterglow::ExtendedEmissionLuminosity(double T) {
    double dum = 0.0;
    if (T < T_a_EE) {
        dum = L_EE * pow(T / T_a_EE, -alpha_fl_EE);
    } else {
        dum = L_EE * pow(T / T_a_EE, -alpha_st_EE);
    }
    return dum;
}

double GRBAfterglow::PlateauEmissionLuminosity(double T) {
    double dum = 0.0;
    if (T < T_a_PL) {
        dum = L_PL * pow(T / T_a_PL, -alpha_fl_PL);
    } else {
        dum = L_PL * pow(T / T_a_PL, -alpha_st_PL);
    }
    return dum;
}

void GRBAfterglow::InitJet(Jet& jet) {
    double sin_v = sin(param.view_angle);
    double cos_v = cos(param.view_angle);
    double sin_theta;
    double theta_min = 0.001 / param.Gamma0;
    double theta_max = param.open_angle;

    jet.theta.resize(jet.angular_num);
    jet.theta_bin.resize(jet.angular_num);

    /*
    for (size_t i = 0; i < jet.angular_num; i++) {
      jet.theta[i] = 0. + param.open_angle * i / (jet.angular_num - 1);
      jet.theta_bin[i] = param.open_angle / jet.angular_num;
    }
    */
    if (param.view_angle > 0) {
        for (size_t i = 0; i < jet.angular_num; i++) {
            jet.theta[i] = 0. + param.open_angle * i / (jet.angular_num - 1);
            jet.theta_bin[i] = param.open_angle / jet.angular_num;
            // jet.theta[i] = theta_min + (1 - theta_min) * i / (jet.angular_num
            // - 1); jet.theta_bin[i] = (1 - theta_min) / jet.angular_num;
        }
    } else {

        // Calculate the logarithmic spacing factor
        double log_factor = (std::log10(theta_max) - std::log10(theta_min)) /
                            (jet.angular_num - 1);

        // Generate the logarithmically spaced grid

        for (int i = 0; i < jet.angular_num; ++i) {
            jet.theta[i] = std::pow(10, std::log10(theta_min) + i * log_factor);
        }
        double theta_l =
            std::pow(10, std::log10(theta_min) + (-1) * log_factor);
        double theta_r =
            std::pow(10, std::log10(theta_min) + jet.angular_num * log_factor);
        for (int i = 0; i < jet.angular_num; ++i) {
            if (i == 0) {
                jet.theta_bin[i] = 0.5 * (jet.theta[i + 1] - theta_l);
            } else if (i == jet.angular_num - 1) {
                jet.theta_bin[i] = 0.5 * (theta_r - jet.theta[i - 1]);
            } else {
                jet.theta_bin[i] = 0.5 * (jet.theta[i + 1] - jet.theta[i - 1]);
            }
        }
    }

    jet.phi.resize(jet.phi_num);
    jet.phi_bin.resize(jet.phi_num);
    for (size_t i = 0; i < jet.phi_num; i++) {
        jet.phi[i] = 0. + 2 * PI * i / (jet.phi_num - 1);
        jet.phi_bin[i] = 2 * PI / jet.phi_num;
    }

    jet.Gamma0.resize(jet.angular_num);
    jet.M_ej.resize(jet.angular_num);
    jet.E_ej.resize(jet.angular_num);
    jet.Gamma.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.radius.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.Gamma_RS.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.radius_RS.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.mass_ejecta.resize(jet.angular_num,
                           std::vector<double>(jet.phi_num, 0.));
    jet.radius_cross.resize(jet.angular_num,
                            std::vector<double>(jet.phi_num, 0.));
    jet.Gamma_cross.resize(jet.angular_num,
                           std::vector<double>(jet.phi_num, 0.));
    jet.mass_cross.resize(jet.angular_num,
                          std::vector<double>(jet.phi_num, 0.));
    jet.mu.resize(jet.angular_num, std::vector<double>(jet.phi_num, 0.));
    jet.dynamical_timescale.resize(jet.angular_num,
                                   std::vector<double>(jet.phi_num, 0.));
    jet.adiabatic_timescale.resize(jet.angular_num,
                                   std::vector<double>(jet.phi_num, 0.));
    jet.dynamical_timescale_cross.resize(jet.angular_num,
                                         std::vector<double>(jet.phi_num, 0.));
    jet.adiabatic_timescale_cross.resize(jet.angular_num,
                                         std::vector<double>(jet.phi_num, 0.));

    for (size_t i = 0; i < jet.theta.size(); i++) {
        for (size_t j = 0; j < jet.phi.size(); j++) {
            jet.mu[i][j] =
                sin(jet.theta[i]) * sin(param.view_angle) * cos(jet.phi[j]) +
                cos(jet.theta[i]) * cos(param.view_angle);
        }
    }
}

void GRBAfterglow::calcDynamic(const std::vector<double>& time_array) {
    const int n = 7;
    std::vector<std::vector<double>> shock_dynamic;
    // double t_max = 1e15;

    std::ofstream mystream;
    mystream.open(output_folder + "/dynamic.dat");
    EvolutionThinShell(shock_dynamic, time_array);
    for (size_t i = 0; i < time_array.size(); i++) {
        double Gamma = shock_dynamic[i][0];
        double beta = sqrt(1. - 1. / (Gamma * Gamma));
        double gamma_hat = (4. * Gamma + 1) / 3. / Gamma;
        double gamma2_sh = ((Gamma + 1) * (gamma_hat * (Gamma - 1) + 1) *
                            (gamma_hat * (Gamma - 1) + 1)) /
                           (gamma_hat * (2 - gamma_hat) * (Gamma - 1) + 2);
        double gamma_sh = sqrt(gamma2_sh);
        double beta_sh = sqrt(1. - 1. / (gamma_sh * gamma_sh));

        mystream << std::setprecision(16);
        mystream << time_array[i] << " " << shock_dynamic[i][0] << " "
                 << shock_dynamic[i][1] << " " << shock_dynamic[i][2] << " "
                 << Gamma * beta << " " << gamma_sh * beta_sh << std::endl;
    }
    mystream.close();
}

void GRBAfterglow::EvolutionThinShell(
    std::vector<std::vector<double>>& shock_dynamic,
    const std::vector<double>& time_array) {
    const int n = 7;
    std::vector<double> y(n);
    std::vector<double> dydx(n);
    std::vector<double> yout(n);

    double epsilon = 0.0;
    double C_BM = epsilon + (9.0 - 2.0 * param.k_ex) /
                                (17.0 - 4.0 * param.k_ex) *
                                (1.0 - epsilon); // correction factor
    double C_ST = 16. / 9 * (3.0 - param.k_ex) * (5.0 - param.k_ex) /
                  (10.0 - 3.0 * param.k_ex);
    double E_ej, M_ej, Gamma, beta;

    E_ej = param.E_ej * (2 * PI * (1 - cos(param.open_angle)) / (4.0 * PI));
    Gamma = param.Gamma0;
    M_ej = E_ej / (Gamma * c_cnst * c_cnst);
    beta = sqrt(1 - 1. / Gamma / Gamma);

    DynamicThinShell dyn(M_ej, C_BM, C_ST, param.n_ex, param.k_ex, epsilon);

    y[0] = 1e10; // assume stellar radius
    y[1] = 0.0;
    y[2] = Gamma;
    y[3] = param.open_angle;
    y[4] = (Gamma - 1) * y[1] * c_cnst * c_cnst;
    y[5] = y[0] / y[2] / c_cnst;
    y[6] = y[0] / y[2] / c_cnst;

    double x = 0;
    dyn(x, y, dydx);
    rk5<DynamicThinShell> evol(dyn, y, dydx, yout, x, n);

    shock_dynamic = evol.integrate(time_array);
}

void GRBAfterglow::EvolutionThinShell(Jet& jet, double T) {
    const int n = 7;
    std::vector<double> y(n);
    std::vector<double> dydx(n);
    std::vector<double> yout(n);
    std::vector<double> shock_dynamic(n + 1);

    double epsilon = 0.0;
    double C_BM = epsilon + (9.0 - 2.0 * param.k_ex) /
                                (17.0 - 4.0 * param.k_ex) *
                                (1.0 - epsilon); // correction factor
    double C_ST = 16. / 9 * (3.0 - param.k_ex) * (5.0 - param.k_ex) /
                  (10.0 - 3.0 * param.k_ex);
    double E_ej, M_ej, Gamma, beta;
    double x = 0;

    // Initial values
    for (size_t i = 0; i < jet.theta.size(); i++) {
        // structured jet
        if ((param.gaussian_cone > 0) && (param.jet_index > 0)) {
            E_ej = param.E_ej *
                   pow(1 + jet.theta[i] / param.gaussian_cone * jet.theta[i] /
                               param.gaussian_cone,
                       -param.jet_index / 2) *
                   (2 * PI * (1 - cos(param.open_angle)) / (4 * PI));
            Gamma =
                param.Gamma0 * pow(1 + jet.theta[i] / param.gaussian_cone *
                                           jet.theta[i] / param.gaussian_cone,
                                   -param.jet_index / 2);
        } else if ((param.gaussian_cone > 0) && (param.jet_index <= 0)) {
            E_ej = param.E_ej *
                   exp(-0.5 * jet.theta[i] * jet.theta[i] /
                       param.gaussian_cone / param.gaussian_cone) *
                   (2 * PI * (1 - cos(param.open_angle)) / (4 * PI));
            Gamma = param.Gamma0;
        } else {
            E_ej = param.E_ej *
                   (2 * PI * (1 - cos(param.open_angle)) / (4.0 * PI));
            Gamma = param.Gamma0;
        }
        M_ej = E_ej / (Gamma * c_cnst * c_cnst);
        beta = sqrt(1 - 1. / Gamma / Gamma);

        jet.Gamma0[i] = Gamma;
        jet.M_ej[i] = M_ej;
        jet.E_ej[i] = E_ej;

        DynamicThinShell dyn(M_ej, C_BM, C_ST, param.n_ex, param.k_ex, epsilon);

        if (param.view_angle == 0) {

            x = 0;
            y[0] = 1e10; // assume stellar radius
            y[1] = 0.;
            y[2] = Gamma;
            y[3] = param.open_angle;
            y[4] = (Gamma - 1) * y[1] * c_cnst * c_cnst;
            y[5] = y[0] / y[2] / c_cnst;
            y[6] = y[0] / y[2] / c_cnst;

            dyn(x, y, dydx);
            rk5<DynamicThinShell> evol(dyn, y, dydx, yout, x, n);
            shock_dynamic = evol.integrateEATS(T / (1 + param.z), jet.mu[i][0]);

            for (size_t j = 0; j < jet.phi.size(); j++) {
                jet.radius[i][j] = shock_dynamic[0];
                jet.Gamma[i][j] = shock_dynamic[2];
                jet.dynamical_timescale[i][j] = shock_dynamic[5];
                jet.adiabatic_timescale[i][j] = shock_dynamic[6];
            }
        } else {
            for (size_t j = 0; j < jet.phi.size(); j++) {
                x = 0;
                y[0] = 1e10; // assume stellar radius
                y[1] = 0.;
                y[2] = Gamma;
                y[3] = param.open_angle;
                y[4] = (Gamma - 1) * y[1] * c_cnst * c_cnst;
                y[5] = y[0] / y[2] / c_cnst;
                y[6] = y[0] / y[2] / c_cnst;

                dyn(x, y, dydx);
                rk5<DynamicThinShell> evol(dyn, y, dydx, yout, x, n);
                shock_dynamic =
                    evol.integrateEATS(T / (1 + param.z), jet.mu[i][j]);

                jet.radius[i][j] = shock_dynamic[0];
                jet.Gamma[i][j] = shock_dynamic[2];
                jet.dynamical_timescale[i][j] = shock_dynamic[5];
                jet.adiabatic_timescale[i][j] = shock_dynamic[6];
            }
        }
    }
}

// dynamics with FS
DynamicThinShell::DynamicThinShell(const double M_ej, const double C_BM,
                                   const double C_ST, const double n_ex,
                                   const double k_ex, const double epsilon)
    : M_ej(M_ej), C_BM(C_BM), C_ST(C_ST), n_ex(n_ex), k_ex(k_ex),
      epsilon(epsilon) {}

void DynamicThinShell::operator()(const double x, std::vector<double>& y,
                                  std::vector<double>& dydx) {

    beta = sqrt(1 - 1. / y[2] / y[2]);
    gamma_hat = (4. + 1. / y[2]) / 3.;
    Gamma_eff = (gamma_hat * y[2] * y[2] - gamma_hat + 1) / y[2];
    dGamma_eff = (gamma_hat * y[2] * y[2] + gamma_hat - 1) / (y[2] * y[2]);
    A = 2. * PI * (1 - cos(y[3])) * y[0] * y[0] * mp *
        (beta * (C_BM - C_ST) + C_ST) * n_ex * pow(y[0], -k_ex);

    // double cs = sqrt(gamma_hat * (gamma_hat - 1) * (y[2] - 1) / (1 +
    // gamma_hat * (y[2] - 1)) *
    //                  c_cnst * c_cnst);

    // dr/dt
    dydx[0] = beta * c_cnst;

    // dm/dt = dm/dr * dr/dt
    dydx[1] = A * dydx[0];

    // dEad / dr
    dEad = -(gamma_hat - 1) * (3. / y[0] - 1 / y[2] * dydx[2] / dydx[0]) * y[4];

    // dEint/dt = dEint / dr * dr / dt
    dydx[4] =
        ((1.0 - epsilon) * (y[2] - 1.0) * A * c_cnst * c_cnst + dEad) * dydx[0];

    // dGamma/dt = dGamma /dr * dr/dt
    dydx[2] =
        -(A * y[2] * (y[2] * y[2] - 1) * (gamma_hat * y[2] - gamma_hat + 1) *
              c_cnst * c_cnst -
          (gamma_hat - 1) * y[2] * (gamma_hat * y[2] * y[2] - gamma_hat + 1) *
              3 * y[4] / y[0]) /
        (y[2] * y[2] * (M_ej + y[1]) * c_cnst * c_cnst +
         (gamma_hat * gamma_hat * y[2] * y[2] - gamma_hat * gamma_hat +
          3 * gamma_hat - 2) *
             y[4]) *
        dydx[0];

    // dtheta/ dt
    // dydx[3] = cs / c_cnst / y[0] / y[2] * dydx[0];
    dydx[3] = 0.0;

    dydx[5] =
        1.0 / c_cnst / beta / y[2] * dydx[0]; // comving dynamical timescale

    dydx[6] = 0.0;
    y[6] = 1.0 / (c_cnst * beta * y[2] *
                  (1.0 / y[0] - 1.0 / 3.0 / y[2] * dydx[2] /
                                    dydx[0])); // comving adiabatic timescale
}

double GRBAfterglow::DecelerationRadius() {
    double M_ej = param.E_ej / (param.Gamma0 * c_cnst * c_cnst);
    return pow((3 - param.k_ex) * M_ej /
                   (4 * PI * param.n_ex * mp * param.Gamma0),
               1 / (3 - param.k_ex));
}

double GRBAfterglow::DecelerationRadius(double E_ej, double Gamma) {
    double M_ej = E_ej / (Gamma * c_cnst * c_cnst);
    return pow((3 - param.k_ex) * M_ej / (4 * PI * param.n_ex * mp * Gamma),
               1 / (3 - param.k_ex));
}

void GRBAfterglow::EvolutionThinShellAnalytic(double T, double& Gamma,
                                              double& radius) {
    Gamma = 43.6 * exp(0.125 * log(param.E_ej / 1e53 / param.n_ex)) *
            exp(-0.375 * log(T / 1e3)) *
            exp(0.375 * log(1 + param.z)); // Sari98

    if (Gamma < 1.00001) {
        Gamma = 1.00001;
    }
    if (Gamma > param.Gamma0) {
        Gamma = param.Gamma0;
    }
    radius = 4. * Gamma * Gamma * c_cnst * T / (1 + param.z); // Sari98
}

double GRBAfterglow::MomentumInjRS(double Gamma) {
    double energy_kin; // kinetic energy
    if (param.spectral_e_RS > 2) {
        energy_kin = proton_mass * param.epsilon_e_RS / param.fraction_e_RS *
                     (param.spectral_e_RS - 2) / (param.spectral_e_RS - 1) *
                     (Gamma - 1);
    } else {
        std::cout << "Warnning: spectral index should larger than 2 !!! \n"
                  << std::endl;
        return 0;
    }
    return s.getElectron().Energy2Momentum(energy_kin + e_mass);
}

double GRBAfterglow::MomentumMaxRS(double mag_strength, double radius,
                                   double beta) {
    // t_syn  == t_acc
    double gamma_max_syn =
        sqrt(6 * PI * e_charge / param.eta_acc_e_RS / sigmaT / mag_strength);

    // t_ad  == t_acc
    double gamma_max_ad = e_charge * mag_strength * radius /
                          (param.eta_acc_e_RS * me * c_cnst * c_cnst * beta);

    return s.getElectron().Energy2Momentum(
        e_mass * std::min(gamma_max_syn, gamma_max_ad));
}

void GRBAfterglow::LorentzBoost(const std::vector<double>& photon_energy,
                                std::vector<double>& photon_spectrum,
                                double Gamma) {
    int num_p = photon_energy.size();
    std::vector<double> photon_temp(num_p, 0.0);
    for (size_t i = 0; i < num_p; i++) {
        photon_temp[i] = photon_energy[i] * photon_spectrum[i];
    }
    for (size_t i = 0; i < num_p; i++) {
        photon_spectrum[i] =
            (1 + param.z) * Gamma *
            utility.Interpolate(photon_energy, photon_temp,
                                (1 + param.z) * photon_energy[i] / Gamma);
    }
}

void GRBAfterglow::FluxNorm(const std::vector<double>& photon_energy,
                            std::vector<double>& photon_spectrum, double radius,
                            double Delta) {
    double flux_norm = eV2erg * (4.0 * PI * radius * radius * Delta) /
                       (4 * PI * param.dl * param.dl);
    int num_p = photon_energy.size();
    for (size_t i = 0; i < num_p; i++) {
        photon_spectrum[i] *= flux_norm * photon_energy[i];
    }
}

void GRBAfterglow::EvolutionThickShell(
    const std::vector<double>& time_array,
    std::vector<std::vector<double>>& jet_dyn) {
    const int n = 11;
    std::vector<double> y(n);
    std::vector<double> dydx(n);
    std::vector<double> yout(n);
    std::vector<double> time_array_corr = time_array;
    jet_dyn.resize(time_array.size(), std::vector<double>(n, 0.));

    double epsilon = 0.0;
    double C_BM = epsilon + (9.0 - 2.0 * param.k_ex) /
                                (17.0 - 4.0 * param.k_ex) *
                                (1.0 - epsilon); // correction factor
    double C_ST = 16. / 9 * (3.0 - param.k_ex) * (5.0 - param.k_ex) /
                  (10.0 - 3.0 * param.k_ex);
    // C_BM76 = 1.0;
    double Delta = param.T_ej / (1 + param.z) * c_cnst;
    double E_ej, M_ej, Gamma, beta;

    double R_dec = DecelerationRadius();
    double t_dec = R_dec / c_cnst / param.Gamma0 / param.Gamma0 * (1 + param.z);
    double RS_index = (7.0 - 2.0 * param.k_ex) / 2.0;
    if (Delta < R_dec / param.Gamma0 / param.Gamma0) {
        Delta = R_dec / param.Gamma0 / param.Gamma0;
        RS_index = 2.0;
    }

    E_ej = param.E_ej * (2 * PI * (1 - cos(param.open_angle)) / (4.0 * PI));
    Gamma = param.Gamma0;
    M_ej = E_ej / (Gamma * c_cnst * c_cnst);
    beta = sqrt(1 - 1. / Gamma / Gamma);
    DynamicThickShell dyn(M_ej, C_BM, C_ST, param.n_ex, param.k_ex, epsilon,
                          Gamma, beta, Delta, RS_index, param.open_angle);

    y[0] = 1e10; // assume stellar radius
    y[1] = 0.0;
    y[2] = Gamma;
    y[3] = param.open_angle;
    y[4] = (Gamma - 1) * y[1] * c_cnst * c_cnst;
    y[5] = M_ej / 1e6;
    y[6] = y[5] * c_cnst * c_cnst;
    y[7] = y[2];
    y[8] = y[0] / y[2] / c_cnst;
    y[9] = y[0] / y[2] / c_cnst;
    y[10] = 1e10;

    double x = 0;
    dyn(x, y, dydx);
    rk5<DynamicThickShell> evol(dyn, y, dydx, yout, x, n);

    for (int i = 0; i < time_array_corr.size(); i++) {
        time_array_corr[i] /= (1 + param.z);
    }
    jet_dyn = evol.integrate(time_array_corr);
}

void GRBAfterglow::EvolutionThickShell(Jet& jet, double T) {
    const int n = 11;
    std::vector<double> y(n);
    std::vector<double> dydx(n);
    std::vector<double> yout(n);
    std::vector<double> shock_dynamic(n);

    double epsilon = 0.0;
    double C_BM = epsilon + (9.0 - 2.0 * param.k_ex) /
                                (17.0 - 4.0 * param.k_ex) *
                                (1.0 - epsilon); // correction factor
    double C_ST = 16. / 9 * (3.0 - param.k_ex) * (5.0 - param.k_ex) /
                  (10.0 - 3.0 * param.k_ex);
    double E_ej, M_ej, Gamma, beta;

    // Initial values
    double RS_index, R_dec, Delta;
    for (size_t i = 0; i < jet.theta.size(); i++) {
        Delta = param.T_ej / (1 + param.z) * c_cnst;
        RS_index = (7.0 - 2.0 * param.k_ex) / 2.0;
        R_dec = DecelerationRadius();
        if (Delta < R_dec / 2. / param.Gamma0 / param.Gamma0) {
            Delta = R_dec / 2. / param.Gamma0 / param.Gamma0;
            RS_index = 2.0;
        }

        // structured jet
        if ((param.gaussian_cone > 0) && (param.jet_index > 0)) {
            E_ej = param.E_ej *
                   pow(1 + jet.theta[i] / param.gaussian_cone * jet.theta[i] /
                               param.gaussian_cone,
                       -param.jet_index / 2) *
                   (2 * PI * (1 - cos(param.open_angle)) / (4 * PI));
            Gamma =
                param.Gamma0 * pow(1 + jet.theta[i] / param.gaussian_cone *
                                           jet.theta[i] / param.gaussian_cone,
                                   -param.jet_index / 2);
        } else if ((param.gaussian_cone > 0) && (param.jet_index <= 0)) {
            E_ej = param.E_ej *
                   exp(-0.5 * jet.theta[i] * jet.theta[i] /
                       param.gaussian_cone / param.gaussian_cone) *
                   (2 * PI * (1 - cos(param.open_angle)) / (4 * PI));
            Gamma = param.Gamma0;
        } else {
            E_ej = param.E_ej *
                   (2 * PI * (1 - cos(param.open_angle)) / (4.0 * PI));
            Gamma = param.Gamma0;
        }
        M_ej = E_ej / (Gamma * c_cnst * c_cnst);
        beta = sqrt(1 - 1. / Gamma / Gamma);

        jet.Gamma0[i] = Gamma;
        jet.M_ej[i] = M_ej;
        jet.E_ej[i] = E_ej;

        DynamicThickShell dyn(M_ej, C_BM, C_ST, param.n_ex, param.k_ex, epsilon,
                              Gamma, beta, Delta, RS_index, param.open_angle);

        y[0] = 1e10; // assume stellar radius
        y[1] = 0.0;
        y[2] = Gamma;
        y[3] = param.open_angle;
        y[4] = (Gamma - 1) * y[1] * c_cnst * c_cnst;
        y[5] = M_ej / 1e6;
        y[6] = y[5] * c_cnst * c_cnst;
        y[7] = y[2];
        y[8] = y[0] / y[2] / c_cnst;
        y[9] = y[0] / y[2] / c_cnst;
        y[10] = y[0];

        double x = 0;
        dyn(x, y, dydx);
        rk5<DynamicThickShell> evol(dyn, y, dydx, yout, x, n);

        shock_dynamic = evol.integrateEATS(T / (1 + param.z), jet.mu[i][0]);

        for (size_t j = 0; j < jet.phi.size(); j++) {

            jet.radius[i][j] = shock_dynamic[0];
            jet.Gamma[i][j] = shock_dynamic[2];
            jet.mass_ejecta[i][j] = shock_dynamic[5];
            jet.dynamical_timescale[i][j] = shock_dynamic[8];
            jet.adiabatic_timescale[i][j] = shock_dynamic[9];
            jet.Gamma_RS[i][j] = shock_dynamic[7];
            jet.radius_RS[i][j] = shock_dynamic[10];
            jet.radius_cross[i][j] = dyn.radius_cross;
            jet.Gamma_cross[i][j] = dyn.Gamma_cross;
            jet.mass_cross[i][j] = dyn.mass_cross;
            jet.dynamical_timescale_cross[i][j] = dyn.dynamical_timescale_cross;
            jet.adiabatic_timescale_cross[i][j] = dyn.adiabatic_timescale_cross;
        }
    }
}

DynamicThickShell::DynamicThickShell(const double M_ej, const double C_BM,
                                     const double C_ST, const double n_ex,
                                     const double k_ex, const double epsilon,
                                     const double Gamma0, const double beta0,
                                     const double Delta, double RS_index,
                                     double theta_c)
    : M_ej(M_ej), C_BM(C_BM), C_ST(C_ST), n_ex(n_ex), k_ex(k_ex),
      epsilon(epsilon), Gamma0(Gamma0), beta0(beta0), Delta(Delta),
      RS_index(RS_index), theta_c(theta_c) {}

void DynamicThickShell::operator()(const double x, std::vector<double>& y,
                                   std::vector<double>& dydx) {
    beta = sqrt(1 - 1. / y[2] / y[2]);
    gamma_hat = (4. + 1. / y[2]) / 3.;
    Gamma_eff = (gamma_hat * y[2] * y[2] - gamma_hat + 1) / y[2];
    dGamma_eff = (gamma_hat * y[2] * y[2] + gamma_hat - 1) / (y[2] * y[2]);
    Gamma34 = Gamma0 * y[7] * (1 - beta0 * beta);

    A = 2. * PI * (1 - cos(y[3])) * y[0] * y[0] * mp *
        (beta * (C_BM - C_ST) + C_ST) * n_ex * pow(y[0], -k_ex);

    double cs = sqrt(gamma_hat * (gamma_hat - 1) * (y[2] - 1) /
                     (1 + gamma_hat * (y[2] - 1)) * c_cnst * c_cnst);

    // if (y[0] < 1e11) {
    // std::cout << "Gamma " << y[0] << " " << y[2] << std::endl;
    // }

    // dr/dt
    dydx[0] = beta * c_cnst;

    // dm/dt = dm/dr * dr/dt
    dydx[1] = A * dydx[0];

    // dEad / dr
    dEad2 =
        -(gamma_hat - 1) * (3. / y[0] - 1 / y[2] * dydx[2] / dydx[0]) * y[4];
    dEad3 = 0;

    // dGamma/dt = dGamma /dr * dr/dt
    dydx[2] =
        -((Gamma_eff + 1) * (y[2] - 1) * c_cnst * c_cnst * dydx[1] / dydx[0] +
          (y[2] - Gamma0 - Gamma_eff + Gamma_eff * Gamma34) * c_cnst * c_cnst *
              dydx[5] / dydx[0] +
          Gamma_eff * (dEad2 + dEad3)) /
        ((y[5] + y[1]) * c_cnst * c_cnst + (y[4] + y[6]) * dGamma_eff) *
        dydx[0];

    // dydx[3] = cs / c_cnst / y[0] / y[2] * dydx[0];
    // dydx[3] = cs * (y[2] + sqrt(y[2] * y[2] - 1)) / y[0];
    // if (beta * y[2] > 1 / 3. / sqrt(2) / theta_c) {
    //  dydx[3] = 0;
    //} else {
    //  dydx[3] = sqrt(3. / 8) * cs / c_cnst / y[2] * c_cnst / y[0];
    //}
    dydx[3] = 0;

    // dEint2/dt = dEint2 / dr * dr / dt
    dydx[4] = ((1.0 - epsilon) * (y[2] - 1.0) * A * c_cnst * c_cnst + dEad2) *
              dydx[0];

    // dM_ej3/dt = dM_ej3/dr * dr/dt
    dydx[5] = M_ej / Delta * (beta0 - beta) / beta * dydx[0];
    // dydx[5] = M_ej / Delta * (beta0 - beta) /
    //           (1 - Gamma0 / y[7] / ((gamma_hat * Gamma34 + 1.0) / (gamma_hat
    //           - 1.0))) * dydx[0];

    // dEint3/dt = dEint3 / dr * dr / dt
    // Eint3 = (Gamma34 - 1) * y[4] * c_cnst * c_cnst;
    dydx[6] = (Gamma34 - 1) * dydx[5] * c_cnst * c_cnst + dEad3 * dydx[0];

    // Gamma RS
    dydx[7] = 0;
    y[7] = y[2];

    if (have_cross) {
        if (y[5] > M_ej) {
            Gamma_cross = y[2];
            beta_cross = sqrt(1. - 1. / y[2] / y[2]);
            radius_cross = y[0];
            mass_cross = y[5];
            dynamical_timescale_cross = y[8];
            adiabatic_timescale_cross = y[9];
            have_cross = false;
        }
    }

    if (y[5] > M_ej) {
        dydx[5] = 0.0;
        dydx[6] = 0.0;
        GammaBeta =
            Gamma_cross * beta_cross * pow(y[10] / radius_cross, -RS_index);
        y[7] = sqrt(1 + GammaBeta * GammaBeta);
    }

    dydx[8] =
        1. / c_cnst / beta / y[2] * dydx[0]; // comving dynamical timescale
    dydx[9] = 0.;
    y[9] = 1. / (c_cnst * beta * y[2] *
                 (1. / y[0] - 1. / 3 / y[2] * dydx[2] / dydx[0]));

    // dr/dt
    dydx[10] = sqrt(1.0 - 1.0 / y[7] / y[7]) * c_cnst;
}

