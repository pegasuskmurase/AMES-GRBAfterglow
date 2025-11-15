//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file photonbackground.hpp
// \brief photonbackground class.

#include "module/photonbackground.hpp"
#include "module/constants.hpp"
#include "module/gammagamma.hpp"
#include "module/source.hpp"

#include <cmath>
#include <fstream>

Photonbackground::Photonbackground(Source& s) : s(s) {}

void Photonbackground::GreyBody(double T_ph, double u_ph) {
    int target_size = s.target.getSize();
    double ueps;
    std::vector<double> energy = s.target.getEnergy();
    for (size_t i = 0; i < target_size; i++) {
        ueps = 15.0 * u_ph / (PI * PI * PI * PI) *
               pow(energy[i] / (k_B * T_ph), 4) /
               (exp(energy[i] / (k_B * T_ph)) - 1.0);
        s.target.setSpectrum(i, ueps / energy[i] / energy[i]);
    }
}

void Photonbackground::BlackBody(double T_ph) {
    int target_size = s.target.getSize();
    double CMBtemp = T_ph * 1.3806503e-16 / e_cnst;
    std::vector<double> energy = s.target.getEnergy();
    for (size_t i = 0; i < target_size; i++) {
        if ((1 - (exp(-energy[i] / T_ph))) < 1e-6) {
            s.target.setSpectrum(i, 1.318e13 * (energy[i] * CMBtemp));
        } else {
            s.target.setSpectrum(i, 1.318e13 * energy[i] * energy[i] /
                                        (exp(energy[i] / CMBtemp) - 1));
        }
    }
}

void Photonbackground::Powerlaw(double eps_min, double eps_max, double index,
                                double u_ph) {
    double ueps;
    int target_size = s.target.getSize();
    std::vector<double> energy = s.target.getEnergy();
    for (size_t i = 0; i < target_size; i++) {
        if (energy[i] < eps_min) {
            s.target.setSpectrum(i, 0.0);
        } else {
            s.target.setSpectrum(i, energy[i] * pow(energy[i], -index) *
                                        exp(-energy[i] / eps_max));
        }
    }
    u_ph /= utility.Integrate(energy, s.target.getSpectrum());
    std::vector<double> spectrum = s.target.getSpectrum();
    for (size_t i = 0; i < target_size; i++) {
        if (energy[i] < eps_min) {
            s.target.setSpectrum(i, 0.0);
        } else {
            s.target.setSpectrum(i, u_ph * pow(energy[i], -index) *
                                        exp(-energy[i] / eps_max));
        }
    }
}

void Photonbackground::BrokenPowerlaw(double eps_min, double eps_break,
                                      double eps_max, double index_1,
                                      double index_2, double u_ph) {
    double ueps;
    int target_size = s.target.getSize();
    std::vector<double> energy = s.target.getEnergy();
    for (size_t i = 0; i < target_size; i++) {
        if (energy[i] < eps_min) {
            s.target.setSpectrum(i, 0.0);
        } else if ((energy[i] >= eps_min) && (energy[i] < eps_break)) {
            s.target.setSpectrum(i, energy[i] *
                                        pow(energy[i] / eps_break, -index_1));
        } else if ((energy[i] >= eps_break) && (energy[i] < eps_max)) {
            s.target.setSpectrum(i, energy[i] *
                                        pow(energy[i] / eps_break, -index_2));
        } else if (energy[i] >= eps_max) {
            s.target.setSpectrum(i, 0.0);
        }
    }
    u_ph /= utility.Integrate(energy, s.target.getSpectrum());
    std::vector<double> spectrum = s.target.getSpectrum();
    for (size_t i = 0; i < target_size; i++) {
        if (energy[i] < eps_min) {
            s.target.setSpectrum(i, 0.0);
        } else if ((energy[i] >= eps_min) && (energy[i] < eps_break)) {
            s.target.setSpectrum(i,
                                 u_ph * pow(energy[i] / eps_break, -index_1));
        } else if ((energy[i] >= eps_break) && (energy[i] < eps_max)) {
            s.target.setSpectrum(i,
                                 u_ph * pow(energy[i] / eps_break, -index_2));
        } else if (energy[i] >= eps_max) {
            s.target.setSpectrum(i, 0.0);
        }
    }
}

void Photonbackground::BandFunction(double eps_min, double eps_peak,
                                    double eps_max, double index_1,
                                    double index_2, double u_ph) {
    double ueps;
    double alpha_B = index_1;
    double beta_B = index_2;
    double x;
    int target_size = s.target.getSize();
    std::vector<double> energy = s.target.getEnergy();
    for (size_t i = 0; i < target_size; i++) {
        x = energy[i] / eps_peak;
        if (energy[i] < eps_min) {
            s.target.setSpectrum(i, 0.0);
        } else if (energy[i] >= eps_max) {
            s.target.setSpectrum(i, 0.0);
        } else {
            s.target.setSpectrum(i,
                                 energy[i] * 1. / eps_peak * pow(x, -alpha_B) *
                                     pow(1 + x * x, (-beta_B + alpha_B) / 2.));
        }
    }
    u_ph /= utility.Integrate(energy, s.target.getSpectrum());
    for (size_t i = 0; i < target_size; i++) {
        x = energy[i] / eps_peak;
        if (energy[i] < eps_min) {
            s.target.setSpectrum(i, 0.0);
        } else if (energy[i] >= eps_max) {
            s.target.setSpectrum(i, 0.0);
        } else {
            s.target.setSpectrum(i,
                                 u_ph * 1. / eps_peak * pow(x, -alpha_B) *
                                     pow(1 + x * x, (-beta_B + alpha_B) / 2.));
        }
    }
}

void Photonbackground::BandFunctionPL(double eps_min, double eps_peak,
                                      double eps_max, double index_1,
                                      double index_2, double u_ph) {
    double ueps;
    double alpha_B = index_1;
    double beta_B = index_2;
    double x;
    int target_size = s.target.getSize();
    std::vector<double> energy = s.target.getEnergy();
    for (size_t i = 0; i < target_size; i++) {
        x = energy[i] / eps_peak;
        s.target.setSpectrum(i, energy[i] * 1. / eps_peak * pow(x, -alpha_B) *
                                    pow(1 + x * x, (-beta_B + alpha_B) / 2.) *
                                    exp(-energy[i] / eps_max));
    }
    u_ph /= utility.Integrate(energy, s.target.getSpectrum());
    for (size_t i = 0; i < target_size; i++) {
        x = energy[i] / eps_peak;
        s.target.setSpectrum(i, u_ph * 1. / eps_peak * pow(x, alpha_B) *
                                    pow(1 + x * x, (-beta_B + alpha_B) / 2.) *
                                    exp(-energy[i] / eps_max));
    }
}

void Photonbackground::BrokenPowerlawSegment(double eps_min, double eps_1,
                                             double eps_2, double eps_max,
                                             double index_1, double index_2,
                                             double index_3, double flux_peak) {
    int target_size = s.target.getSize();
    std::vector<double> energy = s.target.getEnergy();
    for (size_t i = 0; i < target_size; i++) {
        if (energy[i] < eps_min) {
            s.target.setSpectrum(i, 0.0);
        } else if ((energy[i] >= eps_min) && (energy[i] < eps_1)) {
            s.target.setSpectrum(i,
                                 flux_peak * pow(energy[i] / eps_1, -index_1));
        } else if ((energy[i] >= eps_1) && (energy[i] < eps_2)) {
            s.target.setSpectrum(i,
                                 flux_peak * pow(energy[i] / eps_1, -index_2));
        } else if ((energy[i] >= eps_2) && (energy[i] < eps_max)) {
            s.target.setSpectrum(i, flux_peak *
                                        pow(energy[i] / eps_2, -index_3) *
                                        pow(eps_2 / eps_1, index_2));
        } else if (energy[i] >= eps_max) {
            s.target.setSpectrum(i, 0.0);
        }
    }
}

void Photonbackground::CMB(double redshift) {
    int target_size = s.target.getSize();
    double CMBtemp = 2.725 * 1.3806503e-16 / e_cnst;
    std::vector<double> energy = s.target.getEnergy();
    for (size_t i = 0; i < target_size; i++) {
        if ((1 - (exp(-energy[i] / T_CMB / (1 + redshift)))) < 1e-6) {
            s.target.setSpectrum(i, 1.318e13 *
                                        (energy[i] * CMBtemp * (1 + redshift)));
        } else {
            s.target.setSpectrum(
                i, 1.318e13 * energy[i] * energy[i] /
                       (exp(energy[i] / CMBtemp / (1 + redshift)) - 1));
        }
    }
}

void Photonbackground::EBL(double redshift, std::string EBL_model) {
    int target_size = s.target.getSize();
    CMB(redshift);
    std::vector<double> cmb = s.target.getSpectrum();
    if (EBL_model == "EBL_Franceschini08") {
        EBL_Franceschini08(redshift);
    } else if (EBL_model == "EBL_Gilmore12") {
        EBL_Gilmore12(redshift);
    } else if (EBL_model == "EBL_Finke10") {
        EBL_Finke10(redshift);
    } else {
        EBL_Gilmore12(redshift);
    }
    std::vector<double> cib = s.target.getSpectrum();
    for (size_t i = 0; i < target_size; i++) {
        s.target.setSpectrum(i, cmb[i] + cib[i]);
    }
}

void Photonbackground::EBLTau(double redshift, std::string EBL_model) {
    if (EBL_model == "EBL_Franceschini08") {
        EBLTau_Franceschini08(redshift);
    } else if (EBL_model == "EBL_Gilmore12") {
        EBLTau_Gilmore12(redshift);
    } else {
        std::cout << "No such table!" << std::endl;
    }
}

void Photonbackground::EBLAttenuation(const std::vector<double>& photon_energy,
                                      std::vector<double>& photon_spectrum) {
    double tau;
    int num_p = photon_spectrum.size();
    for (size_t i = 0; i < num_p; i++) {
        tau = utility.Interpolate(s.photon.energy, s.photon.optdepth,
                                  photon_energy[i]);
        photon_spectrum[i] *= exp(-tau);
    }
}

double Photonbackground::EBLAttenuation(double photon_energy,
                                        double photon_spectrum) {
    double tau;
    tau =
        utility.Interpolate(s.photon.energy, s.photon.optdepth, photon_energy);
    photon_spectrum *= exp(-tau);
    return photon_spectrum;
}

void Photonbackground::EBL_Franceschini08(double redshift) {

    std::vector<double> raw_redshift;
    std::vector<double> raw_energy;
    std::vector<std::vector<double>> raw_density; // 2d
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Franceschini_2008/target_redshift.dat",
                       raw_redshift);
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Franceschini_2008/target_energies.dat",
                       raw_energy);
    utility.ReadFile2D(std::string(DATA_PATH) +
                           "/data/EBL/Franceschini_2008/target_density.dat",
                       raw_density); // proper density

    std::vector<double> target = s.target.getEnergy();
    std::vector<double> raw_energy_z(raw_energy.size(), 0.0);
    for (size_t i = 0; i < raw_energy_z.size(); i++) {
        raw_energy_z[i] = raw_energy[i] * (1 + redshift);

        const int target_size = target.size();
        for (size_t i = 0; i < target_size; i++) {
            s.target.setSpectrum(
                i, utility.Interpolate2D(raw_redshift, raw_energy_z,
                                         raw_density, redshift, target[i]));
        }
    }
}

void Photonbackground::EBLTau_Franceschini08(double redshift) {
    std::vector<double> opdep_redshift;
    std::vector<double> opdep_energy;
    std::vector<std::vector<double>> raw_tau; // 2d
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Franceschini_2008/opdep_redshift.dat",
                       opdep_redshift);
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Franceschini_2008/opdep_energies.dat",
                       opdep_energy);
    utility.ReadFile2D(std::string(DATA_PATH) +
                           "/data/EBL/Franceschini_2008/opdep.dat",
                       raw_tau);

    std::vector<double> photon_energy = s.photon.getEnergy();
    const int num_p = photon_energy.size();
    for (size_t i = 0; i < num_p; i++) {
        s.photon.setOptdepth(
            i, utility.Interpolate2D(opdep_redshift, opdep_energy, raw_tau,
                                     redshift, photon_energy[i]));
    }
}

void Photonbackground::EBL_Gilmore12(double redshift) {

    std::vector<double> raw_redshift;
    std::vector<double> raw_energy;
    std::vector<std::vector<double>> raw_density; // 2d
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Gilmore_2012/target_redshift.dat",
                       raw_redshift);
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Gilmore_2012/target_energies.dat",
                       raw_energy);
    utility.ReadFile2D(std::string(DATA_PATH) +
                           "/data/EBL/Gilmore_2012/target_density.dat",
                       raw_density); // proper density

    std::vector<double> target = s.target.getEnergy();
    const int target_size = target.size();
    for (size_t i = 0; i < target_size; i++) {
        s.target.setSpectrum(i, utility.Interpolate2D(raw_redshift, raw_energy,
                                                      raw_density, redshift,
                                                      target[i]));
    }
}

void Photonbackground::EBLTau_Gilmore12(double redshift) {
    std::vector<double> opdep_redshift;
    std::vector<double> opdep_energy;
    std::vector<std::vector<double>> raw_tau; // 2d
    /*
    utility.ReadFile1D(std::string(DATA_PATH) +
    "/data/EBL/Gilmore_2012/opdep_redshift.dat", opdep_redshift);
    utility.ReadFile1D(std::string(DATA_PATH) +
    "/data/EBL/Gilmore_2012/opdep_energies.dat", opdep_energy);
    utility.ReadFile2D(std::string(DATA_PATH) +
    "/data/EBL/Gilmore_2012/opdep.dat", raw_tau);
    */
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Gilmore_2012/opdep_redshift_wCMB.dat",
                       opdep_redshift);
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Gilmore_2012/opdep_energies_wCMB.dat",
                       opdep_energy);
    utility.ReadFile2D(std::string(DATA_PATH) +
                           "/data/EBL/Gilmore_2012/opdep_wCMB.dat",
                       raw_tau);

    std::vector<double> photon_energy = s.photon.getEnergy();
    const int num_p = photon_energy.size();
    for (size_t i = 0; i < num_p; i++) {
        s.photon.setOptdepth(
            i, utility.Interpolate2D(opdep_redshift, opdep_energy, raw_tau,
                                     redshift, photon_energy[i]));
    }
}

void Photonbackground::EBL_Finke10(double redshift) {
    std::vector<double> raw_redshift;
    std::vector<double> raw_energy;
    std::vector<std::vector<double>> raw_density; // 2d
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Finke_2010/target_redshift.dat",
                       raw_redshift);
    utility.ReadFile1D(std::string(DATA_PATH) +
                           "/data/EBL/Finke_2010/target_energies.dat",
                       raw_energy);
    utility.ReadFile2D(std::string(DATA_PATH) +
                           "/data/EBL/Finke_2010/target_density.dat",
                       raw_density);

    std::vector<double> target = s.target.getEnergy();
    const int target_size = target.size();
    for (size_t i = 0; i < target_size; i++) {
        s.target.setSpectrum(i, utility.Interpolate2D(raw_redshift, raw_energy,
                                                      raw_density, redshift,
                                                      target[i]));
    }
}

void Photonbackground::calcOpticalDepth(double redshift,
                                        std::string EBL_model) {
    CosmologicalOpticalDepth(redshift, EBL_model);

    std::ofstream mystream;
    if (output_file.size() > 0) {
        mystream.open(output_file);
    } else {
        mystream.open("optical_depth.dat");
    }
    int N = s.photon.optdepth.size();
    for (size_t i = 0; i < N; i++) {
        mystream << s.photon.energy[i] << " " << s.photon.optdepth[i]
                 << std::endl;
    }
    mystream.close();
}

void Photonbackground::CosmologicalOpticalDepth(double redshift,
                                                std::string EBL_model) {
    const int N = s.photon.energy.size();
    s.photon.optdepth.assign(N, 0.0);

    Cosmology cosmology;
    double z, z0, z1, dz, dl, dt; // dl(light travel distance)
    z0 = redshift;
    z = z0;
    z1 = z0;

    EBL(z, EBL_model);
    // EBL_Gilmore12(z);

    double l = cosmology.Redshift2LighttravelDistance(z0);
    dl = 0.1 * Mpc; // TODO
    dt = dl / c_cnst;
    int istep = 0;
    int i = 0;
    bool have_redshift_evolution = true;

    GammaGamma gg(s);
    gg.Intetime(s.target.energy, s.target.spectrum, s.photon.energy,
                s.photon.intetime, z);

    while (z > 0) {
        for (size_t i = 0; i < N; i++) {
            s.photon.optdepth[i] += s.photon.intetime[i] * dl / c_cnst;
        }
        if (have_redshift_evolution) {
            if (z - z1 < 0.00001) {
                z1 -= 0.01;
                EBL(z, EBL_model);
                gg.Intetime(s.target.energy, s.target.spectrum, s.photon.energy,
                            s.photon.intetime, z);
                std::cout << "z " << z << " " << z1 << " " << dz << " "
                          << z + dz << std::endl;
            }
        }
        dz = dl / cosmology.HubbleDistance() * (1 + z) * cosmology.Ez(z);
        z -= dz;
        istep++;
    }

    /*
    GammaGamma gg(s);
    size_t angular_num = 51;
    std::vector<double> angular;
    std::vector<double> angular_bin;
    angular.resize(angular_num);
    angular_bin.resize(angular_num);
    for (size_t i = 0; i < angular_num; i++) {
      angular[i] = -1. + 2. * i / (angular_num - 1.);
      angular_bin[i] = 2. / (angular_num - 1);
    }

    double mandels;
    int N_target = s.target.energy.size();
    int N_photon = s.photon.energy.size();
    std::vector<double> angular_temp(angular_num, 0.0);
    double Emin, Eth;
    for (size_t i = 0; i < N_photon; i++) {
      for (size_t j = 0; j < angular_num - 1; j++) {
        Eth = 2 * e_mass * e_mass / s.photon.energy[i] / (1 - angular[j]);
        for (size_t k = 0; k < N_target; k++) {
          Emin = Eth / (1 + z);
          mandels = 2. * s.photon.energy[i] * (1 + z) * s.target.energy[k] * (1
    - angular[j]) / (e_mass * e_mass); if (s.target.energy[k] > Emin) {
            s.target.temp[k] = s.target.spectrum[k] * gg.GGCross(mandels) /
    c_cnst; } else { s.target.temp[k] = 0.0;
          }
        }
        angular_temp[j] = 0.5 * utility.Integrate(s.target.energy,
    s.target.temp);
      }
      s.photon.temp[i] = utility.Integrate(angular, angular_temp);
    }

    dl = 0.1 * Mpc; // TODO
    dt = dl / c_cnst;
    int istep = 0;
    int i = 0;
    bool have_redshift_evolution = true;
    while (z > 0) {
      for (size_t i = 0; i < N; i++) {
        s.photon.optdepth[i] += s.photon.temp[i] * dl;
      }

      if (have_redshift_evolution) {
        if (z - z1 < 0.00001) {
          z1 -= 0.01;
          //EBL(z, EBL_model);
          EBL_Gilmore12(z);
          for (size_t i = 0; i < N_photon; i++) {
            for (size_t j = 0; j < angular_num - 1; j++) {
              Eth = 2 * e_mass * e_mass / s.photon.energy[i] / (1 - angular[j]);
              for (size_t k = 0; k < N_target; k++) {
                Emin = Eth / (1 + z);
                mandels = 2. * s.photon.energy[i] * (1 + z) * s.target.energy[k]
    * (1 - angular[j]) / (e_mass * e_mass); if (s.target.energy[k] > Emin) {
                  s.target.temp[k] = s.target.spectrum[k] * gg.GGCross(mandels)
    / c_cnst; } else { s.target.temp[k] = 0.0;
                }
              }
              angular_temp[j] = 0.5 * utility.Integrate(s.target.energy,
    s.target.temp);
            }
            s.photon.temp[i] = utility.Integrate(angular, angular_temp);
          }

          std::cout << istep << " " << z << " " << dz << " " << dt << " " << dl
    << " " << std::endl;
        }
      }
      dz = dl / cosmology.HubbleDistance() * (1 + z) * cosmology.Ez(z);
      z -= dz;
      istep++;
    }
    std::cout << "final " << istep << " " << z << " " << dz << " " << dt << " "
    << dl << " "
              << std::endl;
    */
}
