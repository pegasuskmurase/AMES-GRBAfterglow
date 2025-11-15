//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file source.hpp
// \brief This is the source class.
#include "module/source.hpp"
#include "module/constants.hpp"
#include <fstream>

Source::Source() {
    setParticle();
    setParam();
}

Source::Source(Source& s) {
    photon = s.photon;
    photon_syn = s.photon_syn;
    photon_IC = s.photon_IC;
    photon_pm = s.photon_pm;
    photon_pa = s.photon_pa;
    electron_pa = s.electron_pa;
    electron_pm = s.electron_pm;
    target = s.target;
    electron = s.electron;
    neutrino = s.neutrino;
    pionzero = s.pionzero;
    pionminus = s.pionminus;
    pionplus = s.pionplus;
    muonminusL = s.muonminusL;
    muonplusL = s.muonplusL;
    muonminusR = s.muonminusR;
    muonplusR = s.muonplusR;
    proton = s.proton;
    neutron = s.neutron;
    nucleus = s.nucleus;
    param = s.param;

    if (s.syn)
        syn = std::make_unique<Synchrotron>(*s.syn);
    if (s.IC)
        IC = std::make_unique<InverseCompton>(*s.IC);
    if (s.gg)
        gg = std::make_unique<GammaGamma>(*s.gg);
    if (s.ph)
        ph = std::make_unique<Photonbackground>(*s.ph);
    if (s.ED)
        ED = std::make_unique<ElectronDistribution>(*s.ED);

    std::cout << "Calling the copy constructor " << std::endl;
}

void Source::InitSource() {
    syn = std::make_unique<Synchrotron>(*this);
    IC = std::make_unique<InverseCompton>(*this);
    gg = std::make_unique<GammaGamma>(*this);
    ph = std::make_unique<Photonbackground>(*this);
    ED = std::make_unique<ElectronDistribution>(*this);
}

void Source::clear() {
    syn.reset();
    IC.reset();
    gg.reset();
}

void Source::setParticle() {
    int ID;
    std::string name;

    ID = 22;
    name = "target";
    target.setID(ID);
    target.setName(name);
    target.setCharge(0);
    target.setMassNumber(0);
    target.setMass(0.0);
    target.setLife(1e30);
    target.setMomentum(361, 1e-8, 1e10);

    ID = 22;
    name = "photon";
    photon.setID(ID);
    photon.setName(name);
    photon.setCharge(0);
    photon.setMassNumber(0);
    photon.setMass(0.0);
    photon.setLife(1e30);
    photon.setMomentum(561, 1e-8, 1e20);

    ID = 2201;
    name = "photon_syn";
    photon_syn.setID(ID);
    photon_syn.setName(name);
    photon_syn.setCharge(0);
    photon_syn.setMassNumber(0);
    photon_syn.setMass(0.0);
    photon_syn.setLife(1e30);
    photon_syn.setMomentum(561, 1e-8, 1e20);

    ID = 2202;
    name = "photon_IC";
    photon_IC.setID(ID);
    photon_IC.setName(name);
    photon_IC.setCharge(0);
    photon_IC.setMassNumber(0);
    photon_IC.setMass(0.0);
    photon_IC.setLife(1e30);
    photon_IC.setMomentum(561, 1e-8, 1e20);

    ID = 2203;
    name = "photon_pm";
    photon_pm.setID(ID);
    photon_pm.setName(name);
    photon_pm.setCharge(0);
    photon_pm.setMassNumber(0);
    photon_pm.setMass(0.0);
    photon_pm.setLife(1e30);
    photon_pm.setMomentum(561, 1e-8, 1e20);

    ID = 2204;
    name = "photon_pa";
    photon_pa.setID(ID);
    photon_pa.setName(name);
    photon_pa.setCharge(0);
    photon_pa.setMassNumber(0);
    photon_pa.setMass(0.0);
    photon_pa.setLife(1e30);
    photon_pa.setMomentum(561, 1e-8, 1e20);

    ID = 11;
    name = "electron";
    electron.setID(ID);
    electron.setName(name);
    electron.setCharge(1);
    electron.setMassNumber(0);
    electron.setMass(e_mass);
    electron.setLife(1e30);
    electron.setMomentum(381, 1e1, 1e20);

    ID = 1100;
    name = "electron_prim";
    electron_prim.setID(ID);
    electron_prim.setName(name);
    electron_prim.setCharge(1);
    electron_prim.setMassNumber(0);
    electron_prim.setMass(e_mass);
    electron_prim.setLife(1e30);
    electron_prim.setMomentum(381, 1e1, 1e20);

    ID = 1101;
    name = "electron_pm";
    electron_pm.setID(ID);
    electron_pm.setName(name);
    electron_pm.setCharge(1);
    electron_pm.setMassNumber(0);
    electron_pm.setMass(e_mass);
    electron_pm.setLife(1e30);
    electron_pm.setMomentum(381, 1e1, 1e20);

    ID = 1102;
    name = "electron_pa";
    electron_pa.setID(ID);
    electron_pa.setName(name);
    electron_pa.setCharge(1);
    electron_pa.setMassNumber(0);
    electron_pa.setMass(e_mass);
    electron_pa.setLife(1e30);
    electron_pa.setMomentum(381, 1e1, 1e20);

    ID = 12; // represents all flavor neutrinos
    name = "neutrino";
    neutrino.setID(ID);
    neutrino.setName(name);
    neutrino.setCharge(0);
    neutrino.setMassNumber(0);
    neutrino.setMass(0.0);
    neutrino.setLife(1e30);
    neutrino.setMomentum(161, 1e6, 1e22);

    ID = 10006;
    name = "pionzero";
    pionzero.setID(ID);
    pionzero.setName(name);
    pionzero.setCharge(0);
    pionzero.setMassNumber(0);
    pionzero.setMass(piz_Mass);
    pionzero.setLife(piz_Life);
    pionzero.setMomentum(161, 1e6, 1e22);

    ID = 10007;
    name = "pionminus";
    pionminus.setID(ID);
    pionminus.setName(name);
    pionminus.setCharge(1);
    pionminus.setMassNumber(0);
    pionminus.setMass(pic_Mass);
    pionminus.setLife(pic_Life);
    pionminus.setMomentum(161, 1e6, 1e22);

    ID = 10008;
    name = "pionplus";
    pionplus.setID(ID);
    pionplus.setName(name);
    pionplus.setCharge(1);
    pionplus.setMassNumber(0);
    pionplus.setMass(pic_Mass);
    pionplus.setLife(pic_Life);
    pionplus.setMomentum(161, 1e6, 1e22);

    ID = 10009;
    name = "kaonminus";
    kaonminus.setID(ID);
    kaonminus.setName(name);
    kaonminus.setCharge(1);
    kaonminus.setMassNumber(0);
    kaonminus.setMass(kac_Mass);
    kaonminus.setLife(kac_Life);
    kaonminus.setMomentum(161, 1e6, 1e22);

    ID = 10010;
    name = "kaonplus";
    kaonplus.setID(ID);
    kaonplus.setName(name);
    kaonplus.setCharge(1);
    kaonplus.setMassNumber(0);
    kaonplus.setMass(kac_Mass);
    kaonplus.setLife(kac_Life);
    kaonplus.setMomentum(161, 1e6, 1e22);

    ID = 100041;
    name = "muonminusL";
    muonminusL.setID(ID);
    muonminusL.setName(name);
    muonminusL.setCharge(-1);
    muonminusL.setMassNumber(0);
    muonminusL.setMass(mu_Mass);
    muonminusL.setLife(mu_Life);
    muonminusL.setMomentum(161, 1e6, 1e22);

    ID = 100042;
    name = "muonminusR";
    muonminusR.setID(ID);
    muonminusR.setName(name);
    muonminusR.setCharge(-1);
    muonminusR.setMassNumber(0);
    muonminusR.setMass(mu_Mass);
    muonminusR.setLife(mu_Life);
    muonminusR.setMomentum(161, 1e6, 1e22);

    ID = 100051;
    name = "muonplusL";
    muonplusL.setID(ID);
    muonplusL.setName(name);
    muonplusL.setCharge(1);
    muonplusL.setMassNumber(0);
    muonplusL.setMass(mu_Mass);
    muonplusL.setLife(mu_Life);
    muonplusL.setMomentum(161, 1e6, 1e22);

    ID = 100052;
    name = "muonplusR";
    muonplusR.setID(ID);
    muonplusR.setName(name);
    muonplusR.setCharge(1);
    muonplusR.setMassNumber(0);
    muonplusR.setMass(mu_Mass);
    muonplusR.setLife(mu_Life);
    muonplusR.setMomentum(161, 1e6, 1e22);

    ID = 15;
    name = "muon";
    muon.setID(ID);
    muon.setName(name);
    muon.setCharge(1);
    muon.setMassNumber(0);
    muon.setMass(mu_Mass);
    muon.setLife(mu_Life);
    muon.setMomentum(161, 1e6, 1e22);

    ID = 13;
    name = "proton";
    proton.setID(ID);
    proton.setName(name);
    proton.setCharge(1);
    proton.setMassNumber(1);
    proton.setMass(proton_mass);
    proton.setLife(1e30);
    proton.setMomentum(161, 1e6, 1e22);

    ID = 14;
    name = "neutron";
    neutron.setID(ID);
    neutron.setName(name);
    neutron.setCharge(0);
    neutron.setMassNumber(1);
    neutron.setMass(neutron_mass);
    neutron.setLife(neutron_Life);
    neutron.setMomentum(161, 1e6, 1e22);

    ID = 202; // need to be set later
    name = "helium";
    helium.setID(ID);
    helium.setName(name);
    helium.setCharge(2);
    helium.setMassNumber(4);
    helium.setMass(proton_mass * 2.0 + neutron_mass * 2.0);
    helium.setLife(1e30);
    helium.setMomentum(161, 1e6, 1e22);

    ID = 606; // need to be set later
    name = "carbon";
    carbon.setID(ID);
    carbon.setName(name);
    carbon.setCharge(6);
    carbon.setMassNumber(12);
    carbon.setMass(proton_mass * 6.0 + neutron_mass * 6.0);
    carbon.setLife(1e30);
    carbon.setMomentum(161, 1e6, 1e22);

    ID = 808; // need to be set later
    name = "oxygen";
    oxygen.setID(ID);
    oxygen.setName(name);
    oxygen.setCharge(8);
    oxygen.setMassNumber(16);
    oxygen.setMass(proton_mass * 8.0 + neutron_mass * 8.0);
    oxygen.setLife(1e30);
    oxygen.setMomentum(161, 1e6, 1e22);

    ID = 1414; // need to be set later
    name = "silicon";
    silicon.setID(ID);
    silicon.setName(name);
    silicon.setCharge(14);
    silicon.setMassNumber(28);
    silicon.setMass(proton_mass * 14.0 + neutron_mass * 14.0);
    silicon.setLife(1e30);
    silicon.setMomentum(161, 1e6, 1e22);

    ID = 2630; // need to be set later
    name = "iron";
    iron.setID(ID);
    iron.setName(name);
    iron.setCharge(26);
    iron.setMassNumber(56);
    iron.setMass(proton_mass * 26.0 + neutron_mass * 30.0);
    iron.setLife(1e30);
    iron.setMomentum(161, 1e6, 1e22);

    ID = 2630; // need to be set later
    name = "iron";
    nucleus.setID(ID);
    nucleus.setName(name);
    nucleus.setCharge(26);
    nucleus.setMassNumber(56);
    nucleus.setMass(proton_mass * 26.0 + neutron_mass * 30.0);
    nucleus.setLife(1e30);
    nucleus.setMomentum(161, 1e6, 1e22);
}

void Source::setParam() {
    param.setName("CenA");
    param.setGeometry("blob");
    param.setRedshift(0.00785);
    param.setLumDistance(2.7 * Mpc);
    // param.setMagStrength(1e-4);         // G
    param.setEmissionRadius(1.0 * Kpc); // G
    param.setEmissionVolume(1.);
    param.setDopplerFactor(1.);
    param.setAmbientDensity(1.);
    param.setZeff(2.5);
}

void Source::setNucleus(int mass_number, int charge_number) {
    nucleus.setMass(charge_number * proton_mass +
                    (mass_number - charge_number) * neutron_mass);
    nucleus.setCharge(charge_number);
    nucleus.setMassNumber(mass_number);
}

void Source::showParam() {
    std::cout << " name: the name of the source "
              << " \n "
              << "geometry: can be 'blob' or 'shell'"
              << " \n "
              << "redshift: source redshift "
              << " \n "
              << "lum_distance: source luminosity distance [cm]"
              << " \n "
              << "mag_strength: the strength of magnetic field [G]"
              << " \n "
              << "emission_radiu: emission radius for spherical blob [cm]"
              << " \n "
              << "emission_volume: emission volume [cm^3]"
              << " \n "
              << "doppler_factor: source Doppler factor "
              << " \n "
              << "ambient_density: ambient matter density [cm^-3] "
              << " \n " << std::endl;
}

void Source::LorentzBoost(const std::vector<double>& photon_energy,
                          std::vector<double>& photon_spectrum,
                          double doppler_factor) {

    int num_p = photon_energy.size();
    std::vector<double> photon_temp(num_p, 0.0);
    double D3 = doppler_factor * doppler_factor * doppler_factor;
    for (size_t i = 0; i < num_p; i++) {
        photon_temp[i] = photon_spectrum[i] * photon_energy[i];
    }
    for (size_t i = 0; i < num_p; i++) {
        photon_spectrum[i] =
            (1 + param.getRedshift()) * D3 *
            utility.Interpolate(photon_energy, photon_temp,
                                (1 + param.getRedshift()) * photon_energy[i] /
                                    doppler_factor);
    }
}

void Source::LorentzBoost(const std::vector<double>& particle_energy,
                          std::vector<double>& particle_spectrum) {

    int num_p = particle_energy.size();
    std::vector<double> particle_temp(num_p, 0.0);
    double D3 = param.getDopplerFactor() * param.getDopplerFactor() *
                param.getDopplerFactor();
    for (size_t i = 0; i < num_p; i++) {
        particle_temp[i] = particle_spectrum[i] * particle_energy[i];
    }
    for (size_t i = 0; i < num_p; i++) {
        particle_spectrum[i] =
            (1 + param.getRedshift()) * D3 *
            utility.Interpolate(particle_energy, particle_temp,
                                (1 + param.getRedshift()) * particle_energy[i] /
                                    param.getDopplerFactor());
    }
}

void Source::FluxNorm(const std::vector<double>& particle_energy,
                      std::vector<double>& particle_spectrum) {
    double flux_norm =
        (4.0 * PI / 3.0 * param.getEmissionRadius() *
         param.getEmissionRadius() * param.getEmissionRadius()) /
        (4 * PI * param.getLumDistance() * param.getLumDistance());
    int num_p = particle_energy.size();
    for (size_t i = 0; i < num_p; i++) {
        particle_spectrum[i] *= flux_norm / particle_energy[i];
    }
}

void Source::LorentzBoostShell(const std::vector<double>& photon_energy,
                               std::vector<double>& photon_spectrum,
                               double doppler_factor) {

    int num_p = photon_energy.size();
    std::vector<double> photon_temp(num_p, 0.0);
    double D1 = doppler_factor;
    for (size_t i = 0; i < num_p; i++) {
        photon_temp[i] = photon_spectrum[i] * photon_energy[i];
    }
    for (size_t i = 0; i < num_p; i++) {
        photon_spectrum[i] =
            (1 + param.getRedshift()) * D1 *
            utility.Interpolate(photon_energy, photon_temp,
                                (1 + param.getRedshift()) * photon_energy[i] /
                                    doppler_factor);
    }
}

void Source::LorentzBoostShell(const std::vector<double>& photon_energy,
                               std::vector<double>& photon_spectrum) {

    int num_p = photon_energy.size();
    std::vector<double> photon_temp(num_p, 0.0);
    double D1 = param.getDopplerFactor();
    for (size_t i = 0; i < num_p; i++) {
        photon_temp[i] = photon_spectrum[i] * photon_energy[i];
    }
    for (size_t i = 0; i < num_p; i++) {
        photon_spectrum[i] =
            (1 + param.getRedshift()) * D1 *
            utility.Interpolate(photon_energy, photon_temp,
                                (1 + param.getRedshift()) * photon_energy[i] /
                                    param.getDopplerFactor());
    }
}

void Source::FluxNormShell(const std::vector<double>& photon_energy,
                           std::vector<double>& photon_spectrum) {
    double flux_norm =
        (4.0 * PI / 1.0 * param.getEmissionRadius() *
         param.getEmissionRadius() * param.getShellWidth()) /
        (4 * PI * param.getLumDistance() * param.getLumDistance());
    int num_p = photon_energy.size();
    for (size_t i = 0; i < num_p; i++) {
        photon_spectrum[i] *= flux_norm / photon_energy[i];
    }
}

void Source::Output(std::vector<ParticleTransport>& p) {
    std::ofstream mystream_energy;
    mystream_energy.open(output_folder + "/LeptohadronicEnergy.dat");
    for (auto& particle : p) {
        std::cout << particle.ID << " " << particle.Z << " " << particle.N
                  << std::endl;
        mystream_energy << particle.ID << " " << particle.Z << " "
                        << particle.N;
        for (size_t i = 0; i < particle.x.size(); i++) {
            mystream_energy << " " << particle.x_e[i];
        }
        mystream_energy << std::endl;
    }
    mystream_energy.close();

    std::ofstream mystream;
    mystream.open(output_folder + "/LeptohadronicSpectrum.dat");

    for (auto& particle : p) {
        mystream << particle.ID << " " << particle.Z << " " << particle.N;
        for (size_t i = 0; i < particle.x.size(); i++) {
            mystream << " " << particle.u_middle[i];
        }
        mystream << std::endl;
    }
    mystream.close();
}

void Source::Output(std::vector<ParticleTransportChangCooper>& p) {

    std::ofstream mystream_energy;
    mystream_energy.open(output_folder + "/LeptohadronicEnergy.dat");
    for (auto& particle : p) {
        mystream_energy << particle.ID << " " << particle.Z << " "
                        << particle.N;
        for (size_t i = 0; i < particle.x.size(); i++) {
            mystream_energy << " " << particle.x_e[i];
        }
        mystream_energy << std::endl;
    }
    mystream_energy.close();

    std::ofstream mystream;
    mystream.open(output_folder + "/LeptohadronicSpectrum.dat");
    double timescale_escape = param.getEmissionRadius() / c_cnst;
    for (auto& particle : p) {
        for (size_t i = 0; i < particle.x.size(); i++) {
            particle.u_middle[i] = particle.x_e[i] / particle.x[i] *
                                   particle.u[i] / timescale_escape;
        }

        /*
        if (have_SSAAttenuation) {
          std::cout << "have SSAAttenuation" << std::endl;
          syn.SSAOpticalDepth(s.param.getEmissionRadius(),
      s.param.getMagStrength()); syn.Attenuation(p[p_id].u);
        }

        if (have_BremsAttenuation) {
          std::cout << "have BremsAttenuation" << std::endl;
          Br.FFOpticalDepth(s.param.getEmissionRadius(),
      s.param.getAmbientDensity()); Br.Attenuation(p[p_id].u);
        }

      if (have_EBLAttenuation) {
        if (p_id == 22) {
          std::cout << "have EBLAttenuation" << std::endl;
          Photonbackground ph(s);
          ph.EBLTau(s.param.getRedshift(), ph.EBL_Gilmore12());
          ph.EBLAttenuation(p[p_id].x_e, p[p_id].u_middle);
        }
      }
        */

        LorentzBoost(particle.x_e, particle.u_middle);
        FluxNorm(particle.x_e, particle.u_middle);

        mystream << particle.ID << " " << particle.Z << " " << particle.N;
        for (size_t i = 0; i < particle.x.size(); i++) {
            mystream << " " << particle.u_middle[i];
        }
        mystream << std::endl;
    }
    mystream.close();
}

