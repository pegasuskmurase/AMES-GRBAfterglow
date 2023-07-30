// constants.hpp
#ifndef CONSTANTS_hpp
#define CONSTANTS_hpp

const double PI = 3.14159265358979323846;

const double e_mass = 510998;          // electron mass(eV)
const double pion_mass = 0.13957018e9; // pi+-mass(eV)
const double proton_mass = 938.271998e6;
const double neutron_mass = 939.56542052e6;
const double muon_mass = 0.105658357e9;
const double pic_Mass = 0.13957018e9;            // pi+-mass(GeV)
const double piz_Mass = 0.1349766e9;             // pi0
const double mu_Mass = 0.105658357e9;

const double eV2erg = 1.6021e-12; // eV2erg
const double erg2eV = 6.24152e11; // erg2eV
const double sigmaT = 6.652e-25;  // Thomson [cm^2]
const double mubar = 1.0e-30;     // mubar [cm^2]
const double mbar = 1.0e-27;      // mbar [cm^2]

const double u_cnst = 0.931494043;                 // unit atom mass
const double c_cnst = 2.99792458e10;               // light velocity
const double hbar_cnst = 1.054571596e-27;          // Planck constant [erg s]
const double e_cnst = 1.60217462e-12;              // charge
const double e_charge = 4.80320420e-10;            // charge
const double T_c_cnst = 6.652 * 2.99792458 / 1e15; // Thomson*light velocity(cm3/sec)
const double sigma_T = 6.652e-25;                   // Thomson*light velocity(cm3/sec)
const double eMass = 0.000510998902;               // electron mass(GeV)
const double cleradius = e_charge * e_charge / eMass / 1e9 / e_cnst;
const double electron_radius = 2.81794032e-13; // electron radius [cm]
const double pMass = 0.938271998;              // proton mass(GeV)
const double mp = 1.67262e-24;                 // proton mass (g)
const double me = 9.109382e-28;                // proton mass (g)
const double neutron_Life = 879.6;             // free neutron decay (s), mean life time
const double pic_Life = 2.6033 / 1e8;          // pi+-life(sec)
const double piz_Life = 8.4 / 1e17;            //
const double kac_Mass = 0.493677 * 1e9;              // ka+-
const double kac_Life = 1.2384 / 1e8;          //
const double kaz_Mass = 0.497672 * 1e9;              // ka0
const double kas_Life = 0.8935 / 1e10;         // kas
const double kal_Life = 5.17 / 1e8;            // kal
const double mu_Life = 2.19703 / 1e6;
const double me_eV = 511998; // eV

const double Kpc = 3.0856775807 * 1e21;
const double Mpc = 3.0856775807 * 1e24;
const double Gpc = 3.0856775807 * 1e27;
const double pc = 3.0856775807e18;
const double yr = 3.1536 * 1e7;
const double days = 3600. * 25;

const double eV2Hz = 2.4e14; //
const double eV2ang = 12291.472;
const double Hz2eV = 4.16e-15; //
const double erg2Jy = 1e23;
const double Jy2erg = 1e-23;

const double CMBtemp = 2.725 * 1.3806503 * 1e-16 / e_cnst;
const double h_Planck = 4.13e-15;     // eV s
const double k_B = 8.617342294984e-5; // [eV/K] Boltzman constant
const double T_CMB = 2.72548;         // CMB temperature [K]

const double M_sun = 1.988e33; // g
const double R_sun = 6.957e10; // cm

#endif
