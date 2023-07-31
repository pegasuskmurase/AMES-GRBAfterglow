%module AMES
%{
#include <setjmp.h>
#include <signal.h>

static sigjmp_buf timeout;

static void backout(int sig) {
  siglongjmp(timeout, sig);
}

#include "module/particle.hpp"
#include "module/constants.hpp"
#include "module/source.hpp"
#include "module/photonbackground.hpp"
#include "module/synchrotron.hpp"
#include "module/inversecompton.hpp"
#include "module/gammagamma.hpp"
#include "module/electrondistribution.hpp"
#include "GRBAfterglow.hpp"
%}

// Import standard types
%include "std_string.i"
%include "std_vector.i"
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(VecVecdouble) std::vector< std::vector<double> >;


%include <exception.i>

%exception {
  if (!sigsetjmp(timeout, 1)) {
    signal(SIGINT,backout); // Check return?
    $action
  }
  else {
    // raise a Python exception
    SWIG_exception(SWIG_RuntimeError, "Timeout in $decl");
  }
}


#include "module/constants.hpp"
#include "module/particle.hpp"
#include "module/source.hpp"
#include "module/photonbackground.hpp"
#include "module/synchrotron.hpp"
#include "module/inversecompton.hpp"
#include "module/gammagamma.hpp"
#include "module/adiabatic.hpp"
#include "module/electrondistribution.hpp"
#include "GRBAfterglow.hpp"

class Particle {
  friend class RadiativeProcess;
  friend class CRDistribution;
  friend class ElectronDistribution;
  friend class Photonbackground;
  friend class Synchrotron;
  friend class InverseCompton;
  friend class Bremsstrahlung;
  friend class GammaGamma;
  friend class Protonproton;
  friend class Photomeson;
  friend class Photopair;
  friend class Photodisintegration;

public:
  Particle();
  Particle(size_t p_num, double p_min, double p_max)
      : p_num(p_num), p_min(p_min), p_max(p_max) {
    setMomentum();
  }
  Particle(Particle &p);
  ~Particle();

  int getID() { return ID; }
  std::string getName() { return name; }
  double getMass() {return mass;}
  double getChargeNumber() {return Z;}
  double getMassNumber() {return A;}
  std::vector<double> getMomentum() { return momentum; }
  std::vector<double> getMomentumbin() { return momentumbin; }
  std::vector<double> getMomentumghost() { return momentumghost; }
  std::vector<double> getEnergy() { return energy; }
  std::vector<double> getEnergybin() { return energybin; }
  std::vector<double> getSpectrum() { return spectrum; }
  double getSpectrum(double value) {
    return utility.Interpolate(momentum, spectrum, value);
  }

  std::vector<double> getLosstime() { return losstime; }
  std::vector<double> getIntetime() { return intetime; }
  std::vector<double> getOptdepth() { return optdepth; }
  std::vector<double> getTemp() { return temp; }

  double getSpectraInj() { return p_inj; }
  double getSpectraCut() { return p_cut; }
  double getSpectraIndex() { return p_index; }

  int getSize() { return momentum.size(); }

  void setID(int _ID);
  void setName(std::string _name);
  void setMass(double _mass);
  void setCharge(double _Z);
  void setMassNumber(double _A);

  void setMomentum();
  void setMomentum(size_t _p_num, double _p_min, double _p_max);
  void setMomentum(int i, double value) { momentum[i] = value; }
  double getIndex(double value);
  void setEnergy(size_t _p_num, double _p_min, double _p_max);

  void setSpectrum(std::vector<double> vec) { spectrum = vec; }
  void setSpectrum(int i, double value) { spectrum[i] = value; }
  void setSpectrumPL(double p_inj, double p_cut, double index, double norm,
                     bool norm_type);
  void setSpectrumPLwExp(double _p_inj, double _p_cut, double _index, double norm,
                             bool norm_type);
  void setSpectrumPLNoExp(double p_inj, double p_cut, double index, double norm,
                          bool norm_type);
  void setSpectrumBPL(double p_inj, double p_break, double p_cut,
                      double index_l, double index_h, double norm,
                      bool norm_type);
  void setSpectrumBPLNoExp(double p_inj, double p_break, double p_cut,
                      double index_l, double index_h, double norm,
                      bool norm_type);
  void setSpectrumDelta(double mu, double a, double norm);
  void setSpectrumDelta(double mu, double a, double norm, bool norm_type);
  void setSpectrumMono(double value);
  void setEnergySpectrumPL(double _p_inj, double _p_cut, double _index, double norm,
                           bool norm_type);
  void setEnergySpectrumBPL(double p_inj, double p_break, double p_cut, double index_l, double index_h,
                      double norm, bool norm_type);
  void clear();


  void setTemp(int i, double value) { temp[i] = value; }

  void setLosstime(int i, double value) { losstime[i] = value; }
  void setOptdepth(int i, double value) { optdepth[i] = value; }

  double Momentum2Energy(double momentum_x);
  double Energy2Momentum(double energy_x);
  void Momentum2Energy();
  void Energy2Momentum();

  double Interpolate(double value) {
    return utility.Interpolate(energy, spectrum, value);
  }

  double InterpolateOptdepth(double value) {
    return utility.Interpolate(energy, optdepth, value);
  }

  double Interpolate(const std::vector<double> x, const std::vector<double> y,
                     double value) {
    return utility.Interpolate(x, y, value);
  }
};

class Param {
public:
  void setName(std::string _name);
  void setGeometry(std::string _geometry);
  void setRedshift(double _redshift); 
  void setLumDistance(double _lum_distance);
  void setMagStrength(double _mag_strength);
  void setEmissionRadius(double _emission_radius);
  void setEmissionVolume(double _emission_volume);
  void setAmbientDensity(double _ambient_density);
  void setDynamicalTime(double _dynamical_time);
  void setDopplerFactor(double _doppler_factor);
  void setLorentzFactor(double _lorentz_factor);
  void setObsAngle(double _obs_angle);
  void setVelocity(double _velocity);

  std::string getName();
  std::string getGeometry();
  double getRedshift();
  double getLumDistance();
  double getMagStrength();
  double getEmissionRadius();
  double getEmissionVolume();
  double getAmbientDensity();
  double getDopplerFactor();
  double getVelocity();
};

class Source
{
    public:
    Source();
    Source(Source& s);

    void setParticle(); //set particles, default in constructor
    void setParam(); //set parameters, default in constructor
    void showParam();

    Particle& getPhoton();    
    Particle& getTarget(); 
    Particle& getElectron(); 
    Param& getParam(); 
};

class Cosmology
{
    public:
    Cosmology();
    double Redshift2LuminosityDistance(double z);
    double Redshift2ComovingDistance(double z);
    double Redshift2LighttravelDistance(double z);
    double HubbleDistance();
    static double Ez(double z);
};

class Photonbackground
{
    public:
    Photonbackground(Source& s);
    void GreyBody(double T_ph, double u_ph);
    void BlackBody(double T_ph);
    void Powerlaw(double eps_min, double eps_max, double index, double u_ph);
    void BrokenPowerlaw(double eps_min, double eps_break, double eps_max, double index_1, double index_2, double u_ph);

    void CMB(double redshift);
    void EBL(double redshift, std::string EBL_model = "EBL_Gilmore12");
    void EBLTau(double redshift, std::string EBL_model = "EBL_Gilmore12");
    void EBLAttenuation(const std::vector<double> &photon_energy,
                      std::vector<double> &photon_spectrum);
    void EBL_Franceschini08(double redshift);
    void EBL_Gilmore12(double redshift);
    void EBL_Finke10(double redshift);
    std::string EBL_Franceschini08();
    std::string EBL_Gilmore12();
    std::string EBL_Finke10();
    void setOutputFile(std::string _output_file);
};

class Synchrotron {
  /*
  @class: Synchrotron 
  @brief: basic class for Synchrotron emission process
  */
public:
  Synchrotron(Source &s);
  void Losstime();
  void Losstime(const double mag_strength);
  void Losstime(const double mag_strength, const std::vector<double> &electron_energy,
                std::vector<double> &electron_losstime);
  // radiation emissivity in units [eV^-1 s^-1 cm^-3]
  void Emissivity();
  void Emissivity(double mag_strength);
  void Spec(const double mag_strength, const std::vector<double> &primary_energy,
            const std::vector<double> &secondary_energy);
  std::vector<std::vector<double>> spec;
  void Test();
};

class InverseCompton {
  /*
  @class: InverseCompton
  @brief: Class for Inverse Compton emission process
  */
public:
  InverseCompton(Source &s);
  // energy loss time in units [s^-1]
  void Losstime(); // Electron Inverse-Compton energy loss time
  void Losstime(const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &electron_energy, std::vector<double> &electron_losstime);

  // Interaction time in units [s^-1]
  void Intetime();

  void Spec(const std::vector<double> &primary_energy,
              const std::vector<double> &secondary_energy, const std::vector<double> &target_energy,
              const std::vector<double> &target_spectrum);
  void Emissivity(); // Inverse-Compton emissivity
  void Emissivity(const std::vector<double> &target_energy,
                  const std::vector<double> &target_spectrum);
  void Table();
  std::vector<std::vector<double>> spec;

  void Test();
};

class GammaGamma {
  /*
  @class: GammaGamma
  @brief: Class for gamma-gamma pair annihilation process
  */
public:
  GammaGamma(Source &s);

  // energy loss time in units [s^-1]
  void Losstime(); // Electron Inverse-Compton energy loss time
  void Losstime(const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &photon_energy, std::vector<double> &photon_losstime);

  // Interaction time in units [s^-1]
  void Intetime();
  void Intetime(const std::vector<double> &target_energy,
                const std::vector<double> &target_spectrum,
                const std::vector<double> &photon_energy, std::vector<double> &photon_intetime);

  void Test();
};

class ElectronDistribution {
public:
  ElectronDistribution(Source &s);

  void IterationSolutionSteadyState(Synchrotron &syn, InverseCompton &IC, bool have_EIC = false,
                         double t_dyn = 1e10, double t_esc_photon = 1e10);
  void IterationSolution(Synchrotron &syn, InverseCompton &IC, bool have_EIC = false,
                         double t_dyn = 1e10, double t_ad = 1e10, double t_esc_photon = 1e10);
  void Test();
};


class Utility
{
    public:
    Utility() {}
    double Integrate(const std::vector<double>& x, const std::vector<double>& y);
    double Integrate(const std::vector<double> &x, const std::vector<double> &y,
                   const double x_begin);
    double Interpolate(const std::vector<double>& x, const std::vector<double>& y, double value);
};

class GRBAfterglow {
public:
  GRBAfterglow(Source &s);

  void help();
  void setGRBAfterglowParam(std::vector<double> _param);

  void setOutputFolder(std::string _output_folder);
  void haveOneZone(bool _have_onezeon);                 // One zone
  void haveEdgeEffect(bool _have_edge_effect);          // jet break, edge effect
  void haveSSCSpec(bool _have_SSCSpec);                 // SSC
  void haveAttenuSSA(bool _have_attenu_SSA);            // SSA
  void haveAttenuFF(bool _have_attenu_FF);              // Free-free
  void haveAttenuGGSource(bool _have_attenu_GG_source); // gamma-gamma attenuation inside sources
  void haveAttenuGGCosmic(
      bool _have_attenu_GG_cosmic); // gamma-gamma attenuation in the intergalactic space
  void Flux(const std::vector<double> &time_array, const std::vector<double> &energy_array_min,
            const std::vector<double> &energy_array_max);
};


const double Kpc;
const double Mpc;
const double Gpc;
const double pc;

const double PI;
const double c_cnst;
const double e_cnst;
const double e_charge;
const double e_mass;
const double pion_mass;
const double proton_mass;
const double neutron_mass = 939.56542052e6;
const double neutron_Life = 879.6;             // free neutron decay (s), mean life time
const double sigmaT;
const double hbar_cnst;

const double yr;
const double days;

const double mp;
const double me;

const double eV2erg;
const double erg2eV;
const double eV2Hz; // 
const double Hz2eV; // 

const double erg2Jy = 1e23;
const double Jy2erg = 1e-23;

const double CMBtemp=2.725*1.3806503*1e-16/e_cnst;
const double h_Planck = 4.13e-15; // eV s
const double k_B = 8.617342294984e-5; // [eV/K] Boltzman constant
const double T_CMB = 2.72548;  // CMB temperature [K]
