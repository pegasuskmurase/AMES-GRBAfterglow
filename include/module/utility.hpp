//========================================================================================
// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file utility.hpp
// \brief Useful tools.

#ifndef UTILITY_H
#define UTILITY_H

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <string>
#include <vector>

class Utility {
public:
  Utility() {}
  void ResizeVector2d(std::vector<std::vector<double>> &v, int nrow, int ncol);
  int FindCloset(const std::vector<double> &x, double value);
  double Integrate(const std::vector<double> &x, const std::vector<double> &y);
  double Integrate(const std::vector<double> &x, const std::vector<double> &y,
                   const double x_begin);
  double Integrate(const std::vector<double> &x, const std::vector<double> &y, const double x_begin,
                   const double x_end);
  double Interpolate(const std::vector<double> &x, const std::vector<double> &y, double value);
  double Interpolate(const std::vector<double> &x, const std::vector<double> &y, double value,
                     int idx);
  double Interpolate2D(const std::vector<double> &x, const std::vector<double> &y,
                       const std::vector<std::vector<double>> &xy, double value_x, double value_y);
  double Interpolate2D(const std::vector<double> &x, const std::vector<double> &y,
                       const std::vector<std::vector<double>> &xy, double value_x, double value_y,
                       int idx_x, int idx_y);

  bool solveQuadratic(const double &a, const double &b, const double &c, double &x0,
                      double &x1) const {
    double discr = b * b - 4 * a * c;
    if (discr < 0)
      return false;
    else if (discr == 0) {
      x0 = x1 = -0.5 * b / a;
    } else {
      float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
      x0 = q / a;
      x1 = c / q;
    }

    return true;
  }

  void TDMA(const std::vector<double> &a, const std::vector<double> &b,
            const std::vector<double> &c, const std::vector<double> &d, std::vector<double> &f);

  bool ReadFile2D(std::string file_name, std::vector<std::vector<double>> &raw);
  bool ReadFile2D(std::string file_name, std::vector<std::vector<int>> &raw);
  bool ReadFile1D(std::string file_name, std::vector<double> &raw);
  /*
  double PhotohadronicInterp(double E_proton, double E_photon,
                             double E_secondary, int N_proton, int N_photon,
                             int N_secondary, std::vector<double> &E_proton_vec,
                             std::vector<double> &E_photon_vec,
                             std::vector<double> &E_secondary_vec,
                             std::vector<std::vector<double>> &raw_table);
  */

  void Output();
};

struct ParticleTransport {
  int ID;
  int Z;
  int N;
  double mass;
  double life;
  std::vector<double> x;
  std::vector<double> dx;
  std::vector<double> x_e;
  std::vector<double> dx_e;
  std::vector<double> u;
  std::vector<double> u_middle;
  std::vector<double> u_next;
  std::vector<double> A;
  std::vector<double> A_dt;
  std::vector<std::vector<double>> B;
  std::vector<std::vector<std::vector<double>>> C;
  std::vector<int> C_id;
  std::vector<double> Q;
};

struct ParticleTransportChangCooper {
  int ID;
  int Z;
  int N;
  double mass;
  double life;
  std::vector<double> x;
  std::vector<double> dx;
  std::vector<double> x_e;
  std::vector<double> dx_e;
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;
  std::vector<double> d;
  std::vector<double> u;
  std::vector<double> u_middle;
  std::vector<double> a_m;
  std::vector<double> c_m;
  std::vector<double> a_p;
  std::vector<double> c_p;
  std::vector<double> Q; // injection term
  std::vector<double> B; // advective term and cooling term
  std::vector<double> D; // diffusive term
  std::vector<double> T; // escape term
  std::vector<double> A_dt;
  std::vector<double> xghost;
  std::vector<double> xghost_e;
  std::vector<double> Bghost;
  std::vector<double> Dghost;
  std::vector<int> C_id;
  std::vector<std::vector<std::vector<double>>> C;
};

template <typename T> class Vector3d {
public:
  Vector3d() : x(T(0)), y(T(0)), z(T(0)) {}
  Vector3d(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
  Vector3d operator-(const Vector3d &v) const { return Vector3d(x - v.x, y - v.y, z - v.z); }
  Vector3d operator+(const Vector3d &v) const { return Vector3d(x + v.x, y + v.y, z + v.z); }
  Vector3d operator=(const Vector3d &v) const { return Vector3d(v.x, v.y, v.z); }
  T DotProduct(const Vector3d<T> &v) const { return x * v.x + y * v.y + z * v.z; }
  T Norm() const { return x * x + y * y + z * z; }
  Vector3d Product(double f) const { return Vector3d(f * x, f * y, f * z); }
  T Length() const { return sqrt(Norm()); }
  const T &operator[](uint8_t i) const { return (&x)[i]; }
  T &operator[](uint8_t i) { return (&x)[i]; }

  void setValue(T _x, T _y, T _z) {
    x = _x;
    y = _y;
    z = _z;
  }

  Vector3d &normalize() {
    T n = Norm();
    if (n > 0) {
      T factor = 1.0 / sqrt(n);
      x *= factor, y *= factor, z *= factor;
    }

    return *this;
  }

  T x, y, z;
};

class Cosmology {
public:
  Cosmology();
  double Redshift2LuminosityDistance(double z);
  double Redshift2ComovingDistance(double z);
  double Redshift2LighttravelDistance(double z);
  double HubbleDistance();
  static double Ez(double z);

private:
  static double ComovingDistanceInte(double z);
  static double LighttravelDistanceInte(double z);

  double H0 = 67.4; // km s^-1 Mpc^-1
  double Omega_M = 0.315;
  double Omega_L = 0.685;
  // double H0 = 73; // km s^-1 Mpc^-1
  // double Omega_M = 0.3;
  // double Omega_L = 0.7;
};

// methods
class simpson {
public:
  template <typename F, typename Float> double operator()(F f, Float x, Float h) const {
    return (f(x) + 4 * f(x + h / 2) + f(x + h)) / 6;
  }
};

template <typename Method, typename F, typename Float>
double integrate(F f, Float a, Float b, int steps, Method m) {
  double s = 0;
  double h = (b - a) / steps;
  for (int i = 0; i < steps; ++i)
    s += m(f, a + h * i, h);
  return h * s;
}

#endif
