// AMES code
// Copyright(C) 2022, Bing Theodore Zhang <bing.zhang@yukawa.kyoto-u.ac.jp> and
// Kohta Murase <murase@psu.edu> licensed under the GNU GENERAL PUBLIC LICENSE,
// see LICENSE file for details

// \file utility.hpp
// \brief Useful tools
#include <pybind11/pybind11.h>
#include <csignal>

#include "module/utility.hpp"

#include <string.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "module/constants.hpp"

void Utility::ResizeVector2d(std::vector<std::vector<double>> &v, int nrow, int ncol) {
  v.resize(nrow);
  for (int i = 0; i < nrow; ++i)
    v[i].resize(ncol);
}

int Utility::FindCloset(const std::vector<double> &x, double value) {
  auto x_it = std::upper_bound(x.begin(), x.end(), value);
  int idx = x_it - x.begin() - 1;
  return idx;
}

double Utility::Integrate(const std::vector<double>& x, const std::vector<double>& y) {
    double sum = 0.0;
    const size_t n = x.size();
     #pragma omp parallel for reduction(+:sum) 
    for (size_t i = 0; i < n - 1; i++) {
        sum += (y[i+1] + y[i]) * (x[i+1] - x[i]);
    }
    const double ret = sum * 0.5;
    assert(ret >= 0.0 && "Integration result must be non-negative");
    return ret;
}

template<size_t N>
double Utility::Integrate(const std::array<double, N>& x, const std::array<double, N>& y) {
    double sum = 0.0;

    #pragma omp parallel for reduction(+:sum) 
    for (size_t i = 0; i < N - 1; i++) {
        sum += (y[i+1] + y[i]) * (x[i+1] - x[i]);
    }
    const double ret = sum * 0.5;
    assert(ret >= 0.0 && "Integration result must be non-negative");
    return ret;
}

template double Utility::Integrate<121>(const std::array<double, 121>&, const std::array<double, 121>&);

double Utility::Integrate(const std::vector<double> &x, const std::vector<double> &y,
                          const double x_begin) {
  // trapezoid rule for integration of discrete data
  // begining at x_begin
  double sum = 0;

  #pragma omp parallel for reduction(+:sum)  
  for (size_t i = 0; i < x.size() - 1; i++)
    if (x[i] >= x_begin) {
      sum += (y[i + 1] + y[i]) * (x[i + 1] - x[i]);
    }
  sum *= 0.5;
  return sum;
}

/*
double Utility::Integrate(const std::vector<double> &x, const std::vector<double> &y,
                          const double x_begin, const double x_end) {
  // trapezoid rule for integration of discrete data
  // begining at x_begin
  double sum = 0;
  for (size_t i = 0; i < x.size() - 1; i++)
    if ((x[i] >= x_begin) && (x[i] <= x_end)) {
      sum += (y[i + 1] + y[i]) * (x[i + 1] - x[i]);
    }
  sum *= 0.5;
  return sum;
}
*/

double Utility::Integrate(const std::vector<double>& x, 
                         const std::vector<double>& y,
                         double x_begin, 
                         double x_end) 
{
    // Input validation
    if (x.empty() || y.empty() || x.size() != y.size()) {
        return 0.0;
    }

    // Find the integration bounds in the vector
    auto lower = std::lower_bound(x.begin(), x.end(), x_begin);
    auto upper = std::upper_bound(x.begin(), x.end(), x_end);
    
    // Calculate the valid range
    size_t start_idx = std::distance(x.begin(), lower);
    size_t end_idx = std::distance(x.begin(), upper);
    
    // Handle edge cases
    if (start_idx >= end_idx || end_idx >= x.size()) {
        return 0.0;
    }

    double sum = 0.0;
    const size_t n = end_idx - start_idx;

    // Parallel reduction with OpenMP
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (size_t i = start_idx; i < end_idx - 1; ++i) {
        double dx = x[i+1] - x[i];
        double avg_y = (y[i+1] + y[i]) * 0.5;
        sum += dx * avg_y;
    }

    return sum;
}

/*
double Utility::Interpolate(const std::vector<double> &x, const std::vector<double> &y,
                            double value) {
  const unsigned int size = x.size();

  // linear interpolate
  if (value < x[0]) {
    return 0.;
  }

  if (value > x[size - 1]) {
    return 0.;
  }

  // vector<double>::iterator x_it;
  auto x_it = upper_bound(x.begin(), x.end(), value);
  int idx = x_it - x.begin() - 1;
  const double x1 = x[idx];
  const double y1 = y[idx];
  if (x1 == value)
    return y1;
  ++idx;
  const double x2 = x[idx];
  const double y2 = y[idx];
  return y1 + (y2 - y1) * (value - x1) / (x2 - x1);
}*/

double Utility::Interpolate(const std::vector<double> &x, const std::vector<double> &y, double value) {
  const unsigned int size = x.size();

  if (size == 0 || y.size() != size) {
    std::cerr << "Error: x and y sizes do not match or zero size.\n";
    return 0.;
  }

  if (value < x.front()) {
    return 0.;
  }

  if (value > x.back()) {
    return 0.;
  }

  // Use lower_bound instead to avoid negative index issues
  auto x_it = std::lower_bound(x.begin(), x.end(), value);
  int idx = x_it - x.begin();

  // value == x[0] or smaller, return y[0]
  if (x[idx] == value){
      assert(y[idx] >= 0.0 && "Interpolation result must be non-negative");
      return y[idx];}

  // idx points to first x >= value
  // interpolate between x[idx-1] and x[idx]
  double x1 = x[idx - 1];
  double y1 = y[idx - 1];
  double x2 = x[idx];
  double y2 = y[idx];

  double denom = (x2 - x1);
  if (denom == 0) {
    if (value > x.back()) {
        return 0.0;
    }
    //std::cerr << "Warning: zero denominator in interpolation between "
    //          << x1 << " and " << x2 << "\n";
    assert(y2 >= 0.0 && "Interpolation result must be non-negative");
    return y2; // or y2, or average, depending on intent
  }
  double ret = y1 + (y2 - y1) * (value - x1) / denom;

  assert(ret >= 0.0 && "Interpolation result must be non-negative");
  return ret;
}

double Utility::Interpolate(const std::vector<double> &x, const std::vector<double> &y,
                            double value, int idx) {
  const unsigned int size = x.size();

  // linear interpolate
  if (value < x[0]) {
    return 0.;
  }

  if (value > x[size - 1]) {
    return 0.;
  }

  const double x1 = x[idx];
  const double y1 = y[idx];
  if (x1 == value){
      assert(y1 >= 0.0 && "Interpolation result must be non-negative");
      return y1;
  }
  ++idx;
  const double x2 = x[idx];
  const double y2 = y[idx];
  return y1 + (y2 - y1) * (value - x1) / (x2 - x1);
}

double Utility::Interpolate2D(const std::vector<double> &x, const std::vector<double> &y,
                              const std::vector<std::vector<double>> &xy, double value_x,
                              double value_y) {
  const unsigned int size_x = x.size();
  const unsigned int size_y = y.size();

  // linear interpolate
  if ((value_x < x[0])) {
    return 0.;
  }

  if (value_x > x[size_x - 1]) {
    return 0.;
  }

  if (value_y < y[0]) {
    return 0.;
  }

  if (value_y > y[size_y - 1]) {
    return 0.;
  }

  auto x_it = std::upper_bound(x.begin(), x.end(), value_x);
  int idx_x = x_it - x.begin() - 1;
  auto y_it = std::upper_bound(y.begin(), y.end(), value_y);
  int idx_y = y_it - y.begin() - 1;
  const double x1 = x[idx_x];
  const double y1 = y[idx_y];
  /*
  const double xy1 = xy[idx_x][idx_y];
  if ((x1 == value_x) && (y1 == value_y)) {
    return xy1;
  }
  */

  if ((idx_x == size_x - 1) && (idx_y <= size_y - 1)) {
    int idx_y2 = idx_y + 1;
    const double y2 = y[idx_y2];
    return ((y2 - value_y) * xy[idx_x][idx_y] + (value_y - y1) * xy[idx_x][idx_y2]) / (y2 - y1);
  }

  if ((idx_y == size_y - 1) && (idx_x <= size_x - 1)) {
    int idx_x2 = idx_x + 1;
    const double x2 = x[idx_x2];
    return ((x2 - value_x) * xy[idx_x][idx_y] + (value_x - x1) * xy[idx_x2][idx_y]) / (x2 - x1);
  }


  int idx_x2 = idx_x + 1;
  int idx_y2 = idx_y + 1;
  const double x2 = x[idx_x2];
  const double y2 = y[idx_y2];
  double w11 = (x2 - value_x) * (y2 - value_y);
  double w12 = (x2 - value_x) * (value_y - y1);
  double w21 = (value_x - x1) * (y2 - value_y);
  double w22 = (value_x - x1) * (value_y - y1);

  return (w11 * xy[idx_x][idx_y] + w12 * xy[idx_x][idx_y2] + w21 * xy[idx_x2][idx_y] +
          w22 * xy[idx_x2][idx_y2]) /
         (x2 - x1) / (y2 - y1);
}

double Utility::Interpolate2D(const std::vector<double> &x, const std::vector<double> &y,
                              const std::vector<std::vector<double>> &xy, double value_x,
                              double value_y, int idx_x, int idx_y) {
  const unsigned int size = x.size();

  // linear interpolate
  if ((value_x < x[0])) {
    return 0.;
  }

  if (value_x > x[size - 2]) {
    return 0.;
  }

  if (value_y < y[0]) {
    return 0.;
  }

  if (value_y > y[size - 2]) {
    return 0.;
  }


  const double x1 = x[idx_x];
  const double y1 = y[idx_y];
  const double xy1 = xy[idx_x][idx_y];
  if ((x1 == value_x) && (y1 == value_y)) {
    return xy1;
  }
  int idx_x2 = idx_x + 1;
  int idx_y2 = idx_y + 1;
  const double x2 = x[idx_x2];
  const double y2 = y[idx_y2];
  const double xy2 = xy[idx_x2][idx_y];
  double w11 = (x2 - value_x) * (y2 - value_y);
  double w12 = (x2 - value_x) * (value_y - y1);
  double w21 = (value_x - x1) * (y2 - value_y);
  double w22 = (value_x - x1) * (value_y - y1);
  return (w11 * xy[idx_x][idx_y] + w12 * xy[idx_x][idx_y2] + w21 * xy[idx_x2][idx_y] +
          w22 * xy[idx_x2][idx_y2]) /
         (x2 - x1) / (y2 - y1);
  // return xy1 + (xy2 - xy1) * (value_x - x1) / (x2 - x1) +
  //        (xy2 - xy1) * (value_y - y1) / (y2 - y1);
}

void Utility::TDMA(const std::vector<double> &a, const std::vector<double> &b,
                   const std::vector<double> &c, const std::vector<double> &d,
                   std::vector<double> &u) {
  size_t N = d.size();
  std::vector<double> c_star(N, 0.0);
  std::vector<double> d_star(N, 0.0);

  c_star[0] = c[0] / b[0];
  d_star[0] = d[0] / b[0];

  for (int i = 1; i < N; i++) {
    double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
    c_star[i] = c[i] * m;
    d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
  }

  u[N - 1] = d_star[N - 1];
  for (int i = N - 2; i-- > 0;) {
    u[i] = d_star[i] - c_star[i] * d[i + 1];
  }
}

bool Utility::ReadFile1D(std::string file_name, std::vector<double> &raw) {
  raw.clear();
  std::ifstream in(file_name.c_str());
  // check if object is valid
  if (!in) {
    std::cerr << "Cannot open the file :  << " << file_name << std::endl;
    return false;
  }

  std::string temp_str;
  while (std::getline(in, temp_str)) {
    if (temp_str[0] == '#')
      continue;
    raw.push_back(std::stod(temp_str));
  }
  return true;
}

size_t Utility::ReadFile2D(const std::string& file_name, std::vector<double>& flat_data) {
    std::ifstream in(file_name);
    if (!in) {
        std::cerr << "Cannot open the file: " << file_name << std::endl;
        return 0;
    }

    std::string line;
    std::vector<double> temp_row;
    size_t num_cols = 0;
    flat_data.clear();

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        temp_row.clear();
        double value;
        while (iss >> value) {
            temp_row.push_back(value);
        }

        if (num_cols == 0) {
            num_cols = temp_row.size();
        } else if (temp_row.size() != num_cols) {
            std::cerr << "Inconsistent number of columns in file: " << file_name << std::endl;
            return 0;
        }

        flat_data.insert(flat_data.end(), temp_row.begin(), temp_row.end());
    }

    return num_cols;  // Caller can compute row count as flat_data.size() / num_cols
}


bool Utility::ReadFile2D(std::string file_name, std::vector<std::vector<double>> &raw) {
  std::ifstream in(file_name.c_str());
  // check if object is valid
  if (!in) {
    std::cerr << "Cannot open the file :  << " << file_name << std::endl;
    return false;
  }
  std::string temp_str;
  double temp_value;
  while (std::getline(in, temp_str)) {
    if (temp_str[0] == '#')
      continue;
    std::istringstream iss(temp_str);
    std::vector<double> tempv;
    while (iss >> temp_value) {
      tempv.push_back(temp_value);
    }
    raw.push_back(tempv);
    tempv.clear();
  }
  return true;
}

bool Utility::ReadFile2D(std::string file_name, std::vector<std::vector<int>> &raw) {
  std::ifstream in(file_name.c_str());
  // check if object is valid
  if (!in) {
    std::cerr << "Cannot open the file :  << " << file_name << std::endl;
    return false;
  }

  std::string temp_str;
  double temp_value;
  while (std::getline(in, temp_str)) {
    if (temp_str[0] == '#')
      continue;
    std::istringstream iss(temp_str);
    std::vector<int> tempv;
    while (iss >> temp_value) {
      tempv.push_back(temp_value);
    }
    raw.push_back(tempv);
    tempv.clear();
  }
  return true;
}

Cosmology::Cosmology() {}

double Cosmology::Redshift2LuminosityDistance(double z) {
  return (1 + z) * Redshift2ComovingDistance(z);
}

double Cosmology::Redshift2ComovingDistance(double z) {
  return HubbleDistance() * integrate(ComovingDistanceInte, 0.0, z, 30, simpson());
}

double Cosmology::Redshift2LighttravelDistance(double z) {
  return HubbleDistance() * integrate(LighttravelDistanceInte, 0.0, z, 30, simpson());
}

double Cosmology::LighttravelDistanceInte(double z) { return 1. / (1 + z) / Ez(z); }

double Cosmology::ComovingDistanceInte(double z) { return 1. / Ez(z); }

double Cosmology::Ez(double z) { return sqrt(0.315 * (1 + z) * (1 + z) * (1 + z) + 0.685); }

double Cosmology::HubbleDistance() { return 1e-5 * c_cnst / H0 * Mpc; }

