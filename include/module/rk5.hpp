#include <cmath>
#include <vector>

template <typename DerivedFunction> class rk5 {
public:
  rk5(DerivedFunction &derivs, std::vector<double> &y, std::vector<double> &dydx,
      std::vector<double> &yout, double x, const int n)
      : derivs(derivs), y(y), dydx(dydx), yout(yout), x(x), n(n) {
    k2.resize(n);
    k3.resize(n);
    k4.resize(n);
    k5.resize(n);
    k6.resize(n);
    rcont1.resize(n);
    rcont2.resize(n);
    rcont3.resize(n);
    rcont4.resize(n);
    rcont5.resize(n);
    dydxnew.resize(n);
    yerr.resize(n);
    ytemp.resize(n);
    EPS = std::numeric_limits<double>::epsilon();
    atol = 1e-13;
    rtol = 1e-13;
    dense = true;
  }

  std::vector<double> integrate(const double x2);
  std::vector<double> integrateEATS(const double x2, const double mu);
  std::vector<std::vector<double>> integrate(const std::vector<double> x2);
  void step();
  void dy(const double &h);
  void prepare_dense(const double h);

  double error();
  bool success(const double err, double &h);
  double dense_out(int i, const double x, const double h);

private:
  DerivedFunction &derivs;
  bool dense;
  double scale, safe, atol, rtol;
  double EPS;
  double hnext, errold;
  bool reject;

  std::vector<double> yerr, ytemp, dydxnew;
  std::vector<double> rcont1, rcont2, rcont3, rcont4, rcont5;
  std::vector<double> k2, k3, k4, k5, k6;

  double x;
  int n;
  std::vector<double> y;
  std::vector<double> dydx;
  std::vector<double> yout;

  double htry;
  double xold, hold;
};

template <typename DerivedFunction>
std::vector<std::vector<double>> rk5<DerivedFunction>::integrate(const std::vector<double> x2) {
  std::vector<std::vector<double>> yvalue;
  yvalue.resize(x2.size(), std::vector<double>(n + 1, 0.));
  htry = 1e-11;
  double xmax = x2.back();
  int idx = 0;
  double beta;
  double t_lab = 0.0;
  double y_pre = y[0];
  double dr;
  while (x < xmax) {
    step();
    beta = sqrt(1 - 1. / (y[2] * y[2]));
    dr = y[0] - y_pre;
    y_pre = y[0];
    t_lab += dr / (beta * c_cnst);
    while ((x2[idx] < xold + hold) && (x2[idx] > xold)) {
      for (size_t i = 0; i < n; i++) {
        yvalue[idx][i] = dense_out(i, x2[idx], xold);
      }
      yvalue[idx][n] = t_lab - (y[0] - yvalue[idx][0]) / (beta * c_cnst);
      idx++;
    }
  }
  return yvalue;
}

template <typename DerivedFunction>
std::vector<double> rk5<DerivedFunction>::integrate(const double x2) {
  std::vector<double> yvalue(n + 1);
  htry = 1e-11;
  double beta;
  double t_lab = 0.0;
  while (x < x2) {
    step();
    beta = sqrt(1 - 1. / (y[2] * y[2]));
    t_lab += htry / (1 - beta);
    if ((x2 < xold + hold) && (x2 > xold)) {
      for (size_t i = 0; i < n; i++) {
        yvalue[i] = dense_out(i, x2, xold);
      }
      yvalue[n] = t_lab - (x - x2) / (1 - beta);
    }
  }
  return yvalue;
}


template <typename DerivedFunction>
std::vector<double> rk5<DerivedFunction>::integrateEATS(const double T, const double mu) {
  std::vector<double> yvalue(n);
  htry = 1e-11;
  while (x < T + mu * y[0] / c_cnst) {
    step();
  }
  for (size_t i = 0; i < n; i++) {
    yvalue[i] = y[i];
  }

  if (abs(x - mu * y[0] / c_cnst - T) / T > 0.1) {
      std::cout << "Warning: " << T << " " <<  x - mu * y[0] / c_cnst << " " << (x - mu * y[0] / c_cnst - T) / T << " " << y[2] << " " << y[0] << " " << mu << " " << x << std::endl;
  }

  return yvalue;
}

template <class DerivedFunction> void rk5<DerivedFunction>::step() {
  double h = htry;
  for (;;) {
    dy(h);
    double err = error();
    if (success(err, h))
      break;
    if (std::abs(h) < std::abs(x) * EPS) {
      throw("stepsize underflow in rk5");
    }
  }
  if (dense)
    prepare_dense(h);
  dydx = dydxnew;
  y = yout;
  xold = x;
  x += h;
  hold = h;
  htry = hnext;
}

template <class DerivedFunction> void rk5<DerivedFunction>::dy(const double &h) {
  // Dormand-Prince parameters for Runga-kutta Method
  double c2 = 1. / 5., c3 = 3. / 10., c4 = 4. / 5., c5 = 8. / 9., c6 = 1., c7 = 1.;
  double a21 = 1. / 5., a31 = 3. / 40., a32 = 9. / 40., a41 = 44. / 45., a42 = -56. / 15.,
         a43 = 32. / 9., a51 = 19372. / 6561., a52 = -25360. / 2187., a53 = 64448. / 6561.,
         a54 = -212. / 729., a61 = 9017. / 3168, a62 = -355. / 33., a63 = 46732. / 5247.,
         a64 = 49. / 176., a65 = -5103. / 18656., a71 = 35. / 384., a72 = 0.0, a73 = 500. / 1113.,
         a74 = 125. / 192., a75 = -2187. / 6784., a76 = 11. / 84.;
  double e1 = 71.0 / 57600., e2 = 0.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0,
         e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
  for (size_t i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * a21 * dydx[i];
  }
  derivs(x + c2 * h, ytemp, k2); // 2
  for (size_t i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (a31 * dydx[i] + a32 * k2[i]);
  }
  derivs(x + c3 * h, ytemp, k3); // 3
  for (size_t i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (a41 * dydx[i] + a42 * k2[i] + a43 * k3[i]);
  }
  derivs(x + c4 * h, ytemp, k4); // 4
  for (size_t i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (a51 * dydx[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
  }
  derivs(x + c5 * h, ytemp, k5); // 5
  for (size_t i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (a61 * dydx[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
  }
  double xph = x + h;
  derivs(x + c6 * h, ytemp, k6); // 6
  for (size_t i = 0; i < n; i++) {
    yout[i] = y[i] + h * (a71 * dydx[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
  }
  derivs(xph, yout, dydxnew); // next step
  for (size_t i = 0; i < n; i++) {
    yerr[i] =
        h * (e1 * dydx[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * dydxnew[i]);
  }
}

template <class DerivedFunction> double rk5<DerivedFunction>::error() {
  double err = 0.0;
  double sk;
  for (size_t i = 0; i < n; i++) {
    sk = atol + rtol * std::max(std::abs(y[i]), std::abs(yout[i]));
    err += (yerr[i] / sk) * (yerr[i] / sk);
  }
  return sqrt(err / n);
}

template <class DerivedFunction> bool rk5<DerivedFunction>::success(const double err, double &h) {
  double beta = 0.0, alpha = 0.2 - beta * 0.75, safe = 0.9, minscale = 0.1, maxscale = 10.0;

  if (err < 1.0) {
    if (err == 0.0) {
      scale = maxscale;
    } else {
      scale = safe * pow(err, -alpha) * pow(errold, beta);
      if (scale < minscale)
        scale = minscale;
      if (scale > maxscale)
        scale = maxscale;
    }
    if (reject) {
      hnext = h * std::min(scale, 1.0);
    } else {
      hnext = h * scale;
    }
    errold = std::max(err, 1.0e-4);
    reject = false;
    return true;
  } else {
    scale = std::max(safe * pow(err, -alpha), minscale);
    h *= scale;
    reject = true;
    return false;
  }
}

template <class DerivedFunction> void rk5<DerivedFunction>::prepare_dense(const double h) {
  const double d1 = -12715105075.0 / 11282082432.0, d3 = 87487479700.0 / 32700410799.0,
               d4 = -10690763975.0 / 1880347072.0, d5 = 701980252875.0 / 199316789632.0,
               d6 = -1453857185.0 / 822651844.0, d7 = 69997945.0 / 29380423.0;
  for (size_t i = 0; i < n; i++) {
    rcont1[i] = y[i];
    double ydiff = yout[i] - y[i];
    rcont2[i] = ydiff;
    double bsp1 = h * dydx[i] - ydiff;
    rcont3[i] = bsp1;
    rcont4[i] = ydiff - h * dydxnew[i] - bsp1;
    rcont5[i] =
        h * (d1 * dydx[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * dydxnew[i]);
  }
}

template <class DerivedFunction>
double rk5<DerivedFunction>::dense_out(int i, const double x, const double h) {
  double s = (x - xold) / h;
  double s1 = 1.0 - s;
  return rcont1[i] + s * (rcont2[i] + s1 * (rcont3[i] + s * (rcont4[i] + s1 * rcont5[i])));
}
