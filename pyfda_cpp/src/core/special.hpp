// Special functions needed for filter design (ports of the scipy.special
// routines used by scipy.signal)
#pragma once

#include "types.hpp"

#include <functional>

namespace pyfda {

/// Complete elliptic integral of the first kind K(m), 0 <= m < 1
double ellipk(double m);
/// K(1 - p), accurate for small p
double ellipkm1(double p);

/// Jacobian elliptic functions sn, cn, dn and amplitude ph of u with parameter m
void ellipj(double u, double m, double &sn, double &cn, double &dn, double &ph);

/// Modified Bessel function of the first kind, order 0
double bessel_i0(double x);

/// Normalized sinc: sin(pi x) / (pi x)
double sinc(double x);

/// 10**x - 1 without cancellation
double pow10m1(double x);

/// Bounded scalar minimization (port of scipy.optimize.fminbound, xtol = 1e-5)
double fminbound(const std::function<double(double)> &f, double x1, double x2,
                 double xatol = 1e-5, int maxfun = 500);

}  // namespace pyfda
