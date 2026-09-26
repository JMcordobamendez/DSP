// Polynomial helpers (coefficients ordered with the highest power first, like numpy)
#pragma once

#include "types.hpp"

namespace pyfda {

/// Coefficients of the monic polynomial with the given roots (numpy.poly)
CVec poly(const CVec &roots);
/// Real part of poly(roots), for roots in complex conjugate pairs
Vec poly_real(const CVec &roots);

/// Roots of a polynomial (numpy.roots), using the Aberth-Ehrlich iteration
CVec roots(const Vec &coeffs);

/// Polynomial product (numpy.convolve)
Vec polymul(const Vec &a, const Vec &b);

/// Evaluate polynomial at x (Horner)
cplx polyval(const Vec &c, cplx x);

}  // namespace pyfda
