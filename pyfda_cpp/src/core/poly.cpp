#include "poly.hpp"

#include <algorithm>
#include <cmath>

namespace pyfda {

CVec poly(const CVec &r) {
    CVec c{1.0};
    for (const cplx &root : r) {
        c.push_back(0.0);
        for (size_t i = c.size() - 1; i > 0; --i) c[i] -= root * c[i - 1];
    }
    return c;
}

Vec poly_real(const CVec &r) {
    const CVec c = poly(r);
    Vec out(c.size());
    for (size_t i = 0; i < c.size(); ++i) out[i] = c[i].real();
    return out;
}

Vec polymul(const Vec &a, const Vec &b) {
    if (a.empty() || b.empty()) return {};
    Vec c(a.size() + b.size() - 1, 0.0);
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
    return c;
}

cplx polyval(const Vec &c, cplx x) {
    cplx y = 0.0;
    for (double v : c) y = y * x + v;
    return y;
}

namespace {
// p(x) / p'(x), evaluated with the reversed polynomial for |x| > 1 to avoid overflow
cplx newton_ratio(const Vec &c, cplx x) {
    const int n = int(c.size()) - 1;
    if (std::abs(x) <= 1.0) {
        cplx p = c[0], dp = 0.0;
        for (int i = 1; i <= n; ++i) {
            dp = dp * x + p;
            p = p * x + c[i];
        }
        if (dp == 0.0) return p == 0.0 ? 0.0 : cplx(1e-10, 0);
        return p / dp;
    }
    const cplx y = 1.0 / x;
    // q(y) = sum c[i] y^i  (reversed coefficients)
    cplx q = c[n], dq = 0.0;
    for (int i = n - 1; i >= 0; --i) {
        dq = dq * y + q;
        q = q * y + c[i];
    }
    const cplx den = double(n) / x - dq / (x * x * q);
    if (q == 0.0) return 0.0;
    return 1.0 / den;
}
}  // namespace

CVec roots(const Vec &coeffs_in) {
    // strip leading zeros (highest powers) and count trailing zeros (roots at 0)
    size_t first = 0;
    while (first < coeffs_in.size() && coeffs_in[first] == 0.0) ++first;
    if (first == coeffs_in.size()) return {};
    size_t last = coeffs_in.size();
    while (last > first && coeffs_in[last - 1] == 0.0) --last;
    const size_t n_zero_roots = coeffs_in.size() - last;
    Vec c(coeffs_in.begin() + first, coeffs_in.begin() + last);
    const double c0 = c[0];
    for (double &v : c) v /= c0;
    const int n = int(c.size()) - 1;

    CVec x(n);
    if (n == 1) {
        x[0] = -c[1];
    } else if (n == 2) {
        const cplx disc = std::sqrt(cplx(c[1] * c[1] - 4.0 * c[2]));
        const cplx q = -0.5 * (c[1] + (c[1] >= 0 ? disc : -disc));
        x[0] = q;
        x[1] = (q != 0.0) ? cplx(c[2]) / q : cplx(0.0);
    } else if (n > 2) {
        // initial guesses on a circle with the geometric mean radius of the roots
        const double radius = std::pow(std::fabs(c[n]), 1.0 / n);
        const double r0 = radius > 0 ? radius : 1.0;
        for (int k = 0; k < n; ++k) {
            const double ang = 2.0 * PI * k / n + 0.4;
            x[k] = std::polar(r0, ang);
        }
        for (int iter = 0; iter < 1000; ++iter) {
            double max_rel = 0.0;
            for (int k = 0; k < n; ++k) {
                const cplx ratio = newton_ratio(c, x[k]);
                cplx s = 0.0;
                for (int j = 0; j < n; ++j)
                    if (j != k) s += 1.0 / (x[k] - x[j]);
                const cplx w = ratio / (1.0 - ratio * s);
                if (std::isfinite(w.real()) && std::isfinite(w.imag())) {
                    x[k] -= w;
                    max_rel = std::max(max_rel, std::abs(w) / std::max(std::abs(x[k]), 1e-300));
                }
            }
            if (max_rel < 1e-15) break;
        }
        // clean up tiny imaginary parts of real roots
        for (cplx &v : x)
            if (std::fabs(v.imag()) <= 1e-14 * std::abs(v)) v = v.real();
    }
    x.insert(x.end(), n_zero_roots, cplx(0.0));
    return x;
}

}  // namespace pyfda
