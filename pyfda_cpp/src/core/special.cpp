#include "special.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pyfda {

namespace {
double agm(double a, double b) {
    for (int i = 0; i < 100; ++i) {
        const double an = 0.5 * (a + b);
        const double bn = std::sqrt(a * b);
        if (std::fabs(an - bn) <= 1e-17 * an) return an;
        a = an;
        b = bn;
    }
    return 0.5 * (a + b);
}
}  // namespace

double ellipk(double m) {
    if (m >= 1.0) return m == 1.0 ? std::numeric_limits<double>::infinity()
                                   : std::numeric_limits<double>::quiet_NaN();
    return PI / (2.0 * agm(1.0, std::sqrt(1.0 - m)));
}

double ellipkm1(double p) {
    if (p <= 0.0) return p == 0.0 ? std::numeric_limits<double>::infinity()
                                   : std::numeric_limits<double>::quiet_NaN();
    return PI / (2.0 * agm(1.0, std::sqrt(p)));
}

// Port of the Cephes routine ellpj used by scipy.special.ellipj
void ellipj(double u, double m, double &sn, double &cn, double &dn, double &ph) {
    constexpr double MACHEP = 1.11022302462515654042e-16;
    if (m < 0.0 || m > 1.0 || std::isnan(m)) {
        sn = cn = dn = ph = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    if (m < 1.0e-9) {
        const double t = std::sin(u);
        const double b = std::cos(u);
        const double ai = 0.25 * m * (u - t * b);
        sn = t - ai * b;
        cn = b + ai * t;
        ph = u - ai;
        dn = 1.0 - 0.5 * m * t * t;
        return;
    }
    if (m >= 0.9999999999) {
        double ai = 0.25 * (1.0 - m);
        const double b = std::cosh(u);
        const double t = std::tanh(u);
        const double phi = 1.0 / b;
        const double twon = b * std::sinh(u);
        sn = t + ai * (twon - u) / (b * b);
        ph = 2.0 * std::atan(std::exp(u)) - PI / 2 + ai * (twon - u) / b;
        ai *= t * phi;
        cn = phi - ai * (twon - u);
        dn = phi + ai * (twon + u);
        return;
    }
    double a[9], c[9];
    a[0] = 1.0;
    double b = std::sqrt(1.0 - m);
    c[0] = std::sqrt(m);
    double twon = 1.0;
    int i = 0;
    while (std::fabs(c[i] / a[i]) > MACHEP) {
        if (i > 7) break;
        const double ai = a[i];
        ++i;
        c[i] = (ai - b) / 2.0;
        const double t = std::sqrt(ai * b);
        a[i] = (ai + b) / 2.0;
        b = t;
        twon *= 2.0;
    }
    double phi = twon * a[i] * u;
    do {
        const double t = c[i] * std::sin(phi) / a[i];
        b = phi;
        phi = (std::asin(t) + phi) / 2.0;
    } while (--i);
    sn = std::sin(phi);
    const double t = std::cos(phi);
    cn = t;
    dn = t / std::cos(phi - b);
    ph = phi;
}

double bessel_i0(double x) {
    // power series, converges for all x (used for Kaiser windows, beta < ~50)
    const double y = x * x / 4.0;
    double sum = 1.0, term = 1.0;
    for (int k = 1; k < 500; ++k) {
        term *= y / (double(k) * double(k));
        sum += term;
        if (term < 1e-17 * sum) break;
    }
    return sum;
}

double sinc(double x) {
    if (x == 0.0) return 1.0;
    const double y = PI * x;
    return std::sin(y) / y;
}

double pow10m1(double x) { return std::expm1(std::log(10.0) * x); }

double fminbound(const std::function<double(double)> &func, double x1, double x2,
                 double xatol, int maxfun) {
    const double sqrt_eps = std::sqrt(2.2e-16);
    const double golden_mean = 0.5 * (3.0 - std::sqrt(5.0));
    double a = x1, b = x2;
    double fulc = a + golden_mean * (b - a);
    double nfc = fulc, xf = fulc;
    double rat = 0.0, e = 0.0;
    double x = xf;
    double fx = func(x);
    int num = 1;
    double ffulc = fx, fnfc = fx;
    double xm = 0.5 * (a + b);
    double tol1 = sqrt_eps * std::fabs(xf) + xatol / 3.0;
    double tol2 = 2.0 * tol1;
    auto sign = [](double v) { return v > 0 ? 1.0 : (v < 0 ? -1.0 : 0.0); };

    while (std::fabs(xf - xm) > (tol2 - 0.5 * (b - a))) {
        bool golden = true;
        if (std::fabs(e) > tol1) {
            golden = false;
            double r = (xf - nfc) * (fx - ffulc);
            double q = (xf - fulc) * (fx - fnfc);
            double p = (xf - fulc) * q - (xf - nfc) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = std::fabs(q);
            r = e;
            e = rat;
            if ((std::fabs(p) < std::fabs(0.5 * q * r)) && (p > q * (a - xf)) &&
                (p < q * (b - xf))) {
                rat = (p + 0.0) / q;
                x = xf + rat;
                if (((x - a) < tol2) || ((b - x) < tol2)) {
                    const double si = sign(xm - xf) + ((xm - xf) == 0 ? 1.0 : 0.0);
                    rat = tol1 * si;
                }
            } else {
                golden = true;
            }
        }
        if (golden) {
            e = (xf >= xm) ? a - xf : b - xf;
            rat = golden_mean * e;
        }
        const double si = sign(rat) + (rat == 0 ? 1.0 : 0.0);
        x = xf + si * std::max(std::fabs(rat), tol1);
        const double fu = func(x);
        ++num;
        if (fu <= fx) {
            if (x >= xf) a = xf; else b = xf;
            fulc = nfc; ffulc = fnfc;
            nfc = xf; fnfc = fx;
            xf = x; fx = fu;
        } else {
            if (x < xf) a = x; else b = x;
            if ((fu <= fnfc) || (nfc == xf)) {
                fulc = nfc; ffulc = fnfc;
                nfc = x; fnfc = fu;
            } else if ((fu <= ffulc) || (fulc == xf) || (fulc == nfc)) {
                fulc = x; ffulc = fu;
            }
        }
        xm = 0.5 * (a + b);
        tol1 = sqrt_eps * std::fabs(xf) + xatol / 3.0;
        tol2 = 2.0 * tol1;
        if (num >= maxfun) break;
    }
    return xf;
}

}  // namespace pyfda
