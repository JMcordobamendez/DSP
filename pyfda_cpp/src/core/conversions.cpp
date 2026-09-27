#include "conversions.hpp"

#include "poly.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace pyfda {

Ba zpk2tf(const Zpk &zpk) {
    Ba ba;
    ba.b = poly_real(zpk.z);
    for (double &v : ba.b) v *= zpk.k;
    ba.a = poly_real(zpk.p);
    return ba;
}

namespace {

// Port of scipy.signal._filter_design._cplxreal: split into complex values with
// positive imaginary part (one per conjugate pair) and real values
void cplxreal(CVec z, CVec &zc, Vec &zr) {
    zc.clear();
    zr.clear();
    if (z.empty()) return;
    const double tol = 100 * std::numeric_limits<double>::epsilon();
    std::stable_sort(z.begin(), z.end(), [](const cplx &a, const cplx &b) {
        if (a.real() != b.real()) return a.real() < b.real();
        return std::fabs(a.imag()) < std::fabs(b.imag());
    });
    CVec zp, zn;
    for (const cplx &v : z) {
        if (std::fabs(v.imag()) <= tol * std::abs(v)) zr.push_back(v.real());
        else if (v.imag() > 0) zp.push_back(v);
        else zn.push_back(v);
    }
    if (zp.size() != zn.size())
        throw DesignError("Array contains complex value with no matching conjugate.");
    // sort runs with the same real part by |imag|
    size_t i = 0;
    while (i < zp.size()) {
        size_t j = i + 1;
        while (j < zp.size() && zp[j].real() - zp[j - 1].real() <= tol * std::abs(zp[j - 1])) ++j;
        if (j - i > 1) {
            auto by_imag = [](const cplx &a, const cplx &b) {
                return std::fabs(a.imag()) < std::fabs(b.imag());
            };
            std::stable_sort(zp.begin() + i, zp.begin() + j, by_imag);
            std::stable_sort(zn.begin() + i, zn.begin() + j, by_imag);
        }
        i = j;
    }
    for (size_t k = 0; k < zp.size(); ++k) {
        if (std::abs(zp[k] - std::conj(zn[k])) > tol * std::abs(zn[k]))
            throw DesignError("Array contains complex value with no matching conjugate.");
        zc.push_back((zp[k] + std::conj(zn[k])) / 2.0);
    }
}

bool is_real(const cplx &v) { return v.imag() == 0.0; }

std::array<double, 6> single_zpksos(const CVec &z, const CVec &p) {
    std::array<double, 6> s{0, 0, 0, 0, 0, 0};
    const Vec b = poly_real(z);
    const Vec a = poly_real(p);
    for (size_t i = 0; i < b.size(); ++i) s[3 - b.size() + i] = b[i];
    for (size_t i = 0; i < a.size(); ++i) s[6 - a.size() + i] = a[i];
    return s;
}

// which: 0 = real, 1 = complex, 2 = any
size_t nearest_idx(const CVec &fro, const cplx &to, int which) {
    std::vector<size_t> order(fro.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return std::abs(fro[a] - to) < std::abs(fro[b] - to);
    });
    for (size_t idx : order) {
        if (which == 2) return idx;
        if ((which == 0) == is_real(fro[idx])) return idx;
    }
    throw DesignError("zpk2sos: no matching zero found.");
}

size_t idx_worst(const CVec &p) {
    size_t best = 0;
    double best_val = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < p.size(); ++i) {
        const double v = std::fabs(1.0 - std::abs(p[i]));
        if (v < best_val) {
            best_val = v;
            best = i;
        }
    }
    return best;
}

size_t count_real(const CVec &v) {
    return size_t(std::count_if(v.begin(), v.end(), is_real));
}

cplx take(CVec &v, size_t idx) {
    const cplx x = v[idx];
    v.erase(v.begin() + long(idx));
    return x;
}

}  // namespace

Sos zpk2sos(const Zpk &zpk) {
    CVec z = zpk.z, p = zpk.p;
    if (z.empty() && p.empty()) return {{zpk.k, 0., 0., 1., 0., 0.}};
    if (z.size() < p.size()) z.resize(p.size(), 0.0);
    if (p.size() < z.size()) p.resize(z.size(), 0.0);
    const size_t n_sections = (std::max(p.size(), z.size()) + 1) / 2;
    if (p.size() % 2 == 1) {
        p.push_back(0.0);
        z.push_back(0.0);
    }
    {
        CVec zc;
        Vec zr;
        cplxreal(z, zc, zr);
        z = zc;
        for (double v : zr) z.push_back(v);
        cplxreal(p, zc, zr);
        p = zc;
        for (double v : zr) p.push_back(v);
    }
    Sos sos(n_sections);
    for (long si = long(n_sections) - 1; si >= 0; --si) {
        const cplx p1 = take(p, idx_worst(p));
        if (is_real(p1) && count_real(p) == 0) {
            const cplx z1 = take(z, nearest_idx(z, p1, 0));
            sos[si] = single_zpksos({z1, 0.0}, {p1, 0.0});
        } else if (p.size() + 1 == z.size() && !is_real(p1) && count_real(p) == 1 &&
                   count_real(z) == 1) {
            const cplx z1 = take(z, nearest_idx(z, p1, 1));
            sos[si] = single_zpksos({z1, std::conj(z1)}, {p1, std::conj(p1)});
        } else {
            cplx p2;
            if (is_real(p1)) {
                // worst remaining real pole
                size_t best = 0;
                double best_val = std::numeric_limits<double>::infinity();
                for (size_t i = 0; i < p.size(); ++i) {
                    if (!is_real(p[i])) continue;
                    const double v = std::fabs(1.0 - std::abs(p[i]));
                    if (v < best_val) {
                        best_val = v;
                        best = i;
                    }
                }
                p2 = take(p, best);
            } else {
                p2 = std::conj(p1);
            }
            if (!z.empty()) {
                const cplx z1 = take(z, nearest_idx(z, p1, 2));
                if (!is_real(z1)) {
                    sos[si] = single_zpksos({z1, std::conj(z1)}, {p1, p2});
                } else if (!z.empty()) {
                    const cplx z2 = take(z, nearest_idx(z, p1, 0));
                    sos[si] = single_zpksos({z1, z2}, {p1, p2});
                } else {
                    sos[si] = single_zpksos({z1}, {p1, p2});
                }
            } else {
                sos[si] = single_zpksos({}, {p1, p2});
            }
        }
    }
    for (int i = 0; i < 3; ++i) sos[0][i] *= zpk.k;
    return sos;
}

Ba sos2tf(const Sos &sos) {
    Ba ba{{1.0}, {1.0}};
    for (const auto &s : sos) {
        ba.b = polymul(ba.b, {s[0], s[1], s[2]});
        ba.a = polymul(ba.a, {s[3], s[4], s[5]});
    }
    return ba;
}

Zpk tf2zpk(const Ba &ba_in) {
    Vec b = ba_in.b, a = ba_in.a;
    while (a.size() > 1 && a[0] == 0.0) a.erase(a.begin());
    while (b.size() > 1 && b[0] == 0.0) b.erase(b.begin());
    if (a.empty() || a[0] == 0.0) throw DesignError("Denominator must not be zero.");
    Zpk zpk;
    zpk.k = b.empty() ? 0.0 : b[0] / a[0];
    // pad to equal length like scipy (zeros / poles at the origin)
    if (b.size() < a.size()) b.insert(b.end(), a.size() - b.size(), 0.0);
    if (a.size() < b.size()) a.insert(a.end(), b.size() - a.size(), 0.0);
    zpk.z = roots(b);
    zpk.p = roots(a);
    return zpk;
}

Zpk sos2zpk(const Sos &sos) {
    Zpk zpk;
    zpk.k = 1.0;
    for (const auto &s : sos) {
        const Zpk sec = tf2zpk({{s[0], s[1], s[2]}, {s[3], s[4], s[5]}});
        zpk.z.insert(zpk.z.end(), sec.z.begin(), sec.z.end());
        zpk.p.insert(zpk.p.end(), sec.p.begin(), sec.p.end());
        zpk.k *= sec.k;
    }
    return zpk;
}

namespace {
// sum c[n] z^-n
cplx eval_neg(const Vec &c, cplx zinv) {
    cplx y = 0.0;
    for (size_t i = c.size(); i-- > 0;) y = y * zinv + c[i];
    return y;
}
}  // namespace

CVec freqz(const Ba &ba, const Vec &w) {
    CVec H(w.size());
    for (size_t i = 0; i < w.size(); ++i) {
        const cplx zinv = std::polar(1.0, -w[i]);
        H[i] = eval_neg(ba.b, zinv) / eval_neg(ba.a, zinv);
    }
    return H;
}

CVec freqz(const Sos &sos, const Vec &w) {
    CVec H(w.size(), 1.0);
    for (size_t i = 0; i < w.size(); ++i) {
        const cplx zinv = std::polar(1.0, -w[i]);
        for (const auto &s : sos) {
            H[i] *= (s[0] + zinv * (s[1] + zinv * s[2])) / (s[3] + zinv * (s[4] + zinv * s[5]));
        }
    }
    return H;
}

Vec group_delay(const Ba &ba, const Vec &w) {
    // scipy.signal.group_delay: c = b * conj(a[::-1]), gd = Re(sum n c_n z^-n / sum c_n z^-n) - (len(a) - 1)
    Vec ar(ba.a.rbegin(), ba.a.rend());
    const Vec c = polymul(ba.b, ar);
    Vec cr(c.size());
    for (size_t n = 0; n < c.size(); ++n) cr[n] = c[n] * double(n);
    Vec gd(w.size());
    for (size_t i = 0; i < w.size(); ++i) {
        const cplx zinv = std::polar(1.0, -w[i]);
        const cplx num = eval_neg(cr, zinv);
        const cplx den = eval_neg(c, zinv);
        if (std::abs(den) < 10 * std::numeric_limits<double>::epsilon())
            gd[i] = std::numeric_limits<double>::quiet_NaN();
        else
            gd[i] = (num / den).real() - double(ba.a.size() - 1);
    }
    return gd;
}

Vec group_delay(const Sos &sos, const Vec &w) {
    Vec gd(w.size(), 0.0);
    for (const auto &s : sos) {
        const Vec g = group_delay(Ba{{s[0], s[1], s[2]}, {s[3], s[4], s[5]}}, w);
        for (size_t i = 0; i < w.size(); ++i) gd[i] += g[i];
    }
    return gd;
}

int impz_len(const Zpk &zpk, int n_min, int n_max) {
    double r_max = 0.0;
    for (const cplx &p : zpk.p) r_max = std::max(r_max, std::abs(p));
    if (r_max <= 0.0) return std::max(n_min, int(zpk.z.size()) + 1);
    if (r_max >= 1.0) return n_max;  // unstable or marginally stable
    // decay of the slowest pole to -80 dB, plus the settling of multiple poles
    const int n = int(std::ceil(std::log(1e-4) / std::log(r_max))) + int(zpk.p.size());
    return std::clamp(n, n_min, n_max);
}

Vec unwrap(const Vec &phi) {
    Vec out(phi);
    double offset = 0.0;
    for (size_t i = 1; i < phi.size(); ++i) {
        double d = phi[i] - phi[i - 1];
        if (std::isnan(d)) continue;
        double dd = std::fmod(d + PI, 2 * PI);
        if (dd < 0) dd += 2 * PI;
        dd -= PI;
        if (dd == -PI && d > 0) dd = PI;
        offset += dd - d;
        out[i] = phi[i] + offset;
    }
    return out;
}

}  // namespace pyfda
