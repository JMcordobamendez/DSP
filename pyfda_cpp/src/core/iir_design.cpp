#include "iir_design.hpp"

#include "poly.hpp"
#include "special.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace pyfda {

namespace {
cplx prod(const CVec &v, double sign = 1.0) {
    cplx r = 1.0;
    for (const cplx &x : v) r *= sign * x;
    return r;
}
void check_order(int N) {
    if (N < 0) throw DesignError("Filter order must be a nonnegative integer");
}
}  // namespace

Zpk buttap(int N) {
    check_order(N);
    Zpk zpk;
    for (int m = -N + 1; m < N; m += 2) zpk.p.push_back(-std::exp(cplx(0, PI * m / (2.0 * N))));
    zpk.k = 1.0;
    return zpk;
}

Zpk cheb1ap(int N, double rp) {
    check_order(N);
    Zpk zpk;
    if (N == 0) {
        zpk.k = std::pow(10.0, -rp / 20.0);
        return zpk;
    }
    const double eps = std::sqrt(std::pow(10.0, 0.1 * rp) - 1.0);
    const double mu = 1.0 / N * std::asinh(1.0 / eps);
    for (int m = -N + 1; m < N; m += 2) {
        const double theta = PI * m / (2.0 * N);
        zpk.p.push_back(-std::sinh(cplx(mu, theta)));
    }
    zpk.k = prod(zpk.p, -1.0).real();
    if (N % 2 == 0) zpk.k /= std::sqrt(1 + eps * eps);
    return zpk;
}

Zpk cheb2ap(int N, double rs) {
    check_order(N);
    Zpk zpk;
    if (N == 0) return zpk;
    const double de = 1.0 / std::sqrt(std::pow(10.0, 0.1 * rs) - 1.0);
    const double mu = std::asinh(1.0 / de) / N;
    Vec m;
    if (N % 2) {
        for (int i = -N + 1; i < 0; i += 2) m.push_back(i);
        for (int i = 2; i < N; i += 2) m.push_back(i);
    } else {
        for (int i = -N + 1; i < N; i += 2) m.push_back(i);
    }
    for (double mi : m) zpk.z.push_back(cplx(0, 1) / std::sin(mi * PI / (2.0 * N)));
    for (int i = -N + 1; i < N; i += 2) {
        const double theta = PI * i / (2.0 * N);
        zpk.p.push_back(-1.0 / std::sinh(cplx(mu, theta)));
    }
    zpk.k = (prod(zpk.p, -1.0) / prod(zpk.z, -1.0)).real();
    return zpk;
}

namespace {
constexpr int ELLIPDEG_MMAX = 7;
constexpr int ARC_JAC_SN_MAXITER = 10;

double ellipdeg(int n, double m1) {
    const double K1 = ellipk(m1);
    const double K1p = ellipkm1(m1);
    const double q1 = std::exp(-PI * K1p / K1);
    const double q = std::pow(q1, 1.0 / n);
    double num = 0.0, den = 0.0;
    for (int m = 0; m <= ELLIPDEG_MMAX; ++m) num += std::pow(q, double(m * (m + 1)));
    for (int m = 1; m <= ELLIPDEG_MMAX + 1; ++m) den += std::pow(q, double(m * m));
    den = 1 + 2 * den;
    return 16 * q * std::pow(num / den, 4);
}

cplx arc_jac_sn(cplx w, double m) {
    auto complement = [](cplx kx) { return std::sqrt((1.0 - kx) * (1.0 + kx)); };
    const double k = std::sqrt(m);
    if (k > 1) return cplx(std::nan(""), 0);
    if (k == 1) return std::atanh(w);
    Vec ks{k};
    int niter = 0;
    while (ks.back() != 0) {
        const double k_ = ks.back();
        const double k_p = std::sqrt((1 - k_) * (1 + k_));
        ks.push_back((1 - k_p) / (1 + k_p));
        if (++niter > ARC_JAC_SN_MAXITER) throw DesignError("Landen transformation not converging");
    }
    double K = PI / 2;
    for (size_t i = 1; i < ks.size(); ++i) K *= 1 + ks[i];
    cplx wn = w;
    for (size_t i = 0; i + 1 < ks.size(); ++i) {
        const double kn = ks[i], knext = ks[i + 1];
        wn = 2.0 * wn / ((1 + knext) * (1.0 + complement(kn * wn)));
    }
    const cplx u = 2.0 / PI * std::asin(wn);
    return K * u;
}

double arc_jac_sc1(double w, double m) {
    const cplx zc = arc_jac_sn(cplx(0, w), m);
    if (std::fabs(zc.real()) > 1e-14) throw DesignError("arc_jac_sc1 failed");
    return zc.imag();
}
}  // namespace

Zpk ellipap(int N, double rp, double rs) {
    check_order(N);
    constexpr double EPSILON = 2e-16;
    Zpk zpk;
    if (N == 0) {
        zpk.k = std::pow(10.0, -rp / 20.0);
        return zpk;
    }
    if (N == 1) {
        const double p = -std::sqrt(1.0 / pow10m1(0.1 * rp));
        zpk.p = {p};
        zpk.k = -p;
        return zpk;
    }
    const double eps_sq = pow10m1(0.1 * rp);
    const double eps = std::sqrt(eps_sq);
    const double ck1_sq = eps_sq / pow10m1(0.1 * rs);
    if (ck1_sq == 0) throw DesignError("Cannot design a filter with given rp and rs specifications.");
    const double val0 = ellipk(ck1_sq);
    const double m = ellipdeg(N, ck1_sq);
    const double capk = ellipk(m);
    std::vector<int> j;
    for (int i = 1 - N % 2; i < N; i += 2) j.push_back(i);
    Vec s(j.size()), c(j.size()), d(j.size());
    for (size_t i = 0; i < j.size(); ++i) {
        double ph;
        ellipj(j[i] * capk / N, m, s[i], c[i], d[i], ph);
    }
    CVec z;
    for (double si : s)
        if (std::fabs(si) > EPSILON) z.push_back(cplx(0, 1.0 / (std::sqrt(m) * si)));
    const size_t nz = z.size();
    for (size_t i = 0; i < nz; ++i) z.push_back(std::conj(z[i]));

    const double r = arc_jac_sc1(1.0 / eps, ck1_sq);
    const double v0 = capk * r / (N * val0);
    double sv, cv, dv, ph;
    ellipj(v0, 1 - m, sv, cv, dv, ph);
    CVec p(j.size());
    for (size_t i = 0; i < j.size(); ++i)
        p[i] = -(c[i] * d[i] * sv * cv + cplx(0, 1) * s[i] * dv) / (1 - std::pow(d[i] * sv, 2.0));
    if (N % 2) {
        double ss = 0.0;
        for (const cplx &v : p) ss += std::norm(v);
        const size_t np0 = p.size();
        for (size_t i = 0; i < np0; ++i)
            if (std::fabs(p[i].imag()) > EPSILON * std::sqrt(ss)) p.push_back(std::conj(p[i]));
    } else {
        const size_t np0 = p.size();
        for (size_t i = 0; i < np0; ++i) p.push_back(std::conj(p[i]));
    }
    double k = (prod(p, -1.0) / prod(z, -1.0)).real();
    if (N % 2 == 0) k /= std::sqrt(1 + eps_sq);
    zpk.z = z;
    zpk.p = p;
    zpk.k = k;
    return zpk;
}

Zpk besselap(int N) {
    check_order(N);
    Zpk zpk;
    zpk.k = 1.0;
    if (N == 0) return zpk;
    // reverse Bessel polynomial theta_N(s) = sum a_k s^k, a_k = (2N-k)! / (2^(N-k) k! (N-k)!)
    Vec a(N + 1);
    for (int k = 0; k <= N; ++k) {
        double v = 1.0;  // (2N-k)! / (N-k)! / k! / 2^(N-k), computed in logs for large N
        const double lv = std::lgamma(2.0 * N - k + 1) - std::lgamma(N - k + 1.0) - std::lgamma(k + 1.0) -
                          (N - k) * std::log(2.0);
        v = std::exp(lv);
        a[k] = v;
    }
    Vec coeffs(a.rbegin(), a.rend());  // highest power first: a_N s^N + ... + a_0
    CVec p = roots(coeffs);
    // polish each root with Newton iterations on the full polynomial
    for (cplx &x : p) {
        for (int it = 0; it < 20; ++it) {
            cplx f = 0.0, fp = 0.0;
            for (double c : coeffs) {
                fp = fp * x + f;
                f = f * x + c;
            }
            if (fp == 0.0) break;
            const cplx dx = f / fp;
            x -= dx;
            if (std::abs(dx) <= 1e-16 * std::abs(x)) break;
        }
    }
    // symmetrize conjugate pairs
    std::sort(p.begin(), p.end(), [](const cplx &u, const cplx &v) { return u.imag() < v.imag(); });
    for (int i = 0; i < N / 2; ++i) {
        const cplx m = 0.5 * (p[i] + std::conj(p[N - 1 - i]));
        p[i] = m;
        p[N - 1 - i] = std::conj(m);
    }
    if (N % 2) p[N / 2] = p[N / 2].real();
    // norm = 'phase': scale so that the phase is -N*pi/4 at w = 1 (a_last = a_0 = a_N * prod(-p))
    const double a_last = a[0] / a[N];
    const double scale = std::pow(10.0, -std::log10(a_last) / N);
    for (cplx &x : p) x *= scale;
    zpk.p = p;
    return zpk;
}

namespace {
int relative_degree(const Zpk &zpk) {
    const int d = int(zpk.p.size()) - int(zpk.z.size());
    if (d < 0) throw DesignError("Improper transfer function. Must have at least as many poles as zeros.");
    return d;
}
}  // namespace

Zpk lp2lp_zpk(const Zpk &zpk, double wo) {
    const int degree = relative_degree(zpk);
    Zpk out;
    for (const cplx &z : zpk.z) out.z.push_back(wo * z);
    for (const cplx &p : zpk.p) out.p.push_back(wo * p);
    out.k = zpk.k * std::pow(wo, degree);
    return out;
}

Zpk lp2hp_zpk(const Zpk &zpk, double wo) {
    const int degree = relative_degree(zpk);
    Zpk out;
    for (const cplx &z : zpk.z) out.z.push_back(wo / z);
    for (const cplx &p : zpk.p) out.p.push_back(wo / p);
    out.z.insert(out.z.end(), degree, cplx(0.0));
    out.k = zpk.k * (prod(zpk.z, -1.0) / prod(zpk.p, -1.0)).real();
    return out;
}

Zpk lp2bp_zpk(const Zpk &zpk, double wo, double bw) {
    const int degree = relative_degree(zpk);
    Zpk out;
    CVec zl, pl;
    for (const cplx &z : zpk.z) zl.push_back(z * bw / 2.0);
    for (const cplx &p : zpk.p) pl.push_back(p * bw / 2.0);
    for (const cplx &z : zl) out.z.push_back(z + std::sqrt(z * z - wo * wo));
    for (const cplx &z : zl) out.z.push_back(z - std::sqrt(z * z - wo * wo));
    for (const cplx &p : pl) out.p.push_back(p + std::sqrt(p * p - wo * wo));
    for (const cplx &p : pl) out.p.push_back(p - std::sqrt(p * p - wo * wo));
    out.z.insert(out.z.end(), degree, cplx(0.0));
    out.k = zpk.k * std::pow(bw, degree);
    return out;
}

Zpk lp2bs_zpk(const Zpk &zpk, double wo, double bw) {
    const int degree = relative_degree(zpk);
    Zpk out;
    CVec zh, ph;
    for (const cplx &z : zpk.z) zh.push_back((bw / 2.0) / z);
    for (const cplx &p : zpk.p) ph.push_back((bw / 2.0) / p);
    for (const cplx &z : zh) out.z.push_back(z + std::sqrt(z * z - wo * wo));
    for (const cplx &z : zh) out.z.push_back(z - std::sqrt(z * z - wo * wo));
    for (const cplx &p : ph) out.p.push_back(p + std::sqrt(p * p - wo * wo));
    for (const cplx &p : ph) out.p.push_back(p - std::sqrt(p * p - wo * wo));
    out.z.insert(out.z.end(), degree, cplx(0, wo));
    out.z.insert(out.z.end(), degree, cplx(0, -wo));
    out.k = zpk.k * (prod(zpk.z, -1.0) / prod(zpk.p, -1.0)).real();
    return out;
}

Zpk bilinear_zpk(const Zpk &zpk, double fs) {
    const int degree = relative_degree(zpk);
    const double fs2 = 2.0 * fs;
    Zpk out;
    cplx num = 1.0, den = 1.0;
    for (const cplx &z : zpk.z) {
        out.z.push_back((fs2 + z) / (fs2 - z));
        num *= fs2 - z;
    }
    for (const cplx &p : zpk.p) {
        out.p.push_back((fs2 + p) / (fs2 - p));
        den *= fs2 - p;
    }
    out.z.insert(out.z.end(), degree, cplx(-1.0));
    out.k = zpk.k * (num / den).real();
    return out;
}

Zpk iirfilter(int N, const Vec &Wn, double rp, double rs, BType btype, IirType ftype) {
    if (Wn.empty()) throw DesignError("No critical frequency given.");
    for (double w : Wn)
        if (w <= 0 || w >= 1)
            throw DesignError("Digital filter critical frequencies must be 0 < Wn < f_S / 2.");
    if (Wn.size() > 1 && !(Wn[0] < Wn[1])) throw DesignError("Wn[0] must be less than Wn[1].");
    if (rp < 0) throw DesignError("Passband ripple (rp) must be positive.");
    if (rs < 0) throw DesignError("Stopband attenuation (rs) must be positive.");

    Zpk zpk;
    switch (ftype) {
    case IirType::Butter: zpk = buttap(N); break;
    case IirType::Cheby1: zpk = cheb1ap(N, rp); break;
    case IirType::Cheby2: zpk = cheb2ap(N, rs); break;
    case IirType::Ellip: zpk = ellipap(N, rp, rs); break;
    case IirType::Bessel: zpk = besselap(N); break;
    }
    const double fs = 2.0;
    Vec warped;
    for (double w : Wn) warped.push_back(2 * fs * std::tan(PI * w / fs));

    if (btype == BType::Lowpass || btype == BType::Highpass) {
        if (Wn.size() != 1)
            throw DesignError("Must specify a single critical frequency for a lowpass or highpass filter.");
        zpk = (btype == BType::Lowpass) ? lp2lp_zpk(zpk, warped[0]) : lp2hp_zpk(zpk, warped[0]);
    } else {
        if (Wn.size() != 2)
            throw DesignError("Must specify two critical frequencies for a bandpass or bandstop filter.");
        const double bw = warped[1] - warped[0];
        const double wo = std::sqrt(warped[0] * warped[1]);
        zpk = (btype == BType::Bandpass) ? lp2bp_zpk(zpk, wo, bw) : lp2bs_zpk(zpk, wo, bw);
    }
    return bilinear_zpk(zpk, fs);
}

// ---------------------------------------------------------------------------
// Minimum order estimation
// ---------------------------------------------------------------------------
namespace {

enum class OrdKind { Butter, Cheby, Ellip };

struct WpWs {
    Vec passb, stopb;
    int filter_type;  // 1 low, 2 high, 3 stop, 4 pass
};

WpWs validate_prewarp(const Vec &wp, const Vec &ws, double gpass, double gstop) {
    if (wp.empty() || wp.size() != ws.size() || wp.size() > 2)
        throw DesignError("wp and ws must both have one or two elements.");
    if (gpass <= 0) throw DesignError("gpass should be larger than 0.0");
    if (gstop <= 0) throw DesignError("gstop should be larger than 0.0");
    if (gpass > gstop) throw DesignError("gpass should be smaller than gstop");
    for (double w : wp)
        if (w <= 0 || w >= 1) throw DesignError("Pass band frequencies must be 0 < f < f_S / 2.");
    for (double w : ws)
        if (w <= 0 || w >= 1) throw DesignError("Stop band frequencies must be 0 < f < f_S / 2.");
    WpWs r;
    r.filter_type = 2 * (int(wp.size()) - 1) + 1;
    if (wp[0] >= ws[0]) r.filter_type += 1;
    for (double w : wp) r.passb.push_back(std::tan(PI * w / 2.0));
    for (double w : ws) r.stopb.push_back(std::tan(PI * w / 2.0));
    if (wp.size() == 2) {
        if (r.filter_type == 4 && !(ws[0] < wp[0] && wp[1] < ws[1]))
            throw DesignError("Band pass: F_SB < F_PB < F_PB2 < F_SB2 is required.");
        if (r.filter_type == 3 && !(wp[0] < ws[0] && ws[1] < wp[1]))
            throw DesignError("Band stop: F_PB < F_SB < F_SB2 < F_PB2 is required.");
    }
    return r;
}

double band_stop_obj(double wp, int ind, const Vec &passb, const Vec &stopb, double gpass,
                     double gstop, OrdKind kind) {
    Vec pc = passb;
    pc[ind] = wp;
    double nat = std::numeric_limits<double>::infinity();
    for (double s : stopb) nat = std::min(nat, std::fabs(s * (pc[0] - pc[1]) / (s * s - pc[0] * pc[1])));
    const double GSTOP = std::pow(10.0, 0.1 * std::fabs(gstop));
    const double GPASS = std::pow(10.0, 0.1 * std::fabs(gpass));
    switch (kind) {
    case OrdKind::Butter:
        return std::log10((GSTOP - 1.0) / (GPASS - 1.0)) / (2 * std::log10(nat));
    case OrdKind::Cheby:
        return std::acosh(std::sqrt((GSTOP - 1.0) / (GPASS - 1.0))) / std::acosh(nat);
    case OrdKind::Ellip: {
        const double arg1 = std::sqrt((GPASS - 1.0) / (GSTOP - 1.0));
        const double arg0 = 1.0 / nat;
        return ellipk(arg0 * arg0) * ellipk(1 - arg1 * arg1) /
               (ellipk(1 - arg0 * arg0) * ellipk(arg1 * arg1));
    }
    }
    return 0.0;
}

double find_nat_freq(const Vec &stopb, Vec &passb, double gpass, double gstop, int filter_type,
                     OrdKind kind) {
    Vec nat;
    switch (filter_type) {
    case 1: nat = {stopb[0] / passb[0]}; break;
    case 2: nat = {passb[0] / stopb[0]}; break;
    case 3: {
        const Vec pb = passb;
        const double wp0 = fminbound(
            [&](double w) { return band_stop_obj(w, 0, pb, stopb, gpass, gstop, kind); }, pb[0],
            stopb[0] - 1e-12);
        const double wp1 = fminbound(
            [&](double w) { return band_stop_obj(w, 1, pb, stopb, gpass, gstop, kind); },
            stopb[1] + 1e-12, pb[1]);
        passb = {wp0, wp1};
        for (double s : stopb) nat.push_back(s * (passb[0] - passb[1]) / (s * s - passb[0] * passb[1]));
        break;
    }
    case 4:
        for (double s : stopb) nat.push_back((s * s - passb[0] * passb[1]) / (s * (passb[0] - passb[1])));
        break;
    }
    double m = std::numeric_limits<double>::infinity();
    for (double v : nat) m = std::min(m, std::fabs(v));
    return m;
}

Vec postprocess_wn(const Vec &WN) {
    Vec wn;
    for (double w : WN) wn.push_back(std::atan(w) * 2.0 / PI);
    return wn;
}

}  // namespace

OrdResult buttord(const Vec &wp, const Vec &ws, double gpass, double gstop) {
    WpWs v = validate_prewarp(wp, ws, gpass, gstop);
    const double nat = find_nat_freq(v.stopb, v.passb, gpass, gstop, v.filter_type, OrdKind::Butter);
    const double GSTOP = std::pow(10.0, 0.1 * std::fabs(gstop));
    const double GPASS = std::pow(10.0, 0.1 * std::fabs(gpass));
    const int ord = int(std::ceil(std::log10((GSTOP - 1.0) / (GPASS - 1.0)) / (2 * std::log10(nat))));
    const double W0 = ord == 0 ? 1.0 : std::pow(GPASS - 1.0, -1.0 / (2.0 * ord));
    const Vec &pb = v.passb;
    Vec WN;
    switch (v.filter_type) {
    case 1: WN = {W0 * pb[0]}; break;
    case 2: WN = {pb[0] / W0}; break;
    case 3: {
        const double discr = std::sqrt(std::pow(pb[1] - pb[0], 2) + 4 * W0 * W0 * pb[0] * pb[1]);
        WN = {std::fabs(((pb[1] - pb[0]) + discr) / (2 * W0)), std::fabs(((pb[1] - pb[0]) - discr) / (2 * W0))};
        std::sort(WN.begin(), WN.end());
        break;
    }
    case 4: {
        for (double w0 : {-W0, W0})
            WN.push_back(std::fabs(-w0 * (pb[1] - pb[0]) / 2.0 +
                                   std::sqrt(w0 * w0 / 4.0 * std::pow(pb[1] - pb[0], 2) + pb[0] * pb[1])));
        std::sort(WN.begin(), WN.end());
        break;
    }
    }
    return {ord, postprocess_wn(WN)};
}

OrdResult cheb1ord(const Vec &wp, const Vec &ws, double gpass, double gstop) {
    WpWs v = validate_prewarp(wp, ws, gpass, gstop);
    const double nat = find_nat_freq(v.stopb, v.passb, gpass, gstop, v.filter_type, OrdKind::Cheby);
    const double GSTOP = std::pow(10.0, 0.1 * std::fabs(gstop));
    const double GPASS = std::pow(10.0, 0.1 * std::fabs(gpass));
    const double v_pass_stop = std::acosh(std::sqrt((GSTOP - 1.0) / (GPASS - 1.0)));
    const int ord = int(std::ceil(v_pass_stop / std::acosh(nat)));
    return {ord, postprocess_wn(v.passb)};
}

OrdResult cheb2ord(const Vec &wp, const Vec &ws, double gpass, double gstop) {
    WpWs v = validate_prewarp(wp, ws, gpass, gstop);
    const double nat0 = find_nat_freq(v.stopb, v.passb, gpass, gstop, v.filter_type, OrdKind::Cheby);
    const double GSTOP = std::pow(10.0, 0.1 * std::fabs(gstop));
    const double GPASS = std::pow(10.0, 0.1 * std::fabs(gpass));
    const double v_pass_stop = std::acosh(std::sqrt((GSTOP - 1.0) / (GPASS - 1.0)));
    const int ord = int(std::ceil(v_pass_stop / std::acosh(nat0)));
    const double new_freq = 1.0 / std::cosh(1.0 / ord * v_pass_stop);
    const Vec &pb = v.passb;
    Vec nat;
    switch (v.filter_type) {
    case 1: nat = {pb[0] / new_freq}; break;
    case 2: nat = {pb[0] * new_freq}; break;
    case 3: {
        const double n0 = new_freq / 2.0 * (pb[0] - pb[1]) +
                          std::sqrt(new_freq * new_freq * std::pow(pb[1] - pb[0], 2) / 4.0 + pb[1] * pb[0]);
        nat = {n0, pb[1] * pb[0] / n0};
        break;
    }
    case 4: {
        const double n0 = 1.0 / (2.0 * new_freq) * (pb[0] - pb[1]) +
                          std::sqrt(std::pow(pb[1] - pb[0], 2) / (4.0 * new_freq * new_freq) + pb[1] * pb[0]);
        nat = {n0, pb[0] * pb[1] / n0};
        break;
    }
    }
    return {ord, postprocess_wn(nat)};
}

OrdResult ellipord(const Vec &wp, const Vec &ws, double gpass, double gstop) {
    WpWs v = validate_prewarp(wp, ws, gpass, gstop);
    const double nat = find_nat_freq(v.stopb, v.passb, gpass, gstop, v.filter_type, OrdKind::Ellip);
    const double arg1_sq = pow10m1(0.1 * gpass) / pow10m1(0.1 * gstop);
    const double arg0 = 1.0 / nat;
    const double d00 = ellipk(arg0 * arg0), d01 = ellipkm1(arg0 * arg0);
    const double d10 = ellipk(arg1_sq), d11 = ellipkm1(arg1_sq);
    const int ord = int(std::ceil(d00 * d11 / (d01 * d10)));
    return {ord, postprocess_wn(v.passb)};
}

}  // namespace pyfda
