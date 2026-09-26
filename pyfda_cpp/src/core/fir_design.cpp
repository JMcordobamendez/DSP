#include "fir_design.hpp"

#include "special.hpp"

#include <algorithm>
#include <cmath>

namespace pyfda {

const std::vector<WindowInfo> &window_list() {
    static const std::vector<WindowInfo> list = {
        {WindowType::Rectangular, "Rectangular", nullptr, 0},
        {WindowType::Bartlett, "Bartlett", nullptr, 0},
        {WindowType::Hann, "Hann", nullptr, 0},
        {WindowType::Hamming, "Hamming", nullptr, 0},
        {WindowType::Blackman, "Blackman", nullptr, 0},
        {WindowType::BlackmanHarris, "Blackman-Harris", nullptr, 0},
        {WindowType::Nuttall, "Nuttall", nullptr, 0},
        {WindowType::Flattop, "Flattop", nullptr, 0},
        {WindowType::Kaiser, "Kaiser", "beta", 10.0},
        {WindowType::Gaussian, "Gaussian", "std", 0.25},  // relative to (M-1)/2
        {WindowType::Tukey, "Tukey", "alpha", 0.5},
    };
    return list;
}

namespace {
Vec general_cosine(int M, const Vec &a) {
    Vec w(M, 0.0);
    for (int n = 0; n < M; ++n) {
        const double fac = M == 1 ? 0.0 : -PI + 2 * PI * n / (M - 1);
        for (size_t k = 0; k < a.size(); ++k) w[n] += a[k] * std::cos(double(k) * fac);
    }
    return w;
}
}  // namespace

Vec get_window(WindowType type, int M, double par) {
    if (M < 1) return {};
    if (M == 1) return {1.0};
    Vec w(M, 1.0);
    switch (type) {
    case WindowType::Rectangular: break;
    case WindowType::Bartlett:
        for (int n = 0; n < M; ++n)
            w[n] = n <= (M - 1) / 2.0 ? 2.0 * n / (M - 1) : 2.0 - 2.0 * n / (M - 1);
        break;
    case WindowType::Hann: w = general_cosine(M, {0.5, 0.5}); break;
    case WindowType::Hamming: w = general_cosine(M, {0.54, 0.46}); break;
    case WindowType::Blackman: w = general_cosine(M, {0.42, 0.50, 0.08}); break;
    case WindowType::BlackmanHarris: w = general_cosine(M, {0.35875, 0.48829, 0.14128, 0.01168}); break;
    case WindowType::Nuttall: w = general_cosine(M, {0.3635819, 0.4891775, 0.1365995, 0.0106411}); break;
    case WindowType::Flattop:
        w = general_cosine(M, {0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368});
        break;
    case WindowType::Kaiser: {
        const double alpha = (M - 1) / 2.0;
        const double i0b = bessel_i0(par);
        for (int n = 0; n < M; ++n) {
            const double r = (n - alpha) / alpha;
            w[n] = bessel_i0(par * std::sqrt(std::max(0.0, 1.0 - r * r))) / i0b;
        }
        break;
    }
    case WindowType::Gaussian: {
        // std given relative to (M-1)/2 so that it doesn't depend on the length
        const double std_n = par * (M - 1) / 2.0;
        for (int n = 0; n < M; ++n) {
            const double x = n - (M - 1) / 2.0;
            w[n] = std::exp(-0.5 * x * x / (std_n * std_n));
        }
        break;
    }
    case WindowType::Tukey: {
        const double alpha = par;
        if (alpha <= 0) break;
        if (alpha >= 1.0) return get_window(WindowType::Hann, M);
        const int width = int(std::floor(alpha * (M - 1) / 2.0));
        for (int n = 0; n < M; ++n) {
            if (n <= width)
                w[n] = 0.5 * (1 + std::cos(PI * (-1 + 2.0 * n / alpha / (M - 1))));
            else if (n >= M - width - 1)
                w[n] = 0.5 * (1 + std::cos(PI * (-2.0 / alpha + 1 + 2.0 * n / alpha / (M - 1))));
        }
        break;
    }
    }
    return w;
}

Vec firwin(int numtaps, const Vec &cutoff_in, const Vec &window, bool pass_zero, bool scale) {
    if (cutoff_in.empty()) throw DesignError("At least one cutoff frequency must be given.");
    for (double c : cutoff_in)
        if (c <= 0 || c >= 1)
            throw DesignError("Invalid cutoff frequency: frequencies must be greater than 0 and less than f_S / 2.");
    for (size_t i = 1; i < cutoff_in.size(); ++i)
        if (cutoff_in[i] <= cutoff_in[i - 1])
            throw DesignError("Invalid cutoff frequencies: the frequencies must be strictly increasing.");
    const bool pass_nyquist = bool(cutoff_in.size() & 1) ^ pass_zero;
    if (pass_nyquist && numtaps % 2 == 0)
        throw DesignError("A filter with an even number of coefficients (odd order) must have zero "
                          "response at f_S / 2, use an even order.");
    Vec cutoff;
    if (pass_zero) cutoff.push_back(0.0);
    cutoff.insert(cutoff.end(), cutoff_in.begin(), cutoff_in.end());
    if (pass_nyquist) cutoff.push_back(1.0);

    const double alpha = 0.5 * (numtaps - 1);
    Vec h(numtaps, 0.0);
    for (size_t b = 0; b + 1 < cutoff.size(); b += 2) {
        const double left = cutoff[b], right = cutoff[b + 1];
        for (int n = 0; n < numtaps; ++n) {
            const double m = n - alpha;
            h[n] += right * sinc(right * m) - left * sinc(left * m);
        }
    }
    if (int(window.size()) != numtaps) throw DesignError("Window length doesn't match the number of taps.");
    for (int n = 0; n < numtaps; ++n) h[n] *= window[n];

    if (scale) {
        const double left = cutoff[0], right = cutoff[1];
        double sf;
        if (left == 0) sf = 0.0;
        else if (right == 1) sf = 1.0;
        else sf = 0.5 * (left + right);
        double s = 0.0;
        for (int n = 0; n < numtaps; ++n) s += h[n] * std::cos(PI * (n - alpha) * sf);
        for (double &v : h) v /= s;
    }
    return h;
}

double kaiser_beta(double a) {
    if (a > 50) return 0.1102 * (a - 8.7);
    if (a > 21) return 0.5842 * std::pow(a - 21, 0.4) + 0.07886 * (a - 21);
    return 0.0;
}

void kaiserord(double ripple_db, double width, int &numtaps, double &beta) {
    const double A = std::fabs(ripple_db);
    if (A < 8) throw DesignError("Requested maximum ripple attenuation is too small for the Kaiser formula.");
    beta = kaiser_beta(A);
    numtaps = int(std::ceil((A - 7.95) / 2.285 / (PI * width) + 1));
}

namespace {
int remlplen_herrmann(double fp, double fs, double dp, double ds) {
    const double df = fs - fp;
    const double a[] = {5.309e-3, 7.114e-2, -4.761e-1, -2.66e-3, -5.941e-1, -4.278e-1};
    const double b[] = {11.01217, 0.51244};
    const double ldp = std::log10(dp), lds = std::log10(ds);
    const double dinf = lds * (a[0] * ldp * ldp + a[1] * ldp + a[2]) + a[3] * ldp * ldp + a[4] * ldp + a[5];
    const double f = b[0] + b[1] * (ldp - lds);
    return int(dinf / df - f * df + 1);
}

int remlplen_kaiser(double fp, double fs, double dp, double ds) {
    const double df = fs - fp;
    return int((-20 * std::log10(std::sqrt(dp * ds)) - 13.0) / (14.6 * df) + 1.0);
}

int remlplen_ichige(double fp, double fs, double dp, double ds) {
    auto func_v = [&](double df, double dp_) {
        return 2.325 * std::pow(-std::log10(dp_), -0.445) * std::pow(df, -1.39);
    };
    auto func_g = [&](double df, double fp_) {
        return (2.0 / PI) * std::atan(func_v(df, dp) * (1.0 / fp_ - 1.0 / (0.5 - df)));
    };
    auto func_h = [&](double df, double fp_, double c) {
        return (2.0 / PI) * std::atan((c / df) * (1.0 / fp_ - 1.0 / (0.5 - df)));
    };
    const double df = fs - fp;
    const double nc = std::ceil(1.0 + (1.101 / df) * std::pow(-std::log10(2.0 * dp), 1.1));
    const double nm = (0.52 / df) * std::log10(dp / ds) * std::pow(-std::log10(dp), 0.17);
    const double n3 = std::ceil(nc * (func_g(df, fp) + func_g(df, 0.5 - df - fp) + 1.0) / 3.0);
    const double dn = std::ceil(nm * (func_h(df, fp, 1.1) - (func_h(df, 0.5 - df - fp, 0.29) - 1.0) / 2.0));
    return int(n3 + dn);
}
}  // namespace

RemezOrd remezord(const Vec &freqs_in, const Vec &amps, const Vec &rips_in, double fs, RemezAlg alg) {
    Vec freqs(freqs_in), rips(rips_in);
    for (size_t i = 0; i < rips.size(); ++i) rips[i] /= (amps[i] + (amps[i] == 0.0 ? 1.0 : 0.0));
    for (double &f : freqs) f /= fs;
    for (double f : freqs) {
        if (f > 0.5) throw DesignError("Frequency band edges must not exceed the Nyquist frequency.");
        if (f < 0.0) throw DesignError("Frequency band edges must be nonnegative.");
    }
    for (double r : rips)
        if (r <= 0.0) throw DesignError("Ripples must be nonnegative and non-zero.");
    if (amps.size() != rips.size()) throw DesignError("Number of amplitudes must equal number of ripples.");
    if (freqs.size() != 2 * (amps.size() - 1))
        throw DesignError("Number of band edges must equal 2*(number of amplitudes-1)");
    auto remlplen = [&](double fp, double fst, double dp, double ds) {
        switch (alg) {
        case RemezAlg::Herrmann: return remlplen_herrmann(fp, fst, dp, ds);
        case RemezAlg::Kaiser: return remlplen_kaiser(fp, fst, dp, ds);
        default: return remlplen_ichige(fp, fst, dp, ds);
        }
    };
    int fil_len = 0;
    for (size_t i = 0; i + 1 < amps.size(); ++i) {
        const double f1 = freqs[2 * i], f2 = freqs[2 * i + 1];
        fil_len = std::max({fil_len, remlplen(f1, f2, rips[i], rips[i + 1]),
                            remlplen(0.5 - f2, 0.5 - f1, rips[i + 1], rips[i])});
    }
    RemezOrd r;
    r.numtaps = fil_len;
    r.bands.push_back(0.0);
    r.bands.insert(r.bands.end(), freqs.begin(), freqs.end());
    r.bands.push_back(0.5);
    r.desired = amps;
    const double max_rip = *std::max_element(rips.begin(), rips.end());
    for (double rp : rips) r.weight.push_back(max_rip / rp);
    return r;
}

}  // namespace pyfda
