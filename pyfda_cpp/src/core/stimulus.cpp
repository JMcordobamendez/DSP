#include "stimulus.hpp"

#include "filtering.hpp"
#include "special.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <sstream>

namespace pyfda {

const std::vector<StimInfo> &stim_list() {
    static const std::vector<StimInfo> l = {
        {Stim::None, "none", "None", "dc"},
        {Stim::Dirac, "dirac", "Dirac", "a1 t1 dc"},
        {Stim::Sinc, "sinc", "Sinc", "a1 a2 t1 t2 f1 f2 dc"},
        {Stim::Gauss, "gauss", "Gauss", "a1 a2 t1 t2 f1 f2 bw1 bw2 dc"},
        {Stim::Rect, "rect", "Rect", "a1 t1 tw dc"},
        {Stim::Step, "step", "Step", "a1 t1"},
        {Stim::Sine, "sine", "Sine", "a1 a2 phi1 phi2 f1 f2 dc"},
        {Stim::Cos, "cos", "Cos", "a1 a2 phi1 phi2 f1 f2 dc"},
        {Stim::Diric, "diric", "Diric", "a1 t1 n1 f1 dc"},
        {Stim::Chirp, "chirp", "Chirp", "a1 phi1 f1 f2 t2 dc"},
        {Stim::Triang, "triang", "Triangle", "a1 phi1 f1 bl dc"},
        {Stim::Saw, "saw", "Sawtooth", "a1 phi1 f1 bl dc"},
        {Stim::RectPer, "rect_per", "Rect (periodic)", "a1 phi1 f1 bl duty dc"},
        {Stim::Comb, "comb", "Comb", "a1 phi1 f1 dc"},
        {Stim::AM, "am", "AM", "a1 a2 phi1 phi2 f1 f2 dc"},
        {Stim::PMFM, "pmfm", "PM / FM", "a1 a2 phi1 phi2 f1 f2 dc"},
        {Stim::PWM, "pwm", "PWM", "a1 a2 phi1 phi2 f1 f2 bl dc"},
    };
    return l;
}

const StimInfo &stim_info(Stim s) {
    for (const auto &i : stim_list())
        if (i.stim == s) return i;
    return stim_list().front();
}

bool stim_uses(Stim s, const std::string &param) {
    std::istringstream ss(stim_info(s).params);
    std::string p;
    while (ss >> p)
        if (p == param) return true;
    return false;
}

bool is_impulse(Stim s) { return s == Stim::Dirac || s == Stim::Sinc || s == Stim::Gauss || s == Stim::Rect; }

// ---------------------------------------------------------------------------
// waveforms

namespace {
double pymod(double a, double b) {  // numpy.mod: result has the sign of b
    double r = std::fmod(a, b);
    if (r != 0 && ((r < 0) != (b < 0))) r += b;
    return r;
}
}  // namespace

double sawtooth(double t, double w) {
    if (w > 1 || w < 0) return NAN;
    const double tmod = pymod(t, 2 * PI);
    if (tmod < w * 2 * PI) return tmod / (PI * w) - 1;
    return (PI * (w + 1) - tmod) / (PI * (1 - w));
}

double square(double t, double duty) {
    if (duty > 1 || duty < 0) return NAN;
    return pymod(t, 2 * PI) < duty * 2 * PI ? 1.0 : -1.0;
}

double diric(double x, int n) {
    if (n <= 0) return NAN;
    x /= 2;
    const double denom = std::sin(x);
    if (std::fabs(denom) < 1e-7) return std::pow(-1.0, std::nearbyint(x / PI) * (n - 1));
    return std::sin(n * x) / (n * denom);
}

double gausspulse(double t, double fc, double bw, double bwr) {
    if (fc < 0) throw DesignError("Center frequency of the Gaussian pulse must be >= 0.");
    if (bw <= 0) throw DesignError("Fractional bandwidth of the Gaussian pulse must be > 0.");
    if (bwr >= 0) throw DesignError("Reference level for the bandwidth must be < 0 dB.");
    const double ref = std::pow(10.0, bwr / 20.0);
    const double a = -(PI * fc * bw) * (PI * fc * bw) / (4.0 * std::log(ref));
    return std::exp(-a * t * t) * std::cos(2 * PI * fc * t);
}

double chirp(double t, double f0, double t1, double f1, ChirpType method, double phi) {
    double phase = 0;
    switch (method) {
    case ChirpType::Linear: {
        const double beta = (f1 - f0) / t1;
        phase = 2 * PI * (f0 * t + 0.5 * beta * t * t);
        break;
    }
    case ChirpType::Quadratic: {  // vertex_zero = True
        const double beta = (f1 - f0) / (t1 * t1);
        phase = 2 * PI * (f0 * t + beta * t * t * t / 3);
        break;
    }
    case ChirpType::Logarithmic:
        if (f0 * f1 <= 0) throw DesignError("For a logarithmic chirp, f1 and f2 must be != 0 and have the same sign.");
        if (f0 == f1) phase = 2 * PI * f0 * t;
        else {
            const double beta = t1 / std::log(f1 / f0);
            phase = 2 * PI * beta * f0 * (std::pow(f1 / f0, t / t1) - 1.0);
        }
        break;
    case ChirpType::Hyperbolic:
        if (f0 == 0 || f1 == 0) throw DesignError("For a hyperbolic chirp, f1 and f2 must be != 0.");
        if (f0 == f1) phase = 2 * PI * f0 * t;
        else {
            const double sing = -f1 * t1 / (f0 - f1);
            phase = 2 * PI * (-sing * f0) * std::log(std::fabs(1 - t / sing));
        }
        break;
    }
    return std::cos(phase + phi * PI / 180.0);
}

namespace {
// number of harmonics below the Nyquist frequency, int(fs * pi) with fs = 1 / (t[1] - t[0])
int harmonics(const Vec &t) {
    if (t.size() < 2) throw DesignError("Bandlimited signals need at least 2 samples.");
    const double h = (1.0 / (t[1] - t[0])) * PI;
    if (!std::isfinite(h) || h * double(t.size()) > 5e8)
        throw DesignError("The frequency is too low for a bandlimited signal of this length, "
                          "switch off 'BL' or increase f1.");
    return int(h);
}
}  // namespace

Vec sawtooth_bl(const Vec &t) {
    Vec y(t.size(), 0.0);
    const int nh = harmonics(t);
    for (int h = 1; h <= nh; ++h)
        for (size_t i = 0; i < t.size(); ++i) y[i] += 2 / PI * -std::sin(h * t[i]) / h;
    return y;
}

Vec triang_bl(const Vec &t) {
    Vec y(t.size(), 0.0);
    const int nh = harmonics(t);
    for (int h = 1; h <= nh; h += 2)
        for (size_t i = 0; i < t.size(); ++i) y[i] += 8 / (PI * PI) * -std::cos(h * t[i]) / (double(h) * h);
    return y;
}

Vec rect_bl(const Vec &t, const Vec &duty) {
    Vec t2(t.size());
    for (size_t i = 0; i < t.size(); ++i) t2[i] = t[i] - duty[i] * 2 * PI;
    const Vec s1 = sawtooth_bl(t2), s2 = sawtooth_bl(t);
    Vec y(t.size());
    for (size_t i = 0; i < t.size(); ++i) y[i] = s1[i] - s2[i] + 2 * duty[i] - 1;
    return y;
}

Vec comb_bl(const Vec &t) {
    Vec y(t.size(), 0.0);
    const int N = harmonics(t) + 1;
    for (int h = 1; h < N; ++h)
        for (size_t i = 0; i < t.size(); ++i) y[i] += std::cos(h * t[i]);
    for (double &v : y) v /= N;
    return y;
}

std::vector<int> max_len_seq(int nbits, std::vector<int> &state, int length) {
    static const std::map<int, std::vector<int>> taps_map = {
        {2, {1}}, {3, {2}}, {4, {3}}, {5, {3}}, {6, {5}}, {7, {6}}, {8, {7, 6, 1}}, {9, {5}}, {10, {7}},
        {11, {9}}, {12, {11, 10, 4}}, {13, {12, 11, 8}}, {14, {13, 12, 2}}, {15, {14}}, {16, {15, 13, 4}},
        {17, {14}}, {18, {11}}, {19, {18, 17, 14}}, {20, {17}}, {21, {19}}, {22, {21}}, {23, {18}},
        {24, {23, 22, 17}}, {25, {22}}, {26, {25, 24, 20}}, {27, {26, 25, 22}}, {28, {25}}, {29, {27}},
        {30, {29, 28, 7}}, {31, {28}}, {32, {31, 30, 10}}};
    const auto it = taps_map.find(nbits);
    if (it == taps_map.end()) throw DesignError("The number of MLS bits must be between 2 and 32.");
    if (int(state.size()) != nbits) throw DesignError("MLS state must have nbits elements.");
    for (int &s : state) s = s != 0;
    if (std::all_of(state.begin(), state.end(), [](int s) { return s == 0; }))
        throw DesignError("MLS state must not be all zeros.");
    const std::vector<int> &taps = it->second;
    std::vector<int> seq(size_t(std::max(length, 0)));
    int idx = 0;
    for (int i = 0; i < length; ++i) {
        int feedback = state[size_t(idx)];
        seq[size_t(i)] = feedback;
        for (int tap : taps) feedback ^= state[size_t((tap + idx) % nbits)];
        state[size_t(idx)] = feedback;
        idx = (idx + 1) % nbits;
    }
    std::rotate(state.begin(), state.begin() + idx, state.end());  // np.roll(state, -idx)
    return seq;
}

// ---------------------------------------------------------------------------
// stimulus

double impulse_scale(const StimParams &p) {
    switch (p.stim) {
    case Stim::Dirac: return 1.0;
    case Stim::Sinc: return p.f1 * 2;
    case Stim::Gauss: return p.f1 * 2 * p.bw1;
    case Stim::Rect: return 1.0 / p.tw;
    default: return 1.0;
    }
}

Vec calc_stimulus(const StimParams &p, int n_end) {
    if (n_end < 1) throw DesignError("The number of data points must be >= 1.");
    Vec x(static_cast<size_t>(n_end), 0.0);
    const double phi1 = p.phi1 / 180 * PI, phi2 = p.phi2 / 180 * PI;
    const int t1_idx = int(std::nearbyint(p.t1));
    auto need_f1 = [&] {
        if (p.f1 <= 0) throw DesignError("Frequency f1 needs to be > 0.");
    };
    // phase vector 2 pi f1 n + phi1 for the bandlimited signals
    auto phase1 = [&] {
        Vec t(static_cast<size_t>(n_end));
        for (int n = 0; n < n_end; ++n) t[size_t(n)] = 2 * PI * n * p.f1 + phi1;
        return t;
    };
    const bool bl = p.bl && n_end >= 2;

    switch (p.stim) {
    case Stim::None: break;
    case Stim::Dirac:
        if (t1_idx >= 0 && t1_idx < n_end) x[size_t(t1_idx)] = p.a1;
        break;
    case Stim::Sinc:
        for (int n = 0; n < n_end; ++n)
            x[size_t(n)] = p.a1 * sinc(2 * (n - p.t1) * p.f1) + p.a2 * sinc(2 * (n - p.t2) * p.f2);
        break;
    case Stim::Gauss: {
        if ((p.a1 != 0 && p.f1 < 0) || (p.a2 != 0 && p.f2 < 0))
            throw DesignError("Center frequencies f1, f2 need to be >= 0.");
        const double f1 = p.f1 < 0 ? 0.1 : p.f1, f2 = p.f2 < 0 ? 0.1 : p.f2;
        for (int n = 0; n < n_end; ++n)
            x[size_t(n)] = p.a1 * gausspulse(n - p.t1, f1, p.bw1) + p.a2 * gausspulse(n - p.t2, f2, p.bw2);
        break;
    }
    case Stim::Rect: {
        const double n_rise = double(int(t1_idx - std::floor(p.tw / 2)));
        const double n_min = std::max(n_rise, 0.0), n_max = std::min(n_rise + p.tw, double(n_end));
        for (int n = 0; n < n_end; ++n) x[size_t(n)] = (n >= n_min && n < n_max) ? p.a1 : 0.0;
        break;
    }
    case Stim::Step:
        for (int n = std::max(t1_idx, 0); n < n_end; ++n) x[size_t(n)] = p.a1;
        break;
    case Stim::Cos:
        for (int n = 0; n < n_end; ++n)
            x[size_t(n)] = p.a1 * std::cos(2 * PI * n * p.f1 + phi1) + p.a2 * std::cos(2 * PI * n * p.f2 + phi2);
        break;
    case Stim::Sine:
        for (int n = 0; n < n_end; ++n)
            x[size_t(n)] = p.a1 * std::sin(2 * PI * n * p.f1 + phi1) + p.a2 * std::sin(2 * PI * n * p.f2 + phi2);
        break;
    case Stim::Diric:
        if (p.n1 < 1) throw DesignError("N1 needs to be >= 1.");
        for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * diric(2 * PI * (n - p.t1) * p.f1, p.n1);
        break;
    case Stim::Chirp: {
        const double t_end = p.t2 == 0 ? n_end : p.t2;  // sweep over the whole interval or up to t2
        // pyfda passes the phase in radians to scipy.signal.chirp which expects degrees;
        // here phi1 is used in degrees like for all other stimuli
        for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * chirp(n, p.f1, t_end, p.f2, p.chirp, p.phi1);
        break;
    }
    case Stim::Triang:
        if (bl) {
            need_f1();
            const Vec y = triang_bl(phase1());
            for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * y[size_t(n)];
        } else {
            for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * sawtooth(2 * PI * n * p.f1 + phi1, 0.5);
        }
        break;
    case Stim::Saw:
        if (bl) {
            need_f1();
            const Vec y = sawtooth_bl(phase1());
            for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * y[size_t(n)];
        } else {
            for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * sawtooth(2 * PI * n * p.f1 + phi1);
        }
        break;
    case Stim::RectPer:
        if (bl) {
            need_f1();
            const Vec y = rect_bl(phase1(), Vec(size_t(n_end), p.duty));
            for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * y[size_t(n)];
        } else {
            for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * square(2 * PI * n * p.f1 + phi1, p.duty);
        }
        break;
    case Stim::Comb: {
        need_f1();
        if (n_end < 2) throw DesignError("The comb signal needs at least 2 samples.");
        const Vec y = comb_bl(phase1());
        for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * y[size_t(n)];
        break;
    }
    case Stim::AM:
        for (int n = 0; n < n_end; ++n)
            x[size_t(n)] = p.a1 * std::sin(2 * PI * n * p.f1 + phi1) * p.a2 * std::sin(2 * PI * n * p.f2 + phi2);
        break;
    case Stim::PMFM:
        for (int n = 0; n < n_end; ++n)
            x[size_t(n)] = p.a1 * std::sin(2 * PI * n * p.f1 + phi1 + p.a2 * std::sin(2 * PI * n * p.f2 + phi2));
        break;
    case Stim::PWM: {
        Vec duty(static_cast<size_t>(n_end));
        for (int n = 0; n < n_end; ++n) duty[size_t(n)] = 0.5 + p.a2 / 2 * std::sin(2 * PI * n * p.f2 + phi2);
        if (bl) {
            need_f1();
            const Vec y = rect_bl(phase1(), duty);
            for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * y[size_t(n)];
        } else {
            for (int n = 0; n < n_end; ++n) x[size_t(n)] = p.a1 * square(2 * PI * n * p.f1 + phi1, duty[size_t(n)]);
        }
        break;
    }
    }

    // --- noise
    std::mt19937 rng(p.seed);
    switch (p.noise) {
    case Noise::None: break;
    case Noise::Gauss: {
        std::normal_distribution<double> d(0.0, 1.0);
        for (double &v : x) v += p.noi * d(rng);
        break;
    }
    case Noise::Uniform: {
        std::uniform_real_distribution<double> d(0.0, 1.0);
        for (double &v : x) v += p.noi * (d(rng) - 0.5);
        break;
    }
    case Noise::RandInt: {
        std::uniform_int_distribution<long long> d(0, (long long)(std::fabs(p.noi)));
        for (double &v : x) v += double(d(rng));
        break;
    }
    case Noise::MLS: {
        // fixed seed like pyfda, yielding the same sequence at every run
        static const int seed[] = {1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1,
                                   0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0};
        const int b = std::clamp(p.mls_b, 2, 32);
        std::vector<int> state(seed, seed + b);
        const std::vector<int> s = max_len_seq(b, state, n_end);
        for (int n = 0; n < n_end; ++n) x[size_t(n)] += s[size_t(n)] * p.noi;
        break;
    }
    case Noise::Brownian: {
        std::normal_distribution<double> d(0.0, 1.0);
        double acc = 0;
        for (double &v : x) {
            acc += p.noi * d(rng);
            v += acc;
        }
        break;
    }
    }
    if (stim_uses(p.stim, "dc") && p.dc != 0)
        for (double &v : x) v += p.dc;
    return x;
}

std::string stim_title(const StimParams &p, bool step_error) {
    std::string t;
    switch (p.stim) {
    case Stim::None: t = "Zero Input"; break;
    case Stim::Dirac: t = "Impulse Response"; break;
    case Stim::Sinc: t = "Sinc Impulse"; break;
    case Stim::Gauss: t = "Gaussian Impulse"; break;
    case Stim::Rect: t = "Rect Impulse"; break;
    case Stim::Step: t = step_error ? "Settling Error" : "Step Response"; break;
    case Stim::Cos: t = "Cosine Stimulus"; break;
    case Stim::Sine: t = "Sinusoidal Stimulus"; break;
    case Stim::Diric: t = "Periodic Sinc Stimulus"; break;
    case Stim::Chirp: {
        static const char *names[] = {"Linear", "Quadratic", "Logarithmic", "Hyperbolic"};
        t = std::string(names[int(p.chirp)]) + " Chirp Stimulus";
        break;
    }
    case Stim::Triang: t = p.bl ? "Bandlim. Triangular Stimulus" : "Triangular Stimulus"; break;
    case Stim::Saw: t = p.bl ? "Bandlim. Sawtooth Stimulus" : "Sawtooth Stimulus"; break;
    case Stim::RectPer: t = p.bl ? "Bandlimited Rect. Stimulus" : "Rect. Stimulus"; break;
    case Stim::Comb: t = "Bandlim. Comb Stimulus"; break;
    case Stim::AM: t = "AM Stimulus"; break;
    case Stim::PMFM: t = "PM / FM Stimulus"; break;
    case Stim::PWM: t = "PWM Stimulus"; break;
    }
    switch (p.noise) {
    case Noise::None: break;
    case Noise::Gauss: t += " + Gaussian Noise"; break;
    case Noise::Uniform: t += " + Uniform Noise"; break;
    case Noise::RandInt: t += " + Random Int. Sequence"; break;
    case Noise::MLS: t += " + Max. Length Sequence"; break;
    case Noise::Brownian: t += " + Brownian Noise"; break;
    }
    if (stim_uses(p.stim, "dc") && p.dc != 0) t += " + DC";
    return t;
}

// ---------------------------------------------------------------------------
// spectra

Vec fft_window(WindowType type, int N, double par) {
    if (N < 1) return {};
    Vec w = get_window(type, N + 1, par);
    w.pop_back();
    return w;
}

double window_cgain(const Vec &win) {
    double s = 0;
    for (double v : win) s += v;
    return win.empty() ? 1.0 : s / double(win.size());
}

double window_nenbw(const Vec &win) {
    double s = 0, s2 = 0;
    for (double v : win) {
        s += v;
        s2 += v * v;
    }
    return s == 0 ? 1.0 : double(win.size()) * s2 / (s * s);
}

CVec windowed_fft(const Vec &x, const Vec &win) {
    const size_t N = x.size();
    const double cg = window_cgain(win);
    CVec xw(N);
    for (size_t i = 0; i < N; ++i) xw[i] = x[i] * (win[i] / cg);
    CVec X = fft(xw);
    for (cplx &v : X) v /= double(N);
    return X;
}

CVec ssb_spectrum(const CVec &X) {
    const size_t N = X.size();
    if (N == 0) return {};
    CVec s;
    s.push_back(X[0]);
    for (size_t k = 1; k < N / 2; ++k) s.push_back(X[k] * 2.0);
    return s;
}

}  // namespace pyfda
