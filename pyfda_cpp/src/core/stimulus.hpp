// Stimuli for the transient analysis, ported from pyfda
// (pyfda/plot_widgets/tran/plot_tran_stim.py, libs/special_functions.py) and the
// waveforms of scipy.signal (sawtooth, square, chirp, gausspulse, max_len_seq)
// and scipy.special.diric.
//
// All frequencies are normalized to f_S, all times are given in samples.
#pragma once

#include "fir_design.hpp"
#include "types.hpp"

#include <cstdint>
#include <string>

namespace pyfda {

enum class Stim {
    None, Dirac, Sinc, Gauss, Rect, Step,              // impulses / step
    Sine, Cos, Diric,                                  // sinusoids
    Chirp,
    Triang, Saw, RectPer, Comb,                        // periodic
    AM, PMFM, PWM                                      // modulation
};
enum class ChirpType { Linear, Quadratic, Logarithmic, Hyperbolic };
enum class Noise { None, Gauss, Uniform, RandInt, MLS, Brownian };

struct StimInfo {
    Stim stim;
    const char *key;   // pyfda name, e.g. "rect_per"
    const char *name;  // display name
    const char *params;  // parameters used, subset of "a1 a2 f1 f2 phi1 phi2 t1 t2 tw bw1 bw2 n1 bl duty dc"
};
const std::vector<StimInfo> &stim_list();
const StimInfo &stim_info(Stim s);
bool stim_uses(Stim s, const std::string &param);
/// Impulse shaped stimuli (dirac, sinc, gauss, rect): the spectrum of the
/// response can be scaled to show the frequency response
bool is_impulse(Stim s);

struct StimParams {
    Stim stim = Stim::Dirac;
    double a1 = 1.0, a2 = 0.0;
    double f1 = 0.02, f2 = 0.03;      // normalized to f_S
    double phi1 = 0.0, phi2 = 0.0;    // degrees
    double t1 = 0.0, t2 = 0.0;        // samples
    double tw = 10.0;                 // width of the rect impulse in samples
    double bw1 = 0.5, bw2 = 0.5;      // relative bandwidth of Gaussian pulses
    int n1 = 5;                       // Dirichlet function order
    double duty = 0.5;                // duty cycle of periodic rect
    bool bl = true;                   // bandlimited periodic signals
    ChirpType chirp = ChirpType::Linear;
    Noise noise = Noise::None;
    double noi = 0.1;                 // noise amplitude (std. dev., peak-to-peak, max. int)
    int mls_b = 8;                    // bits of the maximum length sequence
    double dc = 0.0;
    uint32_t seed = 1;                // seed for random noise
};

/// Stimulus x[n], n = 0 ... n_end - 1; throws DesignError for invalid parameters
Vec calc_stimulus(const StimParams &p, int n_end);
/// Energy scaling of impulses (pyfda scale_impz): 1 for dirac, 2 f1 for sinc, ...
double impulse_scale(const StimParams &p);
/// Short title for plots, e.g. "Impulse Response" or "Sinusoidal Stimulus + Gaussian Noise"
std::string stim_title(const StimParams &p, bool step_error = false);

// --- waveforms (scipy.signal / scipy.special / pyfda.libs.special_functions)
double sawtooth(double t, double width = 1.0);
double square(double t, double duty = 0.5);
double diric(double x, int n);
double gausspulse(double t, double fc, double bw = 0.5, double bwr = -6.0);
/// scipy.signal.chirp, phi in degrees
double chirp(double t, double f0, double t1, double f1, ChirpType method, double phi = 0.0);
/// Bandlimited waveforms by Fourier synthesis, t = 2 pi f n + phi (at least 2 values,
/// the number of harmonics is derived from t[1] - t[0] like in pyfda)
Vec sawtooth_bl(const Vec &t);
Vec triang_bl(const Vec &t);
Vec rect_bl(const Vec &t, const Vec &duty);
Vec comb_bl(const Vec &t);
/// scipy.signal.max_len_seq with the default taps, state is updated
std::vector<int> max_len_seq(int nbits, std::vector<int> &state, int length);

// --- spectra for the transient analysis
/// Periodic window for spectral analysis (scipy.signal.get_window(..., fftbins=True))
Vec fft_window(WindowType type, int N, double par = 0.0);
/// Coherent gain (mean of window) and normalized equivalent noise bandwidth
double window_cgain(const Vec &win);
double window_nenbw(const Vec &win);
/// fft(x * win / cgain) / N like pyfda's calc_fft
CVec windowed_fft(const Vec &x, const Vec &win);
/// Single-sided spectrum from a double-sided one (pyfda calc_ssb_spectrum, mag = False)
CVec ssb_spectrum(const CVec &X);

}  // namespace pyfda
