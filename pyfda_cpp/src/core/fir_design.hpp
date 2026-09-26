// FIR filter design: windows, window method (firwin), Parks-McClellan (remez)
// and the minimum order estimates used by pyfda.
#pragma once

#include "types.hpp"

#include <string>

namespace pyfda {

enum class WindowType {
    Rectangular, Bartlett, Hann, Hamming, Blackman, BlackmanHarris, Nuttall, Flattop,
    Kaiser, Gaussian, Tukey,
    Barthann, Bohman, Cosine, Parzen, Triang, Chebwin, DPSS,
    GeneralGaussian  // spectral analysis only, like in pyfda
};

struct WindowInfo {
    WindowType type;
    const char *name;       // display name
    const char *par_name;   // name of the parameter or nullptr
    double par_default;
    const char *par2_name = nullptr;  // second parameter or nullptr
    double par2_default = 0;
    bool fir = true;        // available for the FIR window design (else spectral analysis only)
};

/// All available windows, the ones for the FIR design (fir = true) first
const std::vector<WindowInfo> &window_list();
/// Symmetric window of length M (scipy.signal.get_window(..., fftbins=False))
Vec get_window(WindowType type, int M, double par = 0.0, double par2 = 0.0);

/// Window method, cutoff normalized to f_Ny (like scipy.signal.firwin with nyq = 1)
Vec firwin(int numtaps, const Vec &cutoff, const Vec &window, bool pass_zero, bool scale = true);

/// scipy.signal.kaiserord / kaiser_beta
void kaiserord(double ripple_db, double width, int &numtaps, double &beta);
double kaiser_beta(double a);

enum class RemezAlg { Ichige, Kaiser, Herrmann };
/// Minimum filter length for remez (pyfda.filter_widgets.common.remezord), freqs relative to fs
struct RemezOrd {
    int numtaps = 0;
    Vec bands, desired, weight;
};
RemezOrd remezord(const Vec &freqs, const Vec &amps, const Vec &rips, double fs = 1.0,
                  RemezAlg alg = RemezAlg::Ichige);

enum class RemezType { Bandpass = 1, Differentiator = 2, Hilbert = 3 };
/// scipy.signal.remez
Vec remez(int numtaps, const Vec &bands, const Vec &desired, const Vec &weight,
          RemezType type = RemezType::Bandpass, double fs = 1.0, int maxiter = 25,
          int grid_density = 16);

}  // namespace pyfda
