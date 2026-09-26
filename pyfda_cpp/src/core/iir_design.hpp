// IIR filter design, ports of the scipy.signal routines used by pyfda
// (butter, cheby1, cheby2, ellip, bessel and the corresponding *ord functions).
// All frequencies are normalized to the Nyquist frequency (0 < W < 1), digital only.
#pragma once

#include "types.hpp"

namespace pyfda {

enum class BType { Lowpass, Highpass, Bandpass, Bandstop };
enum class IirType { Butter, Cheby1, Cheby2, Ellip, Bessel };

/// Analog lowpass prototypes
Zpk buttap(int N);
Zpk cheb1ap(int N, double rp);
Zpk cheb2ap(int N, double rs);
Zpk ellipap(int N, double rp, double rs);
Zpk besselap(int N);  // norm = 'phase'

/// Frequency transformations of analog zpk systems
Zpk lp2lp_zpk(const Zpk &zpk, double wo);
Zpk lp2hp_zpk(const Zpk &zpk, double wo);
Zpk lp2bp_zpk(const Zpk &zpk, double wo, double bw);
Zpk lp2bs_zpk(const Zpk &zpk, double wo, double bw);
Zpk bilinear_zpk(const Zpk &zpk, double fs);

/// scipy.signal.iirfilter(N, Wn, rp, rs, btype, analog=False, ftype, output='zpk')
Zpk iirfilter(int N, const Vec &Wn, double rp, double rs, BType btype, IirType ftype);

struct OrdResult {
    int N = 0;
    Vec Wn;  // one (LP / HP) or two (BP / BS) critical frequencies
};

/// Minimum order estimation, wp / ws: one or two band edges (normalized to f_Ny)
OrdResult buttord(const Vec &wp, const Vec &ws, double gpass, double gstop);
OrdResult cheb1ord(const Vec &wp, const Vec &ws, double gpass, double gstop);
OrdResult cheb2ord(const Vec &wp, const Vec &ws, double gpass, double gstop);
OrdResult ellipord(const Vec &wp, const Vec &ws, double gpass, double gstop);

}  // namespace pyfda
