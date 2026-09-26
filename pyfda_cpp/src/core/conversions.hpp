// Conversions between filter representations and frequency / time responses
#pragma once

#include "types.hpp"

namespace pyfda {

Ba zpk2tf(const Zpk &zpk);
/// Digital zpk -> second-order sections with scipy's 'nearest' pairing
Sos zpk2sos(const Zpk &zpk);
Ba sos2tf(const Sos &sos);
Zpk tf2zpk(const Ba &ba);
Zpk sos2zpk(const Sos &sos);

/// Complex frequency response at angular frequencies w (rad / sample)
CVec freqz(const Ba &ba, const Vec &w);
CVec freqz(const Sos &sos, const Vec &w);

/// Group delay in samples at angular frequencies w (rad / sample),
/// NaN where the response is (numerically) zero
Vec group_delay(const Ba &ba, const Vec &w);
Vec group_delay(const Sos &sos, const Vec &w);

/// Suitable length for displaying the impulse response of an IIR filter
int impz_len(const Zpk &zpk, int n_min = 20, int n_max = 20000);

/// Unwrap a phase vector (numpy.unwrap)
Vec unwrap(const Vec &phi);

}  // namespace pyfda
