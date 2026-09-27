// Filtering of data (ports of scipy.signal lfilter, sosfilt, filtfilt, sosfiltfilt)
// and spectrum calculation
#pragma once

#include "types.hpp"

namespace pyfda {

/// Direct form II transposed, zi: initial state (size max(len(a), len(b)) - 1) or empty
Vec lfilter(const Vec &b, const Vec &a, const Vec &x, const Vec &zi = {}, Vec *zf = nullptr);
Vec lfilter_zi(const Vec &b, const Vec &a);

/// Cascade of second-order sections, zi: n_sections x 2 initial states or empty
Vec sosfilt(const Sos &sos, const Vec &x, const std::vector<std::array<double, 2>> &zi = {});
std::vector<std::array<double, 2>> sosfilt_zi(const Sos &sos);

/// Zero-phase filtering with odd extension of the signal (scipy default padtype='odd')
Vec filtfilt(const Vec &b, const Vec &a, const Vec &x);
Vec sosfiltfilt(const Sos &sos, const Vec &x);

/// FFT of arbitrary length (radix 2, Bluestein for other lengths)
CVec fft(const CVec &x);
/// Single-sided amplitude spectrum |X(f)|, frequencies k * fs / N for k = 0 ... N/2
Vec amplitude_spectrum(const Vec &x);

}  // namespace pyfda
