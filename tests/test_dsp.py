import os
import sys

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from DSP import DSP  # noqa: E402

FS = 1000
# Same 9-tap low-pass FIR used in DSP_utilizar.py
FILT = (0.5 / 0.537) * np.array([
    -0.02062630863958946,
    -0.010688729607911805,
    0.09988368103096376,
    0.2805270356057702,
    0.3730973942070046,
    0.2805270356057702,
    0.09988368103096376,
    -0.010688729607911805,
    -0.02062630863958946,
])


def example_signal():
    # Same signal as DSP_utilizar.py: 10 Hz (amp 1) + 400 Hz (amp 0.1) + 0.5 DC
    i = np.arange(1000)
    return 0.1 * np.sin(2 * np.pi * 400 / FS * i) + np.sin(2 * np.pi * 10 / FS * i) + 0.5


def amplitude_at(freqs, spectrum, f):
    return spectrum[int(np.argmin(np.abs(np.asarray(freqs) - f)))]


def test_spectrum_peaks_at_10_and_400_hz():
    dsp = DSP(data=example_signal(), fs=FS)
    dsp.calculate()
    spectrum = np.array(dsp.YFFT)

    # The two strongest non-DC bins are 10 Hz and 400 Hz
    top_two = sorted(np.array(dsp.x)[1:][np.argsort(spectrum[1:])[-2:]])
    assert top_two == pytest.approx([10, 400])

    assert amplitude_at(dsp.x, spectrum, 10) == pytest.approx(1.0, abs=1e-6)
    assert amplitude_at(dsp.x, spectrum, 400) == pytest.approx(0.1, abs=1e-6)
    assert spectrum[0] == pytest.approx(0.5, abs=1e-6)


def test_low_pass_filter_attenuates_400_hz():
    dsp = DSP(data=example_signal(), fs=FS, filt="Yes", num=FILT)
    dsp.calculate()
    original = np.array(dsp.YFFT)
    filtered = np.array(dsp.YF_FFT)

    # 10 Hz passes almost unchanged, 400 Hz is strongly attenuated
    assert amplitude_at(dsp.x, filtered, 10) == pytest.approx(amplitude_at(dsp.x, original, 10), rel=0.05)
    assert amplitude_at(dsp.x, filtered, 400) < 0.1 * amplitude_at(dsp.x, original, 400)
