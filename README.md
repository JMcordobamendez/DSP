# DSP

A small Python helper for looking at a sampled signal in the time and frequency domain, and optionally comparing it with a digitally filtered version.

## Installation

```bash
pip install -r requirements.txt
```

## The `DSP` class

`DSP.py` defines one class:

```python
DSP(data, fs, filt=False, num=np.array([1]), den=np.array([1, 0]))
```

| Argument | Meaning |
|----------|---------|
| `data`   | Sequence of samples (list or NumPy array). |
| `fs`     | Sampling frequency in Hz. |
| `filt`   | `True` to also filter the signal, `False` (default) to skip filtering. The legacy strings `'Yes'`/`'No'` are still accepted. |
| `num`    | Numerator coefficients of the filter (for an FIR filter, the taps). |
| `den`    | Denominator coefficients of the filter. |

Methods:

- `calculate()` computes the one-sided FFT amplitude spectrum of `data` into `YFFT` (DC bin as is, other bins doubled so a sine of amplitude A shows a peak of A). The matching frequency axis is in `x`, with a resolution of `fs / len(data)` Hz. If `filt` is true, it also filters the data with zero-phase `scipy.signal.filtfilt(num, den, data)`, storing the result in `data_f` and its spectrum in `YF_FFT`. It returns `(x, YFFT, YF_FFT)`, with `None` in place of `YF_FFT` when no filter is used.
- `plot()` opens a matplotlib figure with two subplots: the signal over time and its amplitude spectrum, showing the original in blue and, when filtering is enabled, the filtered version in red.

Call `calculate()` before `plot()`; calling it again recomputes the spectra from scratch.

## Example

`DSP_utilizar.py` builds one second of signal sampled at 1000 Hz:

```
y = 1.0 * sin(2π·10·t) + 0.1 * sin(2π·400·t) + 0.5
```

It then applies a 9-tap low-pass FIR filter and plots the result:

```bash
python DSP_utilizar.py
```

The spectrum shows a 0.5 DC component, a peak of 1.0 at 10 Hz and a peak of 0.1 at 400 Hz. The low-pass filter keeps the 10 Hz component and removes most of the 400 Hz one.

## Tests

```bash
pytest
```

`tests/test_dsp.py` checks that the example signal's spectrum peaks at 10 Hz and 400 Hz with the expected amplitudes, and that the example filter attenuates the 400 Hz component.
