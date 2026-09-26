#!/usr/bin/env python3
"""
Compare the results of the pyfda C++ core (via build/pyfda_cli) with scipy.signal.

Usage: python tests/verify_scipy.py path/to/pyfda_cli
Exit code 0 when all checks pass.
"""
import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
import scipy.signal as sig

CLI = sys.argv[1] if len(sys.argv) > 1 else os.path.join("build", "pyfda_cli")
# the CLI writes UTF-8 (column names), independent of the console code page
proc = subprocess.Popen([CLI], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, encoding="utf-8")
n_pass = n_fail = 0
rng = np.random.default_rng(1)


def call(line: str) -> dict:
    proc.stdin.write(line + "\n")
    proc.stdin.flush()
    return json.loads(proc.stdout.readline())


def L(v) -> str:
    return ",".join(repr(float(x)) for x in np.ravel(v))


def cplx(v):
    return np.array([complex(a, b) for a, b in v])


def check(name: str, ok: bool, detail: str = "") -> None:
    global n_pass, n_fail
    if ok:
        n_pass += 1
    else:
        n_fail += 1
        print(f"FAIL {name} {detail}")


def close(name, got, ref, rtol=1e-9, atol=0.0):
    got, ref = np.asarray(got), np.asarray(ref)
    if got.shape != ref.shape:
        check(name, False, f"shape {got.shape} != {ref.shape}")
        return
    scale = max(np.max(np.abs(ref)), 1e-300) if ref.size else 1
    err = np.max(np.abs(got - ref)) / scale if ref.size else 0
    check(name, err <= rtol + atol / scale, f"rel. error {err:.3g}")


def resp_close(name, b, a, b_ref, a_ref, tol_db=1e-6, n=2048):
    """ compare magnitude responses in dB (relative to max), limited to -250 dB """
    w = np.linspace(0, np.pi, n, endpoint=False)
    h = np.abs(sig.freqz(b, a, w)[1])
    h_ref = np.abs(sig.freqz(b_ref, a_ref, w)[1])
    db = 20 * np.log10(np.maximum(h, 1e-300))
    db_ref = 20 * np.log10(np.maximum(h_ref, 1e-300))
    mask = db_ref > np.max(db_ref) - 250
    err = np.max(np.abs(db - db_ref)[mask])
    check(name, err <= tol_db, f"max. dB error {err:.3g}")


# ---------------------------------------------------------------------------
# IIR prototypes and transformations
IIR = {'butter': lambda N, rp, rs, wn, bt: sig.butter(N, wn, bt, output='zpk'),
       'cheby1': lambda N, rp, rs, wn, bt: sig.cheby1(N, rp, wn, bt, output='zpk'),
       'cheby2': lambda N, rp, rs, wn, bt: sig.cheby2(N, rs, wn, bt, output='zpk'),
       'ellip': lambda N, rp, rs, wn, bt: sig.ellip(N, rp, rs, wn, bt, output='zpk'),
       'bessel': lambda N, rp, rs, wn, bt: sig.bessel(N, wn, bt, output='zpk')}
for ft in IIR:
    for bt, wn in (('low', [0.2]), ('high', [0.35]), ('bandpass', [0.2, 0.5]),
                   ('bandstop', [0.3, 0.45])):
        for N in (1, 2, 3, 4, 5, 8, 11):
            rp, rs = 0.5, 55
            r = call(f"iir {ft} {bt} {N} {rp} {rs} {L(wn)}")
            if 'error' in r:
                check(f"iir {ft} {bt} N={N}", False, r['error'])
                continue
            z, p, k = IIR[ft](N, rp, rs, wn if len(wn) > 1 else wn[0], bt)
            name = f"iir {ft} {bt} N={N}"
            close(name + " k", r['k'], k, 1e-9)
            close(name + " poles", np.sort_complex(cplx(r['p'])), np.sort_complex(p), 1e-9)
            if len(z):
                close(name + " zeros", np.sort_complex(cplx(r['z'])), np.sort_complex(z), 1e-8)
            sos_ref = sig.zpk2sos(z, p, k)
            close(name + " sos", r['sos'], sos_ref, 1e-8)
            b_ref, a_ref = sig.zpk2tf(z, p, k)
            close(name + " b", r['b'], b_ref, 1e-8)
            close(name + " a", r['a'], a_ref, 1e-8)
            w = np.linspace(0, np.pi, 1024, endpoint=False)
            h = np.abs(sig.sosfreqz(r['sos'], w)[1])
            h_ref = np.abs(sig.sosfreqz(sos_ref, w)[1])
            close(name + " |H(sos)|", h, h_ref, 1e-9, 1e-12)

# minimum order
ORD = {'butter': sig.buttord, 'cheby1': sig.cheb1ord, 'cheby2': sig.cheb2ord,
       'ellip': sig.ellipord}
for ft, fn in ORD.items():
    for wp, ws in (([0.2], [0.3]), ([0.3], [0.2]), ([0.2, 0.5], [0.1, 0.6]),
                   ([0.1, 0.6], [0.2, 0.5]), ([0.05], [0.07]), ([0.4, 0.42], [0.35, 0.47])):
        for gp, gs in ((1, 40), (0.1, 80), (3, 20)):
            r = call(f"ord {ft} {L(wp)} {L(ws)} {gp} {gs}")
            N, Wn = fn(wp if len(wp) > 1 else wp[0], ws if len(ws) > 1 else ws[0], gp, gs)
            name = f"ord {ft} wp={wp} ws={ws} {gp}/{gs}"
            check(name + " N", r.get('N') == N, f"{r} vs {N}")
            close(name + " Wn", r.get('Wn', []), np.atleast_1d(Wn), 1e-7)

# ---------------------------------------------------------------------------
# windows, firwin, kaiserord, remez
WINDOWS = {'rectangular': lambda M, p: sig.get_window('boxcar', M, False),
           'bartlett': lambda M, p: sig.get_window('bartlett', M, False),
           'hann': lambda M, p: sig.get_window('hann', M, False),
           'hamming': lambda M, p: sig.get_window('hamming', M, False),
           'blackman': lambda M, p: sig.get_window('blackman', M, False),
           'blackman-harris': lambda M, p: sig.get_window('blackmanharris', M, False),
           'nuttall': lambda M, p: sig.get_window('nuttall', M, False),
           'flattop': lambda M, p: sig.get_window('flattop', M, False),
           'kaiser': lambda M, p: sig.get_window(('kaiser', p), M, False),
           'gaussian': lambda M, p: sig.get_window(('gaussian', p * (M - 1) / 2), M, False),
           'tukey': lambda M, p: sig.get_window(('tukey', p), M, False),
           'bartlett-hann': lambda M, p: sig.get_window('barthann', M, False),
           'bohman': lambda M, p: sig.get_window('bohman', M, False),
           'cosine': lambda M, p: sig.get_window('cosine', M, False),
           'parzen': lambda M, p: sig.get_window('parzen', M, False),
           'triangular': lambda M, p: sig.get_window('triang', M, False),
           'dolph-chebyshev': lambda M, p: sig.get_window(('chebwin', p), M, False),
           'dpss': lambda M, p: sig.get_window(('dpss', p), M, False)}
for wname, fn in WINDOWS.items():
    for M in (2, 7, 32, 101, 256):
        par = {'kaiser': 6.5, 'gaussian': 0.4, 'tukey': 0.3, 'dolph-chebyshev': 80, 'dpss': 3.0}.get(wname, 0)
        if wname == 'dpss' and M == 2:
            par = 0.5
        r = call(f"window {wname} {M} {par}")
        close(f"window {wname} M={M}", r['w'], fn(M, par), 1e-9 if wname == 'dpss' else 1e-12, 1e-15)

for numtaps, cutoff, pz in ((31, [0.3], 1), (32, [0.3], 1), (41, [0.4], 0),
                            (51, [0.2, 0.5], 0), (61, [0.2, 0.5], 1), (80, [0.1, 0.3], 0)):
    for wname in ('hamming', 'kaiser', 'blackman'):
        par = 8.0 if wname == 'kaiser' else 0
        r = call(f"firwin {numtaps} {L(cutoff)} {pz} {wname} {par}")
        win = ('kaiser', par) if wname == 'kaiser' else wname
        h = sig.firwin(numtaps, cutoff, window=win, pass_zero=bool(pz))
        close(f"firwin {numtaps} {cutoff} {pz} {wname}", r['h'], h, 1e-12)

for ripple, width in ((60, 0.05), (40, 0.1), (21, 0.2), (90, 0.01)):
    r = call(f"kaiserord {ripple} {width}")
    N, beta = sig.kaiserord(ripple, width)
    check(f"kaiserord {ripple} {width}", r['N'] == N and math.isclose(r['beta'], beta, rel_tol=1e-12),
          f"{r} vs {N}, {beta}")

REMEZ = [
    (35, [0, 0.1, 0.15, 0.5], [1, 0], [1, 1], 'bandpass'),
    (36, [0, 0.1, 0.15, 0.5], [1, 0], [1, 10], 'bandpass'),
    (51, [0, 0.1, 0.15, 0.3, 0.35, 0.5], [0, 1, 0], [10, 1, 10], 'bandpass'),
    (61, [0, 0.1, 0.15, 0.3, 0.35, 0.5], [1, 0, 1], [1, 5, 1], 'bandpass'),
    (40, [0, 0.2, 0.25, 0.5], [0, 1], [1, 1], 'hilbert'),
    (31, [0.05, 0.45], [1], [1], 'hilbert'),
    (101, [0, 0.05, 0.07, 0.5], [1, 0], [1, 100], 'bandpass'),
    (9, [0, 0.1, 0.3, 0.5], [1, 0], [1, 1], 'bandpass'),
    (201, [0, 0.2, 0.21, 0.5], [1, 0], [1, 1], 'bandpass'),
]
for numtaps, bands, des, w, t in REMEZ:
    name = f"remez {numtaps} {bands} {t}"
    r = call(f"remez {numtaps} {L(bands)} {L(des)} {L(w)} {t} 16")
    try:
        h = sig.remez(numtaps, bands, des, weight=w, type=t, fs=1)
    except ValueError as e:
        check(name, 'error' in r, f"scipy: {e}, C++: {r}")
        continue
    if 'error' in r:
        check(name, False, r['error'])
        continue
    close(name, r['h'], h, 1e-9)

# ---------------------------------------------------------------------------
# pyfda's remezord (filter_widgets/common.py), reimplemented here from the source
sys.path.insert(0, os.environ.get("PYFDA_SRC", ""))


def remezord_py(freqs, amps, rips, alg):
    def herrmann(fp, fs, dp, ds):
        df = fs - fp
        a = [5.309e-3, 7.114e-2, -4.761e-1, -2.66e-3, -5.941e-1, -4.278e-1]
        b = [11.01217, 0.51244]
        dinf = np.log10(ds) * (a[0] * np.log10(dp)**2 + a[1] * np.log10(dp) + a[2])\
            + a[3] * np.log10(dp)**2 + a[4] * np.log10(dp) + a[5]
        f = b[0] + b[1] * (np.log10(dp) - np.log10(ds))
        return int(dinf / df - f * df + 1)

    def kaiser(fp, fs, dp, ds):
        return int((-20*np.log10(np.sqrt(dp*ds))-13.0)/(14.6*(fs-fp))+1.0)

    def ichige(fp, fs, dp, ds):
        def v(df, dp):
            return 2.325 * ((-np.log10(dp))**-0.445) * df ** (-1.39)

        def g(df, fp):
            return (2.0 / np.pi) * np.arctan(v(df, dp) * (1.0 / fp - 1.0 / (0.5 - df)))

        def h(df, fp, c):
            return (2.0/np.pi) * np.arctan((c/df)*(1.0/fp-1.0/(0.5-df)))
        df = fs-fp
        nc = np.ceil(1.0+(1.101/df) * (-np.log10(2.0*dp)) ** 1.1)
        nm = (0.52/df)*np.log10(dp/ds)*(-np.log10(dp))**0.17
        n3 = np.ceil(nc*(g(df, fp) + g(df, 0.5-df-fp) + 1.0) / 3.0)
        dn = np.ceil(nm*(h(df, fp, 1.1) - (h(df, 0.5-df-fp, 0.29) - 1.0) / 2.0))
        return int(n3 + dn)
    fn = {'ichige': ichige, 'kaiser': kaiser, 'herrmann': herrmann}[alg]
    freqs, amps, rips = (np.asarray(x, 'd') for x in (freqs, amps, rips))
    rips = rips / (amps + (amps == 0.0))
    f1, f2 = freqs[0:-1:2], freqs[1::2]
    n = 0
    for i in range(len(amps) - 1):
        n = max(n, fn(f1[i], f2[i], rips[i], rips[i+1]), fn(0.5-f2[i], 0.5-f1[i], rips[i+1], rips[i]))
    return n, np.hstack((0.0, freqs, 0.5)), max(rips) / rips


for alg in ('ichige', 'kaiser', 'herrmann'):
    for f, a, rp in (([0.1, 0.15], [1, 0], [0.01, 0.001]), ([0.2, 0.25], [0, 1], [0.001, 0.05]),
                     ([0.1, 0.15, 0.3, 0.35], [0, 1, 0], [0.001, 0.01, 0.001])):
        r = call(f"remezord {alg} {L(f)} {L(a)} {L(rp)}")
        n, bands, w = remezord_py(f, a, rp, alg)
        check(f"remezord {alg} {f}", r['N'] == n, f"{r['N']} vs {n}")
        close(f"remezord {alg} {f} weight", r['weight'], w, 1e-12)

# ---------------------------------------------------------------------------
# filtering
x = rng.standard_normal(500) + np.sin(np.arange(500) * 0.05) * 3 + 1.5
for ft in ('butter', 'ellip', 'cheby2'):
    for bt, wn in (('low', 0.1), ('high', 0.3), ('bandpass', [0.1, 0.3]), ('bandstop', [0.2, 0.4])):
        z, p, k = IIR[ft](4, 1, 50, wn, bt)
        sos = sig.zpk2sos(z, p, k)
        b, a = sig.zpk2tf(z, p, k)
        name = f"filter {ft} {bt}"
        close(name + " sosfilt", call(f"sosfilter sosfilt {L(sos)} {L(x)}")['y'], sig.sosfilt(sos, x), 1e-10)
        close(name + " sosfiltfilt", call(f"sosfilter sosfiltfilt {L(sos)} {L(x)}")['y'],
              sig.sosfiltfilt(sos, x), 1e-10)
        close(name + " lfilter", call(f"filter lfilter {L(b)} {L(a)} {L(x)}")['y'], sig.lfilter(b, a, x), 1e-8)
        close(name + " filtfilt", call(f"filter filtfilt {L(b)} {L(a)} {L(x)}")['y'], sig.filtfilt(b, a, x), 1e-8)
h = sig.firwin(41, 0.2)
close("filter FIR lfilter", call(f"filter lfilter {L(h)} 1 {L(x)}")['y'], sig.lfilter(h, 1, x), 1e-12)
close("filter FIR filtfilt", call(f"filter filtfilt {L(h)} 1 {L(x)}")['y'], sig.filtfilt(h, 1, x), 1e-12)

# frequency response, group delay, roots, spectrum
b, a = sig.ellip(6, 1, 60, 0.2)
w = np.linspace(0.01, np.pi - 0.01, 300)
H = cplx(call(f"freqz {L(b)} {L(a)} {L(w)}")['H'])
close("freqz", H, sig.freqz(b, a, w)[1], 1e-9)
sos = sig.ellip(6, 1, 60, 0.2, output='sos')
gd = call(f"gdsos {L(sos)} {L(w)}")['gd']
close("group delay sos", gd, sig.group_delay((b, a), w)[1], 1e-6)
gd = call(f"gd {L(h)} 1 {L(w)}")['gd']
close("group delay FIR", gd, np.full_like(w, 20.0), 1e-9)
# (multiple roots are left out, they are ill-conditioned for any algorithm)
for c in (sig.firwin(31, 0.3), sig.remez(25, [0, 0.2, 0.3, 0.5], [0, 1]), sig.firwin(151, 0.1),
          np.poly([0.5, -0.2 + 0.3j, -0.2 - 0.3j, 2.0]), [2, 0, 0, 1, 0, 0]):
    r = cplx(call(f"roots {L(c)}")['r'])
    ref = np.roots(c)
    close(f"roots deg {len(c) - 1}", np.sort_complex(np.round(r, 9)), np.sort_complex(np.round(ref, 9)), 1e-8)
    if len(c) < 40:  # reconstruction is ill-conditioned for high orders (also with numpy)
        close(f"roots deg {len(c) - 1} poly", np.real(np.poly(r)) * c[0], c, 1e-9)
for n in (256, 300, 1001):
    xs = rng.standard_normal(n)
    A = np.abs(np.fft.rfft(xs)) / n
    A[1:] *= 2
    close(f"spectrum N={n}", call(f"spectrum {L(xs)}")['A'], A, 1e-10)

# ---------------------------------------------------------------------------
# high level design (pyfda filter widgets logic, with N = filter order)
cases = [
    ("rt=LP method=ellip fo=min f_pb=0.1 f_sb=0.15 A_PB=1 A_SB=60",
     lambda: sig.ellip(*sig.ellipord(0.2, 0.3, 1, 60)[:1], 1, 60, sig.ellipord(0.2, 0.3, 1, 60)[1], output='sos')),
    ("rt=HP method=cheby1 fo=min f_pb=0.2 f_sb=0.15 A_PB=0.5 A_SB=50",
     lambda: sig.cheby1(sig.cheb1ord(0.4, 0.3, 0.5, 50)[0], 0.5, sig.cheb1ord(0.4, 0.3, 0.5, 50)[1], 'high', output='sos')),
    ("rt=BP method=butter fo=min f_sb=0.1 f_pb=0.15 f_pb2=0.3 f_sb2=0.35 A_PB=1 A_SB=40",
     lambda: sig.butter(*sig.buttord([0.3, 0.6], [0.2, 0.7], 1, 40), 'bandpass', output='sos')),
    ("rt=BS method=cheby2 fo=min f_pb=0.1 f_sb=0.15 f_sb2=0.3 f_pb2=0.35 A_PB=1 A_SB=40",
     lambda: sig.cheby2(sig.cheb2ord([0.2, 0.7], [0.3, 0.6], 1, 40)[0], 40, sig.cheb2ord([0.2, 0.7], [0.3, 0.6], 1, 40)[1], 'bandstop', output='sos')),
    ("rt=BS method=ellip fo=min f_s=1000 f_pb=100 f_sb=150 f_sb2=300 f_pb2=350 A_PB=1 A_SB=60",
     lambda: sig.ellip(sig.ellipord([0.2, 0.7], [0.3, 0.6], 1, 60)[0], 1, 60, sig.ellipord([0.2, 0.7], [0.3, 0.6], 1, 60)[1], 'bandstop', output='sos')),
    ("rt=LP method=bessel fo=man N=5 f_c=0.1",
     lambda: sig.bessel(5, 0.2, output='sos')),
    ("rt=BP method=ellip fo=man N=8 f_c=0.1 f_c2=0.2 A_PB=1 A_SB=50",
     lambda: sig.ellip(4, 1, 50, [0.2, 0.4], 'bandpass', output='sos')),
]
for spec, ref in cases:
    r = call("design " + spec)
    if 'error' in r:
        check("design " + spec, False, r['error'])
        continue
    sos_ref = ref()
    close("design " + spec, r['sos'], sos_ref, 1e-8)

# FIR designs
fir_cases = [
    ("rt=LP method=equiripple fo=min f_pb=0.1 f_sb=0.15 A_PB=1 A_SB=60", 'LP'),
    ("rt=BP method=equiripple fo=min f_sb=0.1 f_pb=0.15 f_pb2=0.3 f_sb2=0.35 A_PB=1 A_SB=50", 'BP'),
    ("rt=HP method=equiripple fo=man N=40 f_sb=0.2 f_pb=0.25 W_PB=1 W_SB=1", 'HPman'),
    ("rt=LP method=firwin fo=min window=kaiser f_pb=0.1 f_sb=0.15 A_SB=60", 'kaiser'),
    ("rt=HP method=firwin fo=man window=hamming N=30 f_c=0.2", 'hamming'),
]
for spec, kind in fir_cases:
    r = call("design " + spec)
    if 'error' in r:
        check("design " + spec, False, r['error'])
        continue
    ntaps = r['N'] + 1
    a_pb = (10**(1/20) - 1) / (10**(1/20) + 1)
    if kind == 'LP':
        n, bands, w = remezord_py([0.1, 0.15], [1, 0], [a_pb, 10**-3], 'ichige')
        ref = sig.remez(n, bands, [1, 0], weight=w, fs=1)
    elif kind == 'BP':
        n, bands, w = remezord_py([0.1, 0.15, 0.3, 0.35], [0, 1, 0], [10**-2.5, a_pb, 10**-2.5], 'ichige')
        ref = sig.remez(n, bands, [0, 1, 0], weight=w, fs=1)
    elif kind == 'HPman':
        ref = sig.remez(41, [0, 0.2, 0.25, 0.5], [0, 1], fs=1)
    elif kind == 'kaiser':
        n, beta = sig.kaiserord(60, 0.1)
        ref = sig.firwin(n, 0.125, window=('kaiser', beta), fs=1)
    else:
        ref = sig.firwin(31, 0.2, window='hamming', pass_zero=False, fs=1)
    close("design " + spec, r['b'], ref, 1e-9)

# ---------------------------------------------------------------------------
# data import: compare with read_text_table() of pyfda/plot_widgets/plot_data_filt.py
# (extracted from the Python source when PYFDA_DATA_FILT points to it, else built-in copy)
def load_py_reader():
    import ast
    import csv
    import logging
    src_file = os.environ.get("PYFDA_DATA_FILT", "")
    if not os.path.isfile(src_file):
        return None
    src = open(src_file, encoding='utf-8').read()
    tree = ast.parse(src)
    keep = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.Assign))
            and (not isinstance(n, ast.FunctionDef) or n.name in ('_str2num', '_is_num', 'read_text_table'))]
    logging.getLogger('x').setLevel(logging.ERROR)
    ns = {'np': np, 'csv': csv, 'logging': logging, '__name__': 'x', 'os': os}
    exec(compile(ast.Module(body=keep, type_ignores=[]), src_file, 'exec'), ns)
    return ns['read_text_table']


read_text_table_py = load_py_reader()


def _isnum(c):
    c = str(c).strip()
    if ',' in c and '.' not in c:
        c = c.replace(',', '.')
    try:
        float(c)
        return True
    except ValueError:
        return False
if read_text_table_py is None:
    print("PYFDA_DATA_FILT not set, skipping comparison with the Python CSV reader")
csv_cases = [
    "time,ch1\ns,V\n0,1.5\n0.001,2.5\n0.002,-1\n",
    "# comment\n% another\nt;x\n0;1,5\n1;2,5\n2;3\n",
    "Meas 1\n1\n2\n3\n",
    "value\n0,5\n1,5\n2,5\n",
    "a\tb\tc\n1\t2\t3\n4\t5\t6\n",
    "Device: X\nDate: today\nt,y\n0,1\n1,2\n2,3\n3,4\n",
    "Device: X\nt,y\n" + "".join(f"{i},{i*i}\n" for i in range(20)),
    "1 2\n3 4\n5 6\n",
    "x,y,\n1,2,\n3,4,\n",
    "t;u\n0;1\n1;;\n2;x\n3;4\n",
    '"a b","c"\n"1","2"\n3,4\n',
    "// header\nA|B\n1|2\n3|4\n",
]
with tempfile.TemporaryDirectory() as tmp:
    for i, text in enumerate(csv_cases):
        fn = os.path.join(tmp, f"t{i}.csv")
        with open(fn, 'w', encoding='utf-8') as f:
            f.write(text)
        r = call(f"csv {fn}")
        if read_text_table_py is None:
            check(f"csv case {i}", 'error' not in r, str(r))
            continue
        data, names = read_text_table_py(fn)
        ok = 'error' not in r and r['rows'] == data.shape[0] and r['cols'] == data.shape[1]\
            and r['names'] == names
        if ok:
            got = np.array(r['values'], dtype=float).reshape(data.shape)
            conv = np.array([[float(str(c).replace(',', '.')) if _isnum(c) else np.nan for c in row]
                             for row in data])
            ok = np.allclose(got, conv, equal_nan=True)
        check(f"csv case {i}", ok, f"C++ {r} / Python {names} {data.tolist()}")
    # Excel style: cp1252 encoded header with BOM-less umlaut, CRLF
    fn = os.path.join(tmp, "cp1252.csv")
    with open(fn, 'wb') as f:
        f.write("Zeit;Spannung Ä\r\n0;1,0\r\n1;2,0\r\n".encode('cp1252'))
    r = call(f"csv {fn}")
    check("csv cp1252", r.get('names') == ["Zeit", "Spannung Ä"], f"{r}")
    # wav
    import scipy.io.wavfile as wavfile
    fn = os.path.join(tmp, "t.wav")
    data = (np.sin(np.arange(100) * 0.1) * 20000).astype(np.int16)
    wavfile.write(fn, 8000, np.column_stack((data, -data)))
    r = call(f"csv {fn}")
    check("wav", r.get('rows') == 100 and r.get('cols') == 2 and r.get('fs') == 8000
          and np.allclose(np.array(r['values'])[::2], data / 32768), f"{str(r)[:200]}")
    fn = os.path.join(tmp, "t.npy")
    np.save(fn, np.arange(12.).reshape(4, 3))
    r = call(f"csv {fn}")
    check("npy", r.get('rows') == 4 and r.get('cols') == 3 and r['values'] == list(np.arange(12.)), f"{r}")

# ---------------------------------------------------------------------------
# filter files (JSON) and coefficient export: everything must read back exactly
import re
with tempfile.TemporaryDirectory() as tmp:
    specs = ["rt=LP method=ellip", "rt=HP method=cheby2 fo=manual N=7 f_c=0.2",
             "rt=BP method=butter f_pb=0.15 f_pb2=0.3 f_sb=0.1 f_sb2=0.35",
             "rt=BS method=equiripple f_pb=0.1 f_pb2=0.35 f_sb=0.15 f_sb2=0.3",
             "rt=LP method=firwin window=hann fo=manual N=30 f_s=1000 f_c=100 f_pb=80 f_sb=150",
             "rt=HP method=firwin window=kaiser f_pb=0.2 f_sb=0.15"]
    for i, spec in enumerate(specs):
        ref = call("design " + spec)
        j = call("tojson kHz " + spec)
        doc = json.loads(j['json'])
        ok = doc['format'] == "pyfda_cpp filter" and doc['spec']['unit'] == "kHz" \
            and doc['b'] == ref['b'] and doc['a'] == ref['a'] and doc['sos'] == ref['sos']
        check(f"tojson {spec}", ok, str(doc)[:300])
        fn = os.path.join(tmp, f"f{i}.json")
        with open(fn, 'w', encoding='utf-8') as f:
            f.write(j['json'])
        r = call(f"fromjson {fn}")
        ok = 'error' not in r and r['b'] == ref['b'] and r['a'] == ref['a'] and r['sos'] == ref['sos'] \
            and r['file_b'] == ref['b'] and r['file_sos'] == ref['sos'] and r['unit'] == "kHz"
        check(f"fromjson {spec}", ok, str(r)[:300])

        # Python export: exec and compare bit exact
        t = call("export python " + spec)['text']
        ns = {}
        exec(t, ns)
        ok = list(ns['b']) == ref['b'] and list(ns['a']) == ref['a']
        if ref['sos']:
            ok = ok and ns['sos'].tolist() == ref['sos']
        check(f"export python {spec}", ok, t[:300])
        # MATLAB: parse the bracket contents
        t = call("export matlab " + spec)['text']
        def mat(name):
            body = re.search(name + r" = \[(.*?)\];", t, re.S).group(1)
            rows = [r for r in body.replace("...", "").split(";")]
            return [[float(x) for x in r.split(",")] for r in rows]
        ok = mat("b")[0] == ref['b'] and mat("a")[0] == ref['a']
        if ref['sos']:
            ok = ok and mat("sos") == ref['sos']
        check(f"export matlab {spec}", ok, t[:300])
        # C header: compile-free check of the array initializers and sizes
        t = call("export c " + spec)['text']
        def carr_(name):
            body = re.search(r"lp_filter_" + name + r"\[[^=]*= \{(.*?)\};", t, re.S).group(1)
            return [float(x) for x in re.findall(r"[-+0-9.eE]+", body)]
        nb = int(re.search(r"#define LP_FILTER_NB (\d+)", t).group(1))
        ok = carr_("b") == ref['b'] and carr_("a") == ref['a'] and nb == len(ref['b'])
        if ref['sos']:
            ok = ok and carr_("sos") == [v for row in ref['sos'] for v in row]
        check(f"export c {spec}", ok, t[:300])
        # CSV
        t = call("export csv " + spec)['text']
        rows = [l.split(",") for l in t.strip().split("\n")]
        ok = rows[0] == ["b", "a"] and [float(r[0]) for r in rows[1:] if r[0]] == ref['b'] \
            and [float(r[1]) for r in rows[1:] if r[1]] == ref['a']
        check(f"export csv {spec}", ok, t[:300])
    # fixpoint settings (optional)
    doc = json.loads(call("tojson Hz rt=LP method=ellip")['json'])
    doc['fixpoint'] = {"simulate": True, "acc_auto": False, "QACC": {"WI": 3, "WF": 20, "quant": "ceil", "ovfl": "sat"}}
    fn = os.path.join(tmp, "fx.json")
    with open(fn, 'w', encoding='utf-8') as f:
        json.dump(doc, f)
    r = call(f"fromjson {fn}")
    check("fromjson fixpoint", r.get('has_fx') is True and r.get('fx_sim') is True
          and r.get('qacc') == [3, 20, "ceil", "sat"], str(r)[:300])
    doc['fixpoint']['QACC']['quant'] = "nearest"
    with open(fn, 'w', encoding='utf-8') as f:
        json.dump(doc, f)
    check("fromjson fixpoint invalid", 'error' in call(f"fromjson {fn}"))
    # invalid files give a readable error instead of a crash
    for i, text in enumerate(['{"format": "x"}', '{"format": "pyfda_cpp filter", "spec": {"rt": "XX"}}',
                              '{"format": "pyfda_cpp filter", "spec": {"N": 0}}', '[1, 2', '']):
        fn = os.path.join(tmp, f"bad{i}.json")
        with open(fn, 'w', encoding='utf-8') as f:
            f.write(text)
        r = call(f"fromjson {fn}")
        check(f"fromjson invalid {i}", 'error' in r, str(r))

# ---------------------------------------------------------------------------
# stimuli of the transient analysis: reference = pyfda's calc_stimulus_frame()
# with the bandlimited waveforms of pyfda.libs.special_functions (by Endolith)
from scipy.special import diric


def sawtooth_bl(t):
    y = np.zeros(t.shape)
    fs = 1 / (t[1] - t[0])
    for h in range(1, int(fs * np.pi) + 1):
        y += 2 / np.pi * -np.sin(h * t) / h
    return y


def triang_bl(t):
    y = np.zeros(t.shape)
    fs = 1 / (t[1] - t[0])
    for h in range(1, int(fs * np.pi) + 1, 2):
        y += 8 / np.pi**2 * -np.cos(h * t) / h**2
    return y


def rect_bl(t, duty=0.5):
    return sawtooth_bl(t - duty * 2 * np.pi) - sawtooth_bl(t) + 2 * duty - 1


def comb_bl(t):
    y = np.zeros(t.shape)
    fs = 1 / (t[1] - t[0])
    N = int(fs * np.pi) + 1
    for h in range(1, N):
        y += np.cos(h * t)
    return y / N


def stim_py(n_end, stim="dirac", a1=1., a2=0., f1=0.02, f2=0.03, phi1=0., phi2=0., t1=0., t2=0.,
            tw=10., bw1=0.5, bw2=0.5, n1=5, duty=0.5, bl=True, chirp="linear", noise="none", noi=0.1,
            mls_b=8, dc=0.):
    n = np.arange(n_end)
    x = np.zeros(n_end)
    t1_idx = int(np.round(t1))
    r1, r2 = phi1 / 180 * np.pi, phi2 / 180 * np.pi
    if stim == "dirac":
        if 0 <= t1_idx < n_end:
            x[t1_idx] = a1
    elif stim == "sinc":
        x = a1 * np.sinc(2 * (n - t1) * f1) + a2 * np.sinc(2 * (n - t2) * f2)
    elif stim == "gauss":
        x = a1 * sig.gausspulse(n - t1, fc=f1, bw=bw1) + a2 * sig.gausspulse(n - t2, fc=f2, bw=bw2)
    elif stim == "rect":
        n_rise = int(t1_idx - np.floor(tw / 2))
        x = a1 * np.where((n >= max(n_rise, 0)) & (n < min(n_rise + tw, n_end)), 1, 0)
    elif stim == "step":
        x[t1_idx:] = a1
    elif stim == "cos":
        x = a1 * np.cos(2 * np.pi * n * f1 + r1) + a2 * np.cos(2 * np.pi * n * f2 + r2)
    elif stim == "sine":
        x = a1 * np.sin(2 * np.pi * n * f1 + r1) + a2 * np.sin(2 * np.pi * n * f2 + r2)
    elif stim == "diric":
        x = a1 * diric(2 * np.pi * (n - t1) * f1, n1)
    elif stim == "chirp":  # phase in degrees (pyfda passes radians)
        x = a1 * sig.chirp(n, f1, n_end if t2 == 0 else t2, f2, method=chirp, phi=phi1)
    elif stim == "triang":
        x = a1 * (triang_bl(2 * np.pi * n * f1 + r1) if bl else sig.sawtooth(2 * np.pi * n * f1 + r1, width=0.5))
    elif stim == "saw":
        x = a1 * (sawtooth_bl(2 * np.pi * n * f1 + r1) if bl else sig.sawtooth(2 * np.pi * n * f1 + r1))
    elif stim == "rect_per":
        x = a1 * (rect_bl(2 * np.pi * n * f1 + r1, duty=duty) if bl
                  else sig.square(2 * np.pi * n * f1 + r1, duty=duty))
    elif stim == "comb":
        x = a1 * comb_bl(2 * np.pi * n * f1 + r1)
    elif stim == "am":
        x = a1 * np.sin(2 * np.pi * n * f1 + r1) * a2 * np.sin(2 * np.pi * n * f2 + r2)
    elif stim == "pmfm":
        x = a1 * np.sin(2 * np.pi * n * f1 + r1 + a2 * np.sin(2 * np.pi * n * f2 + r2))
    elif stim == "pwm":
        d = 1 / 2 + a2 / 2 * np.sin(2 * np.pi * n * f2 + r2)
        x = a1 * (rect_bl(2 * np.pi * n * f1 + r1, duty=d) if bl else sig.square(2 * np.pi * n * f1 + r1, duty=d))
    if noise == "mls":
        seed = [1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1,
                0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0][:mls_b]
        x = x + sig.max_len_seq(mls_b, length=n_end, state=seed)[0] * noi
    if stim != "step" and dc != 0:
        x = x + dc
    return np.asarray(x, dtype=float)


stim_cases = [
    (100, dict(stim="dirac", t1=3)), (100, dict(stim="dirac", t1=2.5, a1=2)),
    (200, dict(stim="sinc", t1=50, t2=80, f1=0.1, f2=0.05, a2=0.5)),
    (200, dict(stim="gauss", t1=60, t2=120, f1=0.1, f2=0.2, a2=0.7, bw1=0.3, bw2=0.8)),
    (100, dict(stim="rect", t1=20, tw=7)), (100, dict(stim="rect", t1=2, tw=10)),
    (100, dict(stim="step", t1=10, a1=3)),
    (300, dict(stim="cos", f1=0.013, f2=0.21, a2=0.3, phi1=30, phi2=-45)),
    (300, dict(stim="sine", f1=0.013, f2=0.21, a2=0.3, phi1=30, dc=0.5)),
    (300, dict(stim="diric", f1=0.05, t1=10, n1=7)), (300, dict(stim="diric", f1=0.1, n1=4)),
    (500, dict(stim="chirp", f1=0.01, f2=0.4)), (500, dict(stim="chirp", f1=0.01, f2=0.4, t2=300, phi1=20)),
    (500, dict(stim="chirp", chirp="quadratic", f1=0.01, f2=0.3)),
    (500, dict(stim="chirp", chirp="logarithmic", f1=0.01, f2=0.3)),
    (500, dict(stim="chirp", chirp="hyperbolic", f1=0.01, f2=0.3)),
    (400, dict(stim="triang", f1=0.023, phi1=10)), (400, dict(stim="triang", f1=0.023, bl=False)),
    (400, dict(stim="saw", f1=0.031)), (400, dict(stim="saw", f1=0.031, bl=False, phi1=90)),
    (400, dict(stim="rect_per", f1=0.02, duty=0.3)), (400, dict(stim="rect_per", f1=0.02, duty=0.3, bl=False)),
    (400, dict(stim="comb", f1=0.05)),
    (400, dict(stim="am", f1=0.1, f2=0.01, a2=1)), (400, dict(stim="pmfm", f1=0.1, f2=0.01, a2=2)),
    (400, dict(stim="pwm", f1=0.05, f2=0.005, a2=0.8)), (400, dict(stim="pwm", f1=0.05, f2=0.005, a2=0.8, bl=False)),
    (300, dict(stim="dirac", noise="mls", noi=0.5, mls_b=5)), (5000, dict(stim="sine", noise="mls", noi=1, mls_b=12)),
]
for n_end, kw in stim_cases:
    ref = stim_py(n_end, **kw)
    cmd = f"stim {n_end} " + " ".join(f"{k}={int(v) if isinstance(v, bool) else v}" for k, v in kw.items())
    r = call(cmd)
    if 'error' in r:
        check(cmd, False, r['error'])
        continue
    close(cmd, r['x'], ref, 1e-12, 1e-12)
# random noise can't be compared, check the statistics
for noise, mean, std in [("gauss", 0, 0.5), ("uniform", 0, 0.5 / np.sqrt(12)), ("randint", 1.5, np.sqrt(1.25))]:
    noi = 3 if noise == "randint" else 0.5
    x = np.array(call(f"stim 20000 stim=none noise={noise} noi={noi}")['x'])
    check(f"noise {noise}", abs(x.mean() - mean) < 0.05 and abs(x.std() - std) < 0.05 * std + 0.01,
          f"mean {x.mean()}, std {x.std()}")
# windowed FFT like pyfda's calc_fft (window / cgain, scaled by 1/N)
for w, par in [("rectangular", 0), ("hann", 0), ("kaiser", 8), ("flattop", 0)]:
    x = np.sin(2 * np.pi * 0.1 * np.arange(64)) + 0.3
    r = call(f"wfft {w} {par} {L(x)}")
    win = sig.get_window((w, par) if w == "kaiser" else ("boxcar" if w == "rectangular" else w), 64, fftbins=True)
    close(f"fft window {w}", r['win'], win, 1e-13)
    X = np.fft.fft(x * win / np.mean(win)) / 64
    close(f"windowed fft {w}", cplx(r['X']), X, 1e-12, 1e-14)
    close(f"ssb {w}", cplx(r['ssb']), np.insert(X[1:32] * 2, 0, X[0]), 1e-12, 1e-14)
    close(f"nenbw {w}", [r['nenbw']], [64 * np.sum(win**2) / np.sum(win)**2], 1e-12)

# ---------------------------------------------------------------------------
# fixpoint: reference = pyfda.libs.pyfda_fix_lib.Fixed.fixp() ('qfrac'), rewritten
# without the filterbroker; the C++ quantizer was also compared bit exact with
# pyfda itself (285 formats) and with fir_df_pyfixp
def fixp_py(y, WI, WF, quant, ovfl):
    y = np.asarray(y, dtype=float) * 2.0**WF
    yq = {'floor': np.floor, 'round': np.round, 'fix': np.trunc, 'ceil': np.ceil,
          'none': lambda v: v}[quant](y)
    n_over = 0
    if ovfl != 'none':
        MSB = 2.0**(WI + WF - 1)
        MAX, MIN = 2 * MSB - 1, -2 * MSB
        over_neg, over_pos = yq < MIN, yq > MAX
        n_over = int(np.sum(over_neg) + np.sum(over_pos))
        if ovfl == 'sat':
            yq = np.where(over_pos, MAX, np.where(over_neg, MIN, yq))
        else:
            yq = np.where(over_pos | over_neg,
                          yq - 4. * MSB * np.trunc((np.sign(yq) * 2 * MSB + yq) / (4 * MSB)), yq)
    return yq / 2.0**WF, n_over


rng = np.random.default_rng(1)
for WI in [0, 1, 3, 7]:
    for WF in [0, 3, 8, 15, 30]:
        if WI + WF == 0:
            continue
        for quant in ['floor', 'round', 'fix', 'ceil', 'none']:
            for ovfl in ['wrap', 'sat', 'none']:
                x = rng.normal(0, 2**WI * 1.5, 100)
                x[:4] = [0.5 / 2**WF, 1.5 / 2**WF, -2.5 / 2**WF, 2**WI]
                ref, n_over = fixp_py(x, WI, WF, quant, ovfl)
                r = call(f"fixp {WI} {WF} {quant} {ovfl} {L(x)}")
                check(f"fixp {WI}.{WF} {quant} {ovfl}", r['y'] == list(ref) and r['n_over'] == n_over,
                      f"{r['n_over']} vs {n_over}")
for v, W in [(-5, 8), (1234, 16), (0, 4), (-1, 12), (127, 8), (-128, 8), (5, 3)]:
    r = call(f"fxbase {v} {W}")
    u = v & ((1 << W) - 1)
    csd_val = sum({'+': 1, '-': -1, '0': 0}[c] * 2**i for i, c in enumerate(reversed(r['csd'])))
    nonadj = all(not (a != '0' and b != '0') for a, b in zip(r['csd'], r['csd'][1:]))
    check(f"fxbase {v} {W}", r['bin'] == format(u, f'0{W}b') and int(r['hex'], 16) == u
          and int(r['oct'], 8) == u and csd_val == v and nonadj, str(r))


def fx_fir_py(b, x, qi, qcb, qacc, qo):
    bq, _ = fixp_py(b, *qcb)
    xq, _ = fixp_py(x, *qi)
    y = np.zeros(len(x))
    for k in range(len(x)):
        prods = [fixp_py(xq[k - i] * bq[i], *qacc)[0] for i in range(min(len(bq), k + 1))]
        y[k] = fixp_py(fixp_py(np.sum(prods), *qacc)[0], *qo)[0]
    return y


def fx_sos_py(sos, x, qi, qcb, qca, qacc, qo):
    s, _ = fixp_py(x, *qi)
    for sec in sos:
        b, _ = fixp_py(sec[:3], *qcb)
        a, _ = fixp_py(sec[4:], *qca)
        y = np.zeros(len(s))
        for n in range(len(s)):
            xs = [s[n - i] if n - i >= 0 else 0 for i in range(3)]
            ys = [y[n - i] if n - i >= 0 else 0 for i in (1, 2)]
            acc_b = sum(fixp_py(b[i] * xs[i], *qacc)[0] for i in range(3))
            acc_a = sum(fixp_py(a[i] * ys[i], *qacc)[0] for i in range(2))
            y[n] = fixp_py(fixp_py(acc_b - acc_a, *qacc)[0], *qo)[0]
        s = y
    return s


fq = lambda t: ",".join(map(str, t))
fx_cases = [((0, 15, 'round', 'sat'), (0, 15, 'round', 'sat'), (1, 14, 'round', 'sat'), (2, 30, 'floor', 'wrap'),
             (0, 15, 'round', 'sat')),
            ((0, 7, 'round', 'sat'), (0, 7, 'floor', 'wrap'), (1, 6, 'round', 'sat'), (1, 9, 'round', 'wrap'),
             (0, 5, 'fix', 'sat')),
            ((2, 3, 'round', 'sat'), (1, 5, 'round', 'sat'), (1, 5, 'floor', 'sat'), (3, 6, 'floor', 'sat'),
             (1, 4, 'round', 'wrap'))]
x = rng.normal(0, 0.4, 200)
b_fir = sig.firwin(21, 0.3) * np.linspace(0.8, 1.2, 21)  # asymmetric to check the order of the taps
sos_iir = sig.ellip(4, 1, 40, 0.2, output='sos')
for qi, qcb, qca, qacc, qo in fx_cases:
    keys = f"qi={fq(qi)} qcb={fq(qcb)} qacc={fq(qacc)} qo={fq(qo)}"
    r = call(f"fxfir {L(b_fir)} {L(x)} {keys}")
    check(f"fxfir {keys}", r['y'] == list(fx_fir_py(b_fir, x, qi, qcb, qacc, qo)), str(r)[:200])
    r = call(f"fxsos {L(sos_iir)} {L(x)} {keys} qca={fq(qca)}")
    check(f"fxsos {keys}", r['y'] == list(fx_sos_py(sos_iir, x, qi, qcb, qca, qacc, qo)), str(r)[:200])
# with enough bits the fixpoint response is close to the floating point response
r = call(f"fxfir {L(b_fir)} {L(x)} qi=3,31,round,sat qcb=0,31,round,sat qacc=6,62,floor,wrap qo=4,31,round,sat")
close("fxfir 32 bit vs lfilter", r['y'], sig.lfilter(b_fir, 1, x), 0, 1e-8)
r = call(f"fxsos {L(sos_iir)} {L(x)} qi=3,31,round,sat qcb=1,30,round,sat qca=1,30,round,sat "
         "qacc=8,61,floor,wrap qo=4,31,round,sat")
close("fxsos 32 bit vs sosfilt", r['y'], sig.sosfilt(sos_iir, x), 0, 1e-7)
r = call(f"fxauto fir {L(b_fir)} qi=0,15,round,sat qcb=0,15,round,sat qacc=0,0,floor,wrap qo=0,15,round,sat")
check("fxauto fir", r['qcb'] == [0, 15] and r['qacc'] == [0 + 0 + int(np.ceil(np.log2(np.sum(np.abs(b_fir))))), 30],
      str(r))
r = call(f"fxauto sos {L(sos_iir)} qi=0,15,round,sat qcb=0,15,round,sat qca=0,14,round,sat qacc=0,0,floor,wrap "
         "qo=0,15,round,sat")
check("fxauto sos", r["qcb"] == [1, 15] and r["qca"] == [1, 14] and r["qacc"] == [1 + 3, 30], str(r))

# ---- moving average, delay and manual filters (pyfda filter_widgets/ma.py, delay.py, manual.py) ----
def ceil_odd(x):
    """ smallest odd integer >= x (as documented in pyfda, its round_odd(x + 1) gives x + 2 for odd x) """
    c = int(np.ceil(x))
    return c if c % 2 else c + 1


def root_dist(a, b):
    """ max. distance between two sets of roots (greedy matching) """
    a, e = list(a), 0.
    if len(a) != len(b):
        return np.inf
    for x in b:
        i = int(np.argmin([abs(x - y) for y in a]))
        e = max(e, abs(x - a.pop(i)))
    return e


def calc_ma_py(rt, delays, stages, norm_on, fo='man', f_sb=0.2, A_SB=40.):
    """ copy of pyfda's MA.calc_ma() and its lp_min / hp_min / bp / bs order calculation """
    a_sb = 10 ** (-A_SB / 20)
    if rt == 'LP' and fo == 'min':
        delays = int(np.ceil(1 / (a_sb ** (1 / stages) * np.sin(f_sb * np.pi))))
    elif rt == 'HP' and fo == 'min':
        delays = int(np.ceil(1 / (a_sb ** (1 / stages) * np.sin((0.5 - f_sb) * np.pi))))
    elif rt in ('BP', 'BS'):
        delays = ceil_odd(delays)
    k = 1.
    l_taps = delays + 1
    norm = l_taps
    idx = np.arange(1, l_taps)
    b0 = np.ones(l_taps)
    if rt == 'HP':
        b0[::2] = -1.
        idx = np.arange(l_taps)
        if l_taps % 2 == 0:
            idx = np.delete(idx, round(l_taps / 2.))
        else:
            idx = np.delete(idx, int(l_taps / 2.)) + 0.5
    elif rt == 'BP':
        b0[1::2] = 0
        b0[::4] = -1
        l_taps = l_taps + 1
        idx = np.delete(np.arange(l_taps), [0, l_taps // 2]) + l_taps / 4
        norm = np.sum(abs(b0))
    elif rt == 'BS':
        b0[1::2] = 0
        l_taps = l_taps + 1
        idx = np.delete(np.arange(l_taps), [0, l_taps // 2])
        norm = np.sum(b0)
    z0 = np.exp(-2.j * np.pi * idx / l_taps)
    b = 1
    for _ in range(stages):
        b = np.convolve(b0, b)
    z = np.repeat(z0, stages)
    if norm_on:
        b = b / (norm ** stages)
        k = 1. / norm ** stages
    return delays, b, z, k


for rt in ('LP', 'HP', 'BP', 'BS'):
    for delays in (1, 2, 3, 4, 7, 12, 13):
        for stages in (1, 2, 3):
            for nrm in (0, 1):
                name = f"ma {rt} M={delays} stages={stages} norm={nrm}"
                d_ref, b_ref, z_ref, k_ref = calc_ma_py(rt, delays, stages, nrm)
                r = call(f"design rt={rt} method=ma fo=manual N={delays} stages={stages} norm={nrm} f_s=1")
                check(name + " N", r['N'] == d_ref, str(r['N']))
                close(name + " b", r['b'], b_ref, 1e-12)
                check(name + " a", r['a'] == [1], str(r['a']))
                if rt in ('LP', 'HP'):
                    close(name + " z", cplx(r['z']), z_ref, 1e-12)
                    close(name + " k", r['k'], k_ref, 1e-14)
                # the zeros are the roots of b (pyfda's BP / BS zeros are not, the port uses the roots)
                if stages == 1:
                    check(name + " z vs roots", root_dist(cplx(r['z']), np.roots(b_ref)) < 1e-6)
for rt in ('LP', 'HP'):
    for f_sb, A_SB in ((0.05, 20), (0.1, 40), (0.3, 30), (0.2, 60)):
        for stages in (1, 2, 4):
            d_ref, b_ref, z_ref, k_ref = calc_ma_py(rt, 0, stages, 1, 'min', f_sb, A_SB)
            r = call(f"design rt={rt} method=ma fo=min stages={stages} f_s=2 f_sb={2 * f_sb!r} A_SB={A_SB}")
            check(f"ma min {rt} f_sb={f_sb} A_SB={A_SB} stages={stages}", r['N'] == d_ref and
                  np.allclose(r['b'], b_ref, rtol=1e-12, atol=0), f"{r['N']} != {d_ref}")
for N in (1, 5, 20):
    r = call(f"design method=delay N={N} f_s=1")
    check(f"delay N={N}", r['b'] == [0] * N + [1] and r['a'] == [1] and r['fir'] and len(r['p']) == N
          and len(r['z']) == 0, str(r)[:200])
# manual coefficients: normalized to a[0] = 1, FIR / IIR detected from a
b_m, a_m = [2., -1., 0.5], [2., -0.5, 0.25]
r = call(f"design method=manual b={L(b_m)} a={L(a_m)} f_s=1")
close("manual iir b", r['b'], np.array(b_m) / 2, 1e-15)
close("manual iir a", r['a'], np.array(a_m) / 2, 1e-15)
check("manual iir method", r['method'] == 'manual_iir' and not r['fir'] and len(r['sos']) == 1, str(r))
z_ref, p_ref, k_ref = sig.tf2zpk(b_m, a_m)
check("manual iir z", root_dist(cplx(r['z']), z_ref) < 1e-12 and root_dist(cplx(r['p']), p_ref) < 1e-12)
r = call(f"design method=manual b=1,2,3,2,1 a=1 f_s=1")
check("manual fir", r['method'] == 'manual_fir' and r['fir'] and r['N'] == 4 and r['b'] == [1, 2, 3, 2, 1], str(r))
# manual poles / zeros
zm = [0.5 + 0.5j, 0.5 - 0.5j, -1]
pm = [0.8 * np.exp(0.3j), 0.8 * np.exp(-0.3j)]
r = call("design method=manual z=" + ",".join(f"{float(v.real)!r}:{float(v.imag)!r}" for v in zm) + " p=" +
         ",".join(f"{float(v.real)!r}:{float(v.imag)!r}" for v in pm) + " k=0.5 f_s=1")
b_ref, a_ref = sig.zpk2tf(zm, pm + [0], 0.5)  # one pole added at the origin to make it causal
close("manual zpk b", r['b'], b_ref, 1e-12)
close("manual zpk a", r['a'], a_ref, 1e-12)
check("manual zpk iir", r['method'] == 'manual_iir' and len(r['p']) == 3, str(r)[:200])
r = call("design method=manual z=0.5:0,-1:0 p=0:0,0:0 k=2 f_s=1")
check("manual zpk fir", r['method'] == 'manual_fir' and np.allclose(r['b'], [2, 1, -1]), str(r)[:200])
# JSON round trip of the new methods
for spec in ("rt=BP method=ma fo=manual N=8 stages=2 norm=0", "method=delay N=7",
             "method=manual b=1,-0.5,0.25 a=1,0.3", "method=manual z=0.5:0.5,0.5:-0.5 p=0.9:0 k=3"):
    r = call(f"design {spec} f_s=1")
    js = call(f"tojson f_S {spec} f_s=1")['json']
    with open("_rt.json", "w", encoding="utf-8") as f:
        f.write(js)
    r2 = call("fromjson _rt.json")
    check(f"json round trip {spec}", r2['b'] == r['b'] and r2['a'] == r['a'] and r2['method'] == r['method'],
          f"{r2} != {r}"[:300])
os.remove("_rt.json")

# ---- formula stimulus (numexpr syntax like pyfda's "Formula" stimulus) ----
try:
    import numexpr
except ImportError:
    numexpr = None
if os.environ.get("NO_NUMEXPR"):
    numexpr = None
formulas = [
    "A1 * abs(sin(2 * pi * f1 * n))",
    "A1 * sin(2*pi*f1*n + phi1/180*pi) + A2 * cos(2*pi*f2*n)",
    "where(n < T1, 0, A1) - 0.5 * (n % 7 > 3)",
    "exp(-n / 20.) * cos(2 * pi * f1 * n) ** 2",
    "-2 ** 2 + n * 0 + sqrt(abs(sin(n))) * sign(cos(n))",
    "arctan2(sin(n), cos(n)) / pi + log1p(n) - log10(n + 1) + tanh(n / 10 - 3)",
    "(n > 5) & (n < 20) | (n == 40)",
    "t * f_S - n + minimum(n, 10) + maximum(n, 30) + floor(n / 3) + ceil(n / 4)",
    "e ** (-(n - T2) ** 2 / (2 * BW1 ** 2)) * N1",
    "2.5e-1 * n - 1E1 + .5",
]
pars = dict(A1=0.7, A2=0.3, f1=0.03, f2=0.11, phi1=30.0, phi2=0.0, T1=12.0, T2=25.0, N1=5, BW1=3.0, BW2=0.5)
for fo in formulas:
    r = call(f"stim 64 a1={pars['A1']} a2={pars['A2']} f1={pars['f1']} f2={pars['f2']} phi1={pars['phi1']} "
             f"t1={pars['T1']} t2={pars['T2']} n1={pars['N1']} bw1={pars['BW1']} fs=8 formula={fo}")
    n = np.arange(64, dtype=float)
    ld = dict(pars, n=n, t=n / 8, f_S=8.0, pi=np.pi, e=np.e)
    if numexpr is not None:
        ref = numexpr.evaluate(fo, local_dict=ld).astype(float)
    else:
        ref = eval(fo, {k: getattr(np, k) for k in dir(np) if not k.startswith('_')}, ld)
    close(f"formula {fo}", r['x'], np.broadcast_to(ref, n.shape), 1e-12, 1e-14)
for bad in ("sin(n", "foo(n)", "n +* 2", "q * 2", "3 * j", ""):
    r = call(f"stim 8 formula={bad}")
    check(f"formula error '{bad}'", 'error' in r, str(r)[:100])

# ---- spectrogram (scipy.signal.spectrogram, detrend='constant', one-sided) ----
xs = np.sin(2 * np.pi * 0.05 * np.arange(700) ** 1.3 / 30) + 0.2 * rng.normal(size=700) + 0.3
for mode in ("psd", "magnitude", "angle"):
    for dens in (1, 0):
        for wname, par, nper, novl in (("hann", 0, 64, 32), ("kaiser", 6.0, 101, 50), ("rectangular", 0, 128, 0)):
            if mode != "psd" and dens == 0:
                continue
            win = sig.get_window(('kaiser', par) if wname == 'kaiser' else ('boxcar' if wname == 'rectangular' else wname),
                                 nper, fftbins=True)
            f, t, S = sig.spectrogram(xs, 2.0, window=win, nperseg=nper, noverlap=novl, detrend='constant',
                                      scaling='density' if dens else 'spectrum', mode=mode)
            r = call(f"spgr {mode} {dens} 2 {wname} {par} {nper} {novl} {L(xs)}")
            name = f"spectrogram {mode} dens={dens} {wname} {nper}/{novl}"
            close(name + " f", r['f'], f, 1e-12)
            close(name + " t", r['t'], t, 1e-12)
            got = np.array(r['s']).T
            if mode == "angle":  # compare angles where the magnitude is not tiny
                _, _, M = sig.spectrogram(xs, 2.0, window=win, nperseg=nper, noverlap=novl, mode='magnitude')
                m = M > 1e-6 * M.max()
                check(name, got.shape == S.shape and np.allclose(np.exp(1j * got[m]), np.exp(1j * S[m]), atol=1e-8))
            else:
                close(name, got, S, 1e-9)

# ---- amplitude units (pyfda special_functions.unit2lin / lin2unit) ----
def lin2unit(lin, fir, pb, unit):
    if unit == 'dB':
        if pb:
            return -20 * np.log10(1 - lin) if not fir else 20 * np.log10((1 + lin) / (1 - lin))
        return -20 * np.log10(lin)
    return lin * lin if unit == 'W' else lin


def unit2lin(v, fir, pb, unit):
    if unit == 'dB':
        if pb:
            return 1 - 10 ** (-v / 20) if not fir else (10 ** (v / 20) - 1) / (10 ** (v / 20) + 1)
        return 10 ** (-v / 20)
    return np.sqrt(v) if unit == 'W' else v


for fir in (0, 1):
    for pb in (0, 1):
        for db in (0.1, 1, 3, 40, 80):
            for unit in ('V', 'W'):
                ref = lin2unit(unit2lin(db, fir, pb, 'dB'), fir, pb, unit)
                r = call(f"amp from {db} {unit} {fir} {pb}")
                close(f"amp from dB {db} {unit} fir={fir} pb={pb}", r['v'], ref, 1e-12)
                r = call(f"amp to {ref!r} {unit} {fir} {pb}")
                close(f"amp to dB {unit} fir={fir} pb={pb}", r['v'], db, 1e-10)

proc.stdin.close()
proc.wait()
print(f"{n_pass} checks passed, {n_fail} failed")
sys.exit(1 if n_fail else 0)
