import sys, os, json, subprocess
sys.path.insert(0, os.path.dirname(__file__))
from designs import D
exe, outroot = sys.argv[1], sys.argv[2]
for name, d in D:
    out = os.path.join(outroot, name); os.makedirs(out, exist_ok=True)
    spec = {"rt": d['rt'], "method": d['fc'].lower(), "order": d['fo'], "N": d['N'], "unit": "f_S", "f_S": 1,
            "A_PB": d['A_PB'], "A_SB": d['A_SB'], "W_PB": 1, "W_SB": 1, "window": "Kaiser", "win_par": 10,
            "order_alg": "ichige", "grid_density": 16, "ma_stages": 1, "ma_norm": True}
    for k in ('f_pb','f_pb2','f_sb','f_sb2','f_c','f_c2'):
        if k in d: spec[k.upper().replace('F_','F_')] = d[k]
    spec = {k: v for k, v in spec.items()}
    fn = os.path.join(out, 'filter.json')
    json.dump({"format": "pyfda_cpp filter", "version": 1, "spec": spec}, open(fn, 'w'), indent=1)
    subprocess.run([exe, '--config-dir', os.path.join(outroot, 'cfg'), '--load-filter', fn, '--save-filter', os.path.join(out, 'saved.json')], env=dict(os.environ, QT_QPA_PLATFORM='offscreen'), capture_output=True, timeout=300)
    args = [exe, '--config-dir', os.path.join(outroot, 'cfg'), '--load-filter', fn, '--screenshot', out]
    if name == 'd1_ellip_lp': args += ['--data', os.environ['DATA_CSV'], '--filter']
    r = subprocess.run(args, env=dict(os.environ, QT_QPA_PLATFORM='offscreen'), capture_output=True, text=True, timeout=300)
    print(name, r.returncode, r.stderr.strip().splitlines()[-1:] if r.returncode else '')
