import sys, os, json, traceback
sys.path.insert(0, os.path.dirname(__file__))
from designs import D, a_pb_lin, a_sb_lin
outroot = sys.argv[1]
from pyfda import pyfdax as px
from pyfda import filterbroker as fb
from PyQt5.QtWidgets import QApplication, QTabWidget, QSplitter
from PyQt5.QtCore import QTimer
app = QApplication(sys.argv[:1]); app.setStyle('Fusion')
if len(px.QSS.QSS_RC) > 20: app.setStyleSheet(px.QSS.QSS_RC)
px.rc.params['screen'] = {'ref_dpi':96,'ldpi':96,'pdpi':96,'scaling':1,'height':1000,'width':1600}
px.FilterTreeBuilder().build_fil_tree()
w = px.Pyfda(); w.resize(1600, 1000); w.show()
log = open(os.path.join(outroot, 'pyfda_log.txt'), 'w')
def pe(n=30):
    for _ in range(n): app.processEvents()
def find(cls_name):
    from PyQt5.QtWidgets import QWidget
    return [c for c in w.findChildren(QWidget) if type(c).__name__ == cls_name]
def run():
    try:
        tabs = w.findChildren(QTabWidget)
        plot_tabs = [t for t in tabs if t.count() and t.tabText(0) == '|H(f)|'][0]
        in_tabs = [t for t in tabs if t.count() and t.tabText(0) == 'Specs'][0]
        spl = w.findChildren(QSplitter)
        for s in spl:
            if s.count() == 2 and s.orientation() == 1: s.setSizes([480, 1120])
        specs = find('InputSpecs')[0]
        dfw = find('PlotDataFilt')
        for name, d in D:
            out = os.path.join(outroot, name); os.makedirs(out, exist_ok=True)
            f = fb.fil[0]
            fir = d['ft'] == 'FIR'
            f.update({'ft': d['ft'], 'fc': d['fc'], 'rt': d['rt'].lower(), 'fo': d['fo'], 'N': d['N'],
                      'amp_specs_unit': 'dB', 'freq_specs_unit': 'f_S', 'f_s': 1.0})
            for k in ('f_pb','f_pb2','f_sb','f_sb2','f_c','f_c2'):
                if k in d: f[k] = d[k]
            f['a_pb'] = a_pb_lin(d['A_PB'], fir); f['a_sb'] = a_sb_lin(d['A_SB'])
            specs.emit({'data_changed': 'filter_loaded'}); pe()
            specs.start_design_filt(); pe()
            log.write(f"{name}: N={fb.fil[0]['N']} info={fb.fil[0].get('info')}\n")
            import numpy as np
            ba = fb.fil[0]['ba']
            json.dump({'b': np.real(np.asarray(ba[0])).tolist(), 'a': np.real(np.asarray(ba[1])).tolist(), 'N': int(fb.fil[0]['N'])},
                      open(os.path.join(out, 'ba.json'), 'w'))
            if dfw and name == 'd1_ellip_lp':
                dfw[0].load_file(os.environ['DATA_CSV'], 'csv'); pe(); dfw[0].filter_data(); pe()
            for i in range(plot_tabs.count()):
                plot_tabs.setCurrentIndex(i); pe()
                w.grab().save(os.path.join(out, f'plot{i}.png'))
            plot_tabs.setCurrentIndex(0)
            for i in range(in_tabs.count()):
                in_tabs.setCurrentIndex(i); pe()
                w.grab().save(os.path.join(out, f'input{i}.png'))
            in_tabs.setCurrentIndex(0)
    except Exception:
        log.write(traceback.format_exc())
    log.close(); app.quit()
QTimer.singleShot(1500, run)
app.exec_()
