#include "tran_tab.hpp"
#include "settings.hpp"

#include "conversions.hpp"
#include "data_io.hpp"
#include "expr.hpp"
#include "filtering.hpp"
#include "logger.hpp"
#include "plot_widget.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QTabWidget>
#include <QTextStream>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pyfda;

namespace {
const QStringList F_KEYS = {"f1", "f2"};           // frequencies, displayed in the frequency unit
const QStringList T_KEYS = {"t1", "t2", "tw"};     // times, displayed in the time unit
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

double db(double v) { return 20 * std::log10(std::max(std::fabs(v), 1e-15)); }

QString fmt(double v) { return QString::number(v, 'g', 8); }

bool parse(const QLineEdit *e, double &v) {
    QString t = e->text().trimmed();
    if (t.contains(',') && !t.contains('.')) t.replace(',', '.');
    bool ok = false;
    const double d = t.toDouble(&ok);
    if (!ok || !std::isfinite(d)) return false;
    v = d;
    return true;
}

// complex values like in pyfda, e.g. "1 - 3j", "2j", "0.5" or "exp(1j*pi/4)"
bool parseCplx(const QLineEdit *e, double &re, double &im) {
    QString t = e->text().trimmed();
    if (t.contains(',') && !t.contains('.') && !t.contains('(')) t.replace(',', '.');
    try {
        const Expr ex(t.toStdString(), {"pi", "e"});
        const cplx v = ex.evalc({PI, std::exp(1.0)});
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
        re = v.real();
        im = v.imag();
        return true;
    } catch (const std::exception &) {
        return false;
    }
}

QString fmtCplx(double re, double im) {
    if (im == 0) return fmt(re);
    if (re == 0) return fmt(im) + "j";
    return fmt(re) + (im < 0 ? " - " : " + ") + fmt(std::fabs(im)) + "j";
}
}  // namespace

TranView::TranView(QWidget *parent) : DesignView(parent) {
    auto *lay = new QVBoxLayout(this);
    lay->setContentsMargins(4, 4, 4, 4);

    // --- row 1: stimulus selection
    auto *r1 = new QHBoxLayout();
    m_stim = new QComboBox(this);
    for (const auto &s : stim_list()) m_stim->addItem(s.name, int(s.stim));
    m_stim->setCurrentIndex(m_stim->findData(int(Stim::Dirac)));
    m_stim->setToolTip(tr("Stimulus x[n]"));
    m_chirp = new QComboBox(this);
    m_chirp->addItems({tr("linear"), tr("quadratic"), tr("logarithmic"), tr("hyperbolic")});
    m_chirp->setToolTip(tr("Frequency sweep from f1 at n = 0 to f2 at T2 (T2 = 0: end of the interval)"));
    m_bl = new QCheckBox("BL", this);
    m_bl->setChecked(true);
    m_bl->setToolTip(tr("<span>Bandlimit the signal to the Nyquist frequency to avoid aliasing (slower "
                        "for many data points).</span>"));
    m_step_err = new QCheckBox(tr("Error"), this);
    m_step_err->setToolTip(tr("Show the settling error, i.e. the step response minus its final value |H(0)|"));
    r1->addWidget(new QLabel(tr("Stimulus:"), this));
    r1->addWidget(m_stim);
    r1->addWidget(m_chirp);
    r1->addWidget(m_bl);
    r1->addWidget(m_step_err);
    r1->addSpacing(16);
    m_noise = new QComboBox(this);
    m_noise->addItems({tr("No noise"), tr("Gauss"), tr("Uniform"), tr("Random int."), tr("MLS"), tr("Brownian")});
    m_noise->setToolTip(tr("<span>Additive noise: Gaussian (amplitude = standard deviation), uniform (amplitude = "
                           "peak-to-peak), random integers 0 ... amplitude, maximum length sequence (0 / amplitude) "
                           "or Brownian (integrated Gaussian noise)</span>"));
    m_noi = new QLineEdit(fmt(m_p.noi), this);
    m_noi->setMaximumWidth(70);
    m_noi->setToolTip(tr("Noise amplitude; a complex value like 0.1 + 0.2j adds independent noise to the "
                         "real and the imaginary part"));
    m_lbl_mls = new QLabel(tr("bits:"), this);
    m_mls_b = new QSpinBox(this);
    m_mls_b->setRange(2, 32);
    m_mls_b->setValue(m_p.mls_b);
    m_mls_b->setToolTip(tr("Number of bits of the maximum length sequence (period 2^bits - 1)"));
    r1->addWidget(m_noise);
    r1->addWidget(m_noi);
    r1->addWidget(m_lbl_mls);
    r1->addWidget(m_mls_b);
    r1->addStretch(1);
    m_lbl_cmplx = new QLabel(tr("<b>Complex signals</b>"), this);
    m_lbl_cmplx->setStyleSheet("QLabel {color: darkorange}");
    m_lbl_cmplx->setToolTip(tr("<span>The stimulus is complex: real and imaginary parts are shown in two plots "
                               "and the spectra are two-sided (-f_S/2 ... f_S/2).</span>"));
    m_lbl_cmplx->setVisible(false);
    r1->addWidget(m_lbl_cmplx);
    lay->addLayout(r1);

    // --- rows 2 and 3: parameters, only the ones used by the stimulus are shown
    auto *r2 = new QHBoxLayout();
    auto *r3 = new QHBoxLayout();
    auto add = [&](QHBoxLayout *row, const QString &key, const QString &label, const QString &tip) {
        row->addWidget(m_param_widgets[key] = new QWidget(this));
        auto *h = new QHBoxLayout(m_param_widgets[key]);
        h->setContentsMargins(0, 0, 6, 0);
        auto *l = new QLabel(label, this);
        auto *e = new QLineEdit(this);
        e->setMaximumWidth(80);
        e->setToolTip(tip);
        l->setToolTip(tip);
        h->addWidget(l);
        h->addWidget(e);
        m_edits[key] = e;
        connect(e, &QLineEdit::editingFinished, this, [this] { readEdits(); });
    };
    add(r2, "a1", "A1:", tr("Amplitude of the first signal component; complex values like 1 - 3j are allowed"));
    add(r2, "a2", "A2:", tr("Amplitude of the second component (modulation index for PM / FM, "
                            "modulation depth for PWM); complex values like 1 - 3j are allowed"));
    add(r2, "phi1", QString::fromUtf8("φ1 / °:"), tr("Phase of the first component in degrees"));
    add(r2, "phi2", QString::fromUtf8("φ2 / °:"), tr("Phase of the second component in degrees"));
    add(r2, "duty", QString::fromUtf8("α:"), tr("Duty cycle, 0 ... 1"));
    add(r2, "n1", "N1:", tr("Order of the periodic sinc (Dirichlet) function"));
    add(r2, "dc", "DC:", tr("DC offset added to the stimulus, may be complex"));
    r2->addStretch(1);
    add(r3, "f1", "f1:", tr("Frequency of the first component (center frequency for impulses)"));
    add(r3, "f2", "f2:", tr("Frequency of the second component"));
    add(r3, "t1", "T1:", tr("Position of the (first) impulse or step"));
    add(r3, "t2", "T2:", tr("Position of the second impulse, end of the chirp sweep (0 = end of the interval)"));
    add(r3, "tw", "T_W:", tr("Width of the rect impulse"));
    add(r3, "bw1", "BW1:", tr("Relative bandwidth of the first Gaussian pulse"));
    add(r3, "bw2", "BW2:", tr("Relative bandwidth of the second Gaussian pulse"));
    r3->addStretch(1);
    lay->addLayout(r2);
    lay->addLayout(r3);

    // formula (pyfda: numexpr syntax) and file stimulus
    m_wdg_formula = new QWidget(this);
    auto *hf = new QHBoxLayout(m_wdg_formula);
    hf->setContentsMargins(0, 0, 0, 0);
    m_formula = new QLineEdit(QString::fromStdString(m_p.formula), this);
    m_formula->setToolTip(tr("<span>Formula in numexpr syntax with the index <i>n</i>, the time <i>t</i> = n / f_S, "
                             "the parameters A1, A2, f1, f2, phi1, phi2 (degrees), T1, T2, N1, BW1, BW2 (normalized "
                             "frequencies, times in samples), f_S, pi, e and the functions sin, cos, tan, arcsin, "
                             "arccos, arctan, arctan2, sinh, cosh, tanh, exp, log, log10, sqrt, abs, sign, floor, "
                             "ceil, round, where, minimum, maximum, real, imag, conj, complex. Operators: + - * / % **, "
                             "comparisons, &amp; | ~. Complex values with j, e.g. exp(2j * pi * f1 * n).</span>"));
    hf->addWidget(new QLabel("x[n] =", this));
    hf->addWidget(m_formula, 1);
    lay->addWidget(m_wdg_formula);
    m_wdg_file = new QWidget(this);
    auto *hl = new QHBoxLayout(m_wdg_file);
    hl->setContentsMargins(0, 0, 0, 0);
    m_load = new QPushButton(tr("Load file ..."), this);
    m_load->setToolTip(tr("Load the stimulus from a CSV / text, wav or npy file (first data column)"));
    m_file_lbl = new QLabel(tr("no file loaded"), this);
    m_file_norm = new QCheckBox(tr("Normalize"), this);
    m_file_norm->setChecked(true);
    m_file_norm->setToolTip(tr("Normalize the file data to max |x| = 1 (A1 scales the data)"));
    hl->addWidget(m_load);
    hl->addWidget(m_file_lbl, 1);
    hl->addWidget(m_file_norm);
    lay->addWidget(m_wdg_file);

    // --- row 4: number of points, export
    auto *r4 = new QHBoxLayout();
    m_N = new QSpinBox(this);
    m_N->setRange(0, 10000000);
    m_N->setValue(0);
    m_N->setSpecialValueText(tr("auto"));
    m_N->setToolTip(tr("<span>Number of displayed data points N (0 = auto: length of the impulse response, "
                       "at least 100)</span>"));
    m_N_start = new QSpinBox(this);
    m_N_start->setRange(0, 10000000);
    m_N_start->setToolTip(tr("Index of the first displayed data point; the spectra are calculated from "
                             "N_start ... N_start + N - 1"));
    auto *exp = new QPushButton(tr("Export CSV ..."), this);
    exp->setToolTip(tr("Save n, t, x[n] and y[n] as CSV file"));
    r4->addWidget(new QLabel("N =", this));
    r4->addWidget(m_N);
    r4->addWidget(new QLabel("N_start =", this));
    r4->addWidget(m_N_start);
    r4->addStretch(1);
    r4->addWidget(exp);
    lay->addLayout(r4);

    // --- time / frequency tabs
    m_tabs = new QTabWidget(this);
    auto *wt = new QWidget(this);
    auto *vt = new QVBoxLayout(wt);
    vt->setContentsMargins(2, 2, 2, 2);
    auto *ct = new QHBoxLayout();
    m_t_stim = new QCheckBox(tr("Stimulus x[n]"), this);
    m_t_resp = new QCheckBox(tr("Response y[n]"), this);
    m_t_db = new QCheckBox("dB", this);
    m_t_stim->setChecked(true);
    m_t_resp->setChecked(true);
    m_t_db->setToolTip(tr("Show the magnitude in dB"));
    ct->addWidget(m_t_stim);
    ct->addWidget(m_t_resp);
    ct->addWidget(m_t_db);
    ct->addStretch(1);
    vt->addLayout(ct);
    m_plot_t = new PlotWidget(this);
    vt->addWidget(m_plot_t, 1);
    m_plot_ti = new PlotWidget(this);  // imaginary parts of complex signals
    m_plot_ti->setVisible(false);
    vt->addWidget(m_plot_ti, 1);
    m_tabs->addTab(wt, tr("Time"));

    auto *wf = new QWidget(this);
    auto *vf = new QVBoxLayout(wf);
    vf->setContentsMargins(2, 2, 2, 2);
    auto *cf = new QHBoxLayout();
    m_f_stim = new QCheckBox("X(f)", this);
    m_f_resp = new QCheckBox("Y(f)", this);
    m_f_hid = new QCheckBox("H_id(f)", this);
    m_f_db = new QCheckBox("dB", this);
    m_f_norm = new QCheckBox(tr("Freq. response"), this);
    m_f_stim->setChecked(true);
    m_f_resp->setChecked(true);
    m_f_db->setChecked(true);
    m_f_hid->setToolTip(tr("Show the ideal magnitude response |H(f)| of the filter for comparison"));
    m_f_norm->setToolTip(tr("<span>Scale the spectrum of the response to an impulse so that it shows the "
                            "frequency response (use a rectangular window). Only for impulse stimuli without "
                            "noise and DC.</span>"));
    m_win = new QComboBox(this);
    for (const auto &w : window_list()) m_win->addItem(w.name, int(w.type));
    m_win->setCurrentIndex(m_win->findData(int(WindowType::Rectangular)));
    m_win->setToolTip(tr("Window for the spectral analysis"));
    m_lbl_win_par = new QLabel(this);
    m_win_par = new QLineEdit(this);
    m_win_par->setMaximumWidth(60);
    m_lbl_win_par2 = new QLabel(this);
    m_win_par2 = new QLineEdit(this);
    m_win_par2->setMaximumWidth(60);
    cf->addWidget(m_f_stim);
    cf->addWidget(m_f_resp);
    cf->addWidget(m_f_hid);
    cf->addWidget(m_f_db);
    cf->addWidget(m_f_norm);
    cf->addSpacing(12);
    cf->addWidget(new QLabel(tr("Window:"), this));
    cf->addWidget(m_win);
    cf->addWidget(m_lbl_win_par);
    cf->addWidget(m_win_par);
    cf->addWidget(m_lbl_win_par2);
    cf->addWidget(m_win_par2);
    cf->addStretch(1);
    vf->addLayout(cf);
    m_plot_f = new PlotWidget(this);
    vf->addWidget(m_plot_f, 1);
    m_tabs->addTab(wf, tr("Frequency"));

    // spectrogram (pyfda: time tab, scipy.signal.spectrogram)
    auto *ws = new QWidget(this);
    auto *vs = new QVBoxLayout(ws);
    vs->setContentsMargins(2, 2, 2, 2);
    auto *cs = new QHBoxLayout();
    m_s_sig = new QComboBox(this);
    m_s_sig->addItems({"y[n]", "x[n]"});
    m_s_mode = new QComboBox(this);
    m_s_mode->addItems({"PSD", tr("Magnitude"), tr("Phase")});
    m_s_mode->setToolTip(tr("Power spectral density, magnitude or phase of the short time Fourier transform"));
    m_s_db = new QCheckBox("dB", this);
    m_s_db->setChecked(true);
    m_s_nfft = new QSpinBox(this);
    m_s_nfft->setRange(4, 1 << 20);
    m_s_nfft->setValue(256);
    m_s_nfft->setToolTip(tr("Number of points per segment (the window of the Frequency tab is used)"));
    m_s_ovlp = new QSpinBox(this);
    m_s_ovlp->setRange(0, 1 << 20);
    m_s_ovlp->setValue(128);
    m_s_ovlp->setToolTip(tr("Number of overlapping points between segments"));
    cs->addWidget(new QLabel(tr("Signal:"), this));
    cs->addWidget(m_s_sig);
    cs->addWidget(m_s_mode);
    cs->addWidget(m_s_db);
    cs->addSpacing(12);
    cs->addWidget(new QLabel("NFFT =", this));
    cs->addWidget(m_s_nfft);
    cs->addWidget(new QLabel("N_OVLP =", this));
    cs->addWidget(m_s_ovlp);
    cs->addStretch(1);
    vs->addLayout(cs);
    m_plot_s = new PlotWidget(this);
    vs->addWidget(m_plot_s, 1);
    m_tabs->addTab(ws, tr("Spectrogram"));
    lay->addWidget(m_tabs, 1);
    m_info = new QLabel(this);
    m_info->setTextInteractionFlags(Qt::TextSelectableByMouse);
    lay->addWidget(m_info);

    // --- connections
    auto recalc = [this] { changed(); };
    connect(m_stim, &QComboBox::currentIndexChanged, this, [this] {
        m_p.stim = Stim(m_stim->currentData().toInt());
        updateVisibility();
        changed();
    });
    connect(m_chirp, &QComboBox::currentIndexChanged, this, [this](int i) {
        m_p.chirp = ChirpType(i);
        changed();
    });
    connect(m_bl, &QCheckBox::toggled, this, [this](bool on) {
        m_p.bl = on;
        changed();
    });
    connect(m_noise, &QComboBox::currentIndexChanged, this, [this](int i) {
        m_p.noise = Noise(i);
        updateVisibility();
        changed();
    });
    connect(m_noi, &QLineEdit::editingFinished, this, [this] { readEdits(); });
    connect(m_mls_b, &QSpinBox::valueChanged, this, [this](int b) {
        m_p.mls_b = b;
        changed();
    });
    connect(m_step_err, &QCheckBox::toggled, this, recalc);
    connect(m_N, &QSpinBox::valueChanged, this, recalc);
    connect(m_N_start, &QSpinBox::valueChanged, this, recalc);
    for (QCheckBox *c : {m_t_stim, m_t_resp, m_f_stim, m_f_resp, m_f_hid, m_f_norm})
        connect(c, &QCheckBox::toggled, this, [this] { redrawNow(); });
    for (QCheckBox *c : {m_t_db, m_f_db})
        connect(c, &QCheckBox::toggled, this, [this] {
            m_plot_t->keepView(false);
            m_plot_ti->keepView(false);
            m_plot_f->keepView(false);
            redrawNow();
        });
    connect(m_win, &QComboBox::currentIndexChanged, this, [this] {
        const auto &w = window_list()[m_win->currentIndex()];
        if (w.par_name) m_win_par->setText(fmt(w.par_default));
        if (w.par2_name) m_win_par2->setText(fmt(w.par2_default));
        updateVisibility();
        redrawNow();
    });
    connect(m_win_par, &QLineEdit::editingFinished, this, [this] { redrawNow(); });
    connect(m_win_par2, &QLineEdit::editingFinished, this, [this] { redrawNow(); });
    for (QComboBox *c : {m_s_sig, m_s_mode})
        connect(c, &QComboBox::currentIndexChanged, this, [this] {
            m_s_db->setEnabled(m_s_mode->currentIndex() != 2);
            m_plot_s->keepView(false);
            redrawNow();
        });
    connect(m_s_db, &QCheckBox::toggled, this, [this] { redrawNow(); });
    for (QSpinBox *b : {m_s_nfft, m_s_ovlp})
        connect(b, &QSpinBox::valueChanged, this, [this] {
            m_plot_s->keepView(false);
            redrawNow();
        });
    connect(m_formula, &QLineEdit::editingFinished, this, [this] {
        if (m_p.formula == m_formula->text().toStdString()) return;
        m_p.formula = m_formula->text().toStdString();
        changed();
    });
    connect(m_load, &QPushButton::clicked, this, [this] {
        const QString fn = QFileDialog::getOpenFileName(this, tr("Load stimulus"), config::dir("data"),
                                                        tr("Data files (*.csv *.txt *.dat *.wav *.npy);;All files (*)"));
        if (fn.isEmpty()) return;
        config::setDir("data", fn);
        loadStimFile(fn);
    });
    connect(m_file_norm, &QCheckBox::toggled, this, [this] {
        applyFileNorm();
        changed();
    });
    connect(exp, &QPushButton::clicked, this, [this] {
        QString fn = QFileDialog::getSaveFileName(this, tr("Export transient data"), config::dir("export"),
                                          tr("CSV (*.csv)"));
        if (fn.isEmpty()) return;
        config::setDir("export", fn);
        if (QFileInfo(fn).suffix().isEmpty()) fn += ".csv";
        saveCsv(fn);
    });
    updateEdits();
    updateVisibility();
}

void TranView::setParams(const StimParams &p) {
    m_p = p;
    for (QObject *o : std::initializer_list<QObject *>{m_stim, m_chirp, m_bl, m_noise, m_mls_b}) o->blockSignals(true);
    m_stim->setCurrentIndex(m_stim->findData(int(p.stim)));
    m_chirp->setCurrentIndex(int(p.chirp));
    m_bl->setChecked(p.bl);
    m_noise->setCurrentIndex(int(p.noise));
    m_mls_b->setValue(p.mls_b);
    for (QObject *o : std::initializer_list<QObject *>{m_stim, m_chirp, m_bl, m_noise, m_mls_b}) o->blockSignals(false);
    updateEdits();
    updateVisibility();
    changed();
}

bool TranView::loadStimFile(const QString &file_name) {
    DataTable t;
    try {
        t = load_data_file(QFile::encodeName(file_name).toStdString());
    } catch (const std::exception &e) {
        Logger::error(tr("Couldn't load '%1':\n%2").arg(file_name, e.what()));
        return false;
    }
    // first column with numbers that is not a time / index axis (e.g. n, t, x, y)
    int col = -1;
    for (size_t c = 0; c < t.n_cols && col < 0; ++c) {
        if (c + 1 < t.n_cols && is_time_column(t, c)) continue;
        for (size_t r = 0; r < t.n_rows; ++r)
            if (std::isfinite(t.at(r, c))) {
                col = int(c);
                break;
            }
    }
    if (col < 0) {
        Logger::error(tr("'%1' contains no numeric data.").arg(file_name));
        return false;
    }
    m_file_raw = t.column(size_t(col));
    m_file_raw_im = t.column_complex(size_t(col)) ? t.column_imag(size_t(col)) : Vec();
    for (double &v : m_file_raw)
        if (!std::isfinite(v)) v = 0;  // empty / non-numeric cells
    for (double &v : m_file_raw_im)
        if (!std::isfinite(v)) v = 0;
    m_p.file_name = QFileInfo(file_name).fileName().toStdString();
    m_file_lbl->setText(tr("%1: %2 samples, column '%3'")
                            .arg(QFileInfo(file_name).fileName())
                            .arg(m_file_raw.size())
                            .arg(QString::fromStdString(t.names[size_t(col)])) +
                        (m_file_raw_im.empty() ? QString() : tr(", complex")));
    m_file_lbl->setToolTip(file_name);
    Logger::info(tr("Loaded %1 samples from '%2' as stimulus.").arg(m_file_raw.size()).arg(file_name));
    applyFileNorm();
    m_p.stim = Stim::File;
    m_stim->blockSignals(true);
    m_stim->setCurrentIndex(m_stim->findData(int(Stim::File)));
    m_stim->blockSignals(false);
    updateVisibility();
    changed();
    return true;
}

void TranView::applyFileNorm() {
    if (m_file_raw.empty()) return;
    auto x = std::make_shared<Vec>(m_file_raw);
    auto xi = m_file_raw_im.empty() ? nullptr : std::make_shared<Vec>(m_file_raw_im);
    if (m_file_norm->isChecked()) {
        double mx = 0;  // max. |x|, also for complex data
        for (size_t i = 0; i < x->size(); ++i) mx = std::max(mx, std::hypot((*x)[i], xi ? (*xi)[i] : 0.0));
        if (mx > 0) {
            for (double &v : *x) v /= mx;
            if (xi)
                for (double &v : *xi) v /= mx;
        }
    }
    m_p.x_file = x;
    m_p.x_file_im = xi;
}

void TranView::setFormula(const QString &formula) {
    m_formula->setText(formula);
    m_p.formula = formula.toStdString();
    m_p.stim = Stim::Formula;
    m_stim->blockSignals(true);
    m_stim->setCurrentIndex(m_stim->findData(int(Stim::Formula)));
    m_stim->blockSignals(false);
    updateVisibility();
    changed();
}

void TranView::setFixpoint(const FxSpec &spec, bool on) {
    m_fx = spec;
    m_fx_on = on;
    changed();
}

void TranView::updateVisibility() {
    for (auto it = m_param_widgets.begin(); it != m_param_widgets.end(); ++it)
        it.value()->setVisible(stim_uses(m_p.stim, it.key().toStdString()));
    m_chirp->setVisible(m_p.stim == Stim::Chirp);
    m_wdg_formula->setVisible(m_p.stim == Stim::Formula);
    m_wdg_file->setVisible(m_p.stim == Stim::File);
    m_bl->setVisible(stim_uses(m_p.stim, "bl"));
    m_step_err->setVisible(m_p.stim == Stim::Step);
    m_noi->setVisible(m_p.noise != Noise::None);
    m_lbl_mls->setVisible(m_p.noise == Noise::MLS);
    m_mls_b->setVisible(m_p.noise == Noise::MLS);
    m_f_norm->setVisible(is_impulse(m_p.stim));
    const auto &w = window_list()[std::max(0, m_win->currentIndex())];
    m_lbl_win_par->setVisible(w.par_name != nullptr);
    m_win_par->setVisible(w.par_name != nullptr);
    if (w.par_name) m_lbl_win_par->setText(QString(w.par_name) + ":");
    m_lbl_win_par2->setVisible(w.par2_name != nullptr);
    m_win_par2->setVisible(w.par2_name != nullptr);
    if (w.par2_name) m_lbl_win_par2->setText(QString(w.par2_name) + ":");
}

// display value <-> normalized value
void TranView::updateEdits() {
    const double fs = m_ctx.f_s;
    auto set = [&](const QString &k, double v) {
        if (F_KEYS.contains(k)) v *= fs;
        else if (T_KEYS.contains(k)) v /= fs;
        m_edits[k]->setText(fmt(v));
        m_edits[k]->setStyleSheet(QString());
    };
    set("a1", m_p.a1);
    set("a2", m_p.a2);
    set("dc", m_p.dc);
    m_edits["a1"]->setText(fmtCplx(m_p.a1, m_p.a1_im));
    m_edits["a2"]->setText(fmtCplx(m_p.a2, m_p.a2_im));
    m_edits["dc"]->setText(fmtCplx(m_p.dc, m_p.dc_im));
    set("phi1", m_p.phi1);
    set("phi2", m_p.phi2);
    set("f1", m_p.f1);
    set("f2", m_p.f2);
    set("t1", m_p.t1);
    set("t2", m_p.t2);
    set("tw", m_p.tw);
    set("bw1", m_p.bw1);
    set("bw2", m_p.bw2);
    set("duty", m_p.duty);
    set("n1", m_p.n1);
    m_noi->setText(fmtCplx(m_p.noi, m_p.noi_im));
    const QString fu = m_ctx.f_label;
    for (const QString &k : F_KEYS) m_edits[k]->setToolTip(m_edits[k]->toolTip().section(" (", 0, 0) + " (" + fu + ")");
    const QString tu = m_ctx.unit_to_hz == 0 ? tr("samples") : m_ctx.t_label;
    for (const QString &k : T_KEYS) m_edits[k]->setToolTip(m_edits[k]->toolTip().section(" (", 0, 0) + " (" + tu + ")");
}

void TranView::readEdits() {
    const double fs = m_ctx.f_s;
    bool all_ok = true;
    auto get = [&](const QString &k, double &target) {
        QLineEdit *e = m_edits[k];
        double v;
        if (!parse(e, v)) {
            e->setStyleSheet("QLineEdit {background: #ffd0d0}");
            all_ok = false;
            return;
        }
        e->setStyleSheet(QString());
        if (F_KEYS.contains(k)) v /= fs;
        else if (T_KEYS.contains(k)) v *= fs;
        target = v;
    };
    auto getc = [&](const QString &k, double &re, double &im) {
        QLineEdit *e = m_edits[k];
        if (!parseCplx(e, re, im)) {
            e->setStyleSheet("QLineEdit {background: #ffd0d0}");
            all_ok = false;
            return;
        }
        e->setStyleSheet(QString());
    };
    getc("a1", m_p.a1, m_p.a1_im);
    getc("a2", m_p.a2, m_p.a2_im);
    getc("dc", m_p.dc, m_p.dc_im);
    get("phi1", m_p.phi1);
    get("phi2", m_p.phi2);
    get("f1", m_p.f1);
    get("f2", m_p.f2);
    get("t1", m_p.t1);
    get("t2", m_p.t2);
    get("tw", m_p.tw);
    get("bw1", m_p.bw1);
    get("bw2", m_p.bw2);
    get("duty", m_p.duty);
    double n1 = m_p.n1;
    get("n1", n1);
    m_p.n1 = std::max(1, int(std::lround(n1)));
    double noi = m_p.noi, noi_im = m_p.noi_im;
    if (parseCplx(m_noi, noi, noi_im)) {
        m_p.noi = noi;
        m_p.noi_im = noi_im;
    }
    if (!all_ok) Logger::warning(tr("Invalid stimulus parameter, the previous value is used."));
    changed();
}

void TranView::changed() {
    m_dirty = true;
    m_plot_t->keepView(false);
    m_plot_ti->keepView(false);
    m_plot_f->keepView(false);
    m_plot_s->keepView(false);
    if (isVisible()) redrawNow();
    else m_needs_redraw = true;
}

void TranView::designChanged() {
    updateEdits();
    m_dirty = true;
}

void TranView::calc() {
    m_dirty = false;
    m_error.clear();
    m_x.clear();
    m_y.clear();
    m_xi.clear();
    m_yi.clear();
    m_yi_float.clear();
    m_cmplx = false;
    if (!m_design) return;
    int N = m_N->value();
    if (N == 0) N = std::max(100, m_design->fir ? int(m_design->ba.b.size()) + 5 : impz_len(m_design->zpk));
    if (m_N->value() == 0 && m_p.stim == Stim::File && m_p.x_file) N = std::max(1, int(m_p.x_file->size()));
    m_n_start = m_N_start->value();
    const int n_end = m_n_start + N;
    m_p.f_s = m_ctx.f_s;
    m_formula->setStyleSheet(QString());
    try {
        const CVec xc = calc_stimulus_c(m_p, n_end, &m_cmplx);
        m_x.resize(xc.size());
        for (size_t i = 0; i < xc.size(); ++i) m_x[i] = xc[i].real();
        if (m_cmplx) {
            m_xi.resize(xc.size());
            for (size_t i = 0; i < xc.size(); ++i) m_xi[i] = xc[i].imag();
        }
    } catch (const std::exception &e) {
        m_error = e.what();
        Logger::warning(m_error);
        if (m_p.stim == Stim::Formula) m_formula->setStyleSheet("QLineEdit {background: #ffd0d0}");
        m_x.clear();
        m_cmplx = false;
        m_lbl_cmplx->setVisible(false);
        return;
    }
    m_lbl_cmplx->setVisible(m_cmplx);
    // the coefficients are real: real and imaginary parts are filtered separately
    auto filt = [&](const Vec &x) {
        return m_design->sos.empty() ? lfilter(m_design->ba.b, m_design->ba.a, x) : sosfilt(m_design->sos, x);
    };
    m_y = filt(m_x);
    if (m_cmplx) m_yi = filt(m_xi);
    m_y_float.clear();
    m_fx_info.clear();
    if (m_fx_on) {
        try {
            auto fx = [&](const Vec &x) {
                return m_design->fir ? fx_filter_fir(m_design->ba.b, m_fx, x) : fx_filter_sos(m_design->sos, m_fx, x);
            };
            FxResult r = fx(m_x);
            m_y_float = m_y;
            m_y = r.y;
            double err = 0;
            for (size_t i = size_t(m_n_start); i < m_y.size(); ++i) err = std::max(err, std::fabs(m_y[i] - m_y_float[i]));
            if (m_cmplx) {  // two real fixpoint filters for the real and the imaginary part
                const FxResult ri = fx(m_xi);
                m_yi_float = m_yi;
                m_yi = ri.y;
                r.n_over_i += ri.n_over_i;
                r.n_over_acc += ri.n_over_acc;
                r.n_over_o += ri.n_over_o;
                for (size_t i = size_t(m_n_start); i < m_yi.size(); ++i)
                    err = std::max(err, std::fabs(m_yi[i] - m_yi_float[i]));
            }
            m_fx_info = tr("Fixpoint: overflows input %1, accumulator %2, output %3, coefficients %4; "
                           "max. |y - y_float| = %5")
                            .arg(r.n_over_i)
                            .arg(r.n_over_acc)
                            .arg(r.n_over_o)
                            .arg(r.n_over_coeff)
                            .arg(err, 0, 'g', 4);
            if (r.n_over_i + r.n_over_acc + r.n_over_o + r.n_over_coeff > 0)
                Logger::warning(tr("Fixpoint simulation: %1 overflow(s).")
                                    .arg(r.n_over_i + r.n_over_acc + r.n_over_o + r.n_over_coeff));
        } catch (const std::exception &e) {
            m_fx_info = tr("Fixpoint simulation failed: %1").arg(e.what());
            Logger::error(m_fx_info);
        }
    }
    if (m_p.stim == Stim::Step && m_step_err->isChecked()) {
        // settling error: subtract the DC response from the response after the step
        const double dc = std::abs((m_design->sos.empty() ? freqz(m_design->ba, {0.0}) : freqz(m_design->sos, {0.0}))[0]);
        const int n0 = std::max(m_n_start, int(std::nearbyint(m_p.t1)));
        for (int n = std::max(n0, 0); n < n_end; ++n) m_y[size_t(n)] -= dc;
    }
}

void TranView::redraw() {
    if (m_dirty) calc();
    drawTime();
    drawFreq();
    drawSpgr();
}

Vec TranView::analysisWindow(int N) const {
    const auto &wi = window_list()[m_win->currentIndex()];
    double par = wi.par_default, par2 = wi.par2_default;
    if (wi.par_name && !parse(m_win_par, par)) par = wi.par_default;
    if (wi.par2_name && !parse(m_win_par2, par2)) par2 = wi.par2_default;
    return fft_window(wi.type, N, par, par2);
}

void TranView::drawSpgr() {
    m_plot_s->clear();
    m_plot_s->setXLabel(m_ctx.t_label);
    m_plot_s->setYLabel(m_ctx.f_label);
    if (m_x.empty()) {
        m_plot_s->setMessage(m_error.isEmpty() ? tr("No filter designed") : m_error);
        m_plot_s->autoscale();
        return;
    }
    const bool use_y = m_s_sig->currentIndex() == 0;
    const Vec &sr = use_y ? m_y : m_x, &si = use_y ? m_yi : m_xi;
    const Vec s(sr.begin() + m_n_start, sr.end());
    const int N = int(s.size());
    int nfft = m_s_nfft->value(), novl = m_s_ovlp->value();
    // like pyfda: NFFT <= N and N_OVLP < NFFT
    if (nfft > N) nfft = N;
    if (novl >= nfft) novl = 0;
    const auto &wi = window_list()[m_win->currentIndex()];
    const int mode_i = m_s_mode->currentIndex();
    const SpgrMode mode = mode_i == 0 ? SpgrMode::PSD : mode_i == 1 ? SpgrMode::Magnitude : SpgrMode::Angle;
    Spectrogram sp;
    try {
        const Vec win = analysisWindow(nfft);
        if (m_cmplx) {  // two-sided
            CVec c(s.size());
            for (size_t i = 0; i < s.size(); ++i) c[i] = cplx(s[i], si[size_t(m_n_start) + i]);
            sp = spectrogram(c, m_ctx.f_s, win, novl, mode, true);
        } else {
            sp = spectrogram(s, m_ctx.f_s, win, novl, mode, true);
        }
    } catch (const std::exception &e) {
        m_plot_s->setMessage(e.what());
        m_plot_s->autoscale();
        return;
    }
    const bool log = m_s_db->isChecked() && mode != SpgrMode::Angle;
    const double db_scale = mode == SpgrMode::PSD ? 10 : 20;
    const int nt = int(sp.t.size()), nf = int(sp.f.size());
    std::vector<double> z(size_t(nt) * size_t(nf));
    double zmax = -std::numeric_limits<double>::infinity(), zmin = -zmax;
    for (int i = 0; i < nt; ++i)
        for (int k = 0; k < nf; ++k) {
            double v = sp.s[size_t(i)][size_t(k)];
            if (log) v = db_scale * std::log10(std::max(std::fabs(v), 1e-300));
            z[size_t(i) * size_t(nf) + size_t(k)] = v;
            if (std::isfinite(v)) {
                zmax = std::max(zmax, v);
                zmin = std::min(zmin, v);
            }
        }
    if (log) zmin = std::max(zmin, zmax - 120);  // 120 dB dynamic range like the other plots
    if (mode == SpgrMode::Angle) {
        zmin = -PI;
        zmax = PI;
    }
    if (!(zmax > zmin)) zmax = zmin + 1;
    QImage img(nt, nf, QImage::Format_RGB32);
    for (int i = 0; i < nt; ++i)
        for (int k = 0; k < nf; ++k)
            img.setPixel(i, nf - 1 - k, PlotWidget::colormap((z[size_t(i) * size_t(nf) + size_t(k)] - zmin) / (zmax - zmin)));
    // pixel edges: segments are spaced by (nfft - novl) / f_S, bins by f_S / nfft
    const double t_scale = m_ctx.unit_to_hz == 0 ? m_ctx.f_s : 1.0;  // normalized: time in samples
    const double dt = (nfft - novl) / m_ctx.f_s * t_scale, df = m_ctx.f_s / nfft;
    const double t0 = (sp.t.front() + m_n_start / m_ctx.f_s) * t_scale;
    const QString sym = use_y ? "Y" : "X";
    const QString zl = mode == SpgrMode::PSD ? QString("S_%1%1 / %2").arg(sym.toLower(), log ? "dB re W/Hz" : "W/Hz")
                       : mode == SpgrMode::Magnitude ? QString("|%1| / %2").arg(sym, log ? "dBV" : "V")
                                                     : QString::fromUtf8("∠%1 / rad").arg(sym);
    const double f0 = sp.f.front();
    m_plot_s->setImage(img, t0 - dt / 2, t0 + (nt - 0.5) * dt, f0 - df / 2, f0 + (nf - 0.5) * df, zmin, zmax, zl);
    m_plot_s->setTitle(tr("Spectrogram of %1, %2 window, NFFT = %3, N_OVLP = %4")
                           .arg(use_y ? "y[n]" : "x[n]", wi.name)
                           .arg(nfft)
                           .arg(novl));
    m_plot_s->autoscale();
}

void TranView::drawTime() {
    m_plot_t->clear();
    m_plot_ti->clear();
    m_plot_t->setXLabel(m_ctx.t_label);
    m_plot_ti->setXLabel(m_ctx.t_label);
    m_plot_ti->setVisible(m_cmplx && !m_x.empty());
    if (m_x.empty()) {
        m_plot_t->setMessage(m_error.isEmpty() ? tr("No filter designed") : m_error);
        m_plot_t->autoscale();
        return;
    }
    const bool log = m_t_db->isChecked();
    const bool norm = m_ctx.unit_to_hz == 0;
    const int n_end = int(m_x.size());
    QVector<double> t;
    for (int n = m_n_start; n < n_end; ++n) t << (norm ? n : n / m_ctx.f_s);
    auto data = [&](const Vec &v) {
        QVector<double> d;
        for (int n = m_n_start; n < n_end; ++n) d << (log ? db(v[size_t(n)]) : v[size_t(n)]);
        return d;
    };
    const bool step_err = m_p.stim == Stim::Step && m_step_err->isChecked();
    // real part (or real signal) in the upper plot, imaginary part in the lower one like pyfda
    auto draw = [&](PlotWidget *p, const Vec &x, const Vec &y, const Vec &y_float, const QString &part) {
        auto name = [&](const QString &s) { return part.isEmpty() ? s : QString("%1{%2}").arg(part, s); };
        if (m_t_stim->isChecked()) p->addCurve(t, data(x), PlotWidget::color(1), name("x[n]"), PlotWidget::Style::Line);
        if (m_t_resp->isChecked() && !y_float.empty()) {
            PlotWidget::Curve c;
            c.x = t;
            c.y = data(y_float);
            c.color = PlotWidget::color(2);
            c.name = name("y_float[n]");
            c.pen = Qt::DashLine;
            p->addCurve(c);
        }
        if (m_t_resp->isChecked())
            p->addCurve(t, data(y), PlotWidget::color(0),
                        name(step_err ? QString::fromUtf8("ε[n]") : QString("y[n]")) + (y_float.empty() ? "" : " (fixpoint)"),
                        log || t.size() > 300 ? PlotWidget::Style::Line : PlotWidget::Style::Stem);
        const QString unit = log ? tr("Magnitude / dB") : tr("Amplitude");
        p->setYLabel(part.isEmpty() ? unit : QString("%1 %2").arg(part == "Re" ? tr("Real part,") : tr("Imag. part,"), unit));
        p->setYLimits(log ? -80 : NaN, NaN);
    };
    draw(m_plot_t, m_x, m_y, m_y_float, m_cmplx ? "Re" : "");
    m_plot_t->setTitle((m_y_float.empty() ? QString() : tr("Fixpoint ")) + QString::fromStdString(stim_title(m_p, step_err)));
    m_plot_t->autoscale();
    if (m_cmplx) {
        draw(m_plot_ti, m_xi, m_yi, m_yi_float, "Im");
        m_plot_ti->autoscale();
    }
}

void TranView::drawFreq() {
    m_plot_f->clear();
    m_plot_f->setXLabel(m_ctx.f_label);
    m_info->clear();
    if (m_x.empty()) {
        m_plot_f->setMessage(m_error.isEmpty() ? tr("No filter designed") : m_error);
        m_plot_f->autoscale();
        return;
    }
    const int N = int(m_x.size()) - m_n_start;
    const auto &wi = window_list()[m_win->currentIndex()];
    Vec win;
    try {
        win = analysisWindow(N);
    } catch (const std::exception &e) {
        m_plot_f->setMessage(e.what());
        m_plot_f->autoscale();
        return;
    }
    const double cgain = window_cgain(win), nenbw = window_nenbw(win);
    auto spectrum = [&](const Vec &re, const Vec &im) {
        if (!m_cmplx) return windowed_fft(Vec(re.begin() + m_n_start, re.end()), win);
        CVec c(static_cast<size_t>(N));
        for (int i = 0; i < N; ++i) c[size_t(i)] = cplx(re[size_t(m_n_start + i)], im[size_t(m_n_start + i)]);
        return windowed_fft(c, win);
    };
    const CVec X = spectrum(m_x, m_xi), Y = spectrum(m_y, m_yi);

    // frequency response from an impulse: scale by N * cgain / energy of the impulse
    const bool freq_resp =
        m_f_norm->isVisible() && m_f_norm->isChecked() && m_p.noise == Noise::None && m_p.dc == 0 && m_p.dc_im == 0;
    const double s_imp = impulse_scale(m_p);
    const double scale = freq_resp && s_imp > 0 && std::isfinite(s_imp) ? N * cgain / s_imp : 1.0;
    if (freq_resp && wi.type != WindowType::Rectangular)
        Logger::warning(tr("Use a rectangular window for a correctly scaled FFT of an impulse."));
    auto power = [&](const CVec &S) {
        double p = 0;
        for (const cplx &v : S) p += std::norm(v);
        return p * (freq_resp ? scale : 1.0) / nenbw;
    };
    const double p_x = power(X), p_y = power(Y);
    // single-sided display 0 ... f_S / 2, two-sided -f_S / 2 ... f_S / 2 for complex signals
    auto display = [&](const CVec &S) {
        if (m_cmplx) return fftshift(S);
        return freq_resp ? CVec(S.begin(), S.begin() + (N / 2 + 1)) : ssb_spectrum(S);
    };
    const CVec Xd = display(X), Yd = display(Y);
    const int k0 = m_cmplx ? N / 2 : 0;  // index of f = 0
    const bool log = m_f_db->isChecked();
    auto curve = [&](const CVec &S) {
        QVector<double> f, a;
        for (size_t k = 0; k < S.size(); ++k) {
            f << (double(k) - k0) * m_ctx.f_s / N;
            const double m = std::abs(S[k]) * scale;
            a << (log ? db(m) : m);
        }
        return std::make_pair(f, a);
    };
    auto pstr = [&](double p) {
        return log ? QString("P = %1 dB").arg(10 * std::log10(std::max(p, 1e-300)), 0, 'f', 2)
                   : QString("P = %1").arg(p, 0, 'g', 4);
    };
    if (m_f_stim->isChecked()) {
        const auto c = curve(Xd);
        m_plot_f->addCurve(c.first, c.second, PlotWidget::color(1), "X(f): " + pstr(p_x));
    }
    if (m_f_resp->isChecked()) {
        const auto c = curve(Yd);
        m_plot_f->addCurve(c.first, c.second, PlotWidget::color(0), "Y(f): " + pstr(p_y));
    }
    if (m_f_hid->isChecked()) {
        const int n = 2048;
        Vec w(n);
        QVector<double> f(n), a(n);
        const double w0 = m_cmplx ? -PI : 0;  // two-sided for complex signals
        const double span = m_cmplx ? 2 * PI : PI;
        for (int i = 0; i < n; ++i) {
            w[size_t(i)] = w0 + span * i / (n - 1);
            f[i] = w[size_t(i)] / (2 * PI) * m_ctx.f_s;
        }
        const CVec H = m_design->sos.empty() ? freqz(m_design->ba, w) : freqz(m_design->sos, w);
        for (int i = 0; i < n; ++i) a[i] = log ? db(std::abs(H[size_t(i)])) : std::abs(H[size_t(i)]);
        PlotWidget::Curve c;
        c.x = f;
        c.y = a;
        c.color = PlotWidget::color(2);
        c.name = "|H_id(f)|";
        c.pen = Qt::DashLine;
        m_plot_f->addCurve(c);
    }
    m_plot_f->setTitle(freq_resp ? tr("Frequency response from %1").arg(QString::fromStdString(stim_title(m_p)))
                                 : tr("Spectrum, %1 window, N = %2").arg(wi.name).arg(N));
    m_plot_f->setYLabel(log ? tr("|X|, |Y| / dB") : tr("|X|, |Y|"));
    m_plot_f->setYLimits(log ? -120 : NaN, NaN);
    m_plot_f->autoscale();
    m_info->setText(tr("N = %1 points, window: %2 (NENBW = %3 bins, CGAIN = %4)")
                        .arg(N)
                        .arg(wi.name)
                        .arg(nenbw, 0, 'f', 3)
                        .arg(cgain, 0, 'f', 3) +
                    (m_fx_info.isEmpty() ? QString() : "\n" + m_fx_info));
}

bool TranView::saveCsv(const QString &file_name) {
    if (m_dirty) calc();
    if (m_x.empty()) {
        Logger::error(tr("No transient data to export."));
        return false;
    }
    QFile f(file_name);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        Logger::error(tr("Couldn't write '%1'.").arg(file_name));
        return false;
    }
    const CsvFormat fmt = config::csvFormat();
    const QChar d(fmt.delimiter);
    auto num = [&](double v) { return QString::fromStdString(csv_number(v, fmt)); };
    QTextStream s(&f);
    if (m_cmplx) {  // real and imaginary parts in separate columns
        s << "n" << d << "t" << d << "x_re" << d << "x_im" << d << "y_re" << d << "y_im\n";
        for (size_t n = size_t(m_n_start); n < m_x.size(); ++n)
            s << n << d << num(double(n) / m_ctx.f_s) << d << num(m_x[n]) << d << num(m_xi[n]) << d << num(m_y[n]) << d
              << num(m_yi[n]) << "\n";
    } else {
        s << "n" << d << "t" << d << "x" << d << "y\n";
        for (size_t n = size_t(m_n_start); n < m_x.size(); ++n)
            s << n << d << num(double(n) / m_ctx.f_s) << d << num(m_x[n]) << d << num(m_y[n]) << "\n";
    }
    Logger::info(tr("Exported transient data to '%1'.").arg(file_name));
    return true;
}
