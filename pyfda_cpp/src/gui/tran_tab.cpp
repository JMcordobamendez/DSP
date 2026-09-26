#include "tran_tab.hpp"

#include "conversions.hpp"
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
    m_noi->setToolTip(tr("Noise amplitude"));
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
    add(r2, "a1", "A1:", tr("Amplitude of the first signal component"));
    add(r2, "a2", "A2:", tr("Amplitude of the second component (modulation index for PM / FM, "
                            "modulation depth for PWM)"));
    add(r2, "phi1", QString::fromUtf8("φ1 / °:"), tr("Phase of the first component in degrees"));
    add(r2, "phi2", QString::fromUtf8("φ2 / °:"), tr("Phase of the second component in degrees"));
    add(r2, "duty", QString::fromUtf8("α:"), tr("Duty cycle, 0 ... 1"));
    add(r2, "n1", "N1:", tr("Order of the periodic sinc (Dirichlet) function"));
    add(r2, "dc", "DC:", tr("DC offset added to the stimulus"));
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
    cf->addStretch(1);
    vf->addLayout(cf);
    m_plot_f = new PlotWidget(this);
    vf->addWidget(m_plot_f, 1);
    m_tabs->addTab(wf, tr("Frequency"));
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
            m_plot_f->keepView(false);
            redrawNow();
        });
    connect(m_win, &QComboBox::currentIndexChanged, this, [this] {
        const auto &w = window_list()[m_win->currentIndex()];
        if (w.par_name) m_win_par->setText(fmt(w.par_default));
        updateVisibility();
        redrawNow();
    });
    connect(m_win_par, &QLineEdit::editingFinished, this, [this] { redrawNow(); });
    connect(exp, &QPushButton::clicked, this, [this] {
        QString fn = QFileDialog::getSaveFileName(this, tr("Export transient data"), QString(), tr("CSV (*.csv)"));
        if (fn.isEmpty()) return;
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

void TranView::updateVisibility() {
    for (auto it = m_param_widgets.begin(); it != m_param_widgets.end(); ++it)
        it.value()->setVisible(stim_uses(m_p.stim, it.key().toStdString()));
    m_chirp->setVisible(m_p.stim == Stim::Chirp);
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
    set("dc", m_p.dc);
    m_noi->setText(fmt(m_p.noi));
    const QString fu = m_ctx.unit_to_hz == 0 ? tr("normalized to f_S") : m_ctx.f_label;
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
    get("a1", m_p.a1);
    get("a2", m_p.a2);
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
    get("dc", m_p.dc);
    double n1 = m_p.n1;
    get("n1", n1);
    m_p.n1 = std::max(1, int(std::lround(n1)));
    double noi = m_p.noi;
    if (parse(m_noi, noi)) m_p.noi = noi;
    if (!all_ok) Logger::warning(tr("Invalid stimulus parameter, the previous value is used."));
    changed();
}

void TranView::changed() {
    m_dirty = true;
    m_plot_t->keepView(false);
    m_plot_f->keepView(false);
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
    if (!m_design) return;
    int N = m_N->value();
    if (N == 0) N = std::max(100, m_design->fir ? int(m_design->ba.b.size()) + 5 : impz_len(m_design->zpk));
    m_n_start = m_N_start->value();
    const int n_end = m_n_start + N;
    try {
        m_x = calc_stimulus(m_p, n_end);
    } catch (const std::exception &e) {
        m_error = e.what();
        Logger::warning(m_error);
        m_x.clear();
        return;
    }
    m_y = m_design->sos.empty() ? lfilter(m_design->ba.b, m_design->ba.a, m_x) : sosfilt(m_design->sos, m_x);
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
}

void TranView::drawTime() {
    m_plot_t->clear();
    m_plot_t->setXLabel(m_ctx.t_label);
    if (m_x.empty()) {
        m_plot_t->setMessage(m_error.isEmpty() ? tr("No filter designed") : m_error);
        m_plot_t->autoscale();
        return;
    }
    const bool log = m_t_db->isChecked();
    const bool norm = m_ctx.unit_to_hz == 0;
    const int n_end = int(m_x.size());
    QVector<double> t, x, y;
    for (int n = m_n_start; n < n_end; ++n) {
        t << (norm ? n : n / m_ctx.f_s);
        x << (log ? db(m_x[size_t(n)]) : m_x[size_t(n)]);
        y << (log ? db(m_y[size_t(n)]) : m_y[size_t(n)]);
    }
    const bool step_err = m_p.stim == Stim::Step && m_step_err->isChecked();
    if (m_t_stim->isChecked()) m_plot_t->addCurve(t, x, PlotWidget::color(1), "x[n]", PlotWidget::Style::Line);
    if (m_t_resp->isChecked())
        m_plot_t->addCurve(t, y, PlotWidget::color(0), step_err ? QString::fromUtf8("ε[n]") : "y[n]",
                           log || t.size() > 300 ? PlotWidget::Style::Line : PlotWidget::Style::Stem);
    m_plot_t->setTitle(QString::fromStdString(stim_title(m_p, step_err)));
    m_plot_t->setYLabel(log ? tr("Magnitude / dB") : tr("Amplitude"));
    m_plot_t->setYLimits(log ? -80 : NaN, NaN);
    m_plot_t->autoscale();
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
    const Vec xs(m_x.begin() + m_n_start, m_x.end()), ys(m_y.begin() + m_n_start, m_y.end());
    const auto &wi = window_list()[m_win->currentIndex()];
    double par = wi.par_default;
    if (wi.par_name && !parse(m_win_par, par)) par = wi.par_default;
    const Vec win = fft_window(wi.type, N, par);
    const double cgain = window_cgain(win), nenbw = window_nenbw(win);
    const CVec X = windowed_fft(xs, win), Y = windowed_fft(ys, win);

    // frequency response from an impulse: scale by N * cgain / energy of the impulse
    const bool freq_resp = m_f_norm->isVisible() && m_f_norm->isChecked() && m_p.noise == Noise::None && m_p.dc == 0;
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
    // single-sided display: 0 ... f_S / 2
    const CVec Xd = freq_resp ? CVec(X.begin(), X.begin() + (N / 2 + 1)) : ssb_spectrum(X);
    const CVec Yd = freq_resp ? CVec(Y.begin(), Y.begin() + (N / 2 + 1)) : ssb_spectrum(Y);
    const bool log = m_f_db->isChecked();
    auto curve = [&](const CVec &S) {
        QVector<double> f, a;
        for (size_t k = 0; k < S.size(); ++k) {
            f << double(k) * m_ctx.f_s / N;
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
        for (int i = 0; i < n; ++i) {
            w[size_t(i)] = PI * i / (n - 1);
            f[i] = m_ctx.f_s / 2 * i / (n - 1);
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
                        .arg(cgain, 0, 'f', 3));
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
    QTextStream s(&f);
    s << "n,t,x,y\n";
    for (size_t n = size_t(m_n_start); n < m_x.size(); ++n)
        s << n << "," << QString::number(double(n) / m_ctx.f_s, 'g', 17) << "," << QString::number(m_x[n], 'g', 17)
          << "," << QString::number(m_y[n], 'g', 17) << "\n";
    Logger::info(tr("Exported transient data to '%1'.").arg(file_name));
    return true;
}
