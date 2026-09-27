#include "spec_panel.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>
#include <QVBoxLayout>

#include <cmath>

using namespace pyfda;

namespace {
// frequency units like pyfda's freq_units.py
struct UnitDef {
    const char *key;     // stored in filter files
    const char *text;    // combo box text
    double to_hz;        // 0: normalized
    double norm_fs;      // f_S for normalized units
    const char *f_label;
    const char *t_label;
};
const UnitDef UNITS[] = {
    {"f_S", "f_S (norm.)", 0, 1, "F = f / f_S", "n"},
    {"f_Ny", "f_Ny (norm.)", 0, 2, "F = 2f / f_S", "n"},
    {"mHz", "mHz", 1e-3, 0, "f / mHz", "t / ks"},
    {"Hz", "Hz", 1, 0, "f / Hz", "t / s"},
    {"kHz", "kHz", 1e3, 0, "f / kHz", "t / ms"},
    {"MHz", "MHz", 1e6, 0, "f / MHz", "t / \xC2\xB5s"},
    {"GHz", "GHz", 1e9, 0, "f / GHz", "t / ns"},
};
const UnitDef &unit_def(int i) { return UNITS[std::clamp(i, 0, int(std::size(UNITS)) - 1)]; }

const DesignMethod IIR_METHODS[] = {DesignMethod::Butter, DesignMethod::Cheby1, DesignMethod::Cheby2,
                                    DesignMethod::Ellip, DesignMethod::Bessel, DesignMethod::ManualIIR};
const DesignMethod FIR_METHODS[] = {DesignMethod::Equiripple, DesignMethod::Firwin, DesignMethod::MovingAverage,
                                    DesignMethod::Delay, DesignMethod::ManualFIR};
}  // namespace

SpecPanel::SpecPanel(QWidget *parent) : QWidget(parent) {
    auto *lay = new QVBoxLayout(this);

    // --- filter selection -------------------------------------------------
    auto *grpSel = new QGroupBox(tr("Filter"), this);
    auto *fSel = new QFormLayout(grpSel);
    m_rt = new QComboBox(this);
    m_rt->addItems({tr("Lowpass"), tr("Highpass"), tr("Bandpass"), tr("Bandstop")});
    m_rt->setToolTip(tr("Response type"));
    m_ft = new QComboBox(this);
    m_ft->addItems({"IIR", "FIR"});
    m_method = new QComboBox(this);
    m_method->setToolTip(tr("Design method"));
    fSel->addRow(tr("Response:"), m_rt);
    fSel->addRow(tr("Type:"), m_ft);
    fSel->addRow(tr("Method:"), m_method);

    auto *wOrder = m_wOrder = new QWidget(this);
    auto *hOrder = new QHBoxLayout(wOrder);
    hOrder->setContentsMargins(0, 0, 0, 0);
    m_min = new QRadioButton(tr("Minimum"), this);
    m_man = new QRadioButton(tr("Manual"), this);
    m_min->setChecked(true);
    m_min->setToolTip(tr("Calculate the minimum order that meets the specifications"));
    m_man->setToolTip(tr("Specify the order N manually"));
    m_N = new QSpinBox(this);
    m_N->setRange(1, 5000);
    m_N->setValue(10);
    m_N->setToolTip(tr("Filter order N (number of taps - 1 for FIR filters)"));
    hOrder->addWidget(m_min);
    hOrder->addWidget(m_man);
    hOrder->addWidget(new QLabel("N =", this));
    hOrder->addWidget(m_N, 1);
    fSel->addRow(tr("Order:"), wOrder);

    // moving average: stages and normalization (pyfda filter_widgets/ma.py)
    m_lstages = new QLabel(tr("Stages:"), this);
    auto *wMa = new QWidget(this);
    auto *hMa = new QHBoxLayout(wMa);
    hMa->setContentsMargins(0, 0, 0, 0);
    m_stages = new QSpinBox(this);
    m_stages->setRange(1, 100);
    m_stages->setToolTip(tr("Number of cascaded moving average stages; N is the number of delays per stage"));
    m_norm = new QCheckBox(tr("Normalize"), this);
    m_norm->setChecked(true);
    m_norm->setToolTip(tr("Normalize to |H(f)|max = 1"));
    hMa->addWidget(m_stages, 1);
    hMa->addWidget(m_norm);
    fSel->addRow(m_lstages, wMa);
    m_hint = new QLabel(this);
    m_hint->setWordWrap(true);
    fSel->addRow(m_hint);

    m_lwindow = new QLabel(tr("Window:"), this);
    m_window = new QComboBox(this);
    for (const auto &w : window_list())
        if (w.fir) m_window->addItem(w.name, int(w.type));  // the others are for spectral analysis only
    m_window->setCurrentIndex(m_window->findText("Kaiser"));
    fSel->addRow(m_lwindow, m_window);
    m_lwinpar = new QLabel(tr("beta:"), this);
    m_winpar = new QLineEdit("10", this);
    fSel->addRow(m_lwinpar, m_winpar);
    m_lalg = new QLabel(tr("Order est.:"), this);
    m_alg = new QComboBox(this);
    m_alg->addItems({"ichige", "kaiser", "herrmann"});
    m_alg->setToolTip(tr("Algorithm for estimating the minimum FIR filter order"));
    fSel->addRow(m_lalg, m_alg);
    lay->addWidget(grpSel);

    // --- frequency specs ----------------------------------------------------
    auto *grpF = m_grpF = new QGroupBox(tr("Frequencies"), this);
    auto *gF = new QGridLayout(grpF);
    m_unit = new QComboBox(this);
    for (const UnitDef &u : UNITS) m_unit->addItem(u.text, QString(u.key));
    m_unit->setToolTip(tr("<span>Frequency unit. With <i>f_S (norm.)</i> all frequencies are "
                          "normalized to the sampling frequency f_S = 1, with <i>f_Ny (norm.)</i> to the "
                          "Nyquist frequency f_S / 2 = 1.</span>"));
    gF->addWidget(new QLabel(tr("Unit:"), this), 0, 0);
    gF->addWidget(m_unit, 0, 1);
    m_lfs = new QLabel("f_S =", this);
    m_fs = new QLineEdit("1", this);
    m_fs->setToolTip(tr("Sampling frequency"));
    gF->addWidget(m_lfs, 1, 0);
    gF->addWidget(m_fs, 1, 1);
    addRow(gF, 2, "F_PB", m_fpb, m_lfpb, tr("Pass band edge"));
    addRow(gF, 3, "F_SB", m_fsb, m_lfsb, tr("Stop band edge"));
    addRow(gF, 4, "F_PB2", m_fpb2, m_lfpb2, tr("Second pass band edge"));
    addRow(gF, 5, "F_SB2", m_fsb2, m_lfsb2, tr("Second stop band edge"));
    addRow(gF, 6, "F_C", m_fc, m_lfc, tr("Corner frequency (-3 dB for Butterworth, pass band edge for "
                                         "Chebyshev 1 and elliptic, stop band edge for Chebyshev 2)"));
    addRow(gF, 7, "F_C2", m_fc2, m_lfc2, tr("Second corner frequency"));
    lay->addWidget(grpF);

    // --- amplitude specs ----------------------------------------------------
    auto *grpA = m_grpA = new QGroupBox(tr("Amplitudes"), this);
    auto *gA = new QGridLayout(grpA);
    m_amp_unit = new QComboBox(this);
    m_amp_unit->addItems({"dB", "V", "W"});
    m_amp_unit->setToolTip(tr("<span>Unit of the amplitude specs: dB, linear deviation (V) or its square (W). "
                              "Pass band: IIR: A_dB = -20 log10(1 - A_V), FIR: A_dB = 20 log10((1 + A_V) / "
                              "(1 - A_V)); stop band: A_dB = -20 log10(A_V).</span>"));
    gA->addWidget(new QLabel(tr("Unit:"), this), 4, 0);
    gA->addWidget(m_amp_unit, 4, 1);
    addRow(gA, 0, "A_PB / dB", m_apb, m_lapb, tr("Maximum pass band ripple"));
    addRow(gA, 1, "A_SB / dB", m_asb, m_lasb, tr("Minimum stop band attenuation"));
    addRow(gA, 2, "W_PB", m_wpb, m_lwpb, tr("Pass band weight (equiripple, manual order)"));
    addRow(gA, 3, "W_SB", m_wsb, m_lwsb, tr("Stop band weight (equiripple, manual order)"));
    m_apb->setText("1");
    m_asb->setText("60");
    m_wpb->setText("1");
    m_wsb->setText("1");
    lay->addWidget(grpA);

    m_design = new QPushButton(tr("DESIGN FILTER"), this);
    QFont f = m_design->font();
    f.setBold(true);
    m_design->setFont(f);
    m_design->setMinimumHeight(32);
    m_design->setToolTip(tr("Design the filter with the current specifications (Enter)"));
    lay->addWidget(m_design);
    m_status = new QLabel(this);
    m_status->setWordWrap(true);
    m_status->setTextInteractionFlags(Qt::TextSelectableByMouse);
    lay->addWidget(m_status);
    lay->addStretch(1);

    // default frequency specs per response type (normalized to f_S)
    m_freqs[int(RespType::LP)] = {0.1, 0.35, 0.15, 0.3, 0.125, 0.325};
    m_freqs[int(RespType::HP)] = {0.15, 0.35, 0.1, 0.3, 0.125, 0.325};
    m_freqs[int(RespType::BP)] = {0.15, 0.3, 0.1, 0.35, 0.125, 0.325};
    m_freqs[int(RespType::BS)] = {0.1, 0.35, 0.15, 0.3, 0.125, 0.325};
    m_ft->setCurrentIndex(0);
    for (DesignMethod m : IIR_METHODS) m_method->addItem(method_name(m), int(m));
    m_method->setCurrentIndex(3);  // elliptic
    loadFreqs();

    connect(m_rt, &QComboBox::currentIndexChanged, this, &SpecPanel::onRespTypeChanged);
    connect(m_ft, &QComboBox::currentIndexChanged, this, [this](int idx) {
        m_method->blockSignals(true);
        m_method->clear();
        if (idx == 0)
            for (DesignMethod m : IIR_METHODS) m_method->addItem(method_name(m), int(m));
        else
            for (DesignMethod m : FIR_METHODS) m_method->addItem(method_name(m), int(m));
        m_method->setCurrentIndex(idx == 0 ? 3 : 0);
        m_method->blockSignals(false);
        m_manual.manual_from_zpk = false;  // a manual design starts from the current filter
        updateVisibility();
    });
    connect(m_method, &QComboBox::currentIndexChanged, this, [this] {
        updateVisibility();
        const auto m = DesignMethod(m_method->currentData().toInt());
        if (!m_quiet && (m == DesignMethod::ManualFIR || m == DesignMethod::ManualIIR)) emit manualSelected();
    });
    connect(m_min, &QRadioButton::toggled, this, &SpecPanel::updateVisibility);
    connect(m_window, &QComboBox::currentIndexChanged, this, [this](int) {
        const auto &w = window_list()[m_window->currentIndex()];
        if (w.par_name) {
            m_lwinpar->setText(QString(w.par_name) + ":");
            m_winpar->setText(fmt(w.par_default));
        }
        updateVisibility();
    });
    connect(m_unit, &QComboBox::currentIndexChanged, this, &SpecPanel::onUnitChanged);
    connect(m_amp_unit, &QComboBox::currentIndexChanged, this, [this](int u) {
        // convert the displayed values to the new unit
        const bool fir = is_fir(DesignMethod(m_method->currentData().toInt()));
        for (auto [e, pb] : {std::pair{m_apb, true}, std::pair{m_asb, false}}) {
            try {
                const double db = amp_to_db(parse(e, ""), AmpUnit(m_amp_unit_prev), fir, pb);
                e->setText(fmt(amp_from_db(db, AmpUnit(u), fir, pb)));
            } catch (const DesignError &) {
            }
        }
        m_amp_unit_prev = u;
        const QString un = m_amp_unit->currentText();
        m_lapb->setText("A_PB / " + un + ":");
        m_lasb->setText("A_SB / " + un + ":");
    });
    connect(m_fs, &QLineEdit::editingFinished, this, &SpecPanel::onFsEdited);
    connect(m_design, &QPushButton::clicked, this, &SpecPanel::designRequested);
    for (QLineEdit *e : findChildren<QLineEdit *>())
        connect(e, &QLineEdit::returnPressed, this, &SpecPanel::designRequested);
    updateVisibility();
}

void SpecPanel::addRow(QGridLayout *g, int row, const QString &label, QLineEdit *&edit, QLabel *&lbl,
                       const QString &tip) {
    lbl = new QLabel(label + ":", this);
    edit = new QLineEdit(this);
    edit->setToolTip(tip);
    lbl->setToolTip(tip);
    g->addWidget(lbl, row, 0);
    g->addWidget(edit, row, 1);
}

QString SpecPanel::fmt(double v) { return QString::number(v, 'g', 10); }

double SpecPanel::parse(const QLineEdit *e, const QString &name) {
    QString t = e->text().trimmed();
    if (t.contains(',') && !t.contains('.')) t.replace(',', '.');
    bool ok = false;
    const double v = t.toDouble(&ok);
    if (!ok || !std::isfinite(v)) throw DesignError(("Invalid value for " + name + ": '" + e->text() + "'").toStdString());
    return v;
}

void SpecPanel::storeFreqs() {
    Freqs &f = m_freqs[m_cur_rt];
    const double fs = m_fs_prev;
    auto get = [&](QLineEdit *e, double &v) {
        try {
            v = parse(e, "") / fs;
        } catch (const DesignError &) {
        }
    };
    get(m_fpb, f.f_pb);
    get(m_fpb2, f.f_pb2);
    get(m_fsb, f.f_sb);
    get(m_fsb2, f.f_sb2);
    get(m_fc, f.f_c);
    get(m_fc2, f.f_c2);
}

void SpecPanel::loadFreqs() {
    const Freqs &f = m_freqs[m_cur_rt];
    const double fs = m_fs_prev;
    m_fpb->setText(fmt(f.f_pb * fs));
    m_fpb2->setText(fmt(f.f_pb2 * fs));
    m_fsb->setText(fmt(f.f_sb * fs));
    m_fsb2->setText(fmt(f.f_sb2 * fs));
    m_fc->setText(fmt(f.f_c * fs));
    m_fc2->setText(fmt(f.f_c2 * fs));
}

void SpecPanel::onRespTypeChanged() {
    storeFreqs();
    m_cur_rt = m_rt->currentIndex();
    loadFreqs();
    updateVisibility();
}

void SpecPanel::onFsEdited() {
    double fs;
    try {
        fs = parse(m_fs, "f_S");
    } catch (const DesignError &) {
        m_fs->setText(fmt(m_fs_prev));
        return;
    }
    if (!(fs > 0) || fs == m_fs_prev) {
        m_fs->setText(fmt(m_fs_prev));
        return;
    }
    // keep the normalized frequencies, i.e. scale all frequencies with f_S
    storeFreqs();
    m_fs_prev = fs;
    loadFreqs();
    emit unitsChanged();
}

void SpecPanel::onUnitChanged() {
    const UnitDef &u = unit_def(m_unit->currentIndex());
    const bool norm = u.to_hz == 0;
    if (norm) {  // keep the normalized frequencies
        storeFreqs();
        m_fs_prev = u.norm_fs;
        m_fs->setText(fmt(u.norm_fs));
        loadFreqs();
    }
    m_fs->setEnabled(!norm);
    m_lfs->setText("f_S" + (norm ? QString() : " / " + m_unit->currentText()) + " =");
    emit unitsChanged();
}

double SpecPanel::unitToHz() const { return unit_def(m_unit->currentIndex()).to_hz; }

QString SpecPanel::freqLabel() const { return QString::fromUtf8(unit_def(m_unit->currentIndex()).f_label); }

QString SpecPanel::timeLabel() const { return QString::fromUtf8(unit_def(m_unit->currentIndex()).t_label); }

double SpecPanel::timeScale() const {
    // with f_S in kHz, n / f_S is in ms etc., normalized: f_S = 1 -> n
    return 1.0;
}

void SpecPanel::updateVisibility() {
    const auto rt = RespType(m_rt->currentIndex());
    const auto m = DesignMethod(m_method->currentData().toInt());
    const bool min = m_min->isChecked();
    const bool two = rt == RespType::BP || rt == RespType::BS;
    const bool fir = is_fir(m);
    const bool ma = m == DesignMethod::MovingAverage;
    const bool manual = is_manual(m);
    // moving average band pass / band stop filters only have a manual order
    const bool ma_man = ma && (two || !min);
    const bool equi_man = m == DesignMethod::Equiripple && !min;
    const bool edges = (min || equi_man) && !ma;
    m_min->setEnabled(!manual && !(ma && two));
    if (ma && two && min) m_man->setChecked(true);
    m_N->setEnabled(m == DesignMethod::Delay || (!manual && !min) || ma_man);
    m_man->setVisible(!manual);
    m_min->setVisible(!manual);
    m_lstages->setVisible(ma);
    m_stages->parentWidget()->setVisible(ma);
    m_grpF->setVisible(!manual && !ma_man);
    m_grpA->setVisible(!manual && !ma_man);
    m_N->setToolTip(ma ? tr("Number of delays M per stage")
                    : m == DesignMethod::Delay ? tr("Number of delays N")
                                               : tr("Filter order N (number of taps - 1 for FIR filters)"));
    m_hint->setVisible(manual || ma);
    if (m == DesignMethod::Delay) m_hint->setText(tr("<i>N</i> delays, H(z) = z<sup>-N</sup>."));
    else if (manual)
        m_hint->setText(tr("Enter the coefficients or poles / zeros in the <b>Coeffs</b> tab "
                           "(<i>Edit</i>) and press <i>Apply</i>, or drag poles / zeros in the <b>P / Z</b> "
                           "tab. DESIGN FILTER keeps them."));
    else if (ma)
        m_hint->setText(min && !two ? tr("Minimum number of delays M for the stop band specs F_SB, A_SB.")
                                    : tr("Order N = M delays per stage x stages."));
    if (ma) {  // only the stop band edge and attenuation are used
        for (auto *w : {static_cast<QWidget *>(m_fpb), static_cast<QWidget *>(m_lfpb),
                        static_cast<QWidget *>(m_fsb), static_cast<QWidget *>(m_lfsb),
                        static_cast<QWidget *>(m_fpb2), static_cast<QWidget *>(m_lfpb2),
                        static_cast<QWidget *>(m_fsb2), static_cast<QWidget *>(m_lfsb2),
                        static_cast<QWidget *>(m_fc), static_cast<QWidget *>(m_lfc),
                        static_cast<QWidget *>(m_fc2), static_cast<QWidget *>(m_lfc2),
                        static_cast<QWidget *>(m_apb), static_cast<QWidget *>(m_lapb),
                        static_cast<QWidget *>(m_wpb), static_cast<QWidget *>(m_lwpb),
                        static_cast<QWidget *>(m_wsb), static_cast<QWidget *>(m_lwsb),
                        static_cast<QWidget *>(m_window), static_cast<QWidget *>(m_lwindow),
                        static_cast<QWidget *>(m_winpar), static_cast<QWidget *>(m_lwinpar),
                        static_cast<QWidget *>(m_alg), static_cast<QWidget *>(m_lalg)})
            w->setVisible(false);
        m_fsb->setVisible(true);
        m_lfsb->setVisible(true);
        m_asb->setVisible(true);
        m_lasb->setVisible(true);
        return;
    }
    for (auto *w : {static_cast<QWidget *>(m_fpb), static_cast<QWidget *>(m_lfpb), static_cast<QWidget *>(m_fsb),
                    static_cast<QWidget *>(m_lfsb)})
        w->setVisible(edges);
    for (auto *w : {static_cast<QWidget *>(m_fpb2), static_cast<QWidget *>(m_lfpb2),
                    static_cast<QWidget *>(m_fsb2), static_cast<QWidget *>(m_lfsb2)})
        w->setVisible(edges && two);
    m_fc->setVisible(!edges);
    m_lfc->setVisible(!edges);
    m_fc2->setVisible(!edges && two);
    m_lfc2->setVisible(!edges && two);

    const bool apb = min || m == DesignMethod::Cheby1 || m == DesignMethod::Ellip;
    const bool asb = (min && m != DesignMethod::Bessel) || m == DesignMethod::Cheby2 || m == DesignMethod::Ellip;
    m_apb->setVisible(apb && !(fir && !min));
    m_lapb->setVisible(apb && !(fir && !min));
    m_asb->setVisible(asb && !(fir && !min));
    m_lasb->setVisible(asb && !(fir && !min));
    m_wpb->setVisible(equi_man);
    m_lwpb->setVisible(equi_man);
    m_wsb->setVisible(equi_man);
    m_lwsb->setVisible(equi_man);

    const bool firwin = m == DesignMethod::Firwin;
    m_window->setVisible(firwin);
    m_lwindow->setVisible(firwin);
    const auto &w = window_list()[std::max(0, m_window->currentIndex())];
    const bool par = firwin && w.par_name != nullptr;
    m_winpar->setVisible(par);
    m_lwinpar->setVisible(par);
    const bool alg = fir && min && !(firwin && w.type == WindowType::Kaiser);
    m_alg->setVisible(alg);
    m_lalg->setVisible(alg);
    // no amplitude specs (e.g. manual order Butterworth or FIR): hide the group with the unit
    if (!manual && !ma_man)
        m_grpA->setVisible(!m_apb->isHidden() || !m_asb->isHidden() || !m_wpb->isHidden());
}

FilterSpec SpecPanel::spec() const {
    FilterSpec s;
    s.rt = RespType(m_rt->currentIndex());
    s.method = DesignMethod(m_method->currentData().toInt());
    s.fo = m_min->isChecked() ? OrderMode::Min : OrderMode::Manual;
    s.N = m_N->value();
    s.f_s = parse(m_fs, "f_S");
    s.f_pb = parse(m_fpb, "F_PB");
    s.f_sb = parse(m_fsb, "F_SB");
    s.f_pb2 = parse(m_fpb2, "F_PB2");
    s.f_sb2 = parse(m_fsb2, "F_SB2");
    s.f_c = parse(m_fc, "F_C");
    s.f_c2 = parse(m_fc2, "F_C2");
    const AmpUnit au = AmpUnit(m_amp_unit->currentIndex());
    s.A_PB = amp_to_db(parse(m_apb, "A_PB"), au, is_fir(s.method), true);
    s.A_SB = amp_to_db(parse(m_asb, "A_SB"), au, is_fir(s.method), false);
    s.W_PB = parse(m_wpb, "W_PB");
    s.W_SB = parse(m_wsb, "W_SB");
    s.window = WindowType(m_window->currentData().toInt());
    if (m_winpar->isVisible() || window_list()[m_window->currentIndex()].par_name)
        s.win_par = parse(m_winpar, m_lwinpar->text());
    s.order_alg = RemezAlg(m_alg->currentIndex());
    s.ma_stages = m_stages->value();
    s.ma_norm = m_norm->isChecked();
    s.manual_ba = m_manual.manual_ba;
    s.manual_zpk = m_manual.manual_zpk;
    s.manual_from_zpk = m_manual.manual_from_zpk;
    return s;
}

void SpecPanel::updateFromDesign(const FilterSpec &s) {
    m_N->setValue(s.N);
    m_fc->setText(fmt(s.f_c));
    m_fc2->setText(fmt(s.f_c2));
    if (s.method == DesignMethod::Equiripple) {
        m_wpb->setText(fmt(s.W_PB));
        m_wsb->setText(fmt(s.W_SB));
    }
    if (s.method == DesignMethod::Firwin && s.window == WindowType::Kaiser) m_winpar->setText(fmt(s.win_par));
    if (is_manual(s.method)) {
        m_manual.manual_ba = s.manual_ba;
        m_manual.manual_zpk = s.manual_zpk;
        m_manual.manual_from_zpk = s.manual_from_zpk;
    }
    storeFreqs();
}

QString SpecPanel::unitKey() const { return unit_def(m_unit->currentIndex()).key; }

void SpecPanel::setSpec(const FilterSpec &s, const QString &unit) {
    // unit and f_S without rescaling the frequencies
    const int ui = std::max(0, m_unit->findData(unit));
    const UnitDef &u = unit_def(ui);
    const bool norm = u.to_hz == 0;
    m_unit->blockSignals(true);
    m_unit->setCurrentIndex(ui);
    m_unit->blockSignals(false);
    m_fs->setEnabled(!norm);
    m_lfs->setText("f_S" + (norm ? QString() : " / " + m_unit->currentText()) + " =");
    m_fs_prev = norm ? u.norm_fs : s.f_s;
    m_fs->setText(fmt(m_fs_prev));

    m_rt->blockSignals(true);
    m_rt->setCurrentIndex(int(s.rt));
    m_rt->blockSignals(false);
    m_cur_rt = int(s.rt);
    const double fs = s.f_s;
    m_freqs[m_cur_rt] = {s.f_pb / fs, s.f_pb2 / fs, s.f_sb / fs, s.f_sb2 / fs, s.f_c / fs, s.f_c2 / fs};
    loadFreqs();

    m_quiet = true;
    m_ft->setCurrentIndex(is_fir(s.method) ? 1 : 0);  // refills the method combo
    m_method->setCurrentIndex(std::max(0, m_method->findData(int(s.method))));
    m_quiet = false;
    (s.fo == OrderMode::Min ? m_min : m_man)->setChecked(true);
    m_N->setValue(s.N);
    m_window->setCurrentIndex(std::max(0, m_window->findData(int(s.window))));
    m_winpar->setText(fmt(s.win_par));
    m_alg->setCurrentIndex(int(s.order_alg));
    {
        const AmpUnit au = AmpUnit(m_amp_unit->currentIndex());
        m_apb->setText(fmt(amp_from_db(s.A_PB, au, is_fir(s.method), true)));
        m_asb->setText(fmt(amp_from_db(s.A_SB, au, is_fir(s.method), false)));
    }
    m_wpb->setText(fmt(s.W_PB));
    m_wsb->setText(fmt(s.W_SB));
    m_stages->setValue(s.ma_stages);
    m_norm->setChecked(s.ma_norm);
    m_manual.manual_ba = s.manual_ba;
    m_manual.manual_zpk = s.manual_zpk;
    m_manual.manual_from_zpk = s.manual_from_zpk;
    updateVisibility();
    emit unitsChanged();
}

void SpecPanel::setManual(const FilterSpec &manual, bool fir) {
    m_manual.manual_ba = manual.manual_ba;
    m_manual.manual_zpk = manual.manual_zpk;
    m_manual.manual_from_zpk = manual.manual_from_zpk;
    const auto m = fir ? DesignMethod::ManualFIR : DesignMethod::ManualIIR;
    m_quiet = true;
    if (m_ft->currentIndex() != (fir ? 1 : 0)) {
        const FilterSpec keep = m_manual;
        m_ft->setCurrentIndex(fir ? 1 : 0);  // refills the method combo
        m_manual = keep;
    }
    m_method->setCurrentIndex(std::max(0, m_method->findData(int(m))));
    m_quiet = false;
    updateVisibility();
}

void SpecPanel::setStatus(const QString &text, bool error) {
    m_status->setText(text);
    m_status->setStyleSheet(error ? "QLabel {color: #c00000; font-weight: bold}" : "");
}
