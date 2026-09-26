#include "spec_panel.hpp"

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
const DesignMethod IIR_METHODS[] = {DesignMethod::Butter, DesignMethod::Cheby1, DesignMethod::Cheby2,
                                    DesignMethod::Ellip, DesignMethod::Bessel};
const DesignMethod FIR_METHODS[] = {DesignMethod::Equiripple, DesignMethod::Firwin};
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

    auto *wOrder = new QWidget(this);
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

    m_lwindow = new QLabel(tr("Window:"), this);
    m_window = new QComboBox(this);
    for (const auto &w : window_list()) m_window->addItem(w.name, int(w.type));
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
    auto *grpF = new QGroupBox(tr("Frequencies"), this);
    auto *gF = new QGridLayout(grpF);
    m_unit = new QComboBox(this);
    m_unit->addItems({"f_S (norm.)", "Hz", "kHz", "MHz"});
    m_unit->setToolTip(tr("<span>Frequency unit. With <i>f_S (norm.)</i> all frequencies are "
                          "normalized to the sampling frequency f_S = 1.</span>"));
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
    auto *grpA = new QGroupBox(tr("Amplitudes"), this);
    auto *gA = new QGridLayout(grpA);
    addRow(gA, 0, "A_PB / dB", m_apb, m_lapb, tr("Maximum pass band ripple in dB"));
    addRow(gA, 1, "A_SB / dB", m_asb, m_lasb, tr("Minimum stop band attenuation in dB"));
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
        updateVisibility();
    });
    connect(m_method, &QComboBox::currentIndexChanged, this, &SpecPanel::updateVisibility);
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
    const bool norm = m_unit->currentIndex() == 0;
    if (norm) {
        storeFreqs();
        m_fs_prev = 1.0;
        m_fs->setText("1");
        loadFreqs();
    }
    m_fs->setEnabled(!norm);
    const QString u = norm ? QString() : " / " + m_unit->currentText();
    m_lfs->setText("f_S" + u + " =");
    emit unitsChanged();
}

double SpecPanel::unitToHz() const {
    switch (m_unit->currentIndex()) {
    case 1: return 1.0;
    case 2: return 1e3;
    case 3: return 1e6;
    default: return 0.0;
    }
}

QString SpecPanel::freqLabel() const {
    return m_unit->currentIndex() == 0 ? QString("F = f / f_S") : "f / " + m_unit->currentText();
}

QString SpecPanel::timeLabel() const {
    switch (m_unit->currentIndex()) {
    case 1: return "t / s";
    case 2: return "t / ms";
    case 3: return QString("t / %1s").arg(QChar(0x00B5));
    default: return "n";
    }
}

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
    const bool equi_man = m == DesignMethod::Equiripple && !min;
    const bool edges = min || equi_man;
    m_N->setEnabled(!min);
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
    s.A_PB = parse(m_apb, "A_PB");
    s.A_SB = parse(m_asb, "A_SB");
    s.W_PB = parse(m_wpb, "W_PB");
    s.W_SB = parse(m_wsb, "W_SB");
    s.window = WindowType(m_window->currentData().toInt());
    if (m_winpar->isVisible() || window_list()[m_window->currentIndex()].par_name)
        s.win_par = parse(m_winpar, m_lwinpar->text());
    s.order_alg = RemezAlg(m_alg->currentIndex());
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
    storeFreqs();
}

void SpecPanel::setStatus(const QString &text, bool error) {
    m_status->setText(text);
    m_status->setStyleSheet(error ? "QLabel {color: #c00000; font-weight: bold}" : "");
}
