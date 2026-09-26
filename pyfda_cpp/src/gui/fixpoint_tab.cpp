#include "fixpoint_tab.hpp"

#include "conversions.hpp"
#include "hdl_export.hpp"
#include "logger.hpp"
#include "plot_widget.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QTableWidget>
#include <QVBoxLayout>

#include <cmath>
#include <limits>

using namespace pyfda;

namespace {
const Quant QUANTS[] = {Quant::Floor, Quant::Round, Quant::Fix, Quant::Ceil, Quant::None};
const Ovfl OVFLS[] = {Ovfl::Wrap, Ovfl::Sat, Ovfl::None};

// flattened quantized coefficients: FIR b or SOS (b0 b1 b2 1 a1 a2 per section)
Sos quantized_sos(const Sos &sos, const FxSpec &s, long long &n_over) {
    Quantizer Qb(s.qcb), Qa(s.qca);
    Sos q = sos;
    for (auto &sec : q) {
        for (int i : {0, 1, 2}) sec[size_t(i)] = Qb.fixp(sec[size_t(i)]);
        for (int i : {4, 5}) sec[size_t(i)] = Qa.fixp(sec[size_t(i)]);
    }
    n_over = Qb.overflows() + Qa.overflows();
    return q;
}
}  // namespace

FixpointView::FixpointView(QWidget *parent) : DesignView(parent) {
    auto *lay = new QVBoxLayout(this);
    lay->setContentsMargins(4, 4, 4, 4);

    auto *top = new QHBoxLayout();
    m_sim = new QCheckBox(tr("Fixpoint simulation in y[n]"), this);
    m_sim->setToolTip(tr("Calculate the transient response (tab y[n]) with the fixpoint filter"));
    m_structure = new QLabel(this);
    top->addWidget(m_sim);
    top->addSpacing(20);
    top->addWidget(m_structure);
    top->addStretch(1);
    m_coe = new QPushButton(tr("Export COE ..."), this);
    m_coe->setToolTip(tr("Quantized coefficients as Xilinx COE file (FIR Compiler), FIR filters only"));
    m_vhdl = new QPushButton(tr("Export VHDL ..."), this);
    m_vhdl->setToolTip(tr("<span>Synthesizable VHDL (numeric_std) of the fixpoint filter, bit exact with the "
                          "fixpoint simulation</span>"));
    top->addWidget(m_coe);
    top->addWidget(m_vhdl);
    lay->addLayout(top);

    auto *grp = new QGroupBox(tr("Word formats  (W = WI + WF + 1 bits, range -2^WI ... 2^WI - 2^-WF)"), this);
    auto *g = new QGridLayout(grp);
    const QStringList heads = {"", "WI", "WF", "W", tr("Quantization"), tr("Overflow")};
    for (int c = 0; c < heads.size(); ++c) g->addWidget(new QLabel("<b>" + heads[c] + "</b>", this), 0, c);
    auto add = [&](int row, const QString &label, const QString &tip) {
        Row r;
        r.label = new QLabel(label, this);
        r.label->setToolTip(tip);
        r.wi = new QSpinBox(this);
        r.wi->setRange(0, 32);
        r.wf = new QSpinBox(this);
        r.wf->setRange(0, 60);
        r.w = new QLabel(this);
        r.quant = new QComboBox(this);
        for (Quant q : QUANTS) r.quant->addItem(quant_key(q));
        r.quant->setToolTip(tr("floor: truncation, round: rounding (half to even), fix: towards zero, "
                               "ceil: towards +inf, none: no quantization"));
        r.ovfl = new QComboBox(this);
        for (Ovfl o : OVFLS) r.ovfl->addItem(ovfl_key(o));
        r.ovfl->setToolTip(tr("wrap: two's complement wrap around, sat: saturation, none: no overflow handling"));
        g->addWidget(r.label, row, 0);
        g->addWidget(r.wi, row, 1);
        g->addWidget(r.wf, row, 2);
        g->addWidget(r.w, row, 3);
        g->addWidget(r.quant, row, 4);
        g->addWidget(r.ovfl, row, 5);
        connect(r.wi, &QSpinBox::valueChanged, this, &FixpointView::onChanged);
        connect(r.wf, &QSpinBox::valueChanged, this, &FixpointView::onChanged);
        connect(r.quant, &QComboBox::currentIndexChanged, this, &FixpointView::onChanged);
        connect(r.ovfl, &QComboBox::currentIndexChanged, this, &FixpointView::onChanged);
        return r;
    };
    m_qi = add(1, tr("Input x"), tr("Format of the input signal"));
    m_qcb = add(2, tr("Coefficients b"), tr("Format of the transversal (numerator) coefficients"));
    m_qca = add(3, tr("Coefficients a"), tr("Format of the recursive (denominator) coefficients, a0 = 1"));
    m_qacc = add(4, tr("Accumulator"), tr("Format of the products and the accumulator"));
    m_qo = add(5, tr("Output y"), tr("Format of the output (and of the signals between IIR sections)"));
    m_auto_c = new QCheckBox(tr("auto WI"), this);
    m_auto_c->setToolTip(tr("Integer bits of the coefficients from their largest magnitude"));
    m_auto_acc = new QCheckBox(tr("auto"), this);
    m_auto_acc->setToolTip(tr("<span>Accumulator format from input and coefficient formats: no rounding of the "
                              "products and guard bits against overflows</span>"));
    g->addWidget(m_auto_c, 2, 6);
    g->addWidget(m_auto_acc, 4, 6);
    g->setColumnStretch(7, 1);
    connect(m_auto_c, &QCheckBox::toggled, this, &FixpointView::onChanged);
    connect(m_auto_acc, &QCheckBox::toggled, this, &FixpointView::onChanged);
    lay->addWidget(grp);
    m_ovfl_info = new QLabel(this);
    lay->addWidget(m_ovfl_info);

    auto *split = new QSplitter(Qt::Horizontal, this);
    auto *wt = new QWidget(this);
    auto *vt = new QVBoxLayout(wt);
    vt->setContentsMargins(0, 0, 0, 0);
    auto *ht = new QHBoxLayout();
    m_base = new QComboBox(this);
    m_base->addItems({tr("Hex"), tr("Bin"), tr("Oct"), tr("CSD")});
    m_base->setToolTip(tr("Number base for the integer representation (two's complement, CSD = canonical signed "
                          "digit)"));
    ht->addWidget(new QLabel(tr("Quantized coefficients:"), this));
    ht->addStretch(1);
    ht->addWidget(m_base);
    vt->addLayout(ht);
    m_table = new QTableWidget(this);
    m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    m_table->horizontalHeader()->setStretchLastSection(true);
    vt->addWidget(m_table, 1);
    split->addWidget(wt);
    m_plot = new PlotWidget(this);
    split->addWidget(m_plot);
    split->setStretchFactor(0, 2);
    split->setStretchFactor(1, 3);
    lay->addWidget(split, 1);

    connect(m_sim, &QCheckBox::toggled, this, [this] { emit fxChanged(); });
    connect(m_base, &QComboBox::currentIndexChanged, this, [this] { redrawNow(); });
    connect(m_coe, &QPushButton::clicked, this, [this] {
        QString fn = QFileDialog::getSaveFileName(this, tr("Export COE file"), QString(), tr("Xilinx COE (*.coe)"));
        if (fn.isEmpty()) return;
        if (QFileInfo(fn).suffix().isEmpty()) fn += ".coe";
        exportFile(fn);
    });
    connect(m_vhdl, &QPushButton::clicked, this, [this] {
        QString fn = QFileDialog::getSaveFileName(this, tr("Export VHDL"), QString(), tr("VHDL (*.vhd *.vhdl)"));
        if (fn.isEmpty()) return;
        if (QFileInfo(fn).suffix().isEmpty()) fn += ".vhd";
        exportFile(fn);
    });
    updateUi();
}

bool FixpointView::simulate() const { return m_sim->isChecked(); }

QFormat FixpointView::format(const Row &r) const {
    return QFormat{r.wi->value(), r.wf->value(), QUANTS[r.quant->currentIndex()], OVFLS[r.ovfl->currentIndex()]};
}

void FixpointView::setFormat(Row &r, const QFormat &q) {
    r.wi->setValue(q.WI);
    r.wf->setValue(q.WF);
    r.quant->setCurrentIndex(int(q.quant));
    r.ovfl->setCurrentIndex(int(q.ovfl));
    r.w->setText(QString::number(q.W()));
}

void FixpointView::setSpec(const FxSpec &s, bool simulate) {
    m_spec = s;
    if (m_design) update_auto_formats(m_spec, m_design->ba, m_design->sos, m_design->fir);
    updateUi();
    m_sim->blockSignals(true);
    m_sim->setChecked(simulate);
    m_sim->blockSignals(false);
    m_needs_redraw = true;
    if (isVisible()) redrawNow();
    emit fxChanged();
}

void FixpointView::updateUi() {
    m_updating = true;
    setFormat(m_qi, m_spec.qi);
    setFormat(m_qcb, m_spec.qcb);
    setFormat(m_qca, m_spec.qca);
    setFormat(m_qacc, m_spec.qacc);
    setFormat(m_qo, m_spec.qo);
    m_auto_c->setChecked(m_spec.coeff_auto);
    m_auto_acc->setChecked(m_spec.acc_auto);
    m_qcb.wi->setEnabled(!m_spec.coeff_auto);
    m_qca.wi->setEnabled(!m_spec.coeff_auto);
    m_qacc.wi->setEnabled(!m_spec.acc_auto);
    m_qacc.wf->setEnabled(!m_spec.acc_auto);
    const bool iir = m_design && !m_design->fir;
    for (QWidget *w : std::initializer_list<QWidget *>{m_qca.label, m_qca.wi, m_qca.wf, m_qca.w, m_qca.quant, m_qca.ovfl})
        w->setVisible(iir);
    m_coe->setEnabled(m_design && m_design->fir);
    m_vhdl->setEnabled(m_design != nullptr);
    if (!m_design) m_structure->clear();
    else if (m_design->fir)
        m_structure->setText(tr("Structure: FIR direct form, %1 taps").arg(m_design->ba.b.size()));
    else
        m_structure->setText(tr("Structure: cascade of %1 second-order sections (direct form 1)")
                                 .arg(m_design->sos.size()));
    m_updating = false;
}

void FixpointView::readUi() {
    m_spec.qi = format(m_qi);
    m_spec.qcb = format(m_qcb);
    m_spec.qca = format(m_qca);
    m_spec.qacc = format(m_qacc);
    m_spec.qo = format(m_qo);
    m_spec.coeff_auto = m_auto_c->isChecked();
    m_spec.acc_auto = m_auto_acc->isChecked();
}

void FixpointView::onChanged() {
    if (m_updating) return;
    readUi();
    if (m_design) update_auto_formats(m_spec, m_design->ba, m_design->sos, m_design->fir);
    updateUi();
    m_plot->keepView(true);
    redrawNow();
    emit fxChanged();
}

void FixpointView::designChanged() {
    if (m_design) update_auto_formats(m_spec, m_design->ba, m_design->sos, m_design->fir);
    updateUi();
    m_plot->keepView(false);
    emit fxChanged();
}

void FixpointView::redraw() {
    m_table->clear();
    m_table->setRowCount(0);
    m_plot->clear();
    m_ovfl_info->clear();
    if (!m_design) {
        m_plot->autoscale();
        return;
    }
    const int base_idx = m_base->currentIndex();
    auto int_str = [&](long long v, int W) {
        switch (base_idx) {
        case 0: return QString::fromStdString(to_base(v, W, 16));
        case 1: return QString::fromStdString(to_base(v, W, 2));
        case 2: return QString::fromStdString(to_base(v, W, 8));
        default: return QString::fromStdString(to_csd(v, W));
        }
    };
    auto num = [](double v) { return new QTableWidgetItem(QString::number(v, 'g', 12)); };
    long long n_over = 0;
    CVec H_q;
    const int n = 2048;
    Vec w(n);
    QVector<double> f(n);
    for (int i = 0; i < n; ++i) {
        w[size_t(i)] = PI * i / (n - 1);
        f[i] = m_ctx.f_s / 2 * i / (n - 1);
    }
    m_table->setHorizontalHeaderLabels({});
    if (m_design->fir) {
        const Vec &b = m_design->ba.b;
        const Vec bq = quant_coeffs(b, m_spec.qcb, false, &n_over);
        m_table->setColumnCount(5);
        m_table->setHorizontalHeaderLabels({"b", "b_Q", tr("error"), tr("integer"), m_base->currentText()});
        m_table->setRowCount(int(b.size()));
        for (int i = 0; i < int(b.size()); ++i) {
            const long long v = to_int(bq[size_t(i)], m_spec.qcb);
            m_table->setItem(i, 0, num(b[size_t(i)]));
            m_table->setItem(i, 1, num(bq[size_t(i)]));
            m_table->setItem(i, 2, num(bq[size_t(i)] - b[size_t(i)]));
            m_table->setItem(i, 3, new QTableWidgetItem(QString::number(v)));
            m_table->setItem(i, 4, new QTableWidgetItem(int_str(v, m_spec.qcb.W())));
        }
        H_q = freqz(Ba{bq, {1.0}}, w);
    } else {
        const Sos sq = quantized_sos(m_design->sos, m_spec, n_over);
        m_table->setColumnCount(6);
        m_table->setHorizontalHeaderLabels({tr("section"), tr("coeff."), tr("value"), "Q", tr("integer"),
                                            m_base->currentText()});
        static const char *names[] = {"b0", "b1", "b2", "a0", "a1", "a2"};
        int row = 0;
        m_table->setRowCount(int(sq.size()) * 5);
        for (size_t k = 0; k < sq.size(); ++k)
            for (int i : {0, 1, 2, 4, 5}) {
                const QFormat &q = i < 3 ? m_spec.qcb : m_spec.qca;
                const long long v = to_int(sq[k][size_t(i)], q);
                m_table->setItem(row, 0, new QTableWidgetItem(QString::number(k + 1)));
                m_table->setItem(row, 1, new QTableWidgetItem(names[i]));
                m_table->setItem(row, 2, num(m_design->sos[k][size_t(i)]));
                m_table->setItem(row, 3, num(sq[k][size_t(i)]));
                m_table->setItem(row, 4, new QTableWidgetItem(QString::number(v)));
                m_table->setItem(row, 5, new QTableWidgetItem(int_str(v, q.W())));
                ++row;
            }
        H_q = freqz(sq, w);
    }
    const CVec H = m_design->sos.empty() ? freqz(m_design->ba, w) : freqz(m_design->sos, w);
    QVector<double> a(n), aq(n);
    for (int i = 0; i < n; ++i) {
        a[i] = 20 * std::log10(std::max(std::abs(H[size_t(i)]), 1e-15));
        aq[i] = 20 * std::log10(std::max(std::abs(H_q[size_t(i)]), 1e-15));
    }
    m_plot->addCurve(f, a, PlotWidget::color(0), tr("ideal"));
    PlotWidget::Curve c;
    c.x = f;
    c.y = aq;
    c.color = PlotWidget::color(1);
    c.name = tr("quantized coefficients");
    m_plot->addCurve(c);
    m_plot->setTitle(tr("Magnitude response with quantized coefficients"));
    m_plot->setXLabel(m_ctx.f_label);
    m_plot->setYLabel("|H(f)| / dB");
    m_plot->setYLimits(-150, std::numeric_limits<double>::quiet_NaN());
    m_plot->autoscale();
    m_ovfl_info->setText(n_over ? tr("<span style='color:#c00000'><b>%1 coefficient overflow(s)</b>, increase WI of "
                                     "the coefficients.</span>")
                                      .arg(n_over)
                                : tr("No coefficient overflows. Signal overflows are shown in the y[n] tab."));
}

bool FixpointView::exportFile(const QString &file_name) {
    if (!m_design) {
        Logger::error(tr("No filter designed yet."));
        return false;
    }
    std::string text;
    try {
        const QString suffix = QFileInfo(file_name).suffix().toLower();
        const std::string name = QFileInfo(file_name).completeBaseName().toStdString();
        if (suffix == "coe") {
            if (!m_design->fir) throw DesignError("COE files are only available for FIR filters.");
            text = export_coe(m_design->ba.b, m_spec.qcb, 16);
        } else {
            text = m_design->fir ? export_vhdl_fir(m_design->ba.b, m_spec, name)
                                 : export_vhdl_sos(m_design->sos, m_spec, name);
        }
    } catch (const std::exception &e) {
        Logger::error(e.what());
        return false;
    }
    QFile f(file_name);
    if (!f.open(QIODevice::WriteOnly) || f.write(text.data(), qint64(text.size())) != qint64(text.size())) {
        Logger::error(tr("Couldn't write '%1'.").arg(file_name));
        return false;
    }
    Logger::info(tr("Exported the fixpoint filter to '%1'.").arg(file_name));
    return true;
}
