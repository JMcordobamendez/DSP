#include "window_viewer.hpp"

#include "plot_widget.hpp"
#include "stimulus.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QSplitter>
#include <QSpinBox>
#include <QVBoxLayout>

#include <cmath>
#include <limits>

namespace {
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
}

using namespace pyfda;

WindowViewer::WindowViewer(QWidget *parent) : QDialog(parent) {
    setWindowTitle(tr("Window viewer"));
    auto *lay = new QVBoxLayout(this);
    auto *ctl = new QHBoxLayout();
    m_win = new QComboBox(this);
    for (const auto &w : window_list()) m_win->addItem(w.name, int(w.type));
    m_win->setCurrentIndex(m_win->findData(int(WindowType::Hann)));
    m_lpar = new QLabel(this);
    m_par = new QLineEdit(this);
    m_par->setMaximumWidth(70);
    m_lpar2 = new QLabel(this);
    m_par2 = new QLineEdit(this);
    m_par2->setMaximumWidth(70);
    m_N = new QSpinBox(this);
    m_N->setRange(2, 1 << 16);
    m_N->setValue(128);
    m_N->setToolTip(tr("Window length (periodic window for spectral analysis)"));
    m_log = new QCheckBox(tr("dB"), this);
    m_log->setChecked(true);
    ctl->addWidget(new QLabel(tr("Window:"), this));
    ctl->addWidget(m_win);
    ctl->addWidget(m_lpar);
    ctl->addWidget(m_par);
    ctl->addWidget(m_lpar2);
    ctl->addWidget(m_par2);
    ctl->addWidget(new QLabel("N =", this));
    ctl->addWidget(m_N);
    ctl->addWidget(m_log);
    ctl->addStretch(1);
    lay->addLayout(ctl);
    auto *split = new QSplitter(Qt::Horizontal, this);
    m_time = new PlotWidget(this);
    m_freq = new PlotWidget(this);
    split->addWidget(m_time);
    split->addWidget(m_freq);
    lay->addWidget(split, 1);
    m_props = new QLabel(this);
    m_props->setTextInteractionFlags(Qt::TextSelectableByMouse);
    lay->addWidget(m_props);
    resize(1000, 500);

    connect(m_win, &QComboBox::currentIndexChanged, this, [this] {
        const auto &w = window_list()[m_win->currentIndex()];
        if (w.par_name) m_par->setText(QString::number(w.par_default));
        if (w.par2_name) m_par2->setText(QString::number(w.par2_default));
        update();
    });
    connect(m_par, &QLineEdit::editingFinished, this, &WindowViewer::update);
    connect(m_par2, &QLineEdit::editingFinished, this, &WindowViewer::update);
    connect(m_N, &QSpinBox::valueChanged, this, &WindowViewer::update);
    connect(m_log, &QCheckBox::toggled, this, &WindowViewer::update);
    update();
}

void WindowViewer::setWindow(int index, int N, double par) {
    m_win->setCurrentIndex(index);
    m_N->setValue(N);
    m_par->setText(QString::number(par));
    update();
}

void WindowViewer::update() {
    const auto &wi = window_list()[std::max(0, m_win->currentIndex())];
    m_lpar->setVisible(wi.par_name != nullptr);
    m_par->setVisible(wi.par_name != nullptr);
    if (wi.par_name) m_lpar->setText(QString(wi.par_name) + ":");
    m_lpar2->setVisible(wi.par2_name != nullptr);
    m_par2->setVisible(wi.par2_name != nullptr);
    if (wi.par2_name) m_lpar2->setText(QString(wi.par2_name) + ":");
    bool ok = false;
    double par = m_par->text().toDouble(&ok);
    if (!ok) par = wi.par_default;
    double par2 = m_par2->text().toDouble(&ok);
    if (!ok) par2 = wi.par2_default;
    const int N = m_N->value();
    m_time->clear();
    m_freq->clear();
    Vec w;
    try {
        w = fft_window(wi.type, N, par, par2);
    } catch (const std::exception &e) {
        m_time->setMessage(e.what());
        m_time->autoscale();
        m_freq->autoscale();
        m_props->clear();
        return;
    }
    QVector<double> n(N), y(N);
    for (int i = 0; i < N; ++i) {
        n[i] = i;
        y[i] = w[size_t(i)];
    }
    m_time->addCurve(n, y, PlotWidget::color(0), QString(), N <= 64 ? PlotWidget::Style::Stem : PlotWidget::Style::Line);
    m_time->setTitle(tr("%1 window, N = %2").arg(wi.name).arg(N));
    m_time->setXLabel("n");
    m_time->setYLabel("w[n]");
    m_time->autoscale();

    const WindowProps p = window_props(w);
    const bool log = m_log->isChecked();
    QVector<double> f, a;
    for (size_t k = 0; k < p.F_bins.size() && p.F_bins[k] <= N / 2.0; ++k) {
        f << p.F_bins[k];
        a << (log ? p.W_db[k] : std::pow(10.0, p.W_db[k] / 20));
    }
    m_freq->addCurve(f, a, PlotWidget::color(1));
    m_freq->setTitle(tr("Magnitude spectrum (zero padded)"));
    m_freq->setXLabel(tr("f / bins (f_S / N)"));
    m_freq->setYLabel(log ? "|W(f)| / dB" : "|W(f)|");
    m_freq->setYLimits(log ? -160 : NaN, NaN);
    m_freq->setXLimits(0, std::min(N / 2.0, 20.0));
    m_freq->autoscale();
    m_props->setText(tr("Coherent gain %1 (%2 dB), NENBW %3 bins (%4 dB), scallop loss %5 dB, "
                        "3 dB bandwidth %6 bins, 6 dB bandwidth %7 bins, highest side lobe %8 dB")
                         .arg(p.cgain, 0, 'f', 4)
                         .arg(20 * std::log10(p.cgain), 0, 'f', 2)
                         .arg(p.nenbw, 0, 'f', 4)
                         .arg(10 * std::log10(p.nenbw), 0, 'f', 2)
                         .arg(p.scallop_db, 0, 'f', 2)
                         .arg(p.bw3_bins, 0, 'f', 3)
                         .arg(p.bw6_bins, 0, 'f', 3)
                         .arg(p.sidelobe_db, 0, 'f', 1));
}
