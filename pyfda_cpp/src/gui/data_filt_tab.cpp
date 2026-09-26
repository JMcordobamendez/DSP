#include "data_filt_tab.hpp"

#include "filtering.hpp"
#include "logger.hpp"
#include "plot_widget.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSettings>
#include <QSplitter>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>

using namespace pyfda;

DataFiltView::DataFiltView(QWidget *parent) : DesignView(parent) {
    m_but_load = new QPushButton(tr("Load data ..."), this);
    m_but_load->setToolTip(tr("<span>Load data from a file (csv, txt, wav, npy). Delimiter, decimal "
                              "comma, header and encoding of csv / txt files are detected automatically, "
                              "comment lines (#, %, //) and metadata lines are skipped.</span>"));
    m_lbl_file = new QLabel(tr("No file loaded"), this);
    m_lbl_fs_warn = new QLabel(this);
    m_lbl_fs_warn->setStyleSheet("QLabel {color: darkorange; font-weight: bold}");
    m_lbl_fs_warn->setVisible(false);
    m_cmb_col = new QComboBox(this);
    m_cmb_col->setToolTip(tr("<span>Data column to be filtered.</span>"));
    m_cmb_col->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    m_cmb_time = new QComboBox(this);
    m_cmb_time->setToolTip(tr("<span>Column with the time axis, or '<i>n / f_S</i>' to calculate the "
                              "time axis from the sampling frequency of the current design.</span>"));
    m_cmb_time->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    m_but_filter = new QPushButton(tr("Filter data"), this);
    m_but_filter->setToolTip(tr("<span>Apply the current filter design to the selected data column.</span>"));
    m_but_filter->setEnabled(false);
    m_chk_zero_phase = new QCheckBox(tr("Zero phase"), this);
    m_chk_zero_phase->setToolTip(tr("<span>Filter forwards and backwards (<i>filtfilt</i>) for zero phase "
                                    "delay. The effective magnitude response is squared.</span>"));
    m_chk_spectrum = new QCheckBox(tr("Spectrum"), this);
    m_chk_spectrum->setToolTip(tr("<span>Show the magnitude spectra of original and filtered data "
                                  "in a second plot.</span>"));
    m_but_export = new QPushButton(tr("Export ..."), this);
    m_but_export->setToolTip(tr("<span>Save original and filtered data as a CSV file.</span>"));
    m_but_export->setEnabled(false);

    // two rows of controls so that the tab doesn't force a wide window
    auto *ctl = new QHBoxLayout();
    ctl->addWidget(m_but_load);
    ctl->addWidget(m_lbl_file);
    ctl->addSpacing(12);
    ctl->addWidget(new QLabel(tr("Column:"), this));
    ctl->addWidget(m_cmb_col);
    ctl->addWidget(new QLabel(tr("Time:"), this));
    ctl->addWidget(m_cmb_time);
    ctl->addStretch(10);
    auto *ctl2 = new QHBoxLayout();
    ctl2->addWidget(m_chk_zero_phase);
    ctl2->addWidget(m_but_filter);
    ctl2->addWidget(m_chk_spectrum);
    ctl2->addWidget(m_but_export);
    ctl2->addSpacing(12);
    ctl2->addWidget(m_lbl_fs_warn, 1);
    m_lbl_fs_warn->setWordWrap(true);

    auto *splitter = new QSplitter(Qt::Vertical, this);
    m_plot_t = new PlotWidget(splitter);
    m_plot_f = new PlotWidget(splitter);
    m_plot_f->setVisible(false);
    auto *lay = new QVBoxLayout(this);
    lay->setContentsMargins(4, 4, 4, 4);
    lay->addLayout(ctl);
    lay->addLayout(ctl2);
    lay->addWidget(splitter, 1);

    connect(m_but_load, &QPushButton::clicked, this, &DataFiltView::onLoad);
    connect(m_cmb_col, &QComboBox::currentIndexChanged, this, &DataFiltView::selectData);
    connect(m_cmb_time, &QComboBox::currentIndexChanged, this, [this] { redrawNow(); });
    connect(m_but_filter, &QPushButton::clicked, this, &DataFiltView::filterData);
    connect(m_chk_zero_phase, &QCheckBox::clicked, this, &DataFiltView::refilter);
    connect(m_chk_spectrum, &QCheckBox::clicked, this, [this] { redrawNow(); });
    connect(m_but_export, &QPushButton::clicked, this, &DataFiltView::exportData);
    redraw();
}

void DataFiltView::onLoad() {
    QSettings settings;
    const QString dir = settings.value("data_dir").toString();
    const QString fn = QFileDialog::getOpenFileName(
        this, tr("Load data to be filtered"), dir,
        tr("Data files (*.csv *.txt *.wav *.npy);;CSV (*.csv);;Text (*.txt);;Wave (*.wav);;"
           "Numpy (*.npy);;All files (*)"));
    if (fn.isEmpty()) return;
    settings.setValue("data_dir", QFileInfo(fn).absolutePath());
    loadFile(fn);
}

bool DataFiltView::loadFile(const QString &file_name) {
    DataTable t;
    try {
        // never let a malformed file crash the application
        t = load_data_file(QFile::encodeName(file_name).toStdString());
    } catch (const std::exception &e) {
        Logger::error(tr("Couldn't load '%1':\n%2").arg(file_name, e.what()));
        return false;
    }
    for (const auto &w : t.warnings) Logger::warning(QString::fromStdString(w));
    m_data = t;
    m_file_name = file_name;
    m_lbl_file->setText(QFileInfo(file_name).fileName());
    m_lbl_file->setToolTip(file_name);
    Logger::info(tr("Loaded %1 samples x %2 column(s) from '%3'.").arg(t.n_rows).arg(t.n_cols).arg(file_name));
    if (t.fs > 0) Logger::info(tr("Sampling rate stored in the file: %1 Hz.").arg(t.fs));

    QStringList names;
    for (const auto &n : t.names) names << QString::fromStdString(n);
    m_cmb_col->blockSignals(true);
    m_cmb_time->blockSignals(true);
    m_cmb_col->clear();
    m_cmb_col->addItems(names);
    m_cmb_time->clear();
    m_cmb_time->addItem("n / f_S");
    m_cmb_time->addItems(names);
    // When there are several columns and the first one is increasing,
    // assume it's a time column and preselect it
    const bool is_time = is_time_column(t, 0);
    if (is_time) m_cmb_time->setCurrentIndex(1);
    // Preselect the first data column that contains numbers
    for (size_t c = is_time ? 1 : 0; c < t.n_cols; ++c) {
        bool any = false;
        for (size_t r = 0; r < t.n_rows && !any; ++r) any = !std::isnan(t.at(r, c));
        if (any) {
            m_cmb_col->setCurrentIndex(int(c));
            break;
        }
    }
    m_cmb_col->blockSignals(false);
    m_cmb_time->blockSignals(false);
    selectData();
    return true;
}

void DataFiltView::selectData() {
    if (!m_data || m_cmb_col->currentIndex() < 0) return;
    m_x = m_data->column(size_t(m_cmb_col->currentIndex()));
    m_y.reset();
    m_but_filter->setEnabled(true);
    m_but_export->setEnabled(false);
    redrawNow();
}

void DataFiltView::refilter() {
    if (m_y) filterData();
    else redrawNow();
}

void DataFiltView::designChanged() {
    // the filter design has changed, filtered data is outdated
    if (m_y) computeFilter();
}

void DataFiltView::filterData() {
    computeFilter();
    m_plot_t->keepView(false);
    redrawNow();
}

bool DataFiltView::computeFilter() {
    if (m_x.empty()) return false;
    if (!m_design) {
        Logger::warning(tr("Design a filter first."));
        return false;
    }
    Vec x = m_x;
    size_t n_nan = 0;
    for (double &v : x)
        if (std::isnan(v)) {
            v = 0.0;
            ++n_nan;
        }
    if (n_nan) Logger::warning(tr("Replacing %1 non-numeric value(s) by zero.").arg(n_nan));
    const bool zero_phase = m_chk_zero_phase->isChecked();
    try {
        // second-order sections for IIR filters for better numerical stability
        if (!m_design->sos.empty())
            m_y = zero_phase ? sosfiltfilt(m_design->sos, x) : sosfilt(m_design->sos, x);
        else
            m_y = zero_phase ? filtfilt(m_design->ba.b, m_design->ba.a, x)
                             : lfilter(m_design->ba.b, m_design->ba.a, x);
    } catch (const std::exception &e) {
        Logger::error(tr("Filtering failed:\n%1").arg(e.what()));
        m_y.reset();
    }
    m_but_export->setEnabled(m_y.has_value());
    return m_y.has_value();
}

Vec DataFiltView::timeAxis(QString &label) const {
    const int idx = m_cmb_time->currentIndex();
    if (m_data && idx > 0) {
        label = m_cmb_time->currentText();
        return m_data->column(size_t(idx - 1));
    }
    label = m_ctx.t_label;
    Vec t(m_x.size());
    const bool norm = m_ctx.unit_to_hz == 0;
    for (size_t i = 0; i < t.size(); ++i) t[i] = norm ? double(i) : double(i) / m_ctx.f_s;
    return t;
}

void DataFiltView::checkFs() {
    // Estimate the sampling rate from the time column (assumed to be in seconds) and
    // warn when it doesn't match f_S of the design: the filter frequencies are relative to f_S
    const int idx = m_cmb_time->currentIndex();
    QString msg;
    double f_s_data = 0;
    bool non_uniform = false;
    if (m_data && idx > 0 && m_data->n_rows > 1) {
        Vec dt;
        const Vec t = m_data->column(size_t(idx - 1));
        for (size_t i = 1; i < t.size(); ++i) {
            const double d = t[i] - t[i - 1];
            if (std::isfinite(d)) dt.push_back(d);
        }
        if (!dt.empty()) {
            Vec s = dt;
            std::nth_element(s.begin(), s.begin() + long(s.size() / 2), s.end());
            double med = s[s.size() / 2];
            if (s.size() % 2 == 0) {
                const double lo = *std::max_element(s.begin(), s.begin() + long(s.size() / 2));
                med = 0.5 * (med + lo);
            }
            if (med > 0) f_s_data = 1.0 / med;
            const auto mm = std::minmax_element(dt.begin(), dt.end());
            non_uniform = *mm.second > 1.1 * *mm.first;
        }
    } else if (m_data && idx == 0 && m_data->fs > 0) {
        f_s_data = m_data->fs;  // sampling rate stored in a wav file
    }
    if (f_s_data > 0) {
        if (m_ctx.unit_to_hz == 0) {
            msg = tr("Data sampled at %1 Hz, set f_S = %1 Hz in the specs").arg(f_s_data, 0, 'g', 6);
        } else {
            const double f_s_design = m_ctx.f_s * m_ctx.unit_to_hz;
            if (std::fabs(f_s_design / f_s_data - 1) > 0.01)
                msg = tr("Data sampled at %1 Hz but f_S = %2 Hz").arg(f_s_data, 0, 'g', 6).arg(f_s_design, 0, 'g', 6);
        }
        if (!msg.isEmpty() && non_uniform) msg += tr(" (non-uniform sampling)");
    }
    if (!msg.isEmpty()) {
        m_lbl_fs_warn->setToolTip(tr("<span>The sampling rate is estimated from the time column "
                                     "(assumed to be in seconds) or read from the file. The filter "
                                     "frequencies are specified relative to f_S, so f_S must match the "
                                     "sampling rate of the data.</span>"));
        if (msg != m_lbl_fs_warn->text()) Logger::warning(msg);
    }
    m_lbl_fs_warn->setText(msg);
    m_lbl_fs_warn->setVisible(!msg.isEmpty());
}

void DataFiltView::redraw() {
    checkFs();
    const bool spec = m_chk_spectrum->isChecked();
    m_plot_f->setVisible(spec);
    m_plot_t->clear();
    m_plot_f->clear();
    if (m_x.empty()) {
        m_plot_t->setMessage(tr("Load data with 'Load data ...', then press 'Filter data'"));
        m_plot_t->autoscale();
        return;
    }
    QString t_label;
    const Vec t = timeAxis(t_label);
    const QString name = m_cmb_col->currentText();
    const QVector<double> qt(t.begin(), t.end());
    PlotWidget::Curve orig;
    orig.x = qt;
    orig.y = QVector<double>(m_x.begin(), m_x.end());
    orig.color = PlotWidget::color(0);
    orig.name = tr("Original") + " (" + name + ")";
    orig.alpha = m_y ? 0.6 : 1.0;
    orig.width = 1.2;
    m_plot_t->addCurve(orig);
    if (m_y) {
        PlotWidget::Curve filt;
        filt.x = qt;
        filt.y = QVector<double>(m_y->begin(), m_y->end());
        filt.color = PlotWidget::color(1);
        filt.name = tr("Filtered");
        filt.width = 1.2;
        m_plot_t->addCurve(filt);
    }
    m_plot_t->setXLabel(t_label);
    m_plot_t->setYLabel(name);
    m_plot_t->autoscale();
    m_plot_t->keepView(true);

    if (spec) {
        const size_t N = m_x.size();
        int k = 0;
        for (const Vec *d : {&m_x, m_y ? &*m_y : nullptr}) {
            if (!d) continue;
            Vec clean(*d);
            for (double &v : clean)
                if (!std::isfinite(v)) v = 0.0;
            const Vec A = amplitude_spectrum(clean);
            QVector<double> f(int(A.size())), y(int(A.size()));
            for (size_t i = 0; i < A.size(); ++i) {
                f[int(i)] = double(i) * m_ctx.f_s / double(N);
                y[int(i)] = 20 * std::log10(std::max(A[i], 1e-12));
            }
            m_plot_f->addCurve(f, y, PlotWidget::color(k), k == 0 ? tr("Original") : tr("Filtered"));
            ++k;
        }
        m_plot_f->setXLabel(m_ctx.f_label);
        m_plot_f->setYLabel("|X(f)| / dB");
        m_plot_f->autoscale();
    }
}

void DataFiltView::exportData() {
    if (!m_y) return;
    QSettings settings;
    const QString fn = QFileDialog::getSaveFileName(this, tr("Export filtered data"),
                                                    settings.value("data_dir").toString(), tr("CSV (*.csv)"));
    if (fn.isEmpty()) return;
    saveCsv(fn);
}

bool DataFiltView::saveCsv(const QString &file_name) {
    if (!m_y) return false;
    QString t_label;
    const Vec t = timeAxis(t_label);
    const std::string name = m_cmb_col->currentText().toStdString();
    try {
        write_csv(QFile::encodeName(file_name).toStdString(), {"t", name, name + "_filtered"}, {&t, &m_x, &*m_y});
        Logger::info(tr("Exported filtered data to '%1'.").arg(file_name));
        return true;
    } catch (const std::exception &e) {
        Logger::error(tr("Couldn't write '%1':\n%2").arg(file_name, e.what()));
        return false;
    }
}
