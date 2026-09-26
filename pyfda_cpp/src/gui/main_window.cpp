#include "main_window.hpp"

#include "coeffs_tab.hpp"
#include "data_filt_tab.hpp"
#include "fixpoint_tab.hpp"
#include "info_tab.hpp"
#include "plot_3d.hpp"
#include "window_viewer.hpp"
#include "filter_io.hpp"
#include "logger.hpp"
#include "response_tabs.hpp"
#include "spec_panel.hpp"
#include "settings.hpp"
#include "tran_tab.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QLabel>
#include <QDateTime>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QDir>
#include <QSettings>
#include <QMenuBar>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QScrollArea>
#include <QSplitter>
#include <QTabWidget>

#include <cmath>

using namespace pyfda;

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    setWindowTitle("pyfda C++ - Python Filter Design Analysis Tool");

    m_specs = new SpecPanel(this);
    auto *scroll = new QScrollArea(this);
    scroll->setWidget(m_specs);
    scroll->setWidgetResizable(true);
    scroll->setMinimumWidth(280);

    m_tabs = new QTabWidget(this);
    auto add = [&](DesignView *v, const QString &label, const QString &tip) {
        m_views << v;
        const int i = m_tabs->addTab(v, label);
        m_tabs->setTabToolTip(i, tip);
    };
    add(new MagnitudeView(this), "|H(f)|", tr("Magnitude response"));
    add(new PhaseView(this), QString::fromUtf8("φ(f)"), tr("Phase response"));
    add(new GroupDelayView(this), QString::fromUtf8("τ(f)"), tr("Group delay"));
    add(new PoleZeroView(this), "P / Z", tr("Pole / zero plot"));
    add(new ThreeDView(this), "3D", tr("3D magnitude response |H(z)|"));
    add(new ImpulseView(this), "h[n]", tr("Impulse and step response"));
    m_tran = new TranView(this);
    add(m_tran, "y[n]", tr("Transient analysis: stimulus and response in the time and frequency domain"));
    m_coeffs = new CoeffsView(this);
    add(m_coeffs, tr("Coeffs"), tr("Filter coefficients"));
    m_fix = new FixpointView(this);
    add(m_fix, tr("Fixpoint"), tr("Fixpoint formats, quantized coefficients, COE / VHDL export"));
    connect(m_fix, &FixpointView::fxChanged, this, [this] { m_tran->setFixpoint(m_fix->spec(), m_fix->simulate()); });
    add(new InfoView(this), tr("Info"), tr("Filter properties and specifications vs. achieved values"));
    m_data_filt = new DataFiltView(this);
    add(m_data_filt, tr("Data Filt"), tr("Filter data from a file with the current design"));

    m_log = new QPlainTextEdit(this);
    m_log->setReadOnly(true);
    m_log->setMaximumBlockCount(2000);

    auto *right = m_split_right = new QSplitter(Qt::Vertical, this);
    right->addWidget(m_tabs);
    right->addWidget(m_log);
    right->setStretchFactor(0, 5);
    right->setStretchFactor(1, 1);
    auto *main = m_split_main = new QSplitter(Qt::Horizontal, this);
    main->addWidget(scroll);
    main->addWidget(right);
    main->setStretchFactor(1, 1);
    setCentralWidget(main);

    auto *file = menuBar()->addMenu(tr("&File"));
    file->addAction(tr("&Design filter"), QKeySequence(Qt::CTRL | Qt::Key_D), this, [this] { design(); });
    file->addAction(tr("&Open filter ..."), QKeySequence::Open, this, [this] {
        const QString fn = QFileDialog::getOpenFileName(this, tr("Open filter"), config::dir("filter"),
                                                        tr("pyfda_cpp filter (*.json);;All files (*)"));
        if (!fn.isEmpty()) openFilter(fn);
    });
    file->addAction(tr("&Save filter ..."), QKeySequence::Save, this, [this] {
        QString fn = QFileDialog::getSaveFileName(this, tr("Save filter"), config::dir("filter"),
                                                  tr("pyfda_cpp filter (*.json)"));
        if (fn.isEmpty()) return;
        if (QFileInfo(fn).suffix().isEmpty()) fn += ".json";
        saveFilter(fn);
    });
    file->addAction(tr("&Export coefficients ..."), QKeySequence(Qt::CTRL | Qt::Key_E), this,
                    [this] { m_coeffs->exportDialog(); });
    file->addSeparator();
    file->addAction(tr("&Load data ..."), QKeySequence(Qt::CTRL | Qt::Key_L), this, [this] {
        m_tabs->setCurrentWidget(m_data_filt);
        QMetaObject::invokeMethod(m_data_filt, "onLoad");
    });
    file->addSeparator();
    file->addAction(tr("&Preferences ..."), this, &MainWindow::preferences);
    file->addSeparator();
    file->addAction(tr("&Quit"), QKeySequence::Quit, this, &QWidget::close);
    auto *tools = menuBar()->addMenu(tr("&Tools"));
    tools->addAction(tr("&Window viewer ..."), this, [this] {
        auto *v = new WindowViewer(this);
        v->setAttribute(Qt::WA_DeleteOnClose);
        v->show();
    });
    auto *help = menuBar()->addMenu(tr("&Help"));
    help->addAction(tr("&About"), this, [this] {
        QMessageBox::about(this, tr("About pyfda C++"),
                           tr("<b>pyfda C++</b> %1<br><br>C++ / Qt port of <a href='https://github.com/chipmuenk/pyfda'>"
                              "pyfda</a> (Python Filter Design Analysis Tool, MIT license) with the "
                              "'Data Filt' tab.<br><br>Filter design routines ported from scipy.signal "
                              "(BSD license).").arg(QApplication::applicationVersion()));
    });

    connect(&Logger::instance(), &Logger::message, this, &MainWindow::onLog);
    connect(m_specs, &SpecPanel::designRequested, this, [this] { design(); });
    connect(m_specs, &SpecPanel::unitsChanged, this, &MainWindow::updateViews);
    connect(m_coeffs, &CoeffsView::manualDesignRequested, this, &MainWindow::designManual);
    connect(m_specs, &SpecPanel::manualSelected, this, [this] {
        // like pyfda, a manual filter starts with the coefficients of the current design
        if (!m_design) return;
        FilterSpec m;
        m.manual_ba = m_design->ba;
        m.manual_from_zpk = false;
        m_specs->setManual(m, m_design->fir);
        design();
    });
    resize(1280, 820);
    design();
}

ViewContext MainWindow::context() const {
    ViewContext c;
    try {
        c.f_s = m_specs->spec().f_s;
    } catch (const DesignError &) {
    }
    c.unit_to_hz = m_specs->unitToHz();
    c.f_label = m_specs->freqLabel();
    c.t_label = m_specs->timeLabel();
    return c;
}

bool MainWindow::design() {
    try {
        const FilterSpec s = m_specs->spec();
        auto d = std::make_unique<FilterDesign>(design_filter(s));
        m_specs->updateFromDesign(d->spec);
        m_design = std::move(d);
        m_specs->setStatus(QString::fromStdString(m_design->info), false);
        Logger::info(tr("Designed %1.").arg(QString::fromStdString(m_design->info)));
        updateViews();
        return true;
    } catch (const std::exception &e) {
        m_specs->setStatus(e.what(), true);
        Logger::error(e.what());
        return false;
    }
}

bool MainWindow::designManual(const FilterSpec &manual) {
    FilterSpec s = manual;
    s.method = DesignMethod::ManualIIR;  // FIR / IIR is decided from the coefficients
    try {
        const FilterDesign d = design_filter(s);
        m_specs->setManual(manual, d.fir);
    } catch (const std::exception &e) {
        m_specs->setStatus(e.what(), true);
        Logger::error(e.what());
        return false;
    }
    return design();
}

bool MainWindow::openFilter(const QString &file_name, bool session) {
    QFile f(file_name);
    if (!f.open(QIODevice::ReadOnly)) {
        Logger::error(tr("Couldn't open '%1'.").arg(file_name));
        return false;
    }
    const QByteArray data = f.readAll();
    FilterFile ff;
    try {
        ff = filter_from_json(std::string(data.constData(), size_t(data.size())));
    } catch (const std::exception &e) {
        Logger::error(tr("'%1': %2").arg(QFileInfo(file_name).fileName(), e.what()));
        return false;
    }
    if (!session) config::setDir("filter", file_name);
    m_specs->setSpec(ff.spec, QString::fromStdString(ff.unit));
    Logger::info(session ? tr("Restored the design of the last session.") : tr("Loaded filter '%1'.").arg(file_name));
    if (!design()) return false;
    if (ff.has_fx) m_fix->setSpec(ff.fx, ff.fx_sim);
    // the stored coefficients are only for reference, warn if the design differs
    // (e.g. file edited by hand or written by a different version)
    auto differs = [](const Vec &a, const Vec &b) {
        if (a.empty()) return false;  // not stored
        if (a.size() != b.size()) return true;
        double scale = 0, err = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            scale = std::max(scale, std::fabs(b[i]));
            err = std::max(err, std::fabs(a[i] - b[i]));
        }
        return err > 1e-9 * std::max(scale, 1e-300);
    };
    if (differs(ff.ba.b, m_design->ba.b) || differs(ff.ba.a, m_design->ba.a))
        Logger::warning(tr("The coefficients stored in the file differ from the redesigned filter, "
                           "the filter was designed from the specifications in the file."));
    return true;
}

bool MainWindow::saveFilter(const QString &file_name) {
    if (!m_design) {
        Logger::error(tr("No filter designed yet."));
        return false;
    }
    const std::string text = filter_to_json(*m_design, m_specs->unitKey().toStdString(), &m_fix->spec(), m_fix->simulate());
    QFile f(file_name);
    if (!f.open(QIODevice::WriteOnly) || f.write(text.data(), qint64(text.size())) != qint64(text.size())) {
        Logger::error(tr("Couldn't write '%1'.").arg(file_name));
        return false;
    }
    config::setDir("filter", file_name);
    Logger::info(tr("Saved filter to '%1'.").arg(file_name));
    return true;
}

void MainWindow::updateViews() {
    const ViewContext c = context();
    for (DesignView *v : m_views) v->setDesign(m_design.get(), c);
}

void MainWindow::onLog(int level, const QString &text) {
    static const char *names[] = {"INFO", "WARNING", "ERROR"};
    static const char *colors[] = {"", "darkorange", "#c00000"};
    const QString time = QDateTime::currentDateTime().toString("HH:mm:ss");
    QString html = QString("[%1] <b>%2</b>: %3").arg(time, names[level], text.toHtmlEscaped().replace("\n", "<br>"));
    if (level > 0) html = QString("<span style='color:%1'>%2</span>").arg(colors[level], html);
    m_log->appendHtml(html);
}

void MainWindow::restoreSession() {
    QSettings s;
    if (s.contains("window/geometry")) restoreGeometry(s.value("window/geometry").toByteArray());
    m_split_main->restoreState(s.value("window/split_main").toByteArray());
    m_split_right->restoreState(s.value("window/split_right").toByteArray());
    if (config::restoreSession() && QFileInfo::exists(config::sessionFile())) openFilter(config::sessionFile(), true);
    m_tabs->setCurrentIndex(qBound(0, s.value("window/tab", 0).toInt(), m_tabs->count() - 1));
}

void MainWindow::saveSession() {
    QSettings s;
    s.setValue("window/geometry", saveGeometry());
    s.setValue("window/split_main", m_split_main->saveState());
    s.setValue("window/split_right", m_split_right->saveState());
    s.setValue("window/tab", m_tabs->currentIndex());
    s.sync();  // creates the directory of the INI file
    if (m_design && config::restoreSession()) {
        const std::string text =
            filter_to_json(*m_design, m_specs->unitKey().toStdString(), &m_fix->spec(), m_fix->simulate());
        QFile f(config::sessionFile());
        if (!f.open(QIODevice::WriteOnly) || f.write(text.data(), qint64(text.size())) != qint64(text.size()))
            Logger::error(tr("Couldn't write '%1'.").arg(config::sessionFile()));
    }
}

void MainWindow::closeEvent(QCloseEvent *e) {
    saveSession();
    QApplication::closeAllWindows();  // e.g. window viewers
    e->accept();
}

void MainWindow::preferences() {
    QDialog dlg(this);
    dlg.setWindowTitle(tr("Preferences"));
    auto *form = new QFormLayout(&dlg);
    auto *restore = new QCheckBox(tr("Restore the design and the window of the last session"), &dlg);
    restore->setChecked(config::restoreSession());
    form->addRow(tr("Start:"), restore);
    auto *delim = new QComboBox(&dlg);
    delim->addItem(tr("Comma ,"), "comma");
    delim->addItem(tr("Semicolon ;"), "semicolon");
    delim->addItem(tr("Tab"), "tab");
    QSettings s;
    delim->setCurrentIndex(qMax(0, delim->findData(s.value("csv/delimiter", "comma"))));
    form->addRow(tr("CSV delimiter:"), delim);
    auto *comma = new QCheckBox(tr("Decimal comma (e.g. 0,5)"), &dlg);
    comma->setToolTip(tr("For spreadsheets with a comma as decimal separator, needs ';' or tab as delimiter"));
    comma->setChecked(s.value("csv/decimal_comma", false).toBool());
    auto update = [=] { comma->setEnabled(delim->currentData() != "comma"); };
    connect(delim, &QComboBox::currentIndexChanged, &dlg, update);
    update();
    form->addRow(tr("CSV numbers:"), comma);
    auto *file = new QLabel(QDir::toNativeSeparators(config::iniFile()), &dlg);
    file->setTextInteractionFlags(Qt::TextSelectableByMouse);
    form->addRow(tr("Configuration file:"), file);
    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);
    form->addRow(buttons);
    if (dlg.exec() != QDialog::Accepted) return;
    s.setValue("session/restore", restore->isChecked());
    s.setValue("csv/delimiter", delim->currentData());
    s.setValue("csv/decimal_comma", comma->isChecked());
}
