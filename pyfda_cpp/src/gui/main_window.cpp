#include "main_window.hpp"

#include "coeffs_tab.hpp"
#include "data_filt_tab.hpp"
#include "logger.hpp"
#include "response_tabs.hpp"
#include "spec_panel.hpp"

#include <QApplication>
#include <QDateTime>
#include <QMenuBar>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QScrollArea>
#include <QSplitter>
#include <QTabWidget>

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
    add(new ImpulseView(this), "h[n]", tr("Impulse and step response"));
    add(new CoeffsView(this), tr("Coeffs"), tr("Filter coefficients"));
    m_data_filt = new DataFiltView(this);
    add(m_data_filt, tr("Data Filt"), tr("Filter data from a file with the current design"));

    m_log = new QPlainTextEdit(this);
    m_log->setReadOnly(true);
    m_log->setMaximumBlockCount(2000);

    auto *right = new QSplitter(Qt::Vertical, this);
    right->addWidget(m_tabs);
    right->addWidget(m_log);
    right->setStretchFactor(0, 5);
    right->setStretchFactor(1, 1);
    auto *main = new QSplitter(Qt::Horizontal, this);
    main->addWidget(scroll);
    main->addWidget(right);
    main->setStretchFactor(1, 1);
    setCentralWidget(main);

    auto *file = menuBar()->addMenu(tr("&File"));
    file->addAction(tr("&Design filter"), QKeySequence(Qt::CTRL | Qt::Key_D), this, [this] { design(); });
    file->addAction(tr("&Load data ..."), QKeySequence::Open, this, [this] {
        m_tabs->setCurrentWidget(m_data_filt);
        QMetaObject::invokeMethod(m_data_filt, "onLoad");
    });
    file->addSeparator();
    file->addAction(tr("&Quit"), QKeySequence::Quit, qApp, &QApplication::quit);
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
