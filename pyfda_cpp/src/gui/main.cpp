// pyfda C++: C++ / Qt port of pyfda (Python Filter Design Analysis Tool)
//
// Command line options (mainly for automated tests):
//   --load-filter <file>   load a filter file (JSON) and design it
//   --save-filter <file>   save the current design as filter file
//   --data <file>          load a data file into the "Data Filt" tab
//   --filter               filter the loaded data with the default design
//   --export <file.csv>    export original and filtered data
//   --stim <name>          stimulus of the transient analysis (dirac, step, sine, ...)
//   --formula <expr>       formula stimulus (numexpr syntax), e.g. "sin(2*pi*f1*n)"
//   --stim-file <file>     data file as stimulus
//   --tran-export <file>   export stimulus and response of the transient analysis
//   --fixpoint             fixpoint simulation in the transient analysis
//   --export-hdl <file>    export the fixpoint filter (.vhd, .v or .coe)
//   --hdl-testbench        with --export-hdl: also write a self-checking testbench <name>_tb.vhd / .v
//   --manual-ba <b/a>      enter coefficients in the Coeffs tab editor, e.g. "1,2,1/1,-0.5", and apply
//   --manual-zpk <z/p/k>   enter poles / zeros, e.g. "0.5:0.5,0.5:-0.5/0.9/2" (re:im pairs), and apply
//   --screenshot <dir>     save a screenshot of every tab and quit
//   --drag-pz <i:re:im>    move pole / zero i (zeros first, then poles) in the P/Z tab like with the mouse
//   --config-dir <dir>     directory of the configuration (default: the user's config directory),
//                          the INI file is <dir>/pyfda/pyfda_cpp.ini
//   --lang <auto|en|es>    language of the user interface (default: setting ui/language, auto)
//   --quit                 close the window right away like the user (saves the session)
#include "coeffs_tab.hpp"
#include "data_filt_tab.hpp"
#include "fixpoint_tab.hpp"
#include "main_window.hpp"
#include "response_tabs.hpp"
#include "settings.hpp"
#include "tran_tab.hpp"
#include "translator.hpp"
#include "window_viewer.hpp"

#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QSettings>
#include <QTabWidget>
#include <QTimer>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QApplication::setApplicationName("pyfda_cpp");
    QApplication::setOrganizationName("pyfda");
    QApplication::setApplicationVersion("0.1.0");
    // readable configuration file like pyfda's pyfda_user.conf (also on Windows instead of the registry)
    QSettings::setDefaultFormat(QSettings::IniFormat);

    QCommandLineParser parser;
    parser.setApplicationDescription("C++ port of pyfda");
    parser.addHelpOption();
    parser.addVersionOption();
    QCommandLineOption optData("data", "Load a data file into the Data Filt tab.", "file");
    QCommandLineOption optFilter("filter", "Filter the loaded data.");
    QCommandLineOption optExport("export", "Export original and filtered data as CSV.", "file");
    QCommandLineOption optShot("screenshot", "Save screenshots of all tabs to <dir> and quit.", "dir");
    QCommandLineOption optLoadFilt("load-filter", "Load a filter file (JSON).", "file");
    QCommandLineOption optSaveFilt("save-filter", "Save the design as filter file (JSON) and quit.", "file");
    QCommandLineOption optStim("stim", "Stimulus of the transient analysis.", "name");
    QCommandLineOption optFormula("formula", "Formula stimulus of the transient analysis.", "expr");
    QCommandLineOption optStimFile("stim-file", "Data file as stimulus of the transient analysis.", "file");
    QCommandLineOption optTranExp("tran-export", "Export the transient stimulus and response as CSV.", "file");
    QCommandLineOption optFix("fixpoint", "Fixpoint simulation in the transient analysis.");
    QCommandLineOption optHdl("export-hdl", "Export the fixpoint filter as VHDL (.vhd), Verilog (.v) or COE (.coe) and quit.",
                              "file");
    QCommandLineOption optTb("hdl-testbench", "With --export-hdl: also write a self-checking testbench.");
    QCommandLineOption optManBa("manual-ba", "Enter b / a in the coefficient editor and apply, e.g. 1,2,1/1,-0.5.", "b/a");
    QCommandLineOption optManZpk("manual-zpk", "Enter zeros / poles / gain in the P/Z editor and apply, "
                                               "e.g. 0.5:0.5,0.5:-0.5/0.9/2.", "z/p/k");
    QCommandLineOption optConf("config-dir", "Directory of the configuration file.", "dir");
    QCommandLineOption optDrag("drag-pz", "Move pole / zero i to re + j im in the P/Z tab, e.g. 0:0.5:0.6.", "i:re:im");
    QCommandLineOption optLang("lang", "Language of the user interface: auto, en or es.", "lang");
    QCommandLineOption optQuit("quit", "Close the main window right away (saves the session).");
    parser.addOptions({optConf, optLang, optQuit, optDrag, optData, optFilter, optExport, optShot, optLoadFilt, optSaveFilt, optStim, optTranExp, optFix,
                       optHdl, optTb, optManBa, optManZpk, optFormula, optStimFile});
    parser.process(app);
    if (parser.isSet(optConf))
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QDir(parser.value(optConf)).absolutePath());
    i18n::install(parser.isSet(optLang) ? parser.value(optLang) : config::language());

    MainWindow w;
    w.show();
    int rc = 0;
    if (!parser.isSet(optLoadFilt)) w.restoreSession();
    if (parser.isSet(optLoadFilt) && !w.openFilter(parser.value(optLoadFilt))) rc = 1;
    // manual coefficients / poles and zeros go through the editor of the Coeffs tab
    auto manual = [&](const QString &arg, bool zpk) {
        const QStringList parts = arg.split('/');
        CoeffsView *c = w.coeffs();
        w.tabs()->setCurrentWidget(c);
        c->setEditing(true, zpk ? 2 : 0);
        for (int col = 0; col < 2; ++col) {
            const QStringList vals = parts.value(col).split(',', Qt::SkipEmptyParts);
            for (int r = 0; r < 200; ++r) {  // clear the column, then fill it
                const QString v = vals.value(r);
                if (zpk) {
                    c->setCell(r, 2 * col, v.section(':', 0, 0));
                    c->setCell(r, 2 * col + 1, v.isEmpty() ? QString() : v.section(':', 1, 1).isEmpty() ? "0" : v.section(':', 1, 1));
                } else {
                    c->setCell(r, col, v);
                }
            }
        }
        try {
            pyfda::FilterSpec m = c->editedSpec();
            if (zpk) m.manual_zpk.k = parts.value(2, "1").toDouble();
            return w.designManual(m);
        } catch (const std::exception &) {
            return false;
        }
    };
    if (parser.isSet(optManBa) && !manual(parser.value(optManBa), false)) rc = 1;
    if (parser.isSet(optManZpk) && !manual(parser.value(optManZpk), true)) rc = 1;
    if (parser.isSet(optDrag)) {
        const QStringList v = parser.value(optDrag).split(':');
        w.tabs()->setCurrentWidget(w.poleZero());
        w.poleZero()->setDragEnabled(true);
        const QPointF pos(v.value(1).toDouble(), v.value(2).toDouble());
        w.poleZero()->dragTo(v.value(0).toInt(), pos * 0.5, false);  // half way, then to the end point
        w.poleZero()->dragTo(v.value(0).toInt(), pos, true);
    }
    if (parser.isSet(optFix)) w.fixpoint()->setSpec(w.fixpoint()->spec(), true);
    if (parser.isSet(optSaveFilt) && !w.saveFilter(parser.value(optSaveFilt))) rc = 1;
    if (parser.isSet(optHdl) && !w.fixpoint()->exportFile(parser.value(optHdl), parser.isSet(optTb))) rc = 1;
    if (parser.isSet(optStim)) {
        pyfda::StimParams p = w.tran()->params();
        bool found = false;
        for (const auto &i : pyfda::stim_list())
            if (parser.value(optStim) == i.key) {
                p.stim = i.stim;
                found = true;
            }
        if (found) w.tran()->setParams(p);
        else rc = 1;
    }
    if (parser.isSet(optFormula)) w.tran()->setFormula(parser.value(optFormula));
    if (parser.isSet(optStimFile) && !w.tran()->loadStimFile(parser.value(optStimFile))) rc = 1;
    if (parser.isSet(optTranExp) && !w.tran()->saveCsv(parser.value(optTranExp))) rc = 1;
    if (parser.isSet(optData) && !w.dataFilt()->loadFile(parser.value(optData))) rc = 1;
    if (parser.isSet(optFilter)) w.dataFilt()->filterData();
    if (parser.isSet(optExport) && !w.dataFilt()->saveCsv(parser.value(optExport))) rc = 1;
    if (parser.isSet(optShot)) {
        const QString dir = parser.value(optShot);
        QDir().mkpath(dir);
        for (int i = 0; i < w.tabs()->count(); ++i) {
            w.tabs()->setCurrentIndex(i);
            app.processEvents();
            w.grab().save(QString("%1/tab%2.png").arg(dir).arg(i));
            if (w.tabs()->currentWidget() == w.tran()) {  // frequency domain and spectrogram
                const char *names[] = {"", "freq", "spgr"};
                for (int k = 1; k < 3; ++k) {
                    w.tran()->innerTabs()->setCurrentIndex(k);
                    app.processEvents();
                    w.grab().save(QString("%1/tab%2_%3.png").arg(dir).arg(i).arg(names[k]));
                }
                w.tran()->innerTabs()->setCurrentIndex(0);
            }
        }
        WindowViewer wv(&w);  // window viewer with the default window
        wv.show();
        app.processEvents();
        wv.grab().save(QString("%1/window_viewer.png").arg(dir));
        return rc;
    }
    if (parser.isSet(optQuit)) {
        w.close();
        return rc;
    }
    if (parser.isSet(optExport) || parser.isSet(optSaveFilt) || parser.isSet(optTranExp) || parser.isSet(optHdl))
        return rc;
    return app.exec();
}
