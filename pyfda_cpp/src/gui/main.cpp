// pyfda C++: C++ / Qt port of pyfda (Python Filter Design Analysis Tool)
//
// Command line options (mainly for automated tests):
//   --data <file>          load a data file into the "Data Filt" tab
//   --filter               filter the loaded data with the default design
//   --export <file.csv>    export original and filtered data
//   --screenshot <dir>     save a screenshot of every tab and quit
#include "data_filt_tab.hpp"
#include "main_window.hpp"

#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QTabWidget>
#include <QTimer>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QApplication::setApplicationName("pyfda_cpp");
    QApplication::setOrganizationName("pyfda");
    QApplication::setApplicationVersion("0.1.0");

    QCommandLineParser parser;
    parser.setApplicationDescription("C++ port of pyfda");
    parser.addHelpOption();
    parser.addVersionOption();
    QCommandLineOption optData("data", "Load a data file into the Data Filt tab.", "file");
    QCommandLineOption optFilter("filter", "Filter the loaded data.");
    QCommandLineOption optExport("export", "Export original and filtered data as CSV.", "file");
    QCommandLineOption optShot("screenshot", "Save screenshots of all tabs to <dir> and quit.", "dir");
    parser.addOptions({optData, optFilter, optExport, optShot});
    parser.process(app);

    MainWindow w;
    w.show();
    int rc = 0;
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
        }
        return rc;
    }
    if (parser.isSet(optExport)) return rc;
    return app.exec();
}
