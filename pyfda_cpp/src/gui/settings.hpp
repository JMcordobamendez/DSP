// User configuration (like pyfda's pyfda_user.conf), stored with QSettings as
// INI file <config dir>/pyfda/pyfda_cpp.ini (see main.cpp, --config-dir):
//   [session]  restore = true: reload the last design on start (pyfda_cpp_session.json
//              next to the INI file), geometry, splitters and current tab of the window
//   [dirs]     last directories of the file dialogs: filter, data, export
//   [csv]      delimiter = comma | semicolon | tab, decimal_comma = false: format of exported
//              CSV files (filtered data, transient data, coefficient table)
#pragma once

#include "data_io.hpp"

#include <QFileInfo>
#include <QSettings>
#include <QString>

namespace config {

inline QString iniFile() { return QSettings().fileName(); }

/// JSON file with the design of the last session
inline QString sessionFile() {
    const QFileInfo ini(iniFile());
    return ini.absolutePath() + "/" + ini.completeBaseName() + "_session.json";
}

inline bool restoreSession() { return QSettings().value("session/restore", true).toBool(); }

/// Last directory of the file dialogs of kind `key` (filter, data, export)
inline QString dir(const QString &key) { return QSettings().value("dirs/" + key).toString(); }

/// Remember the directory of `file` for the file dialogs of kind `key`
inline void setDir(const QString &key, const QString &file) {
    if (!file.isEmpty()) QSettings().setValue("dirs/" + key, QFileInfo(file).absolutePath());
}

inline pyfda::CsvFormat csvFormat() {
    QSettings s;
    const QString d = s.value("csv/delimiter", "comma").toString();
    pyfda::CsvFormat f;
    f.delimiter = d == "semicolon" ? ';' : d == "tab" ? '\t' : ',';
    f.decimal_comma = s.value("csv/decimal_comma", false).toBool();
    return f;
}

}  // namespace config
