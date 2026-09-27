// Translations of the user interface (pyfda: pyfda_rc.py / translations). The Qt Linguist
// files translations/pyfda_cpp_<lang>.ts are compiled into the executable as resources and
// read directly, without lrelease; tools/update_ts.py keeps them in sync with the sources.
// Strings are looked up by their source text only (the tr() context is ignored).
#pragma once

#include <QHash>
#include <QString>
#include <QStringList>
#include <QTranslator>

class TsTranslator : public QTranslator {
public:
    using QTranslator::QTranslator;
    /// Read the finished translations of a .ts file, false if it can't be read
    bool loadTs(const QString &file);
    QString translate(const char *context, const char *source, const char *disambiguation = nullptr,
                      int n = -1) const override;
    bool isEmpty() const override { return m_tr.isEmpty(); }

private:
    QHash<QString, QString> m_tr;
};

namespace i18n {

/// Available languages, "en" is the language of the sources
QStringList languages();
/// Display name of a language in the language itself
QString languageName(const QString &lang);
/// "auto" -> the system language if it is available, else "en"
QString resolve(const QString &lang);
/// Install the translators (own strings and Qt's standard dialogs) for `lang`
/// ("auto", "en", "es", ...), before any widget is created; returns the language used
QString install(const QString &lang);

}  // namespace i18n
