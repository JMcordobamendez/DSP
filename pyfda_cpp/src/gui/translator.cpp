#include "translator.hpp"

#include <QCoreApplication>
#include <QFile>
#include <QLibraryInfo>
#include <QLocale>
#include <QXmlStreamReader>

bool TsTranslator::loadTs(const QString &file) {
    QFile f(file);
    if (!f.open(QIODevice::ReadOnly)) return false;
    QXmlStreamReader xml(&f);
    QString source, translation;
    bool finished = false;
    while (!xml.atEnd()) {
        xml.readNext();
        if (xml.isStartElement()) {
            if (xml.name() == QLatin1String("message")) {
                source.clear();
                translation.clear();
                finished = false;
            } else if (xml.name() == QLatin1String("source")) {
                source = xml.readElementText();
            } else if (xml.name() == QLatin1String("translation")) {
                const auto type = xml.attributes().value(QLatin1String("type"));
                finished = type != QLatin1String("unfinished") && type != QLatin1String("obsolete");
                translation = xml.readElementText();
            }
        } else if (xml.isEndElement() && xml.name() == QLatin1String("message")) {
            if (finished && !source.isEmpty() && !translation.isEmpty()) m_tr.insert(source, translation);
        }
    }
    return !xml.hasError();
}

QString TsTranslator::translate(const char *, const char *source, const char *, int) const {
    if (!source) return {};
    return m_tr.value(QString::fromUtf8(source));
}

namespace i18n {

QStringList languages() { return {"en", "es"}; }

QString languageName(const QString &lang) {
    if (lang == "es") return QStringLiteral("Español");
    return QStringLiteral("English");
}

QString resolve(const QString &lang) {
    if (languages().contains(lang)) return lang;
    const QString sys = QLocale::system().name().left(2);
    return languages().contains(sys) ? sys : QStringLiteral("en");
}

QString install(const QString &lang) {
    const QString l = resolve(lang);
    if (l == "en") return l;
    auto *own = new TsTranslator(QCoreApplication::instance());
    if (own->loadTs(":/translations/pyfda_cpp_" + l + ".ts")) QCoreApplication::installTranslator(own);
    // Qt's own strings (buttons of the standard dialogs), if the Qt translations are installed
    // (next to the executable: qt_<lang>.qm from windeployqt --translations)
    auto *qt = new QTranslator(QCoreApplication::instance());
    if (qt->load("qt_" + l, QCoreApplication::applicationDirPath() + "/translations") ||
        qt->load("qtbase_" + l, QLibraryInfo::path(QLibraryInfo::TranslationsPath)))
        QCoreApplication::installTranslator(qt);
    return l;
}

}  // namespace i18n
