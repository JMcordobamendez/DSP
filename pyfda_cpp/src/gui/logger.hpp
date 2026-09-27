// Application wide logger, messages are shown in the log pane of the main window
#pragma once

#include <QObject>
#include <QString>

class Logger : public QObject {
    Q_OBJECT
public:
    enum Level { Info, Warning, Error };
    static Logger &instance() {
        static Logger l;
        return l;
    }
    static void info(const QString &s) { emit instance().message(Info, s); }
    static void warning(const QString &s) { emit instance().message(Warning, s); }
    static void error(const QString &s) { emit instance().message(Error, s); }
signals:
    void message(int level, const QString &text);
};
