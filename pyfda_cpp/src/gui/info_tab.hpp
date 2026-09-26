// Filter properties and specifications vs. achieved values (pyfda's "Info" tab)
#pragma once

#include "design_view.hpp"

class QTextBrowser;

class InfoView : public DesignView {
    Q_OBJECT
public:
    explicit InfoView(QWidget *parent = nullptr);
    QString html() const;

protected:
    void redraw() override;

private:
    QTextBrowser *m_text;
};
