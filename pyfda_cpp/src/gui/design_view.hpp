// Common base class of all tabs that display the current filter design
#pragma once

#include "filter_design.hpp"

#include <QString>
#include <QWidget>

/// Display settings shared by all views
struct ViewContext {
    double f_s = 1.0;           // sampling frequency in the display unit
    double unit_to_hz = 0.0;    // factor display unit -> Hz, 0 = normalized
    QString f_label = "F = f / f_S";
    QString t_label = "n";
};

class DesignView : public QWidget {
    Q_OBJECT
public:
    using QWidget::QWidget;
    /// Called whenever the design or the units change, the view redraws when visible
    void setDesign(const pyfda::FilterDesign *d, const ViewContext &ctx) {
        m_design = d;
        m_ctx = ctx;
        designChanged();
        m_needs_redraw = true;
        if (isVisible()) redrawNow();
    }

protected:
    void showEvent(QShowEvent *e) override {
        QWidget::showEvent(e);
        if (m_needs_redraw) redrawNow();
    }
    void redrawNow() {
        m_needs_redraw = false;
        redraw();
    }
    virtual void redraw() = 0;
    /// Hook for views that have to recalculate something when the design changes
    virtual void designChanged() {}

    const pyfda::FilterDesign *m_design = nullptr;
    ViewContext m_ctx;
    bool m_needs_redraw = true;
};
