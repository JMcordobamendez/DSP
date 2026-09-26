// Tabs showing the responses of the current design:
// magnitude, phase, group delay, pole / zero plot and impulse / step response
#pragma once

#include "design_view.hpp"

class PlotWidget;
class QCheckBox;
class QComboBox;
class QHBoxLayout;

/// Base with a control row above a plot
class PlotView : public DesignView {
    Q_OBJECT
public:
    explicit PlotView(QWidget *parent = nullptr);

protected:
    /// frequency axis 0 ... f_S / 2 (display unit) and the corresponding w = 0 ... pi
    void freqAxis(int n, QVector<double> &f, pyfda::Vec &w) const;
    pyfda::CVec response(const pyfda::Vec &w) const;
    QHBoxLayout *m_controls;
    PlotWidget *m_plot;
};

class MagnitudeView : public PlotView {
    Q_OBJECT
public:
    explicit MagnitudeView(QWidget *parent = nullptr);

protected:
    void redraw() override;

private:
    QComboBox *m_unit;
    QCheckBox *m_specs;
};

class PhaseView : public PlotView {
    Q_OBJECT
public:
    explicit PhaseView(QWidget *parent = nullptr);

protected:
    void redraw() override;

private:
    QComboBox *m_unit;
    QCheckBox *m_wrapped;
};

class GroupDelayView : public PlotView {
    Q_OBJECT
public:
    explicit GroupDelayView(QWidget *parent = nullptr);

protected:
    void redraw() override;
};

class PoleZeroView : public PlotView {
    Q_OBJECT
public:
    explicit PoleZeroView(QWidget *parent = nullptr);

protected:
    void redraw() override;
};

class ImpulseView : public PlotView {
    Q_OBJECT
public:
    explicit ImpulseView(QWidget *parent = nullptr);

protected:
    void redraw() override;

private:
    QComboBox *m_kind;
    QCheckBox *m_log;
};
