// Filter specifications (pyfda's "Specs" tab: select filter, order, frequency
// and amplitude specs) and the "DESIGN FILTER" button.
#pragma once

#include "filter_design.hpp"

#include <QMap>
#include <QWidget>

class QComboBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QRadioButton;
class QSpinBox;
class QGridLayout;

class SpecPanel : public QWidget {
    Q_OBJECT
public:
    explicit SpecPanel(QWidget *parent = nullptr);

    /// Current specifications, throws pyfda::DesignError for invalid entries
    pyfda::FilterSpec spec() const;
    /// Take over order, corner frequencies etc. calculated by the design
    void updateFromDesign(const pyfda::FilterSpec &s);
    void setStatus(const QString &text, bool error);

    /// Factor from the frequency unit to Hz, 0 for normalized frequencies
    double unitToHz() const;
    /// Frequency axis label, e.g. "f / kHz" or "F = f / f_S"
    QString freqLabel() const;
    /// Time axis label, e.g. "t / ms" or "n"
    QString timeLabel() const;
    /// Factor that converts n / f_S into the time unit of timeLabel()
    double timeScale() const;

signals:
    void designRequested();
    void unitsChanged();

private slots:
    void updateVisibility();
    void onRespTypeChanged();
    void onUnitChanged();
    void onFsEdited();

private:
    struct Freqs {
        double f_pb, f_pb2, f_sb, f_sb2, f_c, f_c2;
    };
    void storeFreqs();
    void loadFreqs();
    void addRow(QGridLayout *g, int row, const QString &label, QLineEdit *&edit, QLabel *&lbl,
                const QString &tip);
    static double parse(const QLineEdit *e, const QString &name);
    static QString fmt(double v);

    QComboBox *m_rt, *m_ft, *m_method, *m_unit, *m_window, *m_alg;
    QRadioButton *m_min, *m_man;
    QSpinBox *m_N;
    QLineEdit *m_fs, *m_fpb, *m_fpb2, *m_fsb, *m_fsb2, *m_fc, *m_fc2, *m_apb, *m_asb, *m_wpb, *m_wsb,
        *m_winpar;
    QLabel *m_lfpb, *m_lfpb2, *m_lfsb, *m_lfsb2, *m_lfc, *m_lfc2, *m_lapb, *m_lasb, *m_lwpb, *m_lwsb,
        *m_lwinpar, *m_lwindow, *m_lalg, *m_status, *m_lfs;
    QPushButton *m_design;
    QMap<int, Freqs> m_freqs;  // frequency specs per response type
    int m_cur_rt = 0;
    double m_fs_prev = 1.0;
};
