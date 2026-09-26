// Window viewer (pyfda's FFT window viewer, plot_widgets/plot_fft_win.py): time
// and frequency domain of a window with its figures of merit
#pragma once

#include <QDialog>

class PlotWidget;
class QComboBox;
class QLabel;
class QLineEdit;
class QSpinBox;
class QCheckBox;

class WindowViewer : public QDialog {
    Q_OBJECT
public:
    explicit WindowViewer(QWidget *parent = nullptr);
    void setWindow(int index, int N, double par);

private:
    void update();

    QComboBox *m_win;
    QLineEdit *m_par, *m_par2;
    QLabel *m_lpar2;
    QLabel *m_lpar, *m_props;
    QSpinBox *m_N;
    QCheckBox *m_log;
    PlotWidget *m_time, *m_freq;
};
