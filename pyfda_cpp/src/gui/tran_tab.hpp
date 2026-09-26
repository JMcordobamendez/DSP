// Transient analysis (pyfda's "y[n]" tab, plot_widgets/plot_tran.py): stimulus
// x[n] with noise and DC, response y[n] of the current design, both in the time
// domain and as windowed spectra.
#pragma once

#include "design_view.hpp"
#include "fixpoint.hpp"
#include "stimulus.hpp"

#include <QMap>

class PlotWidget;
class QCheckBox;
class QComboBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QTabWidget;

class TranView : public DesignView {
    Q_OBJECT
public:
    explicit TranView(QWidget *parent = nullptr);
    const pyfda::StimParams &params() const { return m_p; }
    void setParams(const pyfda::StimParams &p);
    /// Calculate the response with the fixpoint filter when `on`
    void setFixpoint(const pyfda::FxSpec &spec, bool on);
    /// Stimulus and response as calculated for the current display
    const pyfda::Vec &stimulus() const { return m_x; }
    const pyfda::Vec &response() const { return m_y; }
    bool saveCsv(const QString &file_name);
    QTabWidget *innerTabs() const { return m_tabs; }
    /// Load a data file (CSV, wav, npy ...) as stimulus and select the "File" stimulus
    bool loadStimFile(const QString &file_name);
    /// Formula of the "Formula" stimulus (numexpr syntax)
    void setFormula(const QString &formula);

protected:
    void redraw() override;
    void designChanged() override;

private:
    QLineEdit *addParam(const QString &key, const QString &label, const QString &tip);
    void readEdits();
    void updateEdits();
    void updateVisibility();
    void calc();
    void drawTime();
    void drawFreq();
    void drawSpgr();
    void applyFileNorm();
    /// Periodic window of the spectral analysis with the parameters of the Frequency tab
    pyfda::Vec analysisWindow(int N) const;
    void changed();

    pyfda::StimParams m_p;
    bool m_dirty = true;
    pyfda::Vec m_x, m_y, m_y_float;
    pyfda::Vec m_xi, m_yi, m_yi_float;  // imaginary parts of complex signals
    bool m_cmplx = false;
    pyfda::FxSpec m_fx;
    bool m_fx_on = false;
    QString m_fx_info;
    int m_n_start = 0;
    QString m_error;

    QComboBox *m_stim, *m_chirp, *m_noise, *m_win;
    QCheckBox *m_bl, *m_step_err;
    QMap<QString, QLineEdit *> m_edits;       // parameter edits by key
    QMap<QString, QWidget *> m_param_widgets;  // label + edit containers by key
    QLineEdit *m_noi, *m_win_par, *m_win_par2;
    QSpinBox *m_N, *m_N_start, *m_mls_b;
    QLabel *m_lbl_mls, *m_lbl_win_par, *m_lbl_win_par2, *m_info;
    QTabWidget *m_tabs;
    PlotWidget *m_plot_t, *m_plot_ti, *m_plot_f, *m_plot_s;
    QLabel *m_lbl_cmplx;
    QLineEdit *m_formula;
    QWidget *m_wdg_formula, *m_wdg_file;
    QPushButton *m_load;
    QLabel *m_file_lbl;
    QCheckBox *m_file_norm;
    pyfda::Vec m_file_raw, m_file_raw_im;
    // spectrogram
    QComboBox *m_s_sig, *m_s_mode;
    QCheckBox *m_s_db;
    QSpinBox *m_s_nfft, *m_s_ovlp;
    QCheckBox *m_t_stim, *m_t_resp, *m_t_db, *m_f_stim, *m_f_resp, *m_f_hid, *m_f_db, *m_f_norm;
};
