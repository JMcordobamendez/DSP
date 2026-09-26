// "Data Filt" tab: load a column of data from a file (csv, txt, wav, npy), filter it
// with the current design and plot original and filtered data (and optionally their
// spectra) overlaid. Port of pyfda/plot_widgets/plot_data_filt.py
#pragma once

#include "data_io.hpp"
#include "design_view.hpp"

#include <optional>

class PlotWidget;
class QCheckBox;
class QComboBox;
class QLabel;
class QPushButton;
class QFrame;

class DataFiltView : public DesignView {
    Q_OBJECT
public:
    explicit DataFiltView(QWidget *parent = nullptr);
    /// Load a data file, returns false (and logs the error) on failure
    bool loadFile(const QString &file_name);
    /// Filter the selected column with the current design
    void filterData();
    bool saveCsv(const QString &file_name);

protected:
    void redraw() override;
    void designChanged() override;

private slots:
    void onLoad();
    void selectData();
    void exportData();

private:
    void refilter();
    bool computeFilter();
    void checkFs();
    pyfda::Vec timeAxis(QString &label) const;

    QPushButton *m_but_load, *m_but_filter, *m_but_export;
    QLabel *m_lbl_file, *m_lbl_fs_warn;
    QComboBox *m_cmb_col, *m_cmb_time;
    QCheckBox *m_chk_zero_phase, *m_chk_spectrum;
    PlotWidget *m_plot_t, *m_plot_f;

    std::optional<pyfda::DataTable> m_data;
    QString m_file_name;
    pyfda::Vec m_x;                  // selected data column
    std::optional<pyfda::Vec> m_y;   // filtered data
};
