// Fixpoint settings (pyfda's "Fixpoint" input tab, input_fixpoint_specs.py with the
// fir_df / iir_df1 widgets): word formats of input, coefficients, accumulator and
// output, quantized coefficients, quantized magnitude response and export of the
// fixpoint filter as Xilinx COE file or VHDL.
#pragma once

#include "design_view.hpp"
#include "fixpoint.hpp"

class PlotWidget;
class QCheckBox;
class QComboBox;
class QLabel;
class QPushButton;
class QSpinBox;
class QTableWidget;

class FixpointView : public DesignView {
    Q_OBJECT
public:
    explicit FixpointView(QWidget *parent = nullptr);
    const pyfda::FxSpec &spec() const { return m_spec; }
    bool simulate() const;
    void setSpec(const pyfda::FxSpec &s, bool simulate);
    /// Export as COE (FIR) or VHDL, the format is taken from the suffix (.coe / .vhd)
    /// Export by suffix: .coe, .vhd / .vhdl or .v; with_tb also writes a
    /// self-checking testbench <name>_tb.vhd / <name>_tb.v
    bool exportFile(const QString &file_name, bool with_tb = false);

signals:
    /// Fixpoint settings or the "simulate" switch changed
    void fxChanged();

protected:
    void redraw() override;
    void designChanged() override;

private:
    struct Row {
        QSpinBox *wi, *wf;
        QComboBox *quant, *ovfl;
        QLabel *w, *label;
    };
    Row addRow(int row, const QString &label, const QString &tip);
    void readUi();
    void updateUi();
    void onChanged();
    pyfda::QFormat format(const Row &r) const;
    void setFormat(Row &r, const pyfda::QFormat &q);

    pyfda::FxSpec m_spec;
    bool m_updating = false;
    QCheckBox *m_sim, *m_auto_c, *m_auto_acc;
    Row m_qi, m_qcb, m_qca, m_qacc, m_qo;
    QLabel *m_structure, *m_ovfl_info;
    QComboBox *m_base;
    QTableWidget *m_table;
    PlotWidget *m_plot;
    QPushButton *m_coe, *m_vhdl;
};
