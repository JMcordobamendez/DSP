// Table with the filter coefficients (b, a) or second-order sections,
// copy to the clipboard and export (CSV, MATLAB, C, Python). In edit mode
// b / a or poles / zeros can be entered by hand (pyfda's input_coeffs and
// input_pz widgets), "Apply" turns them into a "Manual" design.
#pragma once

#include "design_view.hpp"

class QCheckBox;
class QComboBox;
class QLineEdit;
class QPushButton;
class QTableWidget;
class QLabel;

class CoeffsView : public DesignView {
    Q_OBJECT
public:
    explicit CoeffsView(QWidget *parent = nullptr);
    /// Ask for a file name and format and export the coefficients
    void exportDialog();
    /// Edit mode on / off; format 0 = b, a; 2 = zeros / poles
    void setEditing(bool on, int format = -1);
    /// Set a table cell (for tests and the command line)
    void setCell(int row, int col, const QString &text);
    /// Parse the edited table into manual coefficients / poles and zeros,
    /// throws pyfda::DesignError for invalid entries
    pyfda::FilterSpec editedSpec() const;

signals:
    /// "Apply" in edit mode: design a manual filter with these coefficients / poles and zeros
    void manualDesignRequested(const pyfda::FilterSpec &manual);

protected:
    void redraw() override;

private:
    QString asText(QChar sep) const;
    void updateEditWidgets();
    void apply();

    QComboBox *m_format;
    QTableWidget *m_table;
    QLabel *m_info, *m_lk;
    QCheckBox *m_edit;
    QLineEdit *m_k;
    QPushButton *m_add, *m_del, *m_apply, *m_reset;
};
