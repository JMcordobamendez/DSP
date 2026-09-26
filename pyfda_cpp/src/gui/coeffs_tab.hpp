// Table with the filter coefficients (b, a) or second-order sections,
// copy to the clipboard and CSV export
#pragma once

#include "design_view.hpp"

class QComboBox;
class QTableWidget;
class QLabel;

class CoeffsView : public DesignView {
    Q_OBJECT
public:
    explicit CoeffsView(QWidget *parent = nullptr);

protected:
    void redraw() override;

private:
    QString asText(QChar sep) const;
    void exportCsv();

    QComboBox *m_format;
    QTableWidget *m_table;
    QLabel *m_info;
};
