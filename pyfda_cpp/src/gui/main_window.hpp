#pragma once

#include "design_view.hpp"
#include "filter_design.hpp"

#include <QMainWindow>
#include <memory>

class SpecPanel;
class QTabWidget;
class QPlainTextEdit;
class DataFiltView;
class CoeffsView;
class TranView;
class FixpointView;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = nullptr);
    /// Design the filter with the current specs, returns false on errors
    bool design();
    SpecPanel *specPanel() const { return m_specs; }
    DataFiltView *dataFilt() const { return m_data_filt; }
    /// Load a filter file (JSON) and design it, returns false on errors
    bool openFilter(const QString &file_name);
    /// Save the current design as a filter file (JSON)
    bool saveFilter(const QString &file_name);
    QTabWidget *tabs() const { return m_tabs; }
    TranView *tran() const { return m_tran; }
    FixpointView *fixpoint() const { return m_fix; }

private slots:
    void onLog(int level, const QString &text);
    void updateViews();

private:
    ViewContext context() const;

    SpecPanel *m_specs;
    QTabWidget *m_tabs;
    QPlainTextEdit *m_log;
    DataFiltView *m_data_filt;
    CoeffsView *m_coeffs;
    TranView *m_tran;
    FixpointView *m_fix;
    QString m_filter_dir;
    QList<DesignView *> m_views;
    std::unique_ptr<pyfda::FilterDesign> m_design;
};
