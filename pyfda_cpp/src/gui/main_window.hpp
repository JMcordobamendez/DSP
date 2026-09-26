#pragma once

#include "design_view.hpp"
#include "filter_design.hpp"

#include <QMainWindow>
#include <memory>

class SpecPanel;
class QTabWidget;
class QPlainTextEdit;
class DataFiltView;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = nullptr);
    /// Design the filter with the current specs, returns false on errors
    bool design();
    SpecPanel *specPanel() const { return m_specs; }
    DataFiltView *dataFilt() const { return m_data_filt; }
    QTabWidget *tabs() const { return m_tabs; }

private slots:
    void onLog(int level, const QString &text);
    void updateViews();

private:
    ViewContext context() const;

    SpecPanel *m_specs;
    QTabWidget *m_tabs;
    QPlainTextEdit *m_log;
    DataFiltView *m_data_filt;
    QList<DesignView *> m_views;
    std::unique_ptr<pyfda::FilterDesign> m_design;
};
