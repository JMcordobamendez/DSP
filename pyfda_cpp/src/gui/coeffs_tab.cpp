#include "coeffs_tab.hpp"

#include "filter_io.hpp"
#include "logger.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QLineEdit>
#include <QLocale>
#include <QClipboard>
#include <QComboBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QStandardItemModel>

#include <algorithm>
#include <cmath>
#include <QTableWidget>
#include <QTextStream>
#include <QVBoxLayout>

using namespace pyfda;

CoeffsView::CoeffsView(QWidget *parent) : DesignView(parent) {
    auto *lay = new QVBoxLayout(this);
    auto *ctl = new QHBoxLayout();
    m_format = new QComboBox(this);
    m_format->addItems({tr("b, a (transfer function)"), tr("Second-order sections"), tr("Zeros / poles")});
    auto *copy = new QPushButton(tr("Copy"), this);
    copy->setToolTip(tr("Copy the table to the clipboard (tab separated)"));
    auto *exp = new QPushButton(tr("Export ..."), this);
    exp->setToolTip(tr("Export the coefficients for MATLAB / Octave, C, Python or as CSV"));
    ctl->addWidget(new QLabel(tr("Format:"), this));
    ctl->addWidget(m_format);
    ctl->addWidget(copy);
    ctl->addWidget(exp);
    ctl->addStretch(1);
    m_edit = new QCheckBox(tr("Edit"), this);
    m_edit->setToolTip(tr("<span>Enter the coefficients b, a or the poles / zeros by hand; "
                          "<i>Apply</i> designs a <i>Manual</i> filter from them.</span>"));
    ctl->addWidget(m_edit);
    lay->addLayout(ctl);
    auto *ed = new QHBoxLayout();
    m_lk = new QLabel("k =", this);
    m_k = new QLineEdit(this);
    m_k->setToolTip(tr("Gain factor k"));
    m_k->setMaximumWidth(180);
    m_add = new QPushButton(tr("Add row"), this);
    m_add->setToolTip(tr("Insert a row above the current row (at the end without selection)"));
    m_del = new QPushButton(tr("Delete rows"), this);
    m_del->setToolTip(tr("Delete the selected rows"));
    m_reset = new QPushButton(tr("Reset"), this);
    m_reset->setToolTip(tr("Discard the changes and show the current filter"));
    m_apply = new QPushButton(tr("Apply"), this);
    m_apply->setToolTip(tr("Design a manual filter from the table (Ctrl+Enter)"));
    m_apply->setShortcut(QKeySequence("Ctrl+Return"));
    QFont bf = m_apply->font();
    bf.setBold(true);
    m_apply->setFont(bf);
    for (QWidget *w : {static_cast<QWidget *>(m_lk), static_cast<QWidget *>(m_k)}) ed->addWidget(w);
    ed->addStretch(1);
    for (QPushButton *b : {m_add, m_del, m_reset, m_apply}) ed->addWidget(b);
    lay->addLayout(ed);
    m_info = new QLabel(this);
    lay->addWidget(m_info);
    m_table = new QTableWidget(this);
    m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    lay->addWidget(m_table, 1);
    connect(m_format, &QComboBox::currentIndexChanged, this, [this] { redrawNow(); });
    connect(copy, &QPushButton::clicked, this, [this] { QApplication::clipboard()->setText(asText('\t')); });
    connect(exp, &QPushButton::clicked, this, &CoeffsView::exportDialog);
    connect(m_edit, &QCheckBox::toggled, this, [this](bool on) {
        if (on && m_format->currentIndex() == 1) m_format->setCurrentIndex(0);  // sections aren't editable
        updateEditWidgets();
        redrawNow();
    });
    connect(m_add, &QPushButton::clicked, this, [this] {
        const int r = m_table->currentRow() >= 0 && !m_table->selectedItems().isEmpty() ? m_table->currentRow()
                                                                                      : m_table->rowCount();
        m_table->insertRow(r);
        for (int c = 0; c < m_table->columnCount(); ++c) m_table->setItem(r, c, new QTableWidgetItem("0"));
    });
    connect(m_del, &QPushButton::clicked, this, [this] {
        QList<int> rows;
        for (const auto &rng : m_table->selectedRanges())
            for (int r = rng.topRow(); r <= rng.bottomRow(); ++r)
                if (!rows.contains(r)) rows << r;
        std::sort(rows.begin(), rows.end(), std::greater<int>());
        for (int r : rows) m_table->removeRow(r);
    });
    connect(m_reset, &QPushButton::clicked, this, [this] { redrawNow(); });
    connect(m_apply, &QPushButton::clicked, this, &CoeffsView::apply);
    updateEditWidgets();
}

void CoeffsView::updateEditWidgets() {
    const bool on = m_edit->isChecked();
    const bool zpk = m_format->currentIndex() == 2;
    for (QWidget *w : {static_cast<QWidget *>(m_add), static_cast<QWidget *>(m_del),
                       static_cast<QWidget *>(m_apply), static_cast<QWidget *>(m_reset)})
        w->setVisible(on);
    m_lk->setVisible(on && zpk);
    m_k->setVisible(on && zpk);
    m_table->setEditTriggers(on ? QAbstractItemView::AllEditTriggers : QAbstractItemView::NoEditTriggers);
    // second-order sections can't be edited
    auto *model = qobject_cast<QStandardItemModel *>(m_format->model());
    if (model) model->item(1)->setEnabled(!on);
}

void CoeffsView::setEditing(bool on, int format) {
    if (format >= 0) m_format->setCurrentIndex(format);
    m_edit->setChecked(on);
    updateEditWidgets();
    redrawNow();
}

void CoeffsView::setCell(int row, int col, const QString &text) {
    if (row >= m_table->rowCount()) m_table->setRowCount(row + 1);
    m_table->setItem(row, col, new QTableWidgetItem(text));
}

FilterSpec CoeffsView::editedSpec() const {
    auto cell = [this](int r, int c, bool &empty) {
        const QTableWidgetItem *it = m_table->item(r, c);
        QString t = it ? it->text().trimmed() : QString();
        empty = t.isEmpty();
        if (empty) return 0.0;
        if (t.contains(',') && !t.contains('.')) t.replace(',', '.');
        bool ok = false;
        const double v = t.toDouble(&ok);
        if (!ok || !std::isfinite(v))
            throw DesignError(tr("Invalid value '%1' in row %2, column '%3'.")
                                  .arg(it->text())
                                  .arg(r + 1)
                                  .arg(m_table->horizontalHeaderItem(c)->text())
                                  .toStdString());
        return v;
    };
    // a column without trailing empty cells, empty cells in between count as 0
    auto column = [&](int c) {
        Vec v;
        int last = -1;
        std::vector<double> vals(m_table->rowCount());
        for (int r = 0; r < m_table->rowCount(); ++r) {
            bool empty;
            vals[r] = cell(r, c, empty);
            if (!empty) last = r;
        }
        v.assign(vals.begin(), vals.begin() + (last + 1));
        return v;
    };
    FilterSpec s;
    if (m_format->currentIndex() == 2) {
        auto roots = [&](int c) {
            CVec out;
            for (int r = 0; r < m_table->rowCount(); ++r) {
                bool e_re, e_im;
                const double re = cell(r, c, e_re), im = cell(r, c + 1, e_im);
                if (!(e_re && e_im)) out.emplace_back(re, im);
            }
            return out;
        };
        s.manual_zpk.z = roots(0);
        s.manual_zpk.p = roots(2);
        QString t = m_k->text().trimmed();
        if (t.contains(',') && !t.contains('.')) t.replace(',', '.');
        bool ok = false;
        s.manual_zpk.k = t.toDouble(&ok);
        if (!ok || !std::isfinite(s.manual_zpk.k))
            throw DesignError(tr("Invalid gain k = '%1'.").arg(m_k->text()).toStdString());
        s.manual_from_zpk = true;
    } else {
        s.manual_ba.b = column(0);
        s.manual_ba.a = column(1);
        if (s.manual_ba.b.empty()) throw DesignError("Enter at least one coefficient b.");
        if (s.manual_ba.a.empty()) s.manual_ba.a = {1.0};
        s.manual_from_zpk = false;
    }
    return s;
}

void CoeffsView::apply() {
    try {
        emit manualDesignRequested(editedSpec());
    } catch (const DesignError &e) {
        Logger::error(e.what());
    }
}

void CoeffsView::redraw() {
    updateEditWidgets();
    m_table->clear();
    m_table->setRowCount(0);
    m_table->setColumnCount(0);
    if (!m_design) return;
    const int f = m_format->currentIndex();
    auto num = [](double v) { return new QTableWidgetItem(QString::number(v, 'g', QLocale::FloatingPointShortest)); };
    if (f == 0) {
        const Vec &b = m_design->ba.b, &a = m_design->ba.a;
        const int n = int(std::max(b.size(), a.size()));
        m_table->setColumnCount(2);
        m_table->setHorizontalHeaderLabels({"b", "a"});
        m_table->setRowCount(n);
        for (int i = 0; i < n; ++i) {
            if (i < int(b.size())) m_table->setItem(i, 0, num(b[i]));
            if (i < int(a.size())) m_table->setItem(i, 1, num(a[i]));
        }
        m_info->setText(tr("%1 coefficients, order N = %2").arg(n).arg(n - 1));
    } else if (f == 1) {
        m_table->setColumnCount(6);
        m_table->setHorizontalHeaderLabels({"b0", "b1", "b2", "a0", "a1", "a2"});
        m_table->setRowCount(int(m_design->sos.size()));
        for (int r = 0; r < int(m_design->sos.size()); ++r)
            for (int c = 0; c < 6; ++c) m_table->setItem(r, c, num(m_design->sos[r][c]));
        m_info->setText(m_design->sos.empty() ? tr("FIR filters are implemented directly with b.")
                                              : tr("%1 sections").arg(m_design->sos.size()));
    } else {
        const CVec &z = m_design->zpk.z, &p = m_design->zpk.p;
        const int n = int(std::max(z.size(), p.size()));
        m_table->setColumnCount(4);
        m_table->setHorizontalHeaderLabels({"Re(z)", "Im(z)", "Re(p)", "Im(p)"});
        m_table->setRowCount(n);
        for (int i = 0; i < n; ++i) {
            if (i < int(z.size())) {
                m_table->setItem(i, 0, num(z[i].real()));
                m_table->setItem(i, 1, num(z[i].imag()));
            }
            if (i < int(p.size())) {
                m_table->setItem(i, 2, num(p[i].real()));
                m_table->setItem(i, 3, num(p[i].imag()));
            }
        }
        m_info->setText(tr("Gain k = %1").arg(m_design->zpk.k, 0, 'g', 17));
        m_k->setText(QString::number(m_design->zpk.k, 'g', QLocale::FloatingPointShortest));
    }
}

QString CoeffsView::asText(QChar sep) const {
    QString s;
    QStringList head;
    for (int c = 0; c < m_table->columnCount(); ++c) head << m_table->horizontalHeaderItem(c)->text();
    s += head.join(sep) + "\n";
    for (int r = 0; r < m_table->rowCount(); ++r) {
        QStringList row;
        for (int c = 0; c < m_table->columnCount(); ++c) {
            const QTableWidgetItem *it = m_table->item(r, c);
            row << (it ? it->text() : QString());
        }
        s += row.join(sep) + "\n";
    }
    return s;
}

void CoeffsView::exportDialog() {
    if (!m_design) {
        Logger::warning(tr("No filter designed yet."));
        return;
    }
    const CoeffFormat formats[] = {CoeffFormat::Matlab, CoeffFormat::CHeader, CoeffFormat::Python, CoeffFormat::Csv};
    QStringList filters;
    for (CoeffFormat f : formats) filters << coeff_format_filter(f);
    const QString table = tr("CSV table of the current view (*.csv)");
    filters << table;
    QString selected = filters[0];
    QString fn = QFileDialog::getSaveFileName(this, tr("Export coefficients"), QString(), filters.join(";;"), &selected);
    if (fn.isEmpty()) return;
    if (selected == table) {
        if (QFileInfo(fn).suffix().isEmpty()) fn += ".csv";
        QFile f(fn);
        if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
            Logger::error(tr("Couldn't write '%1'.").arg(fn));
            return;
        }
        redrawNow();
        QTextStream(&f) << asText(',');
        Logger::info(tr("Exported the coefficient table to '%1'.").arg(fn));
        return;
    }
    CoeffFormat fmt = formats[std::max<qsizetype>(0, filters.indexOf(selected))];
    if (QFileInfo(fn).suffix().isEmpty()) fn += QString(".") + coeff_format_suffix(fmt);
    // write with QFile for Unicode file names on Windows
    QFile f(fn);
    const std::string text = export_coeffs(*m_design, fmt, QFileInfo(fn).completeBaseName().toStdString());
    if (!f.open(QIODevice::WriteOnly) || f.write(text.data(), qint64(text.size())) != qint64(text.size())) {
        Logger::error(tr("Couldn't write '%1'.").arg(fn));
        return;
    }
    Logger::info(tr("Exported coefficients to '%1'.").arg(fn));
}
