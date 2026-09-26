#include "coeffs_tab.hpp"

#include "logger.hpp"

#include <QApplication>
#include <QClipboard>
#include <QComboBox>
#include <QFile>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
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
    auto *exp = new QPushButton(tr("Export CSV ..."), this);
    ctl->addWidget(new QLabel(tr("Format:"), this));
    ctl->addWidget(m_format);
    ctl->addWidget(copy);
    ctl->addWidget(exp);
    ctl->addStretch(1);
    lay->addLayout(ctl);
    m_info = new QLabel(this);
    lay->addWidget(m_info);
    m_table = new QTableWidget(this);
    m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    lay->addWidget(m_table, 1);
    connect(m_format, &QComboBox::currentIndexChanged, this, [this] { redrawNow(); });
    connect(copy, &QPushButton::clicked, this, [this] { QApplication::clipboard()->setText(asText('\t')); });
    connect(exp, &QPushButton::clicked, this, &CoeffsView::exportCsv);
}

void CoeffsView::redraw() {
    m_table->clear();
    m_table->setRowCount(0);
    m_table->setColumnCount(0);
    if (!m_design) return;
    const int f = m_format->currentIndex();
    auto num = [](double v) { return new QTableWidgetItem(QString::number(v, 'g', 17)); };
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
        m_info->setText(tr("%1 coefficients, order N = %2").arg(n).arg(m_design->spec.N));
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

void CoeffsView::exportCsv() {
    const QString fn = QFileDialog::getSaveFileName(this, tr("Export coefficients"), QString(), tr("CSV (*.csv)"));
    if (fn.isEmpty()) return;
    QFile f(fn);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        Logger::error(tr("Couldn't write '%1'.").arg(fn));
        return;
    }
    QTextStream(&f) << asText(',');
    Logger::info(tr("Exported coefficients to '%1'.").arg(fn));
}
