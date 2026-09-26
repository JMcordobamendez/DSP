#include "info_tab.hpp"

#include "filter_info.hpp"

#include <QTextBrowser>
#include <QVBoxLayout>

#include <cmath>

using namespace pyfda;

InfoView::InfoView(QWidget *parent) : DesignView(parent) {
    auto *lay = new QVBoxLayout(this);
    lay->setContentsMargins(4, 4, 4, 4);
    m_text = new QTextBrowser(this);
    lay->addWidget(m_text);
}

QString InfoView::html() const { return m_text->toHtml(); }

void InfoView::redraw() {
    if (!m_design) {
        m_text->setHtml(tr("<p>No filter designed.</p>"));
        return;
    }
    const FilterInfo fi = filter_info(*m_design);
    const FilterSpec &s = m_design->spec;
    auto num = [](double v, int prec = 4) { return QString::number(v, 'g', prec); };
    auto yes = [this](bool b, bool good = true) {
        return QString("<span style='color:%1'>%2</span>").arg(b == good ? "#008000" : "#c00000", b ? tr("yes") : tr("no"));
    };
    const QString fu = m_ctx.unit_to_hz == 0 ? QString() : " " + m_ctx.f_label.section("/ ", 1);
    QString h = "<h3>" + QString::fromStdString(m_design->info).toHtmlEscaped() + "</h3>";
    h += "<table cellspacing='0' cellpadding='3' border='0'>";
    auto row = [&](const QString &k, const QString &v) { h += "<tr><td><b>" + k + "</b></td><td>" + v + "</td></tr>"; };
    row(tr("Type"), m_design->fir ? "FIR" : "IIR");
    row(tr("Order N"), QString::number(fi.order));
    row(tr("Coefficients"), QString("b: %1, a: %2").arg(fi.n_b).arg(fi.n_a) +
                                (fi.n_sos ? tr(", %1 second-order sections").arg(fi.n_sos) : QString()));
    row(tr("f_S"), num(s.f_s, 8) + fu);
    if (!m_design->fir) row(tr("Stable"), yes(fi.stable) + tr(" (max. |p| = %1)").arg(num(fi.max_pole_radius, 8)));
    row(tr("Minimum phase"), yes(fi.min_phase));
    if (m_design->fir) row(tr("Linear phase"), yes(fi.linear_phase));
    row(tr("|H(0)|"), num(fi.gain_dc, 6) + QString(" (%1 dB)").arg(num(20 * std::log10(std::max(fi.gain_dc, 1e-15)), 5)));
    row(tr("|H(f_S / 2)|"), num(fi.gain_ny, 6) + QString(" (%1 dB)").arg(num(20 * std::log10(std::max(fi.gain_ny, 1e-15)), 5)));
    row(tr("max |H|"), num(fi.h_max_db, 5) + " dB");
    h += "</table>";
    if (!fi.bands.empty()) {
        h += "<h4>" + tr("Specifications vs. achieved") + "</h4>";
        h += "<table cellspacing='0' cellpadding='4' border='1' style='border-collapse:collapse'>";
        h += "<tr><th>" + tr("Band") + "</th><th>" + tr("Frequencies") + "</th><th>" + tr("Spec") + "</th><th>" +
             tr("Achieved") + "</th><th></th></tr>";
        for (const BandCheck &b : fi.bands) {
            h += "<tr><td>" + QString::fromStdString(b.name) + "</td><td>" + num(b.f0, 6) + " ... " + num(b.f1, 6) + fu +
                 "</td><td>" + (b.pass ? tr("ripple &le; %1 dB") : tr("attenuation &ge; %1 dB")).arg(num(b.spec_db, 5)) +
                 "</td><td>" + num(b.achieved_db, 5) + " dB</td><td style='color:" + (b.ok ? "#008000'>OK" : "#c00000'>" + tr("not met")) +
                 "</td></tr>";
        }
        h += "</table>";
        bool all = true;
        for (const BandCheck &b : fi.bands) all = all && b.ok;
        if (!all)
            h += "<p>" + tr("The minimum order estimate isn't always sufficient (e.g. for equiripple filters), "
                            "increase the order or relax the specs.") + "</p>";
    } else {
        h += "<p>" + tr("No band edge specifications (manual order / manual filter).") + "</p>";
    }
    m_text->setHtml(h);
}
