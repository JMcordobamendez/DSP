#include "response_tabs.hpp"

#include "conversions.hpp"
#include "filtering.hpp"
#include "plot_widget.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

using namespace pyfda;

namespace {
constexpr int N_FREQ = 2048;
const QColor SPEC_COLOR(255, 80, 80, 40);
}  // namespace

PlotView::PlotView(QWidget *parent) : DesignView(parent) {
    auto *lay = new QVBoxLayout(this);
    lay->setContentsMargins(4, 4, 4, 4);
    m_controls = new QHBoxLayout();
    lay->addLayout(m_controls);
    m_plot = new PlotWidget(this);
    lay->addWidget(m_plot, 1);
}

void PlotView::freqAxis(int n, QVector<double> &f, Vec &w) const {
    f.resize(n);
    w.resize(n);
    for (int i = 0; i < n; ++i) {
        w[i] = PI * i / (n - 1);
        f[i] = m_ctx.f_s / 2 * i / (n - 1);
    }
}

CVec PlotView::response(const Vec &w) const {
    if (!m_design->sos.empty()) return freqz(m_design->sos, w);
    return freqz(m_design->ba, w);
}

// ---------------------------------------------------------------------------
MagnitudeView::MagnitudeView(QWidget *parent) : PlotView(parent) {
    m_unit = new QComboBox(this);
    m_unit->addItems({"dB", "V", "W"});
    m_unit->setToolTip(tr("Unit of the magnitude: dB, linear (V) or squared (W)"));
    m_specs = new QCheckBox(tr("Specs"), this);
    m_specs->setChecked(true);
    m_specs->setToolTip(tr("Show the specifications (forbidden regions are shaded)"));
    m_controls->addWidget(new QLabel(tr("Unit:"), this));
    m_controls->addWidget(m_unit);
    m_controls->addWidget(m_specs);
    m_controls->addStretch(1);
    connect(m_unit, &QComboBox::currentIndexChanged, this, [this] { m_plot->keepView(false); redrawNow(); });
    connect(m_specs, &QCheckBox::toggled, this, [this] { redrawNow(); });
}

void MagnitudeView::redraw() {
    m_plot->clear();
    if (!m_design) {
        m_plot->setMessage(tr("Design a filter"));
        m_plot->autoscale();
        return;
    }
    QVector<double> f;
    Vec w;
    freqAxis(N_FREQ, f, w);
    const CVec H = response(w);
    const int unit = m_unit->currentIndex();
    QVector<double> y(N_FREQ);
    for (int i = 0; i < N_FREQ; ++i) {
        const double a = std::abs(H[i]);
        y[i] = unit == 0 ? 20 * std::log10(std::max(a, 1e-15)) : unit == 1 ? a : a * a;
    }
    const FilterSpec &s = m_design->spec;
    if (m_specs->isChecked() && s.fo == OrderMode::Min) {
        // pass band limits
        const double d_pb = m_design->fir ? fir_a_pb_lin(s.A_PB) : 1 - std::pow(10.0, -s.A_PB / 20);
        const double hi_pb = m_design->fir ? 1 + d_pb : 1.0, lo_pb = 1 - d_pb;
        const double sb = sb_lin(s.A_SB);
        auto conv = [&](double a) { return unit == 0 ? 20 * std::log10(a) : unit == 1 ? a : a * a; };
        const double BIG = 1e9;
        auto pass = [&](double f0, double f1) {
            m_plot->addRegion({f0, f1, conv(hi_pb), BIG, SPEC_COLOR});
            m_plot->addRegion({f0, f1, -BIG, conv(lo_pb), SPEC_COLOR});
        };
        auto stop = [&](double f0, double f1) { m_plot->addRegion({f0, f1, conv(sb), BIG, SPEC_COLOR}); };
        const double fn = s.f_s / 2;
        switch (s.rt) {
        case RespType::LP: pass(0, s.f_pb); stop(s.f_sb, fn); break;
        case RespType::HP: stop(0, s.f_sb); pass(s.f_pb, fn); break;
        case RespType::BP: stop(0, s.f_sb); pass(s.f_pb, s.f_pb2); stop(s.f_sb2, fn); break;
        case RespType::BS: pass(0, s.f_pb); stop(s.f_sb, s.f_sb2); pass(s.f_pb2, fn); break;
        }
    }
    m_plot->addCurve(f, y, PlotWidget::color(0), QString());
    m_plot->setTitle(QString::fromStdString(m_design->info));
    m_plot->setXLabel(m_ctx.f_label);
    m_plot->setYLabel(unit == 0 ? "|H(f)| / dB" : unit == 1 ? "|H(f)|" : "|H(f)|²");
    m_plot->setXLimits(0, m_ctx.f_s / 2);
    if (unit == 0) {
        const double mx = *std::max_element(y.begin(), y.end());
        const double floor_db = std::max(-std::max(s.A_SB, 40.0) - 40, mx - 200);
        m_plot->setYLimits(floor_db, std::numeric_limits<double>::quiet_NaN());
    }
    m_plot->autoscale();
}

// ---------------------------------------------------------------------------
PhaseView::PhaseView(QWidget *parent) : PlotView(parent) {
    m_unit = new QComboBox(this);
    m_unit->addItems({"rad", "rad / pi", "deg"});
    m_wrapped = new QCheckBox(tr("Wrapped"), this);
    m_wrapped->setToolTip(tr("Show the phase wrapped to -pi ... pi"));
    m_controls->addWidget(new QLabel(tr("Unit:"), this));
    m_controls->addWidget(m_unit);
    m_controls->addWidget(m_wrapped);
    m_controls->addStretch(1);
    connect(m_unit, &QComboBox::currentIndexChanged, this, [this] { redrawNow(); });
    connect(m_wrapped, &QCheckBox::toggled, this, [this] { redrawNow(); });
}

void PhaseView::redraw() {
    m_plot->clear();
    if (!m_design) {
        m_plot->autoscale();
        return;
    }
    QVector<double> f;
    Vec w;
    freqAxis(N_FREQ, f, w);
    const CVec H = response(w);
    Vec phi(N_FREQ);
    for (int i = 0; i < N_FREQ; ++i) phi[i] = std::arg(H[i]);
    if (!m_wrapped->isChecked()) phi = unwrap(phi);
    const int u = m_unit->currentIndex();
    const double scale = u == 0 ? 1.0 : u == 1 ? 1.0 / PI : 180.0 / PI;
    QVector<double> y(N_FREQ);
    for (int i = 0; i < N_FREQ; ++i) y[i] = phi[i] * scale;
    m_plot->addCurve(f, y, PlotWidget::color(0));
    m_plot->setTitle(tr("Phase response"));
    m_plot->setXLabel(m_ctx.f_label);
    m_plot->setYLabel(QString("phi(f) / ") + m_unit->currentText());
    m_plot->setXLimits(0, m_ctx.f_s / 2);
    m_plot->autoscale();
}

// ---------------------------------------------------------------------------
GroupDelayView::GroupDelayView(QWidget *parent) : PlotView(parent) { m_controls->addStretch(1); }

void GroupDelayView::redraw() {
    m_plot->clear();
    if (!m_design) {
        m_plot->autoscale();
        return;
    }
    QVector<double> f;
    Vec w;
    freqAxis(N_FREQ, f, w);
    const Vec gd = m_design->sos.empty() ? group_delay(m_design->ba, w) : group_delay(m_design->sos, w);
    const bool norm = m_ctx.unit_to_hz == 0;
    QVector<double> y(N_FREQ);
    for (int i = 0; i < N_FREQ; ++i) y[i] = norm ? gd[i] : gd[i] / m_ctx.f_s;
    m_plot->addCurve(f, y, PlotWidget::color(0));
    m_plot->setTitle(tr("Group delay"));
    m_plot->setXLabel(m_ctx.f_label);
    m_plot->setYLabel(norm ? QString("tau_g(f) / samples") : "tau_g(f) / " + m_ctx.t_label.mid(4));
    m_plot->setXLimits(0, m_ctx.f_s / 2);
    m_plot->autoscale();
}

// ---------------------------------------------------------------------------
PoleZeroView::PoleZeroView(QWidget *parent) : PlotView(parent) {
    m_controls->addWidget(new QLabel(tr("o: zeros, x: poles"), this));
    m_drag = new QCheckBox(tr("Drag P/Z"), this);
    m_drag->setToolTip(tr("<span>Move poles and zeros with the mouse. The filter becomes a 'Manual' design, "
                          "complex poles / zeros keep their complex conjugate partner, real ones stay on the "
                          "real axis. Ctrl + drag pans the view.</span>"));
    m_controls->addSpacing(12);
    m_controls->addWidget(m_drag);
    m_controls->addStretch(1);
    connect(m_drag, &QCheckBox::toggled, this, [this] { redrawNow(); });
    connect(m_plot, &PlotWidget::pointDragged, this, &PoleZeroView::dragTo);
}

void PoleZeroView::setDragEnabled(bool on) { m_drag->setChecked(on); }

void PoleZeroView::dragTo(int index, QPointF pos, bool finished) {
    if (!m_design) return;
    if (!m_editing) {  // start of the drag: work on a copy of the current poles and zeros
        m_edit = m_design->zpk;
        if (m_design->fir) m_edit.p.clear();  // FIR: poles at the origin are added by the design
        m_editing = true;
        m_plot->keepView(true);
    }
    const bool zero = index < int(m_edit.z.size());
    CVec &v = zero ? m_edit.z : m_edit.p;
    const size_t i = size_t(zero ? index : index - int(m_edit.z.size()));
    if (i >= v.size()) return;
    const cplx old = v[i];
    if (std::fabs(old.imag()) < 1e-10) {
        v[i] = cplx(pos.x(), 0.0);  // real roots stay on the real axis
    } else {
        // complex conjugate partner: the root closest to conj(old)
        size_t partner = i;
        double best = std::numeric_limits<double>::infinity();
        for (size_t k = 0; k < v.size(); ++k)
            if (k != i && std::abs(v[k] - std::conj(old)) < best) {
                best = std::abs(v[k] - std::conj(old));
                partner = k;
            }
        cplx n(pos.x(), pos.y());
        if (std::fabs(n.imag()) < 1e-6) n.imag(old.imag() > 0 ? 1e-6 : -1e-6);  // don't merge on the real axis
        if ((n.imag() > 0) != (old.imag() > 0)) n.imag(-n.imag());  // stay in the half plane
        v[i] = n;
        if (partner != i) v[partner] = std::conj(n);
    }
    emit zpkEdited(m_edit, finished);
    if (finished) {
        m_editing = false;
        m_plot->keepView(false);
    }
}

void PoleZeroView::redraw() {
    m_plot->clear();
    if (!m_design) {
        m_plot->autoscale();
        return;
    }
    m_plot->setEqualAspect(true);
    m_plot->addCircle(0, 0, 1, palette().color(QPalette::Text));
    auto add = [&](const CVec &v, PlotWidget::Marker mk, const QColor &c, const QString &name) {
        PlotWidget::Curve cu;
        for (const cplx &x : v) {
            if (!std::isfinite(x.real()) || !std::isfinite(x.imag())) continue;
            cu.x.append(x.real());
            cu.y.append(x.imag());
        }
        cu.style = PlotWidget::Style::Markers;
        cu.marker = mk;
        cu.color = c;
        cu.name = name;
        m_plot->addCurve(cu);
        // label multiplicities
        std::map<std::pair<long long, long long>, int> cnt;
        for (const cplx &x : v)
            ++cnt[{std::llround(x.real() * 1e6), std::llround(x.imag() * 1e6)}];
        for (const auto &kv : cnt)
            if (kv.second > 1)
                m_plot->addLabel({kv.first.first * 1e-6, kv.first.second * 1e-6, QString::number(kv.second), c});
    };
    CVec p = m_design->zpk.p;
    if (m_design->fir) p.assign(m_design->ba.b.size() - 1, cplx(0.0));  // FIR: poles at the origin
    add(m_design->zpk.z, PlotWidget::Marker::Circle, PlotWidget::color(0), tr("Zeros"));
    add(p, PlotWidget::Marker::Cross, PlotWidget::color(3), tr("Poles"));
    // drag points: zeros, then poles (not the fixed poles at the origin of FIR filters)
    QVector<QPointF> drag;
    if (m_drag->isChecked()) {
        const Zpk &src = m_editing ? m_edit : m_design->zpk;
        for (const cplx &x : src.z) drag << QPointF(x.real(), x.imag());
        if (!m_design->fir || m_editing)
            for (const cplx &x : src.p) drag << QPointF(x.real(), x.imag());
    }
    m_plot->setDragPoints(drag);
    double pmax = 0;
    for (const cplx &x : m_design->zpk.p) pmax = std::max(pmax, std::abs(x));
    m_plot->setTitle(m_design->fir ? tr("Pole / zero plot") :
                     tr("Pole / zero plot, max. |p| = %1%2").arg(pmax, 0, 'g', 6)
                         .arg(pmax >= 1 ? tr(" (unstable!)") : QString()));
    m_plot->setXLabel("Re(z)");
    m_plot->setYLabel("Im(z)");
    m_plot->autoscale();
}

// ---------------------------------------------------------------------------
ImpulseView::ImpulseView(QWidget *parent) : PlotView(parent) {
    m_kind = new QComboBox(this);
    m_kind->addItems({tr("Impulse response"), tr("Step response")});
    m_log = new QCheckBox(tr("dB"), this);
    m_log->setToolTip(tr("Show the magnitude in dB"));
    m_controls->addWidget(m_kind);
    m_controls->addWidget(m_log);
    m_controls->addStretch(1);
    connect(m_kind, &QComboBox::currentIndexChanged, this, [this] { redrawNow(); });
    connect(m_log, &QCheckBox::toggled, this, [this] { redrawNow(); });
}

void ImpulseView::redraw() {
    m_plot->clear();
    if (!m_design) {
        m_plot->autoscale();
        return;
    }
    const int n = m_design->fir ? int(m_design->ba.b.size()) + 5 : impz_len(m_design->zpk);
    Vec x(n, 0.0);
    const bool step = m_kind->currentIndex() == 1;
    if (step) std::fill(x.begin(), x.end(), 1.0);
    else x[0] = 1.0;
    const Vec h = m_design->sos.empty() ? lfilter(m_design->ba.b, m_design->ba.a, x) : sosfilt(m_design->sos, x);
    const bool norm = m_ctx.unit_to_hz == 0;
    QVector<double> t(n), y(n);
    for (int i = 0; i < n; ++i) {
        t[i] = norm ? i : i / m_ctx.f_s;
        y[i] = m_log->isChecked() ? 20 * std::log10(std::max(std::fabs(h[i]), 1e-15)) : h[i];
    }
    m_plot->addCurve(t, y, PlotWidget::color(0), QString(),
                     m_log->isChecked() ? PlotWidget::Style::Line : PlotWidget::Style::Stem);
    m_plot->setTitle(step ? tr("Step response") : tr("Impulse response"));
    m_plot->setXLabel(m_ctx.t_label);
    m_plot->setYLabel(m_log->isChecked() ? (step ? "|s[n]| / dB" : "|h[n]| / dB") : (step ? "s[n]" : "h[n]"));
    if (m_log->isChecked()) m_plot->setYLimits(-150, std::numeric_limits<double>::quiet_NaN());
    m_plot->autoscale();
}
