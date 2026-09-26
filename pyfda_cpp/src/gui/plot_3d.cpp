#include "plot_3d.hpp"

#include "filter_info.hpp"
#include "plot_widget.hpp"
#include "settings.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QVBoxLayout>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pyfda;

namespace {
constexpr double DEG = PI / 180;
constexpr double Z_HEIGHT = 1.3;  // height of the box relative to its width (normalized units)

// "nice" tick positions like the 2D plots
QVector<double> ticks(double lo, double hi, int n) {
    QVector<double> t;
    if (!(hi > lo)) return t;
    const double raw = (hi - lo) / n, mag = std::pow(10.0, std::floor(std::log10(raw)));
    double step = mag;
    for (double m : {1.0, 2.0, 2.5, 5.0, 10.0})
        if (m * mag >= raw) {
            step = m * mag;
            break;
        }
    for (double v = std::ceil(lo / step - 1e-9) * step; v <= hi + 1e-9 * step; v += step)
        t << (std::fabs(v) < 1e-12 * step ? 0.0 : v);
    return t;
}

QColor shade(QColor c, double f) {
    return QColor::fromRgbF(float(std::clamp(c.redF() * f, 0.0, 1.0)), float(std::clamp(c.greenF() * f, 0.0, 1.0)),
                            float(std::clamp(c.blueF() * f, 0.0, 1.0)));
}
}  // namespace

// ---------------------------------------------------------------------------
// Plot3D

Plot3D::Plot3D(QWidget *parent) : QWidget(parent) {
    setMinimumSize(300, 250);
    setMouseTracking(false);
    setToolTip(tr("Drag to rotate, wheel to zoom, double click to reset the view"));
}

void Plot3D::setSurface(const QVector<P3> &grid, int nu, int nv, bool filled, bool wrap_v) {
    m_grid = grid;
    m_nu = nu;
    m_nv = nv;
    m_filled = filled;
    m_wrap = wrap_v;
}

void Plot3D::clearSurface() {
    m_grid.clear();
    m_nu = m_nv = 0;
}

void Plot3D::setRange(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                      const QString &zlabel) {
    m_xmin = xmin;
    m_xmax = xmax;
    m_ymin = ymin;
    m_ymax = ymax;
    m_zmin = zmin;
    m_zmax = zmax > zmin ? zmax : zmin + 1;
    m_zlabel = zlabel;
}

void Plot3D::resetView() {
    m_azim = -60;
    m_elev = 30;
    m_zoom = 1;
    update();
}

Plot3D::R3 Plot3D::rotate(const P3 &p) const {
    const double xn = (2 * p.x - m_xmin - m_xmax) / (m_xmax - m_xmin);
    const double yn = (2 * p.y - m_ymin - m_ymax) / (m_ymax - m_ymin);
    const double zn = ((p.z - m_zmin) / (m_zmax - m_zmin) - 0.5) * Z_HEIGHT;
    const double a = m_azim * DEG, e = m_elev * DEG;
    // azimuth around the z axis, then tilt by the elevation (camera above the xy plane)
    const double x1 = xn * std::cos(a) + yn * std::sin(a);
    const double y1 = -xn * std::sin(a) + yn * std::cos(a);  // away from the viewer
    return {x1, zn * std::cos(e) + y1 * std::sin(e), -y1 * std::cos(e) + zn * std::sin(e)};
}

void Plot3D::paintEvent(QPaintEvent *) {
    QPainter p(this);
    const QPalette pal = palette();
    p.fillRect(rect(), pal.color(QPalette::Base));
    const QColor fg = pal.color(QPalette::Text);
    const QFontMetrics fm = p.fontMetrics();
    const bool cbar = m_colorbar && m_filled && !m_grid.isEmpty();
    // plot area without the title and the color bar
    const QRectF area(0, fm.height() + 8, width() - (cbar ? 90 : 0), height() - fm.height() - 8);
    const double s = std::min(area.width(), area.height()) * 0.36 * m_zoom;
    const QPointF c = area.center();
    auto scr = [&](const R3 &r) { return QPointF(c.x() + r.sx * s, c.y() - r.sy * s); };
    auto proj = [&](const P3 &q) { return scr(rotate(q)); };

    if (!m_title.isEmpty()) {
        QFont f = p.font();
        f.setBold(true);
        p.setFont(f);
        p.setPen(fg);
        p.drawText(QRectF(0, 4, width(), fm.height() + 4), Qt::AlignHCenter, m_title);
        p.setFont(font());
    }

    // --- axes box: floor with grid, back walls, ticks
    const P3 corners[4] = {{m_xmin, m_ymin, 0}, {m_xmax, m_ymin, 0}, {m_xmax, m_ymax, 0}, {m_xmin, m_ymax, 0}};
    auto at = [](P3 q, double z) {
        q.z = z;
        return q;
    };
    QColor grid = fg;
    grid.setAlpha(40);
    QColor box = fg;
    box.setAlpha(110);
    p.setRenderHint(QPainter::Antialiasing, true);
    // back walls: the two vertical sides whose center is farther away
    for (int i = 0; i < 4; ++i) {
        const P3 a = corners[i], b = corners[(i + 1) % 4];
        const P3 mid{(a.x + b.x) / 2, (a.y + b.y) / 2, (m_zmin + m_zmax) / 2};
        if (rotate(mid).depth > rotate({(m_xmin + m_xmax) / 2, (m_ymin + m_ymax) / 2, mid.z}).depth) continue;
        p.setPen(QPen(grid, 1));
        for (double z : ticks(m_zmin, m_zmax, 5)) p.drawLine(proj(at(a, z)), proj(at(b, z)));
        p.setPen(QPen(box, 1));
        p.drawLine(proj(at(a, m_zmin)), proj(at(a, m_zmax)));
        p.drawLine(proj(at(b, m_zmin)), proj(at(b, m_zmax)));
        p.drawLine(proj(at(a, m_zmax)), proj(at(b, m_zmax)));
    }
    // floor
    p.setPen(QPen(grid, 1));
    for (double x : ticks(m_xmin, m_xmax, 6)) p.drawLine(proj({x, m_ymin, m_zmin}), proj({x, m_ymax, m_zmin}));
    for (double y : ticks(m_ymin, m_ymax, 6)) p.drawLine(proj({m_xmin, y, m_zmin}), proj({m_xmax, y, m_zmin}));
    p.setPen(QPen(box, 1));
    for (int i = 0; i < 4; ++i) p.drawLine(proj(at(corners[i], m_zmin)), proj(at(corners[(i + 1) % 4], m_zmin)));
    // tick labels: x along the front y edge, y along the front x edge, z at the leftmost corner
    const double y_front = rotate({0, m_ymin, m_zmin}).depth > rotate({0, m_ymax, m_zmin}).depth ? m_ymin : m_ymax;
    const double x_front = rotate({m_xmin, 0, m_zmin}).depth > rotate({m_xmax, 0, m_zmin}).depth ? m_xmin : m_xmax;
    p.setPen(fg);
    auto label = [&](const QPointF &pt, const QPointF &out, const QString &t) {
        const double w = fm.horizontalAdvance(t);
        p.drawText(QPointF(pt.x() + out.x() - w / 2, pt.y() + out.y() + fm.ascent() / 2.0 - 1), t);
    };
    auto outward = [&](const P3 &q, const P3 &center) {
        QPointF d = proj(q) - proj(center);
        const double l = std::hypot(d.x(), d.y());
        return l > 0 ? d * (18.0 / l) : QPointF(0, 18);
    };
    for (double x : ticks(m_xmin, m_xmax, 6)) {
        const P3 q{x, y_front, m_zmin};
        label(proj(q), outward(q, {x, (m_ymin + m_ymax) / 2, m_zmin}), QString::number(x, 'g', 3));
    }
    for (double y : ticks(m_ymin, m_ymax, 6)) {
        const P3 q{x_front, y, m_zmin};
        label(proj(q), outward(q, {(m_xmin + m_xmax) / 2, y, m_zmin}), QString::number(y, 'g', 3));
    }
    {
        const P3 qx{(m_xmin + m_xmax) / 2, y_front, m_zmin}, qy{x_front, (m_ymin + m_ymax) / 2, m_zmin};
        label(proj(qx), outward(qx, {qx.x, (m_ymin + m_ymax) / 2, m_zmin}) * 2.3, "Re");
        label(proj(qy), outward(qy, {(m_xmin + m_xmax) / 2, qy.y, m_zmin}) * 2.3, "Im");
    }
    int left = 0;
    for (int i = 1; i < 4; ++i)
        if (proj(corners[i]).x() < proj(corners[left]).x()) left = i;
    for (double z : ticks(m_zmin, m_zmax, 5)) {
        const QPointF pt = proj(at(corners[left], z));
        const QString t = QString::number(z, 'g', 4);
        p.drawText(QPointF(pt.x() - fm.horizontalAdvance(t) - 8, pt.y() + fm.ascent() / 2.0 - 1), t);
    }
    {
        const QPointF top = proj(at(corners[left], m_zmax));
        p.drawText(QPointF(top.x() - fm.horizontalAdvance(m_zlabel) / 2.0, top.y() - 8), m_zlabel);
    }

    // --- depth sorted primitives: surface quads, line segments, markers
    struct Prim {
        double depth;
        int kind;  // 0 quad, 1 segment, 2 marker
        int a, b;  // quad: u, v; segment: line, index; marker: index
    };
    std::vector<Prim> prims;
    QVector<R3> rg(m_grid.size());
    for (int i = 0; i < m_grid.size(); ++i) rg[i] = rotate(m_grid[i]);
    const int nvq = m_wrap ? m_nv : m_nv - 1;
    for (int u = 0; u + 1 < m_nu; ++u)
        for (int v = 0; v < nvq; ++v) {
            const int v1 = (v + 1) % m_nv;
            const double d = (rg[u * m_nv + v].depth + rg[(u + 1) * m_nv + v].depth + rg[(u + 1) * m_nv + v1].depth +
                              rg[u * m_nv + v1].depth) / 4;
            prims.push_back({d, 0, u, v});
        }
    for (int l = 0; l < m_lines.size(); ++l)
        for (int i = 0; i + 1 < m_lines[l].pts.size(); ++i) {
            // lines lie on or above the surface: small bias towards the viewer
            const double d = (rotate(m_lines[l].pts[i]).depth + rotate(m_lines[l].pts[i + 1]).depth) / 2 + 0.01;
            prims.push_back({d, 1, l, i});
        }
    for (int i = 0; i < m_markers.size(); ++i) prims.push_back({rotate(m_markers[i].p).depth + 0.02, 2, i, 0});
    std::sort(prims.begin(), prims.end(), [](const Prim &x, const Prim &y) { return x.depth < y.depth; });

    // light from the viewer's upper left for the shading of the surface
    const double la = (m_azim + 45) * DEG;
    const double lx = -std::sin(la) * 0.5, ly = -std::cos(la) * 0.5, lz = 0.8;
    const double ll = std::sqrt(lx * lx + ly * ly + lz * lz);
    QColor mesh = fg;
    mesh.setAlpha(90);
    for (const Prim &pr : prims) {
        if (pr.kind == 0) {
            const int u = pr.a, v = pr.b, v1 = (v + 1) % m_nv;
            const int idx[4] = {u * m_nv + v, (u + 1) * m_nv + v, (u + 1) * m_nv + v1, u * m_nv + v1};
            QPointF poly[4];
            double zs = 0;
            for (int k = 0; k < 4; ++k) {
                poly[k] = scr(rg[idx[k]]);
                zs += m_grid[idx[k]].z;
            }
            if (m_filled) {
                // normal in normalized coordinates for a simple Lambert shading
                auto nrm = [&](int k) {
                    const P3 &q = m_grid[idx[k]];
                    return P3{(2 * q.x - m_xmin - m_xmax) / (m_xmax - m_xmin), (2 * q.y - m_ymin - m_ymax) / (m_ymax - m_ymin),
                              (q.z - m_zmin) / (m_zmax - m_zmin) * Z_HEIGHT};
                };
                const P3 p0 = nrm(0), p1 = nrm(1), p2 = nrm(2), p3 = nrm(3);
                const P3 d1{p2.x - p0.x, p2.y - p0.y, p2.z - p0.z}, d2{p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};
                P3 n{d1.y * d2.z - d1.z * d2.y, d1.z * d2.x - d1.x * d2.z, d1.x * d2.y - d1.y * d2.x};
                const double nl = std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
                double f = 1.0;
                if (nl > 0) f = 0.55 + 0.5 * std::fabs(n.x * lx + n.y * ly + n.z * lz) / (nl * ll);
                const QColor col = shade(QColor::fromRgb(PlotWidget::colormap((zs / 4 - m_zmin) / (m_zmax - m_zmin))), f);
                p.setRenderHint(QPainter::Antialiasing, false);  // no seams between the quads
                p.setPen(QPen(col, 1));
                p.setBrush(col);
                p.drawPolygon(poly, 4);
            } else {
                p.setRenderHint(QPainter::Antialiasing, true);
                p.setPen(QPen(mesh, 1));
                p.setBrush(Qt::NoBrush);
                p.drawPolygon(poly, 4);
            }
        } else if (pr.kind == 1) {
            const Line &l = m_lines[pr.a];
            p.setRenderHint(QPainter::Antialiasing, true);
            p.setPen(QPen(l.color, l.width, l.pen, Qt::RoundCap));
            p.drawLine(proj(l.pts[pr.b]), proj(l.pts[pr.b + 1]));
        } else {
            const Marker &mk = m_markers[pr.a];
            const QPointF q = proj(mk.p);
            p.setRenderHint(QPainter::Antialiasing, true);
            p.setPen(QPen(mk.color, 2));
            p.setBrush(Qt::NoBrush);
            if (mk.cross) {
                p.drawLine(q + QPointF(-5, -5), q + QPointF(5, 5));
                p.drawLine(q + QPointF(-5, 5), q + QPointF(5, -5));
            } else {
                p.drawEllipse(q, 5, 5);
            }
        }
    }

    // --- color bar
    if (cbar) {
        p.setRenderHint(QPainter::Antialiasing, false);
        const QRectF bar(width() - 80, area.top() + area.height() * 0.15, 16, area.height() * 0.7);
        for (int i = 0; i < int(bar.height()); ++i)
            p.fillRect(QRectF(bar.left(), bar.bottom() - i - 1, bar.width(), 1),
                       QColor::fromRgb(PlotWidget::colormap(i / (bar.height() - 1))));
        p.setPen(QPen(box, 1));
        p.setBrush(Qt::NoBrush);
        p.drawRect(bar);
        p.setPen(fg);
        for (double z : ticks(m_zmin, m_zmax, 5)) {
            const double y = bar.bottom() - (z - m_zmin) / (m_zmax - m_zmin) * bar.height();
            p.drawLine(QPointF(bar.right(), y), QPointF(bar.right() + 4, y));
            p.drawText(QPointF(bar.right() + 7, y + fm.ascent() / 2.0 - 1), QString::number(z, 'g', 4));
        }
    }
}

void Plot3D::mousePressEvent(QMouseEvent *e) { m_last = e->pos(); }

void Plot3D::mouseMoveEvent(QMouseEvent *e) {
    if (!(e->buttons() & Qt::LeftButton)) return;
    const QPoint d = e->pos() - m_last;
    m_last = e->pos();
    m_azim -= d.x() * 0.5;
    m_elev = std::clamp(m_elev + d.y() * 0.5, -90.0, 90.0);
    update();
}

void Plot3D::mouseDoubleClickEvent(QMouseEvent *) { resetView(); }

void Plot3D::wheelEvent(QWheelEvent *e) {
    m_zoom = std::clamp(m_zoom * std::pow(1.0015, e->angleDelta().y()), 0.3, 5.0);
    update();
}

void Plot3D::contextMenuEvent(QContextMenuEvent *e) {
    QMenu menu(this);
    QAction *reset = menu.addAction(tr("Reset view (double click)"));
    QAction *copy = menu.addAction(tr("Copy image"));
    QAction *save = menu.addAction(tr("Save image ..."));
    QAction *a = menu.exec(e->globalPos());
    if (a == reset) {
        resetView();
    } else if (a == copy) {
        QApplication::clipboard()->setPixmap(grab());
    } else if (a == save) {
        const QString f = QFileDialog::getSaveFileName(this, tr("Save image"), config::dir("export"),
                                                       tr("PNG image (*.png)"));
        if (f.isEmpty()) return;
        config::setDir("export", f);
        grab().save(f, "PNG");
    }
}

// ---------------------------------------------------------------------------
// ThreeDView

ThreeDView::ThreeDView(QWidget *parent) : DesignView(parent) {
    auto *lay = new QVBoxLayout(this);
    lay->setContentsMargins(4, 4, 4, 4);
    auto *row = new QHBoxLayout();
    m_mode = new QComboBox(this);
    m_mode->addItems({tr("Surface"), tr("Mesh"), tr("None")});
    m_mode->setToolTip(tr("Display |H(z)| as colored surface, as wire mesh or not at all"));
    m_log = new QCheckBox("dB", this);
    m_log->setToolTip(tr("Logarithmic scale"));
    m_bottom = new QLineEdit("0", this);
    m_top = new QLineEdit("4", this);
    for (QLineEdit *e : {m_bottom, m_top}) e->setMaximumWidth(60);
    m_bottom->setToolTip(tr("Lower limit of the display, the surface is clipped there"));
    m_top->setToolTip(tr("Upper limit of the display, the surface is clipped there"));
    m_polar = new QCheckBox("|z| < 1", this);
    m_polar->setChecked(true);
    m_polar->setToolTip(tr("Plot |H(z)| only inside the unit circle (polar grid), else for |Re|, |Im| < 1.5"));
    m_uc = new QCheckBox("UC", this);
    m_uc->setChecked(true);
    m_uc->setToolTip(tr("Show the unit circle"));
    m_hf = new QCheckBox("H(f)", this);
    m_hf->setChecked(true);
    m_hf->setToolTip(tr("Show |H(f)| along the unit circle"));
    m_pz = new QCheckBox("P/Z", this);
    m_pz->setChecked(true);
    m_pz->setToolTip(tr("Show poles (x) and zeros (o)"));
    m_cbar = new QCheckBox(tr("Colorbar"), this);
    m_cbar->setChecked(true);
    row->addWidget(m_mode);
    row->addWidget(m_log);
    row->addWidget(new QLabel(tr("Bottom:"), this));
    row->addWidget(m_bottom);
    row->addWidget(new QLabel(tr("Top:"), this));
    row->addWidget(m_top);
    row->addSpacing(8);
    row->addWidget(m_polar);
    row->addWidget(m_uc);
    row->addWidget(m_hf);
    row->addWidget(m_pz);
    row->addWidget(m_cbar);
    row->addStretch(1);
    lay->addLayout(row);
    m_plot = new Plot3D(this);
    lay->addWidget(m_plot, 1);

    for (QCheckBox *c : {m_polar, m_uc, m_hf, m_pz, m_cbar})
        connect(c, &QCheckBox::toggled, this, [this] { redrawNow(); });
    connect(m_mode, &QComboBox::currentIndexChanged, this, [this] { redrawNow(); });
    for (QLineEdit *e : {m_bottom, m_top}) connect(e, &QLineEdit::editingFinished, this, [this] { redrawNow(); });
    // like pyfda: switch the limits between linear and dB
    connect(m_log, &QCheckBox::toggled, this, [this](bool on) {
        bool ok_t;
        const double t = m_top->text().toDouble(&ok_t);
        if (on) {
            m_bottom->setText("-80");
            m_top->setText(QString::number(ok_t && t > 0 ? std::round(2000 * std::log10(t)) / 100 : 12.04));
        } else {
            m_bottom->setText("0");
            m_top->setText(QString::number(ok_t ? std::round(100 * std::pow(10.0, t / 20)) / 100 : 4));
        }
        redrawNow();
    });
}

void ThreeDView::redraw() {
    m_plot->clearSurface();
    m_plot->setLines({});
    m_plot->setMarkers({});
    m_plot->setTitle(tr("|H(z)| and |H(f)| along the unit circle"));
    if (!m_design) {
        m_plot->update();
        return;
    }
    const Ba &ba = m_design->ba;
    const bool log = m_log->isChecked();
    auto num = [](const QLineEdit *e, double def) {
        bool ok;
        const double v = e->text().toDouble(&ok);
        return ok && std::isfinite(v) ? v : def;
    };
    double top = num(m_top, log ? 12.04 : 4), bottom_set = num(m_bottom, log ? -80 : 0);

    // |H| along the unit circle
    const int n_uc = 400;
    QVector<cplx> uc;
    double h_min = std::numeric_limits<double>::infinity();
    Vec h_uc;
    for (int i = 0; i <= n_uc; ++i) {
        uc << std::polar(1.0, 2 * PI * i / n_uc);
        h_uc.push_back(h_mag_z(ba, uc.back()));
        if (std::isfinite(h_uc.back())) h_min = std::min(h_min, h_uc.back());
    }
    // display limits like pyfda: the bottom is raised to the minimum of |H(f)|
    double bottom;
    auto val = [&](double m) { return log ? 20 * std::log10(std::max(m, 1e-300)) : m; };
    if (log) bottom = std::floor(std::max(bottom_set, std::isfinite(h_min) ? val(h_min) : bottom_set) / 10) * 10;
    else bottom = std::max(bottom_set, std::isfinite(h_min) ? h_min : bottom_set);
    if (!(top > bottom)) top = bottom + (log ? 10 : 1);
    auto clip = [&](double m) { return std::clamp(std::isfinite(m) ? val(m) : top, bottom, top); };

    // surface
    const int mode = m_mode->currentIndex();
    const bool polar = m_polar->isChecked();
    if (mode < 2) {
        QVector<Plot3D::P3> grid;
        int nu, nv;
        if (polar) {  // radius x angle
            nu = 41;
            nv = 120;
            for (int u = 0; u < nu; ++u)
                for (int v = 0; v < nv; ++v) {
                    const cplx z = std::polar(double(u) / (nu - 1), 2 * PI * v / nv);
                    grid << Plot3D::P3{z.real(), z.imag(), clip(h_mag_z(ba, z))};
                }
        } else {
            nu = nv = 91;
            for (int u = 0; u < nu; ++u)
                for (int v = 0; v < nv; ++v) {
                    const cplx z(-1.5 + 3.0 * v / (nv - 1), -1.5 + 3.0 * u / (nu - 1));
                    grid << Plot3D::P3{z.real(), z.imag(), clip(h_mag_z(ba, z))};
                }
        }
        m_plot->setSurface(grid, nu, nv, mode == 0, polar);
    }

    QVector<Plot3D::Line> lines;
    if (m_uc->isChecked()) {
        Plot3D::Line l;
        for (const cplx &z : uc) l.pts << Plot3D::P3{z.real(), z.imag(), bottom};
        l.color = palette().color(QPalette::Text);
        l.width = 1.5;
        lines << l;
    }
    if (m_hf->isChecked()) {
        Plot3D::Line l;
        for (int i = 0; i <= n_uc; ++i) l.pts << Plot3D::P3{uc[i].real(), uc[i].imag(), clip(h_uc[size_t(i)])};
        l.color = PlotWidget::color(0);
        l.width = 3;
        lines << l;
        // thin vertical lines from the bottom to |H(f)|
        for (int i = 0; i < n_uc; i += 10) {
            Plot3D::Line v;
            v.pts << Plot3D::P3{uc[i].real(), uc[i].imag(), bottom} << Plot3D::P3{uc[i].real(), uc[i].imag(), clip(h_uc[size_t(i)])};
            v.color = QColor(128, 128, 128);
            lines << v;
        }
    }
    QVector<Plot3D::Marker> markers;
    if (m_pz->isChecked()) {
        const double span = top - bottom;
        const double zlevel = bottom + 0.1 * span;                  // height of the zero markers
        const double plevel_btm = mode == 2 ? bottom : top;         // pole stems start at the top
        const double plevel_top = mode == 2 ? bottom + 0.1 * span : top + 0.05 * span;
        auto inside = [](const cplx &z) { return std::fabs(z.real()) <= 1.5 && std::fabs(z.imag()) <= 1.5; };
        for (const cplx &z : m_design->zpk.z) {
            if (!inside(z)) continue;
            Plot3D::Line s;
            s.pts << Plot3D::P3{z.real(), z.imag(), bottom} << Plot3D::P3{z.real(), z.imag(), zlevel};
            s.color = QColor(0, 0, 255);
            lines << s;
            markers << Plot3D::Marker{{z.real(), z.imag(), zlevel}, QColor(0, 0, 255), false};
        }
        for (const cplx &z : m_design->zpk.p) {
            if (!inside(z)) continue;
            Plot3D::Line s;
            s.pts << Plot3D::P3{z.real(), z.imag(), plevel_btm} << Plot3D::P3{z.real(), z.imag(), plevel_top};
            s.color = QColor(220, 0, 0);
            lines << s;
            markers << Plot3D::Marker{{z.real(), z.imag(), plevel_top}, QColor(220, 0, 0), true};
        }
    }
    m_plot->setLines(lines);
    m_plot->setMarkers(markers);
    m_plot->setColorbar(m_cbar->isChecked());
    // the box includes the pole markers above the top
    const double zmax = m_pz->isChecked() && mode != 2 && !m_design->zpk.p.empty() ? top + 0.05 * (top - bottom) : top;
    m_plot->setRange(-1.5, 1.5, -1.5, 1.5, bottom, zmax, log ? "|H| / dB" : "|H|");
    m_plot->update();
}
