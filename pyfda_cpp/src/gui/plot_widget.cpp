#include "plot_widget.hpp"
#include "settings.hpp"

#include <QApplication>
#include <QClipboard>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QMenu>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
const double NaN = std::numeric_limits<double>::quiet_NaN();

// "nice" tick spacing for the range [lo, hi] with about n ticks
double nice_step(double range, int n) {
    if (!(range > 0)) return 1.0;
    const double raw = range / std::max(n, 1);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double f = raw / mag;
    double nf;
    if (f < 1.5) nf = 1;
    else if (f < 3) nf = 2;
    else if (f < 7) nf = 5;
    else nf = 10;
    return nf * mag;
}

QString fmt_tick(double v, double step) {
    if (std::fabs(v) < step * 1e-9) v = 0.0;
    const double a = std::fabs(v);
    if (a != 0 && (a >= 1e5 || a < 1e-3)) return QString::number(v, 'g', 4);
    const int dec = std::max(0, int(-std::floor(std::log10(step) + 1e-9)));
    return QString::number(v, 'f', std::min(dec, 8));
}
}  // namespace

PlotWidget::PlotWidget(QWidget *parent) : QWidget(parent) {
    m_xlo = m_xhi = m_ylo = m_yhi = NaN;
    setMouseTracking(true);
    setMinimumSize(200, 150);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    setBackgroundRole(QPalette::Base);
    setAutoFillBackground(true);
}

QColor PlotWidget::color(int i) {
    static const QColor cycle[] = {QColor(31, 119, 180), QColor(255, 127, 14), QColor(44, 160, 44),
                                   QColor(214, 39, 40),  QColor(148, 103, 189), QColor(140, 86, 75),
                                   QColor(227, 119, 194), QColor(127, 127, 127)};
    return cycle[((i % 8) + 8) % 8];
}

void PlotWidget::setImage(const QImage &img, double x0, double x1, double y0, double y1, double zlo, double zhi,
                          const QString &zlabel) {
    m_image = img;
    m_ix0 = x0;
    m_ix1 = x1;
    m_iy0 = y0;
    m_iy1 = y1;
    m_zlo = zlo;
    m_zhi = zhi;
    m_zlabel = zlabel;
}

QRgb PlotWidget::colormap(double t) {
    // samples of matplotlib's viridis
    static const double c[][3] = {{0.267, 0.005, 0.329}, {0.283, 0.141, 0.458}, {0.254, 0.265, 0.530},
                                  {0.207, 0.372, 0.553}, {0.164, 0.471, 0.558}, {0.128, 0.567, 0.551},
                                  {0.135, 0.659, 0.518}, {0.267, 0.749, 0.441}, {0.478, 0.821, 0.318},
                                  {0.741, 0.873, 0.150}, {0.993, 0.906, 0.144}};
    if (!std::isfinite(t)) t = 0;
    t = std::clamp(t, 0.0, 1.0) * 10;
    const int i = std::min(9, int(t));
    const double f = t - i;
    auto ch = [&](int k) { return int(255 * (c[i][k] + f * (c[i + 1][k] - c[i][k])) + 0.5); };
    return qRgb(ch(0), ch(1), ch(2));
}

void PlotWidget::clear() {
    m_curves.clear();
    m_regions.clear();
    m_labels.clear();
    m_image = QImage();
    m_title.clear();
    m_message.clear();
    m_xlo = m_xhi = m_ylo = m_yhi = NaN;
    m_equal = false;
}

void PlotWidget::addCurve(const Curve &c) { m_curves.append(c); }

void PlotWidget::addCurve(const QVector<double> &x, const QVector<double> &y, const QColor &color,
                          const QString &name, Style style) {
    Curve c;
    c.x = x;
    c.y = y;
    c.color = color;
    c.name = name;
    c.style = style;
    m_curves.append(c);
}

void PlotWidget::addCircle(double cx, double cy, double r, const QColor &color) {
    Curve c;
    for (int i = 0; i <= 360; ++i) {
        const double a = i * M_PI / 180.0;
        c.x.append(cx + r * std::cos(a));
        c.y.append(cy + r * std::sin(a));
    }
    c.color = color;
    c.width = 1.0;
    c.pen = Qt::DashLine;
    m_curves.append(c);
}

void PlotWidget::autoscale() {
    if (m_keep_view && m_has_view) {
        update();
        return;
    }
    double x0 = std::numeric_limits<double>::infinity(), x1 = -x0, y0 = x0, y1 = -x0;
    for (const Curve &c : m_curves) {
        for (int i = 0; i < c.x.size() && i < c.y.size(); ++i) {
            const double x = c.x[i], y = c.y[i];
            if (!std::isfinite(x) || !std::isfinite(y)) continue;
            if (!std::isnan(m_ylo) && y < m_ylo) continue;  // limited range, e.g. -inf dB
            x0 = std::min(x0, x);
            x1 = std::max(x1, x);
            y0 = std::min(y0, y);
            y1 = std::max(y1, y);
        }
        if (c.style == Style::Stem) {
            y0 = std::min(y0, 0.0);
            y1 = std::max(y1, 0.0);
        }
    }
    if (!m_image.isNull()) {
        x0 = std::min(x0, m_ix0);
        x1 = std::max(x1, m_ix1);
        y0 = std::min(y0, m_iy0);
        y1 = std::max(y1, m_iy1);
    }
    if (!std::isfinite(x0)) {
        x0 = 0;
        x1 = 1;
        y0 = 0;
        y1 = 1;
    }
    if (!std::isnan(m_xlo)) x0 = m_xlo;
    if (!std::isnan(m_xhi)) x1 = m_xhi;
    if (!std::isnan(m_ylo)) y0 = std::max(y0, m_ylo);
    if (!std::isnan(m_yhi)) y1 = std::min(y1, m_yhi);
    if (x1 <= x0) {
        x0 -= 0.5;
        x1 += 0.5;
    }
    if (y1 <= y0) {
        const double d = std::max(std::fabs(y0) * 0.1, 1e-12);
        y0 -= d;
        y1 += d;
    }
    const double dy = m_image.isNull() ? (y1 - y0) * 0.05 : 0.0;
    m_x0 = x0;
    m_x1 = x1;
    m_y0 = y0 - dy;
    m_y1 = y1 + dy;
    if (m_equal) {
        const double dx = (x1 - x0) * 0.05;
        m_x0 -= dx;
        m_x1 += dx;
        applyEqualAspect();
    }
    m_has_view = true;
    update();
}

QRectF PlotWidget::plotRect() const {
    const QFontMetrics fm(font());
    const double left = fm.horizontalAdvance("-0.00000") + fm.height() + 8;
    const double top = m_title.isEmpty() ? 10 : fm.height() + 14;
    const double bottom = 2 * fm.height() + 12;
    const double right = m_image.isNull() ? 16 : 16 + 20 + fm.horizontalAdvance("-000.0") + fm.height() + 10;
    return QRectF(left, top, std::max(10.0, width() - left - right), std::max(10.0, height() - top - bottom));
}

void PlotWidget::applyEqualAspect() {
    const QRectF r = plotRect();
    const double sx = (m_x1 - m_x0) / r.width(), sy = (m_y1 - m_y0) / r.height();
    if (sx > sy) {
        const double c = 0.5 * (m_y0 + m_y1), h = sx * r.height() / 2;
        m_y0 = c - h;
        m_y1 = c + h;
    } else {
        const double c = 0.5 * (m_x0 + m_x1), w = sy * r.width() / 2;
        m_x0 = c - w;
        m_x1 = c + w;
    }
}

QPointF PlotWidget::toPixel(double x, double y) const {
    const QRectF r = plotRect();
    return QPointF(r.left() + (x - m_x0) / (m_x1 - m_x0) * r.width(),
                   r.bottom() - (y - m_y0) / (m_y1 - m_y0) * r.height());
}

QPointF PlotWidget::toData(const QPointF &p) const {
    const QRectF r = plotRect();
    return QPointF(m_x0 + (p.x() - r.left()) / r.width() * (m_x1 - m_x0),
                   m_y0 + (r.bottom() - p.y()) / r.height() * (m_y1 - m_y0));
}

void PlotWidget::drawCurve(QPainter &p, const Curve &c) const {
    QColor col = c.color;
    col.setAlphaF(c.alpha);
    QPen pen(col, c.width, c.pen);
    pen.setCapStyle(Qt::RoundCap);
    p.setPen(pen);
    const int n = std::min(c.x.size(), c.y.size());
    const QRectF r = plotRect();

    if (c.style == Style::Line || c.style == Style::Steps) {
        // decimation: for monotonic x with many points per pixel draw min / max per column
        bool monotonic = n > 4 * r.width();
        for (int i = 1; monotonic && i < n; ++i)
            if (c.x[i] < c.x[i - 1]) monotonic = false;
        QPainterPath path;
        bool pen_down = false;
        auto line_to = [&](const QPointF &pt) {
            if (!pen_down) path.moveTo(pt);
            else path.lineTo(pt);
            pen_down = true;
        };
        if (monotonic) {
            int col_prev = std::numeric_limits<int>::min();
            double ymin = 0, ymax = 0, yfirst = 0, ylast = 0;
            auto flush = [&](int column) {
                if (column == std::numeric_limits<int>::min()) return;
                const double px = column + 0.5;
                line_to(QPointF(px, yfirst));
                line_to(QPointF(px, ymin));
                line_to(QPointF(px, ymax));
                line_to(QPointF(px, ylast));
            };
            for (int i = 0; i < n; ++i) {
                if (!std::isfinite(c.y[i]) || !std::isfinite(c.x[i])) {
                    flush(col_prev);
                    col_prev = std::numeric_limits<int>::min();
                    pen_down = false;
                    continue;
                }
                const QPointF pt = toPixel(c.x[i], c.y[i]);
                if (pt.x() < r.left() - 2 || pt.x() > r.right() + 2) continue;
                const int column = int(std::floor(pt.x()));
                if (column != col_prev) {
                    flush(col_prev);
                    col_prev = column;
                    ymin = ymax = yfirst = ylast = pt.y();
                } else {
                    ymin = std::min(ymin, pt.y());
                    ymax = std::max(ymax, pt.y());
                    ylast = pt.y();
                }
            }
            flush(col_prev);
        } else {
            QPointF prev;
            for (int i = 0; i < n; ++i) {
                if (!std::isfinite(c.y[i]) || !std::isfinite(c.x[i])) {
                    pen_down = false;
                    continue;
                }
                const QPointF pt = toPixel(c.x[i], c.y[i]);
                if (c.style == Style::Steps && pen_down) line_to(QPointF(pt.x(), prev.y()));
                line_to(pt);
                prev = pt;
            }
        }
        p.drawPath(path);
    } else if (c.style == Style::Stem) {
        const double y0 = toPixel(0, 0).y();
        const bool draw_markers = n < r.width() / 3;
        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(c.y[i])) continue;
            const QPointF pt = toPixel(c.x[i], c.y[i]);
            p.drawLine(QPointF(pt.x(), y0), pt);
            if (draw_markers) {
                p.setBrush(col);
                p.drawEllipse(pt, 2.5, 2.5);
                p.setBrush(Qt::NoBrush);
            }
        }
    }
    if (c.style == Style::Markers || c.marker != Marker::None) {
        const Marker m = c.marker == Marker::None ? Marker::Dot : c.marker;
        p.setPen(QPen(col, 1.8));
        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(c.y[i]) || !std::isfinite(c.x[i])) continue;
            const QPointF pt = toPixel(c.x[i], c.y[i]);
            const double s = 5;
            switch (m) {
            case Marker::Circle: p.drawEllipse(pt, s, s); break;
            case Marker::Cross:
                p.drawLine(pt + QPointF(-s, -s), pt + QPointF(s, s));
                p.drawLine(pt + QPointF(-s, s), pt + QPointF(s, -s));
                break;
            default:
                p.setBrush(col);
                p.drawEllipse(pt, 2.5, 2.5);
                p.setBrush(Qt::NoBrush);
            }
        }
    }
}

void PlotWidget::paintEvent(QPaintEvent *) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    const QPalette pal = palette();
    const QColor fg = pal.color(QPalette::Text);
    QColor grid = fg;
    grid.setAlpha(45);
    const QRectF r = plotRect();
    const QFontMetrics fm(font());

    if (!m_title.isEmpty()) {
        QFont bf = font();
        bf.setBold(true);
        p.setFont(bf);
        p.setPen(fg);
        p.drawText(QRectF(0, 2, width(), fm.height() + 8), Qt::AlignCenter, m_title);
        p.setFont(font());
    }
    if (m_curves.isEmpty() && m_image.isNull() && !m_message.isEmpty()) {
        p.setPen(fg);
        p.drawRect(r);
        p.drawText(r, Qt::AlignCenter | Qt::TextWordWrap, m_message);
        return;
    }

    // ticks and grid
    const double xs = nice_step(m_x1 - m_x0, std::max(2, int(r.width() / 90)));
    const double ys = nice_step(m_y1 - m_y0, std::max(2, int(r.height() / 45)));
    p.setPen(QPen(grid, 1));
    for (double x = std::ceil(m_x0 / xs) * xs; x <= m_x1 + xs * 1e-9; x += xs) {
        const double px = toPixel(x, 0).x();
        p.setPen(QPen(grid, 1));
        p.drawLine(QPointF(px, r.top()), QPointF(px, r.bottom()));
        p.setPen(fg);
        const QString t = fmt_tick(x, xs);
        p.drawText(QPointF(px - fm.horizontalAdvance(t) / 2.0, r.bottom() + fm.ascent() + 4), t);
    }
    for (double y = std::ceil(m_y0 / ys) * ys; y <= m_y1 + ys * 1e-9; y += ys) {
        const double py = toPixel(0, y).y();
        p.setPen(QPen(grid, 1));
        p.drawLine(QPointF(r.left(), py), QPointF(r.right(), py));
        p.setPen(fg);
        const QString t = fmt_tick(y, ys);
        p.drawText(QPointF(r.left() - fm.horizontalAdvance(t) - 5, py + fm.ascent() / 2.5), t);
    }
    // axis labels
    p.setPen(fg);
    p.drawText(QRectF(r.left(), r.bottom() + fm.height() + 6, r.width(), fm.height() + 4), Qt::AlignCenter,
               m_xlabel);
    p.save();
    p.translate(fm.height() / 2.0 + 2, r.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-r.height() / 2, -fm.height() / 2.0, r.height(), fm.height() + 2), Qt::AlignCenter,
               m_ylabel);
    p.restore();

    // data
    p.save();
    p.setClipRect(r.adjusted(-1, -1, 1, 1));
    if (!m_image.isNull()) {
        const QPointF a = toPixel(m_ix0, m_iy1), b = toPixel(m_ix1, m_iy0);
        p.save();
        p.setRenderHint(QPainter::SmoothPixmapTransform, false);
        p.drawImage(QRectF(a, b).normalized(), m_image);
        p.restore();
    }
    for (const Region &g : m_regions) {
        const QPointF a = toPixel(g.x0, g.y1), b = toPixel(g.x1, g.y0);
        p.fillRect(QRectF(a, b).normalized(), g.color);
    }
    for (const Curve &c : m_curves) drawCurve(p, c);
    for (const Label &l : m_labels) {
        p.setPen(l.color);
        p.drawText(toPixel(l.x, l.y) + QPointF(6, -6), l.text);
    }
    p.restore();
    p.setPen(QPen(fg, 1));
    p.drawRect(r);

    // color bar of the image
    if (!m_image.isNull()) {
        const QRectF cb(r.right() + 12, r.top(), 16, r.height());
        QImage bar(1, 256, QImage::Format_RGB32);
        for (int i = 0; i < 256; ++i) bar.setPixel(0, 255 - i, colormap(i / 255.0));
        p.drawImage(cb, bar);
        p.setPen(QPen(fg, 1));
        p.drawRect(cb);
        const double zs = nice_step(m_zhi - m_zlo, std::max(2, int(r.height() / 45)));
        if (m_zhi > m_zlo && zs > 0)
            for (double z = std::ceil(m_zlo / zs) * zs; z <= m_zhi + zs * 1e-9; z += zs) {
                const double py = cb.bottom() - (z - m_zlo) / (m_zhi - m_zlo) * cb.height();
                p.drawLine(QPointF(cb.right(), py), QPointF(cb.right() + 3, py));
                p.drawText(QPointF(cb.right() + 5, py + fm.ascent() / 2.5), fmt_tick(z, zs));
            }
        p.save();
        p.translate(width() - fm.height() / 2.0 - 2, cb.center().y());
        p.rotate(-90);
        p.drawText(QRectF(-cb.height() / 2, -fm.height() / 2.0, cb.height(), fm.height() + 2), Qt::AlignCenter,
                   m_zlabel);
        p.restore();
    }

    // legend
    QStringList names;
    for (const Curve &c : m_curves)
        if (!c.name.isEmpty()) names << c.name;
    if (!names.isEmpty()) {
        int w = 0;
        for (const QString &s : names) w = std::max(w, fm.horizontalAdvance(s));
        const QRectF lr(r.right() - w - 44, r.top() + 6, w + 38, names.size() * (fm.height() + 2) + 8);
        QColor bg = pal.color(QPalette::Base);
        bg.setAlpha(220);
        p.fillRect(lr, bg);
        p.setPen(QPen(grid, 1));
        p.drawRect(lr);
        double y = lr.top() + 4;
        for (const Curve &c : m_curves) {
            if (c.name.isEmpty()) continue;
            const double cy = y + fm.height() / 2.0;
            p.setPen(QPen(c.color, 2.0, c.pen));  // same line style as the curve (e.g. dashed)
            p.drawLine(QPointF(lr.left() + 4, cy), QPointF(lr.left() + 28, cy));
            p.setPen(fg);
            p.drawText(QPointF(lr.left() + 32, y + fm.ascent()), c.name);
            y += fm.height() + 2;
        }
    }

    // rubber band
    if (m_rubber) {
        QColor sel = pal.color(QPalette::Highlight);
        sel.setAlpha(60);
        p.fillRect(QRectF(m_press, m_current).normalized(), sel);
        p.setPen(QPen(pal.color(QPalette::Highlight), 1, Qt::DashLine));
        p.drawRect(QRectF(m_press, m_current).normalized());
    }
    // cursor coordinates
    if (m_mouse_inside && r.contains(m_mouse)) {
        const QPointF d = toData(m_mouse);
        const QString t = QString("x = %1   y = %2").arg(d.x(), 0, 'g', 6).arg(d.y(), 0, 'g', 6);
        p.setPen(fg);
        p.drawText(QPointF(r.left() + 6, r.bottom() - 6), t);
    }
}

void PlotWidget::mousePressEvent(QMouseEvent *e) {
    if (e->button() == Qt::LeftButton && plotRect().contains(e->position())) {
        m_rubber = true;
        m_press = m_current = e->position();
    } else if (e->button() == Qt::MiddleButton ||
               (e->button() == Qt::LeftButton && (e->modifiers() & Qt::ControlModifier))) {
        m_panning = true;
        m_press = e->position();
        m_px0 = m_x0;
        m_px1 = m_x1;
        m_py0 = m_y0;
        m_py1 = m_y1;
    }
}

void PlotWidget::mouseMoveEvent(QMouseEvent *e) {
    m_mouse = e->position();
    m_mouse_inside = true;
    if (m_rubber) {
        m_current = e->position();
    } else if (m_panning) {
        const QRectF r = plotRect();
        const double dx = (e->position().x() - m_press.x()) / r.width() * (m_px1 - m_px0);
        const double dy = (e->position().y() - m_press.y()) / r.height() * (m_py1 - m_py0);
        m_x0 = m_px0 - dx;
        m_x1 = m_px1 - dx;
        m_y0 = m_py0 + dy;
        m_y1 = m_py1 + dy;
        emit viewChanged();
    }
    update();
}

void PlotWidget::mouseReleaseEvent(QMouseEvent *e) {
    if (m_rubber && e->button() == Qt::LeftButton) {
        m_rubber = false;
        const QRectF sel = QRectF(m_press, e->position()).normalized();
        if (sel.width() > 5 && sel.height() > 5) {
            const QPointF a = toData(sel.bottomLeft()), b = toData(sel.topRight());
            m_x0 = a.x();
            m_x1 = b.x();
            m_y0 = a.y();
            m_y1 = b.y();
            if (m_equal) applyEqualAspect();
            emit viewChanged();
        }
        update();
    }
    m_panning = false;
}

void PlotWidget::mouseDoubleClickEvent(QMouseEvent *) {
    const bool keep = m_keep_view;
    m_keep_view = false;
    autoscale();
    m_keep_view = keep;
    emit viewChanged();
}

void PlotWidget::wheelEvent(QWheelEvent *e) {
    const QRectF r = plotRect();
    if (!r.contains(e->position())) return;
    const double f = std::pow(0.85, e->angleDelta().y() / 120.0);
    const QPointF c = toData(e->position());
    const bool only_x = e->modifiers() & Qt::ShiftModifier;
    const bool only_y = e->modifiers() & Qt::ControlModifier;
    if (!only_y) {
        m_x0 = c.x() + (m_x0 - c.x()) * f;
        m_x1 = c.x() + (m_x1 - c.x()) * f;
    }
    if (!only_x) {
        m_y0 = c.y() + (m_y0 - c.y()) * f;
        m_y1 = c.y() + (m_y1 - c.y()) * f;
    }
    emit viewChanged();
    update();
}

void PlotWidget::leaveEvent(QEvent *) {
    m_mouse_inside = false;
    update();
}

void PlotWidget::contextMenuEvent(QContextMenuEvent *e) {
    QMenu menu(this);
    QAction *reset = menu.addAction(tr("Reset view (double click)"));
    QAction *copy = menu.addAction(tr("Copy image"));
    QAction *save = menu.addAction(tr("Save image ..."));
    QAction *a = menu.exec(e->globalPos());
    if (a == reset) {
        mouseDoubleClickEvent(nullptr);
    } else if (a == copy) {
        QApplication::clipboard()->setPixmap(grab());
    } else if (a == save) {
        const QString f = QFileDialog::getSaveFileName(this, tr("Save image"), config::dir("export"),
                                                       tr("PNG image (*.png)"));
        if (f.isEmpty()) return;
        config::setDir("export", f);
        savePng(f);
    }
}

bool PlotWidget::savePng(const QString &file) { return grab().save(file, "PNG"); }
