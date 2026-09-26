// Lightweight 2D plot widget drawn with QPainter (replaces matplotlib in pyfda).
// Line / stem / marker curves, grid, legend, rubber band and wheel zoom,
// panning with the middle mouse button or Ctrl + drag, double click to reset the view,
// equal aspect ratio for pole / zero plots. Long curves are decimated per pixel.
#pragma once

#include <QColor>
#include <QImage>
#include <QPointF>
#include <QString>
#include <QVector>
#include <QWidget>

class PlotWidget : public QWidget {
    Q_OBJECT
public:
    enum class Style { Line, Stem, Markers, Steps };
    enum class Marker { None, Circle, Cross, Dot };

    struct Curve {
        QVector<double> x, y;
        QColor color;
        QString name;
        Style style = Style::Line;
        Marker marker = Marker::None;
        double width = 1.5;
        Qt::PenStyle pen = Qt::SolidLine;
        double alpha = 1.0;
    };
    struct Region {  // shaded rectangle, e.g. filter specs
        double x0, x1, y0, y1;
        QColor color;
    };
    struct Label {  // text at a data position
        double x, y;
        QString text;
        QColor color;
    };

    explicit PlotWidget(QWidget *parent = nullptr);

    void clear();
    void addCurve(const Curve &c);
    void addCurve(const QVector<double> &x, const QVector<double> &y, const QColor &color,
                  const QString &name = QString(), Style style = Style::Line);
    void addRegion(const Region &r) { m_regions.append(r); }
    void addLabel(const Label &l) { m_labels.append(l); }
    void addCircle(double cx, double cy, double r, const QColor &color);  // drawn as a curve
    /// Image in data coordinates (e.g. spectrogram), row 0 at y1, with a color bar
    /// for the values zlo ... zhi
    void setImage(const QImage &img, double x0, double x1, double y0, double y1, double zlo, double zhi,
                  const QString &zlabel);
    /// Viridis like color map, t = 0 ... 1
    static QRgb colormap(double t);
    void setTitle(const QString &t) { m_title = t; }
    void setXLabel(const QString &t) { m_xlabel = t; }
    void setYLabel(const QString &t) { m_ylabel = t; }
    void setMessage(const QString &t) { m_message = t; }  // shown when there are no curves
    void setEqualAspect(bool on) { m_equal = on; }
    /// Limits for autoscaling, NaN = automatic
    void setYLimits(double lo, double hi) { m_ylo = lo; m_yhi = hi; }
    void setXLimits(double lo, double hi) { m_xlo = lo; m_xhi = hi; }
    /// Recalculate the autoscaled view and repaint
    void autoscale();
    /// Keep the current zoom when replacing curves (e.g. after refiltering)
    void keepView(bool keep) { m_keep_view = keep; }
    bool savePng(const QString &file);

    static QColor color(int i);  // default color cycle

signals:
    void viewChanged();

protected:
    void paintEvent(QPaintEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;
    void mouseDoubleClickEvent(QMouseEvent *) override;
    void wheelEvent(QWheelEvent *) override;
    void leaveEvent(QEvent *) override;
    void contextMenuEvent(QContextMenuEvent *) override;

private:
    QRectF plotRect() const;
    QPointF toPixel(double x, double y) const;
    QPointF toData(const QPointF &p) const;
    void applyEqualAspect();
    void drawCurve(class QPainter &p, const Curve &c) const;

    QVector<Curve> m_curves;
    QVector<Region> m_regions;
    QVector<Label> m_labels;
    QImage m_image;
    double m_ix0 = 0, m_ix1 = 1, m_iy0 = 0, m_iy1 = 1, m_zlo = 0, m_zhi = 1;
    QString m_zlabel;
    QString m_title, m_xlabel, m_ylabel, m_message;
    bool m_equal = false;
    bool m_keep_view = false;
    bool m_has_view = false;
    double m_xlo, m_xhi, m_ylo, m_yhi;  // autoscale limits (NaN = auto)
    double m_x0 = 0, m_x1 = 1, m_y0 = 0, m_y1 = 1;  // current view
    // mouse interaction
    bool m_rubber = false, m_panning = false;
    QPointF m_press, m_current, m_mouse;
    double m_px0, m_px1, m_py0, m_py1;
    bool m_mouse_inside = false;
};
