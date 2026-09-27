// 3D plot of |H(z)| over the z-plane (pyfda's "3D" tab, plot_widgets/plot_3d.py):
// surface or mesh of the magnitude, |H(f)| along the unit circle, poles and zeros
// with stems. Plot3D is a small QPainter based 3D widget (orthographic projection,
// painter's algorithm), rotated with the mouse and zoomed with the wheel.
#pragma once

#include "design_view.hpp"

#include <QColor>
#include <QVector>
#include <QWidget>

class QCheckBox;
class QComboBox;
class QLineEdit;

class Plot3D : public QWidget {
    Q_OBJECT
public:
    struct P3 {
        double x, y, z;
    };
    struct Line {
        QVector<P3> pts;
        QColor color;
        double width = 1.0;
        Qt::PenStyle pen = Qt::SolidLine;
    };
    struct Marker {
        P3 p;
        QColor color;
        bool cross = false;  // 'x' (pole) or 'o' (zero)
    };

    explicit Plot3D(QWidget *parent = nullptr);
    /// Grid of nu x nv points (row-major, index u * nv + v), z already clipped;
    /// `filled`: colored surface, else a gray wire mesh; wrap_v: v is periodic (polar grid)
    void setSurface(const QVector<P3> &grid, int nu, int nv, bool filled, bool wrap_v);
    void clearSurface();
    void setLines(const QVector<Line> &lines) { m_lines = lines; }
    void setMarkers(const QVector<Marker> &markers) { m_markers = markers; }
    /// Axis ranges used for scaling and the box, z label
    void setRange(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, const QString &zlabel);
    void setTitle(const QString &t) { m_title = t; }
    void setColorbar(bool on) { m_colorbar = on; }
    void resetView();

protected:
    void paintEvent(QPaintEvent *) override;
    void mousePressEvent(QMouseEvent *e) override;
    void mouseMoveEvent(QMouseEvent *e) override;
    void mouseDoubleClickEvent(QMouseEvent *e) override;
    void wheelEvent(QWheelEvent *e) override;
    void contextMenuEvent(QContextMenuEvent *e) override;

private:
    // rotated coordinates: sx, sy on the screen (before scaling), depth (larger = nearer)
    struct R3 {
        double sx, sy, depth;
    };
    R3 rotate(const P3 &p) const;

    QVector<P3> m_grid;
    int m_nu = 0, m_nv = 0;
    bool m_filled = true, m_wrap = false, m_colorbar = true;
    QVector<Line> m_lines;
    QVector<Marker> m_markers;
    double m_xmin = -1, m_xmax = 1, m_ymin = -1, m_ymax = 1, m_zmin = 0, m_zmax = 1;
    QString m_zlabel, m_title;
    double m_azim = -65, m_elev = 30, m_zoom = 1;  // degrees, pyfda's default view
    QPoint m_last;
};

class ThreeDView : public DesignView {
    Q_OBJECT
public:
    explicit ThreeDView(QWidget *parent = nullptr);

protected:
    void redraw() override;

private:
    Plot3D *m_plot;
    QComboBox *m_mode;
    QCheckBox *m_log, *m_polar, *m_hf, *m_pz, *m_uc, *m_cbar;
    QLineEdit *m_bottom, *m_top;
};
