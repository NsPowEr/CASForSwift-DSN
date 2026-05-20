#include "PlotCanvas.h"
#include <QPainterPath>
#include <QVariantMap>

namespace cas::gui {

PlotCanvas::PlotCanvas(QQuickItem* parent) : QQuickPaintedItem(parent) {
    setAntialiasing(true);
}

void PlotCanvas::setSeriesList(const QVariantList& list) {
    m_seriesList = list;
    emit seriesChanged();
    update();
}

void PlotCanvas::setLineColor(const QColor& c) {
    m_lineColor = c;
    emit lineColorChanged();
    update();
}

void PlotCanvas::setXMin(double value) {
    if (m_xMin == value) return;
    m_xMin = value;
    emit rangesChanged();
    update();
}

void PlotCanvas::setXMax(double value) {
    if (m_xMax == value) return;
    m_xMax = value;
    emit rangesChanged();
    update();
}

void PlotCanvas::setYMin(double value) {
    if (m_yMin == value) return;
    m_yMin = value;
    emit rangesChanged();
    update();
}

void PlotCanvas::setYMax(double value) {
    if (m_yMax == value) return;
    m_yMax = value;
    emit rangesChanged();
    update();
}

void PlotCanvas::paint(QPainter* p) {
    p->setRenderHint(QPainter::Antialiasing);

    qreal w = width();
    qreal h = height();
    const double xSpan = (m_xMax - m_xMin) == 0.0 ? 1.0 : (m_xMax - m_xMin);
    const double ySpan = (m_yMax - m_yMin) == 0.0 ? 1.0 : (m_yMax - m_yMin);

    // Draw grid
    p->setPen(QPen(QColor(220, 220, 220, 100), 1, Qt::DashLine));
    const int gridSteps = 10;
    for (int i = 0; i <= gridSteps; ++i) {
        qreal px = w * i / gridSteps;
        p->drawLine(QPointF(px, 0), QPointF(px, h));

        qreal py = h * i / gridSteps;
        p->drawLine(QPointF(0, py), QPointF(w, py));
    }

    // Draw axes
    p->setPen(QPen(QColor(128, 128, 128, 150), 1, Qt::SolidLine));
    
    // X-axis (y = 0)
    if (m_yMin <= 0 && m_yMax >= 0) {
        qreal pyZero = h - static_cast<qreal>((0 - m_yMin) / ySpan) * h;
        p->drawLine(QPointF(0, pyZero), QPointF(w, pyZero));
    }
    
    // Y-axis (x = 0)
    if (m_xMin <= 0 && m_xMax >= 0) {
        qreal pxZero = static_cast<qreal>((0 - m_xMin) / xSpan) * w;
        p->drawLine(QPointF(pxZero, 0), QPointF(pxZero, h));
    }

    if (m_seriesList.isEmpty()) return;

    for (const auto& sVar : m_seriesList) {
        QVariantMap s = sVar.toMap();
        QVariantList points = s["points"].toList();
        if (points.isEmpty()) continue;

        QColor color = s["color"].value<QColor>();
        p->setPen(QPen(color, 2, Qt::SolidLine));

        QPainterPath path;
        bool first = true;
        
        for (const auto& v : points) {
            QPointF pt = v.toPointF();
            qreal px = static_cast<qreal>((pt.x() - m_xMin) / xSpan) * w;
            qreal py = h - static_cast<qreal>((pt.y() - m_yMin) / ySpan) * h;
            
            if (first) {
                path.moveTo(px, py);
                first = false;
            } else {
                path.lineTo(px, py);
            }
        }
        p->drawPath(path);
    }
}

} // namespace cas::gui
