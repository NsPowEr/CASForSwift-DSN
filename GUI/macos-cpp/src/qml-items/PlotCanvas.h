#pragma once
#include <QQuickPaintedItem>
#include <QPainter>
#include <QVariantList>

namespace cas::gui {

class PlotCanvas : public QQuickPaintedItem {
    Q_OBJECT
    QML_ELEMENT
    Q_PROPERTY(QVariantList seriesList READ seriesList WRITE setSeriesList NOTIFY seriesChanged)
    Q_PROPERTY(QColor lineColor READ lineColor WRITE setLineColor NOTIFY lineColorChanged)
    Q_PROPERTY(double xMin READ xMin WRITE setXMin NOTIFY rangesChanged)
    Q_PROPERTY(double xMax READ xMax WRITE setXMax NOTIFY rangesChanged)
    Q_PROPERTY(double yMin READ yMin WRITE setYMin NOTIFY rangesChanged)
    Q_PROPERTY(double yMax READ yMax WRITE setYMax NOTIFY rangesChanged)

public:
    explicit PlotCanvas(QQuickItem* parent = nullptr);
    void paint(QPainter* p) override;

    QVariantList seriesList() const { return m_seriesList; }
    void setSeriesList(const QVariantList& list);

    QColor lineColor() const { return m_lineColor; }
    void setLineColor(const QColor& c);
    double xMin() const { return m_xMin; }
    double xMax() const { return m_xMax; }
    double yMin() const { return m_yMin; }
    double yMax() const { return m_yMax; }
    void setXMin(double value);
    void setXMax(double value);
    void setYMin(double value);
    void setYMax(double value);

signals:
    void seriesChanged();
    void lineColorChanged();
    void rangesChanged();

private:
    QVariantList m_seriesList;
    QColor m_lineColor = Qt::blue;
    double m_xMin = -10.0;
    double m_xMax = 10.0;
    double m_yMin = -10.0;
    double m_yMax = 10.0;
};

} // namespace cas::gui
