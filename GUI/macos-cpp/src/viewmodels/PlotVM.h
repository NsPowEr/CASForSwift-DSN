#pragma once
#include <QObject>
#include <QVariantList>
#include <QPointF>
#include <QString>
#include <QColor>
#include <vector>
#include <cmath>

namespace cas::gui {

class PlotVM : public QObject {
    Q_OBJECT
    Q_PROPERTY(QVariantList seriesList READ seriesList NOTIFY seriesChanged)
    Q_PROPERTY(bool hasPoints READ hasPoints NOTIFY seriesChanged)
    Q_PROPERTY(QString state READ state NOTIFY stateChanged)
    Q_PROPERTY(QString status READ status NOTIFY statusChanged)
    Q_PROPERTY(int sampleCount READ sampleCount NOTIFY seriesChanged)
    Q_PROPERTY(QString variable READ variable WRITE setVariable NOTIFY settingsChanged)
    Q_PROPERTY(double xMin READ xMin WRITE setXMin NOTIFY rangeChanged)
    Q_PROPERTY(double xMax READ xMax WRITE setXMax NOTIFY rangeChanged)
    Q_PROPERTY(double yMin READ yMin NOTIFY rangeChanged)
    Q_PROPERTY(double yMax READ yMax NOTIFY rangeChanged)

public:
    struct Series {
        QString name;
        QVariantList points;
        QColor color;
    };

    explicit PlotVM(QObject* parent = nullptr) : QObject(parent) {}

    QVariantList seriesList() const {
        QVariantList result;
        for (const auto& s : m_series) {
            QVariantMap m;
            m["name"] = s.name;
            m["points"] = s.points;
            m["color"] = s.color;
            result.append(m);
        }
        return result;
    }

    bool hasPoints() const {
        for (const auto& s : m_series) if (!s.points.isEmpty()) return true;
        return false;
    }

    QVariantList points() const {
        if (m_series.empty()) return {};
        return m_series.front().points;
    }

    QString state() const { return m_state; }
    QString status() const { return m_status; }
    int sampleCount() const { 
        int total = 0;
        for (const auto& s : m_series) total += s.points.size();
        return total;
    }
    QString variable() const { return m_variable; }
    double xMin() const { return m_xMin; }
    double xMax() const { return m_xMax; }
    double yMin() const { return m_yMin; }
    double yMax() const { return m_yMax; }

    void setVariable(const QString& value) {
        if (m_variable == value) return;
        emit settingsAboutToChange();
        m_variable = value;
        emit settingsChanged();
    }

    Q_INVOKABLE void setXMin(double v) {
        if (v == m_xMin) return;
        setRange(v, m_xMax);
    }
    Q_INVOKABLE void setXMax(double v) {
        if (v == m_xMax) return;
        setRange(m_xMin, v);
    }
    Q_INVOKABLE void setRange(double minValue, double maxValue) {
        if (!std::isfinite(minValue) || !std::isfinite(maxValue)) return;
        if (minValue == m_xMin && maxValue == m_xMax) return;
        emit settingsAboutToChange();
        m_xMin = minValue;
        m_xMax = maxValue;
        emit rangeChanged();
        emit settingsChanged();
    }

    Q_INVOKABLE void zoomIn() { scaleRange(0.5); }
    Q_INVOKABLE void zoomOut() { scaleRange(2.0); }
    Q_INVOKABLE void panLeft() { shiftRange(-0.2); }
    Q_INVOKABLE void panRight() { shiftRange(0.2); }
    Q_INVOKABLE void resetRange() {
        setRange(m_defaultXMin, m_defaultXMax);
    }
    Q_INVOKABLE void applyPreset(const QString& presetId) {
        if (presetId == QStringLiteral("symmetric_10")) setRange(-10.0, 10.0);
        else if (presetId == QStringLiteral("positive_10")) setRange(0.0, 10.0);
        else if (presetId == QStringLiteral("trig_pi")) setRange(-M_PI, M_PI);
    }

    Q_INVOKABLE void addSeries(const QString& name, const QVariantList& points, const QColor& color = QColor()) {
        clearSeries();
        if (points.isEmpty()) {
            setEmpty(QStringLiteral("Plot empty: no samples"));
            return;
        }
        m_series.push_back({name, points, color.isValid() ? color : nextColor()});
        emit seriesChanged();
        updateYRange();
        setState(QStringLiteral("ready"));
        setStatus(QStringLiteral("Plot ready: %1 samples").arg(sampleCount()));
    }

    Q_INVOKABLE void clearSeries() {
        m_series.clear();
        m_colorIndex = 0;
        emit seriesChanged();
    }

    Q_INVOKABLE void clear() {
        m_variable = "x";
        setRange(m_defaultXMin, m_defaultXMax);
        clearSeries();
        m_yMin = -10.0;
        m_yMax = 10.0;
        setState(QStringLiteral("idle"));
        setStatus(QStringLiteral("Plot cleared"));
        emit rangeChanged();
    }

    void beginLoading() {
        clearSeries();
        m_yMin = -10.0;
        m_yMax = 10.0;
        setState(QStringLiteral("loading"));
        setStatus(QStringLiteral("Sampling plot..."));
        emit rangeChanged();
    }
    void setIdle(const QString& status = QStringLiteral("Plot idle")) {
        clearSeries();
        m_yMin = -10.0;
        m_yMax = 10.0;
        setState(QStringLiteral("idle"));
        setStatus(status);
        emit rangeChanged();
    }
    void setEmpty(const QString& status = QStringLiteral("No plot data")) {
        clearSeries();
        m_yMin = -10.0;
        m_yMax = 10.0;
        setState(QStringLiteral("empty"));
        setStatus(status);
        emit rangeChanged();
    }

    void setError(const QString& status) {
        clearSeries();
        setState(QStringLiteral("error"));
        setStatus(status);
    }

signals:
    void settingsAboutToChange();
    void seriesChanged();
    void rangeChanged();
    void stateChanged();
    void statusChanged();
    void settingsChanged();

private:
    void setState(const QString& value) {
        if (m_state == value) return;
        m_state = value;
        emit stateChanged();
    }
    void setStatus(const QString& value) {
        if (m_status == value) return;
        m_status = value;
        emit statusChanged();
    }
    void scaleRange(double factor) {
        if (!(factor > 0.0)) return;
        const double center = (m_xMin + m_xMax) / 2.0;
        const double halfSpan = ((m_xMax - m_xMin) / 2.0) * factor;
        setRange(center - halfSpan, center + halfSpan);
    }
    void shiftRange(double fraction) {
        const double span = m_xMax - m_xMin;
        const double delta = span * fraction;
        setRange(m_xMin + delta, m_xMax + delta);
    }
    QColor nextColor() {
        static const QVector<QColor> palette = {
            "#007aff", "#34c759", "#ff9500", "#ff3b30", "#af52de", "#5856d6"
        };
        return palette[m_colorIndex++ % palette.size()];
    }
    void updateYRange() {
        bool first = true;
        for (const auto& s : m_series) {
            for (const auto& v : s.points) {
                QPointF p = v.toPointF();
                if (first) { m_yMin = m_yMax = p.y(); first = false; }
                else { m_yMin = std::min(m_yMin, p.y()); m_yMax = std::max(m_yMax, p.y()); }
            }
        }
        if (m_yMin == m_yMax) { m_yMin -= 1.0; m_yMax += 1.0; }
        emit rangeChanged();
    }

    std::vector<Series> m_series;
    QString m_state = "idle";
    QString m_status = "Plot idle";
    QString m_variable = "x";
    double m_xMin = -10;
    double m_xMax = 10;
    double m_yMin = -10;
    double m_yMax = 10;
    double m_defaultXMin = -10;
    double m_defaultXMax = 10;
    int m_colorIndex = 0;
};

} // namespace cas::gui
