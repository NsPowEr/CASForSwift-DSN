// ThemeVM.h — espone i design tokens (estratti dal prototipo HTML) a QML.
// Tutti i colori sono derivati da `accent` + `dark` come nella funzione palette() del prototipo.

#pragma once
#include <QObject>
#include <QColor>
#include <QQmlEngine>

class ThemeVM : public QObject {
    Q_OBJECT

    Q_PROPERTY(bool    dark    READ dark    WRITE setDark    NOTIFY changed)
    Q_PROPERTY(QColor  accent  READ accent  WRITE setAccent  NOTIFY changed)

    // Palette derivata
    Q_PROPERTY(QColor bg         READ bg         NOTIFY changed)
    Q_PROPERTY(QColor bgRaised   READ bgRaised   NOTIFY changed)
    Q_PROPERTY(QColor bgCard     READ bgCard     NOTIFY changed)
    Q_PROPERTY(QColor stroke     READ stroke     NOTIFY changed)
    Q_PROPERTY(QColor text       READ text       NOTIFY changed)
    Q_PROPERTY(QColor textMuted  READ textMuted  NOTIFY changed)
    Q_PROPERTY(QColor textFaint  READ textFaint  NOTIFY changed)
    Q_PROPERTY(QColor accentSoft READ accentSoft NOTIFY changed)
    Q_PROPERTY(QColor ok         READ ok         CONSTANT)
    Q_PROPERTY(QColor error      READ error      CONSTANT)

    // Type
    Q_PROPERTY(QString fontUI      READ fontUI      CONSTANT)
    Q_PROPERTY(QString fontDisplay READ fontDisplay CONSTANT)
    Q_PROPERTY(QString fontMono    READ fontMono    CONSTANT)
    Q_PROPERTY(QString fontMath    READ fontMath    CONSTANT)

    // Spacing tokens
    Q_PROPERTY(int radiusS  READ radiusS  CONSTANT)
    Q_PROPERTY(int radiusM  READ radiusM  CONSTANT)
    Q_PROPERTY(int radiusL  READ radiusL  CONSTANT)

public:
    explicit ThemeVM(QObject* parent = nullptr);

    bool   dark()   const { return m_dark; }
    QColor accent() const { return m_accent; }
    void setDark(bool v);
    void setAccent(QColor c);

    QColor bg() const;
    QColor bgRaised() const;
    QColor bgCard() const;
    QColor stroke() const;
    QColor text() const;
    QColor textMuted() const;
    QColor textFaint() const;
    QColor accentSoft() const;
    QColor ok()    const { return QColor("#30D158"); }
    QColor error() const { return QColor("#FF453A"); }

    QString fontUI()      const { return "Helvetica Neue"; }
    QString fontDisplay() const { return "Helvetica Neue"; }
    QString fontMono()    const { return "Menlo"; }
    QString fontMath()    const { return "STIX Two Math"; }

    int radiusS() const { return 6; }
    int radiusM() const { return 10; }
    int radiusL() const { return 16; }

    Q_INVOKABLE void resetAccent() { setAccent(QColor("#C77DFF")); }
    Q_INVOKABLE void persist();   // salva su QSettings
    Q_INVOKABLE void load();

signals:
    void changed();

private:
    bool   m_dark   = true;
    QColor m_accent = QColor("#C77DFF");
};
