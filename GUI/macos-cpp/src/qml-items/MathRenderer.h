// MathRenderer.h — Custom Qt Quick item che disegna espressioni matematiche LaTeX-like.
// Opzione A (raccomandata): integra microtex (https://github.com/NanoMichael/MicroTeX) — TeX rendering puro C++/OpenGL.
// Opzione B: porting del parser "Tex-lite" del prototipo HTML (sufficiente per il 90% dei casi).
//
// L'item è guidato da Q_PROPERTY e ridisegna su `update()`.

#pragma once
#include <QQuickPaintedItem>
#include <QPainter>
#include <memory>

class MathRenderer : public QQuickPaintedItem {
    Q_OBJECT
    QML_ELEMENT

    Q_PROPERTY(QString latex    READ latex    WRITE setLatex    NOTIFY changed)
    Q_PROPERTY(qreal   fontSize READ fontSize WRITE setFontSize NOTIFY changed)
    Q_PROPERTY(QColor  color    READ color    WRITE setColor    NOTIFY changed)
    Q_PROPERTY(bool    block    READ block    WRITE setBlock    NOTIFY changed)
    // dimensione naturale calcolata dopo il layout — utile per content-sizing in QML
    Q_PROPERTY(qreal   contentWidth  READ contentWidth  NOTIFY layoutChanged)
    Q_PROPERTY(qreal   contentHeight READ contentHeight NOTIFY layoutChanged)

public:
    explicit MathRenderer(QQuickItem* parent = nullptr);
    ~MathRenderer() override;

    void paint(QPainter* p) override;

    QString latex()    const { return m_latex; }
    qreal   fontSize() const { return m_fontSize; }
    QColor  color()    const { return m_color; }
    bool    block()    const { return m_block; }
    qreal   contentWidth()  const { return m_contentSize.width();  }
    qreal   contentHeight() const { return m_contentSize.height(); }

    void setLatex(const QString& s);
    void setFontSize(qreal v);
    void setColor(QColor c);
    void setBlock(bool b);

signals:
    void changed();
    void layoutChanged();

private:
    void rebuildLayout();

    QString  m_latex;
    qreal    m_fontSize = 18.0;
    QColor   m_color    = Qt::black;
    bool     m_block    = false;
    QSizeF   m_contentSize;

    // PIMPL: tiene il box-model della formula (fraction, sup, sub, sqrt, sum, int, matrix, ...)
    struct LayoutBox;
    std::unique_ptr<LayoutBox> m_root;
};
