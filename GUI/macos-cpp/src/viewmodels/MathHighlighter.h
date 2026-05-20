#pragma once

#include <QSyntaxHighlighter>
#include <QTextCharFormat>
#include <QRegularExpression>
#include <QtQml>

namespace cas::gui {

class MathHighlighter : public QSyntaxHighlighter {
    Q_OBJECT
    QML_ELEMENT
    Q_PROPERTY(QObject* target READ target WRITE setTarget NOTIFY targetChanged)

public:
    explicit MathHighlighter(QObject* parent = nullptr);

    QObject* target() const { return m_target; }
    void setTarget(QObject* target);

signals:
    void targetChanged();

protected:
    void highlightBlock(const QString& text) override;

private:
    struct HighlightingRule {
        QRegularExpression pattern;
        QTextCharFormat format;
    };
    QVector<HighlightingRule> m_highlightingRules;

    QTextCharFormat m_functionFormat;
    QTextCharFormat m_variableFormat;
    QTextCharFormat m_numberFormat;
    QTextCharFormat m_operatorFormat;
    QTextCharFormat m_bracketFormat;
    QTextCharFormat m_latexCommandFormat;

    QPointer<QObject> m_target;
    void setupRules();
};

} // namespace cas::gui
