#include "MathHighlighter.h"
#include <QQuickTextDocument>
#include <QTextDocument>

namespace cas::gui {

MathHighlighter::MathHighlighter(QObject* parent) : QSyntaxHighlighter(parent) {
    setupRules();
}

void MathHighlighter::setTarget(QObject* target) {
    if (m_target == target) return;
    m_target = target;
    
    if (m_target) {
        // QML TextArea has a 'textDocument' property that is a QQuickTextDocument wrapper
        QVariant docVar = m_target->property("textDocument");
        if (docVar.isValid()) {
            QQuickTextDocument* quickDoc = qvariant_cast<QQuickTextDocument*>(docVar);
            if (quickDoc) {
                setDocument(quickDoc->textDocument());
            }
        }
    } else {
        setDocument(nullptr);
    }
    
    emit targetChanged();
}

void MathHighlighter::setupRules() {
    // Basic CAS syntax rules
    
    // LaTeX Commands (\sin, \frac, etc) - stop at non-letters
    m_latexCommandFormat.setForeground(QColor("#7b1fa2")); // Purple
    m_latexCommandFormat.setFontWeight(QFont::Bold);
    m_highlightingRules.append({QRegularExpression(R"(\\[a-zA-Z]+)"), m_latexCommandFormat});

    // Valid Math functions (recognized by CAS) - use word boundaries
    m_functionFormat.setForeground(QColor("#007aff")); // macOS Blue
    m_functionFormat.setFontWeight(QFont::Medium);
    m_highlightingRules.append({QRegularExpression(R"(\b(sin|cos|tan|arcsin|arccos|arctan|exp|log|ln|diff|integrate|solve|limit|simplify|expand|factor|assume|forget|plot|roots|gcd|lcm|derivative|integral|sum|product|min|max|abs|sqrt)\b)"), m_functionFormat});

    // Numbers - more robust pattern
    m_numberFormat.setForeground(QColor("#d35400")); // Orange
    m_highlightingRules.append({QRegularExpression(R"(\b\d+(\.\d+)?([eE][+-]?\d+)?\b)"), m_numberFormat});

    // Variables - ensure they don't consume parts of functions or commands
    // and don't match if preceded by backslash (LaTeX command)
    m_variableFormat.setForeground(QColor("#2c3e50")); // Dark grey/blue
    m_variableFormat.setFontItalic(true);
    m_highlightingRules.append({QRegularExpression(R"((?<!\\[a-zA-Z]*)\b[a-zA-Z][a-zA-Z0-9]*\b)"), m_variableFormat});

    // Operators
    m_operatorFormat.setForeground(QColor("#c0392b")); // Red
    m_highlightingRules.append({QRegularExpression(R"([\+\-\*\/\^=])"), m_operatorFormat});

    // Brackets
    m_bracketFormat.setForeground(QColor("#7f8c8d")); // Grey
    m_highlightingRules.append({QRegularExpression(R"([\{\}\(\)\[\]])"), m_bracketFormat});
}

void MathHighlighter::highlightBlock(const QString& text) {
    // 1. First, mark everything as default
    setFormat(0, text.length(), QTextCharFormat());

    // 2. Apply rules in order of specificity
    // We use a bitmask or similar if we wanted zero overlap, 
    // but for CAS standard rules, we can just apply from most specific to least.
    
    // Commands and Functions first (must not be sub-matched by variable rule)
    for (const HighlightingRule& rule : m_highlightingRules) {
        QRegularExpressionMatchIterator matchIterator = rule.pattern.globalMatch(text);
        while (matchIterator.hasNext()) {
            QRegularExpressionMatch match = matchIterator.next();
            setFormat(match.capturedStart(), match.capturedLength(), rule.format);
        }
    }
}

} // namespace cas::gui
