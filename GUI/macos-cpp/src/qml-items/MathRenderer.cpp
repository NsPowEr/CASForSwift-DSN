#include "MathRenderer.h"
#include <QFontDatabase>
#include <QFont>
#include <QFontMetricsF>
#include <QPainter>
#include <cmath>
#include <vector>

namespace {

QString math_font_family() {
    const auto families = QFontDatabase::families();
    if (families.contains("STIX Two Math")) {
        return "STIX Two Math";
    }
    if (families.contains("STIX Two Text")) {
        return "STIX Two Text";
    }
    if (families.contains("Times New Roman")) {
        return "Times New Roman";
    }
    return "Helvetica Neue";
}

QString strip_outer_braces(QString text) {
    if (text.size() < 2 || !text.startsWith('{') || !text.endsWith('}')) {
        return text;
    }

    int depth = 0;
    for (int i = 0; i < text.size(); ++i) {
        if (text[i] == '{') {
            ++depth;
        } else if (text[i] == '}') {
            --depth;
            if (depth == 0 && i != text.size() - 1) {
                return text;
            }
        }
    }
    return depth == 0 ? text.mid(1, text.size() - 2) : text;
}

int find_matching_brace(const QString& text, int open_index) {
    if (open_index < 0 || open_index >= text.size() || text[open_index] != '{') {
        return -1;
    }
    int depth = 0;
    for (int i = open_index; i < text.size(); ++i) {
        if (text[i] == '{') {
            ++depth;
        } else if (text[i] == '}') {
            --depth;
            if (depth == 0) {
                return i;
            }
        }
    }
    return -1;
}

int find_top_level_script_op(const QString& text) {
    int depth = 0;
    for (int i = 0; i < text.size(); ++i) {
        if (text[i] == '{') {
            ++depth;
        } else if (text[i] == '}') {
            --depth;
        } else if (depth == 0 && (text[i] == '^' || text[i] == '_')) {
            return i;
        }
    }
    return -1;
}

bool extract_script_token(const QString& text, QString* script, QString* tail) {
    const QString trimmed = text.trimmed();
    if (trimmed.isEmpty()) {
        script->clear();
        tail->clear();
        return false;
    }

    if (trimmed[0] == '{') {
        const int close = find_matching_brace(trimmed, 0);
        if (close == -1) {
            *script = trimmed.mid(1);
            tail->clear();
            return true;
        }
        *script = trimmed.mid(1, close - 1);
        *tail = trimmed.mid(close + 1).trimmed();
        return true;
    }

    *script = QString(trimmed[0]);
    *tail = trimmed.mid(1).trimmed();
    return true;
}

QString normalize_text(QString text) {
    return text
        .replace("\\cdot", "·")
        .replace("\\pi", "π")
        .replace("\\theta", "θ")
        .replace("\\alpha", "α")
        .replace("\\beta", "β")
        .replace("\\gamma", "γ")
        .replace("\\lambda", "λ")
        .replace("\\phi", "φ")
        .replace("\\omega", "ω")
        .replace("\\sigma", "σ")
        .replace("\\tau", "τ")
        .replace("\\delta", "δ")
        .replace("\\epsilon", "ε")
        .replace("\\eta", "η")
        .replace("\\rho", "ρ")
        .replace("\\int", "∫")
        .replace("\\sum", "Σ")
        .replace("\\infty", "∞")
        .replace("\\to", "→")
        .replace("\\rightarrow", "→")
        .replace("\\partial", "∂")
        .replace("\\nabla", "∇")
        .replace("\\pm", "±")
        .replace("\\neq", "≠")
        .replace("\\le", "≤")
        .replace("\\ge", "≥")
        .replace("\\approx", "≈")
        .replace("\\equiv", "≡")
        .replace("\\sin", "sin")
        .replace("\\cos", "cos")
        .replace("\\tan", "tan")
        .replace("\\exp", "exp")
        .replace("\\log", "log")
        .replace("\\ln", "ln")
        .replace("\\sqrt", "√");
}

}

struct MathRenderer::LayoutBox {
    enum Kind { Text, Frac, Sup, Sub, Group, Sqrt };
    Kind kind = Text;
    QString content;
    qreal x = 0, y = 0, w = 0, h = 0, ascent = 0;
    std::vector<std::unique_ptr<LayoutBox>> children;

    void draw(QPainter* p, qreal ox, qreal oy, qreal fontSize) {
        if (kind == Text) {
            p->drawText(QPointF(ox + x, oy + y + ascent), content);
        }
        for (auto& child : children) {
            qreal childScale = 1.0;
            if (kind == Sup || kind == Sub) childScale = 0.7;
            
            QFont f = p->font();
            f.setPointSizeF(fontSize * childScale);
            p->setFont(f);
            child->draw(p, ox + x, oy + y, fontSize * childScale);
            f.setPointSizeF(fontSize);
            p->setFont(f);
        }
        if (kind == Frac) {
            p->setPen(p->pen().color());
            // Draw fraction line between numerator and denominator
            qreal lineY = oy + y + children[0]->h + 1;
            p->drawLine(QPointF(ox + x, lineY), QPointF(ox + x + w, lineY));
        } else if (kind == Sqrt && !children.empty()) {
            const qreal baseX = ox + x;
            const qreal baseY = oy + y;
            const qreal rootWidth = std::max<qreal>(8.0, fontSize * 0.55);
            const qreal tickX = baseX + rootWidth * 0.2;
            const qreal tickY = baseY + children[0]->h * 0.75;
            const qreal midX = baseX + rootWidth * 0.45;
            const qreal midY = baseY + children[0]->h;
            const qreal armX = baseX + rootWidth * 0.72;
            const qreal armTopY = baseY + children[0]->y + 2;
            p->drawLine(QPointF(tickX, tickY), QPointF(midX, midY));
            p->drawLine(QPointF(midX, midY), QPointF(armX, armTopY));
            p->drawLine(QPointF(armX, armTopY), QPointF(baseX + w, armTopY));
        }
    }
};

MathRenderer::MathRenderer(QQuickItem* parent) : QQuickPaintedItem(parent) {
    setAntialiasing(true);
}

MathRenderer::~MathRenderer() = default;

void MathRenderer::setLatex(const QString& l) {
    if (m_latex == l) return;
    m_latex = l;
    rebuildLayout();
    emit changed();
}

void MathRenderer::setFontSize(qreal v) {
    if (m_fontSize == v) return;
    m_fontSize = v;
    rebuildLayout();
    emit changed();
}

void MathRenderer::setColor(QColor c) {
    if (m_color == c) return;
    m_color = c;
    emit changed();
    update();
}

void MathRenderer::setBlock(bool b) {
    if (m_block == b) return;
    m_block = b;
    emit changed();
    update();
}

void MathRenderer::rebuildLayout() {
    QFont f(math_font_family(), m_fontSize);
    QFontMetricsF fm(f);

    auto parse = [&](auto self, const QString& src, bool small) -> std::unique_ptr<LayoutBox> {
        auto box = std::make_unique<LayoutBox>();
        const QString trimmed = strip_outer_braces(src.trimmed());
        qreal currentFontSize = small ? m_fontSize * 0.7 : m_fontSize;
        QFont currentFont(math_font_family(), currentFontSize);
        QFontMetricsF cfm(currentFont);

        if (trimmed.isEmpty()) {
            box->kind = LayoutBox::Text;
            box->content.clear();
            box->w = 0;
            box->h = cfm.height();
            box->ascent = cfm.ascent();
        } else if (trimmed.startsWith("\\frac{")) {
            box->kind = LayoutBox::Frac;
            const int num_open = trimmed.indexOf('{', 5);
            const int num_close = find_matching_brace(trimmed, num_open);
            const int den_open = num_close + 1 < trimmed.size() ? trimmed.indexOf('{', num_close + 1) : -1;
            const int den_close = find_matching_brace(trimmed, den_open);

            if (num_open != -1 && num_close != -1 && den_open == num_close + 1 && den_close != -1) {
                const QString num = trimmed.mid(num_open + 1, num_close - num_open - 1);
                const QString den = trimmed.mid(den_open + 1, den_close - den_open - 1);
                const QString tail = trimmed.mid(den_close + 1).trimmed();

                auto numerator = self(self, num, small);
                auto denominator = self(self, den, small);
                numerator->x = 0;
                numerator->y = 0;
                denominator->x = 0;
                denominator->y = numerator->h + 6;

                const qreal frac_width = std::max(numerator->w, denominator->w) + 8;
                numerator->x = (frac_width - numerator->w) / 2;
                denominator->x = (frac_width - denominator->w) / 2;

                box->children.push_back(std::move(numerator));
                box->children.push_back(std::move(denominator));
                box->w = frac_width;
                box->h = box->children[0]->h + box->children[1]->h + 6;
                box->ascent = box->children[0]->h;

                if (!tail.isEmpty()) {
                    auto trailing = self(self, tail, small);
                    trailing->x = box->w + 6;
                    trailing->y = std::max<qreal>(0, box->ascent - trailing->ascent);
                    box->w += 6 + trailing->w;
                    box->h = std::max(box->h, trailing->h);
                    box->children.push_back(std::move(trailing));
                }
            }
        } else if (trimmed.startsWith("\\sqrt{")) {
            box->kind = LayoutBox::Sqrt;
            const int radicand_open = trimmed.indexOf('{', 5);
            const int radicand_close = find_matching_brace(trimmed, radicand_open);

            if (radicand_open != -1 && radicand_close != -1) {
                const QString radicand = trimmed.mid(radicand_open + 1, radicand_close - radicand_open - 1);
                const QString tail = trimmed.mid(radicand_close + 1).trimmed();
                auto radicandBox = self(self, radicand, small);
                const qreal rootWidth = std::max<qreal>(8.0, currentFontSize * 0.55);
                const qreal topPadding = 4.0;
                radicandBox->x = rootWidth + 4;
                radicandBox->y = topPadding;

                box->children.push_back(std::move(radicandBox));
                box->w = box->children[0]->x + box->children[0]->w + 2;
                box->h = box->children[0]->h + topPadding + 2;
                box->ascent = box->children[0]->y + box->children[0]->ascent;

                if (!tail.isEmpty()) {
                    auto trailing = self(self, tail, small);
                    trailing->x = box->w + 4;
                    trailing->y = std::max<qreal>(0, box->ascent - trailing->ascent);
                    box->w += 4 + trailing->w;
                    box->h = std::max(box->h, trailing->h);
                    box->children.push_back(std::move(trailing));
                }
            }
        } else if (const int opIdx = find_top_level_script_op(trimmed); opIdx != -1) {
            box->kind = LayoutBox::Group;
            const QString base = trimmed.left(opIdx).trimmed();
            const QChar firstOp = trimmed[opIdx];
            QString firstScript;
            QString tail;
            extract_script_token(trimmed.mid(opIdx + 1), &firstScript, &tail);

            auto baseBox = self(self, base, small);
            auto firstScriptBox = self(self, firstScript, true);
            std::unique_ptr<LayoutBox> secondScriptBox;
            QChar secondOp;
            QString trailingTail = tail;

            if (!tail.isEmpty() && (tail[0] == '^' || tail[0] == '_')) {
                secondOp = tail[0];
                QString secondScript;
                extract_script_token(tail.mid(1), &secondScript, &trailingTail);
                secondScriptBox = self(self, secondScript, true);
            }

            baseBox->x = 0;
            baseBox->y = 0;
            firstScriptBox->x = baseBox->w + 2;
            firstScriptBox->y = firstOp == '^' ? -baseBox->ascent * 0.45 : baseBox->ascent * 0.35;

            qreal scriptColumnWidth = firstScriptBox->w;
            qreal top = std::min<qreal>(0.0, firstScriptBox->y);
            qreal bottom = std::max(baseBox->h, firstScriptBox->y + firstScriptBox->h);

            if (secondScriptBox) {
                secondScriptBox->x = baseBox->w + 2;
                secondScriptBox->y = secondOp == '^' ? -baseBox->ascent * 0.45 : baseBox->ascent * 0.35;
                scriptColumnWidth = std::max(scriptColumnWidth, secondScriptBox->w);
                top = std::min(top, secondScriptBox->y);
                bottom = std::max(bottom, secondScriptBox->y + secondScriptBox->h);
            }

            if (top < 0) {
                baseBox->y -= top;
                firstScriptBox->y -= top;
                if (secondScriptBox) {
                    secondScriptBox->y -= top;
                }
            }

            box->w = baseBox->w + 2 + scriptColumnWidth;
            box->h = bottom - top;
            box->ascent = baseBox->ascent - top;
            box->children.push_back(std::move(baseBox));
            box->children.push_back(std::move(firstScriptBox));
            if (secondScriptBox) {
                box->children.push_back(std::move(secondScriptBox));
            }

            if (!trailingTail.isEmpty()) {
                auto trailing = self(self, trailingTail, small);
                trailing->x = box->w + 4;
                trailing->y = std::max<qreal>(0, box->ascent - trailing->ascent);
                box->w += 4 + trailing->w;
                box->h = std::max(box->h, trailing->h);
                box->children.push_back(std::move(trailing));
            }
        } else {
            box->kind = LayoutBox::Text;
            box->content = normalize_text(trimmed);
            box->w = cfm.horizontalAdvance(box->content);
            box->h = cfm.height();
            box->ascent = cfm.ascent();
        }
        return box;
    };

    m_root = parse(parse, m_latex, false);
    if (m_root) {
        m_contentSize = QSizeF(m_root->w, m_root->h);
        setImplicitWidth(m_root->w);
        setImplicitHeight(m_root->h);
    }
    emit layoutChanged();
    update();
}

void MathRenderer::paint(QPainter* p) {
    if (!m_root) return;
    p->setPen(m_color);
    p->setFont(QFont(math_font_family(), m_fontSize));
    m_root->draw(p, 0, 0, m_fontSize);
}
