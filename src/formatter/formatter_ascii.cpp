#include "cas/formatter.hpp"
#include <sstream>
#include <algorithm>

namespace cas::formatter {

namespace {
std::string spaces(int n) { return n > 0 ? std::string(n, ' ') : ""; }
std::string dashes(int n) { return n > 0 ? std::string(n, '-') : ""; }

/**
 * Unisce due box orizzontalmente allineandoli sulla baseline.
 */
LayoutBox join_horizontal(const LayoutBox& left, const LayoutBox& right, const std::string& op = "") {
    LayoutBox res;
    int op_w = static_cast<int>(op.length());
    res.width = left.width + right.width + op_w;
    
    // Altezza sopra la baseline e sotto
    int above_left = left.baseline;
    int below_left = left.height - left.baseline - 1;
    int above_right = right.baseline;
    int below_right = right.height - right.baseline - 1;
    
    int max_above = std::max(above_left, above_right);
    int max_below = std::max(below_left, below_right);
    
    res.height = max_above + max_below + 1;
    res.baseline = max_above;
    
    for (int i = 0; i < res.height; ++i) {
        std::string line;
        
        // Parte sinistra
        int l_idx = i - (max_above - above_left);
        if (l_idx >= 0 && l_idx < left.height) line += left.lines[l_idx];
        else line += spaces(left.width);
        
        // Operatore
        if (i == res.baseline) line += op;
        else line += spaces(op_w);
        
        // Parte destra
        int r_idx = i - (max_above - above_right);
        if (r_idx >= 0 && r_idx < right.height) line += right.lines[r_idx];
        else line += spaces(right.width);
        
        res.lines.push_back(line);
    }
    return res;
}

/**
 * Crea un box per una frazione.
 */
LayoutBox make_fraction(const LayoutBox& num, const LayoutBox& den) {
    LayoutBox res;
    res.width = std::max(num.width, den.width) + 2;
    res.height = num.height + den.height + 1;
    res.baseline = num.height; // La linea della frazione è la nostra baseline
    
    // Numeratore
    for (int i = 0; i < num.height; ++i) {
        int pad = (res.width - num.width) / 2;
        res.lines.push_back(spaces(pad) + num.lines[i] + spaces(res.width - num.width - pad));
    }
    
    // Linea
    res.lines.push_back(dashes(res.width));
    
    // Denominatore
    for (int i = 0; i < den.height; ++i) {
        int pad = (res.width - den.width) / 2;
        res.lines.push_back(spaces(pad) + den.lines[i] + spaces(res.width - den.width - pad));
    }
    
    return res;
}
} // namespace

LayoutBox Ascii2DFormatter::layout(ExprPtr expr) {
    if (!expr) return {{"?"}, 1, 1, 0};

    return visit_expr(expr, [this](const auto& node) -> LayoutBox {
        using NodeT = std::decay_t<decltype(node)>;

        if constexpr (std::is_same_v<NodeT, IntegerLit>) {
            std::string s = node.value.decimal();
            return {{s}, static_cast<int>(s.length()), 1, 0};
        } else if constexpr (std::is_same_v<NodeT, RationalLit>) {
            return make_fraction(
                {{node.numerator.decimal()}, static_cast<int>(node.numerator.decimal().length()), 1, 0},
                {{node.denominator.decimal()}, static_cast<int>(node.denominator.decimal().length()), 1, 0}
            );
        } else if constexpr (std::is_same_v<NodeT, DecimalLit>) {
            return {{node.text}, static_cast<int>(node.text.length()), 1, 0};
        } else if constexpr (std::is_same_v<NodeT, Symbol>) {
            return {{node.name}, static_cast<int>(node.name.length()), 1, 0};
        } else if constexpr (std::is_same_v<NodeT, Constant>) {
            std::string s = (node.value == MathConstant::Pi ? "pi" : "e");
            return {{s}, static_cast<int>(s.length()), 1, 0};
        } else if constexpr (std::is_same_v<NodeT, Unary>) {
            LayoutBox op = layout(node.operand);
            LayoutBox sign = {{"-"}, 1, 1, 0};
            return join_horizontal(sign, op);
        } else if constexpr (std::is_same_v<NodeT, Binary>) {
            std::string op_str;
            switch (node.op) {
                case BinaryOp::Add: op_str = " + "; break;
                case BinaryOp::Sub: op_str = " - "; break;
                case BinaryOp::Mul: op_str = " * "; break;
                case BinaryOp::Mod: op_str = " % "; break;
                case BinaryOp::Pow: op_str = "^"; break;
                case BinaryOp::Equal: op_str = " = "; break;
                case BinaryOp::Div: return make_fraction(layout(node.left), layout(node.right));
            }
            return join_horizontal(layout(node.left), layout(node.right), op_str);
        } else if constexpr (std::is_same_v<NodeT, Sum>) {
            LayoutBox res = layout(node.terms[0]);
            for (size_t i = 1; i < node.terms.size(); ++i) {
                res = join_horizontal(res, layout(node.terms[i]), " + ");
            }
            return res;
        } else if constexpr (std::is_same_v<NodeT, Product>) {
            LayoutBox res = layout(node.factors[0]);
            for (size_t i = 1; i < node.factors.size(); ++i) {
                res = join_horizontal(res, layout(node.factors[i]), " * ");
            }
            return res;
        } else if constexpr (std::is_same_v<NodeT, FuncCall>) {
            if (node.func_id == BuiltinOp::Sqrt) {
                LayoutBox inner = layout(node.args[0]);
                LayoutBox res;
                res.width = inner.width + 3;
                res.height = inner.height + 1;
                res.baseline = inner.baseline + 1;
                res.lines.push_back("  " + dashes(inner.width + 1));
                std::string first_line = "\\/";
                for (int i = 0; i < inner.height; ++i) {
                    if (i == 0) res.lines.push_back("\\/" + inner.lines[i] + " ");
                    else res.lines.push_back("| " + inner.lines[i] + " ");
                }
                return res;
            }
            if (node.func_id == BuiltinOp::BesselZero) {
                std::string prefix = "BesselZero(";
                LayoutBox res = {{prefix}, static_cast<int>(prefix.length()), 1, 0};
                for (size_t i = 0; i < node.args.size(); ++i) {
                    if (i > 0) res = join_horizontal(res, {{", "}, 2, 1, 0});
                    res = join_horizontal(res, layout(node.args[i]));
                }
                res = join_horizontal(res, {{")"}, 1, 1, 0});
                return res;
            }
            // Fallback per altre funzioni
            std::string prefix = node.name + "(";
            LayoutBox res = {{prefix}, static_cast<int>(prefix.length()), 1, 0};
            for (size_t i = 0; i < node.args.size(); ++i) {
                if (i > 0) res = join_horizontal(res, {{", "}, 2, 1, 0});
                res = join_horizontal(res, layout(node.args[i]));
            }
            res = join_horizontal(res, {{")"}, 1, 1, 0});
            return res;
        }
        
        return {{"?"}, 1, 1, 0};
    });
}

std::string Ascii2DFormatter::format(ExprPtr expr) {
    LayoutBox box = layout(expr);
    std::string result;
    for (const auto& line : box.lines) {
        result += line + "\n";
    }
    return result;
}

} // namespace cas::formatter
