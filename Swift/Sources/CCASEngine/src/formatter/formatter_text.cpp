#include "cas/formatter.hpp"
#include <sstream>
#include <algorithm>

namespace cas::formatter {

int TextFormatter::precedence(ExprPtr expr) {
    if (!expr) return 0;
    return visit_expr(expr, [](const auto& node) -> int {
        using NodeT = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<NodeT, IntegerLit> || std::is_same_v<NodeT, DecimalLit> || 
                      std::is_same_v<NodeT, Symbol> || std::is_same_v<NodeT, Constant> || 
                      std::is_same_v<NodeT, FuncCall> || std::is_same_v<NodeT, Matrix>) {
            return 100; // Atomici
        } else if constexpr (std::is_same_v<NodeT, Unary>) {
            return 80;
        } else if constexpr (std::is_same_v<NodeT, Binary>) {
            switch (node.op) {
                case BinaryOp::Pow: return 90;
                case BinaryOp::Mul:
                case BinaryOp::Div: return 70;
                case BinaryOp::Add:
                case BinaryOp::Sub: return 60;
                default: return 0;
            }
        } else if constexpr (std::is_same_v<NodeT, Sum>) {
            return 60;
        } else if constexpr (std::is_same_v<NodeT, Product>) {
            return 70;
        }
        return 0;
    });
}

std::string TextFormatter::format_child(ExprPtr parent, ExprPtr child, bool is_right) {
    int p_prec = precedence(parent);
    int c_prec = precedence(child);
    
    std::string s = format(child);
    
    // Parentesi se la precedenza del figlio è minore.
    // Se uguale, dipende dall'associatività (per semplicità mettiamo parentesi se is_right e non commutativo)
    if (c_prec < p_prec) return "(" + s + ")";
    
    // Caso speciale: a - (b + c)
    if (c_prec == p_prec && is_right) {
        if (const auto* b = expr_cast<Binary>(parent)) {
            if (b->op == BinaryOp::Sub || b->op == BinaryOp::Div) return "(" + s + ")";
        }
    }
    
    return s;
}

std::string TextFormatter::format(ExprPtr expr) {
    if (!expr) return "";

    return visit_expr(expr, [this, expr](const auto& node) -> std::string {
        using NodeT = std::decay_t<decltype(node)>;

        if constexpr (std::is_same_v<NodeT, IntegerLit>) {
            return node.value.decimal();
        } else if constexpr (std::is_same_v<NodeT, RationalLit>) {
            return node.numerator.decimal() + "/" + node.denominator.decimal();
        } else if constexpr (std::is_same_v<NodeT, DecimalLit>) {
            return node.text;
        } else if constexpr (std::is_same_v<NodeT, Symbol>) {
            return node.name;
        } else if constexpr (std::is_same_v<NodeT, Constant>) {
            switch (node.value) {
                case MathConstant::Pi: return "pi";
                case MathConstant::E: return "e";
                case MathConstant::Infinity: return "inf";
                default: return "const";
            }
        } else if constexpr (std::is_same_v<NodeT, Unary>) {
            return "-" + format_child(expr, node.operand);
        } else if constexpr (std::is_same_v<NodeT, Binary>) {
            std::string op_str;
            switch (node.op) {
                case BinaryOp::Add: op_str = " + "; break;
                case BinaryOp::Sub: op_str = " - "; break;
                case BinaryOp::Mul: op_str = " * "; break;
                case BinaryOp::Div: op_str = " / "; break;
                case BinaryOp::Mod: op_str = " % "; break;
                case BinaryOp::Pow: return format_child(expr, node.left) + "^" + format_child(expr, node.right, true);
            }
            return format_child(expr, node.left) + op_str + format_child(expr, node.right, true);
        } else if constexpr (std::is_same_v<NodeT, Sum>) {
            std::string s;
            for (size_t i = 0; i < node.terms.size(); ++i) {
                if (i > 0) s += " + ";
                s += format_child(expr, node.terms[i]);
            }
            return s;
        } else if constexpr (std::is_same_v<NodeT, Product>) {
            std::string s;
            for (size_t i = 0; i < node.factors.size(); ++i) {
                if (i > 0) s += " * ";
                s += format_child(expr, node.factors[i]);
            }
            return s;
        } else if constexpr (std::is_same_v<NodeT, FuncCall>) {
            std::string s = node.name + "(";
            for (size_t i = 0; i < node.args.size(); ++i) {
                if (i > 0) s += ", ";
                s += format(node.args[i]);
            }
            return s + ")";
        } else if constexpr (std::is_same_v<NodeT, Matrix>) {
            std::string s = "[";
            for (size_t r = 0; r < node.rows; ++r) {
                s += "[";
                for (size_t c = 0; c < node.cols; ++c) {
                    if (c > 0) s += ", ";
                    s += format(node.elements[r * node.cols + c]);
                }
                s += "]";
                if (r < node.rows - 1) s += ", ";
            }
            return s + "]";
        }
        return "?";
    });
}

} // namespace cas::formatter
