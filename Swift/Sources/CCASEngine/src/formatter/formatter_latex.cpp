#include "cas/formatter.hpp"
#include <sstream>

namespace cas::formatter {

std::string LaTeXFormatter::format(ExprPtr expr) {
    if (!expr) return "";

    return visit_expr(expr, [this](const auto& node) -> std::string {
        using NodeT = std::decay_t<decltype(node)>;

        if constexpr (std::is_same_v<NodeT, IntegerLit>) {
            return node.value.decimal();
        } else if constexpr (std::is_same_v<NodeT, RationalLit>) {
            return "\\frac{" + node.numerator.decimal() + "}{" + node.denominator.decimal() + "}";
        } else if constexpr (std::is_same_v<NodeT, DecimalLit>) {
            return node.text;
        } else if constexpr (std::is_same_v<NodeT, Symbol>) {
            if (node.name == "alpha") return "\\alpha";
            if (node.name == "beta") return "\\beta";
            if (node.name == "gamma") return "\\gamma";
            if (node.name == "theta") return "\\theta";
            return node.name;
        } else if constexpr (std::is_same_v<NodeT, Constant>) {
            switch (node.value) {
                case MathConstant::Pi: return "\\pi";
                case MathConstant::E: return "e";
                case MathConstant::Infinity: return "\\infty";
                default: return "const";
            }
        } else if constexpr (std::is_same_v<NodeT, Unary>) {
            return "-" + format(node.operand);
        } else if constexpr (std::is_same_v<NodeT, Binary>) {
            switch (node.op) {
                case BinaryOp::Add: return format(node.left) + " + " + format(node.right);
                case BinaryOp::Sub: return format(node.left) + " - " + format(node.right);
                case BinaryOp::Mul: return format(node.left) + " \\cdot " + format(node.right);
                case BinaryOp::Div: return "\\frac{" + format(node.left) + "}{" + format(node.right) + "}";
                case BinaryOp::Mod: return format(node.left) + " \\pmod{" + format(node.right) + "}";
                case BinaryOp::Pow: return "{" + format(node.left) + "}^{" + format(node.right) + "}";
            }
            return "?";
        } else if constexpr (std::is_same_v<NodeT, Sum>) {
            std::string s;
            for (size_t i = 0; i < node.terms.size(); ++i) {
                if (i > 0) s += " + ";
                s += format(node.terms[i]);
            }
            return s;
        } else if constexpr (std::is_same_v<NodeT, Product>) {
            std::string s;
            for (size_t i = 0; i < node.factors.size(); ++i) {
                if (i > 0) s += " \\cdot ";
                s += format(node.factors[i]);
            }
            return s;
        } else if constexpr (std::is_same_v<NodeT, FuncCall>) {
            if (node.name == "sqrt") return "\\sqrt{" + format(node.args[0]) + "}";
            if (node.name == "sin" || node.name == "cos" || node.name == "tan" || node.name == "ln" || node.name == "exp") {
                return "\\" + node.name + "{" + format(node.args[0]) + "}";
            }
            std::string s = "\\text{" + node.name + "}(";
            for (size_t i = 0; i < node.args.size(); ++i) {
                if (i > 0) s += ", ";
                s += format(node.args[i]);
            }
            return s + ")";
        } else if constexpr (std::is_same_v<NodeT, Matrix>) {
            std::string s = "\\begin{pmatrix} ";
            for (size_t r = 0; r < node.rows; ++r) {
                for (size_t c = 0; c < node.cols; ++c) {
                    if (c > 0) s += " & ";
                    s += format(node.elements[r * node.cols + c]);
                }
                if (r < node.rows - 1) s += " \\\\ ";
            }
            return s + " \\end{pmatrix}";
        }
        return "?";
    });
}

} // namespace cas::formatter
