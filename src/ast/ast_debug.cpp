#include "cas/ast_debug.hpp"

#include <sstream>
#include <utility>

namespace cas {

namespace {

[[nodiscard]] std::string unary_op_name(UnaryOp op) {
    switch (op) {
    case UnaryOp::Neg:
        return "Neg";
    case UnaryOp::Factorial:
        return "Factorial";
    }

    return "UnknownUnary";
}

[[nodiscard]] std::string binary_op_name(BinaryOp op) {
    switch (op) {
    case BinaryOp::Add:
        return "Add";
    case BinaryOp::Sub:
        return "Sub";
    case BinaryOp::Mul:
        return "Mul";
    case BinaryOp::Div:
        return "Div";
    case BinaryOp::Pow:
        return "Pow";
    case BinaryOp::Mod:
        return "Mod";
    case BinaryOp::Equal:
        return "Equal";
    case BinaryOp::Less:
        return "Less";
    case BinaryOp::Greater:
        return "Greater";
    case BinaryOp::LessEqual:
        return "LessEqual";
    case BinaryOp::GreaterEqual:
        return "GreaterEqual";
    }

    return "UnknownBinary";
}

[[nodiscard]] std::string constant_name(MathConstant constant) {
    switch (constant) {
    case MathConstant::Pi:
        return "Pi";
    case MathConstant::E:
        return "E";
    case MathConstant::I:
        return "I";
    case MathConstant::Infinity:
        return "Infinity";
    case MathConstant::NegInfinity:
        return "NegInfinity";
    case MathConstant::ComplexInfinity:
        return "ComplexInfinity";
    case MathConstant::Indeterminate:
        return "Indeterminate";
    case MathConstant::NaN:
        return "NaN";
    case MathConstant::EulerGamma:
        return "EulerGamma";
    }

    return "UnknownConstant";
}

[[nodiscard]] std::string limit_direction_name(LimitDirection direction) {
    switch (direction) {
    case LimitDirection::Left:
        return "Left";
    case LimitDirection::Right:
        return "Right";
    case LimitDirection::Both:
        return "Both";
    }

    return "UnknownDirection";
}

[[nodiscard]] std::string join_expr_list(const std::vector<ExprPtr>& exprs);
[[nodiscard]] std::string debug_print_impl(ExprPtr expr);

[[nodiscard]] std::string optional_expr_string(const std::optional<ExprPtr>& expr) {
    return expr.has_value() ? debug_print_impl(*expr) : "None";
}

[[nodiscard]] std::string join_expr_list(const std::vector<ExprPtr>& exprs) {
    std::ostringstream out;
    for (std::size_t i = 0; i < exprs.size(); ++i) {
        if (i != 0U) {
            out << ", ";
        }
        out << debug_print_impl(exprs[i]);
    }
    return std::move(out).str();
}

struct DebugPrintVisitor {
    [[nodiscard]] std::string operator()(const IntegerLit& value) const {
        if (value.value.is_negative()) {
            return "IntegerLit(-" + value.value.decimal() + ")";
        }
        return "IntegerLit(" + value.value.decimal() + ")";
    }

    [[nodiscard]] std::string operator()(const RationalLit& value) const {
        std::ostringstream out;
        out << "RationalLit(";
        if (value.numerator.is_negative()) {
            out << '-';
        }
        out << value.numerator.decimal() << ", ";
        if (value.denominator.is_negative()) {
            out << '-';
        }
        out << value.denominator.decimal() << ')';
        return std::move(out).str();
    }
    [[nodiscard]] std::string operator()(const ComplexLit& value) const {
        return "ComplexLit(" + value.re_num.decimal() + "/" + value.re_den.decimal() + " + " +
               value.im_num.decimal() + "/" + value.im_den.decimal() + "i)";
    }


    [[nodiscard]] std::string operator()(const DecimalLit& value) const {
        return "DecimalLit(" + value.text + ")";
    }

    [[nodiscard]] std::string operator()(const Symbol& value) const {
        return "Symbol(" + value.name + ")";
    }

    [[nodiscard]] std::string operator()(const Constant& value) const {
        return "Constant(" + constant_name(value.value) + ")";
    }

    [[nodiscard]] std::string operator()(const Unary& value) const {
        return "Unary(" + unary_op_name(value.op) + ", " + debug_print_impl(value.operand) + ")";
    }

    [[nodiscard]] std::string operator()(const Binary& value) const {
        return "Binary(" + binary_op_name(value.op) + ", " + debug_print_impl(value.left) + ", " +
            debug_print_impl(value.right) + ")";
    }

    [[nodiscard]] std::string operator()(const FuncCall& value) const {
        return "FuncCall(" + value.name + ", [" + join_expr_list(value.args) + "])";
    }

    [[nodiscard]] std::string operator()(const Sum& value) const {
        return "Sum([" + join_expr_list(value.terms) + "])";
    }

    [[nodiscard]] std::string operator()(const Product& value) const {
        return "Product([" + join_expr_list(value.factors) + "])";
    }

    [[nodiscard]] std::string operator()(const Integral& value) const {
        return "Integral(" + debug_print_impl(value.integrand) + ", Symbol(" + value.variable.name + "), " +
            optional_expr_string(value.lower) + ", " + optional_expr_string(value.upper) + ")";
    }

    [[nodiscard]] std::string operator()(const Derivative& value) const {
        return "Derivative(" + debug_print_impl(value.expression) + ", Symbol(" + value.variable.name + "), " +
            std::to_string(value.order) + ")";
    }

    [[nodiscard]] std::string operator()(const Limit& value) const {
        return "Limit(" + debug_print_impl(value.expression) + ", Symbol(" + value.variable.name + "), " +
            debug_print_impl(value.point) + ", " + limit_direction_name(value.direction) + ")";
    }

    [[nodiscard]] std::string operator()(const RootOf& value) const {
        return "RootOf(" + debug_print_impl(value.polynomial) + ", Symbol(" + value.variable.name + "), " +
            (value.root_index.has_value() ? std::to_string(*value.root_index) : std::string("None")) + ")";
    }

    [[nodiscard]] std::string operator()(const Matrix& value) const {
        std::ostringstream out;
        out << "Matrix(" << value.rows << ", " << value.cols << ", [" << join_expr_list(value.elements) << "])";
        return std::move(out).str();
    }

    [[nodiscard]] std::string operator()(const SeriesExp& value) const {
        std::ostringstream out;
        out << "SeriesExp(" << value.var.name << ", point=" << debug_print_impl(value.point) << ", terms=[";
        for (std::size_t i = 0; i < value.terms.size(); ++i) {
            if (i > 0) out << ", ";
            out << "(" << value.terms[i].first << ", " << debug_print_impl(value.terms[i].second) << ")";
        }
        out << "], order=" << value.order << ")";
        return std::move(out).str();
    }

    [[nodiscard]] std::string operator()(const Quantity& value) const {
        std::ostringstream out;
        out << "Quantity(" << debug_print_impl(value.value) << ", SI(";
        bool first = true;
        auto append = [&](const char* name, int16_t exp) {
            if (exp == 0) return;
            if (!first) out << ", ";
            out << name << ":" << exp;
            first = false;
        };
        append("m", value.dimensions.m);
        append("kg", value.dimensions.kg);
        append("s", value.dimensions.s);
        append("A", value.dimensions.A);
        append("K", value.dimensions.K);
        append("mol", value.dimensions.mol);
        append("cd", value.dimensions.cd);
        out << "))";
        return std::move(out).str();
    }
};

[[nodiscard]] std::string debug_print_impl(ExprPtr expr) {
    if (!expr) {
        return "Null";
    }

    return visit_expr(expr, DebugPrintVisitor{});
}

}  // namespace

std::string debug_print(ExprPtr expr) {
    return debug_print_impl(expr);
}

}  // namespace cas
