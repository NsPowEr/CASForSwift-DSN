#include "cas/parser.hpp"

#include "cas/error.hpp"

#include <sstream>

namespace cas {

namespace {

enum class ExprPrecedence : int {
    Lowest = 0,
    Sum = 10,
    Product = 20,
    Power = 30,
    Prefix = 40,
    Postfix = 50,
    Atom = 60,
};

[[nodiscard]] int precedence_of(ExprPtr expr) {
    if (!expr) {
        return static_cast<int>(ExprPrecedence::Atom);
    }

    switch (expr_kind(expr)) {
    case ExprKind::Sum:
        return static_cast<int>(ExprPrecedence::Sum);
    case ExprKind::Product:
        return static_cast<int>(ExprPrecedence::Product);
    case ExprKind::Binary: {
        const auto& binary = expr_ref<Binary>(expr);
        switch (binary.op) {
        case BinaryOp::Add:
        case BinaryOp::Sub:
            return static_cast<int>(ExprPrecedence::Sum);
        case BinaryOp::Mul:
        case BinaryOp::Div:
        case BinaryOp::Mod:
            return static_cast<int>(ExprPrecedence::Product);
        case BinaryOp::Pow:
            return static_cast<int>(ExprPrecedence::Power);
        case BinaryOp::Equal:
            return 10; // Lowest precedence for equality
        }
    }
    case ExprKind::Unary: {
        const auto& unary = expr_ref<Unary>(expr);
        return unary.op == UnaryOp::Factorial
            ? static_cast<int>(ExprPrecedence::Postfix)
            : static_cast<int>(ExprPrecedence::Prefix);
    }
    default:
        return static_cast<int>(ExprPrecedence::Atom);
    }
}

[[nodiscard]] CASError make_round_trip_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] Result<std::string> constant_text(MathConstant value) {
    switch (value) {
    case MathConstant::Pi:
        return ok(std::string("pi"));
    case MathConstant::E:
        return ok(std::string("e"));
    case MathConstant::I:
        return ok(std::string("i"));
    case MathConstant::Infinity:
        return ok(std::string("inf"));
    case MathConstant::NaN:
        return ok(std::string("nan"));
    }
    return fail<std::string>(make_round_trip_error(CASErrorKind::InternalError, "unknown constant"));
}

[[nodiscard]] Result<std::string> limit_direction_text(LimitDirection direction) {
    switch (direction) {
    case LimitDirection::Left:
        return ok(std::string("left"));
    case LimitDirection::Right:
        return ok(std::string("right"));
    case LimitDirection::Both:
        return ok(std::string("both"));
    }
    return fail<std::string>(make_round_trip_error(CASErrorKind::InternalError, "unknown limit direction"));
}

class RoundTripPrinter {
public:
    [[nodiscard]] Result<std::string> print(ExprPtr expr) const {
        return print_expr(expr, static_cast<int>(ExprPrecedence::Lowest));
    }

private:
    [[nodiscard]] Result<std::string> print_expr(ExprPtr expr, int parent_precedence) const {
        if (!expr) {
            return fail<std::string>(
                make_round_trip_error(CASErrorKind::InvalidArgument, "cannot print null expression"));
        }

        const auto self_precedence = precedence_of(expr);
        auto rendered = visit_expr(
            expr,
            [this, self_precedence](const auto& node) -> Result<std::string> {
                return print_node(node, self_precedence);
            });
        if (rendered.is_error()) {
            return rendered;
        }

        if (self_precedence < parent_precedence) {
            return ok(std::string("(") + rendered.value() + ")");
        }

        return rendered;
    }

    [[nodiscard]] Result<std::string> print_node(const IntegerLit& value, int) const {
        return ok(value.value.is_negative() ? "-" + value.value.decimal() : value.value.decimal());
    }

    [[nodiscard]] Result<std::string> print_node(const RationalLit& value, int) const {
        std::ostringstream out;
        if (value.numerator.is_negative()) {
            out << '-';
        }
        out << value.numerator.decimal() << '/';
        if (value.denominator.is_negative()) {
            out << '-';
        }
        out << value.denominator.decimal();
        return ok(std::move(out).str());
    }

    [[nodiscard]] Result<std::string> print_node(const DecimalLit& value, int) const {
        return ok(value.text);
    }

    [[nodiscard]] Result<std::string> print_node(const Symbol& value, int) const {
        return ok(value.name);
    }

    [[nodiscard]] Result<std::string> print_node(const Constant& value, int) const {
        return constant_text(value.value);
    }

    [[nodiscard]] Result<std::string> print_node(const Unary& value, int self_precedence) const {
        auto operand = print_expr(value.operand, self_precedence);
        if (operand.is_error()) {
            return operand;
        }
        if (value.op == UnaryOp::Neg) {
            return ok(std::string("-") + operand.value());
        }

        return ok(operand.value() + "!");
    }

    [[nodiscard]] Result<std::string> print_node(const Binary& value, int self_precedence) const {
        const auto left_precedence = value.op == BinaryOp::Pow ? self_precedence + 1 : self_precedence;
        const auto right_precedence =
            (value.op == BinaryOp::Sub || value.op == BinaryOp::Div || value.op == BinaryOp::Mod || value.op == BinaryOp::Pow)
            ? self_precedence + 1
            : self_precedence;

        const char* op = "";
        switch (value.op) {
        case BinaryOp::Add:
            op = "+";
            break;
        case BinaryOp::Sub:
            op = "-";
            break;
        case BinaryOp::Mul:
            op = "*";
            break;
        case BinaryOp::Div:
            op = "/";
            break;
        case BinaryOp::Pow:
            op = "^";
            break;
        case BinaryOp::Mod:
            op = "%";
            break;
        case BinaryOp::Equal:
            op = "=";
            break;
        }

        auto left = print_expr(value.left, left_precedence);
        if (left.is_error()) {
            return left;
        }
        auto right = print_expr(value.right, right_precedence);
        if (right.is_error()) {
            return right;
        }

        return ok(left.value() + op + right.value());
    }

    [[nodiscard]] Result<std::string> print_node(const FuncCall& value, int) const {
        std::ostringstream out;
        out << value.name << '(';
        for (std::size_t i = 0; i < value.args.size(); ++i) {
            if (i != 0U) {
                out << ',';
            }
            auto arg = print_expr(value.args[i], static_cast<int>(ExprPrecedence::Lowest));
            if (arg.is_error()) {
                return arg;
            }
            out << arg.value();
        }
        out << ')';
        return ok(std::move(out).str());
    }

    [[nodiscard]] Result<std::string> print_node(const Sum& value, int self_precedence) const {
        std::ostringstream out;
        for (std::size_t i = 0; i < value.terms.size(); ++i) {
            if (i != 0U) {
                out << '+';
            }
            auto term = print_expr(value.terms[i], self_precedence);
            if (term.is_error()) {
                return term;
            }
            out << term.value();
        }
        return ok(std::move(out).str());
    }

    [[nodiscard]] Result<std::string> print_node(const Product& value, int self_precedence) const {
        std::ostringstream out;
        for (std::size_t i = 0; i < value.factors.size(); ++i) {
            if (i != 0U) {
                out << '*';
            }
            auto factor = print_expr(value.factors[i], self_precedence);
            if (factor.is_error()) {
                return factor;
            }
            out << factor.value();
        }
        return ok(std::move(out).str());
    }

    [[nodiscard]] Result<std::string> print_node(const Integral& value, int) const {
        auto integrand = print_expr(value.integrand, static_cast<int>(ExprPrecedence::Lowest));
        if (integrand.is_error()) {
            return integrand;
        }
        std::ostringstream out;
        out << "int(" << integrand.value() << ',' << value.variable.name;
        if (value.lower.has_value()) {
            auto lower = print_expr(*value.lower, static_cast<int>(ExprPrecedence::Lowest));
            if (lower.is_error()) {
                return lower;
            }
            auto upper = print_expr(*value.upper, static_cast<int>(ExprPrecedence::Lowest));
            if (upper.is_error()) {
                return upper;
            }
            out << ',' << lower.value() << ',' << upper.value();
        }
        out << ')';
        return ok(std::move(out).str());
    }

    [[nodiscard]] Result<std::string> print_node(const Derivative& value, int) const {
        auto expression = print_expr(value.expression, static_cast<int>(ExprPrecedence::Lowest));
        if (expression.is_error()) {
            return expression;
        }
        std::ostringstream out;
        out << "diff(" << expression.value() << ',' << value.variable.name;
        if (value.order != 1U) {
            out << ',' << value.order;
        }
        out << ')';
        return ok(std::move(out).str());
    }

    [[nodiscard]] Result<std::string> print_node(const Limit& value, int) const {
        auto expression = print_expr(value.expression, static_cast<int>(ExprPrecedence::Lowest));
        if (expression.is_error()) {
            return expression;
        }
        auto point = print_expr(value.point, static_cast<int>(ExprPrecedence::Lowest));
        if (point.is_error()) {
            return point;
        }
        std::ostringstream out;
        out << "lim(" << expression.value() << ',' << value.variable.name << ',' << point.value();
        if (value.direction != LimitDirection::Both) {
            auto direction = limit_direction_text(value.direction);
            if (direction.is_error()) {
                return direction;
            }
            out << ',' << direction.value();
        }
        out << ')';
        return ok(std::move(out).str());
    }

    [[nodiscard]] Result<std::string> print_node(const RootOf& value, int) const {
        auto polynomial = print_expr(value.polynomial, static_cast<int>(ExprPrecedence::Lowest));
        if (polynomial.is_error()) {
            return polynomial;
        }
        std::ostringstream out;
        out << "RootOf(" << polynomial.value() << ',' << value.variable.name;
        if (value.root_index.has_value()) {
            out << ',' << *value.root_index;
        }
        out << ')';
        return ok(std::move(out).str());
    }

    [[nodiscard]] Result<std::string> print_node(const Matrix& value, int) const {
        std::ostringstream out;
        out << '[';
        for (std::size_t row = 0; row < value.rows; ++row) {
            if (row != 0U) {
                out << ',';
            }
            out << '[';
            for (std::size_t col = 0; col < value.cols; ++col) {
                if (col != 0U) {
                    out << ',';
                }
                auto element = print_expr(value.elements[row * value.cols + col], static_cast<int>(ExprPrecedence::Lowest));
                if (element.is_error()) {
                    return element;
                }
                out << element.value();
            }
            out << ']';
        }
        out << ']';
        return ok(std::move(out).str());
    }
};

}  // namespace

Result<std::string> to_round_trip_text(ExprPtr expr) {
    if (!expr) {
        return Result<std::string>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "round-trip printer requires a non-null expression",
            .hint = "pass a parsed AST node",
        });
    }

    return RoundTripPrinter{}.print(expr);
}

}  // namespace cas
