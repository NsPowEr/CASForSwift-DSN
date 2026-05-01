#pragma once

#include "cas/ast.hpp"

namespace cas::detail {

enum class Prec : int {
    Lowest  = 0,
    Sum     = 10,
    Product = 20,
    Power   = 30,
    Prefix  = 40,
    Postfix = 50,
    Atom    = 60,
};

[[nodiscard]] inline int expr_prec(ExprPtr expr) noexcept {
    if (!expr) {
        return static_cast<int>(Prec::Atom);
    }
    switch (expr_kind(expr)) {
    case ExprKind::Sum:
        return static_cast<int>(Prec::Sum);
    case ExprKind::Product:
        return static_cast<int>(Prec::Product);
    case ExprKind::Binary: {
        const auto& b = expr_ref<Binary>(expr);
        switch (b.op) {
        case BinaryOp::Add: case BinaryOp::Sub:
            return static_cast<int>(Prec::Sum);
        case BinaryOp::Mul: case BinaryOp::Div: case BinaryOp::Mod:
            return static_cast<int>(Prec::Product);
        case BinaryOp::Pow:
            return static_cast<int>(Prec::Power);
        }
        return static_cast<int>(Prec::Atom);
    }
    case ExprKind::Unary: {
        const auto& u = expr_ref<Unary>(expr);
        return u.op == UnaryOp::Factorial
            ? static_cast<int>(Prec::Postfix)
            : static_cast<int>(Prec::Prefix);
    }
    default:
        return static_cast<int>(Prec::Atom);
    }
}

}  // namespace cas::detail
