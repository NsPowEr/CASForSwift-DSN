/// @file ast_kinds.hpp
/// @brief ExprKind enum, MathConstant, UnaryOp, BinaryOp, LimitDirection.
///
/// Split from ast.hpp (F8.0 / Task 1.1) — do NOT include ast.hpp from here;
/// this header must remain dependency-free (only <cstdint>).
#pragma once

#include <cstdint>

namespace cas {

// ---------------------------------------------------------------------------
// Expression kind discriminator
// ---------------------------------------------------------------------------

enum class ExprKind : std::uint8_t {
    Null,
    IntegerLit,
    RationalLit,
    DecimalLit,
    Symbol,
    Constant,
    Unary,
    Binary,
    FuncCall,
    Sum,
    Product,
    Integral,
    Derivative,
    Limit,
    RootOf,
    ComplexLit,

    Matrix,
    SeriesExp,
    Quantity,
};

// ---------------------------------------------------------------------------
// Mathematical constants
// ---------------------------------------------------------------------------

enum class MathConstant : std::uint8_t {
    Pi,
    E,
    I,
    Infinity,         // +∞ (canonical positive infinity)
    NegInfinity,      // -∞ (canonical signed negative infinity, F7.5.F1)
    ComplexInfinity,  // ∞̃  unsigned/directionless infinity, e.g. 1/0, Γ(0) (F7.5.F1)
    Indeterminate,    // 0·∞, ∞-∞, 0/0 indeterminate form (F7.5.F1)
    NaN,
    EulerGamma,
};

// ---------------------------------------------------------------------------
// Operator enumerations
// ---------------------------------------------------------------------------

enum class UnaryOp : std::uint8_t {
    Neg,
    Factorial,
};

enum class BinaryOp : std::uint8_t {
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Mod,
    Equal,
    Less,
    Greater,
    LessEqual,
    GreaterEqual,
};

enum class LimitDirection : std::uint8_t {
    Left,
    Right,
    Both,
};

}  // namespace cas
