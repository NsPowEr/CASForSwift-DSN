// Shared primitives for the parametric Risch DE solvers over Q(x)
// (risch_parametric.cpp and risch_parametric_rational.cpp).  Kept here so the
// two translation units share one definition (anti-monolith split of the
// originally single file; no duplicated logic).
#pragma once

#include "cas/ast.hpp"
#include "cas/bigint.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus::detail {

// Extract a Rational from an integer/rational literal or a Neg thereof.
[[nodiscard]] inline std::optional<Rational> as_rational(ExprPtr e) {
    if (const auto* lit = expr_cast<IntegerLit>(e)) return Rational(lit->value);
    if (const auto* lit = expr_cast<RationalLit>(e))
        return Rational(lit->numerator, lit->denominator);
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (auto inner = as_rational(un->operand)) return -*inner;
    }
    return std::nullopt;
}

// Reduce M to (non-reduced) row-echelon form in place; return pivot columns.
// rank == pivots.size().
[[nodiscard]] inline std::vector<std::size_t> row_echelon(
    std::vector<std::vector<Rational>>& M, std::size_t n_cols) {
    std::vector<std::size_t> pivots;
    std::size_t row = 0;
    const std::size_t n_rows = M.size();
    for (std::size_t col = 0; col < n_cols && row < n_rows; ++col) {
        std::size_t best = row;
        while (best < n_rows && M[best][col].numerator().is_zero()) ++best;
        if (best == n_rows) continue;  // free column (no pivot)
        if (best != row) std::swap(M[best], M[row]);
        Rational pivot = M[row][col];
        for (std::size_t c = col; c < n_cols; ++c) M[row][c] = M[row][c] / pivot;
        for (std::size_t r = 0; r < n_rows; ++r) {
            if (r == row) continue;
            Rational f = M[r][col];
            if (f.numerator().is_zero()) continue;
            for (std::size_t c = col; c < n_cols; ++c)
                M[r][c] = M[r][c] - f * M[row][c];
        }
        pivots.push_back(col);
        ++row;
    }
    return pivots;
}

// Null-space basis of the homogeneous system M·x = 0 given the pivot columns
// from row_echelon.  Each free column yields one basis vector.
[[nodiscard]] inline std::vector<std::vector<Rational>> null_space_basis(
    const std::vector<std::vector<Rational>>& M,
    const std::vector<std::size_t>& pivots,
    std::size_t n_cols) {
    std::vector<bool> is_pivot(n_cols, false);
    for (auto p : pivots) is_pivot[p] = true;
    std::vector<std::size_t> free_cols;
    for (std::size_t c = 0; c < n_cols; ++c) if (!is_pivot[c]) free_cols.push_back(c);

    std::vector<std::vector<Rational>> basis;
    basis.reserve(free_cols.size());
    for (std::size_t f : free_cols) {
        std::vector<Rational> v(n_cols, Rational(BigInt(0)));
        v[f] = Rational(BigInt(1));
        for (std::size_t r = 0; r < pivots.size(); ++r) v[pivots[r]] = -M[r][f];
        basis.push_back(std::move(v));
    }
    return basis;
}

// Rational constant as an ExprPtr literal (IntegerLit when denominator is 1).
[[nodiscard]] inline ExprPtr rational_to_expr(const Rational& c, AstArena& arena) {
    return (c.denominator() == BigInt(1))
        ? static_cast<ExprPtr>(arena.make<IntegerLit>(c.numerator()))
        : static_cast<ExprPtr>(arena.make<RationalLit>(c.numerator(), c.denominator()));
}

// Q-coefficient vector of a polynomial-in-var ExprPtr (nullopt if not a
// polynomial over Q in var).
[[nodiscard]] inline std::optional<std::vector<Rational>>
poly_coeffs_q(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    auto pr = algebra::parse_polynomial(e, var, ctx);
    if (pr.is_error()) return std::nullopt;
    std::vector<Rational> out;
    out.reserve(pr.value().size());
    for (ExprPtr c : pr.value().coefficients()) {
        auto r = as_rational(c);
        if (!r) { if (auto s = ctx.simplify(c); s.is_ok()) r = as_rational(s.value()); }
        if (!r) return std::nullopt;
        out.push_back(*r);
    }
    return out;
}

}  // namespace cas::calculus::detail
