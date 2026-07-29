// Bronstein, "Symbolic Integration I", §7.1 — ConstantSystem (Fig 7.1) and its
// correctness (Lemma 7.1.2).  Reduces a homogeneous linear system A·c = 0 whose
// entries live in a differential field K to the equivalent system for its
// CONSTANT solutions (c ∈ Const(K)), returning a null-space basis over Q.
//
// This is the general replacement for the Q(x)-only "clear denominators + equate
// x-power coefficients" shortcut that previously guarded the residual step of the
// parametric PolyRischDE non-cancellation solver (risch_param_nocancel.cpp).  With
// a proper tower derivation (DifferentialField::derive), the same algorithm now
// handles deep towers K = k(t_1,…,t_j) ⊋ Q(x): the derivation D itself supplies
// the extra rows R_{m+1} = D(R_i)/D(a_ij) that turn a non-constant column into a
// constant one, so no assumption K = Q(x) is baked in (closes HC-F8-PENDING-17).
//
// Algorithm (verbatim Fig 7.1, homogeneous case u = 0 ⇒ v = 0):
//   1. RowEchelon(A) over K.
//   2. while A is not all-constant:
//        j ← minimal column with a non-constant entry
//        i ← any row with a_ij ∉ Const(K)   (then D(a_ij) ≠ 0)
//        R_{m+1} ← D(R_i)/D(a_ij)            (first nonzero entry is 1 in col j)
//        for every row s:  R_s ← R_s − a_sj·R_{m+1}   (clears column j)
//        append R_{m+1}
//   Each pass makes columns 1..j constant with j strictly increasing, so it
//   terminates in ≤ m passes.  Every returned candidate is independently verified
//   downstream by field back-substitution, so this can only affect completeness.

#include "calculus_internal.hpp"
#include "risch_parametric_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

using detail::as_rational;
using detail::null_space_basis;
using detail::row_echelon;

// Rational-normalise an entry of K (common denominator + simplify) so the zero
// test and constant test are reliable.
[[nodiscard]] ExprPtr norm(ExprPtr e, symbolic::CASContext& ctx) {
    auto tog = algebra::together(e, ctx);
    ExprPtr r = tog.is_ok() ? tog.value() : e;
    if (auto s = ctx.simplify(r); s.is_ok()) r = s.value();
    return r;
}

[[nodiscard]] bool is_zero_expr(ExprPtr e) {
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
    return false;
}

[[nodiscard]] ExprPtr k_sub(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    return norm(ctx.arena().make<Binary>(BinaryOp::Sub, a, b), ctx);
}
[[nodiscard]] ExprPtr k_mul(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    return norm(ctx.arena().make<Binary>(BinaryOp::Mul, a, b), ctx);
}
[[nodiscard]] ExprPtr k_div(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    return norm(ctx.arena().make<Binary>(BinaryOp::Div, a, b), ctx);
}

// D(e) ≡ 0 in the field ⇔ e ∈ Const(K).
[[nodiscard]] Result<bool> is_constant(ExprPtr e, const DifferentialField& field,
                                       symbolic::CASContext& ctx) {
    auto d = field.derive(e, ctx);
    if (d.is_error()) return fail<bool>(d.error());
    return ok(is_zero_expr(norm(d.value(), ctx)));
}

// Reduced row-echelon form of A (m columns) over K, in place.
void field_row_echelon(std::vector<std::vector<ExprPtr>>& A, std::size_t m,
                       symbolic::CASContext& ctx) {
    std::size_t row = 0;
    for (std::size_t col = 0; col < m && row < A.size(); ++col) {
        std::size_t best = row;
        while (best < A.size() && is_zero_expr(A[best][col])) ++best;
        if (best == A.size()) continue;  // no pivot in this column
        if (best != row) std::swap(A[best], A[row]);
        ExprPtr pivot = A[row][col];
        for (std::size_t c = col; c < m; ++c) A[row][c] = k_div(A[row][c], pivot, ctx);
        for (std::size_t r = 0; r < A.size(); ++r) {
            if (r == row) continue;
            ExprPtr f = A[r][col];
            if (is_zero_expr(f)) continue;
            for (std::size_t c = col; c < m; ++c)
                A[r][c] = k_sub(A[r][c], k_mul(f, A[row][c], ctx), ctx);
        }
        ++row;
    }
}

}  // namespace

Result<std::vector<std::vector<Rational>>>
constant_system_nullspace(std::vector<std::vector<ExprPtr>> A, std::size_t m,
                          const DifferentialField& field, symbolic::CASContext& ctx) {
    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<std::vector<std::vector<Rational>>>(
            "calculus", "constant_system_nullspace", msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Bronstein Symbolic Integration I §7.1 ConstantSystem / Lemma 7.1.2",
            "HC-F8-PENDING-17");
    };

    for (auto& r : A)
        for (auto& e : r) e = norm(e, ctx);

    // Step 1: row-echelon over K.
    field_row_echelon(A, m, ctx);

    // Step 2: clear non-constant columns using derivation-generated rows.  Bounded
    // by ≤ m passes (each pass makes one more leading column constant).
    for (std::size_t pass = 0; pass <= m + 1U; ++pass) {
        // Locate the minimal non-constant column j and a witness row i.
        std::size_t j = m, iw = 0;
        bool found = false;
        for (std::size_t col = 0; col < m && !found; ++col) {
            for (std::size_t r = 0; r < A.size(); ++r) {
                auto c = is_constant(A[r][col], field, ctx);
                if (c.is_error()) return fail<std::vector<std::vector<Rational>>>(c.error());
                if (!c.value()) { j = col; iw = r; found = true; break; }
            }
        }
        if (!found) break;  // A is all-constant

        // R_{m+1} = D(R_i)/D(a_ij).  a_ij non-constant ⇒ D(a_ij) ≠ 0.
        auto dpiv_res = field.derive(A[iw][j], ctx);
        if (dpiv_res.is_error()) return fail<std::vector<std::vector<Rational>>>(dpiv_res.error());
        ExprPtr dpiv = norm(dpiv_res.value(), ctx);
        if (is_zero_expr(dpiv)) return fail_unimpl("degenerate derivation of a non-constant entry");

        std::vector<ExprPtr> Rnew(m);
        for (std::size_t c = 0; c < m; ++c) {
            auto dc = field.derive(A[iw][c], ctx);
            if (dc.is_error()) return fail<std::vector<std::vector<Rational>>>(dc.error());
            Rnew[c] = k_div(norm(dc.value(), ctx), dpiv, ctx);
        }

        // Eliminate column j from every existing row: R_s ← R_s − a_sj·R_{m+1}.
        for (auto& Rs : A) {
            ExprPtr asj = Rs[j];
            if (is_zero_expr(asj)) continue;
            for (std::size_t c = 0; c < m; ++c)
                Rs[c] = k_sub(Rs[c], k_mul(asj, Rnew[c], ctx), ctx);
        }
        A.push_back(std::move(Rnew));

        if (pass == m + 1U) return fail_unimpl("ConstantSystem did not terminate");
    }

    // All entries are now constant: convert to a rational matrix, dropping zero
    // rows.  A surviving non-rational constant would mean Const(K) ⊋ Q.
    std::vector<std::vector<Rational>> B;
    for (auto& r : A) {
        std::vector<Rational> row(m, Rational(BigInt(0)));
        bool nonzero = false;
        for (std::size_t c = 0; c < m; ++c) {
            auto q = as_rational(r[c]);
            if (!q) { if (auto s = ctx.simplify(r[c]); s.is_ok()) q = as_rational(s.value()); }
            if (!q) return fail_unimpl("reduced constant is not rational over Q (Const(K) ⊋ Q)");
            row[c] = *q;
            if (!q->numerator().is_zero()) nonzero = true;
        }
        if (nonzero) B.push_back(std::move(row));
    }

    // Empty constraint set ⇒ every c_i is free.
    if (B.empty()) {
        std::vector<std::vector<Rational>> basis;
        basis.reserve(m);
        for (std::size_t c = 0; c < m; ++c) {
            std::vector<Rational> e(m, Rational(BigInt(0)));
            e[c] = Rational(BigInt(1));
            basis.push_back(std::move(e));
        }
        return ok(std::move(basis));
    }

    auto pivots = row_echelon(B, m);
    return ok(null_space_basis(B, pivots, m));
}

}  // namespace cas::calculus
