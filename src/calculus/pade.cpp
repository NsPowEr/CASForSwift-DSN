// L2-12 / F5.5: Pade approximant [m/n] of an analytic function at a finite
// centre.  Symbolic-coefficient version (B5 closure): Taylor coefficients
// remain in the symbolic algebra so closed forms involving algebraic
// constants (√2, π, e, RootOf, …) survive end-to-end without forcing the
// caller to pre-normalise to Q.
//
// Given f(x) = c_0 + c_1·(x − a) + c_2·(x − a)² + …  (Taylor series at a),
// the [m/n] Pade approximant is the rational function P(x)/Q(x) with
//   deg P ≤ m,  deg Q ≤ n,  Q(0) = 1,  P/Q ≡ f mod (x − a)^{m+n+1}.
//
// Writing P − f·Q ≡ 0 mod (x − a)^{m+n+1} yields:
//   degree k = 0..m   →   p_k = Σ_{j=0..min(k,n)} c_{k−j} · q_j
//   degree k = m+1..m+n → Σ_{j=0..min(k,n)} c_{k−j} · q_j = 0
//
// The second block, after substituting q_0 = 1, is a linear system in
// q_1..q_n with a Toeplitz coefficient matrix.  Solved exactly over the
// symbolic field via Gauss-Jordan elimination — pivot non-vanishing is
// decided by the simplifier; if the simplifier cannot prove pivot ≠ 0 a
// diagnostic Unimplemented is returned (Cat 4: explicit, with guidance).

#include "cas/calculus.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] BigInt factorial_bigint(unsigned int n) {
    BigInt acc(1);
    for (unsigned int k = 2; k <= n; ++k) acc = acc * BigInt(k);
    return acc;
}

[[nodiscard]] ExprPtr zero_expr(AstArena& arena) {
    return arena.make<IntegerLit>(BigInt(0));
}

[[nodiscard]] ExprPtr one_expr(AstArena& arena) {
    return arena.make<IntegerLit>(BigInt(1));
}

[[nodiscard]] bool is_literal_zero(ExprPtr e) {
    if (!e) return false;
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
    return false;
}

[[nodiscard]] Result<bool> is_known_zero(ExprPtr e, symbolic::CASContext& ctx) {
    if (is_literal_zero(e)) return ok(true);
    auto s = ctx.simplify(e);
    if (s.is_error()) return fail<bool>(s.error());
    return ok(is_literal_zero(s.value()));
}

[[nodiscard]] Result<ExprPtr> simp_mul(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    auto raw = ctx.arena().make<Binary>(BinaryOp::Mul, a, b);
    return ctx.simplify(raw);
}

[[nodiscard]] Result<ExprPtr> simp_sub(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    auto raw = ctx.arena().make<Binary>(BinaryOp::Sub, a, b);
    return ctx.simplify(raw);
}

[[nodiscard]] Result<ExprPtr> simp_div(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    auto raw = ctx.arena().make<Binary>(BinaryOp::Div, a, b);
    return ctx.simplify(raw);
}

[[nodiscard]] Result<ExprPtr> simp_neg(ExprPtr a, symbolic::CASContext& ctx) {
    auto raw = ctx.arena().make<Unary>(UnaryOp::Neg, a);
    return ctx.simplify(raw);
}

// Compute c_k = (1/k!) · d^k f / dx^k |_{x = centre}, kept fully symbolic.
[[nodiscard]] Result<std::vector<ExprPtr>> taylor_coefficients_symbolic(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int order,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> coeffs;
    coeffs.reserve(static_cast<std::size_t>(order) + 1U);

    for (unsigned int k = 0; k <= order; ++k) {
        Result<ExprPtr> deriv = (k == 0U)
            ? ok(expr)
            : diff(expr, var, k, ctx);
        if (deriv.is_error()) return fail<std::vector<ExprPtr>>(deriv.error());
        auto sub = ctx.substitute(deriv.value(), var, center);
        if (sub.is_error()) return fail<std::vector<ExprPtr>>(sub.error());
        auto simp = ctx.simplify(sub.value());
        if (simp.is_error()) return fail<std::vector<ExprPtr>>(simp.error());

        ExprPtr coeff_expr = simp.value();
        if (k >= 2U) {
            ExprPtr fact_expr = arena.make<IntegerLit>(factorial_bigint(k));
            auto divided = simp_div(coeff_expr, fact_expr, ctx);
            if (divided.is_error()) return fail<std::vector<ExprPtr>>(divided.error());
            coeff_expr = divided.value();
        }
        coeffs.push_back(coeff_expr);
    }
    return ok(std::move(coeffs));
}

// Symbolic Gauss-Jordan over the field of expressions.  Pivot non-vanishing
// is decided by ctx.simplify; if the simplifier cannot decide we treat the
// candidate as zero (skip), and emit a clear Unimplemented if no pivot is
// available in the column.  Returns the solution vector on success.
[[nodiscard]] Result<std::vector<ExprPtr>> solve_linear_symbolic(
    std::vector<std::vector<ExprPtr>> M,
    std::vector<ExprPtr> rhs,
    symbolic::CASContext& ctx) {
    const std::size_t n = M.size();
    if (n == 0U) return ok(std::vector<ExprPtr>{});
    if (M[0].size() != n || rhs.size() != n) {
        return fail<std::vector<ExprPtr>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Pade linear solve: matrix shape mismatch",
        });
    }
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot_row = col;
        while (pivot_row < n) {
            auto zr = is_known_zero(M[pivot_row][col], ctx);
            if (zr.is_error()) return fail<std::vector<ExprPtr>>(zr.error());
            if (!zr.value()) break;
            ++pivot_row;
        }
        if (pivot_row == n) {
            return fail<std::vector<ExprPtr>>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Pade: Toeplitz system singular or simplifier cannot "
                           "decide pivot non-vanishing; provide assumptions or "
                           "reduce the requested order",
            });
        }
        if (pivot_row != col) {
            std::swap(M[col], M[pivot_row]);
            std::swap(rhs[col], rhs[pivot_row]);
        }
        const ExprPtr pivot = M[col][col];
        for (std::size_t j = col; j < n; ++j) {
            auto r = simp_div(M[col][j], pivot, ctx);
            if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
            M[col][j] = r.value();
        }
        {
            auto r = simp_div(rhs[col], pivot, ctx);
            if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
            rhs[col] = r.value();
        }
        for (std::size_t r = 0; r < n; ++r) {
            if (r == col) continue;
            const ExprPtr factor = M[r][col];
            auto zr = is_known_zero(factor, ctx);
            if (zr.is_error()) return fail<std::vector<ExprPtr>>(zr.error());
            if (zr.value()) continue;
            for (std::size_t j = col; j < n; ++j) {
                auto prod = simp_mul(factor, M[col][j], ctx);
                if (prod.is_error()) return fail<std::vector<ExprPtr>>(prod.error());
                auto diff_e = simp_sub(M[r][j], prod.value(), ctx);
                if (diff_e.is_error()) return fail<std::vector<ExprPtr>>(diff_e.error());
                M[r][j] = diff_e.value();
            }
            auto prod = simp_mul(factor, rhs[col], ctx);
            if (prod.is_error()) return fail<std::vector<ExprPtr>>(prod.error());
            auto diff_e = simp_sub(rhs[r], prod.value(), ctx);
            if (diff_e.is_error()) return fail<std::vector<ExprPtr>>(diff_e.error());
            rhs[r] = diff_e.value();
        }
    }
    return ok(std::move(rhs));
}

// Build  Σ_{k=0..coeffs.size()-1} coeffs[k] · (var − centre)^k  as an ExprPtr.
[[nodiscard]] ExprPtr build_poly_in_shifted_var(
    const std::vector<ExprPtr>& coeffs,
    const Symbol& var,
    ExprPtr center,
    AstArena& arena) {
    ExprPtr x = arena.make<Symbol>(var);
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, x, center);
    std::vector<ExprPtr> terms;
    terms.reserve(coeffs.size());
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        if (is_literal_zero(coeffs[k])) continue;
        ExprPtr term;
        if (k == 0U) {
            term = coeffs[k];
        } else {
            ExprPtr power = (k == 1U)
                ? delta
                : arena.make<Binary>(BinaryOp::Pow, delta,
                    arena.make<IntegerLit>(BigInt(static_cast<long long>(k))));
            term = arena.make<Binary>(BinaryOp::Mul, coeffs[k], power);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return zero_expr(arena);
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

}  // namespace

Result<PadeApproximant> pade_approximant(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int numerator_order,
    unsigned int denominator_order,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    const unsigned int taylor_order = numerator_order + denominator_order;
    auto coeffs_res = taylor_coefficients_symbolic(expr, var, center, taylor_order, ctx);
    if (coeffs_res.is_error()) return fail<PadeApproximant>(coeffs_res.error());
    const std::vector<ExprPtr>& c = coeffs_res.value();

    auto coeff = [&](long long i) -> ExprPtr {
        if (i < 0 || static_cast<std::size_t>(i) >= c.size()) return zero_expr(arena);
        return c[static_cast<std::size_t>(i)];
    };

    const long long m = numerator_order;
    const long long n = denominator_order;

    std::vector<ExprPtr> q_coeffs(static_cast<std::size_t>(n) + 1U,
                                  zero_expr(arena));
    q_coeffs[0] = one_expr(arena);

    if (n > 0) {
        // Toeplitz system: for k = m+1..m+n,
        //   Σ_{j=1..n} c_{k − j} · q_j = − c_k.
        std::vector<std::vector<ExprPtr>> M(static_cast<std::size_t>(n),
            std::vector<ExprPtr>(static_cast<std::size_t>(n), zero_expr(arena)));
        std::vector<ExprPtr> rhs(static_cast<std::size_t>(n), zero_expr(arena));
        for (long long i = 0; i < n; ++i) {
            const long long k = m + 1 + i;
            for (long long j = 1; j <= n; ++j) {
                M[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)] = coeff(k - j);
            }
            auto neg = simp_neg(coeff(k), ctx);
            if (neg.is_error()) return fail<PadeApproximant>(neg.error());
            rhs[static_cast<std::size_t>(i)] = neg.value();
        }
        auto q_solve = solve_linear_symbolic(std::move(M), std::move(rhs), ctx);
        if (q_solve.is_error()) return fail<PadeApproximant>(q_solve.error());
        for (long long j = 0; j < n; ++j) {
            q_coeffs[static_cast<std::size_t>(j) + 1U] =
                q_solve.value()[static_cast<std::size_t>(j)];
        }
    }

    std::vector<ExprPtr> p_coeffs(static_cast<std::size_t>(m) + 1U,
                                  zero_expr(arena));
    for (long long k = 0; k <= m; ++k) {
        ExprPtr acc = zero_expr(arena);
        for (long long j = 0; j <= std::min<long long>(k, n); ++j) {
            auto prod = simp_mul(coeff(k - j), q_coeffs[static_cast<std::size_t>(j)], ctx);
            if (prod.is_error()) return fail<PadeApproximant>(prod.error());
            auto raw_sum = arena.make<Binary>(BinaryOp::Add, acc, prod.value());
            auto sum_simp = ctx.simplify(raw_sum);
            if (sum_simp.is_error()) return fail<PadeApproximant>(sum_simp.error());
            acc = sum_simp.value();
        }
        p_coeffs[static_cast<std::size_t>(k)] = acc;
    }

    PadeApproximant out;
    out.center = center;
    out.numerator_order = numerator_order;
    out.denominator_order = denominator_order;
    out.numerator = build_poly_in_shifted_var(p_coeffs, var, center, arena);
    out.denominator = build_poly_in_shifted_var(q_coeffs, var, center, arena);
    return ok(std::move(out));
}

}  // namespace cas::calculus
