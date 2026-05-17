// L2-12: Pade approximant [m/n] of an analytic function at a finite centre.
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
// q_1..q_n with a Toeplitz coefficient matrix.  Solved exactly over Q via
// Gauss-Jordan elimination — the system has rational right-hand side
// because the Taylor coefficients we generate are reduced to Rational form
// before assembly.

#include "cas/calculus.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] BigInt factorial_bigint(unsigned int n) {
    BigInt acc(1);
    for (unsigned int k = 2; k <= n; ++k) acc = acc * BigInt(k);
    return acc;
}

[[nodiscard]] std::optional<Rational> exact_rational(ExprPtr expr) {
    if (!expr) return std::nullopt;
    if (const auto* il = expr_cast<IntegerLit>(expr)) return Rational(il->value);
    if (const auto* rl = expr_cast<RationalLit>(expr)) {
        return Rational(rl->numerator, rl->denominator);
    }
    if (const auto* un = expr_cast<Unary>(expr); un && un->op == UnaryOp::Neg) {
        auto inner = exact_rational(un->operand);
        if (inner.has_value()) return -inner.value();
    }
    return std::nullopt;
}

[[nodiscard]] Result<std::vector<Rational>> taylor_coefficients_rational(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int order,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    std::vector<Rational> coeffs;
    coeffs.reserve(static_cast<std::size_t>(order) + 1U);

    for (unsigned int k = 0; k <= order; ++k) {
        Result<ExprPtr> deriv = (k == 0U)
            ? ok(expr)
            : diff(expr, var, k, ctx);
        if (deriv.is_error()) return fail<std::vector<Rational>>(deriv.error());
        auto sub = ctx.substitute(deriv.value(), var, center);
        if (sub.is_error()) return fail<std::vector<Rational>>(sub.error());
        auto simp = ctx.simplify(sub.value());
        if (simp.is_error()) return fail<std::vector<Rational>>(simp.error());

        ExprPtr coeff_expr = simp.value();
        if (k >= 2U) {
            ExprPtr fact_expr = arena.make<IntegerLit>(factorial_bigint(k));
            ExprPtr divided = arena.make<Binary>(BinaryOp::Div, coeff_expr, fact_expr);
            auto divided_simp = ctx.simplify(divided);
            if (divided_simp.is_error()) return fail<std::vector<Rational>>(divided_simp.error());
            coeff_expr = divided_simp.value();
        }
        auto rational = exact_rational(coeff_expr);
        if (!rational.has_value()) {
            return fail<std::vector<Rational>>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Pade: Taylor coefficient at order " + std::to_string(k)
                         + " is not an exact rational; symbolic Pade requires rational data",
            });
        }
        coeffs.push_back(rational.value());
    }
    return ok(std::move(coeffs));
}

// Gauss-Jordan over Rational.  M is (rows × cols), rhs is rows.  On success
// returns x with M·x = rhs.  Returns Unimplemented on singular matrix.
[[nodiscard]] Result<std::vector<Rational>> solve_linear_rational(
    std::vector<std::vector<Rational>> M,
    std::vector<Rational> rhs) {
    const std::size_t n = M.size();
    if (n == 0U) return ok(std::vector<Rational>{});
    if (M[0].size() != n || rhs.size() != n) {
        return fail<std::vector<Rational>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Pade linear solve: matrix shape mismatch",
        });
    }
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot_row = col;
        while (pivot_row < n && M[pivot_row][col].numerator().is_zero()) ++pivot_row;
        if (pivot_row == n) {
            return fail<std::vector<Rational>>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Pade: Toeplitz system is singular (defective table entry)",
            });
        }
        if (pivot_row != col) {
            std::swap(M[col], M[pivot_row]);
            std::swap(rhs[col], rhs[pivot_row]);
        }
        const Rational pivot = M[col][col];
        for (std::size_t j = col; j < n; ++j) M[col][j] = M[col][j] / pivot;
        rhs[col] = rhs[col] / pivot;
        for (std::size_t r = 0; r < n; ++r) {
            if (r == col) continue;
            const Rational factor = M[r][col];
            if (factor.numerator().is_zero()) continue;
            for (std::size_t j = col; j < n; ++j) {
                M[r][j] = M[r][j] - factor * M[col][j];
            }
            rhs[r] = rhs[r] - factor * rhs[col];
        }
    }
    return ok(std::move(rhs));
}

[[nodiscard]] ExprPtr rational_to_expr(const Rational& r, AstArena& arena) {
    if (r.denominator() == BigInt(1)) return arena.make<IntegerLit>(r.numerator());
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

// Build  Σ_{k=0..coeffs.size()-1} coeffs[k] · (var − centre)^k  as an ExprPtr.
[[nodiscard]] ExprPtr build_poly_in_shifted_var(
    const std::vector<Rational>& coeffs,
    const Symbol& var,
    ExprPtr center,
    AstArena& arena) {
    ExprPtr x = arena.make<Symbol>(var);
    ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, x, center);
    std::vector<ExprPtr> terms;
    terms.reserve(coeffs.size());
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        if (coeffs[k].numerator().is_zero()) continue;
        ExprPtr coeff_expr = rational_to_expr(coeffs[k], arena);
        ExprPtr term;
        if (k == 0U) {
            term = coeff_expr;
        } else {
            ExprPtr power = (k == 1U)
                ? delta
                : arena.make<Binary>(BinaryOp::Pow, delta,
                    arena.make<IntegerLit>(BigInt(static_cast<long long>(k))));
            term = arena.make<Binary>(BinaryOp::Mul, coeff_expr, power);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return arena.make<IntegerLit>(BigInt(0));
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
    auto coeffs_res = taylor_coefficients_rational(expr, var, center, taylor_order, ctx);
    if (coeffs_res.is_error()) return fail<PadeApproximant>(coeffs_res.error());
    const std::vector<Rational>& c = coeffs_res.value();

    // Auxiliary indexing.  c_i for i > taylor_order is treated as 0 — never
    // needed in the system below because k − j ≤ k ≤ taylor_order = m + n
    // throughout.
    auto coeff = [&](long long i) -> Rational {
        if (i < 0 || static_cast<std::size_t>(i) >= c.size()) return Rational(BigInt(0));
        return c[static_cast<std::size_t>(i)];
    };

    const long long m = numerator_order;
    const long long n = denominator_order;

    std::vector<Rational> q_coeffs(static_cast<std::size_t>(n) + 1U,
                                   Rational(BigInt(0)));
    q_coeffs[0] = Rational(BigInt(1));

    if (n > 0) {
        // Toeplitz system: for k = m+1..m+n,
        //   Σ_{j=1..n} c_{k − j} · q_j = − c_k.
        std::vector<std::vector<Rational>> M(static_cast<std::size_t>(n),
            std::vector<Rational>(static_cast<std::size_t>(n), Rational(BigInt(0))));
        std::vector<Rational> rhs(static_cast<std::size_t>(n), Rational(BigInt(0)));
        for (long long i = 0; i < n; ++i) {
            const long long k = m + 1 + i;
            for (long long j = 1; j <= n; ++j) {
                M[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)] = coeff(k - j);
            }
            rhs[static_cast<std::size_t>(i)] = -coeff(k);
        }
        auto q_solve = solve_linear_rational(std::move(M), std::move(rhs));
        if (q_solve.is_error()) return fail<PadeApproximant>(q_solve.error());
        for (long long j = 0; j < n; ++j) {
            q_coeffs[static_cast<std::size_t>(j) + 1U] =
                q_solve.value()[static_cast<std::size_t>(j)];
        }
    }

    std::vector<Rational> p_coeffs(static_cast<std::size_t>(m) + 1U,
                                   Rational(BigInt(0)));
    for (long long k = 0; k <= m; ++k) {
        Rational acc(BigInt(0));
        for (long long j = 0; j <= std::min<long long>(k, n); ++j) {
            acc = acc + coeff(k - j) * q_coeffs[static_cast<std::size_t>(j)];
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
