#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "../symbolic/summation_gosper.hpp"
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace cas::calculus {

[[nodiscard]] static bool is_one(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == BigInt(1);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == BigInt(1) && rational->denominator == BigInt(1);
    }
    return false;
}

[[nodiscard]] static bool is_positive_infinity(ExprPtr expr) {
    const auto* constant = expr_cast<Constant>(expr);
    return constant != nullptr && constant->value == MathConstant::Infinity;
}

[[nodiscard]] static std::optional<unsigned int> positive_integer_u32(ExprPtr expr) {
    const auto* integer = expr_cast<IntegerLit>(expr);
    if (integer == nullptr || integer->value <= BigInt(0)) {
        return std::nullopt;
    }
    if (integer->value.bit_length() > std::numeric_limits<unsigned int>::digits) {
        return std::nullopt;
    }
    return static_cast<unsigned int>(integer->value.to_u64());
}

[[nodiscard]] static BigInt factorial(unsigned int value) {
    BigInt result(1);
    for (unsigned int factor = 2U; factor <= value; ++factor) {
        result *= BigInt(static_cast<std::int64_t>(factor));
    }
    return result;
}

[[nodiscard]] static BigInt pow_bigint_nonnegative(BigInt base, unsigned int exponent) {
    BigInt result(1);
    while (exponent > 0U) {
        if ((exponent % 2U) == 1U) {
            result *= base;
        }
        exponent /= 2U;
        if (exponent > 0U) {
            base *= base;
        }
    }
    return result;
}

[[nodiscard]] static ExprPtr rational_expr(AstArena& arena, const Rational& value) {
    if (value.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] static ExprPtr pow_expr(AstArena& arena, ExprPtr base, unsigned int exponent) {
    if (exponent == 0U) return arena.make<IntegerLit>(BigInt(1));
    if (exponent == 1U) return base;
    return arena.make<Binary>(BinaryOp::Pow, base, arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(exponent))));
}

[[nodiscard]] static std::optional<unsigned int> reciprocal_power_exponent(ExprPtr term, const Symbol& var) {
    const auto* division = expr_cast<Binary>(term);
    if (division == nullptr || division->op != BinaryOp::Div || !is_one(division->left)) {
        return std::nullopt;
    }

    const auto* power = expr_cast<Binary>(division->right);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return std::nullopt;
    }

    const auto* symbol = expr_cast<Symbol>(power->left);
    if (symbol == nullptr || symbol->name != var.name) {
        return std::nullopt;
    }

    return positive_integer_u32(power->right);
}

[[nodiscard]] static Result<ExprPtr> zeta_even_value(unsigned int exponent, symbolic::CASContext& ctx) {
    if (exponent == 0U || (exponent % 2U) != 0U) {
        // F0.8-MIGRATED
        return make_unimplemented<ExprPtr>(
            "calculus", "zeta_even_value",
            "p-series with odd or zero exponent",
            cas::error::reason_codes::SUMMATION_GENERAL,
            "Implement zeta at odd arguments (open problem) or return unevaluated Zeta(n) form",
            "F0.8");
    }

    const unsigned int m = exponent / 2U;
    const auto bernoulli = cas::numtheory::bernoulli_numbers(exponent);
    Rational coefficient = bernoulli[exponent] * Rational(pow_bigint_nonnegative(BigInt(2), exponent - 1U), factorial(exponent));
    if ((m % 2U) == 0U) {
        coefficient = -coefficient;
    }

    AstArena& arena = ctx.arena();
    ExprPtr pi_power = pow_expr(arena, arena.make<Constant>(MathConstant::Pi), exponent);
    ExprPtr coeff = rational_expr(arena, coefficient);
    if (const auto* integer = expr_cast<IntegerLit>(coeff); integer != nullptr && integer->value == BigInt(1)) {
        return ok(pi_power);
    }
    return ctx.simplify(arena.make<Binary>(BinaryOp::Mul, coeff, pi_power));
}

// F5.7 sub-block 0 — definite hypergeometric summation via Gosper.
//
// gosper_sum(t, k) returns (when it succeeds) an antidifference S(k) such
// that S(k+1) − S(k) = t(k).  The Newton-Leibniz analogue for the finite
// calculus then gives  Σ_{k=a}^{b} t(k) = S(b+1) − S(a).
//
// Indeterminate bounds (free symbols, RootOf, …) are passed verbatim into
// the substitution machinery; the simplifier folds them when possible and
// leaves them symbolic otherwise.
[[nodiscard]] static Result<ExprPtr> try_gosper_definite(
    ExprPtr term,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto antidiff = symbolic::gosper_sum(term, var, ctx);
    if (antidiff.is_error()) return fail<ExprPtr>(antidiff.error());
    if (!antidiff.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Gosper: term is not Gosper-summable",
        });
    }
    AstArena& arena = ctx.arena();
    ExprPtr S = antidiff.value().value();
    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));
    auto S_upper = ctx.substitute(S, var, upper_plus_one);
    if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
    auto S_lower = ctx.substitute(S, var, lower);
    if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
    ExprPtr diff = arena.make<Binary>(BinaryOp::Sub, S_upper.value(), S_lower.value());
    return ctx.simplify(diff);
}

// F5.7 sub-block 1 — Abramov rational summation via polygamma antidifference.
//
// Polygamma identity (Abramowitz & Stegun 6.4.6):
//   ψ^(n)(k+1) − ψ^(n)(k) = (−1)^n · n! / k^(n+1).
// So for a single rational atom  A / (k + a)^m  (m ≥ 1) the antidifference
// is
//   S(k) = A · (−1)^(m−1) / (m − 1)! · ψ^(m−1)(k + a),
// where ψ^0 = digamma.  The definite sum follows by  S(b+1) − S(a).
//
// We pattern-match the input term against this shape post-`together`,
// covering all univariate-fraction rational summands of the form
// `A / (linear(k))^m`.  Multi-atom partial fractions reduce to a finite sum
// of these shapes; full Abramov decomposition over Q[k] is the natural
// extension and is left as a follow-on (the current path correctly delegates
// to the diagnostic Unimplemented when the term is not of single-atom form).
[[nodiscard]] static bool extract_rational_const(ExprPtr e, Rational& out) {
    if (const auto* il = expr_cast<IntegerLit>(e)) {
        out = Rational(il->value);
        return true;
    }
    if (const auto* rl = expr_cast<RationalLit>(e)) {
        out = Rational(rl->numerator, rl->denominator);
        return true;
    }
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (extract_rational_const(un->operand, out)) {
            out = -out;
            return true;
        }
    }
    return false;
}

// Parse `c0 + c1·k` (any constant c0, c1) into (c0, c1).  Returns false if
// the expression depends on `k` in a non-linear way.
[[nodiscard]] static bool parse_affine_in_k(
    ExprPtr e, const Symbol& k, ExprPtr& c0_out, ExprPtr& c1_out,
    symbolic::CASContext& ctx) {
    auto coeffs_res = algebra::univariate_coefficients(e, k, ctx);
    if (coeffs_res.is_error()) return false;
    const auto& cs = coeffs_res.value();
    if (cs.empty() || cs.size() > 2U) return false;
    c0_out = cs[0];
    c1_out = cs.size() == 2U ? cs[1] : ctx.arena().make<IntegerLit>(BigInt(0));
    // c1 must be a non-zero constant (no k dependence already enforced by
    // univariate_coefficients).
    auto c1_simp = ctx.simplify(c1_out);
    if (c1_simp.is_error()) return false;
    Rational c1_rat;
    if (!extract_rational_const(c1_simp.value(), c1_rat)) return false;
    if (c1_rat.numerator().is_zero()) return false;
    c1_out = c1_simp.value();
    return true;
}

// Build the antidifference of A / (k + a)^m as
//   A · (−1)^(m−1) / (m − 1)! · ψ^(m−1)(k + a).
// All arguments are ExprPtrs so symbolic coefficients pass through.
[[nodiscard]] static ExprPtr polygamma_antidiff(
    ExprPtr A, ExprPtr k_plus_a, unsigned int m, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    // (−1)^(m−1) / (m − 1)!
    BigInt sign = (m % 2U == 1U) ? BigInt(1) : BigInt(-1);
    BigInt fact(1);
    for (unsigned int i = 2; i < m; ++i) fact *= BigInt(static_cast<long long>(i));
    Rational scale(sign, fact);

    ExprPtr scale_expr;
    if (scale.denominator() == BigInt(1)) {
        scale_expr = arena.make<IntegerLit>(scale.numerator());
    } else {
        scale_expr = arena.make<RationalLit>(scale.numerator(), scale.denominator());
    }
    ExprPtr coefficient = arena.make<Binary>(BinaryOp::Mul, A, scale_expr);

    // Build polygamma(m−1, k+a)  (m=1 → digamma).
    ExprPtr ant;
    if (m == 1U) {
        std::vector<ExprPtr> args{k_plus_a};
        ant = arena.make<FuncCall>(BuiltinOp::Digamma, std::move(args));
    } else {
        std::vector<ExprPtr> args{
            arena.make<IntegerLit>(BigInt(static_cast<long long>(m - 1))),
            k_plus_a,
        };
        ant = arena.make<FuncCall>(BuiltinOp::Polygamma, std::move(args));
    }
    return arena.make<Binary>(BinaryOp::Mul, coefficient, ant);
}

// Fast structural check: does `expr` contain a node of the form Div or
// Pow with `k` in the divisor / negative exponent?  If not, the term is a
// polynomial in `k` (or doesn't depend on `k` at all) and the polygamma
// path is guaranteed to return nullopt — short-circuit before paying the
// `together` / `apart` overhead.
[[nodiscard]] static bool has_rational_dependency(
    ExprPtr expr, const std::string& k_name) {
    if (!expr) return false;
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Div) {
            // Any subtree (including the divisor) may carry k.
            return has_rational_dependency(bin->left, k_name) ||
                   has_rational_dependency(bin->right, k_name) ||
                   true;  // a Div node itself counts as candidate.
        }
        if (bin->op == BinaryOp::Pow) {
            if (const auto* exp_lit = expr_cast<IntegerLit>(bin->right);
                exp_lit && exp_lit->value.is_negative()) {
                return true;
            }
        }
        return has_rational_dependency(bin->left, k_name) ||
               has_rational_dependency(bin->right, k_name);
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr t : sum->terms) {
            if (has_rational_dependency(t, k_name)) return true;
        }
        return false;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr t : prod->factors) {
            if (has_rational_dependency(t, k_name)) return true;
        }
        return false;
    }
    if (const auto* un = expr_cast<Unary>(expr)) {
        return has_rational_dependency(un->operand, k_name);
    }
    return false;
}

// Try to interpret `term` as A / (linear(k))^m and return the polygamma
// antidifference.  Returns nullopt if the shape does not match.
[[nodiscard]] static std::optional<ExprPtr> try_polygamma_antidiff(
    ExprPtr term, const Symbol& k, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    if (!has_rational_dependency(term, k.name)) return std::nullopt;

    // Normalise: bring to single fraction.
    auto together_res = algebra::together(term, ctx);
    if (together_res.is_error()) return std::nullopt;
    auto parts = algebra::apart_num_den(together_res.value(), ctx);
    if (parts.is_error()) return std::nullopt;
    ExprPtr N = parts.value().numerator;
    ExprPtr D = parts.value().denominator;

    auto N_simp = ctx.simplify(N);
    auto D_simp = ctx.simplify(D);
    if (N_simp.is_error() || D_simp.is_error()) return std::nullopt;
    N = N_simp.value();
    D = D_simp.value();

    // Numerator must be a constant in k.
    auto N_coeffs_res = algebra::univariate_coefficients(N, k, ctx);
    if (N_coeffs_res.is_error()) return std::nullopt;
    if (N_coeffs_res.value().size() != 1U) return std::nullopt;
    ExprPtr A_const = N_coeffs_res.value()[0];

    // Denominator must be (linear)^m for some m ≥ 1.
    ExprPtr base = D;
    unsigned int m = 1U;
    if (const auto* bin = expr_cast<Binary>(D); bin && bin->op == BinaryOp::Pow) {
        const auto* exp_lit = expr_cast<IntegerLit>(bin->right);
        if (!exp_lit) return std::nullopt;
        if (exp_lit->value <= BigInt(0)) return std::nullopt;
        if (exp_lit->value.bit_length() > 16U) return std::nullopt;
        m = static_cast<unsigned int>(exp_lit->value.to_u64());
        base = bin->left;
    }

    ExprPtr c0;
    ExprPtr c1;
    if (!parse_affine_in_k(base, k, c0, c1, ctx)) return std::nullopt;

    // Normalise base to (k + a):  base = c1·k + c0 = c1·(k + c0/c1).
    // Pull out c1^m from the denominator into the numerator coefficient.
    ExprPtr a = arena.make<Binary>(BinaryOp::Div, c0, c1);
    auto a_simp = ctx.simplify(a);
    if (a_simp.is_error()) return std::nullopt;
    a = a_simp.value();
    ExprPtr k_sym = arena.make<Symbol>(k);
    ExprPtr k_plus_a = arena.make<Binary>(BinaryOp::Add, k_sym, a);
    auto k_plus_a_simp = ctx.simplify(k_plus_a);
    if (k_plus_a_simp.is_error()) return std::nullopt;
    k_plus_a = k_plus_a_simp.value();

    // Effective coefficient: A_const / c1^m.
    ExprPtr c1_m = (m == 1U)
        ? c1
        : arena.make<Binary>(BinaryOp::Pow, c1,
            arena.make<IntegerLit>(BigInt(static_cast<long long>(m))));
    ExprPtr A_eff_raw = arena.make<Binary>(BinaryOp::Div, A_const, c1_m);
    auto A_eff = ctx.simplify(A_eff_raw);
    if (A_eff.is_error()) return std::nullopt;

    return polygamma_antidiff(A_eff.value(), k_plus_a, m, ctx);
}

// Definite-sum variant of try_polygamma_antidiff: builds S, then returns
// S(b+1) − S(a).
[[nodiscard]] static Result<ExprPtr> try_polygamma_definite(
    ExprPtr term, const Symbol& var, ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto antidiff_opt = try_polygamma_antidiff(term, var, ctx);
    if (!antidiff_opt.has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Polygamma summation: term is not of the form "
                       "A/(linear(k))^m",
        });
    }
    AstArena& arena = ctx.arena();
    ExprPtr S = antidiff_opt.value();
    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));
    auto S_upper = ctx.substitute(S, var, upper_plus_one);
    if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
    auto S_lower = ctx.substitute(S, var, lower);
    if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
    ExprPtr diff = arena.make<Binary>(BinaryOp::Sub, S_upper.value(), S_lower.value());
    return ctx.simplify(diff);
}

[[nodiscard]] Result<ExprPtr> symbolic_sum(
    ExprPtr term,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    if (is_positive_infinity(upper)) {
        if (const auto* l = expr_cast<IntegerLit>(lower); l != nullptr && l->value == BigInt(1)) {
            if (auto exponent = reciprocal_power_exponent(term, var); exponent.has_value()) {
                return zeta_even_value(*exponent, ctx);
            }
        }
        // Infinite Gosper-summable series:  if S(k) → 0 as k → ∞ we cannot
        // decide this without an asymptotic analysis, so fall through to
        // diagnostic Unimplemented below rather than risk silent error.
    } else {
        // Definite finite-bound sum — try Gosper first (hypergeometric class
        // including all polynomial / rational-with-rational-antidifference
        // summands), then fall back to the polygamma path which closes the
        // remaining rational atoms A/(linear(k))^m via the digamma/polygamma
        // identity ψ^(n)(k+1) − ψ^(n)(k) = (−1)^n n!/k^(n+1).
        auto gosper_res = try_gosper_definite(term, var, lower, upper, ctx);
        if (gosper_res.is_ok()) return gosper_res;
        auto poly_res = try_polygamma_definite(term, var, lower, upper, ctx);
        if (poly_res.is_ok()) return poly_res;
        // Neither path succeeded: fall through to diagnostic.
    }

    return make_unimplemented<ExprPtr>(
        "calculus", "sum_closed_form",
        "general summand not in closed-form table",
        cas::error::reason_codes::SUMMATION_GENERAL,
        "Implement Petkovšek-WZ / Zeilberger creative telescoping or Abramov "
        "rational summation for terms outside the Gosper hypergeometric class",
        "F0.8");
}

Result<ExprPtr> sum(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    return symbolic_sum(expr, var, lower, upper, ctx);
}

} // namespace cas::calculus
