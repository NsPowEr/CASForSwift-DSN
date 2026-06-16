// F7.5 — ∫ 1/(x² ± a²)^n dx (integer n ≥ 1) via the Apostol/Bronstein
// closed-form recursion:
//
//   For (x² + a²)^n :   I_n = x / (2(n−1)·a²·(x²+a²)^(n−1))
//                            + (2n−3)/(2(n−1)·a²) · I_{n−1},
//                       I_1 = (1/a)·arctan(x/a).
//
//   For (x² − a²)^n :   I_n = − x / (2(n−1)·a²·(x²−a²)^(n−1))
//                            − (2n−3)/(2(n−1)·a²) · I_{n−1},
//                       I_1 = (1/(2a))·ln|(x−a)/(x+a)|.
//
// Both are universal in n; the matcher only verifies `radicand = x² ± a²`.
// Routes around the (currently buggy) Lazard-Rioboo-Trager path for
// `c/(x²+a²)^n` (see open task #5 / #6 — LRT garbage on irreducible
// quadratic denominators) while remaining mathematically sound on every
// integer multiplicity.

#include "integrate_engine.hpp"
#include "calculus_internal.hpp"
#include "cas/symbolic.hpp"

#include <cstdint>

namespace cas::calculus::integrate_detail {

namespace {

[[nodiscard]] ExprPtr make_integer_local(AstArena& arena, std::int64_t v) {
    return arena.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] ExprPtr make_pow(AstArena& arena, ExprPtr base, ExprPtr exp) {
    return arena.make<Binary>(BinaryOp::Pow, base, exp);
}

[[nodiscard]] ExprPtr make_div(AstArena& arena, ExprPtr num, ExprPtr den) {
    return arena.make<Binary>(BinaryOp::Div, num, den);
}

[[nodiscard]] ExprPtr make_neg(AstArena& arena, ExprPtr v) {
    return arena.make<Unary>(UnaryOp::Neg, v);
}

[[nodiscard]] ExprPtr make_sum_l(AstArena& arena, std::vector<ExprPtr> ts) {
    return arena.make<Sum>(std::move(ts));
}

[[nodiscard]] ExprPtr make_product_l(AstArena& arena, std::vector<ExprPtr> fs) {
    return arena.make<Product>(std::move(fs));
}

[[nodiscard]] ExprPtr make_func(AstArena& arena, const std::string& name,
                                std::vector<ExprPtr> args) {
    return arena.make<FuncCall>(name, std::move(args));
}

[[nodiscard]] ExprPtr make_rational_l(AstArena& arena, std::int64_t n,
                                      std::int64_t d) {
    return arena.make<RationalLit>(BigInt(n), BigInt(d));
}

}  // namespace

Result<ExprPtr> Integrator::integrate_inverse_quadratic_power(
    ExprPtr radicand, const BigInt& negative_exponent, const Symbol& var) {
    if (!negative_exponent.is_negative()) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "integrate_inverse_quadratic_power: non-negative exponent"));
    }
    const BigInt neg = -negative_exponent;
    if (neg.is_zero() || !(neg <= BigInt(64))) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "integrate_inverse_quadratic_power: exponent out of range"));
    }
    // Decompose radicand as x² ± C where C is constant in `var`. C represents
    // a² in the standard formulation; the recursion runs entirely in terms
    // of a² and never needs to extract a itself, except to build I_1.
    AstArena& a = arena_;
    auto match_plus_constant = [&](ExprPtr e, bool& is_plus, ExprPtr& a_sq) -> bool {
        if (const auto* sum = expr_cast<Sum>(e); sum && sum->terms.size() == 2U) {
            for (int i = 0; i < 2; ++i) {
                if (!matches_square_of_variable(sum->terms[i], var)) continue;
                ExprPtr other = sum->terms[1 - i];
                if (depends_on(other, var)) return false;
                if (const auto* un = expr_cast<Unary>(other); un && un->op == UnaryOp::Neg) {
                    is_plus = false; a_sq = un->operand; return true;
                }
                is_plus = true; a_sq = other; return true;
            }
        }
        if (const auto* bn = expr_cast<Binary>(e); bn && (bn->op == BinaryOp::Add || bn->op == BinaryOp::Sub)) {
            if (matches_square_of_variable(bn->left, var) && !depends_on(bn->right, var)) {
                is_plus = (bn->op == BinaryOp::Add); a_sq = bn->right; return true;
            }
            if (matches_square_of_variable(bn->right, var) && !depends_on(bn->left, var) && bn->op == BinaryOp::Add) {
                is_plus = true; a_sq = bn->left; return true;
            }
        }
        return false;
    };
    bool plus = false;
    ExprPtr a_sq;
    if (!match_plus_constant(radicand, plus, a_sq)) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "integrate_inverse_quadratic_power: radicand not x²±a²"));
    }
    const std::int64_t n = neg.to_u64();
    ExprPtr x = a.make<Symbol>(var);
    ExprPtr cb = make_func(a, "sqrt", std::vector<ExprPtr>{a_sq});

    ExprPtr I_curr;
    if (plus) {
        I_curr = make_div(a,
            make_func(a, "arctan",
                std::vector<ExprPtr>{make_div(a, x, cb)}),
            cb);
    } else {
        ExprPtr num_log = make_div(a,
            make_sum_l(a, {x, make_neg(a, cb)}),
            make_sum_l(a, {x, cb}));
        I_curr = make_product_l(a, {
            make_div(a, make_rational_l(a, 1, 2), cb),
            make_func(a, "ln",
                std::vector<ExprPtr>{make_func(a, "abs",
                    std::vector<ExprPtr>{num_log})})});
    }

    for (std::int64_t k = 2; k <= n; ++k) {
        const std::int64_t denom_factor = 2 * (k - 1);
        ExprPtr scale = make_product_l(a, {
            make_integer_local(a, denom_factor),
            a_sq});
        ExprPtr pow_km1 = make_pow(a, radicand,
            make_integer_local(a, k - 1));
        ExprPtr term_rational = make_div(a, x,
            make_product_l(a, {scale, pow_km1}));
        ExprPtr coeff_rec = make_div(a,
            make_integer_local(a, 2 * k - 3), scale);
        ExprPtr term_rec = make_product_l(a, {coeff_rec, I_curr});
        if (plus) {
            I_curr = make_sum_l(a, {term_rational, term_rec});
        } else {
            // (x² − a²) reduction: both terms carry a minus sign.
            I_curr = make_sum_l(a, {
                make_neg(a, term_rational),
                make_neg(a, term_rec)});
        }
    }

    auto simp = context_.simplify(I_curr);
    return simp.is_ok() ? ok(simp.value()) : ok(I_curr);
}

}  // namespace cas::calculus::integrate_detail
