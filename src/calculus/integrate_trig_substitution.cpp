#include "integrate_engine.hpp"
#include "cas/algebra.hpp"

#include <optional>
#include <vector>

namespace cas::calculus::integrate_detail {
namespace {

struct ConstantSquareMinusVariableSquare {
    ExprPtr constant_base;
};

struct VariableSquareAndConstantSquare {
    ExprPtr constant_base;
};

[[nodiscard]] bool extract_square_base(ExprPtr expr, ExprPtr& base) {
    const auto* power = expr_cast<Binary>(expr);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return false;
    }
    if (!is_rational_value(power->right, 2, 1)) {
        return false;
    }
    base = power->left;
    return true;
}

[[nodiscard]] bool extract_negative_square_base(ExprPtr expr, ExprPtr& base) {
    if (const auto* unary = expr_cast<Unary>(expr)) {
        return unary->op == UnaryOp::Neg && extract_square_base(unary->operand, base);
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        if (product->factors.size() != 2U) {
            return false;
        }
        if (is_negative_one(product->factors[0]) && extract_square_base(product->factors[1], base)) {
            return true;
        }
        if (is_negative_one(product->factors[1]) && extract_square_base(product->factors[0], base)) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] bool extract_constant_square_base(AstArena& arena, ExprPtr expr, const Symbol& var, ExprPtr& base) {
    if (is_one(expr)) {
        base = make_integer(arena, 1);
        return true;
    }
    if (!extract_square_base(expr, base)) {
        return false;
    }
    return !depends_on(base, var);
}

[[nodiscard]] bool matches_negative_square_of_variable(ExprPtr expr, const Symbol& var) {
    ExprPtr base{};
    return extract_negative_square_base(expr, base) && is_same_symbol(base, var);
}

[[nodiscard]] std::optional<ConstantSquareMinusVariableSquare> match_constant_square_minus_variable_square(
    AstArena& arena,
    ExprPtr expr,
    const Symbol& var) {
    auto pair_matches = [&](ExprPtr constant_term, ExprPtr variable_term) -> std::optional<ConstantSquareMinusVariableSquare> {
        ExprPtr constant_base{};
        if (!extract_constant_square_base(arena, constant_term, var, constant_base)) {
            return std::nullopt;
        }
        if (!matches_negative_square_of_variable(variable_term, var)) {
            return std::nullopt;
        }
        return ConstantSquareMinusVariableSquare{constant_base};
    };

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Sub) {
            ExprPtr constant_base{};
            if (extract_constant_square_base(arena, binary->left, var, constant_base) &&
                matches_square_of_variable(binary->right, var)) {
                return ConstantSquareMinusVariableSquare{constant_base};
            }
        }
        if (binary->op == BinaryOp::Add) {
            if (auto matched = pair_matches(binary->left, binary->right); matched.has_value()) {
                return matched;
            }
            return pair_matches(binary->right, binary->left);
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() != 2U) {
            return std::nullopt;
        }
        if (auto matched = pair_matches(sum->terms[0], sum->terms[1]); matched.has_value()) {
            return matched;
        }
        return pair_matches(sum->terms[1], sum->terms[0]);
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<VariableSquareAndConstantSquare> match_variable_square_plus_constant_square(
    AstArena& arena,
    ExprPtr expr,
    const Symbol& var) {
    auto pair_matches = [&](ExprPtr variable_term, ExprPtr constant_term) -> std::optional<VariableSquareAndConstantSquare> {
        ExprPtr constant_base{};
        if (!matches_square_of_variable(variable_term, var) ||
            !extract_constant_square_base(arena, constant_term, var, constant_base)) {
            return std::nullopt;
        }
        return VariableSquareAndConstantSquare{constant_base};
    };

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Add) {
            if (auto matched = pair_matches(binary->left, binary->right); matched.has_value()) {
                return matched;
            }
            return pair_matches(binary->right, binary->left);
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() != 2U) {
            return std::nullopt;
        }
        if (auto matched = pair_matches(sum->terms[0], sum->terms[1]); matched.has_value()) {
            return matched;
        }
        return pair_matches(sum->terms[1], sum->terms[0]);
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<VariableSquareAndConstantSquare> match_variable_square_minus_constant_square(
    AstArena& arena,
    ExprPtr expr,
    const Symbol& var) {
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Sub) {
            ExprPtr constant_base{};
            if (matches_square_of_variable(binary->left, var) &&
                extract_constant_square_base(arena, binary->right, var, constant_base)) {
                return VariableSquareAndConstantSquare{constant_base};
            }
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() != 2U) {
            return std::nullopt;
        }
        for (std::size_t variable_index = 0; variable_index < 2U; ++variable_index) {
            const std::size_t constant_index = 1U - variable_index;
            ExprPtr constant_base{};
            if (matches_square_of_variable(sum->terms[variable_index], var) &&
                extract_negative_square_base(sum->terms[constant_index], constant_base) &&
                !depends_on(constant_base, var)) {
                return VariableSquareAndConstantSquare{constant_base};
            }
        }
    }

    return std::nullopt;
}

[[nodiscard]] ExprPtr square_expr(AstArena& arena, ExprPtr base) {
    return make_binary(arena, BinaryOp::Pow, base, make_integer(arena, 2));
}

[[nodiscard]] ExprPtr half_times(AstArena& arena, ExprPtr expr) {
    return make_product(arena, {make_rational(arena, 1, 2), expr});
}

}  // namespace

Result<ExprPtr> Integrator::integrate_sqrt_quadratic_impl(ExprPtr radicand, const Symbol& var) {
    ExprPtr x = arena_.make<Symbol>(var);
    ExprPtr sqrt_radicand = make_function(arena_, "sqrt", {radicand});
    ExprPtr x_sqrt = make_product(arena_, {x, sqrt_radicand});

    if (auto matched = match_constant_square_minus_variable_square(arena_, radicand, var); matched.has_value()) {
        ExprPtr a_squared = square_expr(arena_, matched->constant_base);
        ExprPtr arcsin_arg = make_binary(arena_, BinaryOp::Div, x, matched->constant_base);
        return ok(half_times(arena_, make_sum(arena_, {
            x_sqrt,
            make_product(arena_, {a_squared, make_function(arena_, "arcsin", {arcsin_arg})}),
        })));
    }

    if (auto matched = match_variable_square_plus_constant_square(arena_, radicand, var); matched.has_value()) {
        // ∫ √(x² + a²) dx = ½(x·√(x²+a²) + a²·asinh(x/a))
        // Emit asinh directly so the equivalence checker doesn't have to
        // bridge ln|x+√(x²+a²)| ↔ asinh(x/a) for every comparison.
        ExprPtr a_squared = square_expr(arena_, matched->constant_base);
        ExprPtr arg = make_binary(arena_, BinaryOp::Div, x, matched->constant_base);
        return ok(half_times(arena_, make_sum(arena_, {
            x_sqrt,
            make_product(arena_, {
                a_squared,
                make_function(arena_, "asinh", {arg}),
            }),
        })));
    }

    if (auto matched = match_variable_square_minus_constant_square(arena_, radicand, var); matched.has_value()) {
        // ∫ √(x² − a²) dx = ½(x·√(x²−a²) − a²·acosh(x/a))
        ExprPtr a_squared = square_expr(arena_, matched->constant_base);
        ExprPtr arg = make_binary(arena_, BinaryOp::Div, x, matched->constant_base);
        return ok(half_times(arena_, make_sum(arena_, {
            x_sqrt,
            make_unary(arena_, UnaryOp::Neg, make_product(arena_, {
                a_squared,
                make_function(arena_, "acosh", {arg}),
            })),
        })));
    }

    return fail<ExprPtr>(make_error(
        CASErrorKind::Unimplemented,
        "No supported trigonometric substitution pattern found for sqrt radicand"));
}

// HC-IBP-VDU: ∫ x² / √R(x) dx for R(x) ∈ {a² − x², a² + x², x² − a²}.
//   √(a² − x²) → (a²/2)·arcsin(x/a) − (x/2)·√(a² − x²)
//   √(a² + x²) → (x/2)·√(a² + x²) − (a²/2)·arcsinh(x/a)
//   √(x² − a²) → (x/2)·√(x² − a²) − (a²/2)·arccosh(x/a)
Result<ExprPtr> Integrator::integrate_xsq_over_sqrt_quadratic_impl(
    ExprPtr radicand, const Symbol& var) {
    ExprPtr x = arena_.make<Symbol>(var);
    ExprPtr sqrt_r = make_function(arena_, "sqrt", {radicand});
    ExprPtr half = make_rational(arena_, 1, 2);

    if (auto m = match_constant_square_minus_variable_square(arena_, radicand, var);
        m.has_value()) {
        ExprPtr a_sq = square_expr(arena_, m->constant_base);
        ExprPtr arc = make_function(arena_, "arcsin",
            {make_binary(arena_, BinaryOp::Div, x, m->constant_base)});
        return ok(make_sum(arena_, {
            make_product(arena_, {half, a_sq, arc}),
            make_unary(arena_, UnaryOp::Neg,
                make_product(arena_, {half, x, sqrt_r})),
        }));
    }
    if (auto m = match_variable_square_plus_constant_square(arena_, radicand, var);
        m.has_value()) {
        ExprPtr a_sq = square_expr(arena_, m->constant_base);
        ExprPtr arc = make_function(arena_, "asinh",
            {make_binary(arena_, BinaryOp::Div, x, m->constant_base)});
        return ok(make_sum(arena_, {
            make_product(arena_, {half, x, sqrt_r}),
            make_unary(arena_, UnaryOp::Neg,
                make_product(arena_, {half, a_sq, arc})),
        }));
    }
    if (auto m = match_variable_square_minus_constant_square(arena_, radicand, var);
        m.has_value()) {
        ExprPtr a_sq = square_expr(arena_, m->constant_base);
        ExprPtr arc = make_function(arena_, "acosh",
            {make_binary(arena_, BinaryOp::Div, x, m->constant_base)});
        return ok(make_sum(arena_, {
            make_product(arena_, {half, x, sqrt_r}),
            make_unary(arena_, UnaryOp::Neg,
                make_product(arena_, {half, a_sq, arc})),
        }));
    }
    return fail<ExprPtr>(make_error(
        CASErrorKind::Unimplemented,
        "No supported trig-substitution pattern for x² / sqrt(quadratic)"));
}

// ∫ dx / √(A x² + B x + C) for rational A≠0, B, C — general completing-the-square.
// Generalizes the specific a²−x² / x²±a² matchers (no closed pattern table).
//
//   Let R = A x² + B x + C = A·(x + h)² + k  with h = B/(2A), k = C − B²/(4A).
//   • A > 0:  ∫ = (1/√A)·ln| √A·(x+h) + √R |          (d/dx ⇒ 1/√R, verified)
//   • A < 0 (needs k > 0, else no real domain):
//             ∫ = (1/√|A|)·arcsin( √|A|·(x+h) / √k )
//   √A, √|A|, √k are built via sqrt() and stay exact for any rational argument.
Result<ExprPtr> Integrator::integrate_inverse_sqrt_quadratic_impl(ExprPtr radicand, const Symbol& var) {
    auto coeffs_res = algebra::univariate_coefficients(radicand, var, context_);
    if (coeffs_res.is_error())
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "1/sqrt(quadratic): radicand not a univariate polynomial in var"));
    const auto& coeffs = coeffs_res.value();              // ascending: [C, B, A]
    if (coeffs.size() != 3U)                              // degree exactly 2
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "1/sqrt(quadratic): radicand is not quadratic in var"));

    auto C = exact_scalar_from_expr(coeffs[0]);
    auto B = exact_scalar_from_expr(coeffs[1]);
    auto A = exact_scalar_from_expr(coeffs[2]);
    if (!A.has_value() || !B.has_value() || !C.has_value() || A->numerator().is_zero())
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "1/sqrt(quadratic): non-rational or degenerate coefficients"));

    const Rational two(BigInt(2)), four(BigInt(4));
    const Rational h = *B / (two * *A);                   // B/(2A)
    const Rational k = *C - (*B * *B) / (four * *A);      // C − B²/(4A)

    ExprPtr x = arena_.make<Symbol>(var);
    // u = x + h  (drop the +0 term when h == 0 to keep the form minimal)
    ExprPtr u = h.numerator().is_zero()
        ? x
        : make_sum(arena_, {x, make_rational(arena_, h)});
    ExprPtr sqrt_R = make_function(arena_, "sqrt", {radicand});

    if (A->numerator().is_negative()) {
        if (!(k.numerator().is_positive()))              // A<0 needs k>0 for real values
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                "1/sqrt(quadratic): A<0 with k≤0 has no real domain"));
        const Rational absA = -*A;
        ExprPtr sqrt_absA = make_function(arena_, "sqrt", {make_rational(arena_, absA)});
        ExprPtr sqrt_k = make_function(arena_, "sqrt", {make_rational(arena_, k)});
        ExprPtr arg = make_binary(arena_, BinaryOp::Div,
            make_product(arena_, {sqrt_absA, u}), sqrt_k);
        ExprPtr inv = make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), sqrt_absA);
        return ok(make_product(arena_, {inv, make_function(arena_, "arcsin", {arg})}));
    }

    // A > 0
    ExprPtr sqrt_A = make_function(arena_, "sqrt", {make_rational(arena_, *A)});
    ExprPtr ln_arg = make_sum(arena_, {make_product(arena_, {sqrt_A, u}), sqrt_R});
    ExprPtr inv = make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), sqrt_A);
    return ok(make_product(arena_, {inv,
        make_function(arena_, "ln", {make_function(arena_, "abs", {ln_arg})})}));
}

// ∫ xᵏ / √(c − d·x²) dx for integer k ≥ 0, rational c>0, d>0 — the radical-monomial
// family the asin/acos integration-by-parts chain reduces to. Driven by the reduction
//   I_k = [(k−1)·c·I_{k−2} − xᵏ⁻¹·√(c−dx²)] / (k·d)              (k ≥ 2),
// derived from d/dx[xᵏ⁻¹√(c−dx²)] = (k−1)c·xᵏ⁻²/√R − k·d·xᵏ/√R. Bases:
//   I₀ = ∫ dx/√(c−dx²)          (delegated to integrate_inverse_sqrt_quadratic)
//   I₁ = ∫ x dx/√(c−dx²) = −(1/d)·√(c−dx²).
// The descent steps by 2, so only the single base of matching parity is needed.
// No closed pattern table; k is the actual integer power, not a hardcoded set.
Result<ExprPtr> Integrator::integrate_monomial_over_sqrt_quadratic_impl(
    long long k, ExprPtr radicand, const Symbol& var) {
    if (k < 0)
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "x^k/sqrt(quadratic): negative power"));

    // Validate the radicand is c − d·x²  (ascending coeffs [C, B, A], B=0, A<0, C>0).
    auto coeffs_res = algebra::univariate_coefficients(radicand, var, context_);
    if (coeffs_res.is_error())
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "x^k/sqrt(quadratic): radicand not a univariate polynomial in var"));
    const auto& coeffs = coeffs_res.value();
    if (coeffs.size() != 3U)
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "x^k/sqrt(quadratic): radicand is not quadratic in var"));
    auto C = exact_scalar_from_expr(coeffs[0]);
    auto B = exact_scalar_from_expr(coeffs[1]);
    auto A = exact_scalar_from_expr(coeffs[2]);
    if (!A.has_value() || !B.has_value() || !C.has_value())
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "x^k/sqrt(quadratic): non-rational coefficients"));
    if (!B->numerator().is_zero() || !A->numerator().is_negative() ||
        !C->numerator().is_positive())
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "x^k/sqrt(quadratic): only c − d·x² (B=0, A<0, C>0) supported"));

    const Rational& c = *C;
    const Rational d = -*A;                       // d > 0
    const Rational one(BigInt(1));

    ExprPtr x = arena_.make<Symbol>(var);
    ExprPtr sqrt_R = make_function(arena_, "sqrt", {radicand});

    // Base of matching parity.
    ExprPtr cur;
    const long long start = k % 2;                // 0 (even) or 1 (odd)
    if (start == 0) {
        auto base = integrate_inverse_sqrt_quadratic(radicand, var);   // I₀
        if (base.is_error())
            return base;
        cur = base.value();
    } else {
        ExprPtr inv_d = make_rational(arena_, one / d);                // I₁ = −(1/d)√R
        cur = make_unary(arena_, UnaryOp::Neg,
            make_product(arena_, {inv_d, sqrt_R}));
    }

    // Climb by 2 up to k.
    for (long long m = start + 2; m <= k; m += 2) {
        const Rational coef_prev = Rational(BigInt(m - 1)) * c;        // (k−1)·c
        const Rational inv_md = one / (Rational(BigInt(m)) * d);       // 1/(k·d)
        ExprPtr xpow = (m - 1 == 1)
            ? x
            : make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, m - 1));
        ExprPtr reduced = make_product(arena_, {make_rational(arena_, coef_prev), cur});
        ExprPtr boundary = make_unary(arena_, UnaryOp::Neg,
            make_product(arena_, {xpow, sqrt_R}));
        cur = make_product(arena_, {
            make_rational(arena_, inv_md),
            make_sum(arena_, {reduced, boundary}),
        });
    }
    return ok(cur);
}

}  // namespace cas::calculus::integrate_detail
