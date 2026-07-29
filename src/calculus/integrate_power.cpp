#include "integrate_engine.hpp"
#include "cas/algebra.hpp"

#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

// integrate_power — extracted from integrate_product_power.cpp (anti-monolith, 2026-06-20).
// Handles ∫ base^exponent: power rule, reciprocal-quadratic-radical (arcsin/ln),
// 1/(quadratic)^n, ln^n recurrence, trig-square shortcuts. Shared helpers via integrate_engine.hpp.
Result<ExprPtr> Integrator::integrate_power_impl(const Binary& power, const Symbol& var) {
    if (is_same_symbol(power.right, var) && !depends_on(power.left, var)) {
        if (expr_is<DecimalLit>(power.left)) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Decimal literals are not supported in symbolic integration"));
        }
        if (is_one(power.left)) {
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "The integral of 1^x is undefined in this symbolic form"));
        }
        return ok(make_binary(arena_, BinaryOp::Div, make_binary(arena_, BinaryOp::Pow, power.left, arena_.make<Symbol>(var)), make_function(arena_, "ln", {power.left})));
    }

    if (is_rational_value(power.right, -1, 1)) {
        if (const auto* sqrt_call = expr_cast<FuncCall>(power.left);
            sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U) {
            if (matches_one_minus_square(sqrt_call->args.front(), var)) {
                return ok(make_function(arena_, "arcsin", {arena_.make<Symbol>(var)}));
            }
            ExprPtr constant_base{};
            if (matches_square_plus_constant_square(sqrt_call->args.front(), var, constant_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {x, power.left})})}));
            }
            // General fallback: ∫1/√(Ax²+Bx+C) via completing the square (Pow(Sqrt,−1) form).
            if (auto g = integrate_inverse_sqrt_quadratic(sqrt_call->args.front(), var); g.is_ok())
                return g;
        }

        if (matches_one_plus_square(power.left, var)) {
            return ok(make_function(arena_, "arctan", {arena_.make<Symbol>(var)}));
        }

        if (auto affine = extract_affine_argument(power.left, var);
            affine.has_value() && affine->coefficient.numerator() != BigInt(0)) {
            ExprPtr primitive = make_function(arena_, "ln", {make_function(arena_, "abs", {power.left})});
            if (affine->coefficient == Rational(BigInt(1))) {
                return ok(primitive);
            }
            return ok(make_product(arena_, {make_rational(arena_, Rational(BigInt(1)) / affine->coefficient), primitive}));
        }

        ExprPtr constant_base{};
        if (matches_square_minus_constant_square(power.left, var, constant_base)) {
            ExprPtr x = arena_.make<Symbol>(var);
            return ok(make_product(arena_, {
                make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_product(arena_, {make_integer(arena_, 2), constant_base})),
                make_function(arena_, "ln", {make_function(arena_, "abs", {make_binary(arena_, BinaryOp::Div, make_sum(arena_, {x, make_unary(arena_, UnaryOp::Neg, constant_base)}), make_sum(arena_, {x, constant_base}))})}),
            }));
        }

        if (auto r = integrate_via_partial_fractions(make_binary(arena_, BinaryOp::Pow, power.left, power.right), var); r.is_ok()) return r;
    }
    // F7.5: ∫ c/(x²±a²)^n via Apostol/Bronstein recursion (helper file).
    if (const auto* ie = expr_cast<IntegerLit>(power.right); ie && ie->value.is_negative())
        if (auto r = integrate_inverse_quadratic_power(power.left, ie->value, var); r.is_ok()) return r;

    // A46: una potenza NEGATIVA di un polinomio e' una funzione razionale, e la
    // via generale per le razionali e' Hermite + Lazard-Rioboo-Trager. Senza
    // questo passo la stessa integranda scritta `Pow(Q,−1)` (la forma che
    // producono partial_fractions e il simplifier) saltava la via razionale —
    // che il ramo `Div` imbocca subito — e finiva nelle strategie generiche a
    // valle. Misurato su 1/(x⁴+1): 40s come `Pow` contro 105ms come `Div`, e la
    // seconda frazione parziale di 1/(x⁶+1) sbagliava per 81s prima di cedere.
    // La forma sintattica non deve decidere ne' il costo ne' l'esito.
    if (const auto* ie = expr_cast<IntegerLit>(power.right); ie && ie->value.is_negative()) {
        if (auto base_degree = algebra::polynomial_degree(power.left, var, context_);
            base_degree.is_ok() && base_degree.value() > 0U) {
            ExprPtr positive_power = make_binary(arena_, BinaryOp::Pow, power.left,
                arena_.make<IntegerLit>(-ie->value));
            if (auto r = integrate_rational(
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), positive_power), var);
                r.is_ok()) {
                return r;
            }
        }
    }

    if (is_rational_value(power.right, -1, 2)) {
        if (matches_one_minus_square(power.left, var)) {
            return ok(make_function(arena_, "arcsin", {arena_.make<Symbol>(var)}));
        }
        ExprPtr constant_base{};
        if (matches_square_plus_constant_square(power.left, var, constant_base)) {
            ExprPtr x = arena_.make<Symbol>(var);
            return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {x, make_function(arena_, "sqrt", {power.left})})})}));
        }
        if (matches_square_minus_constant_square(power.left, var, constant_base)) {  // F7.5: ∫1/√(x²−a²) = ln|x+√…|
            ExprPtr x = arena_.make<Symbol>(var);
            return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {x, make_function(arena_, "sqrt", {power.left})})})})); }
        // General fallback: ∫1/√(Ax²+Bx+C) via completing the square (Pow(R,−1/2) form).
        if (auto g = integrate_inverse_sqrt_quadratic(power.left, var); g.is_ok()) return g;
    }

    if (const auto* integer = expr_cast<IntegerLit>(power.right)) {
        if (auto affine = extract_affine_argument(power.left, var);
            affine.has_value() && affine->coefficient.numerator() != BigInt(0)) {
            const BigInt next = integer->value + BigInt(1);
            if (!next.is_zero()) {
                return ok(make_product(arena_, {make_rational(arena_, Rational(BigInt(1), next) / affine->coefficient), make_binary(arena_, BinaryOp::Pow, power.left, arena_.make<IntegerLit>(next))}));
            }
        }
    }

    if (const auto* call = expr_cast<FuncCall>(power.left); call != nullptr && expr_is<IntegerLit>(power.right)) {
        const auto& exponent = expr_ref<IntegerLit>(power.right);
        BuiltinOp func_id = call->func_id;
        if (exponent.value == BigInt(2) && call->args.size() == 1U && is_same_symbol(call->args.front(), var)) {
            if (func_id == BuiltinOp::Sec) {
                return ok(make_function(arena_, "tan", {arena_.make<Symbol>(var)}));
            }
            if (func_id == BuiltinOp::Csc) {
                return ok(make_unary(arena_, UnaryOp::Neg, make_function(arena_, "cot", {arena_.make<Symbol>(var)})));
            }
            if (func_id == BuiltinOp::Tan) return ok(make_sum(arena_, {make_function(arena_, "tan", {arena_.make<Symbol>(var)}), make_unary(arena_, UnaryOp::Neg, arena_.make<Symbol>(var))}));
            if (func_id == BuiltinOp::Cot) return ok(make_sum(arena_, {make_unary(arena_, UnaryOp::Neg, make_function(arena_, "cot", {arena_.make<Symbol>(var)})), make_unary(arena_, UnaryOp::Neg, arena_.make<Symbol>(var))}));
        }
    }

    if (const auto* call = expr_cast<FuncCall>(power.left); call != nullptr && (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) && call->args.size() == 1U && is_same_symbol(call->args.front(), var) && expr_is<IntegerLit>(power.right)) {
        const auto& exponent = expr_ref<IntegerLit>(power.right);
        if (exponent.value > BigInt(0)) {
            // int ln(x)^n dx = x*ln(x)^n - n * int ln(x)^(n-1) dx
            ExprPtr x = arena_.make<Symbol>(var);
            ExprPtr term1 = make_product(arena_, {x, arena_.make<Binary>(BinaryOp::Pow, power.left, power.right)});
            
            ExprPtr next_exp = arena_.make<IntegerLit>(exponent.value - BigInt(1));
            ExprPtr next_pow = (exponent.value == BigInt(1)) ? make_integer(arena_, 1) : arena_.make<Binary>(BinaryOp::Pow, power.left, next_exp);
            
            auto next_integral = integrate_once(next_pow, var);
            if (next_integral.is_ok()) {
                ExprPtr term2 = make_product(arena_, {arena_.make<IntegerLit>(exponent.value), next_integral.value()});
                return ok(make_sum(arena_, {term1, make_unary(arena_, UnaryOp::Neg, term2)}));
            }
        }
    }

    if (!is_same_symbol(power.left, var) || depends_on(power.right, var)) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Only powers of the integration variable with constant exponent are implemented"));
    }

    if (is_negative_one(power.right)) {
        return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {arena_.make<Symbol>(var)})}));
    }

    if (const auto* integer = expr_cast<IntegerLit>(power.right)) {
        const BigInt next = integer->value + BigInt(1);
        return ok(make_product(arena_, {arena_.make<RationalLit>(BigInt(1), next), make_binary(arena_, BinaryOp::Pow, arena_.make<Symbol>(var), arena_.make<IntegerLit>(next))}));
    }

    ExprPtr exponent_plus_one = make_sum(arena_, {power.right, make_integer(arena_, 1)});
    return ok(make_binary(arena_, BinaryOp::Div, make_binary(arena_, BinaryOp::Pow, arena_.make<Symbol>(var), exponent_plus_one), exponent_plus_one));
}

}  // namespace cas::calculus::integrate_detail
