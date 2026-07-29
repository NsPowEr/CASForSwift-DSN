// Differentiation rules for builtin functions of arity > 1 — split out of
// differentiate_rules.cpp (anti-monolith, 500-line limit). Pure relocation:
// the dispatch order and every rule body are unchanged; the single-argument
// path stays in differentiate_rules.cpp and tail-calls into here.
//
// Covers: Bessel J/Y/I/K and BesselZero, incomplete gamma (upper/lower),
// Hypergeometric2F1, and Meijer G (variable arity).

#include "differentiate_internal.hpp"
#include "calculus_internal.hpp"

#include "cas/meijerg.hpp"

#include <string>

namespace cas::calculus {

Result<ExprPtr> Differentiator::differentiate_function_multiarg(
    const FuncCall& call, const Symbol& var) {
    if (call.args.size() == 2U) {
        if (call.func_id == BuiltinOp::BesselZero) {
            return ok(make_integer(arena_, 0));
        }
        if (call.func_id == BuiltinOp::BesselJ || call.func_id == BuiltinOp::BesselY || call.func_id == BuiltinOp::BesselI || call.func_id == BuiltinOp::BesselK) {
            ExprPtr nu = call.args[0];
            ExprPtr x = call.args[1];

            if (depends_on(nu, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of Bessel functions w.r.t order is not implemented"));
            }

            auto x_derivative = differentiate_once(x, var);
            if (x_derivative.is_error()) return x_derivative;

            ExprPtr outer;
            if (call.func_id == BuiltinOp::BesselJ) {
                outer = make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_integer(arena_, 2)),
                    make_sum(arena_, {
                        make_function(arena_, "BesselJ", {make_sum(arena_, {nu, make_integer(arena_, -1)}), x}),
                        make_unary(arena_, UnaryOp::Neg, make_function(arena_, "BesselJ", {make_sum(arena_, {nu, make_integer(arena_, 1)}), x}))
                    })
                });
            } else if (call.func_id == BuiltinOp::BesselY) {
                outer = make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_integer(arena_, 2)),
                    make_sum(arena_, {
                        make_function(arena_, "BesselY", {make_sum(arena_, {nu, make_integer(arena_, -1)}), x}),
                        make_unary(arena_, UnaryOp::Neg, make_function(arena_, "BesselY", {make_sum(arena_, {nu, make_integer(arena_, 1)}), x}))
                    })
                });
            } else if (call.func_id == BuiltinOp::BesselI) {
                outer = make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_integer(arena_, 2)),
                    make_sum(arena_, {
                        make_function(arena_, "BesselI", {make_sum(arena_, {nu, make_integer(arena_, -1)}), x}),
                        make_function(arena_, "BesselI", {make_sum(arena_, {nu, make_integer(arena_, 1)}), x})
                    })
                });
            } else if (call.func_id == BuiltinOp::BesselK) {
                outer = make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_unary(arena_, UnaryOp::Neg, make_integer(arena_, 1)), make_integer(arena_, 2)),
                    make_sum(arena_, {
                        make_function(arena_, "BesselK", {make_sum(arena_, {nu, make_integer(arena_, -1)}), x}),
                        make_function(arena_, "BesselK", {make_sum(arena_, {nu, make_integer(arena_, 1)}), x})
                    })
                });
            }

            return ok(make_product(arena_, {outer, x_derivative.value()}));
        }
        if (call.func_id == BuiltinOp::GammaIncomplete
            || call.func_id == BuiltinOp::GammaIncompleteLower) {
            // A7 §5.9 (DLMF 8.8.13):  d/dz Γ(a,z) = −z^{a−1} e^{−z},
            //                         d/dz γ(a,z) = +z^{a−1} e^{−z}.
            // Differentiation w.r.t. the order a needs the incomplete-gamma
            // derivative w.r.t. its first parameter (a Meijer-G / G-function
            // series) — not implemented: structured error, never a wrong value.
            ExprPtr a = call.args[0];
            ExprPtr z = call.args[1];
            if (depends_on(a, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of incomplete gamma w.r.t. order a not implemented"));
            }
            auto z_d = differentiate_once(z, var);
            if (z_d.is_error()) return z_d;
            ExprPtr kernel = make_product(arena_, {
                make_power(arena_, z, make_sum(arena_, {a, make_integer(arena_, -1)})),
                make_function(arena_, "exp",
                    {make_unary(arena_, UnaryOp::Neg, z)})});
            ExprPtr outer = call.func_id == BuiltinOp::GammaIncomplete
                ? make_unary(arena_, UnaryOp::Neg, kernel)
                : kernel;
            return ok(make_product(arena_, {outer, z_d.value()}));
        }
        if (call.func_id == BuiltinOp::Hypergeometric0F1) {
            ExprPtr b = call.args[0];
            ExprPtr z = call.args[1];
            if (depends_on(b, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of Hypergeometric0F1 w.r.t parametri non supportata"));
            }
            auto z_d = differentiate_once(z, var);
            if (z_d.is_error()) return z_d;
            ExprPtr b_plus_1 = make_sum(arena_, {b, make_integer(arena_, 1)});
            ExprPtr new_fc = make_function(arena_, "Hypergeometric0F1", {b_plus_1, z});
            ExprPtr div_expr = arena_.make<Binary>(BinaryOp::Div, make_integer(arena_, 1), b);
            ExprPtr result = make_product(arena_, {div_expr, new_fc, z_d.value()});
            return ok(result);
        }
    } else if (call.args.size() == 3U) {
        if (call.func_id == BuiltinOp::Hypergeometric1F1) {
            ExprPtr a = call.args[0];
            ExprPtr b = call.args[1];
            ExprPtr z = call.args[2];
            if (depends_on(a, var) || depends_on(b, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of Hypergeometric1F1 w.r.t parametri non supportata"));
            }
            auto z_d = differentiate_once(z, var);
            if (z_d.is_error()) return z_d;
            ExprPtr a_plus_1 = make_sum(arena_, {a, make_integer(arena_, 1)});
            ExprPtr b_plus_1 = make_sum(arena_, {b, make_integer(arena_, 1)});
            ExprPtr new_fc = make_function(arena_, "Hypergeometric1F1",
                {a_plus_1, b_plus_1, z});
            ExprPtr a_over_b = arena_.make<Binary>(BinaryOp::Div, a, b);
            return ok(make_product(arena_, {a_over_b, new_fc, z_d.value()}));
        }
    } else if (call.args.size() == 4U) {
        if (call.func_id == BuiltinOp::Hypergeometric2F1) {
            ExprPtr a = call.args[0];
            ExprPtr b = call.args[1];
            ExprPtr c = call.args[2];
            ExprPtr z = call.args[3];
            if (depends_on(a, var) || depends_on(b, var) || depends_on(c, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of Hypergeometric2F1 w.r.t parametri non supportata"));
            }
            auto z_d = differentiate_once(z, var);
            if (z_d.is_error()) return z_d;
            ExprPtr a1 = make_sum(arena_, {a, make_integer(arena_, 1)});
            ExprPtr b1 = make_sum(arena_, {b, make_integer(arena_, 1)});
            ExprPtr c1 = make_sum(arena_, {c, make_integer(arena_, 1)});
            ExprPtr new_fc = make_function(arena_, "Hypergeometric2F1",
                {a1, b1, c1, z});
            ExprPtr ab = make_product(arena_, {a, b});
            ExprPtr ab_over_c = arena_.make<Binary>(BinaryOp::Div, ab, c);
            return ok(make_product(arena_, {ab_over_c, new_fc, z_d.value()}));
        }
    }

    // Meijer G (variable arity): §6.5 h=+1 theta-shift (DLMF 16.19.5 family,
    // numerically certified — see meijerg_derivative_shift):
    //   d/dx G(u|a;b) = (u'/u) * G^{m,n+1}_{p+1,q+1}(u | 0,a ; b, +1).
    if (call.func_id == BuiltinOp::MeijerG) {
        auto view = symbolic::view_meijerg(call);
        if (view.is_error()) return fail<ExprPtr>(view.error());
        for (ExprPtr p : view.value().a) {
            if (depends_on(p, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of MeijerG w.r.t. parametri non supportata"));
            }
        }
        for (ExprPtr p : view.value().b) {
            if (depends_on(p, var)) {
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Differentiation of MeijerG w.r.t. parametri non supportata"));
            }
        }
        ExprPtr z = view.value().z;
        auto z_d = differentiate_once(z, var);
        if (z_d.is_error()) return z_d;
        auto shifted = symbolic::meijerg_derivative_shift(context_, call);
        if (shifted.is_error()) return shifted;
        ExprPtr chain = arena_.make<Binary>(BinaryOp::Div, z_d.value(), z);
        return ok(make_product(arena_, {chain, shifted.value()}));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation is not implemented for function '" + call.name + "' with " + std::to_string(call.args.size()) + " arguments"));
}

}  // namespace cas::calculus
