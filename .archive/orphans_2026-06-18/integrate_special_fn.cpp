#include "integrate_engine.hpp"
#include "cas/algebra.hpp"
#include <optional>
#include <vector>

namespace cas::calculus::integrate_detail {

namespace {

/**
 * @brief Recursively flattens an integrand expression into its product factors,
 * representing division X / Y as a product of X and Y^-1.
 * @param expr The expression to flatten.
 * @param factors Out-parameter to collect the flattened factors.
 * @param arena The AST arena to allocate the inverse power nodes.
 */
void flatten_factors(ExprPtr expr, std::vector<ExprPtr>& factors, AstArena& arena) {
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) {
            flatten_factors(f, factors, arena);
        }
    } else if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Mul) {
            flatten_factors(bin->left, factors, arena);
            flatten_factors(bin->right, factors, arena);
        } else if (bin->op == BinaryOp::Div) {
            flatten_factors(bin->left, factors, arena);
            // Represent right as right^-1
            ExprPtr inv_right = arena.make<Binary>(BinaryOp::Pow, bin->right, arena.make<IntegerLit>(BigInt(-1)));
            flatten_factors(inv_right, factors, arena);
        } else {
            factors.push_back(expr);
        }
    } else {
        factors.push_back(expr);
    }
}

/**
 * @brief Computes the integral of exp(a*x + b) * (c*x + d)^-k recursively.
 * @param a Coefficient of the exponent linear term.
 * @param b Constant of the exponent linear term.
 * @param c Coefficient of the denominator linear term.
 * @param d Constant of the denominator linear term.
 * @param k Exponent of the denominator (k >= 1).
 * @param var The integration variable.
 * @param arena The AST arena.
 * @return The integrated expression or a CASError on failure.
 */
Result<ExprPtr> integrate_exp_linear_power(
    Rational a, Rational b, Rational c, Rational d, int k,
    const Symbol& var, AstArena& arena)
{
    if (a.numerator().is_zero() || c.numerator().is_zero()) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Invalid linear coefficients"));
    }
    if (k < 1) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Power must be >= 1"));
    }

    if (k == 1) {
        // ans = (1/c) * exp((b*c - a*d)/c) * Ei( (a/c)*(c*x + d) )
        ExprPtr var_expr = arena.make<Symbol>(var);
        ExprPtr arg = make_product(arena, {
            make_rational(arena, a / c),
            make_sum(arena, {make_product(arena, {make_rational(arena, c), var_expr}), make_rational(arena, d)})
        });
        ExprPtr ei = make_function(arena, "Ei", {arg});
        ExprPtr coeff = make_product(arena, {
            make_rational(arena, Rational(BigInt(1)) / c),
            make_function(arena, "exp", {make_rational(arena, (b * c - a * d) / c)})
        });
        return ok(make_product(arena, {coeff, ei}));
    } else {
        // IBP: ∫ e^{ax+b} (cx+d)^-k dx = e^{ax+b} (cx+d)^{1-k} / ((1-k)*c) - a/((1-k)*c) * ∫ e^{ax+b} (cx+d)^{1-k} dx
        ExprPtr var_expr = arena.make<Symbol>(var);
        ExprPtr term1 = make_product(arena, {
            make_rational(arena, Rational(BigInt(1), BigInt(1 - k)) / c),
            make_function(arena, "exp", {make_sum(arena, {make_product(arena, {make_rational(arena, a), var_expr}), make_rational(arena, b)})}),
            make_binary(arena, BinaryOp::Pow, make_sum(arena, {make_product(arena, {make_rational(arena, c), var_expr}), make_rational(arena, d)}), make_integer(arena, 1 - k))
        });
        auto sub_res = integrate_exp_linear_power(a, b, c, d, k - 1, var, arena);
        if (sub_res.is_error()) return sub_res;
        ExprPtr term2 = make_product(arena, {
            make_rational(arena, Rational(-a) / (Rational(BigInt(1 - k)) * c)),
            sub_res.value()
        });
        return ok(make_sum(arena, {term1, term2}));
    }
}

} // namespace

Result<ExprPtr> Integrator::try_integrate_special_function(ExprPtr expr, const Symbol& var) {
    // 1. Single Exp(quadratic) case: exp(a*x^2 + b*x + c)
    if (const auto* fc = expr_cast<FuncCall>(expr)) {
        if (fc->func_id == BuiltinOp::Exp && fc->args.size() == 1U) {
            ExprPtr arg = fc->args[0];
            if (auto quad = extract_quadratic_argument(arg, var)) {
                if (!quad->quadratic.numerator().is_zero()) {
                    Rational a = quad->quadratic;
                    Rational b = quad->linear;
                    Rational c = quad->constant;
                    ExprPtr var_expr = arena_.make<Symbol>(var);
                    
                    if (a > Rational(0)) {
                        // C * erfi(u)
                        // C = sqrt(pi) / (2 * sqrt(a)) * exp(c - b^2 / (4*a))
                        // u = sqrt(a) * (x + b / (2*a))
                        ExprPtr u = make_product(arena_, {
                            make_function(arena_, "sqrt", {make_rational(arena_, a)}),
                            make_sum(arena_, {var_expr, make_rational(arena_, b / (a * Rational(2)))})
                        });
                        ExprPtr C = make_product(arena_, {
                            make_binary(arena_, BinaryOp::Div,
                                make_function(arena_, "sqrt", {arena_.make<Constant>(MathConstant::Pi)}),
                                make_product(arena_, {make_integer(arena_, 2), make_function(arena_, "sqrt", {make_rational(arena_, a)})})
                            ),
                            make_function(arena_, "exp", {make_rational(arena_, c - b * b / (a * Rational(4)))})
                        });
                        return ok(make_product(arena_, {C, make_function(arena_, "erfi", {u})}));
                    } else {
                        // a < 0
                        // C * erf(u)
                        // C = sqrt(pi) / (2 * sqrt(-a)) * exp(c - b^2 / (4*a))
                        // u = sqrt(-a) * (x + b / (2*a))
                        Rational neg_a = -a;
                        ExprPtr u = make_product(arena_, {
                            make_function(arena_, "sqrt", {make_rational(arena_, neg_a)}),
                            make_sum(arena_, {var_expr, make_rational(arena_, b / (a * Rational(2)))})
                        });
                        ExprPtr C = make_product(arena_, {
                            make_binary(arena_, BinaryOp::Div,
                                make_function(arena_, "sqrt", {arena_.make<Constant>(MathConstant::Pi)}),
                                make_product(arena_, {make_integer(arena_, 2), make_function(arena_, "sqrt", {make_rational(arena_, neg_a)})})
                            ),
                            make_function(arena_, "exp", {make_rational(arena_, c - b * b / (a * Rational(4)))})
                        });
                        return ok(make_product(arena_, {C, make_function(arena_, "erf", {u})}));
                    }
                }
            }
        }
    }

    // 2. Product/Div cases
    std::vector<ExprPtr> factors;
    flatten_factors(expr, factors, arena_);
    
    ExprPtr exp_factor = nullptr;
    ExprPtr ln_factor = nullptr;
    ExprPtr trig_factor = nullptr;
    ExprPtr poly_denom = nullptr;
    int denom_power = 0;
    std::vector<ExprPtr> const_factors;
    bool match_failed = false;

    for (ExprPtr f : factors) {
        if (!depends_on(f, var)) {
            const_factors.push_back(f);
            continue;
        }
        
        // Check exp(affine)
        if (const auto* fc = expr_cast<FuncCall>(f)) {
            if (fc->func_id == BuiltinOp::Exp && fc->args.size() == 1U) {
                if (auto aff = extract_affine_argument(fc->args[0], var)) {
                    if (!aff->coefficient.numerator().is_zero()) {
                        if (!exp_factor) {
                            exp_factor = fc->args[0];
                            continue;
                        }
                    }
                }
            }
            if ((fc->func_id == BuiltinOp::Ln || fc->func_id == BuiltinOp::Log) && fc->args.size() == 1U) {
                if (auto aff = extract_affine_argument(fc->args[0], var)) {
                    if (!aff->coefficient.numerator().is_zero()) {
                        if (!ln_factor) {
                            ln_factor = fc->args[0];
                            continue;
                        }
                    }
                }
            }
            if ((fc->func_id == BuiltinOp::Sin || fc->func_id == BuiltinOp::Cos) && fc->args.size() == 1U) {
                if (auto aff = extract_affine_argument(fc->args[0], var)) {
                    if (!aff->coefficient.numerator().is_zero()) {
                        if (!trig_factor) {
                            trig_factor = f;
                            continue;
                        }
                    }
                }
            }
        }
        
        // Check B^-k
        if (const auto* bin = expr_cast<Binary>(f)) {
            if (bin->op == BinaryOp::Pow) {
                if (const auto* il = expr_cast<IntegerLit>(bin->right)) {
                    if (il->value.is_negative()) {
                        if (auto aff = extract_affine_argument(bin->left, var)) {
                            if (!aff->coefficient.numerator().is_zero()) {
                                if (!poly_denom) {
                                    poly_denom = bin->left;
                                    denom_power = static_cast<int>((-il->value).to_u64());
                                    continue;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        match_failed = true;
        break;
    }

    if (match_failed) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Not a recognized special function product"));
    }

    // Now check which specific combination matched!
    // A) exp(A) / B^k
    if (exp_factor && poly_denom && denom_power >= 1 && !ln_factor && !trig_factor) {
        auto aff_e = extract_affine_argument(exp_factor, var);
        auto aff_d = extract_affine_argument(poly_denom, var);
        if (aff_e && aff_d) {
            auto prim_res = integrate_exp_linear_power(
                aff_e->coefficient, aff_e->constant,
                aff_d->coefficient, aff_d->constant,
                denom_power, var, arena_
            );
            if (prim_res.is_ok()) {
                if (const_factors.empty()) return prim_res;
                const_factors.push_back(prim_res.value());
                return ok(make_product(arena_, std::move(const_factors)));
            }
        }
    }

    // B) sin/cos(W) / B (power 1 only)
    if (trig_factor && poly_denom && denom_power == 1 && !exp_factor && !ln_factor) {
        const auto* fc = expr_cast<FuncCall>(trig_factor);
        auto aff_t = extract_affine_argument(fc->args[0], var);
        auto aff_d = extract_affine_argument(poly_denom, var);
        if (aff_t && aff_d) {
            Rational a = aff_t->coefficient;
            Rational b = aff_t->constant;
            Rational c = aff_d->coefficient;
            Rational d = aff_d->constant;
            Rational theta = b - (a * d) / c;
            
            ExprPtr var_expr = arena_.make<Symbol>(var);
            ExprPtr arg = make_product(arena_, {
                make_rational(arena_, a / c),
                make_sum(arena_, {make_product(arena_, {make_rational(arena_, c), var_expr}), make_rational(arena_, d)})
            });
            ExprPtr si = make_function(arena_, "Si", {arg});
            ExprPtr ci = make_function(arena_, "Ci", {arg});
            ExprPtr cos_theta = make_function(arena_, "cos", {make_rational(arena_, theta)});
            ExprPtr sin_theta = make_function(arena_, "sin", {make_rational(arena_, theta)});
            
            ExprPtr prim;
            if (fc->func_id == BuiltinOp::Sin) {
                ExprPtr term1 = make_product(arena_, {make_rational(arena_, Rational(BigInt(1)) / c), cos_theta, si});
                ExprPtr term2 = make_product(arena_, {make_rational(arena_, Rational(BigInt(1)) / c), sin_theta, ci});
                prim = make_sum(arena_, {term1, term2});
            } else {
                ExprPtr term1 = make_product(arena_, {make_rational(arena_, Rational(BigInt(1)) / c), cos_theta, ci});
                ExprPtr term2 = make_unary(arena_, UnaryOp::Neg,
                    make_product(arena_, {make_rational(arena_, Rational(BigInt(1)) / c), sin_theta, si})
                );
                prim = make_sum(arena_, {term1, term2});
            }
            
            if (const_factors.empty()) return ok(prim);
            const_factors.push_back(prim);
            return ok(make_product(arena_, std::move(const_factors)));
        }
    }

    // C) exp(A) * ln(B)
    if (exp_factor && ln_factor && !poly_denom && !trig_factor) {
        auto aff_e = extract_affine_argument(exp_factor, var);
        auto aff_l = extract_affine_argument(ln_factor, var);
        if (aff_e && aff_l) {
            Rational a = aff_e->coefficient;
            Rational b = aff_e->constant;
            Rational c = aff_l->coefficient;
            Rational d = aff_l->constant;
            
            ExprPtr sub_integral = nullptr;
            if (c == Rational(1) && d == Rational(0)) {
                ExprPtr var_expr = arena_.make<Symbol>(var);
                if (a == Rational(-1)) {
                    // sub = -gamma_incomplete(0, x)
                    sub_integral = make_unary(arena_, UnaryOp::Neg,
                        make_function(arena_, "gamma_incomplete", {make_integer(arena_, 0), var_expr})
                    );
                } else if (a == Rational(1)) {
                    // sub = -gamma_incomplete(0, -x)
                    sub_integral = make_unary(arena_, UnaryOp::Neg,
                        make_function(arena_, "gamma_incomplete", {make_integer(arena_, 0), make_unary(arena_, UnaryOp::Neg, var_expr)})
                    );
                }
            }
            
            if (!sub_integral) {
                // Fallback to Ei representation
                auto sub_res = integrate_exp_linear_power(a, b, c, d, 1, var, arena_);
                if (sub_res.is_error()) return sub_res;
                sub_integral = sub_res.value();
            }
            
            ExprPtr var_expr = arena_.make<Symbol>(var);
            ExprPtr term1 = make_product(arena_, {
                make_rational(arena_, Rational(1) / a),
                make_function(arena_, "exp", {make_sum(arena_, {make_product(arena_, {make_rational(arena_, a), var_expr}), make_rational(arena_, b)})}),
                make_function(arena_, "ln", {make_sum(arena_, {make_product(arena_, {make_rational(arena_, c), var_expr}), make_rational(arena_, d)})})
            });
            ExprPtr term2 = make_product(arena_, {
                make_rational(arena_, -c / a),
                sub_integral
            });
            ExprPtr prim = make_sum(arena_, {term1, term2});
            
            if (const_factors.empty()) return ok(prim);
            const_factors.push_back(prim);
            return ok(make_product(arena_, std::move(const_factors)));
        }
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "No matched special function pattern"));
}

} // namespace cas::calculus::integrate_detail
