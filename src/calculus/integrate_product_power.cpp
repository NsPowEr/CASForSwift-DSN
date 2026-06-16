#include "integrate_engine.hpp"
#include "cas/algebra.hpp"

#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

namespace {

// Returns true if `expr` is exp(arg) where `arg` has negative or symbolic
// powers of `var` (i.e. `arg` is not a polynomial in `var`).  Such products
// cannot be closed by integration-by-parts in finite steps; the Risch
// algorithm must handle them.  This guard prevents the growing-denominator
// IBP loop that caused BUG-HANG-001.
bool has_exp_rational_non_poly_factor(const std::vector<ExprPtr>& factors, const Symbol& var) {
    for (ExprPtr f : factors) {
        const auto* call = expr_cast<FuncCall>(f);
        if (!call || call->func_id != BuiltinOp::Exp || call->args.size() != 1U) continue;
        ExprPtr arg = call->args[0];
        // Check: does the arg depend on var and contain any negative or
        // fractional power of var (e.g. 1/x, 1/x², √x)?
        // Walk the arg; if we find x^n with n < 0 or non-integer, or Div by var,
        // declare it non-polynomial.
        struct NonPolyWalker {
            const Symbol& v;
            bool found{false};
            void walk(ExprPtr e) {
                if (!e || found) return;
                if (const auto* bin = expr_cast<Binary>(e)) {
                    if (bin->op == BinaryOp::Div) {
                        // Check if denominator depends on v.
                        if (depends_on(bin->right, v)) { found = true; return; }
                    }
                    if (bin->op == BinaryOp::Pow && depends_on(bin->left, v)) {
                        // Negative or non-integer exponent → not polynomial.
                        if (const auto* el = expr_cast<IntegerLit>(bin->right)) {
                            if (el->value < BigInt(0)) { found = true; return; }
                        } else if (const auto* rl = expr_cast<RationalLit>(bin->right)) {
                            // fractional → not polynomial
                            (void)rl;
                            found = true; return;
                        } else {
                            // symbolic exponent depending on v
                            if (depends_on(bin->right, v)) { found = true; return; }
                        }
                    }
                    walk(bin->left); walk(bin->right); return;
                }
                if (const auto* un = expr_cast<Unary>(e)) { walk(un->operand); return; }
                if (const auto* s = expr_cast<Sum>(e)) {
                    for (ExprPtr t : s->terms) walk(t);
                    return;
                }
                if (const auto* p = expr_cast<Product>(e)) {
                    for (ExprPtr t : p->factors) walk(t);
                    return;
                }
                if (const auto* fc = expr_cast<FuncCall>(e)) {
                    for (ExprPtr a : fc->args) walk(a);
                    return;
                }
            }
        };
        NonPolyWalker w{var};
        w.walk(arg);
        if (w.found) return true;
    }
    return false;
}

}  // namespace

Result<ExprPtr> Integrator::integrate_product(const Product& product, const Symbol& var) {
    {  // F7.5: ∫e^{ax}·sin/cos(bx) — Laplace-style closed form.
        ExprPtr exp_arg = nullptr;
        bool exp_seen = false;
        const FuncCall* trig = nullptr;
        bool trig_seen = false;
        std::vector<ExprPtr> consts;
        bool shape_ok = true;
        for (ExprPtr f : product.factors) {
            if (!depends_on(f, var)) { consts.push_back(f); continue; }
            if (const auto* fc = expr_cast<FuncCall>(f); fc && fc->args.size() == 1U) {
                if (fc->func_id == BuiltinOp::Exp) {
                    if (exp_seen) { shape_ok = false; break; }
                    exp_seen = true;
                    exp_arg = fc->args[0];
                    continue;
                }
                if (fc->func_id == BuiltinOp::Sin || fc->func_id == BuiltinOp::Cos) {
                    if (trig_seen) { shape_ok = false; break; }
                    trig_seen = true;
                    trig = fc;
                    continue;
                }
            }
            shape_ok = false; break;
        }
        if (shape_ok && exp_seen && trig_seen) {
            auto aff_e = extract_affine_argument(exp_arg, var);
            auto aff_t = extract_affine_argument(trig->args[0], var);
            if (aff_e.has_value() && aff_t.has_value()
                && !aff_e->coefficient.numerator().is_zero()
                && !aff_t->coefficient.numerator().is_zero()) {
                const Rational a = aff_e->coefficient;
                const Rational b = aff_t->coefficient;
                const Rational denom = a * a + b * b;
                // F = e^{ax+p} · (a·trig(bx+q) ∓ b·other(bx+q)) / (a² + b²)
                ExprPtr e_factor = make_function(arena_, "exp", {exp_arg});
                ExprPtr s_fn = make_function(arena_, "sin", {trig->args[0]});
                ExprPtr c_fn = make_function(arena_, "cos", {trig->args[0]});
                ExprPtr a_expr = make_rational(arena_, a);
                ExprPtr b_expr = make_rational(arena_, b);
                ExprPtr denom_expr = make_rational(arena_, denom);
                ExprPtr inv_denom = make_binary(arena_, BinaryOp::Div,
                    make_integer(arena_, 1), denom_expr);
                ExprPtr bracket;
                if (trig->func_id == BuiltinOp::Sin) {
                    bracket = make_sum(arena_, {
                        make_product(arena_, {a_expr, s_fn}),
                        make_unary(arena_, UnaryOp::Neg,
                            make_product(arena_, {b_expr, c_fn})),
                    });
                } else {  // Cos
                    bracket = make_sum(arena_, {
                        make_product(arena_, {a_expr, c_fn}),
                        make_product(arena_, {b_expr, s_fn}),
                    });
                }
                ExprPtr primitive = make_product(arena_, {
                    inv_denom, e_factor, bracket});
                if (consts.empty()) return ok(primitive);
                consts.push_back(primitive);
                return ok(make_product(arena_, std::move(consts)));
            }
        }
    }

    // HC-IBP-VDU: detect rational-function shape  Π c_i · Π p_j(x) · Π q_k(x)^{-1}
    // upfront and route through Div(N, D) → poly-divide + PFD.
    // Without this gate, the recursive u-substitution path can return a
    // wrong-but-internally-self-consistent stub (e.g. (1/2)x for
    // x²/(2(1+x²))) on shapes IBP produces as v·du.
    {
        auto is_poly = [&](ExprPtr e) {
            return algebra::polynomial_degree(e, var, context_).is_ok();
        };
        std::vector<ExprPtr> num_factors;
        std::vector<ExprPtr> den_factors;
        bool all_rational = true;
        for (ExprPtr f : product.factors) {
            if (!depends_on(f, var)) { num_factors.push_back(f); continue; }
            if (const auto* pw = expr_cast<Binary>(f);
                pw && pw->op == BinaryOp::Pow) {
                if (const auto* el = expr_cast<IntegerLit>(pw->right)) {
                    // Restrict to single-power denominators (Pow exp = -1).
                    // Repeated factors (exp ≤ -2) are best served by the
                    // Hermite reduction path; bypassing it here regresses
                    // CalculusIntegrateTest.P1_HermiteReduction_*.
                    if (el->value == BigInt(-1) && is_poly(pw->left)) {
                        den_factors.push_back(pw->left);
                        continue;
                    }
                    if (!el->value.is_negative() && is_poly(f)) {
                        num_factors.push_back(f);
                        continue;
                    }
                }
                all_rational = false; break;
            }
            if (is_poly(f)) { num_factors.push_back(f); continue; }
            all_rational = false; break;
        }
        if (all_rational && !den_factors.empty()) {
            ExprPtr N = num_factors.empty()
                ? arena_.make<IntegerLit>(BigInt(1))
                : (num_factors.size() == 1U ? num_factors[0]
                    : arena_.make<Product>(std::move(num_factors)));
            ExprPtr D = den_factors.size() == 1U ? den_factors[0]
                : arena_.make<Product>(std::move(den_factors));
            ExprPtr div = arena_.make<Binary>(BinaryOp::Div, N, D);
            // Only trigger when numerator degree ≥ denominator degree.
            // For deg(N) < deg(D) the existing dispatch (PFD via the Binary(Div)
            // path or Hermite) already handles the integrand correctly;
            // preempting here causes regressions on ln/log integrands whose
            // IBP path produces Product([N, Pow(D,-1)]) sub-integrands that
            // need to flow through the standard rational pipeline.
            auto dn = algebra::polynomial_degree(N, var, context_);
            auto dd = algebra::polynomial_degree(D, var, context_);
            if (dn.is_ok() && dd.is_ok() && dn.value() >= dd.value()) {
                auto rat = integrate_via_partial_fractions(div, var);
                if (rat.is_ok()) return rat;
                // Polynomial-divide: N = D·Q + R → ∫N/D = ∫Q + ∫R/D.
                auto dm = algebra::polynomial_divmod(N, D, var, context_);
                if (dm.is_ok()) {
                    ExprPtr leftover = arena_.make<Binary>(BinaryOp::Div,
                        dm.value().remainder, D);
                    auto q_int = integrate(dm.value().quotient, var);
                    auto r_int = integrate(leftover, var);
                    if (q_int.is_ok() && r_int.is_ok()) {
                        return ok(arena_.make<Sum>(std::vector<ExprPtr>{
                            q_int.value(), r_int.value()}));
                    }
                }
            }
        }
    }

    // HC-IBP-VDU: detect shape  c · x² · sqrt(R)^{-1}  →
    //   c · integrate_xsq_over_sqrt_quadratic(R).
    // R(x) recognised by match_constant_square_minus_variable_square etc.
    {
        std::vector<ExprPtr> consts;
        ExprPtr xsq = nullptr;
        ExprPtr inv_sqrt_arg = nullptr;
        bool shape_ok = true;
        for (ExprPtr f : product.factors) {
            if (!depends_on(f, var)) { consts.push_back(f); continue; }
            // x² factor?
            if (const auto* pw = expr_cast<Binary>(f);
                pw && pw->op == BinaryOp::Pow && !xsq) {
                if (auto* sym = expr_cast<Symbol>(pw->left);
                    sym && sym->name == var.name) {
                    if (const auto* el = expr_cast<IntegerLit>(pw->right);
                        el && el->value == BigInt(2)) {
                        xsq = f; continue;
                    }
                }
            }
            // 1/sqrt(R)?  Accept both Pow(sqrt(R), -1) and Div(1, sqrt(R)).
            if (const auto* pw = expr_cast<Binary>(f);
                pw && pw->op == BinaryOp::Pow && !inv_sqrt_arg) {
                if (const auto* el = expr_cast<IntegerLit>(pw->right);
                    el && el->value == BigInt(-1)) {
                    if (const auto* fc = expr_cast<FuncCall>(pw->left);
                        fc && fc->func_id == BuiltinOp::Sqrt
                        && fc->args.size() == 1U) {
                        inv_sqrt_arg = fc->args[0]; continue;
                    }
                }
            }
            if (const auto* dv = expr_cast<Binary>(f);
                dv && dv->op == BinaryOp::Div && !inv_sqrt_arg) {
                if (const auto* nu = expr_cast<IntegerLit>(dv->left);
                    nu && nu->value == BigInt(1)) {
                    if (const auto* fc = expr_cast<FuncCall>(dv->right);
                        fc && fc->func_id == BuiltinOp::Sqrt
                        && fc->args.size() == 1U) {
                        inv_sqrt_arg = fc->args[0]; continue;
                    }
                }
            }
            shape_ok = false; break;
        }
        if (shape_ok && xsq && inv_sqrt_arg) {
            auto primitive = integrate_xsq_over_sqrt_quadratic(inv_sqrt_arg, var);
            if (primitive.is_ok()) {
                if (consts.empty()) return primitive;
                consts.push_back(primitive.value());
                return ok(make_product(arena_, std::move(consts)));
            }
        }
    }

    auto substitution = try_u_substitution_for_product(product, var);
    if (substitution.is_ok()) {
        return substitution;
    }

    std::vector<ExprPtr> constant_factors;
    std::vector<ExprPtr> variable_factors;
    constant_factors.reserve(product.factors.size());
    variable_factors.reserve(product.factors.size());

    for (ExprPtr factor : product.factors) {
        if (depends_on(factor, var)) {
            variable_factors.push_back(factor);
        } else {
            constant_factors.push_back(factor);
        }
    }

    if (variable_factors.empty()) {
        constant_factors.push_back(arena_.make<Symbol>(var));
        return ok(make_product(arena_, std::move(constant_factors)));
    }

    if (constant_factors.size() == 1U && is_negative_one(constant_factors.front()) &&
        variable_factors.size() == 1U && matches_reciprocal_sqrt_one_minus_square(variable_factors.front(), var)) {
        return ok(make_function(arena_, "arccos", {arena_.make<Symbol>(var)}));
    }

    for (std::size_t index = 0; index < variable_factors.size(); ++index) {
        const auto* power = expr_cast<Binary>(variable_factors[index]);
        if (power == nullptr || power->op != BinaryOp::Pow || !is_rational_value(power->right, -1, 1)) {
            continue;
        }

        std::vector<ExprPtr> numerator_factors = constant_factors;
        for (std::size_t factor_index = 0; factor_index < variable_factors.size(); ++factor_index) {
            if (factor_index != index) {
                numerator_factors.push_back(variable_factors[factor_index]);
            }
        }
        ExprPtr numerator = numerator_factors.empty()
            ? make_integer(arena_, 1)
            : make_product(arena_, std::move(numerator_factors));
        Binary quotient{BinaryOp::Div, numerator, power->left};
        if (auto quadratic_integral = integrate_linear_over_quadratic(quotient, var);
            quadratic_integral.is_ok()) {
            return quadratic_integral;
        }
    }

    if (variable_factors.size() >= 2U) {
        if (auto rational_integral = integrate_via_partial_fractions(arena_.make<Product>(product.factors), var);
            rational_integral.is_ok()) {
            if (constant_factors.empty()) {
                return rational_integral;
            }
            std::vector<ExprPtr> result_factors = constant_factors;
            result_factors.push_back(rational_integral.value());
            return ok(make_product(arena_, std::move(result_factors)));
        }

        // BUG-HANG-001: skip IBP when a factor is exp(g) with g non-polynomial
        // in var (e.g. exp(1/x)).  IBP generates integrands with growing
        // denominator degree and never terminates.  Risch DE solver handles
        // these cases in step 3 of Integrator::integrate.
        if (has_exp_rational_non_poly_factor(variable_factors, var)) {
            return fail<ExprPtr>(make_error(
                CASErrorKind::Unimplemented,
                "integrate_product: IBP skipped for exp(non-poly) factor; defer to Risch DE"));
        }

        auto ibp_res = integrate_by_parts(arena_.make<Product>(variable_factors), var, context_);
        if (ibp_res.is_ok()) {
            if (constant_factors.empty()) {
                return ibp_res;
            }
            std::vector<ExprPtr> result_factors = constant_factors;
            result_factors.push_back(ibp_res.value());
            return ok(make_product(arena_, std::move(result_factors)));
        }
    }

    if (variable_factors.size() == 1U) {
        auto inner = integrate_once(variable_factors.front(), var);
        if (inner.is_error()) {
            auto ibp_res = integrate_by_parts(make_product(arena_, {make_integer(arena_, 1), variable_factors.front()}), var, context_);
            if (ibp_res.is_ok()) {
                if (constant_factors.empty()) {
                    return ibp_res;
                }
                std::vector<ExprPtr> result_factors = constant_factors;
                result_factors.push_back(ibp_res.value());
                return ok(make_product(arena_, std::move(result_factors)));
            }
            return inner;
        }

        if (constant_factors.empty()) {
            return inner;
        }

        constant_factors.push_back(inner.value());
        return ok(make_product(arena_, std::move(constant_factors)));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Integration by parts or substitution is not implemented for this product"));
}

Result<ExprPtr> Integrator::integrate_power(const Binary& power, const Symbol& var) {
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
