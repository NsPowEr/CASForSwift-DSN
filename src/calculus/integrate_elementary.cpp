#include "integrate_engine.hpp"
#include "cas/algebra.hpp"

#include <string>
#include <utility>

namespace cas::calculus::integrate_detail {

Result<ExprPtr> Integrator::integrate_function_direct(const std::string& name, ExprPtr argument) {
    BuiltinOp func_id = get_builtin_op(name);
    if (func_id == BuiltinOp::Sin) {
        return ok(make_product(arena_, {make_integer(arena_, -1), make_function(arena_, "cos", {argument})}));
    }
    if (func_id == BuiltinOp::Cos) {
        return ok(make_function(arena_, "sin", {argument}));
    }
    if (func_id == BuiltinOp::Tan) {
        return ok(make_product(arena_, {make_integer(arena_, -1), make_function(arena_, "ln", {make_function(arena_, "abs", {make_function(arena_, "cos", {argument})})})}));
    }
    if (func_id == BuiltinOp::Cot) {
        return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_function(arena_, "sin", {argument})})}));
    }
    if (func_id == BuiltinOp::Sec) {
        return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {make_function(arena_, "sec", {argument}), make_function(arena_, "tan", {argument})})})}));
    }
    if (func_id == BuiltinOp::Csc) {
        return ok(make_product(arena_, {make_integer(arena_, -1), make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {make_function(arena_, "csc", {argument}), make_function(arena_, "cot", {argument})})})})}));
    }
    if (func_id == BuiltinOp::Exp) {
        return ok(make_function(arena_, "exp", {argument}));
    }
    if (func_id == BuiltinOp::Sinh) {
        return ok(make_function(arena_, "cosh", {argument}));
    }
    if (func_id == BuiltinOp::Cosh) {
        return ok(make_function(arena_, "sinh", {argument}));
    }
    if (func_id == BuiltinOp::Tanh) {
        // ∫tanh(x) dx = ln(cosh(x))
        return ok(make_function(arena_, "ln",
            {make_function(arena_, "cosh", {argument})}));
    }
    if (func_id == BuiltinOp::Ln || func_id == BuiltinOp::Log) {
        // Both Ln and Log are natural log in this engine (see
        // differentiate.cpp same fix). ∫ln(x) dx = x·ln(x) - x.
        ExprPtr x = argument;
        return ok(make_sum(arena_, {make_product(arena_, {x, make_function(arena_, "ln", {x})}), make_unary(arena_, UnaryOp::Neg, x)}));
    }
    if (func_id == BuiltinOp::Atan) {
        ExprPtr x = argument;
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "arctan", {x})}),
            make_product(arena_, {
                make_rational(arena_, -1, 2),
                make_function(arena_, "ln", {make_sum(arena_, {
                    make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2)),
                    make_integer(arena_, 1),
                })}),
            }),
        }));
    }
    if (func_id == BuiltinOp::Asin) {
        // ∫asin(x) dx = x·asin(x) + sqrt(1 - x²)
        ExprPtr x = argument;
        ExprPtr one_minus_x2 = make_sum(arena_, {
            make_integer(arena_, 1),
            make_unary(arena_, UnaryOp::Neg,
                make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2))),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "arcsin", {x})}),
            make_function(arena_, "sqrt", {one_minus_x2}),
        }));
    }
    if (func_id == BuiltinOp::Acos) {
        // ∫acos(x) dx = x·acos(x) - sqrt(1 - x²)
        ExprPtr x = argument;
        ExprPtr one_minus_x2 = make_sum(arena_, {
            make_integer(arena_, 1),
            make_unary(arena_, UnaryOp::Neg,
                make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2))),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "arccos", {x})}),
            make_unary(arena_, UnaryOp::Neg,
                make_function(arena_, "sqrt", {one_minus_x2})),
        }));
    }
    // F7.5.B1: inverse hyperbolic standalone integrals.
    // These functions are not in BuiltinOp (would touch 76 switch
    // statements under -Wswitch -Werror, deferred to Fase 8); they are
    // parsed as FuncCall(name, …) with BuiltinOp::Custom. Match by
    // canonical name (and 'arc' aliases Maxima sometimes emits).
    auto matches_name = [&](std::initializer_list<const char*> names) {
        for (auto* n : names) if (name == n) return true;
        return false;
    };
    if (matches_name({"asinh", "arcsinh"})) {
        // ∫asinh(x) dx = x·asinh(x) - sqrt(x² + 1)
        ExprPtr x = argument;
        ExprPtr x2_plus_1 = make_sum(arena_, {
            make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2)),
            make_integer(arena_, 1),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "asinh", {x})}),
            make_unary(arena_, UnaryOp::Neg, make_function(arena_, "sqrt", {x2_plus_1})),
        }));
    }
    if (matches_name({"acosh", "arccosh"})) {
        // ∫acosh(x) dx = x·acosh(x) - sqrt(x² - 1)
        ExprPtr x = argument;
        ExprPtr x2_minus_1 = make_sum(arena_, {
            make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2)),
            make_integer(arena_, -1),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "acosh", {x})}),
            make_unary(arena_, UnaryOp::Neg, make_function(arena_, "sqrt", {x2_minus_1})),
        }));
    }
    if (matches_name({"atanh", "arctanh"})) {
        // ∫atanh(x) dx = x·atanh(x) + ½·ln(1 - x²)
        ExprPtr x = argument;
        ExprPtr one_minus_x2 = make_sum(arena_, {
            make_integer(arena_, 1),
            make_unary(arena_, UnaryOp::Neg,
                make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2))),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "atanh", {x})}),
            make_product(arena_, {
                make_rational(arena_, 1, 2),
                make_function(arena_, "ln", {one_minus_x2}),
            }),
        }));
    }
    if (matches_name({"acoth", "arccoth"})) {
        // ∫acoth(x) dx = x·acoth(x) + ½·ln(x² - 1)
        ExprPtr x = argument;
        ExprPtr x2_minus_1 = make_sum(arena_, {
            make_binary(arena_, BinaryOp::Pow, x, make_integer(arena_, 2)),
            make_integer(arena_, -1),
        });
        return ok(make_sum(arena_, {
            make_product(arena_, {x, make_function(arena_, "acoth", {x})}),
            make_product(arena_, {
                make_rational(arena_, 1, 2),
                make_function(arena_, "ln", {x2_minus_1}),
            }),
        }));
    }
    if (func_id == BuiltinOp::Sqrt) {
        return integrate_power_direct(argument, make_rational(arena_, 1, 2), Symbol("_u_"));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "This elementary function integral is not implemented"));
}

Result<ExprPtr> Integrator::integrate_power_direct(ExprPtr base, ExprPtr exponent, const Symbol& var) {
    if (depends_on(exponent, var)) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Only powers with constant exponent support direct substitution"));
    }

    if (is_negative_one(exponent)) {
        return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {base})}));
    }

    if (const auto* integer = expr_cast<IntegerLit>(exponent)) {
        const BigInt next = integer->value + BigInt(1);
        return ok(make_product(arena_, {arena_.make<RationalLit>(BigInt(1), next), make_binary(arena_, BinaryOp::Pow, base, arena_.make<IntegerLit>(next))}));
    }

    ExprPtr exponent_plus_one = make_sum(arena_, {exponent, make_integer(arena_, 1)});
    return ok(make_binary(arena_, BinaryOp::Div, make_binary(arena_, BinaryOp::Pow, base, exponent_plus_one), exponent_plus_one));
}

namespace {

// True if `e` carries a genuine radical in `var`: sqrt of a var-dependent
// argument, or a fractional power of a var-dependent base.
[[nodiscard]] bool has_var_radical(ExprPtr e, const Symbol& var) {
    if (!e) return false;
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        if (fc->func_id == BuiltinOp::Sqrt)
            for (ExprPtr a : fc->args) if (depends_on(a, var)) return true;
        for (ExprPtr a : fc->args) if (has_var_radical(a, var)) return true;
        return false;
    }
    if (const auto* b = expr_cast<Binary>(e)) {
        if (b->op == BinaryOp::Pow && depends_on(b->left, var)
            && expr_is<RationalLit>(b->right)) return true;
        return has_var_radical(b->left, var) || has_var_radical(b->right, var);
    }
    if (const auto* u = expr_cast<Unary>(e)) return has_var_radical(u->operand, var);
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms) if (has_var_radical(t, var)) return true;
        return false;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr f : p->factors) if (has_var_radical(f, var)) return true;
        return false;
    }
    return false;
}

// Conjugate of a radical-sum: negate exactly the terms carrying a var-radical.
// For S = A + B√q (A radical-free, B√q the radical part) returns A − B√q.
[[nodiscard]] ExprPtr radical_sum_conjugate(const Sum& s, const Symbol& var, AstArena& arena) {
    std::vector<ExprPtr> terms;
    terms.reserve(s.terms.size());
    for (ExprPtr t : s.terms)
        terms.push_back(has_var_radical(t, var) ? make_unary(arena, UnaryOp::Neg, t) : t);
    return make_sum(arena, std::move(terms));
}

// Canonical negation of an explicit integer/rational exponent literal, so the
// rewritten power never carries a non-canonical Neg(IntegerLit(neg)) node.
[[nodiscard]] ExprPtr negate_exponent(ExprPtr e, AstArena& arena) {
    if (const auto* il = expr_cast<IntegerLit>(e)) return arena.make<IntegerLit>(-il->value);
    if (const auto* rl = expr_cast<RationalLit>(e)) return arena.make<RationalLit>(-rl->numerator, rl->denominator);
    return make_unary(arena, UnaryOp::Neg, e);
}

// General conjugate rationalisation of radical-sum denominators (Regola Zero:
// algorithm, not pattern table).  Recursively rewrites every factor Pow(S, n<0)
// whose base S is a Sum carrying a var-radical:
//   Pow(S, n) → conj^{-n} · Pow(S·conj, n),   conj = radical terms of S negated,
// since for S = A + B√q the product S·conj = A² − B²q is radical-free.  The
// rewrite is applied ONLY when S·conj actually loses the radical (guard), so a
// multi-radical denominator like √2+√3 — which a single conjugation does not
// clear — is left untouched.  The caller re-simplifies and verifies
// D(result) ≡ integrand, so a wrong rewrite can never leak a silent wrong answer.
[[nodiscard]] ExprPtr rationalize_radical_denominators(
    ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    if (!e) return e;
    AstArena& arena = ctx.arena();
    if (const auto* b = expr_cast<Binary>(e)) {
        if (b->op == BinaryOp::Pow) {
            bool neg = false;
            if (const auto* il = expr_cast<IntegerLit>(b->right)) neg = il->value.is_negative();
            else if (const auto* rl = expr_cast<RationalLit>(b->right)) neg = rl->numerator.is_negative();
            const auto* base_sum = expr_cast<Sum>(b->left);
            if (neg && base_sum != nullptr && has_var_radical(b->left, var)) {
                ExprPtr conj = radical_sum_conjugate(*base_sum, var, arena);
                auto rs = ctx.simplify(make_product(arena, {b->left, conj}));
                ExprPtr R = rs.is_ok() ? rs.value() : nullptr;
                if (R != nullptr && !has_var_radical(R, var)) {
                    ExprPtr conj_pow = make_binary(arena, BinaryOp::Pow, conj,
                        negate_exponent(b->right, arena));
                    ExprPtr R_pow = make_binary(arena, BinaryOp::Pow, R, b->right);
                    return make_product(arena, {conj_pow, R_pow});
                }
            }
        }
        ExprPtr l = rationalize_radical_denominators(b->left, var, ctx);
        ExprPtr r = rationalize_radical_denominators(b->right, var, ctx);
        return (l == b->left && r == b->right) ? e : make_binary(arena, b->op, l, r);
    }
    if (const auto* u = expr_cast<Unary>(e)) {
        ExprPtr o = rationalize_radical_denominators(u->operand, var, ctx);
        return (o == u->operand) ? e : make_unary(arena, u->op, o);
    }
    if (const auto* s = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> terms; terms.reserve(s->terms.size());
        bool changed = false;
        for (ExprPtr t : s->terms) {
            ExprPtr nt = rationalize_radical_denominators(t, var, ctx);
            changed = changed || (nt != t);
            terms.push_back(nt);
        }
        return changed ? make_sum(arena, std::move(terms)) : e;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        std::vector<ExprPtr> factors; factors.reserve(p->factors.size());
        bool changed = false;
        for (ExprPtr f : p->factors) {
            ExprPtr nf = rationalize_radical_denominators(f, var, ctx);
            changed = changed || (nf != f);
            factors.push_back(nf);
        }
        return changed ? make_product(arena, std::move(factors)) : e;
    }
    return e;
}

// Distribute a numeric literal coefficient out of a power whose base is a
// Product containing a var-radical: (c·rest)^n → c^n · Pow(rest, n).  After
// conjugate rationalisation, `together` can leave the parts remainder as
// c·x·(c·√(x²+a))⁻¹ (c = (S·conj)² , e.g. 16 when S·conj=−4) whose constant the
// universal simplifier does not pull out of the reciprocal (T-054, deferred
// there because doing it in the hot power/product path is costly and regresses
// Gruntz).  Lifting it HERE — once, on the integration remainder — lets the
// constant cancel with its sibling so the remainder integrates, with no
// hot-path cost and no Gruntz interaction.  No-op when the radical base has no
// numeric coefficient (e.g. the a=1 asinh remainder x·√(x²+1)⁻¹).
[[nodiscard]] ExprPtr lift_numeric_from_radical_powers(
    ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    if (!e) return e;
    AstArena& arena = ctx.arena();
    if (const auto* b = expr_cast<Binary>(e)) {
        if (b->op == BinaryOp::Pow && expr_is<IntegerLit>(b->right)) {
            if (const auto* bp = expr_cast<Product>(b->left)) {
                std::vector<ExprPtr> nums, rest;
                bool rest_radical = false;
                for (ExprPtr bf : bp->factors) {
                    if (expr_is<IntegerLit>(bf) || expr_is<RationalLit>(bf)) {
                        nums.push_back(bf);
                    } else {
                        rest.push_back(bf);
                        if (has_var_radical(bf, var)) rest_radical = true;
                    }
                }
                if (!nums.empty() && !rest.empty() && rest_radical) {
                    ExprPtr num = nums.size() == 1 ? nums[0] : make_product(arena, std::move(nums));
                    ExprPtr rst = rest.size() == 1 ? rest[0] : make_product(arena, std::move(rest));
                    return make_product(arena, {
                        make_binary(arena, BinaryOp::Pow, num, b->right),
                        make_binary(arena, BinaryOp::Pow, rst, b->right)});
                }
            }
        }
        ExprPtr l = lift_numeric_from_radical_powers(b->left, var, ctx);
        ExprPtr r = lift_numeric_from_radical_powers(b->right, var, ctx);
        return (l == b->left && r == b->right) ? e : make_binary(arena, b->op, l, r);
    }
    if (const auto* u = expr_cast<Unary>(e)) {
        ExprPtr o = lift_numeric_from_radical_powers(u->operand, var, ctx);
        return (o == u->operand) ? e : make_unary(arena, u->op, o);
    }
    if (const auto* s = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> terms; terms.reserve(s->terms.size());
        bool changed = false;
        for (ExprPtr t : s->terms) {
            ExprPtr nt = lift_numeric_from_radical_powers(t, var, ctx);
            changed = changed || (nt != t);
            terms.push_back(nt);
        }
        return changed ? make_sum(arena, std::move(terms)) : e;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        std::vector<ExprPtr> factors; factors.reserve(p->factors.size());
        bool changed = false;
        for (ExprPtr f : p->factors) {
            ExprPtr nf = lift_numeric_from_radical_powers(f, var, ctx);
            changed = changed || (nf != f);
            factors.push_back(nf);
        }
        return changed ? make_product(arena, std::move(factors)) : e;
    }
    return e;
}

}  // namespace

Result<ExprPtr> Integrator::integrate_function(const FuncCall& call, const Symbol& var) {
    if (call.args.size() != 1U) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Only unary function integration is implemented"));
    }

    ExprPtr argument = call.args.front();
    BuiltinOp func_id = call.func_id;
    if (func_id == BuiltinOp::Sqrt) {
        // F7.5.B1: try quadratic argument first; fall back to power_direct
        // for sqrt(affine in x) — `sqrt(x)` = x^(1/2), affine path scales.
        auto quad = integrate_sqrt_quadratic(argument, var);
        if (quad.is_ok()) return quad;
        auto affine = extract_affine_argument(argument, var);
        if (affine.has_value() && affine->coefficient.numerator() != BigInt(0)) {
            // ∫sqrt(a·x + b) dx = (1/a) · (2/3) · (a·x + b)^(3/2)
            auto primitive = integrate_power_direct(argument,
                make_rational(arena_, 1, 2), var);
            if (primitive.is_ok()) {
                if (affine->coefficient == Rational(BigInt(1))) {
                    return primitive;
                }
                return ok(make_product(arena_, {
                    make_rational(arena_, Rational(BigInt(1)) / affine->coefficient),
                    primitive.value(),
                }));
            }
        }
        return quad;  // propagate original failure
    }

    auto affine = extract_affine_argument(argument, var);
    if (!affine.has_value() || affine->coefficient.numerator() == BigInt(0)) {
        // F7.5: IBP fallback for non-affine argument.
        //   ∫ f(g(x)) dx  with u = f(g(x)), v = x  →
        //     x·f(g(x)) − ∫ x · f'(g(x)) · g'(x) dx
        // Effective on log(poly), asin/acos/atan/atanh(poly), etc.
        // Reconstruct the call as ExprPtr.  Direct address-of would be safer
        // but is not exposed; rebuild using known constructor.
        ExprPtr u = arena_.make<FuncCall>(call.func_id, std::vector<ExprPtr>{argument});
        ExprPtr v = arena_.make<Symbol>(var);
        auto du_res = diff(u, var, 1U, context_);
        if (du_res.is_ok()) {
            ExprPtr du_simp = du_res.value();
            if (auto ds = context_.simplify(du_simp); ds.is_ok())
                du_simp = ds.value();
            ExprPtr v_du = arena_.make<Product>(std::vector<ExprPtr>{v, du_simp});
            auto vdu_simp = context_.simplify(v_du);
            ExprPtr vdu_in = vdu_simp.is_ok() ? vdu_simp.value() : v_du;
            // HC-IBP-RADSUM-RATIONALIZE: the parts remainder x·g'/g of
            // ∫f(g(x))dx can carry a radical-sum denominator (e.g.
            // ∫log(x+√(x²+1)) leaves x·(x+√(x²+1))^{-1}·…) that simplify alone
            // does not cancel.  Conjugate rationalisation clears the radical;
            // together() then combines over a common denominator (x·√(x²+1) −
            // x³/√(x²+1) → x/√(x²+1)) so the remainder integrates.
            ExprPtr vdu_rat = rationalize_radical_denominators(vdu_in, var, context_);
            if (vdu_rat != vdu_in) {
                if (auto tg = algebra::together(vdu_rat, context_); tg.is_ok())
                    vdu_rat = tg.value();
                // together can strand a numeric coefficient inside a radical
                // reciprocal (c·x·(c·√(x²+a))⁻¹ for a≠1); lift it so it cancels.
                vdu_rat = lift_numeric_from_radical_powers(vdu_rat, var, context_);
                auto rs = context_.simplify(vdu_rat);
                vdu_in = rs.is_ok() ? rs.value() : vdu_rat;
            }
            auto int_vdu = ::cas::calculus::integrate(vdu_in, var, context_);
            if (int_vdu.is_ok()) {
                ExprPtr candidate = arena_.make<Sum>(std::vector<ExprPtr>{
                    arena_.make<Product>(std::vector<ExprPtr>{v, u}),
                    arena_.make<Unary>(UnaryOp::Neg, int_vdu.value())});
                // Verify D(candidate) ≡ original integrand (call.func_id(arg)).
                auto D_res = diff(candidate, var, 1U, context_);
                if (D_res.is_ok()) {
                    ExprPtr orig = arena_.make<FuncCall>(call.func_id,
                        std::vector<ExprPtr>{argument});
                    ExprPtr delta = arena_.make<Binary>(BinaryOp::Sub,
                        D_res.value(), orig);
                    auto check_zero_e = [&](ExprPtr e) {
                        return expr_is<IntegerLit>(e)
                            && expr_ref<IntegerLit>(e).value.is_zero();
                    };
                    // Rationalise any radical-sum denominator in the residual
                    // first (no-op when absent): for a log-of-radical-sum the
                    // residual still carries the uncancellable g'/g fraction
                    // (x·g'/g − x/√(x²+1)); conjugate rationalisation lets it
                    // collapse to 0 instead of the simplifier churning on it.
                    ExprPtr delta_r = rationalize_radical_denominators(delta, var, context_);
                    auto delta_simp = context_.simplify(delta_r);
                    if (delta_simp.is_ok() && check_zero_e(delta_simp.value()))
                        return ok(candidate);
                    // Fallback: combine fractions then re-simplify; covers the
                    // common log(p(x)) IBP shape where the delta is a sum of
                    // rational terms whose numerator collapses to zero but the
                    // canonicaliser does not group them automatically.
                    if (auto tg = algebra::together(delta_r, context_); tg.is_ok())
                        if (auto ts = context_.simplify(tg.value()); ts.is_ok())
                            if (check_zero_e(ts.value())) return ok(candidate);
                }
            }
        }
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
            "Function integration: non-affine arg, IBP fallback failed"));
    }

    auto primitive = integrate_function_direct(call.name, argument);
    if (primitive.is_error()) {
        return primitive;
    }

    if (affine->coefficient == Rational(BigInt(1))) {
        return primitive;
    }

    return ok(make_product(arena_, {make_rational(arena_, Rational(BigInt(1)) / affine->coefficient), primitive.value()}));
}

}  // namespace cas::calculus::integrate_detail
