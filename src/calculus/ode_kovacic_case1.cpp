// Kovacic Case 1 — rational Liouvillian solutions of z'' = r·z.
// Ref: Kovacic J.J. (1986) §3, "An Algorithm for Solving Second Order Linear
// Homogeneous Differential Equations" J. Symbolic Computation 2, 3–43.
//
// Seeks ω ∈ Q(x) s.t. ω' + ω² = r, giving z₁ = exp(∫ω dx).
// Returns OmegaPair{ω₊, ω₋} for the two sign choices at each order-2 pole.
// Fails with Unimplemented if any condition for Case 1 is violated.

#include "ode_kovacic_internal.hpp"
#include <optional>
#include <string>

namespace cas::calculus {
namespace kovacic_impl {

// Helper functions moved to ode_kovacic_laurent.cpp

Result<ExprPtr> compute_r(
    ExprPtr a2, ExprPtr a1, ExprPtr a0,
    const Symbol& x, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();
    ExprPtr p = kv_div(a, a1, a2);
    ExprPtr q = kv_div(a, a0, a2);

    auto dp_res = diff(p, x, 1U, ctx);
    if (dp_res.is_error()) return dp_res;

    ExprPtr r_raw = kv_sub(a,
        kv_add(a,
            kv_div(a, dp_res.value(), kv_int(a, 2)),
            kv_div(a, kv_mul(a, p, p), kv_int(a, 4))),
        q);

    return ctx.simplify(r_raw);
}

[[nodiscard]] static std::optional<ExprPtr> factor_pole(
    ExprPtr factor, const Symbol& x, AstArena& a) {

    if (auto* sym = expr_cast<Symbol>(factor))
        if (sym->name == x.name) return kv_int(a, 0);

    if (auto* bin = expr_cast<Binary>(factor)) {
        auto* lsym = expr_cast<Symbol>(bin->left);
        if (!lsym || lsym->name != x.name) return std::nullopt;
        if (bin->op == BinaryOp::Sub) return bin->right;
        if (bin->op == BinaryOp::Add) return kv_neg(a, bin->right);
    }
    return std::nullopt;
}

// A pole factor (x − c)^{−k} as emitted by simplify/partial_fractions, i.e.
// Pow(base, IntegerLit(neg)).  Returns the pole base and the (positive) order k.
[[nodiscard]] static std::optional<std::pair<ExprPtr, unsigned>> match_neg_power(
    ExprPtr f, const Symbol& x, AstArena& a) {

    auto* pw = expr_cast<Binary>(f);
    if (!pw || pw->op != BinaryOp::Pow) return std::nullopt;
    auto* el = expr_cast<IntegerLit>(pw->right);
    if (!el || !el->value.is_negative()) return std::nullopt;

    const long long ev = el->value.abs().to_double();
    if (ev < 1 || ev > 1024) return std::nullopt;

    auto pole_opt = factor_pole(pw->left, x, a);
    if (!pole_opt) return std::nullopt;
    return std::make_pair(*pole_opt, static_cast<unsigned>(ev));
}

[[nodiscard]] static std::optional<PFPole> try_extract_pole(
    ExprPtr term, const Symbol& x, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();

    // Form 1 — explicit division  num / (x−c)^k  (legacy representation).
    if (auto* dv = expr_cast<Binary>(term); dv && dv->op == BinaryOp::Div) {
        ExprPtr num  = dv->left;
        ExprPtr den  = dv->right;
        unsigned pow = 1;
        ExprPtr base = den;

        if (auto* pw = expr_cast<Binary>(den)) {
            if (pw->op == BinaryOp::Pow) {
                if (auto* el = expr_cast<IntegerLit>(pw->right)) {
                    if (el->value.is_negative() || el->value.is_zero())
                        return std::nullopt;
                    long long ev = el->value.to_double();
                    if (ev < 1 || ev > 1024) return std::nullopt;
                    pow  = static_cast<unsigned>(ev);
                    base = pw->left;
                } else { return std::nullopt; }
            }
        }

        auto pole_opt = factor_pole(base, x, a);
        if (!pole_opt) return std::nullopt;
        auto ps = ctx.simplify(*pole_opt);
        return PFPole{ps.is_ok() ? ps.value() : *pole_opt, pow, num};
    }

    // Form 2 — coeff·(x−c)^{−k}: the negative-power Pow representation that
    // simplify/partial_fractions actually emit.  The (optional) coefficient is
    // the product of all non-pole factors.
    ExprPtr coeff = kv_int(a, 1);
    ExprPtr pole_factor = nullptr;

    if (auto* prod = expr_cast<Product>(term)) {
        for (ExprPtr f : prod->factors) {
            if (match_neg_power(f, x, a)) {
                if (pole_factor) return std::nullopt;  // ≥2 pole factors — unsupported
                pole_factor = f;
            } else {
                coeff = kv_mul(a, coeff, f);
            }
        }
    } else {
        pole_factor = term;
    }
    if (!pole_factor) return std::nullopt;

    auto pp = match_neg_power(pole_factor, x, a);
    if (!pp) return std::nullopt;

    auto ps = ctx.simplify(pp->first);
    auto cs = ctx.simplify(coeff);
    return PFPole{ps.is_ok() ? ps.value() : pp->first, pp->second,
                  cs.is_ok() ? cs.value() : coeff};
}

Result<OmegaPair> case1_omega(
    ExprPtr r, const Symbol& x, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();

    auto parts = algebra::apart_num_den(r, ctx);
    if (parts.is_error()) return fail<OmegaPair>(parts.error());

    ExprPtr num = parts.value().numerator;
    ExprPtr den = parts.value().denominator;

    auto dn_res = algebra::polynomial_degree(num, x, ctx);
    auto dd_res = algebra::polynomial_degree(den, x, ctx);
    if (dn_res.is_error()) return fail<OmegaPair>(dn_res.error());
    if (dd_res.is_error()) return fail<OmegaPair>(dd_res.error());

    ExprPtr poly_op = kv_int(a, 0);
    ExprPtr poly_om = kv_int(a, 0);
    ExprPtr r_proper = r;

    auto build_laurent_inf_poly = [&](const std::vector<Rational>& s_prime, unsigned v, bool positive_branch) -> ExprPtr {
        ExprPtr res = kv_int(a, 0);
        for (unsigned k = 0; k <= v; ++k) {
            Rational coeff = s_prime[k];
            if (!positive_branch) coeff = coeff * Rational(BigInt(-1));
            if (coeff.numerator().is_zero()) continue;
            
            ExprPtr coeff_expr = coeff.is_integer()
                ? static_cast<ExprPtr>(a.make<IntegerLit>(coeff.numerator()))
                : static_cast<ExprPtr>(a.make<RationalLit>(coeff.numerator(), coeff.denominator()));
                
            ExprPtr term = coeff_expr;
            if (k < v) {
                ExprPtr xpow = (v - k == 1)
                    ? static_cast<ExprPtr>(a.make<Symbol>(x.name))
                    : a.make<Binary>(BinaryOp::Pow, a.make<Symbol>(x.name), kv_int(a, v - k));
                term = kv_mul(a, term, xpow);
            }
            res = kv_is_zero(res, ctx) ? term : kv_add(a, res, term);
        }
        return res;
    };

    if (dn_res.value() >= dd_res.value()) {
        auto dm = algebra::polynomial_divmod(num, den, x, ctx);
        if (dm.is_error()) return fail<OmegaPair>(dm.error());

        ExprPtr quotient  = dm.value().quotient;
        ExprPtr remainder = dm.value().remainder;

        auto dq_res = algebra::polynomial_degree(quotient, x, ctx);
        if (dq_res.is_error()) return fail<OmegaPair>(dq_res.error());

        if (dq_res.value() == 1U) {
            return fail<OmegaPair>(kv_unimpl(
                "Kovacic Case 1: r has polynomial part of degree 1 (e.g. Airy y''=xy). "
                "Case 1 inapplicable; solution is non-Liouvillian (Airy functions). "
                "Case 2/3 not applicable either."));
        }
        if (dq_res.value() >= 2U && (dq_res.value() % 2U == 1U)) {
            return fail<OmegaPair>(kv_unimpl(
                "Kovacic Case 1: r has polynomial part of odd degree " +
                std::to_string(dq_res.value()) +
                " — Case 1 inapplicable."));
        }
        
        if (dq_res.value() >= 2U) {
            unsigned v = dq_res.value() / 2U;
            Symbol y(ctx.make_fresh_symbol("y"));
            
            auto num_rev_res = reverse_polynomial(num, x, y, a, ctx);
            auto den_rev_res = reverse_polynomial(den, x, y, a, ctx);
            if (num_rev_res.is_error() || den_rev_res.is_error()) {
                return fail<OmegaPair>(kv_unimpl("Kovacic Case 1: reversing polynomial failed."));
            }
            
            auto u_coeffs_opt = compute_taylor_rational(num_rev_res.value(), den_rev_res.value(), y, kv_int(a, 0), v + 2, ctx);
            if (!u_coeffs_opt) {
                return fail<OmegaPair>(kv_unimpl("Kovacic Case 1: Taylor coefficient calculation at infinity failed."));
            }
            
            auto s_coeffs_opt = compute_laurent_sqrt(*u_coeffs_opt, v + 2);
            if (!s_coeffs_opt) {
                return fail<OmegaPair>(kv_unimpl("Kovacic Case 1: Laurent coefficient calculation at infinity failed (sqrt is not rational)."));
            }
            
            poly_op = build_laurent_inf_poly(*s_coeffs_opt, v, true);
            poly_om = build_laurent_inf_poly(*s_coeffs_opt, v, false);
        } else {
            auto qs = ctx.simplify(quotient);
            if (qs.is_error()) return fail<OmegaPair>(qs.error());
            auto c_opt = as_rational(qs.value());
            if (!c_opt) {
                return fail<OmegaPair>(kv_unimpl(
                    "Kovacic Case 1: polynomial part of r is not a rational literal."));
            }
            auto sc = rational_sqrt(*c_opt);
            if (!sc) {
                return fail<OmegaPair>(kv_unimpl(
                    "Kovacic Case 1: √c ∉ ℚ for the polynomial part of r."));
            }
            ExprPtr sc_expr = sc->is_integer()
                ? static_cast<ExprPtr>(a.make<IntegerLit>(sc->numerator()))
                : static_cast<ExprPtr>(a.make<RationalLit>(sc->numerator(), sc->denominator()));
            poly_op = sc_expr;
            poly_om = kv_neg(a, sc_expr);
        }

        auto rp_s = ctx.simplify(kv_div(a, remainder, den));
        r_proper = rp_s.is_ok() ? rp_s.value() : kv_div(a, remainder, den);
    }

    if (kv_is_zero(r_proper, ctx)) {
        return ok(OmegaPair{poly_op, poly_om});
    }

    auto pf_res = algebra::partial_fractions(r_proper, x, ctx);
    if (pf_res.is_error()) {
        return fail<OmegaPair>(kv_unimpl(
            "Kovacic Case 1: partial_fractions failed: " + pf_res.error().message));
    }

    auto r_prop_parts = algebra::apart_num_den(r_proper, ctx);
    if (r_prop_parts.is_error()) return fail<OmegaPair>(r_prop_parts.error());
    ExprPtr num_prop = r_prop_parts.value().numerator;
    ExprPtr den_prop = r_prop_parts.value().denominator;

    // apart_num_den returns num/den over a common denominator without reducing
    // to lowest terms (e.g. x⁻⁴−2x⁻³ → (x³−2x⁴)/x⁷ rather than (1−2x)/x⁴).
    // The pole-expansion below needs den_prop to carry each pole at exactly its
    // multiplicity, so cancel the polynomial gcd first.
    {
        auto g_res = algebra::polynomial_gcd(num_prop, den_prop, x, ctx);
        if (g_res.is_ok()) {
            auto dg_res = algebra::polynomial_degree(g_res.value(), x, ctx);
            if (dg_res.is_ok() && dg_res.value() >= 1U) {
                auto nq = algebra::polynomial_divmod(num_prop, g_res.value(), x, ctx);
                auto dq = algebra::polynomial_divmod(den_prop, g_res.value(), x, ctx);
                if (nq.is_ok() && dq.is_ok()) {
                    num_prop = nq.value().quotient;
                    den_prop = dq.value().quotient;
                }
            }
        }
    }

    std::vector<PFPole> collected_poles;
    for (ExprPtr term : pf_res.value()) {
        auto pole_opt = try_extract_pole(term, x, ctx);
        if (!pole_opt) {
            continue;
        }
        
        bool found = false;
        for (auto& cp : collected_poles) {
            auto diff_res = ctx.simplify(kv_sub(a, cp.pole, pole_opt->pole));
            if (diff_res.is_ok() && kv_is_zero(diff_res.value(), ctx)) {
                if (pole_opt->power > cp.power) {
                    cp.power = pole_opt->power;
                }
                found = true;
                break;
            }
        }
        if (!found) {
            collected_poles.push_back(*pole_opt);
        }
    }

    std::vector<ExprPtr> terms_plus  = {poly_op};
    std::vector<ExprPtr> terms_minus = {poly_om};

    auto build_laurent_pole_term = [&](const std::vector<Rational>& s, unsigned v, ExprPtr pole, bool positive_branch) -> ExprPtr {
        std::vector<ExprPtr> terms;
        ExprPtr xmc = kv_is_zero(pole, ctx)
            ? static_cast<ExprPtr>(a.make<Symbol>(x.name))
            : kv_sub(a, a.make<Symbol>(x.name), pole);
            
        for (unsigned j = 2; j <= v; ++j) {
            Rational coeff = s[v - j];
            if (!positive_branch) coeff = coeff * Rational(BigInt(-1));
            if (coeff.numerator().is_zero()) continue;
            
            ExprPtr coeff_expr = coeff.is_integer()
                ? static_cast<ExprPtr>(a.make<IntegerLit>(coeff.numerator()))
                : static_cast<ExprPtr>(a.make<RationalLit>(coeff.numerator(), coeff.denominator()));
                
            ExprPtr denom = a.make<Binary>(BinaryOp::Pow, xmc, kv_int(a, j));
            terms.push_back(kv_div(a, coeff_expr, denom));
        }
        
        Rational alpha = s[v - 1];
        if (!positive_branch) alpha = alpha * Rational(BigInt(-1));
        alpha = alpha + Rational(BigInt(v), BigInt(2));
        
        if (!alpha.numerator().is_zero()) {
            ExprPtr alpha_expr = alpha.is_integer()
                ? static_cast<ExprPtr>(a.make<IntegerLit>(alpha.numerator()))
                : static_cast<ExprPtr>(a.make<RationalLit>(alpha.numerator(), alpha.denominator()));
            terms.push_back(kv_div(a, alpha_expr, xmc));
        }
        
        if (terms.empty()) return kv_int(a, 0);
        ExprPtr res = terms[0];
        for (std::size_t i = 1; i < terms.size(); ++i) {
            res = kv_add(a, res, terms[i]);
        }
        return res;
    };

    for (const auto& cp : collected_poles) {
        if (cp.power % 2U == 1U) {
            return fail<OmegaPair>(kv_unimpl(
                "Kovacic Case 1: r has an odd-order pole of order " +
                std::to_string(cp.power) +
                " — Case 1 inapplicable."));
        }
        
        if (cp.power == 2U) {
            ExprPtr A_expr = nullptr;
            for (ExprPtr term : pf_res.value()) {
                auto pole_opt = try_extract_pole(term, x, ctx);
                if (pole_opt && pole_opt->power == 2U) {
                    auto diff_res = ctx.simplify(kv_sub(a, pole_opt->pole, cp.pole));
                    if (diff_res.is_ok() && kv_is_zero(diff_res.value(), ctx)) {
                        A_expr = pole_opt->coeff;
                        break;
                    }
                }
            }
            if (!A_expr) {
                continue;
            }
            
            auto as = ctx.simplify(A_expr);
            if (as.is_error()) return fail<OmegaPair>(as.error());
            auto A_opt = as_rational(as.value());
            if (!A_opt) {
                return fail<OmegaPair>(kv_unimpl(
                    "Kovacic Case 1: order-2 pole coeff is not a rational literal."));
            }

            Rational one(BigInt(1), BigInt(1));
            Rational four(BigInt(4));
            Rational disc = one + four * (*A_opt);
            auto sd = rational_sqrt(disc);
            if (!sd) {
                return fail<OmegaPair>(kv_unimpl(
                    "Kovacic Case 1: √(1+4A) ∉ ℚ — Case 2/3 needed."));
            }

            Rational half(BigInt(1), BigInt(2));
            Rational kp = half + half * (*sd);
            Rational km = half - half * (*sd);

            ExprPtr xmc = kv_is_zero(cp.pole, ctx)
                ? static_cast<ExprPtr>(a.make<Symbol>(x.name))
                : kv_sub(a, a.make<Symbol>(x.name), cp.pole);

            auto make_term = [&](const Rational& k) -> ExprPtr {
                if (k.numerator().is_zero()) return kv_int(a, 0);
                ExprPtr ke = k.is_integer()
                    ? static_cast<ExprPtr>(a.make<IntegerLit>(k.numerator()))
                    : static_cast<ExprPtr>(a.make<RationalLit>(k.numerator(), k.denominator()));
                return kv_div(a, ke, xmc);
            };

            terms_plus.push_back(make_term(kp));
            terms_minus.push_back(make_term(km));
        } else {
            unsigned v = cp.power / 2U;
            ExprPtr base = kv_sub(a, a.make<Symbol>(x.name), cp.pole);
            ExprPtr base_pow = a.make<Binary>(BinaryOp::Pow, base, kv_int(a, cp.power));
            auto div_res = algebra::polynomial_divmod(den_prop, base_pow, x, ctx);
            if (div_res.is_error()) {
                return fail<OmegaPair>(kv_unimpl("Kovacic Case 1: polynomial division failed."));
            }
            ExprPtr d_other = div_res.value().quotient;

            auto u_coeffs_opt = compute_taylor_rational(num_prop, d_other, x, cp.pole, v, ctx);
            if (!u_coeffs_opt) {
                return fail<OmegaPair>(kv_unimpl("Kovacic Case 1: Taylor coefficient calculation failed."));
            }
            
            auto s_coeffs_opt = compute_laurent_sqrt(*u_coeffs_opt, v);
            if (!s_coeffs_opt) {
                return fail<OmegaPair>(kv_unimpl("Kovacic Case 1: Laurent coefficient calculation failed (sqrt is not rational)."));
            }
            
            terms_plus.push_back(build_laurent_pole_term(*s_coeffs_opt, v, cp.pole, true));
            terms_minus.push_back(build_laurent_pole_term(*s_coeffs_opt, v, cp.pole, false));
        }
    }

    auto build_sum = [&](std::vector<ExprPtr>& ts) -> Result<ExprPtr> {
        std::vector<ExprPtr> nz;
        for (auto& t : ts) if (!kv_is_zero(t, ctx)) nz.push_back(t);
        if (nz.empty()) return ok(kv_int(a, 0));
        if (nz.size() == 1) return ok(nz[0]);
        ExprPtr s = nz[0];
        for (std::size_t i = 1; i < nz.size(); ++i) s = kv_add(a, s, nz[i]);
        return ctx.simplify(s);
    };

    auto op = build_sum(terms_plus);
    auto om = build_sum(terms_minus);
    if (op.is_error()) return fail<OmegaPair>(op.error());
    if (om.is_error()) return fail<OmegaPair>(om.error());

    return ok(OmegaPair{op.value(), om.value()});
}

} // namespace kovacic_impl
} // namespace cas::calculus
