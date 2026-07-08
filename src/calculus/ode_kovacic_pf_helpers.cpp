// ode_kovacic_pf_helpers.cpp — implementation of shared PF/∞ extractors.
// See ode_kovacic_pf_helpers.hpp.

#include "ode_kovacic_pf_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"

namespace cas::calculus::kovacic_impl {

std::optional<Rational> try_get_rational(
    ExprPtr e, symbolic::CASContext& ctx) {
    auto s = ctx.simplify(e);
    if (s.is_error()) return std::nullopt;
    return as_rational(s.value());
}

std::optional<long long> rational_to_int(const Rational& r) {
    if (!r.is_integer()) return std::nullopt;
    const BigInt& n = r.numerator();
    if (n.bit_length() > 32U) return std::nullopt;
    long long v = static_cast<long long>(n.to_u64());
    if (n.is_negative()) v = -v;
    return v;
}

std::optional<ExprPtr> extract_pole_loc(
    ExprPtr base, const Symbol& x, AstArena& a, symbolic::CASContext& ctx) {
    // General linear-factor root extraction: base = c0 + c1*x  (any AST shape,
    // e.g. Binary(Sub,...) or the canonical n-ary Sum node produced by arena
    // normalisation) has pole location  c = −c0 / c1.  Parsing through
    // `parse_polynomial` avoids pattern-matching on a specific tree shape
    // (REGOLA ZERO categoria 8): any representation of a degree-1 polynomial
    // in x is recognised uniformly.
    auto poly_res = algebra::parse_polynomial(base, x, ctx);
    if (poly_res.is_error()) return std::nullopt;
    const auto& poly = poly_res.value();
    if (poly.size() == 1U) return std::nullopt;   // constant: no root in x.
    if (poly.size() != 2U) return std::nullopt;    // not a linear factor.
    ExprPtr c0 = poly[0];
    ExprPtr c1 = poly[1];
    ExprPtr loc = a.make<Binary>(BinaryOp::Div,
        a.make<Unary>(UnaryOp::Neg, c0), c1);
    auto s = ctx.simplify(loc);
    return s.is_ok() ? s.value() : loc;
}

std::optional<std::pair<ExprPtr, unsigned>> as_neg_pow(ExprPtr e) {
    auto* pw = expr_cast<Binary>(e);
    if (!pw || pw->op != BinaryOp::Pow) return std::nullopt;
    auto* el = expr_cast<IntegerLit>(pw->right);
    if (!el || !el->value.is_negative()) return std::nullopt;
    BigInt n = el->value.abs();
    if (n.bit_length() > 16U) return std::nullopt;
    unsigned u = static_cast<unsigned>(n.to_u64());
    if (u < 1U || u > 1024U) return std::nullopt;
    return std::make_pair(pw->left, u);
}

std::vector<PFPole> collect_pf_poles(
    const std::vector<ExprPtr>& pf_terms,
    const Symbol& x, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    std::vector<PFPole> out;

    auto record = [&](ExprPtr pole_loc, unsigned pw, ExprPtr coeff) {
        auto ps = ctx.simplify(pole_loc);
        if (ps.is_ok()) pole_loc = ps.value();
        auto cs = ctx.simplify(coeff);
        if (cs.is_ok()) coeff = cs.value();
        out.push_back(PFPole{pole_loc, pw, coeff});
    };

    for (ExprPtr term : pf_terms) {
        // (a) Div(N, D).
        if (auto* dv = expr_cast<Binary>(term);
            dv && dv->op == BinaryOp::Div) {
            ExprPtr den = dv->right;
            unsigned pw = 1U;
            ExprPtr base = den;
            if (const auto* pw_node = expr_cast<Binary>(den);
                pw_node && pw_node->op == BinaryOp::Pow) {
                const auto* el = expr_cast<IntegerLit>(pw_node->right);
                if (el && !el->value.is_negative() && !el->value.is_zero()) {
                    long long ev = static_cast<long long>(el->value.to_u64());
                    if (ev >= 1 && ev <= 1024) {
                        pw = static_cast<unsigned>(ev);
                        base = pw_node->left;
                    }
                }
            }
            auto loc = extract_pole_loc(base, x, a, ctx);
            if (loc) record(*loc, pw, dv->left);
            continue;
        }

        // (b) Pow(base, −k) on its own → coeff = 1.
        if (auto np = as_neg_pow(term)) {
            auto loc = extract_pole_loc(np->first, x, a, ctx);
            if (loc) record(*loc, np->second, kv_int(a, 1));
            continue;
        }

        // (c) Product with a single negative-power factor in x.
        if (auto* prod = expr_cast<Product>(term)) {
            std::vector<ExprPtr> coeffs;
            coeffs.reserve(prod->factors.size());
            std::optional<ExprPtr> pole_loc;
            unsigned pw = 0U;
            bool ambiguous = false;
            for (ExprPtr f : prod->factors) {
                auto np = as_neg_pow(f);
                if (!np) { coeffs.push_back(f); continue; }
                auto loc = extract_pole_loc(np->first, x, a, ctx);
                if (!loc) { ambiguous = true; break; }
                if (pole_loc.has_value()) { ambiguous = true; break; }
                pole_loc = *loc;
                pw = np->second;
            }
            if (ambiguous || !pole_loc.has_value()) continue;
            ExprPtr coeff = coeffs.empty()
                ? kv_int(a, 1)
                : (coeffs.size() == 1U
                    ? coeffs[0]
                    : static_cast<ExprPtr>(a.make<Product>(std::move(coeffs))));
            record(*pole_loc, pw, coeff);
            continue;
        }
    }
    return out;
}

std::optional<InfinityData> compute_infinity_data(
    ExprPtr r, const Symbol& x, symbolic::CASContext& ctx) {
    auto parts = algebra::apart_num_den(r, ctx);
    if (parts.is_error()) return std::nullopt;
    ExprPtr num = parts.value().numerator;
    ExprPtr den = parts.value().denominator;
    auto dn_res = algebra::polynomial_degree(num, x, ctx);
    auto dd_res = algebra::polynomial_degree(den, x, ctx);
    if (dn_res.is_error() || dd_res.is_error()) return std::nullopt;
    long long deg_num = static_cast<long long>(dn_res.value());
    long long deg_den = static_cast<long long>(dd_res.value());
    InfinityData info;
    info.ord = deg_den - deg_num;
    info.leading_set = false;
    // Leading Laurent coefficient at ∞ for the branch (∞₂)/(γ in §5):
    //   b = leading_coeff(num) / leading_coeff(den).
    auto nc_res = algebra::univariate_coefficients(num, x, ctx);
    auto dc_res = algebra::univariate_coefficients(den, x, ctx);
    if (nc_res.is_ok() && dc_res.is_ok()
        && !nc_res.value().empty() && !dc_res.value().empty()) {
        auto lc_num_q = try_get_rational(nc_res.value().back(), ctx);
        auto lc_den_q = try_get_rational(dc_res.value().back(), ctx);
        if (lc_num_q && lc_den_q && !lc_den_q->numerator().is_zero()) {
            info.leading_b = *lc_num_q / *lc_den_q;
            info.leading_set = true;
        }
    }
    return info;
}

} // namespace cas::calculus::kovacic_impl
