// L2-05 closure: Laurent expansion of N(x)/D(x) around a point where both
// numerator and denominator are analytic but the denominator may vanish to
// finite order m, producing a pole of order ≤ m in the quotient.
//
// Algorithm (independent of the rational fast path used elsewhere):
//   1. Compute Taylor coefficients of N and D at `center` to order N + m_max,
//      where m_max is the requested negative-exponent budget plus the
//      positive expansion order.
//   2. Detect the leading exponent k of D as the smallest index whose
//      coefficient is non-zero.  Require N's leading exponent ≥ 0 (the
//      caller is responsible for stripping a polynomial part of N if any
//      finer factor cancellation is desired; here a Laurent expansion of
//      pole order at most k is produced when N has no zero at the centre).
//   3. Invert D's series via the geometric‑series identity
//        1/D(x) = (x − c)^{−k} · (1/c_k) · 1/(1 + u),
//        u = Σ_{i ≥ 1} (c_{k+i}/c_k) · (x − c)^i ,
//        1/(1+u) = Σ_{j = 0..N+k} (−u)^j   (truncated).
//   4. Multiply by N's series, truncate, and emit a LaurentExpansion struct.

#include "cas/calculus.hpp"
#include "cas/error.hpp"

#include "calculus_internal.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

struct GeneralSeries {
    // (exponent, coefficient) pairs, exponents may be negative.
    std::vector<std::pair<long long, ExprPtr>> terms;
    // Hard truncation cap: terms with exp >= order are unknown.  Order can
    // be negative when the leading exponent itself is below zero.
    long long order;
};

[[nodiscard]] bool is_zero_lit_expr(ExprPtr expr) {
    if (!expr) return true;
    const auto* il = expr_cast<IntegerLit>(expr);
    if (il && il->value.is_zero()) return true;
    return false;
}

[[nodiscard]] BigInt factorial_bigint(unsigned int n) {
    BigInt acc(1);
    for (unsigned int k = 2; k <= n; ++k) acc = acc * BigInt(k);
    return acc;
}

// Compute Taylor coefficients [c_0, c_1, ..., c_order] of `expr` at point.
// c_k = (d^k expr / dx^k)(point) / k!.
[[nodiscard]] Result<std::vector<ExprPtr>> taylor_coefficients(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    unsigned int order,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> coeffs;
    coeffs.reserve(static_cast<std::size_t>(order) + 1U);

    for (unsigned int k = 0; k <= order; ++k) {
        Result<ExprPtr> deriv = (k == 0U)
            ? ok(expr)
            : diff(expr, var, k, ctx);
        if (deriv.is_error()) return fail<std::vector<ExprPtr>>(deriv.error());

        auto sub = ctx.substitute(deriv.value(), var, point);
        if (sub.is_error()) return fail<std::vector<ExprPtr>>(sub.error());

        auto simplified = ctx.simplify(sub.value());
        if (simplified.is_error()) return fail<std::vector<ExprPtr>>(simplified.error());

        ExprPtr coeff = simplified.value();
        if (k >= 2U) {
            const BigInt fact = factorial_bigint(k);
            ExprPtr fact_expr = arena.make<IntegerLit>(fact);
            ExprPtr divided = arena.make<Binary>(BinaryOp::Div, coeff, fact_expr);
            auto simp = ctx.simplify(divided);
            if (simp.is_error()) return fail<std::vector<ExprPtr>>(simp.error());
            coeff = simp.value();
        }
        coeffs.push_back(coeff);
    }
    return ok(std::move(coeffs));
}

[[nodiscard]] GeneralSeries series_from_taylor(
    const std::vector<ExprPtr>& coeffs,
    long long order_cap) {
    GeneralSeries out;
    out.order = order_cap;
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        if (is_zero_lit_expr(coeffs[k])) continue;
        out.terms.emplace_back(static_cast<long long>(k), coeffs[k]);
    }
    return out;
}

[[nodiscard]] Result<GeneralSeries> series_mul(
    const GeneralSeries& a,
    const GeneralSeries& b,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    const long long order = std::min(a.order, b.order);
    std::vector<std::pair<long long, ExprPtr>> out_terms;

    for (const auto& [ea, ca] : a.terms) {
        for (const auto& [eb, cb] : b.terms) {
            const long long exp = ea + eb;
            if (exp >= order) continue;
            ExprPtr product = arena.make<Product>(std::vector<ExprPtr>{ca, cb});
            auto simp = ctx.simplify(product);
            if (simp.is_error()) return fail<GeneralSeries>(simp.error());
            if (is_zero_lit_expr(simp.value())) continue;
            bool merged = false;
            for (auto& [ee, cc] : out_terms) {
                if (ee == exp) {
                    ExprPtr sum_expr = arena.make<Sum>(std::vector<ExprPtr>{cc, simp.value()});
                    auto s = ctx.simplify(sum_expr);
                    if (s.is_error()) return fail<GeneralSeries>(s.error());
                    cc = s.value();
                    merged = true;
                    break;
                }
            }
            if (!merged) out_terms.emplace_back(exp, simp.value());
        }
    }
    std::sort(out_terms.begin(), out_terms.end(),
              [](const auto& l, const auto& r) { return l.first < r.first; });
    GeneralSeries result;
    result.order = order;
    for (auto& [e, c] : out_terms) {
        if (is_zero_lit_expr(c)) continue;
        result.terms.emplace_back(e, c);
    }
    return ok(std::move(result));
}

[[nodiscard]] Result<GeneralSeries> series_add(
    const GeneralSeries& a,
    const GeneralSeries& b,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    const long long order = std::min(a.order, b.order);
    std::vector<std::pair<long long, ExprPtr>> out_terms;
    for (const auto& [ea, ca] : a.terms) {
        if (ea >= order) continue;
        out_terms.emplace_back(ea, ca);
    }
    for (const auto& [eb, cb] : b.terms) {
        if (eb >= order) continue;
        bool merged = false;
        for (auto& [ee, cc] : out_terms) {
            if (ee == eb) {
                ExprPtr sum_expr = arena.make<Sum>(std::vector<ExprPtr>{cc, cb});
                auto s = ctx.simplify(sum_expr);
                if (s.is_error()) return fail<GeneralSeries>(s.error());
                cc = s.value();
                merged = true;
                break;
            }
        }
        if (!merged) out_terms.emplace_back(eb, cb);
    }
    std::sort(out_terms.begin(), out_terms.end(),
              [](const auto& l, const auto& r) { return l.first < r.first; });
    GeneralSeries result;
    result.order = order;
    for (auto& [e, c] : out_terms) {
        if (is_zero_lit_expr(c)) continue;
        result.terms.emplace_back(e, c);
    }
    return ok(std::move(result));
}

[[nodiscard]] Result<GeneralSeries> series_negate(
    const GeneralSeries& a,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    GeneralSeries out;
    out.order = a.order;
    for (const auto& [e, c] : a.terms) {
        ExprPtr neg = arena.make<Unary>(UnaryOp::Neg, c);
        auto simp = ctx.simplify(neg);
        if (simp.is_error()) return fail<GeneralSeries>(simp.error());
        out.terms.emplace_back(e, simp.value());
    }
    return ok(std::move(out));
}

[[nodiscard]] Result<GeneralSeries> series_invert(
    const GeneralSeries& s,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (s.terms.empty()) {
        return fail<GeneralSeries>(CASError{
            .kind = CASErrorKind::Undefined,
            .message = "Laurent general: cannot invert zero series"});
    }
    // Leading exponent k with nonzero coefficient.
    long long k = s.terms.front().first;
    ExprPtr c_k = s.terms.front().second;

    const long long target_order = s.order - k;
    if (target_order <= 0) {
        ExprPtr inv_ck_raw = arena.make<Binary>(BinaryOp::Div,
            arena.make<IntegerLit>(BigInt(1)), c_k);
        auto inv_ck = ctx.simplify(inv_ck_raw);
        if (inv_ck.is_error()) return fail<GeneralSeries>(inv_ck.error());
        GeneralSeries out;
        out.order = -k + 1;
        out.terms.emplace_back(-k, inv_ck.value());
        return ok(std::move(out));
    }

    // u = Σ_{e > k} (c_e / c_k) · (x − pt)^{e − k}
    GeneralSeries u;
    u.order = target_order;
    for (const auto& [exp, coeff] : s.terms) {
        if (exp == k) continue;
        const long long rel = exp - k;
        if (rel >= target_order) continue;
        ExprPtr ratio = arena.make<Binary>(BinaryOp::Div, coeff, c_k);
        auto rs = ctx.simplify(ratio);
        if (rs.is_error()) return fail<GeneralSeries>(rs.error());
        u.terms.emplace_back(rel, rs.value());
    }
    std::sort(u.terms.begin(), u.terms.end(),
              [](const auto& l, const auto& r) { return l.first < r.first; });

    auto neg_u_res = series_negate(u, ctx);
    if (neg_u_res.is_error()) return neg_u_res;
    GeneralSeries neg_u = neg_u_res.value();

    // Iterative geometric series: result = Σ_{j ≥ 0} (−u)^j up to target_order.
    GeneralSeries accumulator;
    accumulator.order = target_order;
    accumulator.terms.emplace_back(0, arena.make<IntegerLit>(BigInt(1)));
    GeneralSeries power = accumulator;

    for (long long step = 1; step <= target_order; ++step) {
        auto next_power = series_mul(power, neg_u, ctx);
        if (next_power.is_error()) return next_power;
        if (next_power.value().terms.empty()) break;
        auto next_acc = series_add(accumulator, next_power.value(), ctx);
        if (next_acc.is_error()) return next_acc;
        accumulator = next_acc.value();
        power = next_power.value();
    }

    // 1/s = (1/c_k) · (x − pt)^{−k} · accumulator.
    ExprPtr inv_ck_raw = arena.make<Binary>(BinaryOp::Div,
        arena.make<IntegerLit>(BigInt(1)), c_k);
    auto inv_ck = ctx.simplify(inv_ck_raw);
    if (inv_ck.is_error()) return fail<GeneralSeries>(inv_ck.error());

    GeneralSeries result;
    result.order = accumulator.order - k;
    for (const auto& [exp, coeff] : accumulator.terms) {
        if (is_zero_lit_expr(coeff)) continue;
        ExprPtr scaled = arena.make<Product>(std::vector<ExprPtr>{coeff, inv_ck.value()});
        auto simp = ctx.simplify(scaled);
        if (simp.is_error()) return fail<GeneralSeries>(simp.error());
        if (is_zero_lit_expr(simp.value())) continue;
        result.terms.emplace_back(exp - k, simp.value());
    }
    return ok(std::move(result));
}

// Split expr into numerator and denominator candidates without invoking
// algebra::apart_num_den (which assumes a polynomial denominator).  Walks
// Binary(Div) and Product nesting, leaving non‑quotient expressions as
// N/1.  Conservative: factors out only top‑level Div nodes.
struct QuotientPair {
    ExprPtr numerator;
    ExprPtr denominator;
};

[[nodiscard]] QuotientPair split_quotient(ExprPtr expr, AstArena& arena) {
    if (const auto* binary = expr_cast<Binary>(expr); binary && binary->op == BinaryOp::Div) {
        return {binary->left, binary->right};
    }
    return {expr, arena.make<IntegerLit>(BigInt(1))};
}

}  // namespace

[[nodiscard]] Result<LaurentExpansion> laurent_series_general(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr center,
    unsigned int positive_order,
    unsigned int pole_budget,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    QuotientPair parts = split_quotient(expr, arena);
    const unsigned int taylor_order = positive_order + pole_budget;

    auto num_coeffs = taylor_coefficients(parts.numerator, var, center, taylor_order, ctx);
    if (num_coeffs.is_error()) return fail<LaurentExpansion>(num_coeffs.error());
    auto den_coeffs = taylor_coefficients(parts.denominator, var, center, taylor_order, ctx);
    if (den_coeffs.is_error()) return fail<LaurentExpansion>(den_coeffs.error());

    GeneralSeries num_series = series_from_taylor(num_coeffs.value(),
        static_cast<long long>(taylor_order) + 1);
    GeneralSeries den_series = series_from_taylor(den_coeffs.value(),
        static_cast<long long>(taylor_order) + 1);

    if (den_series.terms.empty()) {
        return fail<LaurentExpansion>(CASError{
            .kind = CASErrorKind::Undefined,
            .message = "Laurent general: denominator is identically zero at the expansion point"});
    }

    auto inv_den = series_invert(den_series, ctx);
    if (inv_den.is_error()) return fail<LaurentExpansion>(inv_den.error());

    auto quotient_series = series_mul(num_series, inv_den.value(), ctx);
    if (quotient_series.is_error()) return fail<LaurentExpansion>(quotient_series.error());

    const auto& qs = quotient_series.value();
    if (qs.terms.empty()) {
        // Quotient is identically zero (up to truncation).  Report a flat
        // Laurent at leading_order = 0 with no coefficients.
        LaurentExpansion exp_out;
        exp_out.center = center;
        exp_out.leading_order = 0;
        exp_out.coefficients = {};
        exp_out.positive_order = positive_order;
        exp_out.remainder = arena.make<IntegerLit>(BigInt(0));
        return ok(std::move(exp_out));
    }

    const long long leading_exp = qs.terms.front().first;
    const long long upper_exp = static_cast<long long>(positive_order);

    std::vector<ExprPtr> coefficients;
    for (long long e = leading_exp; e <= upper_exp; ++e) {
        ExprPtr c = arena.make<IntegerLit>(BigInt(0));
        for (const auto& [te, tc] : qs.terms) {
            if (te == e) { c = tc; break; }
        }
        coefficients.push_back(c);
    }

    LaurentExpansion out;
    out.center = center;
    out.leading_order = static_cast<int>(leading_exp);
    out.coefficients = std::move(coefficients);
    out.positive_order = positive_order;
    out.remainder = arena.make<IntegerLit>(BigInt(0));
    return ok(std::move(out));
}

}  // namespace cas::calculus
