#pragma once

#include "cas/error.hpp"
#include "cas/symbolic.hpp"

#include <type_traits>
#include <vector>

namespace cas {
namespace algebra {

namespace resultant_detail {

[[nodiscard]] inline CASError make_resultant_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}

template <typename Coeff>
[[nodiscard]] bool coeff_is_zero(const Coeff& value, symbolic::CASContext* /*ctx*/) {
    return CoeffOps<Coeff>::is_zero(value);
}

template <>
[[nodiscard]] inline bool coeff_is_zero<ExprPtr>(const ExprPtr& value, symbolic::CASContext* /*ctx*/) {
    return !value || poly_is_zero_expr(value);
}

template <typename Coeff>
[[nodiscard]] Coeff coeff_zero_like(const Coeff& sample, symbolic::CASContext* /*ctx*/) {
    return CoeffOps<Coeff>::zero_like(sample);
}

template <>
[[nodiscard]] inline ExprPtr coeff_zero_like<ExprPtr>(const ExprPtr& /*sample*/, symbolic::CASContext* ctx) {
    return poly_make_integer(ctx->arena(), 0);
}

template <typename Coeff>
[[nodiscard]] Coeff coeff_one_like(const Coeff& sample, symbolic::CASContext* /*ctx*/) {
    return CoeffOps<Coeff>::one_like(sample);
}

template <>
[[nodiscard]] inline ExprPtr coeff_one_like<ExprPtr>(const ExprPtr& /*sample*/, symbolic::CASContext* ctx) {
    return poly_make_integer(ctx->arena(), 1);
}

template <typename Coeff>
[[nodiscard]] Result<Coeff> coeff_negate(const Coeff& value, symbolic::CASContext* /*ctx*/) {
    return ok(-value);
}

template <>
[[nodiscard]] inline Result<ExprPtr> coeff_negate<ExprPtr>(const ExprPtr& value, symbolic::CASContext* ctx) {
    return poly_simplify_expr(ctx->arena().make<Unary>(UnaryOp::Neg, value), *ctx);
}

template <typename Coeff>
[[nodiscard]] Result<Coeff> coeff_add(const Coeff& lhs, const Coeff& rhs, symbolic::CASContext* /*ctx*/) {
    return ok(lhs + rhs);
}

template <>
[[nodiscard]] inline Result<ExprPtr> coeff_add<ExprPtr>(
    const ExprPtr& lhs,
    const ExprPtr& rhs,
    symbolic::CASContext* ctx) {
    return poly_simplify_expr(ctx->arena().make<Binary>(BinaryOp::Add, lhs, rhs), *ctx);
}

template <typename Coeff>
[[nodiscard]] Result<Coeff> coeff_sub(const Coeff& lhs, const Coeff& rhs, symbolic::CASContext* /*ctx*/) {
    return ok(lhs - rhs);
}

template <>
[[nodiscard]] inline Result<ExprPtr> coeff_sub<ExprPtr>(
    const ExprPtr& lhs,
    const ExprPtr& rhs,
    symbolic::CASContext* ctx) {
    return poly_simplify_expr(ctx->arena().make<Binary>(BinaryOp::Sub, lhs, rhs), *ctx);
}

template <typename Coeff>
[[nodiscard]] Result<Coeff> coeff_mul(const Coeff& lhs, const Coeff& rhs, symbolic::CASContext* /*ctx*/) {
    return ok(lhs * rhs);
}

template <>
[[nodiscard]] inline Result<ExprPtr> coeff_mul<ExprPtr>(
    const ExprPtr& lhs,
    const ExprPtr& rhs,
    symbolic::CASContext* ctx) {
    return poly_simplify_expr(ctx->arena().make<Binary>(BinaryOp::Mul, lhs, rhs), *ctx);
}

template <typename Coeff>
[[nodiscard]] Result<Coeff> coeff_div(const Coeff& lhs, const Coeff& rhs, symbolic::CASContext* /*ctx*/) {
    auto inv = CoeffOps<Coeff>::inverse(rhs);
    if (inv.is_error()) return fail<Coeff>(inv.error());
    return ok(lhs * inv.value());
}

template <>
[[nodiscard]] inline Result<ExprPtr> coeff_div<ExprPtr>(
    const ExprPtr& lhs,
    const ExprPtr& rhs,
    symbolic::CASContext* ctx) {
    return poly_simplify_expr(ctx->arena().make<Binary>(BinaryOp::Div, lhs, rhs), *ctx);
}

template <typename Coeff>
[[nodiscard]] Result<Coeff> coeff_pow(Coeff base, std::size_t exponent, symbolic::CASContext* ctx) {
    const Coeff one = coeff_one_like(base, ctx);
    Coeff result = one;
    while (exponent > 0U) {
        if ((exponent & 1U) != 0U) {
            auto mul = coeff_mul(result, base, ctx);
            if (mul.is_error()) return fail<Coeff>(mul.error());
            result = mul.value();
        }
        exponent >>= 1U;
        if (exponent == 0U) break;
        auto square = coeff_mul(base, base, ctx);
        if (square.is_error()) return fail<Coeff>(square.error());
        base = square.value();
    }
    return ok(std::move(result));
}

template <typename Coeff>
inline void strip_trailing(std::vector<Coeff>& poly, symbolic::CASContext* ctx) {
    while (!poly.empty() && coeff_is_zero(poly.back(), ctx)) poly.pop_back();
}

template <typename Coeff>
[[nodiscard]] inline std::size_t degree(const std::vector<Coeff>& poly) {
    return poly.empty() ? 0U : poly.size() - 1U;
}

template <typename Coeff>
[[nodiscard]] inline bool is_zero_poly(const std::vector<Coeff>& poly, symbolic::CASContext* ctx) {
    if (poly.empty()) return true;
    for (const auto& coeff : poly) {
        if (!coeff_is_zero(coeff, ctx)) return false;
    }
    return true;
}

template <typename Coeff>
[[nodiscard]] Result<std::vector<Coeff>> multiply_by_scalar(
    const std::vector<Coeff>& poly,
    const Coeff& scalar,
    symbolic::CASContext* ctx) {
    std::vector<Coeff> out;
    out.reserve(poly.size());
    for (const auto& coeff : poly) {
        auto scaled = coeff_mul(coeff, scalar, ctx);
        if (scaled.is_error()) return fail<std::vector<Coeff>>(scaled.error());
        out.push_back(std::move(scaled.value()));
    }
    strip_trailing(out, ctx);
    return ok(std::move(out));
}

}  // namespace resultant_detail

template <typename Coeff>
[[nodiscard]] Result<std::vector<Coeff>> pseudo_remainder_generic(
    std::vector<Coeff> a,
    const std::vector<Coeff>& b,
    symbolic::CASContext* ctx) {
    using namespace resultant_detail;
    if constexpr (std::is_same_v<Coeff, ExprPtr>) {
        if (!ctx) {
            return fail<std::vector<Coeff>>(make_resultant_error(
                CASErrorKind::InvalidArgument,
                "Expr pseudo-remainder requires CASContext"));
        }
    }
    strip_trailing(a, ctx);
    std::vector<Coeff> divisor = b;
    strip_trailing(divisor, ctx);
    if (divisor.empty()) {
        return fail<std::vector<Coeff>>(make_resultant_error(
            CASErrorKind::InvalidArgument,
            "Generic pseudo-remainder: divisor cannot be zero"));
    }

    const std::size_t m = degree(a);
    const std::size_t n = degree(divisor);
    if (a.empty() || m < n) return ok(std::move(a));

    const Coeff& lead_b = divisor.back();
    std::vector<Coeff> remainder = std::move(a);
    for (std::size_t step = 0; step <= m - n; ++step) {
        strip_trailing(remainder, ctx);
        if (remainder.empty()) break;

        const std::size_t deg_r = degree(remainder);
        if (deg_r == m - step) {
            const Coeff lead_r = remainder.back();
            auto scaled = multiply_by_scalar(remainder, lead_b, ctx);
            if (scaled.is_error()) return fail<std::vector<Coeff>>(scaled.error());
            remainder = std::move(scaled.value());

            const std::size_t shift = deg_r - n;
            if (remainder.size() < divisor.size() + shift) {
                const Coeff zero = coeff_zero_like(divisor.front(), ctx);
                remainder.resize(divisor.size() + shift, zero);
            }
            for (std::size_t i = 0; i < divisor.size(); ++i) {
                auto term = coeff_mul(divisor[i], lead_r, ctx);
                if (term.is_error()) return fail<std::vector<Coeff>>(term.error());
                auto next = coeff_sub(remainder[shift + i], term.value(), ctx);
                if (next.is_error()) return fail<std::vector<Coeff>>(next.error());
                remainder[shift + i] = std::move(next.value());
            }
        } else {
            auto scaled = multiply_by_scalar(remainder, lead_b, ctx);
            if (scaled.is_error()) return fail<std::vector<Coeff>>(scaled.error());
            remainder = std::move(scaled.value());
        }
    }
    strip_trailing(remainder, ctx);
    return ok(std::move(remainder));
}

// Subresultant algorithm of Collins/Brown (Bronstein, Symbolic_Integration_I.md
// :1044-1071). The gamma/beta recursion of spec lines 1019-1024 is tracked in
// `c`/`b_scale` below; because `c` is carried across steps, the final `s_last`
// already carries the tau_k correction factor of Theorem 1.5.3 (spec line 1035),
// so this routine is correct on *defective* sequences too — the ones where
// deg(R_{k-1}) > 1 and "resultant = last chain element" silently under-counts.
// (Verified against sympy/Maxima on 21 cases, defective included:
//  scripts/a45_prs_simulation.py.)
template <typename Coeff>
[[nodiscard]] Result<Coeff> resultant_generic(
    std::vector<Coeff> a,
    std::vector<Coeff> b,
    symbolic::CASContext* ctx,
    const ResultantDeadline& deadline,
    std::vector<std::vector<Coeff>>* chain_out) {
    using namespace resultant_detail;
    if (chain_out != nullptr) chain_out->clear();
    auto deadline_check = [&]() -> bool {
        return deadline.has_value() && std::chrono::steady_clock::now() >= *deadline;
    };
    if constexpr (std::is_same_v<Coeff, ExprPtr>) {
        if (!ctx) {
            return fail<Coeff>(make_resultant_error(
                CASErrorKind::InvalidArgument,
                "Expr resultant requires CASContext"));
        }
    }
    strip_trailing(a, ctx);
    strip_trailing(b, ctx);
    if (a.empty() || b.empty()) {
        if constexpr (std::is_same_v<Coeff, ExprPtr>) {
            return ok(poly_make_integer(ctx->arena(), 0));
        } else {
            if (a.empty() && b.empty()) {
                return fail<Coeff>(make_resultant_error(
                    CASErrorKind::InvalidArgument,
                    "Generic resultant is undefined for two zero polynomials without coefficient context"));
            }
            const Coeff sample = !a.empty() ? a.front() : b.front();
            return ok(coeff_zero_like(sample, ctx));
        }
    }

    std::size_t n = degree(a);
    std::size_t m = degree(b);
    int sign_correction = 1;
    if (n < m) {
        std::swap(a, b);
        std::swap(n, m);
        if ((n * m) % 2U != 0U) sign_correction = -1;
    }

    const Coeff minus_one = [&]() -> Coeff {
        if constexpr (std::is_same_v<Coeff, ExprPtr>) {
            return poly_make_integer(ctx->arena(), -1);
        } else {
            Coeff one = coeff_one_like(a.front(), ctx);
            return -one;
        }
    }();

    Coeff b_scale = (((n - m + 1U) % 2U) != 0U) ? minus_one : coeff_one_like(a.front(), ctx);

    // R_0 and R_1 of the PRS. Note these are post-swap: if the caller passed
    // deg(a) < deg(b) the chain is that of (b, a), matching Theorem 1.4.1.
    if (chain_out != nullptr) {
        chain_out->push_back(a);
        chain_out->push_back(b);
    }

    auto prem = pseudo_remainder_generic(a, b, ctx);
    if (prem.is_error()) return fail<Coeff>(prem.error());
    std::vector<Coeff> h = prem.value();

    for (auto& coeff : h) {
        auto scaled = coeff_mul(coeff, b_scale, ctx);
        if (scaled.is_error()) return fail<Coeff>(scaled.error());
        coeff = std::move(scaled.value());
    }
    strip_trailing(h, ctx);
    if (chain_out != nullptr && !is_zero_poly(h, ctx)) chain_out->push_back(h);

    Coeff lc = b.back();
    auto s_last_res = [&]() -> Result<Coeff> {
        auto lc_pow = coeff_pow(lc, n - m, ctx);
        if (lc_pow.is_error()) return fail<Coeff>(lc_pow.error());
        auto neg = coeff_mul(minus_one, lc_pow.value(), ctx);
        if (neg.is_error()) return fail<Coeff>(neg.error());
        return ok(neg.value());
    }();
    if (s_last_res.is_error()) return fail<Coeff>(s_last_res.error());
    Coeff s_last = std::move(s_last_res.value());
    Coeff c = s_last;

    while (!is_zero_poly(h, ctx)) {
        if (deadline_check()) {
            return fail<Coeff>(make_resultant_error(
                CASErrorKind::Unimplemented,
                "resultant_generic: ctx.timeout() exceeded during subresultant chain"));
        }
        const std::size_t k = degree(h);
        a = std::move(b);
        b = std::move(h);
        m = k;
        const std::size_t d_next = degree(a) - k;

        auto c_pow = coeff_pow(c, d_next, ctx);
        if (c_pow.is_error()) return fail<Coeff>(c_pow.error());
        auto neg_lc = coeff_mul(minus_one, lc, ctx);
        if (neg_lc.is_error()) return fail<Coeff>(neg_lc.error());
        auto next_b = coeff_mul(neg_lc.value(), c_pow.value(), ctx);
        if (next_b.is_error()) return fail<Coeff>(next_b.error());
        b_scale = next_b.value();

        auto next_prem = pseudo_remainder_generic(a, b, ctx);
        if (next_prem.is_error()) return fail<Coeff>(next_prem.error());
        h = next_prem.value();

        for (auto& coeff : h) {
            auto scaled = coeff_div(coeff, b_scale, ctx);
            if (scaled.is_error()) return fail<Coeff>(scaled.error());
            coeff = std::move(scaled.value());
        }
        strip_trailing(h, ctx);
        if (chain_out != nullptr && !is_zero_poly(h, ctx)) chain_out->push_back(h);

        lc = b.back();
        if (d_next > 1U) {
            auto q = coeff_pow(c, d_next - 1U, ctx);
            if (q.is_error()) return fail<Coeff>(q.error());
            auto neg_lc_next = coeff_mul(minus_one, lc, ctx);
            if (neg_lc_next.is_error()) return fail<Coeff>(neg_lc_next.error());
            auto neg_lc_pow = coeff_pow(neg_lc_next.value(), d_next, ctx);
            if (neg_lc_pow.is_error()) return fail<Coeff>(neg_lc_pow.error());
            auto next_c = coeff_div(neg_lc_pow.value(), q.value(), ctx);
            if (next_c.is_error()) return fail<Coeff>(next_c.error());
            c = next_c.value();
        } else {
            auto next_c = coeff_mul(minus_one, lc, ctx);
            if (next_c.is_error()) return fail<Coeff>(next_c.error());
            c = next_c.value();
        }

        auto next_s = coeff_mul(minus_one, c, ctx);
        if (next_s.is_error()) return fail<Coeff>(next_s.error());
        s_last = next_s.value();
    }

    if (degree(b) > 0U) {
        if constexpr (std::is_same_v<Coeff, ExprPtr>) {
            return ok(poly_make_integer(ctx->arena(), 0));
        } else {
            return ok(coeff_zero_like(b.front(), ctx));
        }
    }

    if (sign_correction == -1) {
        auto neg = coeff_mul(minus_one, s_last, ctx);
        if (neg.is_error()) return fail<Coeff>(neg.error());
        return ok(neg.value());
    }
    return ok(s_last);
}

}  // namespace algebra
}  // namespace cas
