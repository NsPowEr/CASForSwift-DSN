#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "cas/algebra.hpp"
#include "calculus_internal.hpp"
#include <vector>

namespace cas::calculus {

[[nodiscard]] Result<std::vector<Asymptote>> find_asymptotes(ExprPtr f, const Symbol& x, symbolic::CASContext& ctx) {
    std::vector<Asymptote> results;

    auto make_pos_inf = [&]() -> ExprPtr {
        return ctx.arena().make<Constant>(MathConstant::Infinity);
    };
    auto make_neg_inf = [&]() -> ExprPtr {
        return ctx.arena().make<Unary>(UnaryOp::Neg, ctx.arena().make<Constant>(MathConstant::Infinity));
    };
    auto is_inf = [](ExprPtr e) {
        return expr_is<Constant>(e) && expr_cast<Constant>(e)->value == MathConstant::Infinity;
    };
    auto is_neg_inf = [](ExprPtr e) {
        const auto* u = expr_cast<Unary>(e);
        return u && u->op == UnaryOp::Neg &&
               expr_is<Constant>(u->operand) &&
               expr_cast<Constant>(u->operand)->value == MathConstant::Infinity;
    };
    auto is_nan = [](ExprPtr e) {
        return expr_is<Constant>(e) && expr_cast<Constant>(e)->value == MathConstant::NaN;
    };
    auto is_zero_int = [](ExprPtr e) {
        return expr_is<IntegerLit>(e) && expr_cast<IntegerLit>(e)->value.is_zero();
    };

    auto try_oblique = [&](ExprPtr pt) {
        auto f_div_x = ctx.arena().make<Binary>(BinaryOp::Div, f, ctx.arena().make<Symbol>(x.name));
        auto lm = limit(f_div_x, x, pt, LimitDirection::Both, ctx);
        if (!lm.is_ok()) return;
        auto m = lm.value();
        if (is_zero_int(m) || is_inf(m) || is_neg_inf(m) || is_nan(m)) return;
        auto mx = ctx.arena().make<Binary>(BinaryOp::Mul, m, ctx.arena().make<Symbol>(x.name));
        auto lq = limit(ctx.arena().make<Binary>(BinaryOp::Sub, f, mx), x, pt, LimitDirection::Both, ctx);
        if (!lq.is_ok()) return;
        auto q = lq.value();
        if (is_inf(q) || is_neg_inf(q) || is_nan(q)) return;
        ExprPtr slant = ctx.arena().make<Binary>(BinaryOp::Add, mx, q);
        results.push_back({Asymptote::Type::Slant, slant});
    };

    // 1a. Orizzontali / Obliqui (x -> +inf)
    auto limit_pos = limit(f, x, make_pos_inf(), LimitDirection::Both, ctx);
    if (limit_pos.is_ok()) {
        auto val = limit_pos.value();
        if (is_inf(val) || is_neg_inf(val)) {
            try_oblique(make_pos_inf());
        } else if (!is_nan(val)) {
            results.push_back({Asymptote::Type::Horizontal, val});
        }
    }

    // 1b. Orizzontali / Obliqui (x -> -inf)
    auto limit_neg = limit(f, x, make_neg_inf(), LimitDirection::Both, ctx);
    if (limit_neg.is_ok()) {
        auto val = limit_neg.value();
        if (is_inf(val) || is_neg_inf(val)) {
            try_oblique(make_neg_inf());
        } else if (!is_nan(val)) {
            bool already = false;
            for (const auto& a : results) {
                if (a.type == Asymptote::Type::Horizontal && structural_equal(a.expression, val)) {
                    already = true; break;
                }
            }
            if (!already) results.push_back({Asymptote::Type::Horizontal, val});
        }
    }

    // 2. Verticali (poli: denominatore = 0)
    auto qv = extract_quotient_view(f, ctx.arena());
    if (qv.has_value()) {
        auto roots = algebra::solve_polynomial(qv->denominator, x, ctx);
        if (roots.is_ok()) {
            for (auto root : roots.value()) {
                results.push_back({Asymptote::Type::Vertical, root});
            }
        }
    }

    return ok(results);
}

} // namespace cas::calculus
