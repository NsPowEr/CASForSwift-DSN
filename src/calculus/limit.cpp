#include "cas/calculus.hpp"
#include "cas/algebra.hpp"

#include "limit_internal.hpp"
#include "calculus_internal.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <utility>

namespace cas::calculus {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

}  // namespace

Result<ExprPtr> LimitEngine::compute(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir) {
    const unsigned int tower_h = transcendental_tower_depth(expr, var);
    max_depth_budget_ = std::max<unsigned int>(8U, 2U * tower_h + 4U);
    auto together_res = algebra::together(expr, context_);
    auto simplified_expr = context_.simplify(together_res.is_ok() ? together_res.value() : expr);
    if (simplified_expr.is_error()) {
        return simplified_expr;
    }

    auto simplified_point = context_.simplify(point);
    if (simplified_point.is_error()) {
        return simplified_point;
    }

    const bool point_is_pos_inf = limit_is_infinity(simplified_point.value()) &&
                                  !expr_is<Unary>(simplified_point.value());
    const bool point_is_neg_inf = limit_is_infinity(simplified_point.value()) &&
                                  expr_is<Unary>(simplified_point.value());

    if (point_is_pos_inf || point_is_neg_inf) {
        auto infinite_limit = try_infinite_limit(simplified_expr.value(), var, simplified_point.value(), context_);
        if (infinite_limit.is_ok() || infinite_limit.error().kind != CASErrorKind::Unimplemented) {
            return infinite_limit;
        }
    }

    ExprPtr pre_mrv_expr = simplified_expr.value();
    if (auto cancelled = try_cancel_product_pow_inverse(pre_mrv_expr, context_);
        cancelled.has_value()) {
        auto re_simp = context_.simplify(cancelled.value());
        if (re_simp.is_ok()) pre_mrv_expr = re_simp.value();
    }

    if (point_is_pos_inf || point_is_neg_inf) {
        auto q_inf = extract_quotient_view(pre_mrv_expr, arena_);
        if (q_inf.has_value()) {
            LimitDirection inf_dir = point_is_pos_inf
                ? LimitDirection::Right : LimitDirection::Left;
            auto ll = try_log_log_limit(q_inf.value(), var,
                                             simplified_point.value(), inf_dir, 0U);
            if (ll.has_value()) {
                return ll.value();
            }
        }
    }

    if ((point_is_pos_inf || point_is_neg_inf)) {
        if (auto sum_res = try_limit_sum_termwise(
                pre_mrv_expr, var, simplified_point.value(), dir, context_);
            sum_res.has_value()) {
            return *sum_res;
        }
    }

    if (point_is_pos_inf || point_is_neg_inf) {
        auto mrv_res = compute_limit_mrv(pre_mrv_expr, var, simplified_point.value(), context_);
        if (mrv_res.is_ok()) {
            return mrv_res;
        }
        if (mrv_res.error().kind != CASErrorKind::Unimplemented) {
            return fail<ExprPtr>(mrv_res.error());
        }
    }

    if (point_is_pos_inf || point_is_neg_inf) {
        const std::string t_name = var.name + "_inv";
        const Symbol t_var(t_name);
        ExprPtr t_sym = arena_.make<Symbol>(t_name);

        ExprPtr one_over_t = arena_.make<Binary>(BinaryOp::Div, limit_make_integer(arena_, 1), t_sym);
        ExprPtr x_sub = point_is_pos_inf
            ? one_over_t
            : arena_.make<Unary>(UnaryOp::Neg, one_over_t);
        auto expr_sub = context_.substitute(pre_mrv_expr, var, x_sub);
        if (expr_sub.is_error()) return expr_sub;
        
        auto together_t = algebra::together(expr_sub.value(), context_);
        auto expr_t = context_.simplify(together_t.is_ok() ? together_t.value() : expr_sub.value());
        if (expr_t.is_error()) return expr_t;

        ExprPtr zero_pt = limit_make_integer(arena_, 0);
        auto res = compute_recursive(expr_t.value(), t_var, zero_pt, LimitDirection::Right, 0U);
        if (res.is_ok()) return res;

        const auto kind = res.error().kind;
        if (kind != CASErrorKind::Unimplemented &&
            kind != CASErrorKind::Undefined &&
            kind != CASErrorKind::DivisionByZero) {
            return res;
        }
    }

    return compute_recursive(pre_mrv_expr, var, simplified_point.value(), dir, 0U);
}

Result<ExprPtr> LimitEngine::substitute_and_simplify(ExprPtr expr, const Symbol& var, ExprPtr point) {
    auto substituted = context_.substitute(expr, var, point);
    if (substituted.is_error()) {
        return substituted;
    }
    auto res = context_.simplify(substituted.value());
    return res;
}

Result<ExprPtr> LimitEngine::compute_recursive(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    unsigned int depth) {
    if (depth >= max_depth_budget_) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::Unimplemented,
            "Il limite richiede piu' iterazioni di quelle supportate"));
    }

    auto direct = substitute_and_simplify(expr, var, point);

    // HC-F8-PENDING-20 §3.2: direction-limit table at branch-cut edges.
    // Plain substitution is edge-blind (always the principal/top edge); when
    // the argument of sqrt/ln lands on the negative real axis with a
    // var-dependent imaginary part, the side of the approach decides the
    // value. nullopt → no cut involvement, fall through unchanged.
    if (auto cut = try_branch_cut_directional(expr, var, point, dir, depth);
        cut.has_value()) {
        return cut.value();
    }

    auto quotient = extract_quotient_view(expr, arena_);
    const bool has_quotient = quotient.has_value();

    if (!limit_is_infinity(point)) {
        auto logarithmic_limit = try_logarithmic_root_limit(expr, var, point, dir, arena_);
        if (logarithmic_limit.has_value()) {
            return logarithmic_limit.value();
        }
    }

    if (limit_is_infinity(point)) {
        auto infinite_limit = try_infinite_limit(expr, var, point, context_);
        if (infinite_limit.is_ok() || infinite_limit.error().kind != CASErrorKind::Unimplemented) {
            return infinite_limit;
        }
    }

    if (direct.is_ok() && (!expr_is<Binary>(direct.value()) || direct.value() == expr)) {
        if (!expr_is<Binary>(expr) || (!limit_is_zero(direct.value()) && !limit_is_infinity(direct.value()))) {
            return direct;
        }
    }

    if (const auto* prod = expr_cast<Product>(expr)) {
        bool has_zero = false;
        bool has_infinity = false;
        bool all_bounded = true;
        std::vector<ExprPtr> zero_factors;
        std::vector<ExprPtr> inf_factors;
        std::vector<ExprPtr> non_zero_factors;

        for (auto factor : prod->factors) {
            auto lim_f = compute_recursive(factor, var, point, dir, depth + 1U);
            if (lim_f.is_ok() && limit_is_zero(lim_f.value())) {
                has_zero = true;
                zero_factors.push_back(factor);
            } else if (lim_f.is_ok() && limit_is_infinity(lim_f.value())) {
                has_infinity = true;
                all_bounded = false;
                inf_factors.push_back(factor);
            } else {
                if (!is_bounded(factor, var)) all_bounded = false;
                non_zero_factors.push_back(factor);
            }
        }
        if (has_zero && all_bounded) {
            return ok(limit_make_integer(arena_, 0));
        }
        // A41: prefer the NATURAL quotient view of the whole expression when
        // one already exists (e.g. (1/t)*sin(t) — one factor is itself a
        // fraction — already reduces via extract_quotient_view's Mul rule to
        // the clean sin(t)/t, analytic at the point) over building a fresh
        // one from reciprocals below. Reciprocating a side that is ALREADY a
        // quotient stacks a second pole on top of the first instead of
        // cancelling it — confirmed to blow up compute_lhopital_or_taylor
        // into runaway derivatives of an ever-more-singular quotient and time
        // out (x*sin(1/x), x->inf, reached via LimitEngine::compute()'s 1/t
        // substitution to a FINITE point where this exact shape appears).
        if (has_zero && has_infinity && has_quotient) {
            return compute_quotient_limit(expr, var, point, dir, depth + 1U, quotient.value(), direct);
        }
        // Only rewrite 0*infinity into a fresh quotient at a FINITE point,
        // and only when no natural quotient view exists (the case above).
        // At an infinite point this branch is normally unreachable anyway
        // (LimitEngine::compute() intercepts point-at-infinity cases earlier
        // via the MRV/Gruntz machinery); routing an infinite-point case
        // through here instead sends compute_quotient_limit's pole/sign
        // analysis (try_signed_pole_via_reciprocal) down a path it was never
        // designed for. Finite points are unaffected and get the real fix;
        // infinite points keep whatever this branch did before (effectively
        // inert, since simplify() collapses the old Pow(-1) construction back
        // to the input product, so control falls through to other
        // strategies).
        if (has_zero && has_infinity && !limit_is_infinity(point)) {
            // This rewrites a 0*infinity product into num/den, one side being
            // the reciprocal of the other, so that L'Hopital/Taylor (which
            // need an actual 0/0 or inf/inf quotient) can take over. The
            // construction used to build `Div(num, Pow(other_side,-1))` and
            // hand THAT to context_.simplify() before extracting a quotient
            // view from the result — but A/(B^-1) simplifies right back to
            // A*B (a correct, general simplification rule, and exactly the
            // double-reciprocal cancellation that undoes what this rewrite
            // was trying to build), so `extract_quotient_view` found nothing
            // and the code recursed on the SAME original product it started
            // from. Fix: build the QuotientView struct directly (numerator,
            // denominator as the two sides — no Pow(-1) wrapping to then
            // discard, no simplify() round-trip to undo it) and hand it
            // straight to compute_quotient_limit.
            auto get_priority = [](ExprPtr f) -> int {
                if (const auto* call = expr_cast<FuncCall>(f)) {
                    if (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) return 3;
                    if (call->func_id == BuiltinOp::Exp) return 1;
                }
                if (expr_cast<Symbol>(f)) return 2;
                if (const auto* bin = expr_cast<Binary>(f)) {
                    if (bin->op == BinaryOp::Pow) return 2;
                }
                return 0;
            };

            int max_zero_pri = -1;
            for (auto f : zero_factors) max_zero_pri = std::max(max_zero_pri, get_priority(f));
            int max_inf_pri = -1;
            for (auto f : inf_factors) max_inf_pri = std::max(max_inf_pri, get_priority(f));

            ExprPtr num, den;
            if (max_zero_pri >= max_inf_pri) {
                // 0/0 shape: numerator = zero side, denominator = 1/(infinity side).
                num = zero_factors.size() == 1 ? zero_factors[0] : arena_.make<Product>(std::move(zero_factors));
                for (auto factor : non_zero_factors) inf_factors.push_back(factor);
                ExprPtr inf_side = inf_factors.size() == 1 ? inf_factors[0] : arena_.make<Product>(std::move(inf_factors));
                den = arena_.make<Binary>(BinaryOp::Pow, inf_side, limit_make_integer(arena_, -1));
            } else {
                // infinity/infinity shape: numerator = infinity side, denominator = 1/(zero side).
                for (auto factor : non_zero_factors) inf_factors.push_back(factor);
                num = inf_factors.size() == 1 ? inf_factors[0] : arena_.make<Product>(std::move(inf_factors));
                ExprPtr zero_side = zero_factors.size() == 1 ? zero_factors[0] : arena_.make<Product>(std::move(zero_factors));
                den = arena_.make<Binary>(BinaryOp::Pow, zero_side, limit_make_integer(arena_, -1));
            }

            return compute_quotient_limit(expr, var, point, dir, depth + 1U,
                                          QuotientView{.numerator = num, .denominator = den}, direct);
        }
    }

    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Mul) {
            auto lim_l = compute_recursive(bin->left, var, point, dir, depth + 1U);
            auto lim_r = compute_recursive(bin->right, var, point, dir, depth + 1U);

            if (lim_l.is_ok() && lim_r.is_ok()) {
                bool l_zero = limit_is_zero(lim_l.value());
                bool r_zero = limit_is_zero(lim_r.value());
                bool l_inf = limit_is_infinity(lim_l.value());
                bool r_inf = limit_is_infinity(lim_r.value());

                if ((l_zero && is_bounded(bin->right, var)) || (r_zero && is_bounded(bin->left, var))) {
                    return ok(limit_make_integer(arena_, 0));
                }

                if ((l_zero && r_inf) || (r_zero && l_inf)) {
                    // A41: prefer the natural quotient view when one already
                    // exists — see the Product branch above for why.
                    if (has_quotient) {
                        return compute_quotient_limit(expr, var, point, dir, depth + 1U, quotient.value(), direct);
                    }
                }
                // Finite points only — see the Product branch above for why
                // an infinite point must not be routed through here.
                if (((l_zero && r_inf) || (r_zero && l_inf)) && !limit_is_infinity(point)) {
                    // Same fix as the Product 0*infinity branch above — build
                    // the QuotientView directly instead of routing through
                    // context_.simplify(), which cancels A/(B^-1) back to A*B
                    // and silently undoes the rewrite.
                    auto get_priority = [](ExprPtr f) -> int {
                        if (const auto* call = expr_cast<FuncCall>(f)) {
                            if (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) return 3;
                            if (call->func_id == BuiltinOp::Exp) return 1;
                        }
                        if (expr_cast<Symbol>(f)) return 2;
                        if (const auto* b = expr_cast<Binary>(f)) {
                            if (b->op == BinaryOp::Pow) return 2;
                        }
                        return 0;
                    };

                    ExprPtr zero_side = l_zero ? bin->left : bin->right;
                    ExprPtr inf_side = l_zero ? bin->right : bin->left;

                    ExprPtr num, den;
                    if (get_priority(zero_side) >= get_priority(inf_side)) {
                        num = zero_side;
                        den = arena_.make<Binary>(BinaryOp::Pow, inf_side, limit_make_integer(arena_, -1));
                    } else {
                        num = inf_side;
                        den = arena_.make<Binary>(BinaryOp::Pow, zero_side, limit_make_integer(arena_, -1));
                    }

                    return compute_quotient_limit(expr, var, point, dir, depth + 1U,
                                                  QuotientView{.numerator = num, .denominator = den}, direct);
                }
            }
        }
    }

    if (const auto* power = expr_cast<Binary>(expr); power && power->op == BinaryOp::Pow) {
        // A41: the old guard required BOTH base_lim and exp_lim to be
        // successfully classified before attempting exp(exponent*ln(base)),
        // and only for the three combinations spelled out explicitly
        // (1^inf, 0^0, inf^0). That missed 0^inf (e.g. sin(x)^(1/x), x->0+ —
        // not even indeterminate, the answer is plain 0) and broke down
        // whenever exp_lim itself failed to resolve — which it does for
        // (1+x)^(1/x) at x->0 BOTH-sided, because the raw exponent 1/x has no
        // bilateral limit (it diverges with opposite signs) even though the
        // exponent TIMES ln(base) does (ln(1+x) ~ x cancels the pole).
        //
        // General fix: the log-rewrite is valid whenever the BASE tends to
        // 1, 0, or infinity — regardless of what the raw exponent alone does,
        // because exp(exponent*ln(base)) reduces the whole power limit to a
        // single recursive limit of a product, and the Mul/Product branches
        // above already classify 0*infinity generally. No need to
        // pre-classify the exponent at all.
        //
        // Finite points only (same reasoning as the 0*infinity branches
        // above): at an infinite point this is normally unreachable anyway
        // (MRV intercepts it earlier), and routing an infinite-point power
        // through this rewrite risks the same expensive fallback machinery
        // that made x*sin(1/x), x->inf time out during verification.
        if (!limit_is_infinity(point)) {
        auto base_lim = compute_recursive(power->left, var, point, dir, depth + 1U);
        if (base_lim.is_ok()) {
            const bool base_is_one = limit_is_one(base_lim.value());
            const bool base_is_zero = limit_is_zero(base_lim.value());
            const bool base_is_inf = limit_is_infinity(base_lim.value());

            if (base_is_one || base_is_zero || base_is_inf) {
                auto ln_f = arena_.make<FuncCall>("ln", std::vector<ExprPtr>{power->left});
                auto inner_expr = arena_.make<Binary>(BinaryOp::Mul, power->right, ln_f);
                auto inner_lim = compute_recursive(inner_expr, var, point, dir, depth + 1U);
                if (inner_lim.is_ok()) {
                    return ok(arena_.make<FuncCall>("exp", std::vector<ExprPtr>{inner_lim.value()}));
                }
            }
        }
        }
    }

    // A41: last-resort composition continuity for a single-argument function
    // whose argument diverges or otherwise fails direct substitution (e.g.
    // atan(1/x) as x->0+/-: 1/x itself has no value AT x=0, but its LIMIT is
    // +-infinity, and atan is continuous there — simplify already knows
    // atan(+-infinity) = +-pi/2, exp(+-infinity), etc.). This is the general
    // rule lim g(h(x)) = g(lim h(x)) for g continuous at the limit point;
    // it only fires when the inner limit is itself fully determinate (no
    // error, no leftover dependency on var) and simplify manages to reduce
    // the composed value to something that no longer mentions var — so it can
    // never emit an unevaluated placeholder as if it were an answer.
    if (const auto* call = expr_cast<FuncCall>(expr);
        call && call->args.size() == 1U && depends_on(call->args.front(), var)) {
        auto inner_lim = compute_recursive(call->args.front(), var, point, dir, depth + 1U);
        if (inner_lim.is_ok() && !depends_on(inner_lim.value(), var)) {
            ExprPtr candidate = arena_.make<FuncCall>(call->name, std::vector<ExprPtr>{inner_lim.value()});
            auto candidate_simplified = context_.simplify(candidate);
            if (candidate_simplified.is_ok() && !depends_on(candidate_simplified.value(), var)) {
                return candidate_simplified;
            }
        }
    }

    if (direct.is_error() && !has_quotient) {
        return direct;
    }

    if (limit_is_zero(point) && dir == LimitDirection::Right) {
        if (const auto* power = expr_cast<Binary>(expr);
            power != nullptr &&
            power->op == BinaryOp::Pow &&
            expr_cast<Symbol>(power->left) != nullptr &&
            expr_cast<Symbol>(power->right) != nullptr &&
            expr_ref<Symbol>(power->left).name == var.name &&
            expr_ref<Symbol>(power->right).name == var.name) {
            return ok(limit_make_integer(arena_, 1));
        }
    }

    if (quotient.has_value()) {
        if (auto log_log = try_log_log_limit(quotient.value(), var, point, dir, depth)) {
            return log_log.value();
        }
    }

    if (!quotient.has_value()) {
        return direct;
    }

    return compute_quotient_limit(expr, var, point, dir, depth, quotient.value(), direct);
}

Result<ExprPtr> limit(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir, symbolic::CASContext& ctx) {
    auto value = LimitEngine(ctx).compute(expr, var, point, dir);
    if (value.is_error()) {
        return value;
    }
    return symbolic::materialize_expr(value.value(), ctx.arena());
}

}  // namespace cas::calculus
