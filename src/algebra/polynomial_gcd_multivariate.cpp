// polynomial_gcd_multivariate.cpp
// Public entry point polynomial_gcd_multivariate() and the core recursive
// evaluation/interpolation GCD algorithm (gcd_multivariate_recursive +
// gcd_multivariate_eval_interp).
//
// Split companions:
//   polynomial_gcd_multivariate_helpers.cpp  — FactorKey, normalize, conversions
//   polynomial_gcd_multivariate_interp.cpp   — Lagrange + interpolate_polynomial_values
//   polynomial_gcd_multivariate_sparse.cpp   — SparsePoly arithmetic + exact_quotient
//   polynomial_gcd_multivariate_linear.cpp   — primitive_linear_candidates

#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include "polynomial_gcd_multivariate_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <vector>

namespace cas::algebra {

// Forward declaration: try_certified_linear_gcd calls gcd_multivariate_recursive
// and vice-versa; both are defined in this translation unit.
[[nodiscard]] static Result<std::optional<MultivariatePolynomial>> try_certified_linear_gcd(
    const MultivariatePolynomial& p,
    const MultivariatePolynomial& q,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx,
    std::size_t depth,
    std::size_t& call_count);

// ---------------------------------------------------------------------------
// Core recursive GCD: evaluation/interpolation over Z
// ---------------------------------------------------------------------------

[[nodiscard]] Result<MultivariatePolynomial> gcd_multivariate_recursive(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx,
    std::size_t depth,
    std::size_t& call_count) {
    // Configurable total-call budget (CLAUDE.md Cat 1 — no fixed computational
    // limit).  Bound justification: by Schwartz-Zippel, each nesting level
    // expands fan-out by at most (2*D+3) evaluations; default 4096 safely covers
    // all currently-tested inputs (≤4 vars, deg ≤2) while bounding pathological
    // inputs.
    if (++call_count > ctx.max_gcd_total_calls()) {
        return fail<MultivariatePolynomial>(make_error(
            CASErrorKind::Unimplemented,
            "polynomial_gcd_multivariate [module=algebra, function=gcd_multivariate, "
            "reason_code=GCD_MULTIVARIATE_BUDGET_EXCEEDED, ticket=L1-08/F3.1]: "
            "total call budget exceeded — increase ctx.max_gcd_total_calls()"));
    }
    if (depth > ctx.max_gcd_recursion_depth()) {
        return fail<MultivariatePolynomial>(make_error(
            CASErrorKind::Unimplemented,
            "polynomial_gcd_multivariate [module=algebra, function=gcd_multivariate, "
            "reason_code=GCD_MULTIVARIATE_BUDGET_EXCEEDED, ticket=L1-08/F3.1]: "
            "recursion depth limit reached — increase ctx.max_gcd_recursion_depth()"));
    }

    MultivariatePolynomial p = normalize_multivariate_gcd(P);
    MultivariatePolynomial q = normalize_multivariate_gcd(Q);

    if (p.is_zero()) {
        return ok(q);
    }
    if (q.is_zero()) {
        return ok(p);
    }

    const std::vector<Symbol> vars = collect_all_variables(p, q);
    if (vars.empty()) {
        BigInt gcd_value = gcd(p.terms().front().coefficient.abs(), q.terms().front().coefficient.abs());
        return ok(MultivariatePolynomial({MultivariateTerm{.coefficient = gcd_value, .factors = {}}}));
    }

    if (vars.size() == 1U) {
        auto p_int = multivariate_single_var_to_intpoly(p, vars[0]);
        if (p_int.is_error()) {
            return fail<MultivariatePolynomial>(p_int.error());
        }
        auto q_int = multivariate_single_var_to_intpoly(q, vars[0]);
        if (q_int.is_error()) {
            return fail<MultivariatePolynomial>(q_int.error());
        }

        IntPoly gcd_int = gcd_integer_poly_primitive(std::move(p_int.value()), std::move(q_int.value()));
        if (!gcd_int.empty() && gcd_int.leading_coeff().is_negative()) {
            for (auto& coefficient : gcd_int.coefficients()) {
                coefficient = -coefficient;
            }
        }
        return ok(intpoly_to_multivariate(gcd_int, vars[0]));
    }

    // Identical-polynomial shortcut: covers gcd(p,p)=p for all variable counts,
    // including linear trivariate cases where the degree≤1 guard below is absent.
    if (same_polynomial(p, q)) {
        return ok(normalize_multivariate_gcd(p));
    }

    // Linear-candidate fast-path: only worthwhile for ≤2 variables.
    // For 3+ variables the candidate count grows as 5^n which makes the scan
    // more expensive than the eval-interpolation path it is meant to accelerate.
    if (vars.size() <= 2U) {
        auto linear_gcd = try_certified_linear_gcd(p, q, vars, ctx, depth, call_count);
        if (linear_gcd.is_error()) {
            return fail<MultivariatePolynomial>(linear_gcd.error());
        }
        if (linear_gcd.value().has_value()) {
            return ok(linear_gcd.value().value());
        }
        // Bivariate termination: try_certified_linear_gcd calls gcd_multivariate_recursive
        // for cofactor GCDs, which would loop forever if we fell through to interpolation
        // on degree-1 coprime inputs (e.g., gcd(x-y, x+y) = 1).  The identical-polynomial
        // check above already handled p==q; here p≠q so returning 1 is correct.
        if (std::min(p.total_degree(), q.total_degree()) <= 1U) {
            return ok(MultivariatePolynomial({MultivariateTerm{.coefficient = BigInt(1), .factors = {}}}));
        }
    }

    const Symbol& main_var = vars.front();
    const Symbol& interpolation_var = vars.back();

    const std::size_t interpolation_degree_bound =
        std::min(degree_in_var(p, interpolation_var), degree_in_var(q, interpolation_var));
    const std::size_t required_samples = std::max<std::size_t>(interpolation_degree_bound + 1U, 2U);

    // Dynamic samples calculation via Schwartz-Zippel (ZP-1):
    // k = ceil(log(delta) / log(p_fail)) with p_fail = d / S.
    // We choose a conservative evaluation box of size S = 10000.
    const double delta = ctx.zippel_error_probability();
    const double S = 10000.0;
    const double d = std::max(1.0, static_cast<double>(interpolation_degree_bound));
    const double p_fail = std::min(0.5, d / S);
    const std::size_t extra_samples = (delta > 0.0 && delta < 1.0)
        ? static_cast<std::size_t>(std::ceil(std::log(delta) / std::log(p_fail)))
        : 8U;

    const std::size_t max_samples = required_samples + std::max<std::size_t>(2U, extra_samples);

    struct SamplePoint {
        BigInt value;
        MultivariatePolynomial gcd;
        std::size_t main_degree{0U};
    };

    const std::array<long long, 3> offsets = {0LL, 17LL, 43LL};

    for (long long offset : offsets) {
        std::vector<SamplePoint> samples;
        samples.reserve(max_samples);

        bool sampling_ok = true;
        for (std::size_t i = 0U; i < max_samples; ++i) {
            const BigInt point(static_cast<long long>(i) + 1LL + offset);
            ExprPtr value_expr = ctx.arena().make<IntegerLit>(point);

            auto p_eval = p.evaluate_at(interpolation_var, value_expr);
            if (p_eval.is_error()) {
                sampling_ok = false;
                break;
            }
            auto q_eval = q.evaluate_at(interpolation_var, value_expr);
            if (q_eval.is_error()) {
                sampling_ok = false;
                break;
            }

            auto gcd_eval = gcd_multivariate_recursive(p_eval.value(), q_eval.value(), ctx, depth + 1U, call_count);
            if (gcd_eval.is_error()) {
                sampling_ok = false;
                break;
            }

            MultivariatePolynomial normalized_eval = normalize_multivariate_gcd(gcd_eval.value());
            samples.push_back(SamplePoint{
                .value = point,
                .gcd = std::move(normalized_eval),
                .main_degree = degree_in_var(samples.empty() ? MultivariatePolynomial{} : samples.back().gcd, main_var),
            });
            samples.back().main_degree = degree_in_var(samples.back().gcd, main_var);
        }

        if (!sampling_ok) {
            continue;
        }

        std::map<std::size_t, std::vector<SamplePoint>> buckets;
        for (const auto& sample : samples) {
            buckets[sample.main_degree].push_back(sample);
        }

        std::optional<std::size_t> target_degree;
        for (const auto& [degree, bucket] : buckets) {
            if (bucket.size() >= required_samples) {
                if (!target_degree.has_value() || degree > target_degree.value()) {
                    target_degree = degree;
                }
            }
        }

        if (!target_degree.has_value()) {
            continue;
        }

        const auto& chosen_bucket = buckets[target_degree.value()];
        std::vector<BigInt> chosen_points;
        std::vector<MultivariatePolynomial> chosen_polys;
        chosen_points.reserve(required_samples);
        chosen_polys.reserve(required_samples);

        for (std::size_t i = 0U; i < required_samples; ++i) {
            chosen_points.push_back(chosen_bucket[i].value);
            chosen_polys.push_back(chosen_bucket[i].gcd);
        }

        std::vector<MultivariateTerm> candidate_terms;

        bool interpolation_ok = true;
        for (std::size_t degree = 0U; degree <= target_degree.value(); ++degree) {
            std::vector<MultivariatePolynomial> coefficient_values;
            coefficient_values.reserve(chosen_polys.size());
            for (const auto& sample_poly : chosen_polys) {
                coefficient_values.push_back(coefficient_poly_in_var(sample_poly, main_var, degree));
            }

            auto interpolation = interpolate_polynomial_values(coefficient_values, chosen_points, interpolation_var);
            if (interpolation.is_error()) {
                interpolation_ok = false;
                break;
            }

            MultivariatePolynomial lifted = multiply_by_variable_power(interpolation.value(), main_var, degree);
            for (const auto& term : lifted.terms()) {
                candidate_terms.push_back(term);
            }
        }

        if (!interpolation_ok) {
            continue;
        }

        MultivariatePolynomial candidate = normalize_multivariate_gcd(MultivariatePolynomial(std::move(candidate_terms)));
        if (candidate.is_zero()) {
            continue;
        }

        // Unit-candidate short-circuit: if interpolation produced a unit (1),
        // all sample GCDs were 1 at ≥ required_samples evaluation points.
        // By the Schwartz-Zippel uniqueness argument (required_samples ≥
        // interpolation_degree_bound + 1), the unique polynomial consistent with
        // those samples IS 1.  No recursive certification is needed — calling
        // gcd(p/1, q/1) = gcd(p, q) would be circular and causes the
        // 3^depth exponential blowup that constitutes BUG-HANG-002.
        if (is_unit_polynomial(candidate)) {
            return ok(candidate);
        }

        auto quotient_p = exact_quotient(p, candidate, vars, ctx);
        if (quotient_p.is_error() || !quotient_p.value().has_value()) {
            continue;
        }
        auto quotient_q = exact_quotient(q, candidate, vars, ctx);
        if (quotient_q.is_error() || !quotient_q.value().has_value()) {
            continue;
        }

        MultivariatePolynomial certified_gcd = candidate;
        MultivariatePolynomial cofactor_p = quotient_p.value().value();
        MultivariatePolynomial cofactor_q = quotient_q.value().value();

        bool certified = false;
        for (std::size_t refine = 0U; refine < 3U; ++refine) {
            auto cofactor_gcd = gcd_multivariate_recursive(cofactor_p, cofactor_q, ctx, depth + 1U, call_count);
            if (cofactor_gcd.is_error()) {
                certified = false;
                break;
            }

            MultivariatePolynomial extra_factor = normalize_multivariate_gcd(cofactor_gcd.value());
            if (is_unit_polynomial(extra_factor)) {
                certified = true;
                break;
            }

            MultivariatePolynomial next_gcd = normalize_multivariate_gcd(certified_gcd * extra_factor);
            auto next_quotient_p = exact_quotient(p, next_gcd, vars, ctx);
            auto next_quotient_q = exact_quotient(q, next_gcd, vars, ctx);
            if (next_quotient_p.is_error() || next_quotient_q.is_error() ||
                !next_quotient_p.value().has_value() || !next_quotient_q.value().has_value()) {
                certified = false;
                break;
            }

            if (same_polynomial(cofactor_p, next_quotient_p.value().value()) &&
                same_polynomial(cofactor_q, next_quotient_q.value().value())) {
                certified = false;
                break;
            }

            certified_gcd = std::move(next_gcd);
            cofactor_p = std::move(next_quotient_p.value().value());
            cofactor_q = std::move(next_quotient_q.value().value());
        }

        if (certified) {
            return ok(normalize_multivariate_gcd(certified_gcd));
        }
    }

    return fail<MultivariatePolynomial>(make_error(
        CASErrorKind::Unimplemented,
        "polynomial_gcd_multivariate: unable to certify multivariate gcd for this input"));
}

// ---------------------------------------------------------------------------
// Linear GCD certification (uses primitive_linear_candidates from _linear.cpp)
// ---------------------------------------------------------------------------

[[nodiscard]] static Result<std::optional<MultivariatePolynomial>> try_certified_linear_gcd(
    const MultivariatePolynomial& p,
    const MultivariatePolynomial& q,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx,
    std::size_t depth,
    std::size_t& call_count) {
    for (const auto& candidate : primitive_linear_candidates(vars)) {
        auto quotient_p = exact_quotient(p, candidate, vars, ctx);
        if (quotient_p.is_error()) {
            continue;
        }
        if (!quotient_p.value().has_value()) {
            continue;
        }

        auto quotient_q = exact_quotient(q, candidate, vars, ctx);
        if (quotient_q.is_error()) {
            continue;
        }
        if (!quotient_q.value().has_value()) {
            continue;
        }

        auto cofactor_gcd = gcd_multivariate_recursive(
            quotient_p.value().value(),
            quotient_q.value().value(),
            ctx,
            depth + 1U,
            call_count);
        if (cofactor_gcd.is_error()) {
            return fail<std::optional<MultivariatePolynomial>>(cofactor_gcd.error());
        }

        MultivariatePolynomial certified = normalize_multivariate_gcd(candidate * cofactor_gcd.value());
        return ok(std::optional<MultivariatePolynomial>{std::move(certified)});
    }

    return ok(std::optional<MultivariatePolynomial>{std::nullopt});
}

// ---------------------------------------------------------------------------
// Public wrappers (declared in algebra_internal.hpp)
// ---------------------------------------------------------------------------

Result<MultivariatePolynomial> gcd_multivariate_eval_interp(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx) {
    std::size_t call_count = 0U;
    return gcd_multivariate_recursive(P, Q, ctx, 0U, call_count);
}

Result<ExprPtr> polynomial_gcd_multivariate(ExprPtr p, ExprPtr q, symbolic::CASContext& ctx) {
    auto expanded_p = expand_expr_impl(p, ctx);
    if (expanded_p.is_error()) {
        return fail<ExprPtr>(expanded_p.error());
    }

    auto expanded_q = expand_expr_impl(q, ctx);
    if (expanded_q.is_error()) {
        return fail<ExprPtr>(expanded_q.error());
    }

    auto P = parse_multivariate_polynomial(expanded_p.value(), ctx);
    if (P.is_error()) {
        return fail<ExprPtr>(P.error());
    }

    auto Q = parse_multivariate_polynomial(expanded_q.value(), ctx);
    if (Q.is_error()) {
        return fail<ExprPtr>(Q.error());
    }

    auto G = gcd_zippel_sparse(P.value(), Q.value(), ctx);
    if (G.is_error()) {
        return fail<ExprPtr>(G.error());
    }

    auto gcd_expr = multivariate_to_expr(G.value(), ctx);
    if (gcd_expr.is_error()) {
        return fail<ExprPtr>(gcd_expr.error());
    }

    return ok(gcd_expr.value());
}

} // namespace cas::algebra
