#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <functional>
#include <utility>
#include <vector>

namespace cas::symbolic {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr) {
    if (!expr) {
        return std::nullopt;
    }

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }

    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            auto inner = exact_scalar_from_expr(unary->operand);
            if (inner.has_value()) {
                return -(*inner);
            }
        }
    }

    return std::nullopt;
}

[[nodiscard]] ExprPtr negate_expr(ExprPtr expr, AstArena& arena) {
    if (!expr) {
        return expr;
    }

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return arena.make<IntegerLit>(-integer->value);
    }

    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return arena.make<RationalLit>(-rational->numerator, rational->denominator);
    }

    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) {
        return unary->operand;
    }

    return arena.make<Unary>(UnaryOp::Neg, expr);
}

[[nodiscard]] int compare_exact_scalars(const Rational& lhs, const Rational& rhs) {
    const BigInt left_cross = lhs.numerator() * rhs.denominator();
    const BigInt right_cross = rhs.numerator() * lhs.denominator();
    if (left_cross < right_cross) {
        return -1;
    }
    if (left_cross > right_cross) {
        return 1;
    }
    return 0;
}

[[nodiscard]] std::size_t expr_weight(ExprPtr expr) {
    if (!expr) {
        return 0U;
    }

    return 1U + visit_expr(
        expr,
        [](const auto& node) -> std::size_t {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, Symbol> ||
                std::is_same_v<Node, Constant>) {
                return 0U;
            } else if constexpr (std::is_same_v<Node, Unary>) {
                return expr_weight(node.operand);
            } else if constexpr (std::is_same_v<Node, Binary>) {
                return expr_weight(node.left) + expr_weight(node.right);
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::size_t weight = 0U;
                for (ExprPtr arg : node.args) {
                    weight += expr_weight(arg);
                }
                return weight;
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::size_t weight = 0U;
                for (ExprPtr term : node.terms) {
                    weight += expr_weight(term);
                }
                return weight;
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::size_t weight = 0U;
                for (ExprPtr factor : node.factors) {
                    weight += expr_weight(factor);
                }
                return weight;
            } else if constexpr (std::is_same_v<Node, Integral>) {
                return expr_weight(node.integrand) +
                    (node.lower.has_value() ? expr_weight(*node.lower) : 0U) +
                    (node.upper.has_value() ? expr_weight(*node.upper) : 0U);
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                return expr_weight(node.expression);
            } else if constexpr (std::is_same_v<Node, Limit>) {
                return expr_weight(node.expression) + expr_weight(node.point);
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                return expr_weight(node.polynomial);
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::size_t weight = 0U;
                for (ExprPtr element : node.elements) {
                    weight += expr_weight(element);
                }
                return weight;
            } else if constexpr (std::is_same_v<Node, SeriesExp>) {
                std::size_t weight = expr_weight(node.point);
                for (const auto& [exp, coeff] : node.terms) weight += expr_weight(coeff);
                return weight;
            } else {
                return 0U;
            }
        });
}

[[nodiscard]] const ComputationTrace& empty_trace() noexcept {
    static const ComputationTrace trace;
    return trace;
}

[[nodiscard]] Result<ExprPtr> materialize_expr_impl(ExprPtr expr, AstArena& arena) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot clone null expression"));
    }

    return ok(visit_expr(
        expr,
        [&](const auto& node) -> ExprPtr {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, ComplexLit> ||
                std::is_same_v<Node, Symbol> ||
                std::is_same_v<Node, Constant>) {
                return arena.make<Node>(node);
            } else if constexpr (std::is_same_v<Node, Unary>) {
                return arena.make<Unary>(node.op, materialize_expr_impl(node.operand, arena).value());
            } else if constexpr (std::is_same_v<Node, Binary>) {
                return arena.make<Binary>(
                    node.op,
                    materialize_expr_impl(node.left, arena).value(),
                    materialize_expr_impl(node.right, arena).value());
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> args;
                args.reserve(node.args.size());
                for (ExprPtr arg : node.args) {
                    args.push_back(materialize_expr_impl(arg, arena).value());
                }
                return arena.make<FuncCall>(node.name, std::move(args));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> terms;
                terms.reserve(node.terms.size());
                for (ExprPtr term : node.terms) {
                    terms.push_back(materialize_expr_impl(term, arena).value());
                }
                return arena.make<Sum>(std::move(terms));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> factors;
                factors.reserve(node.factors.size());
                for (ExprPtr factor : node.factors) {
                    factors.push_back(materialize_expr_impl(factor, arena).value());
                }
                return arena.make<Product>(std::move(factors));
            } else if constexpr (std::is_same_v<Node, Integral>) {
                return arena.make<Integral>(
                    materialize_expr_impl(node.integrand, arena).value(),
                    node.variable,
                    node.lower.has_value()
                        ? std::optional<ExprPtr>(materialize_expr_impl(*node.lower, arena).value())
                        : std::nullopt,
                    node.upper.has_value()
                        ? std::optional<ExprPtr>(materialize_expr_impl(*node.upper, arena).value())
                        : std::nullopt);
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                return arena.make<Derivative>(
                    materialize_expr_impl(node.expression, arena).value(),
                    node.variable,
                    node.order);
            } else if constexpr (std::is_same_v<Node, Limit>) {
                return arena.make<Limit>(
                    materialize_expr_impl(node.expression, arena).value(),
                    node.variable,
                    materialize_expr_impl(node.point, arena).value(),
                    node.direction);
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                return arena.make<RootOf>(
                    materialize_expr_impl(node.polynomial, arena).value(),
                    node.variable,
                    node.root_index);
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::vector<ExprPtr> elements;
                elements.reserve(node.elements.size());
                for (ExprPtr element : node.elements) {
                    elements.push_back(materialize_expr_impl(element, arena).value());
                }
                return arena.make<Matrix>(node.rows, node.cols, std::move(elements));
            } else if constexpr (std::is_same_v<Node, SeriesExp>) {
                std::vector<std::pair<long long, ExprPtr>> terms;
                terms.reserve(node.terms.size());
                for (const auto& [exp, coeff] : node.terms) {
                    terms.push_back({exp, materialize_expr_impl(coeff, arena).value()});
                }
                return arena.make<SeriesExp>(node.var, materialize_expr_impl(node.point, arena).value(), std::move(terms), node.order);
            } else {
                return ExprPtr{};
            }
        }));
}

[[nodiscard]] ExprPtr instantiate_pattern(ExprPtr pattern, const MatchMap& matches, AstArena& arena) {
    if (!pattern) {
        return ExprPtr{};
    }

    if (const auto* wildcard = expr_cast<Symbol>(pattern)) {
        if (is_wildcard_name(wildcard->name)) {
            const auto found = matches.find(wildcard->name);
            if (found != matches.end()) {
                return found->second;
            }
        }
    }

    return visit_expr(
        pattern,
        [&](const auto& node) -> ExprPtr {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, ComplexLit> ||
                std::is_same_v<Node, Symbol> ||
                std::is_same_v<Node, Constant>) {
                return arena.make<Node>(node);
            } else if constexpr (std::is_same_v<Node, Unary>) {
                return arena.make<Unary>(node.op, instantiate_pattern(node.operand, matches, arena));
            } else if constexpr (std::is_same_v<Node, Binary>) {
                return arena.make<Binary>(
                    node.op,
                    instantiate_pattern(node.left, matches, arena),
                    instantiate_pattern(node.right, matches, arena));
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> args;
                args.reserve(node.args.size());
                for (ExprPtr arg : node.args) {
                    args.push_back(instantiate_pattern(arg, matches, arena));
                }
                return arena.make<FuncCall>(node.name, std::move(args));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> terms;
                terms.reserve(node.terms.size());
                for (ExprPtr term : node.terms) {
                    terms.push_back(instantiate_pattern(term, matches, arena));
                }
                return arena.make<Sum>(std::move(terms));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> factors;
                factors.reserve(node.factors.size());
                for (ExprPtr factor : node.factors) {
                    factors.push_back(instantiate_pattern(factor, matches, arena));
                }
                return arena.make<Product>(std::move(factors));
            } else if constexpr (std::is_same_v<Node, Integral>) {
                return arena.make<Integral>(
                    instantiate_pattern(node.integrand, matches, arena),
                    node.variable,
                    node.lower.has_value()
                        ? std::optional<ExprPtr>(instantiate_pattern(*node.lower, matches, arena))
                        : std::nullopt,
                    node.upper.has_value()
                        ? std::optional<ExprPtr>(instantiate_pattern(*node.upper, matches, arena))
                        : std::nullopt);
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                return arena.make<Derivative>(
                    instantiate_pattern(node.expression, matches, arena),
                    node.variable,
                    node.order);
            } else if constexpr (std::is_same_v<Node, Limit>) {
                return arena.make<Limit>(
                    instantiate_pattern(node.expression, matches, arena),
                    node.variable,
                    instantiate_pattern(node.point, matches, arena),
                    node.direction);
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                return arena.make<RootOf>(
                    instantiate_pattern(node.polynomial, matches, arena),
                    node.variable,
                    node.root_index);
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::vector<ExprPtr> elements;
                elements.reserve(node.elements.size());
                for (ExprPtr element : node.elements) {
                    elements.push_back(instantiate_pattern(element, matches, arena));
                }
                return arena.make<Matrix>(node.rows, node.cols, std::move(elements));
            } else if constexpr (std::is_same_v<Node, SeriesExp>) {
                std::vector<std::pair<long long, ExprPtr>> terms;
                terms.reserve(node.terms.size());
                for (const auto& [exp, coeff] : node.terms) {
                    terms.push_back({exp, instantiate_pattern(coeff, matches, arena)});
                }
                return arena.make<SeriesExp>(node.var, instantiate_pattern(node.point, matches, arena), std::move(terms), node.order);
            } else {
                return ExprPtr{};
            }
        });
}

CASContext::CASContext() : rewrite_provider_(&default_rewrite_provider()) {
    ::cas::algebra::register_algebraic_simplify_hook(*this);
}

void CASContext::define(const Symbol& symbol, ExprPtr value) {
    variables_[symbol.name] = value;
}

void CASContext::clear_variables() noexcept {
    variables_.clear();
}

std::optional<ExprPtr> CASContext::lookup(const Symbol& symbol) const {
    const auto found = variables_.find(symbol.name);
    if (found == variables_.end()) {
        return std::nullopt;
    }
    return found->second;
}

Assumptions& CASContext::assumptions() noexcept {
    return assumptions_;
}

const Assumptions& CASContext::assumptions() const noexcept {
    return assumptions_;
}

AstArena& CASContext::arena() noexcept {
    return arena_;
}

const AstArena& CASContext::arena() const noexcept {
    return arena_;
}

void CASContext::set_rewrite_provider(const RewriteProvider* provider) noexcept {
    rewrite_provider_ = provider;
}

const RewriteProvider* CASContext::rewrite_provider() const noexcept {
    return rewrite_provider_;
}

void CASContext::enable_trace(bool enabled) noexcept {
    trace_enabled_ = enabled;
    if (!enabled) {
        trace_.clear();
    }
}

const ComputationTrace& CASContext::get_trace() const noexcept {
    return trace_enabled_ ? trace_ : empty_trace();
}

void CASContext::set_timeout(std::chrono::milliseconds timeout) noexcept {
    timeout_ = timeout;
}

void CASContext::set_timeout_check_interval(std::uint64_t interval) noexcept {
    timeout_check_interval_ = (interval < 64U) ? 64U : interval;
}

void CASContext::set_max_simplification_depth(int depth) noexcept {
    max_simplification_depth_ = (depth < 10) ? 10 : depth;
}

void CASContext::set_max_integration_depth(std::size_t depth) noexcept {
    max_integration_depth_ = (depth < 1U) ? 1U : (depth > 128U) ? 128U : depth;
}

void CASContext::set_gcd_error_probability(double prob) noexcept {
    if (prob < 1e-6) prob = 1e-6;
    if (prob > 0.1)  prob = 0.1;
    gcd_error_probability_ = prob;
}

void CASContext::set_numeric_precision_digits(unsigned int digits) noexcept {
    // Clamp: minimum 6 digits (≈ 20 bits MPFR), maximum 10000 (~ 33k bits).
    if (digits < 6U) digits = 6U;
    if (digits > 10000U) digits = 10000U;
    numeric_precision_digits_ = digits;
}

void CASContext::set_max_rootof_explicit_degree(std::size_t deg) noexcept {
    max_rootof_explicit_degree_ = (deg < 1U) ? 1U : deg;
}

void CASContext::set_max_gcd_recursion_depth(std::size_t depth) noexcept {
    max_gcd_recursion_depth_ = (depth < 4U) ? 4U : depth;
}

void CASContext::set_min_gcd_division_steps(std::size_t steps) noexcept {
    min_gcd_division_steps_ = (steps < 1U) ? 1U : steps;
}

void CASContext::set_max_gcd_total_calls(std::size_t n) noexcept {
    max_gcd_total_calls_ = (n < 16U) ? 16U : n;
}

void CASContext::set_max_cyclotomic_n(int n) noexcept {
    max_cyclotomic_n_ = n;
}

void CASContext::set_max_q_alpha_bridge_depth(std::size_t depth) noexcept {
    max_q_alpha_bridge_depth_ = (depth < 8U) ? 8U : depth;
}

void CASContext::set_max_gamma_recursion(std::size_t iters) noexcept {
    max_gamma_recursion_ = (iters < 16U) ? 16U : iters;
}

void CASContext::set_improper_leading_order_scan(std::size_t window) noexcept {
    improper_leading_order_scan_ = (window < 1U) ? 1U : window;
}

void CASContext::set_expand_bessel_recurrence(bool enabled) noexcept {
    expand_bessel_recurrence_ = enabled;
}

void CASContext::set_max_trager_tower_shift_attempts(std::size_t attempts) noexcept {
    max_trager_tower_shift_attempts_ = attempts;
}

Symbol CASContext::make_fresh_symbol(const std::string& prefix) {
    // Probe candidate names of the form "<prefix>_<n>" with monotonically
    // increasing n until the name is not present in the user-defined
    // variable map.  The internal counter guarantees we never return
    // the same name twice from one context.
    std::string base = prefix.empty() ? std::string("_g") : prefix;
    while (true) {
        ++fresh_symbol_counter_;
        std::string candidate = base + "_" + std::to_string(fresh_symbol_counter_);
        if (variables_.find(candidate) == variables_.end()) {
            return Symbol(candidate);
        }
    }
}

void CASContext::clear_caches() noexcept {
    simplify_cache_.clear();
    diff_cache_.clear();
    integrate_cache_.clear();
}

void CASContext::set_caching_enabled(bool enabled) noexcept {
    caching_enabled_ = enabled;
    if (!enabled) {
        clear_caches();
    }
}

void CASContext::set_cache_limit(std::size_t limit) noexcept {
    simplify_cache_.set_max_size(limit);
    diff_cache_.set_max_size(limit);
    integrate_cache_.set_max_size(limit);
}

std::size_t CASContext::get_cache_limit() const noexcept {
    return simplify_cache_.max_size();
}

CacheMetrics CASContext::get_simplify_metrics() const noexcept {
    return simplify_cache_.metrics();
}

CacheMetrics CASContext::get_diff_metrics() const noexcept {
    return diff_cache_.metrics();
}

CacheMetrics CASContext::get_integrate_metrics() const noexcept {
    return integrate_cache_.metrics();
}

Result<ExprPtr> CASContext::simplify(ExprPtr expr) {
    if (!expr) return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot simplify null expression"));

    // F7.0-A4.1: invalidate caches if assumptions have changed since last fill.
    // simplify(abs(x)) result depends on whether x is known positive/negative;
    // a stale cache entry from a different assumption regime would corrupt
    // the result. Cheap check: one std::uint64_t comparison per simplify call.
    {
        const std::uint64_t cur_rev = assumptions_.revision();
        if (cur_rev != last_assumptions_revision_) {
            clear_caches();
            last_assumptions_revision_ = cur_rev;
        }
    }

    if (caching_enabled_ && !trace_enabled_) {
        if (auto cached = simplify_cache_.get(expr)) {
            return ok(*cached);
        }
    }

    const bool owns_operation = !operation_active_;
    if (owns_operation) {
        operation_active_ = true;
        trace_capture_active_ = trace_enabled_;
        trace_.clear();
        ops_count_ = 0;
        operation_started_at_ = std::chrono::steady_clock::now();
    }
    auto result = symbolic::simplify(expr, *this);
    if (owns_operation) {
        operation_active_ = false;
        trace_capture_active_ = false;
        ops_count_ = 0;

        // Post-simplify hook: algebraic extension reduction (e.g. Q(alpha) via bridge).
        // Only runs at the top-level call to avoid O(N) overhead on sub-expressions.
        if (result.is_ok() && post_simplify_hook_) {
            auto hook_res = post_simplify_hook_(result.value(), *this);
            if (hook_res.is_ok() && hook_res.value() != result.value()) {
                // Hook produced a distinct expression; re-simplify to canonical form.
                auto resimplified = symbolic::simplify(hook_res.value(), *this);
                result = resimplified.is_ok() ? resimplified : hook_res;
            }
        }
    }

    if (caching_enabled_ && !trace_enabled_ && result.is_ok()) {
        simplify_cache_.put(expr, result.value());
    }

#ifndef NDEBUG
    // F7.0-A4.2: debug-only canonical invariant check.
    // Reports invariant violations on stderr without aborting — emits one
    // warning per top-level simplify call. Helps catch invariant regressions
    // close to their source during dev/test runs.
    if (owns_operation && result.is_ok()
        && !is_strictly_canonical(result.value())) {
        std::fprintf(stderr,
            "[F7.0-A4.2 WARN] simplify result is not strictly canonical\n");
    }
#endif

    return result;
}

Result<ExprPtr> materialize_expr(ExprPtr expr, AstArena& arena) {
    return materialize_expr_impl(expr, arena);
}

// NOTE: mathematically_equal is defined in cas_algebra (algebraic_equal.cpp)
// so it can use split_num_den for rational cross-multiplication.
// The declaration in cas/symbolic.hpp is satisfied by the algebra library at link time.

void CASContext::collect_garbage(const std::vector<ExprPtr*>& external_roots) {
    AstArena new_arena;
    std::unordered_map<ExprPtr, ExprPtr> cache;

    // 1. Update variables
    for (auto& [name, expr] : variables_) {
        expr = clone_into_arena(expr, new_arena, cache);
    }

    // 2. Update assumptions
    assumptions_.update_roots(new_arena, cache);

    // 3. Update trace
    if (trace_enabled_) {
        for (auto& step : trace_) {
            step.target_before = clone_into_arena(step.target_before, new_arena, cache);
            step.target_after = clone_into_arena(step.target_after, new_arena, cache);
            step.root_after = clone_into_arena(step.root_after, new_arena, cache);
        }
    }

    // 4. Update external roots
    for (auto* root_ptr : external_roots) {
        if (root_ptr && *root_ptr) {
            *root_ptr = clone_into_arena(*root_ptr, new_arena, cache);
        }
    }

    // 5. Update caches
    if (!simplify_cache_.empty()) {
        CacheContainer<ExprPtr, ExprPtr, ExprHash, ExprEqual> new_simplify_cache(simplify_cache_.max_size());
        for (auto& it : simplify_cache_) {
            new_simplify_cache.put(clone_into_arena(it.first, new_arena, cache), clone_into_arena(it.second.first, new_arena, cache));
        }
        new_simplify_cache.metrics() = simplify_cache_.metrics();
        simplify_cache_ = std::move(new_simplify_cache);
    }

    if (!diff_cache_.empty()) {
        CacheContainer<DiffKey, ExprPtr, DiffHash> new_diff_cache(diff_cache_.max_size());
        for (auto& it : diff_cache_) {
            const auto& key = it.first;
            DiffKey new_key = {
                .expr = clone_into_arena(key.expr, new_arena, cache),
                .var_name = key.var_name,
                .order = key.order
            };
            new_diff_cache.put(new_key, clone_into_arena(it.second.first, new_arena, cache));
        }
        new_diff_cache.metrics() = diff_cache_.metrics();
        diff_cache_ = std::move(new_diff_cache);
    }

    if (!integrate_cache_.empty()) {
        CacheContainer<IntegrateKey, ExprPtr, IntegrateHash> new_integrate_cache(integrate_cache_.max_size());
        for (auto& it : integrate_cache_) {
            const auto& key = it.first;
            IntegrateKey new_key = {
                .expr = clone_into_arena(key.expr, new_arena, cache),
                .var_name = key.var_name
            };
            new_integrate_cache.put(new_key, clone_into_arena(it.second.first, new_arena, cache));
        }
        new_integrate_cache.metrics() = integrate_cache_.metrics();
        integrate_cache_ = std::move(new_integrate_cache);
    }

    // 6. Swap arena
    arena_ = std::move(new_arena);
}

} // namespace cas::symbolic
