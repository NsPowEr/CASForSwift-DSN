#pragma once

#include "cas/ast.hpp"
#include "cas/cas_cache_keys.hpp"
#include "cas/cas_context_params.hpp"
#include "cas/result.hpp"
#include "cas/trace.hpp"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <iterator>
#include <random>

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra {
[[nodiscard]] Result<ExprPtr> polynomial_gcd(ExprPtr p, ExprPtr q, const Symbol& var, symbolic::CASContext& ctx);
}

namespace cas::symbolic {

using MatchMap = std::unordered_map<std::string, ExprPtr>;


class Assumptions;

enum class TraversalStrategy {
    BottomUp,
    TopDown,
    FixPoint,
};

enum class TermOrderRelation {
    Less,
    Equivalent,
    Greater,
    Incomparable,
};

struct RewriteRule {
    ExprPtr pattern;
    ExprPtr replacement;
    std::function<bool(const MatchMap&)> condition;
};

class RewriteProvider {
public:
    virtual ~RewriteProvider() = default;

    [[nodiscard]] virtual Result<ExprPtr> try_rewrite(
        ExprPtr expr,
        AstArena& arena,
        const Assumptions* assumptions,
        CASContext* context = nullptr) const = 0;
};

struct RangeAssumption {
    ExprPtr lower;
    ExprPtr upper;
};

enum class RelType : std::uint8_t {
    Less,          // <
    LessEqual      // <=
};

struct Relation {
    ExprPtr target;
    RelType type;
};

enum class Domain : std::uint8_t {
    Complex,
    Real,
    Integer,
    Rational,
    Natural,     // >= 0
    Positive,    // > 0
    Negative,    // < 0
    NonZero,
};

class Assumptions {
public:
    void assume_real(const Symbol& symbol);
    void assume_positive(const Symbol& symbol);
    void assume_integer(const Symbol& symbol);
    void assume_nonzero(const Symbol& symbol);
    void assume_in_range(const Symbol& symbol, ExprPtr lower, ExprPtr upper);
    void assume_domain(const Symbol& symbol, Domain domain);

    // Advanced Assumptions (F9-A4)
    void assume_greater(ExprPtr lhs, ExprPtr rhs);
    void assume_greater_equal(ExprPtr lhs, ExprPtr rhs);
    void assume(ExprPtr condition);

    [[nodiscard]] bool is_real(const Symbol& symbol) const;
    [[nodiscard]] bool is_real(ExprPtr expr) const;
    [[nodiscard]] bool is_positive(const Symbol& symbol) const;
    [[nodiscard]] bool is_positive(ExprPtr expr) const;
    [[nodiscard]] bool is_nonnegative(ExprPtr expr) const;
    [[nodiscard]] bool is_negative(ExprPtr expr) const;
    [[nodiscard]] bool is_greater(ExprPtr lhs, ExprPtr rhs) const;
    [[nodiscard]] bool is_greater_equal(ExprPtr lhs, ExprPtr rhs) const;
    [[nodiscard]] bool is_nonzero(const Symbol& symbol) const;
    [[nodiscard]] bool is_nonzero(ExprPtr expr) const;
    [[nodiscard]] bool could_be_zero(const Symbol& symbol) const;
    [[nodiscard]] bool could_be_zero(ExprPtr expr) const;
    [[nodiscard]] bool is_integer(const Symbol& symbol) const;
    [[nodiscard]] bool is_integer(ExprPtr expr) const;
    [[nodiscard]] Domain get_domain(const Symbol& symbol) const;
    [[nodiscard]] std::optional<RangeAssumption> get_range(const Symbol& symbol) const;
    [[nodiscard]] Result<void> check_consistency() const;

    void update_roots(AstArena& target, std::unordered_map<ExprPtr, ExprPtr>& cache);

    // F7.0-A4.1: monotonically increasing revision counter, incremented by
    // every mutator (assume_*, assume()). Used by CASContext::simplify() to
    // detect assumption changes since the last cache fill, and invalidate
    // simplify_cache_ / integrate_cache_ on mismatch.
    //
    // Mathematical correctness: simplify(abs(x)) returns x when x is known
    // positive, but -x when x is known negative — without cache invalidation
    // on assumption change, a stale entry would corrupt the session.
    [[nodiscard]] std::uint64_t revision() const noexcept { return revision_; }

private:
    // F7.0-A4.1: bump in every mutator above (assume_*, assume).
    std::uint64_t revision_{0};

    [[nodiscard]] bool prove_relation(ExprPtr start, ExprPtr end, bool strict, std::unordered_set<const ExprNode*>& visited) const;
    [[nodiscard]] bool prove_positive_linear(ExprPtr expr) const;
    [[nodiscard]] bool prove_positive_product(const Product& prod) const;

    std::unordered_set<std::string> real_symbols_;
    std::unordered_set<std::string> positive_symbols_;
    std::unordered_set<std::string> negative_symbols_;
    std::unordered_set<std::string> integer_symbols_;
    std::unordered_set<std::string> nonzero_symbols_;
    std::unordered_map<std::string, RangeAssumption> range_symbols_;
    std::unordered_map<std::string, Domain> symbol_domains_;

    // Relation graph for deduction chains
    std::unordered_map<ExprPtr, std::vector<Relation>, ExprHash, ExprEqual> relations_;
};

class Substituter;

struct SimplifyHints {
    bool fold_exp_products = false;
};

class CASContext;

class [[nodiscard]] ScopedSimplifyHint {
    CASContext& ctx_;
    SimplifyHints old_hints_;
public:
    ScopedSimplifyHint(CASContext& ctx, SimplifyHints hints) noexcept;
    ~ScopedSimplifyHint();
    ScopedSimplifyHint(const ScopedSimplifyHint&) = delete;
    ScopedSimplifyHint& operator=(const ScopedSimplifyHint&) = delete;
};

// CASContext inherits CASContextParams to expose all algorithm-tuning
// getter/setter pairs directly on ctx without any call-site changes.
// Simple inline setters/getters live in cas_context_params.hpp.
// Setters with clamping/validation are declared here and defined in
// context_core.cpp (they access protected fields via inheritance).
class CASContext : public CASContextParams {
public:
    CASContext();

    [[nodiscard]] SimplifyHints hints() const noexcept { return hints_; }
    void set_hints(SimplifyHints hints) noexcept { hints_ = hints; }
    [[nodiscard]] ScopedSimplifyHint with_hint(SimplifyHints hints) noexcept;

    void define(const Symbol& symbol, ExprPtr value);
    void clear_variables() noexcept;
    [[nodiscard]] std::optional<ExprPtr> lookup(const Symbol& symbol) const;
    [[nodiscard]] const std::unordered_map<std::string, ExprPtr>& variables() const noexcept { return variables_; }

    [[nodiscard]] Assumptions& assumptions() noexcept;
    [[nodiscard]] const Assumptions& assumptions() const noexcept;

    [[nodiscard]] AstArena& arena() noexcept;
    [[nodiscard]] const AstArena& arena() const noexcept;

    void set_rewrite_provider(const RewriteProvider* provider) noexcept;
    [[nodiscard]] const RewriteProvider* rewrite_provider() const noexcept;

    void enable_trace(bool enabled) noexcept;
    [[nodiscard]] const ComputationTrace& get_trace() const noexcept;
    void set_timeout(std::chrono::milliseconds timeout) noexcept;
    // Read-only access to the per-operation timeout budget.  Used by long-running
    // outer loops (e.g. factor_polynomial_tower shift search) that own their wall-clock
    // deadline independently of the CASContext operation-start mechanism.
    [[nodiscard]] std::chrono::milliseconds timeout() const noexcept { return timeout_; }

    // Setters with clamping/validation (implemented in context_core.cpp).
    // Getters for these are inherited from CASContextParams.
    void set_timeout_check_interval(std::uint64_t interval) noexcept;
    void set_max_simplification_depth(int depth) noexcept;
    void set_max_integration_depth(std::size_t depth) noexcept;
    void set_gcd_error_probability(double prob) noexcept;
    void set_zippel_error_probability(double prob) noexcept;
    void set_zippel_density_threshold(double t) noexcept;
    void set_numeric_precision_digits(unsigned int digits) noexcept;
    void set_max_rootof_explicit_degree(std::size_t deg) noexcept;
    void set_max_gcd_recursion_depth(std::size_t depth) noexcept;
    void set_min_gcd_division_steps(std::size_t steps) noexcept;
    void set_max_gcd_total_calls(std::size_t n) noexcept;
    void set_max_cyclotomic_n(int n) noexcept;
    void set_max_q_alpha_bridge_depth(std::size_t depth) noexcept;
    void set_max_gamma_recursion(std::size_t iters) noexcept;
    void set_improper_leading_order_scan(std::size_t window) noexcept;
    void set_expand_bessel_recurrence(bool enabled) noexcept;
    void set_max_trager_tower_shift_attempts(std::size_t attempts) noexcept;

    // HC-004: fresh symbol generator.  Returns a Symbol whose name is unique
    // within this context across all previous make_fresh_symbol calls AND
    // does not collide with any name currently registered through `define`.
    [[nodiscard]] Symbol make_fresh_symbol(const std::string& prefix);

    void interrupt() noexcept { interrupted_ = true; }
    void clear_interrupt() noexcept { interrupted_ = false; }
    [[nodiscard]] bool is_interrupted() const noexcept { return interrupted_; }

    // F7.0-A3.3: poll-point helper for heavy non-simplify loops
    // (polynomial GCD, Groebner, factorization, matrix Bareiss, etc.).
    // Returns a Timeout error with explicit "cancelled" message if the
    // interrupt flag has been set asynchronously by an external controller.
    // Cheap (atomic load + branch); call inside inner loops at safe points.
    [[nodiscard]] Result<void> check_interrupt() const noexcept {
        if (interrupted_) {
            return Result<void>(CASError{
                .kind = CASErrorKind::Timeout,
                .message = "Operation cancelled by interrupt request",
                .hint = std::nullopt,
            });
        }
        return ok();
    }

    [[nodiscard]] std::mt19937& rng() noexcept { return rng_; }

    [[nodiscard]] Result<ExprPtr> simplify(ExprPtr expr);
    [[nodiscard]] Result<ExprPtr> substitute(ExprPtr expr, const Symbol& variable, ExprPtr value);

    // Optional post-simplify hook invoked at the top-level simplify call.
    // The hook receives the already-simplified expression and may return a
    // further-reduced form (e.g. algebraic extension reduction via Q(alpha)).
    // If the hook returns a different pointer, a second standard simplify pass
    // is run on the result.  Set to nullptr to disable.
    using PostSimplifyHook = std::function<Result<ExprPtr>(ExprPtr, CASContext&)>;
    void set_post_simplify_hook(PostSimplifyHook hook) {
        post_simplify_hook_ = std::move(hook);
        clear_caches();
    }
    void clear_post_simplify_hook() noexcept {
        post_simplify_hook_ = nullptr;
        clear_caches();
    }
    [[nodiscard]] bool has_post_simplify_hook() const noexcept {
        return static_cast<bool>(post_simplify_hook_);
    }

    void collect_garbage(const std::vector<ExprPtr*>& external_roots = {});

    // Caching methods
    void clear_caches() noexcept;
    void set_caching_enabled(bool enabled) noexcept;
    [[nodiscard]] bool is_caching_enabled() const noexcept { return caching_enabled_; }

    void set_cache_limit(std::size_t limit) noexcept;
    [[nodiscard]] std::size_t get_cache_limit() const noexcept;

    [[nodiscard]] CacheMetrics get_simplify_metrics() const noexcept;
    [[nodiscard]] CacheMetrics get_diff_metrics() const noexcept;
    [[nodiscard]] CacheMetrics get_integrate_metrics() const noexcept;

    private:
    friend class Substituter;
    friend Result<ExprPtr> simplify(ExprPtr expr, CASContext& context);
    friend Result<ExprPtr> substitute(ExprPtr expr, const Symbol& variable, ExprPtr value, CASContext& context);
    friend Result<bool> mathematically_equal(ExprPtr lhs, ExprPtr rhs, CASContext& context);
    friend Result<ExprPtr> cas::algebra::polynomial_gcd(ExprPtr p, ExprPtr q, const Symbol& var, CASContext& ctx);

    std::unordered_map<std::string, ExprPtr> variables_;
    Assumptions assumptions_;
    AstArena arena_;
    const RewriteProvider* rewrite_provider_{nullptr};
    bool trace_enabled_{false};
    bool trace_capture_active_{false};
    bool operation_active_{false};
    bool caching_enabled_{true};
    ComputationTrace trace_;
    std::chrono::milliseconds timeout_{1000};
    std::chrono::steady_clock::time_point operation_started_at_{};
    std::uint64_t ops_count_{0};
    std::uint64_t fresh_symbol_counter_{0U};
    PostSimplifyHook post_simplify_hook_{nullptr};
    std::atomic_bool interrupted_{false};
    std::mt19937 rng_{0xCAFEBABEULL};

public:
// Performance Caches
struct DiffKey {
    ExprPtr expr;
    std::string var_name;
    unsigned int order;

    bool operator==(const DiffKey& other) const {
        return order == other.order && var_name == other.var_name && ExprEqual{}(expr, other.expr);
    }
};
struct DiffHash {
    std::size_t operator()(const DiffKey& key) const {
        std::size_t h = ExprHash{}(key.expr);
        h ^= std::hash<std::string>{}(key.var_name) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<unsigned int>{}(key.order) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

struct IntegrateKey {
    ExprPtr expr;
    std::string var_name;

    bool operator==(const IntegrateKey& other) const {
        return var_name == other.var_name && ExprEqual{}(expr, other.expr);
    }
};
struct IntegrateHash {
    std::size_t operator()(const IntegrateKey& key) const {
        std::size_t h = ExprHash{}(key.expr);
        h ^= std::hash<std::string>{}(key.var_name) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

CacheContainer<ExprPtr, ExprPtr, ExprHash, ExprEqual> simplify_cache_;
CacheContainer<DiffKey, ExprPtr, DiffHash> diff_cache_;
CacheContainer<IntegrateKey, ExprPtr, IntegrateHash> integrate_cache_;

// F7.0-A4.1: last Assumptions::revision() observed by the cache layer.
// On simplify(), if assumptions_.revision() != this, clear all caches and
// update. Prevents stale-cache corruption when user changes assumptions
// mid-session (e.g. assume(x>0); simplify(abs(x)); assume(x<0); ...).
mutable std::uint64_t last_assumptions_revision_{0};
SimplifyHints hints_;
};

inline ScopedSimplifyHint::ScopedSimplifyHint(CASContext& ctx, SimplifyHints hints) noexcept
    : ctx_(ctx), old_hints_(ctx.hints()) {
    ctx.set_hints(hints);
}

inline ScopedSimplifyHint::~ScopedSimplifyHint() {
    ctx_.set_hints(old_hints_);
}

inline ScopedSimplifyHint CASContext::with_hint(SimplifyHints hints) noexcept {
    return ScopedSimplifyHint(*this, hints);
}

[[nodiscard]] int canonical_compare(ExprPtr lhs, ExprPtr rhs) noexcept;

// F7.0-A4.2: post-simplify canonical-form invariant check.
// Returns true if `expr` and every reachable sub-expression respect the
// invariants that simplify() is supposed to maintain:
//   - Sum::terms sorted by polynomial degree descending then canonical_compare,
//     no nested Sum, no exact-zero IntegerLit / RationalLit summand.
//   - Product::factors sorted by canonical_compare, no nested Product, no
//     exact-one IntegerLit / RationalLit factor.
//   - All Sum/Product nodes have ≥ 2 operands (singletons collapsed).
//
// Used in DEBUG builds via an assert at the end of CASContext::simplify()
// to catch invariant violations close to their source. In release builds
// the function is still available for explicit verification but the assert
// is compiled out.
[[nodiscard]] bool is_strictly_canonical(ExprPtr expr) noexcept;
[[nodiscard]] TermOrderRelation compare_rewrite_terms(ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] bool rewrite_rule_is_oriented(const RewriteRule& rule);
[[nodiscard]] bool is_strongly_normalizing(const std::vector<RewriteRule>& rules);
[[nodiscard]] bool match_pattern(ExprPtr expr, ExprPtr pattern, MatchMap& out_matches);
[[nodiscard]] bool match_ac_pattern(ExprPtr expr, ExprPtr pattern, MatchMap& out_matches);
[[nodiscard]] Result<ExprPtr> apply_rule(ExprPtr expr, const RewriteRule& rule, TraversalStrategy strategy, AstArena& arena);
[[nodiscard]] Result<ExprPtr> apply_rule_set(ExprPtr expr, const std::vector<RewriteRule>& rules, AstArena& arena);
[[nodiscard]] Result<ExprPtr> materialize_expr(ExprPtr expr, AstArena& arena);
[[nodiscard]] const RewriteProvider& default_rewrite_provider();
[[nodiscard]] Result<ExprPtr> simplify(ExprPtr expr, AstArena& arena);
[[nodiscard]] Result<ExprPtr> simplify(ExprPtr expr, CASContext& context);
[[nodiscard]] Result<ExprPtr> substitute(ExprPtr expr, const Symbol& variable, ExprPtr value, CASContext& context);
[[nodiscard]] Result<bool> mathematically_equal(ExprPtr lhs, ExprPtr rhs, CASContext& context);

// L2-19: decidable subset of transcendental equivalence via Risch-style
// log/exp/trig normalisation.  Applies opt-in expansion rules to lhs/rhs
// (log(x*y) -> log(x)+log(y) under x>0, y>0; exp(x+y) -> exp(x)*exp(y);
// exp(n*ln(x)) -> x^n under x>0; sin^2+cos^2 collapse) before delegating
// to mathematically_equal.  Returns false (not Unimplemented) for cases
// outside the decidable subset — Richardson's theorem precludes a total
// decision procedure.
[[nodiscard]] Result<bool> mathematically_equal_subset_risch(
    ExprPtr lhs, ExprPtr rhs, CASContext& context);

}  // namespace cas::symbolic
