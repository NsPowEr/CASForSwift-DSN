#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/trace.hpp"

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

namespace cas::symbolic {
class CASContext;
}

namespace cas::algebra {
[[nodiscard]] Result<ExprPtr> polynomial_gcd(ExprPtr p, ExprPtr q, const Symbol& var, symbolic::CASContext& ctx);
}

namespace cas::symbolic {

using MatchMap = std::unordered_map<std::string, ExprPtr>;

struct CacheMetrics {
    std::uint64_t hits{0};
    std::uint64_t misses{0};
    std::uint64_t evictions{0};
};

template <typename Key, typename Value, typename Hash = std::hash<Key>, typename Equal = std::equal_to<Key>>
class CacheContainer {
public:
    using ListType = std::list<Key>;
    using MapType = std::unordered_map<Key, std::pair<Value, typename ListType::iterator>, Hash, Equal>;

    explicit CacheContainer(std::size_t max_size = 1000) : max_size_(max_size) {}

    void set_max_size(std::size_t size) {
        max_size_ = size;
        evict_if_needed();
    }

    [[nodiscard]] std::size_t max_size() const noexcept { return max_size_; }
    [[nodiscard]] std::size_t size() const noexcept { return map_.size(); }

    [[nodiscard]] std::optional<Value> get(const Key& key) {
        auto it = map_.find(key);
        if (it == map_.end()) {
            metrics_.misses++;
            return std::nullopt;
        }
        metrics_.hits++;
        list_.splice(list_.begin(), list_, it->second.second);
        return it->second.first;
    }

    void put(const Key& key, Value value) {
        if (max_size_ == 0) return;

        auto it = map_.find(key);
        if (it != map_.end()) {
            it->second.first = value;
            list_.splice(list_.begin(), list_, it->second.second);
            return;
        }

        list_.push_front(key);
        map_[key] = {value, list_.begin()};
        evict_if_needed();
    }

    void clear() noexcept {
        map_.clear();
        list_.clear();
    }

    [[nodiscard]] CacheMetrics& metrics() noexcept { return metrics_; }
    [[nodiscard]] const CacheMetrics& metrics() const noexcept { return metrics_; }
    void reset_metrics() noexcept { metrics_ = {}; }

    [[nodiscard]] auto begin() { return map_.begin(); }
    [[nodiscard]] auto end() { return map_.end(); }
    [[nodiscard]] auto begin() const { return map_.begin(); }
    [[nodiscard]] auto end() const { return map_.end(); }
    [[nodiscard]] bool empty() const noexcept { return map_.empty(); }

private:
    void evict_if_needed() {
        while (map_.size() > max_size_ && !list_.empty()) {
            Key last = list_.back();
            list_.pop_back();
            map_.erase(last);
            metrics_.evictions++;
        }
    }

    std::size_t max_size_;
    MapType map_;
    ListType list_;
    CacheMetrics metrics_;
};

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

private:
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

class CASContext {
public:
    CASContext();

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
    void set_timeout_check_interval(std::uint64_t interval) noexcept;
    [[nodiscard]] std::uint64_t timeout_check_interval() const noexcept { return timeout_check_interval_; }

    void set_max_simplification_depth(int depth) noexcept;
    [[nodiscard]] int max_simplification_depth() const noexcept { return max_simplification_depth_; }

    void set_max_integration_depth(std::size_t depth) noexcept;
    [[nodiscard]] std::size_t max_integration_depth() const noexcept { return max_integration_depth_; }

    void set_gcd_error_probability(double prob) noexcept;
    [[nodiscard]] double gcd_error_probability() const noexcept { return gcd_error_probability_; }

    void set_max_trig_power_reduction(long long n) noexcept { max_trig_power_reduction_ = n; }
    [[nodiscard]] long long max_trig_power_reduction() const noexcept { return max_trig_power_reduction_; }

    void set_max_rootof_explicit_degree(std::size_t deg) noexcept;
    [[nodiscard]] std::size_t max_rootof_explicit_degree() const noexcept { return max_rootof_explicit_degree_; }

    void set_max_gcd_recursion_depth(std::size_t depth) noexcept;
    [[nodiscard]] std::size_t max_gcd_recursion_depth() const noexcept { return max_gcd_recursion_depth_; }

    void set_min_gcd_division_steps(std::size_t steps) noexcept;
    [[nodiscard]] std::size_t min_gcd_division_steps() const noexcept { return min_gcd_division_steps_; }

    // -1 = auto-derive from polynomial degree (2*(deg+1)); set to higher value for unusual composites
    void set_max_cyclotomic_n(int n) noexcept;
    [[nodiscard]] int max_cyclotomic_n() const noexcept { return max_cyclotomic_n_; }

    // HC-001..006 configurable knobs.
    void set_max_q_alpha_bridge_depth(std::size_t depth) noexcept;
    [[nodiscard]] std::size_t max_q_alpha_bridge_depth() const noexcept { return max_q_alpha_bridge_depth_; }
    void set_max_gamma_recursion(std::size_t iters) noexcept;
    [[nodiscard]] std::size_t max_gamma_recursion() const noexcept { return max_gamma_recursion_; }
    void set_improper_leading_order_scan(std::size_t window) noexcept;
    [[nodiscard]] std::size_t improper_leading_order_scan() const noexcept { return improper_leading_order_scan_; }

    // Opt-in: expand BesselJ(n, x) and BesselY(n, x) with integer n >= 2 via the
    // three-term recurrence  J_{n}(x) = (2(n-1)/x) J_{n-1}(x) - J_{n-2}(x)
    // until the order reaches {0, 1}.  Default false because the expanded form
    // is rarely what the user wants in symbolic answers.
    void set_expand_bessel_recurrence(bool enabled) noexcept;
    [[nodiscard]] bool expand_bessel_recurrence() const noexcept { return expand_bessel_recurrence_; }

    // L3-06: maximum shift attempts (s1, s2) when searching for a square-free
    // composite Trager norm Res_y2(m2, Res_y1(m1, f(x - s1*y1 - s2*y2))).
    // Default 0 → auto-derive from discriminant collision bound.
    void set_max_trager_tower_shift_attempts(std::size_t attempts) noexcept;
    [[nodiscard]] std::size_t max_trager_tower_shift_attempts() const noexcept { return max_trager_tower_shift_attempts_; }

    // HC-004: fresh symbol generator.  Returns a Symbol whose name is unique
    // within this context across all previous make_fresh_symbol calls AND
    // does not collide with any name currently registered through `define`.
    // The returned Symbol is guaranteed unique by construction: an internal
    // counter is monotonically incremented and the candidate name is
    // probed against the user-defined variable map until a free slot is
    // found.
    [[nodiscard]] Symbol make_fresh_symbol(const std::string& prefix);

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
    std::uint64_t timeout_check_interval_{1024U};
    int max_simplification_depth_{300};
    std::size_t max_integration_depth_{16U};
    double gcd_error_probability_{0.001};
    std::size_t max_rootof_explicit_degree_{2U};
    std::size_t max_gcd_recursion_depth_{16U};
    std::size_t min_gcd_division_steps_{8U};
    int max_cyclotomic_n_{-1};
    std::size_t max_q_alpha_bridge_depth_{256U};
    std::size_t max_gamma_recursion_{1024U};
    std::size_t improper_leading_order_scan_{8U};
    bool expand_bessel_recurrence_{false};
    long long max_trig_power_reduction_{32LL};
    std::size_t max_trager_tower_shift_attempts_{0U};
    std::uint64_t fresh_symbol_counter_{0U};
    PostSimplifyHook post_simplify_hook_{nullptr};

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
};

[[nodiscard]] int canonical_compare(ExprPtr lhs, ExprPtr rhs) noexcept;
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
