#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/trace.hpp"

#include <chrono>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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

class Assumptions {
public:
    void assume_real(const Symbol& symbol);
    void assume_positive(const Symbol& symbol);
    void assume_integer(const Symbol& symbol);
    void assume_nonzero(const Symbol& symbol);
    void assume_in_range(const Symbol& symbol, ExprPtr lower, ExprPtr upper);

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
    [[nodiscard]] std::optional<RangeAssumption> get_range(const Symbol& symbol) const;
    [[nodiscard]] Result<void> check_consistency() const;

    void update_roots(AstArena& target, std::unordered_map<ExprPtr, ExprPtr>& cache);

private:
    [[nodiscard]] bool prove_relation(ExprPtr start, ExprPtr end, bool strict, std::unordered_set<const ExprNode*>& visited) const;
    [[nodiscard]] bool prove_positive_linear(ExprPtr expr) const;

    std::unordered_set<std::string> real_symbols_;
    std::unordered_set<std::string> positive_symbols_;
    std::unordered_set<std::string> integer_symbols_;
    std::unordered_set<std::string> nonzero_symbols_;
    std::unordered_map<std::string, RangeAssumption> range_symbols_;

    // Relation graph for deduction chains
    std::unordered_map<ExprPtr, std::vector<Relation>, ExprHash, ExprEqual> relations_;
};

class Substituter;

class CASContext {
public:
    CASContext();

    void define(const Symbol& symbol, ExprPtr value);
    [[nodiscard]] std::optional<ExprPtr> lookup(const Symbol& symbol) const;

    [[nodiscard]] Assumptions& assumptions() noexcept;
    [[nodiscard]] const Assumptions& assumptions() const noexcept;

    [[nodiscard]] AstArena& arena() noexcept;
    [[nodiscard]] const AstArena& arena() const noexcept;

    void set_rewrite_provider(const RewriteProvider* provider) noexcept;
    [[nodiscard]] const RewriteProvider* rewrite_provider() const noexcept;

    void enable_trace(bool enabled) noexcept;
    [[nodiscard]] const ComputationTrace& get_trace() const noexcept;
    void set_timeout(std::chrono::milliseconds timeout) noexcept;

    [[nodiscard]] Result<ExprPtr> simplify(ExprPtr expr);
    [[nodiscard]] Result<ExprPtr> substitute(ExprPtr expr, const Symbol& variable, ExprPtr value);

    void collect_garbage(const std::vector<ExprPtr*>& external_roots = {});

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
    ComputationTrace trace_;
    std::chrono::milliseconds timeout_{1000};
    std::chrono::steady_clock::time_point operation_started_at_{};
    std::uint64_t ops_count_{0};
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

}  // namespace cas::symbolic
