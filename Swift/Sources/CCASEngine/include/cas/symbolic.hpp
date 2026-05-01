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
        const Assumptions* assumptions) const = 0;
};

struct RangeAssumption {
    ExprPtr lower;
    ExprPtr upper;
};

class Assumptions {
public:
    void assume_real(const Symbol& symbol);
    void assume_positive(const Symbol& symbol);
    void assume_integer(const Symbol& symbol);
    void assume_nonzero(const Symbol& symbol);
    void assume_in_range(const Symbol& symbol, ExprPtr lower, ExprPtr upper);

    [[nodiscard]] bool is_real(const Symbol& symbol) const;
    [[nodiscard]] bool is_positive(const Symbol& symbol) const;
    [[nodiscard]] bool is_nonzero(const Symbol& symbol) const;
    [[nodiscard]] bool could_be_zero(const Symbol& symbol) const;
    [[nodiscard]] bool is_integer(const Symbol& symbol) const;
    [[nodiscard]] std::optional<RangeAssumption> get_range(const Symbol& symbol) const;
    [[nodiscard]] Result<void> check_consistency() const;

private:
    std::unordered_set<std::string> real_symbols_;
    std::unordered_set<std::string> positive_symbols_;
    std::unordered_set<std::string> integer_symbols_;
    std::unordered_set<std::string> nonzero_symbols_;
    std::unordered_map<std::string, RangeAssumption> range_symbols_;
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
