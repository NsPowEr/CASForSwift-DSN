#pragma once

#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include <chrono>
#include <functional>
#include <string>
#include <unordered_set>
#include <vector>
#include <optional>

namespace cas::symbolic {

// Helper functions for simplification logic
namespace detail {

thread_local extern int simplification_depth;
constexpr int MAX_SIMPLIFICATION_DEPTH = 300;

// Fingerprint set for cycle detection: tracks ExprPtr nodes currently
// being simplified on this thread's call stack. If a node is re-entered
// before its first simplification completes, a rewrite cycle is detected.
thread_local extern std::unordered_set<ExprPtr, ExprHash> active_simplify_nodes;

struct DepthGuard {
    DepthGuard();
    ~DepthGuard();
    bool exceeded() const;
};

// RAII guard that inserts expr into active_simplify_nodes on construction
// and removes it on destruction. cycle_detected() returns true if the
// node was already present (i.e., we are inside a cycle).
struct CycleGuard {
    ExprPtr expr_;
    bool cycle_{false};

    explicit CycleGuard(ExprPtr expr) : expr_(expr) {
        cycle_ = !active_simplify_nodes.insert(expr).second;
    }
    ~CycleGuard() {
        if (!cycle_) {
            active_simplify_nodes.erase(expr_);
        }
    }
    [[nodiscard]] bool cycle_detected() const noexcept { return cycle_; }
};

struct LiteralRational {
    Rational value;
    bool exact{false};
};

struct MonomialKey {
    std::vector<std::pair<ExprPtr, BigInt>> factors;
    
    bool operator<(const MonomialKey& other) const {
        if (factors.size() != other.factors.size()) return factors.size() < other.factors.size();
        for (std::size_t i = 0; i < factors.size(); ++i) {
            int cmp = canonical_compare(factors[i].first, other.factors[i].first);
            if (cmp != 0) return cmp < 0;
            if (factors[i].second != other.factors[i].second) {
                return factors[i].second < other.factors[i].second;
            }
        }
        return false;
    }
};

struct MonomialTerm {
    Rational coefficient;
    MonomialKey key;
};

[[nodiscard]] bool is_odd_parity_function(BuiltinOp op);
[[nodiscard]] bool is_even_parity_function(BuiltinOp op);
[[nodiscard]] bool is_parity_rewrite_function(BuiltinOp op);
[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message);
[[nodiscard]] ExprPtr make_integer(AstArena& arena, BigInt value);
[[nodiscard]] ExprPtr make_rational(AstArena& arena, Rational value);
[[nodiscard]] Result<ExprPtr> make_rational_result(AstArena& arena, Rational value);
[[nodiscard]] ExprPtr make_constant(AstArena& arena, MathConstant value);
[[nodiscard]] bool is_zero_expr(ExprPtr expr);
[[nodiscard]] bool is_one_expr(ExprPtr expr);
[[nodiscard]] bool is_constant_expr(ExprPtr expr, MathConstant constant);
[[nodiscard]] Result<bool> try_get_exact_rational(ExprPtr expr, LiteralRational& out);
[[nodiscard]] Rational decimal_to_rational(const DecimalLit& node);
[[nodiscard]] std::optional<BigInt> try_get_integer_exponent(ExprPtr expr);
[[nodiscard]] Rational pow_rational_nonnegative(Rational base, BigInt exponent);
[[nodiscard]] bool is_known_positive_constant(MathConstant value) noexcept;
[[nodiscard]] bool is_known_nonnegative_constant(MathConstant value) noexcept;
[[nodiscard]] bool expr_ptr_sequence_identical(const std::vector<ExprPtr>& lhs, const std::vector<ExprPtr>& rhs) noexcept;
[[nodiscard]] int saturating_add_degree(int lhs, int rhs) noexcept;
[[nodiscard]] int polynomial_degree(ExprPtr expr) noexcept;

class Simplifier {
public:
    explicit Simplifier(
        AstArena& arena,
        const Assumptions* assumptions = nullptr,
        const RewriteProvider* rewrite_provider = nullptr,
        ComputationTrace* trace = nullptr,
        bool trace_enabled = false,
        std::chrono::milliseconds timeout = std::chrono::milliseconds(0),
        std::chrono::steady_clock::time_point* operation_started_at = nullptr,
        std::uint64_t* ops_count = nullptr,
        CASContext* context = nullptr);

    [[nodiscard]] Result<ExprPtr> simplify_expr(ExprPtr expr);

private:
    class ScopedFrame {
    public:
        ScopedFrame(Simplifier& owner, std::function<ExprPtr(ExprPtr)> builder);
        ~ScopedFrame();

    private:
        Simplifier& owner_;
        bool active_{false};
    };

    // Forward declarations of simplification methods
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const IntegerLit& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const RationalLit& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const DecimalLit& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Symbol& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Constant& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Unary& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Binary& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const FuncCall& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Sum& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Product& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Integral& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Derivative& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Limit& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const RootOf& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Matrix& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const SeriesExp& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Quantity& node);

    // Template and generic node handlers
    template <typename Node>
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const Node& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(const DecimalLit& node);
    [[nodiscard]] Result<ExprPtr> simplify_node(const ExprNode& node);
    template <typename Node>
    [[nodiscard]] Result<ExprPtr> simplify_passthrough(ExprPtr original, const Node& node);

    // Specialized simplification methods
    [[nodiscard]] Result<ExprPtr> simplify_additive_chain_fast(ExprPtr original);
    [[nodiscard]] Result<ExprPtr> simplify_power(ExprPtr base, ExprPtr exponent, ExprPtr target_before = ExprPtr{});
    [[nodiscard]] Result<ExprPtr> simplify_sum_terms(const std::vector<ExprPtr>& terms, ExprPtr target_before = ExprPtr{}, bool inputs_are_simplified = false);
    [[nodiscard]] Result<ExprPtr> simplify_product_factors(const std::vector<ExprPtr>& factors, ExprPtr target_before = ExprPtr{}, bool inputs_are_simplified = false);

    // Rewrite and trace helpers
    [[nodiscard]] bool may_rewrite_function_call(BuiltinOp op, const std::vector<ExprPtr>& args) const;
    [[nodiscard]] bool may_rewrite_sum_terms(const std::vector<ExprPtr>& terms) const;
    [[nodiscard]] bool may_rewrite_power(ExprPtr base, ExprPtr exponent) const;
    [[nodiscard]] Result<void> check_timeout();
    [[nodiscard]] ExprPtr build_root_after(ExprPtr target_after);
    void append_trace(RuleId rule_id, ExprPtr before, ExprPtr after, bool allow_identity = false);
    [[nodiscard]] Result<ExprPtr> traced_result(RuleId rule_id, ExprPtr before, ExprPtr after);
    void append_assumption(ExprPtr target);

    // Expression building helpers
    [[nodiscard]] ExprPtr make_sum_target(const std::vector<ExprPtr>& terms);
    [[nodiscard]] ExprPtr make_product_target(const std::vector<ExprPtr>& factors);
    void collect_additive_operands(ExprPtr expr, bool negate, std::vector<std::pair<ExprPtr, bool>>& operands);
    
    // Predicates and property checks
    [[nodiscard]] bool is_assumed_nonzero(ExprPtr expr) const;
    [[nodiscard]] bool is_known_positive(ExprPtr expr) const;
    [[nodiscard]] bool is_known_nonnegative(ExprPtr expr) const;
    [[nodiscard]] bool is_known_negative(ExprPtr expr) const;

    // Monomial helpers
    [[nodiscard]] static bool monomial_keys_equal(const MonomialKey& lhs, const MonomialKey& rhs);
    [[nodiscard]] Result<std::optional<MonomialTerm>> extract_monomial(ExprPtr expr);
    [[nodiscard]] ExprPtr build_monomial(const MonomialKey& key, const Rational& coefficient);
    static void merge_symbolic_factors(std::vector<std::pair<ExprPtr, BigInt>>& factors);

    AstArena& arena_;
    const Assumptions* assumptions_;
    const RewriteProvider* rewrite_provider_;
    ComputationTrace* trace_{nullptr};
    bool trace_enabled_{false};
    std::chrono::milliseconds timeout_{0};
    std::chrono::steady_clock::time_point* operation_started_at_{nullptr};
    std::uint64_t* ops_count_{nullptr};
    CASContext* context_{nullptr};
    std::vector<std::function<ExprPtr(ExprPtr)>> root_frames_;
};

} // namespace detail
} // namespace cas::symbolic
