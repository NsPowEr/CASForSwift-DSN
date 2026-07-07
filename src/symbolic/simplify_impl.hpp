#pragma once

#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "cas/complex_rational.hpp"
#include <chrono>
#include <functional>
#include <string>

#include <vector>
#include <optional>

namespace cas::symbolic {

// Helper functions for simplification logic
namespace detail {

thread_local extern int simplification_depth;
constexpr int MAX_SIMPLIFICATION_DEPTH = 300;

struct DepthGuard {
    explicit DepthGuard(int max_depth = MAX_SIMPLIFICATION_DEPTH);
    ~DepthGuard();
    bool exceeded() const;
private:
    int max_depth_;
};

// F7.0-A3.2: async-aware depth propagation primitives.
// AsyncDepthScope propagates the parent thread's recursion budget to a worker
// thread if async simplification is ever introduced.
[[nodiscard]] int current_simplify_depth() noexcept;

class AsyncDepthScope {
public:
    explicit AsyncDepthScope(int inherited_depth) noexcept;
    ~AsyncDepthScope() noexcept;
    AsyncDepthScope(const AsyncDepthScope&)            = delete;
    AsyncDepthScope& operator=(const AsyncDepthScope&) = delete;
private:
    int prev_;
};

// Path-based cycle guard (A20).  Tracks the current recursion *path* (the
// direct ancestor chain) rather than the global set of all active nodes.
// A node is a cycle only if it appears as a direct ancestor on the current
// descent — i.e. the same pointer is already being simplified in an enclosing
// stack frame on this thread.  Structurally-shared nodes (same ExprPtr
// appearing as a valid child in multiple parents) are NOT false-positives
// because they are inserted+erased independently per call frame.
// Complexity: O(depth) linear scan, depth ≤ MAX_SIMPLIFICATION_DEPTH.
extern thread_local std::vector<ExprPtr> simplify_ancestor_path;

struct CycleGuard {
    ExprPtr expr_;
    bool cycle_{false};

    explicit CycleGuard(ExprPtr expr) : expr_(expr) {
        for (ExprPtr anc : simplify_ancestor_path) {
            if (anc == expr) { cycle_ = true; return; }
        }
        simplify_ancestor_path.push_back(expr);
    }
    ~CycleGuard() {
        if (!cycle_) {
            simplify_ancestor_path.pop_back();
        }
    }
    [[nodiscard]] bool cycle_detected() const noexcept { return cycle_; }
};


struct LiteralRational {
    Rational value;
    bool exact{false};
};

struct LiteralComplex {
    ComplexRational value;
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
[[nodiscard]] ExprPtr make_complex(AstArena& arena, const ComplexRational& value);
[[nodiscard]] bool is_zero_expr(ExprPtr expr);
[[nodiscard]] bool is_one_expr(ExprPtr expr);
[[nodiscard]] bool is_constant_expr(ExprPtr expr, MathConstant constant);

// F7.5.F1 Phase 2 — extended-real arithmetic helpers.
// Return std::nullopt when no operand is extended-real (±∞, ComplexInfinity,
// Indeterminate). Return the propagated extended-real result otherwise.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Extended_Real_AST.md
[[nodiscard]] std::optional<ExprPtr> try_simplify_sum_extended_real(
    const std::vector<ExprPtr>& terms, AstArena& arena);
[[nodiscard]] std::optional<ExprPtr> try_simplify_product_extended_real(
    const std::vector<ExprPtr>& factors, AstArena& arena);
[[nodiscard]] std::optional<ExprPtr> try_simplify_pow_extended_real(
    ExprPtr base, ExprPtr exponent, AstArena& arena);
[[nodiscard]] Result<bool> try_get_exact_rational(ExprPtr expr, LiteralRational& out);
[[nodiscard]] Result<bool> try_get_exact_complex(ExprPtr expr, LiteralComplex& out);
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
    [[nodiscard]] Result<ExprPtr> simplify_node(ExprPtr original, const ComplexLit& node);
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
    // Chebyshev/DeMoivre linearization of sin^n / cos^n (integer n ≥ 2).
    // Returns std::nullopt when the rule does not apply (base not sin/cos, n out
    // of range, etc.) so the caller falls through. Defined in
    // simplify_arithmetic_power_trig.cpp.
    [[nodiscard]] std::optional<Result<ExprPtr>> try_linearize_trig_power(ExprPtr base, const BigInt& n);
    [[nodiscard]] Result<ExprPtr> simplify_sum_terms(const std::vector<ExprPtr>& terms, ExprPtr target_before = ExprPtr{}, bool inputs_are_simplified = false);
    // T-054b: combine sum terms over a common denominator and cancel. Non-recursive
    // (local univariate Rational-polynomial arithmetic, no simplify_expr). Rewrites
    // `terms` in place and returns true when a denominator fully cancels; otherwise
    // leaves the sum untouched. Defined in simplify_combine_fractions.cpp.
    [[nodiscard]] bool try_combine_common_denominator(std::vector<ExprPtr>& terms);
    [[nodiscard]] ExprPtr poly_to_expr_for_combine(const std::vector<Rational>& coeffs, ExprPtr var);
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

    // FuncCall domain sub-dispatchers (called by simplify_node(FuncCall))
    [[nodiscard]] Result<ExprPtr> simplify_funcall_trig(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    [[nodiscard]] Result<ExprPtr> simplify_funcall_arc_trig(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    [[nodiscard]] Result<ExprPtr> simplify_funcall_exp_log_sqrt(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    // W9.3 split: sqrt branch extracted from simplify_funcall_exp_log_sqrt to keep
    // both translation units under the 500-line anti-monolith limit.
    [[nodiscard]] Result<ExprPtr> simplify_funcall_sqrt(ExprPtr original, std::vector<ExprPtr> args, ExprPtr target_before);
    [[nodiscard]] Result<ExprPtr> simplify_funcall_special(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    // F7.5.E1: hypergeometric/elliptic split out from simplify_funcall_special
    // to keep simplify_special_fn.cpp under the 500-line anti-monolith limit.
    [[nodiscard]] Result<ExprPtr> simplify_funcall_hyper_elliptic(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    [[nodiscard]] Result<ExprPtr> simplify_funcall_bessel_orthogonal(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    // F7.5.E2: half-integer Bessel reduction via three-term recurrence from
    // ±1/2 base values. Returns Unimplemented when |p_num|/2 exceeds the
    // configured bound (ctx.max_bessel_half_integer_order()).
    [[nodiscard]] Result<ExprPtr> bessel_half_integer_reduce(BuiltinOp op, const BigInt& p_num, ExprPtr x_arg);
    // F7.5.E2: orthogonal polynomial recurrences (Chebyshev, Hermite,
    // Legendre, Laguerre, Jacobi) split out from simplify_bessel_orthogonal.cpp
    // to keep that file under the 500-line anti-monolith limit.
    [[nodiscard]] Result<ExprPtr> simplify_funcall_orthogonal_polys(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    [[nodiscard]] Result<ExprPtr> simplify_funcall_complex(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    // Exact-value identities for hyperbolic builtins (sinh/cosh/tanh/coth).
    [[nodiscard]] Result<ExprPtr> simplify_funcall_hyperbolic(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);
    // Combinatorial + erfc exact-value identities.
    [[nodiscard]] Result<ExprPtr> simplify_funcall_combinatorial(ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before);

    // Monomial helpers
    [[nodiscard]] static bool monomial_keys_equal(const MonomialKey& lhs, const MonomialKey& rhs);
    [[nodiscard]] Result<std::optional<MonomialTerm>> extract_monomial(ExprPtr expr);
    [[nodiscard]] ExprPtr build_monomial(const MonomialKey& key, const Rational& coefficient);
    static void merge_symbolic_factors(std::vector<std::pair<ExprPtr, BigInt>>& factors);

    // L3-04 Gamma reflection identity applied to symbolic factor pairs.
    // Implementation in simplify_arithmetic_chain_gamma.cpp.
    [[nodiscard]] Result<void> apply_gamma_reflection_pairs(
        std::vector<std::pair<ExprPtr, BigInt>>& symbolic,
        ComplexRational& coefficient);

    // sqrt(a)·sqrt(a) → a and sqrt(a)·sqrt(b) → sqrt(ab) collapses.
    // Implementation in simplify_arithmetic_chain_sqrt.cpp.
    [[nodiscard]] Result<void> collapse_sqrt_pairs(
        std::vector<std::pair<ExprPtr, BigInt>>& symbolic,
        ComplexRational& coefficient);

    // Risch IBP exp-fold.
    // Implementation in simplify_arithmetic_chain_exp.cpp.
    [[nodiscard]] Result<void> fold_exponential_products(
        std::vector<std::pair<ExprPtr, BigInt>>& symbolic,
        ComplexRational& coefficient);

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
