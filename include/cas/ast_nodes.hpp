/// @file ast_nodes.hpp
/// @brief All concrete ExprNode subclasses: ExprNode base, ExprPtr smart handle,
///        SIDimensions, and all node types (IntegerLit … SeriesExp).
///
/// Split from ast.hpp (F8.0 / Task 1.1).
/// Depends on: ast_kinds.hpp, cas/bigint.hpp, cas/builtin_functions.hpp.
#pragma once

#include "cas/ast_kinds.hpp"
#include "cas/bigint.hpp"
#include "cas/builtin_functions.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cas {

// ---------------------------------------------------------------------------
// Base node and smart pointer
// ---------------------------------------------------------------------------

struct ExprNode {
    explicit constexpr ExprNode(ExprKind node_kind) noexcept : kind(node_kind) {}
    ExprNode(const ExprNode& other) noexcept : kind(other.kind), hash_cache_(0) {}
    ExprNode& operator=(const ExprNode& other) noexcept { kind = other.kind; hash_cache_.store(0, std::memory_order_relaxed); return *this; }
    ExprNode(ExprNode&& other) noexcept : kind(other.kind), hash_cache_(0) {}
    ExprNode& operator=(ExprNode&& other) noexcept { kind = other.kind; hash_cache_.store(0, std::memory_order_relaxed); return *this; }

    ExprKind kind;
    mutable std::atomic<std::size_t> hash_cache_{0};
};

class ExprPtr {
public:
    constexpr ExprPtr() noexcept = default;
    constexpr ExprPtr(std::nullptr_t) noexcept : node_(nullptr) {}
    explicit constexpr ExprPtr(const ExprNode* node) noexcept : node_(node) {}

    ExprPtr& operator=(std::nullptr_t) noexcept {
        node_ = nullptr;
        return *this;
    }

    [[nodiscard]] constexpr const ExprNode* get() const noexcept { return node_; }
    [[nodiscard]] constexpr const ExprNode& operator*() const noexcept { return *node_; }
    [[nodiscard]] constexpr const ExprNode* operator->() const noexcept { return node_; }
    [[nodiscard]] constexpr explicit operator bool() const noexcept { return node_ != nullptr; }

    [[nodiscard]] friend constexpr bool operator==(ExprPtr lhs, ExprPtr rhs) noexcept {
        return lhs.node_ == rhs.node_;
    }

    [[nodiscard]] friend constexpr bool operator!=(ExprPtr lhs, ExprPtr rhs) noexcept {
        return lhs.node_ != rhs.node_;
    }

    // WARNING (F7.0-A2.3): pointer-based ordering, NOT mathematically canonical.
    // ASLR randomizes addresses, so the ordering changes between runs.
    //
    // Use this operator ONLY for:
    //   - hash maps / unordered_set keys (key identity only)
    //   - internal set/map deduplication where iteration order is invisible
    //     to the user
    //
    // For any user-facing output (factor lists, root lists, simplified sums/
    // products, serialized AST, test goldens), sort with `canonical_compare`
    // — see `include/cas/symbolic.hpp:307`. Existing call sites verified:
    // normal_form.cpp:71,136; simplify_arithmetic.cpp:813; simplify_arithmetic_chain_*.cpp;
    // integrate_substitution.cpp ExprLess; gaussian_factor.cpp sorts by mathematical norm.
    [[nodiscard]] friend constexpr bool operator<(ExprPtr lhs, ExprPtr rhs) noexcept {
        return lhs.node_ < rhs.node_;
    }

private:
    const ExprNode* node_{nullptr};
};

// ---------------------------------------------------------------------------
// Physical-quantity dimension carrier
// ---------------------------------------------------------------------------

struct SIDimensions {
    int16_t m = 0;
    int16_t kg = 0;
    int16_t s = 0;
    int16_t A = 0;
    int16_t K = 0;
    int16_t mol = 0;
    int16_t cd = 0;

    [[nodiscard]] constexpr bool operator==(const SIDimensions& other) const noexcept {
        return m == other.m && kg == other.kg && s == other.s && A == other.A &&
               K == other.K && mol == other.mol && cd == other.cd;
    }

    [[nodiscard]] constexpr bool operator!=(const SIDimensions& other) const noexcept {
        return !(*this == other);
    }

    [[nodiscard]] constexpr bool operator<(const SIDimensions& other) const noexcept {
        if (m != other.m) return m < other.m;
        if (kg != other.kg) return kg < other.kg;
        if (s != other.s) return s < other.s;
        if (A != other.A) return A < other.A;
        if (K != other.K) return K < other.K;
        if (mol != other.mol) return mol < other.mol;
        return cd < other.cd;
    }

    [[nodiscard]] constexpr bool is_dimensionless() const noexcept {
        return m == 0 && kg == 0 && s == 0 && A == 0 && K == 0 && mol == 0 && cd == 0;
    }
};

// ---------------------------------------------------------------------------
// Concrete node types
// ---------------------------------------------------------------------------

struct Quantity : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Quantity;

    Quantity(ExprPtr val, SIDimensions dims) : ExprNode(KIND), value(val), dimensions(dims) {}

    ExprPtr value;
    SIDimensions dimensions;
};

struct IntegerLit : ExprNode {
    static constexpr ExprKind KIND = ExprKind::IntegerLit;

    explicit IntegerLit(BigInt literal_value) : ExprNode(KIND), value(std::move(literal_value)) {}

    BigInt value;
};

struct RationalLit : ExprNode {
    static constexpr ExprKind KIND = ExprKind::RationalLit;

    RationalLit(BigInt literal_numerator, BigInt literal_denominator)
        : ExprNode(KIND), numerator(std::move(literal_numerator)), denominator(std::move(literal_denominator)) {}

    BigInt numerator;
    BigInt denominator;
};

/// Front-end decimal literal preserved exactly as written by the user.
/// This keeps F1 lossless without introducing floating-point semantics in the symbolic core.
struct DecimalLit : ExprNode {
    static constexpr ExprKind KIND = ExprKind::DecimalLit;

    explicit DecimalLit(std::string literal_text) : ExprNode(KIND), text(std::move(literal_text)) {}
    explicit DecimalLit(double val) : ExprNode(KIND) {
        std::ostringstream oss;
        oss << std::setprecision(17) << val;
        text = oss.str();
    }

    [[nodiscard]] double to_double() const { return std::stod(text); }

    std::string text;
};

struct ComplexLit : ExprNode {
    static constexpr ExprKind KIND = ExprKind::ComplexLit;

    ComplexLit(BigInt re_num, BigInt re_den, BigInt im_num, BigInt im_den)
        : ExprNode(KIND), re_num(std::move(re_num)), re_den(std::move(re_den)),
          im_num(std::move(im_num)), im_den(std::move(im_den)) {}

    BigInt re_num;
    BigInt re_den;
    BigInt im_num;
    BigInt im_den;
};

struct Symbol : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Symbol;

    explicit Symbol(std::string symbol_name) : ExprNode(KIND), name(std::move(symbol_name)) {}

    std::string name;
};

struct Constant : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Constant;

    explicit constexpr Constant(MathConstant constant_value) noexcept : ExprNode(KIND), value(constant_value) {}

    MathConstant value;
};

struct Unary : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Unary;

    constexpr Unary(UnaryOp unary_op, ExprPtr unary_operand) noexcept
        : ExprNode(KIND), op(unary_op), operand(unary_operand) {}

    UnaryOp op;
    ExprPtr operand;
};

struct Binary : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Binary;

    constexpr Binary(BinaryOp binary_op, ExprPtr binary_left, ExprPtr binary_right) noexcept
        : ExprNode(KIND), op(binary_op), left(binary_left), right(binary_right) {}

    BinaryOp op;
    ExprPtr left;
    ExprPtr right;
};

struct FuncCall : ExprNode {
    static constexpr ExprKind KIND = ExprKind::FuncCall;

    FuncCall(std::string function_name, std::vector<ExprPtr> function_args)
        : ExprNode(KIND), name(std::move(function_name)), func_id(get_builtin_op(name)), args(std::move(function_args)) {}

    FuncCall(BuiltinOp function_id, std::vector<ExprPtr> function_args)
        : ExprNode(KIND), name(std::string(builtin_op_name(function_id))), func_id(function_id), args(std::move(function_args)) {}

    std::string name;
    BuiltinOp func_id;
    std::vector<ExprPtr> args;
};

struct Sum : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Sum;

    explicit Sum(std::vector<ExprPtr> sum_terms) : ExprNode(KIND), terms(std::move(sum_terms)) {}

    std::vector<ExprPtr> terms;
};

struct Product : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Product;

    explicit Product(std::vector<ExprPtr> product_factors) : ExprNode(KIND), factors(std::move(product_factors)) {}

    std::vector<ExprPtr> factors;
};

struct Integral : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Integral;

    Integral(ExprPtr integral_integrand, Symbol integral_variable, std::optional<ExprPtr> integral_lower, std::optional<ExprPtr> integral_upper)
        : ExprNode(KIND),
          integrand(integral_integrand),
          variable(std::move(integral_variable)),
          lower(std::move(integral_lower)),
          upper(std::move(integral_upper)) {}

    ExprPtr integrand;
    Symbol variable;
    std::optional<ExprPtr> lower;
    std::optional<ExprPtr> upper;
};

struct Derivative : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Derivative;

    Derivative(ExprPtr derivative_expression, Symbol derivative_variable, unsigned int derivative_order = 1U)
        : ExprNode(KIND),
          expression(derivative_expression),
          variable(std::move(derivative_variable)),
          order(derivative_order) {}

    ExprPtr expression;
    Symbol variable;
    unsigned int order{1U};
};

struct Limit : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Limit;

    Limit(ExprPtr limit_expression, Symbol limit_variable, ExprPtr limit_point, LimitDirection limit_direction = LimitDirection::Both)
        : ExprNode(KIND),
          expression(limit_expression),
          variable(std::move(limit_variable)),
          point(limit_point),
          direction(limit_direction) {}

    ExprPtr expression;
    Symbol variable;
    ExprPtr point;
    LimitDirection direction{LimitDirection::Both};
};

// Rational isolating bound — used inside RootOf to pin a single real root of
// the minimal polynomial. Storing as (num,den) BigInt pairs avoids an
// include-cycle with rational.hpp (Rational forward-decl insufficient since
// the struct holds the value by-value, not by-pointer).
struct IsolatingBound {
    BigInt low_num;
    BigInt low_den;
    BigInt high_num;
    BigInt high_den;
};

struct RootOf : ExprNode {
    static constexpr ExprKind KIND = ExprKind::RootOf;

    // Legacy constructor — root_index references the i-th real root in the
    // canonical Sturm-sorted order. Kept for backward compatibility.
    RootOf(ExprPtr root_polynomial, Symbol root_variable,
           std::optional<std::size_t> root_position = std::nullopt)
        : ExprNode(KIND),
          polynomial(root_polynomial),
          variable(std::move(root_variable)),
          root_index(std::move(root_position)) {}

    // Rigorous constructor — bound is a half-open rational interval
    // [low, high] guaranteed to contain exactly one real root of `polynomial`.
    // The Sturm sequence is run at construction time to verify isolation;
    // see make_rootof_isolated() in algebra.hpp.
    RootOf(ExprPtr root_polynomial, Symbol root_variable,
           IsolatingBound bound,
           std::optional<std::size_t> root_position = std::nullopt)
        : ExprNode(KIND),
          polynomial(root_polynomial),
          variable(std::move(root_variable)),
          root_index(std::move(root_position)),
          isolating_bound(std::move(bound)) {}

    ExprPtr polynomial;
    Symbol variable;
    std::optional<std::size_t> root_index;
    // F8.0-5.4: rational isolating interval for the represented real root.
    // When present, equality / ordering between RootOf instances can be
    // decided without re-running Sturm.
    std::optional<IsolatingBound> isolating_bound;
};

struct Matrix : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Matrix;

    Matrix(std::size_t matrix_rows, std::size_t matrix_cols, std::vector<ExprPtr> matrix_elements)
        : ExprNode(KIND), rows(matrix_rows), cols(matrix_cols), elements(std::move(matrix_elements)) {}

    std::size_t rows{0U};
    std::size_t cols{0U};
    std::vector<ExprPtr> elements;
};

struct SeriesExp : ExprNode {
    static constexpr ExprKind KIND = ExprKind::SeriesExp;

    SeriesExp(Symbol series_var, ExprPtr series_point, std::vector<std::pair<long long, ExprPtr>> series_terms, long long series_order)
        : ExprNode(KIND), var(std::move(series_var)), point(series_point), terms(std::move(series_terms)), order(series_order) {}

    Symbol var;
    ExprPtr point;
    std::vector<std::pair<long long, ExprPtr>> terms; // (exp, coeff)
    long long order;
};

// ---------------------------------------------------------------------------
// Helpers: equality, hash, cast utilities
// ---------------------------------------------------------------------------

[[nodiscard]] inline bool is_equal(ExprPtr lhs, ExprPtr rhs) noexcept {
    return lhs == rhs;
}

[[nodiscard]] ExprKind expr_kind(ExprPtr expr) noexcept;
[[nodiscard]] bool structural_equal(ExprPtr lhs, ExprPtr rhs) noexcept;
[[nodiscard]] std::size_t expr_hash(ExprPtr expr) noexcept;
[[nodiscard]] std::string_view expr_kind_name(ExprKind kind) noexcept;

struct ExprHash {
    std::size_t operator()(ExprPtr expr) const noexcept { return expr_hash(expr); }
};

struct ExprEqual {
    bool operator()(ExprPtr lhs, ExprPtr rhs) const noexcept { return structural_equal(lhs, rhs); }
};

template <typename T>
[[nodiscard]] constexpr bool expr_is(ExprPtr expr) noexcept {
    static_assert(std::is_base_of_v<ExprNode, T>, "expr_is requires an ExprNode-derived type");
    return expr && expr->kind == T::KIND;
}

template <typename T>
[[nodiscard]] constexpr const T* expr_cast(ExprPtr expr) noexcept {
    static_assert(std::is_base_of_v<ExprNode, T>, "expr_cast requires an ExprNode-derived type");
    return expr_is<T>(expr) ? static_cast<const T*>(expr.get()) : nullptr;
}

template <typename T>
[[nodiscard]] constexpr const T& expr_ref(ExprPtr expr) noexcept {
    static_assert(std::is_base_of_v<ExprNode, T>, "expr_ref requires an ExprNode-derived type");
    return *static_cast<const T*>(expr.get());
}

[[nodiscard]] ExprPtr clone_into_arena(ExprPtr expr, class AstArena& target, std::unordered_map<ExprPtr, ExprPtr>& cache);

}  // namespace cas

namespace std {
template <>
struct hash<cas::ExprPtr> {
    std::size_t operator()(cas::ExprPtr expr) const noexcept {
        return cas::expr_hash(expr);
    }
};
}  // namespace std
