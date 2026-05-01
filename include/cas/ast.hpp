#pragma once

#include "cas/bigint.hpp"
#include "cas/builtin_functions.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <array>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <iomanip>
#include <sstream>

namespace cas {

enum class ExprKind : std::uint8_t {
    Null,
    IntegerLit,
    RationalLit,
    DecimalLit,
    Symbol,
    Constant,
    Unary,
    Binary,
    FuncCall,
    Sum,
    Product,
    Integral,
    Derivative,
    Limit,
    RootOf,
    Matrix,
};

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

    [[nodiscard]] friend constexpr bool operator<(ExprPtr lhs, ExprPtr rhs) noexcept {
        return lhs.node_ < rhs.node_;
    }

private:
    const ExprNode* node_{nullptr};
};

[[nodiscard]] inline bool is_equal(ExprPtr lhs, ExprPtr rhs) noexcept {
    return lhs == rhs;
}

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

struct Symbol : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Symbol;

    explicit Symbol(std::string symbol_name) : ExprNode(KIND), name(std::move(symbol_name)) {}

    std::string name;
};

enum class MathConstant : std::uint8_t {
    Pi,
    E,
    I,
    Infinity,
    NaN,
};

struct Constant : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Constant;

    explicit constexpr Constant(MathConstant constant_value) noexcept : ExprNode(KIND), value(constant_value) {}

    MathConstant value;
};

enum class UnaryOp : std::uint8_t {
    Neg,
    Factorial,
};

struct Unary : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Unary;

    constexpr Unary(UnaryOp unary_op, ExprPtr unary_operand) noexcept
        : ExprNode(KIND), op(unary_op), operand(unary_operand) {}

    UnaryOp op;
    ExprPtr operand;
};

enum class BinaryOp : std::uint8_t {
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Mod,
    Equal,
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

enum class LimitDirection : std::uint8_t {
    Left,
    Right,
    Both,
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

struct RootOf : ExprNode {
    static constexpr ExprKind KIND = ExprKind::RootOf;

    RootOf(ExprPtr root_polynomial, Symbol root_variable, std::optional<std::size_t> root_position = std::nullopt)
        : ExprNode(KIND),
          polynomial(root_polynomial),
          variable(std::move(root_variable)),
          root_index(std::move(root_position)) {}

    ExprPtr polynomial;
    Symbol variable;
    std::optional<std::size_t> root_index;
};

struct Matrix : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Matrix;

    Matrix(std::size_t matrix_rows, std::size_t matrix_cols, std::vector<ExprPtr> matrix_elements)
        : ExprNode(KIND), rows(matrix_rows), cols(matrix_cols), elements(std::move(matrix_elements)) {}

    std::size_t rows{0U};
    std::size_t cols{0U};
    std::vector<ExprPtr> elements;
};

[[nodiscard]] ExprKind expr_kind(ExprPtr expr) noexcept;
[[nodiscard]] bool structural_equal(ExprPtr lhs, ExprPtr rhs) noexcept;
[[nodiscard]] std::size_t expr_hash(ExprPtr expr) noexcept;

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

class AstArena {
public:
    static constexpr std::uint32_t API_VERSION = 1;
    static constexpr std::size_t DEFAULT_BLOCK_BYTES = 64U * 1024U;

    AstArena() = default;
    ~AstArena();
    AstArena(const AstArena&) = delete;
    AstArena& operator=(const AstArena&) = delete;
    AstArena(AstArena&& other) noexcept;
    AstArena& operator=(AstArena&& other) noexcept;

    template <typename T, typename... Args>
    [[nodiscard]] ExprPtr make(Args&&... args) {
        static_assert(std::is_base_of_v<ExprNode, T>, "AstArena::make requires an ExprNode-derived type");

        T candidate(std::forward<Args>(args)...);

        std::lock_guard<std::mutex> lock(mutex_);

        // Path veloci per costanti calde (evitano lookup se già presenti)
        if constexpr (std::is_same_v<T, Constant>) {
            const std::size_t index = static_cast<std::size_t>(candidate.value);
            if (index < interned_constants_.size() && interned_constants_[index]) return interned_constants_[index];
        } else if constexpr (std::is_same_v<T, IntegerLit>) {
            if (candidate.value.is_zero() && interned_zero_) return *interned_zero_;
            if (candidate.value == BigInt(1) && interned_one_) return *interned_one_;
            if (candidate.value == BigInt(-1) && interned_negative_one_) return *interned_negative_one_;
        }

        ExprPtr candidate_ptr(&candidate);
        if (auto it = interning_table_.find(candidate_ptr); it != interning_table_.end()) {
            return *it;
        }

        // Se non trovato, allochiamo realmente in arena
        ExprPtr created = make_uncached<T>(std::move(candidate));

        // Aggiorna cache calde se necessario
        if constexpr (std::is_same_v<T, Constant>) {
            const std::size_t index = static_cast<std::size_t>(expr_ref<Constant>(created).value);
            if (index < interned_constants_.size()) interned_constants_[index] = created;
        } else if constexpr (std::is_same_v<T, IntegerLit>) {
            const auto& val = expr_ref<IntegerLit>(created).value;
            if (val.is_zero()) interned_zero_ = created;
            else if (val == BigInt(1)) interned_one_ = created;
            else if (val == BigInt(-1)) interned_negative_one_ = created;
        }

        interning_table_.insert(created);
        return created;
    }

    [[nodiscard]] std::size_t size() const noexcept;

private:
    template <typename T, typename... Args>
    [[nodiscard]] ExprPtr make_uncached(Args&&... args) {
        T* node = static_cast<T*>(allocate(sizeof(T), alignof(T)));
        node = ::new (node) T(std::forward<Args>(args)...);

        if (node_chunks_.empty() || node_chunks_.back().size() == CHUNK_SIZE) {
            node_chunks_.emplace_back();
            node_chunks_.back().reserve(CHUNK_SIZE);
        }
        node_chunks_.back().push_back(node);
        ++total_nodes_;

        return ExprPtr(node);
    }
    struct Block {
        std::unique_ptr<std::byte[]> data;
        std::size_t capacity{0U};
        std::size_t used{0U};
    };

    static constexpr std::size_t CHUNK_SIZE = 8192U;

    [[nodiscard]] void* allocate(std::size_t size, std::size_t alignment);
    void append_block(std::size_t minimum_bytes);

    mutable std::mutex mutex_;
    std::vector<Block> blocks_;
    std::vector<std::vector<ExprNode*>> node_chunks_;
    std::size_t total_nodes_{0U};

    // Tabella di interning universale
    std::unordered_set<ExprPtr, ExprHash, ExprEqual> interning_table_;

    // Cache per costanti matematiche calde (indicizzate per MathConstant)
    std::array<ExprPtr, 5U> interned_constants_{};
    // Cache per literal interi comuni
    std::optional<ExprPtr> interned_zero_{};
    std::optional<ExprPtr> interned_one_{};
    std::optional<ExprPtr> interned_negative_one_{};
};

template <typename Visitor>
decltype(auto) visit_expr(ExprPtr expr, Visitor&& visitor) {
    switch (expr ? expr->kind : ExprKind::Null) {
    case ExprKind::IntegerLit:
        return std::forward<Visitor>(visitor)(expr_ref<IntegerLit>(expr));
    case ExprKind::RationalLit:
        return std::forward<Visitor>(visitor)(expr_ref<RationalLit>(expr));
    case ExprKind::DecimalLit:
        return std::forward<Visitor>(visitor)(expr_ref<DecimalLit>(expr));
    case ExprKind::Symbol:
        return std::forward<Visitor>(visitor)(expr_ref<Symbol>(expr));
    case ExprKind::Constant:
        return std::forward<Visitor>(visitor)(expr_ref<Constant>(expr));
    case ExprKind::Unary:
        return std::forward<Visitor>(visitor)(expr_ref<Unary>(expr));
    case ExprKind::Binary:
        return std::forward<Visitor>(visitor)(expr_ref<Binary>(expr));
    case ExprKind::FuncCall:
        return std::forward<Visitor>(visitor)(expr_ref<FuncCall>(expr));
    case ExprKind::Sum:
        return std::forward<Visitor>(visitor)(expr_ref<Sum>(expr));
    case ExprKind::Product:
        return std::forward<Visitor>(visitor)(expr_ref<Product>(expr));
    case ExprKind::Integral:
        return std::forward<Visitor>(visitor)(expr_ref<Integral>(expr));
    case ExprKind::Derivative:
        return std::forward<Visitor>(visitor)(expr_ref<Derivative>(expr));
    case ExprKind::Limit:
        return std::forward<Visitor>(visitor)(expr_ref<Limit>(expr));
    case ExprKind::RootOf:
        return std::forward<Visitor>(visitor)(expr_ref<RootOf>(expr));
    case ExprKind::Matrix:
        return std::forward<Visitor>(visitor)(expr_ref<Matrix>(expr));
    case ExprKind::Null:
        std::abort();
    }

    std::abort();
}

[[nodiscard]] ExprKind expr_kind(ExprPtr expr) noexcept;
[[nodiscard]] bool structural_equal(ExprPtr lhs, ExprPtr rhs) noexcept;
[[nodiscard]] std::size_t expr_hash(ExprPtr expr) noexcept;
[[nodiscard]] std::string_view expr_kind_name(ExprKind kind) noexcept;

[[nodiscard]] ExprPtr clone_into_arena(ExprPtr expr, AstArena& target, std::unordered_map<ExprPtr, ExprPtr>& cache);

}  // namespace cas

namespace std {
template <>
struct hash<cas::ExprPtr> {
    std::size_t operator()(cas::ExprPtr expr) const noexcept {
        return cas::expr_hash(expr);
    }
};
} // namespace std
