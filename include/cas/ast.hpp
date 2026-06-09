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
    ComplexLit,

    Matrix,
    SeriesExp,
    Quantity,
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

struct Quantity : ExprNode {
    static constexpr ExprKind KIND = ExprKind::Quantity;

    Quantity(ExprPtr val, SIDimensions dims) : ExprNode(KIND), value(val), dimensions(dims) {}

    ExprPtr value;
    SIDimensions dimensions;
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

enum class MathConstant : std::uint8_t {
    Pi,
    E,
    I,
    Infinity,
    NaN,
    EulerGamma,
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
    Less,
    Greater,
    LessEqual,
    GreaterEqual,
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

struct SeriesExp : ExprNode {
    static constexpr ExprKind KIND = ExprKind::SeriesExp;

    SeriesExp(Symbol series_var, ExprPtr series_point, std::vector<std::pair<long long, ExprPtr>> series_terms, long long series_order)
        : ExprNode(KIND), var(std::move(series_var)), point(series_point), terms(std::move(series_terms)), order(series_order) {}

    Symbol var;
    ExprPtr point;
    std::vector<std::pair<long long, ExprPtr>> terms; // (exp, coeff)
    long long order;
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

    // F1.3-NEW: number of interning shards. Power-of-2 for cheap modulo.
    static constexpr std::size_t N_INTERN_SHARDS = 16U;

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

        // F1.3-NEW: fast-path for hot constants (no shard lookup needed).
        // These are protected by the alloc_mutex_ when written; reads after
        // first write are safe because ExprPtr is a raw-pointer wrapper —
        // we check under the shard lock below if the fast-path misses.
        if constexpr (std::is_same_v<T, Constant>) {
            const std::size_t idx = static_cast<std::size_t>(candidate.value);
            if (idx < interned_constants_.size()) {
                // Load with acquire so we see any previously stored value.
                ExprPtr cached = interned_constants_[idx].load(std::memory_order_acquire);
                if (cached) return cached;
            }
        } else if constexpr (std::is_same_v<T, IntegerLit>) {
            if (candidate.value.is_zero()) {
                ExprPtr cached = interned_zero_.load(std::memory_order_acquire);
                if (cached) return cached;
            } else if (candidate.value == BigInt(1)) {
                ExprPtr cached = interned_one_.load(std::memory_order_acquire);
                if (cached) return cached;
            } else if (candidate.value == BigInt(-1)) {
                ExprPtr cached = interned_neg_one_.load(std::memory_order_acquire);
                if (cached) return cached;
            }
        }

        // F1.3-NEW: compute shard from structural hash of the candidate node.
        ExprPtr candidate_ptr(&candidate);
        const std::size_t h = expr_hash(candidate_ptr);
        const std::size_t shard_idx = h & (N_INTERN_SHARDS - 1U);

        std::lock_guard<std::mutex> shard_lock(intern_shards_[shard_idx]);

        // Each shard owns its own unordered_set — no shared container across
        // shard locks (eliminates the data-race on a single interning_table_).
        auto& shard_table = intern_shard_tables_[shard_idx];
        if (auto it = shard_table.find(candidate_ptr); it != shard_table.end()) {
            return *it;
        }

        // Not found: allocate under the global alloc mutex (nested inside the
        // shard lock — always in shard → alloc order, never reversed).
        ExprPtr created;
        {
            std::lock_guard<std::mutex> alloc_lock(alloc_mutex_);
            created = make_uncached<T>(std::move(candidate));
        }

        // Update atomic hot caches (visible to fast-path without lock).
        if constexpr (std::is_same_v<T, Constant>) {
            const std::size_t idx = static_cast<std::size_t>(expr_ref<Constant>(created).value);
            if (idx < interned_constants_.size()) {
                interned_constants_[idx].store(created, std::memory_order_release);
            }
        } else if constexpr (std::is_same_v<T, IntegerLit>) {
            const auto& val = expr_ref<IntegerLit>(created).value;
            if (val.is_zero()) interned_zero_.store(created, std::memory_order_release);
            else if (val == BigInt(1)) interned_one_.store(created, std::memory_order_release);
            else if (val == BigInt(-1)) interned_neg_one_.store(created, std::memory_order_release);
        }

        shard_table.insert(created);
        return created;
    }

    [[nodiscard]] std::size_t size() const noexcept;

    // F7.0-A3.1: controlled reset for long-lived REPL/server processes.
    // Destroys every node ever allocated, clears interning tables, releases
    // all memory blocks, resets hot caches. Required to prevent unbounded
    // AstArena growth in a REPL/web-server loop where each query allocates
    // transient intermediate nodes.
    //
    // ⚠ INVALIDATION: every `ExprPtr` previously vended by this arena becomes
    // a dangling pointer after `reset()`. The caller is responsible for
    // discarding all stale `ExprPtr` values BEFORE calling reset.
    //
    // Thread-safety: acquires every shard mutex AND the alloc mutex (shard
    // → alloc order). Concurrent `make<T>` calls will block until reset
    // completes. NEVER call `reset()` while another thread is using the arena.
    //
    // Root migration (deep-copy preserve) is NOT provided here — see
    // `migrate_into(...)` follow-up task (HC-F70-A31-MIGRATION-TODO).
    void reset();

private:
    template <typename T, typename... Args>
    [[nodiscard]] ExprPtr make_uncached(Args&&... args) {
        // Caller holds alloc_mutex_.
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

    // F1.3-NEW: bump allocator state protected by its own mutex.
    mutable std::mutex alloc_mutex_;
    std::vector<Block> blocks_;
    std::vector<std::vector<ExprNode*>> node_chunks_;
    std::size_t total_nodes_{0U};

    // F1.3-NEW / HPP-016: interning table sharded by hash % N_INTERN_SHARDS.
    //
    // N_INTERN_SHARDS = 16: power-of-2 so shard selection is a cheap bitwise
    // AND.  16 shards reduce lock contention to ~1/16 of a global lock in the
    // common case (distinct hash buckets map to distinct shards).
    //
    // THREAD-SAFETY INVARIANT: each intern_shard_tables_[i] is accessed
    // exclusively while holding intern_shards_[i].  No operation may touch
    // intern_shard_tables_[i] without first acquiring intern_shards_[i].
    // The alloc_mutex_ is a NESTED lock (always acquired AFTER a shard lock,
    // never before) — locking order: shard → alloc, never alloc → shard.
    mutable std::array<std::mutex, N_INTERN_SHARDS> intern_shards_;
    // Per-shard interning sets: intern_shard_tables_[i] is protected by
    // intern_shards_[i].  Using per-shard containers eliminates the data
    // race that would arise from a single shared unordered_set accessed under
    // different shard locks.
    std::array<std::unordered_set<ExprPtr, ExprHash, ExprEqual>, N_INTERN_SHARDS>
        intern_shard_tables_;

    // Compile-time check: N_INTERN_SHARDS must be a power of two so that
    // (hash & (N_INTERN_SHARDS - 1)) == (hash % N_INTERN_SHARDS).
    static_assert((N_INTERN_SHARDS & (N_INTERN_SHARDS - 1U)) == 0U,
                  "N_INTERN_SHARDS must be a power of two");

    // F1.3-NEW: hot-constant caches as atomics for lock-free fast-path reads.
    std::array<std::atomic<ExprPtr>, 5U> interned_constants_{};
    std::atomic<ExprPtr> interned_zero_{ExprPtr{}};
    std::atomic<ExprPtr> interned_one_{ExprPtr{}};
    std::atomic<ExprPtr> interned_neg_one_{ExprPtr{}};
};

template <typename Visitor>
decltype(auto) visit_expr(ExprPtr expr, Visitor&& visitor) {
    switch (expr ? expr->kind : ExprKind::Null) {
    case ExprKind::IntegerLit:
        return std::forward<Visitor>(visitor)(expr_ref<IntegerLit>(expr));
    case ExprKind::RationalLit:
        return std::forward<Visitor>(visitor)(expr_ref<RationalLit>(expr));
    case ExprKind::ComplexLit:
        return std::forward<Visitor>(visitor)(expr_ref<ComplexLit>(expr));

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
    case ExprKind::SeriesExp:
        return std::forward<Visitor>(visitor)(expr_ref<SeriesExp>(expr));
    case ExprKind::Quantity:
        return std::forward<Visitor>(visitor)(expr_ref<Quantity>(expr));
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
