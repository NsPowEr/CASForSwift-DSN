/// @file ast_arena.hpp
/// @brief AstArena: the thread-safe, bump-allocating, interning arena for AST nodes.
///        Also includes visit_expr<Visitor> dispatch template.
///
/// Split from ast.hpp (F8.0 / Task 1.1).
/// Depends on: ast_nodes.hpp (ExprNode, ExprPtr, all node types).
#pragma once

#include "cas/ast_nodes.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cas {

// ---------------------------------------------------------------------------
// AstArena — thread-safe interning allocator
// ---------------------------------------------------------------------------

class AstArena {
public:
    static constexpr std::uint32_t API_VERSION = 1;
    static constexpr std::size_t DEFAULT_BLOCK_BYTES = 64U * 1024U;

    explicit AstArena(std::size_t num_shards = 0);
    ~AstArena();
    void configure_shards(std::size_t num_shards);
    [[nodiscard]] std::size_t num_shards() const noexcept { return num_shards_; }
    AstArena(const AstArena&) = delete;
    AstArena& operator=(const AstArena&) = delete;
    AstArena(AstArena&& other) noexcept;
    AstArena& operator=(AstArena&& other) noexcept;

    template <typename T, typename... Args>
    [[nodiscard]] ExprPtr make(Args&&... args) {
        static_assert(std::is_base_of_v<ExprNode, T>, "AstArena::make requires an ExprNode-derived type");

        T candidate(std::forward<Args>(args)...);

        // F1.3-NEW: fast-path for hot constants (no shard lookup needed).
        if constexpr (std::is_same_v<T, Constant>) {
            const std::size_t idx = static_cast<std::size_t>(candidate.value);
            if (idx < interned_constants_.size()) {
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
        const std::size_t shard_idx = h & (num_shards_ - 1U);

        std::lock_guard<std::mutex> shard_lock(intern_shards_[shard_idx]);

        auto& shard_table = intern_shard_tables_[shard_idx];
        if (auto it = shard_table.find(candidate_ptr); it != shard_table.end()) {
            return *it;
        }

        // Not found: allocate under the global alloc mutex.
        ExprPtr created;
        {
            std::lock_guard<std::mutex> alloc_lock(alloc_mutex_);
            created = make_uncached<T>(std::move(candidate));
        }

        if (!created || hash_dos_detected_) {
            return ExprPtr{};
        }

        // Update atomic hot caches.
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
        // F7.0-A3.7: detect pathological hash-collision chains.
        const std::size_t bucket = shard_table.bucket(created);
        if (shard_table.bucket_size(bucket) > MAX_COLLISION_CHAIN) {
            hash_dos_detected_ = true;
        }
        return created;
    }

    [[nodiscard]] std::size_t size() const noexcept;

    // F7.0-A3.5: total bytes allocated by append_block since last reset.
    [[nodiscard]] std::size_t bytes_allocated() const noexcept;

    // F7.0-A3.7: Hash-DoS defence.
    [[nodiscard]] bool hash_dos_detected() const noexcept;
    void clear_hash_dos_flag() noexcept;
    static constexpr std::size_t MAX_COLLISION_CHAIN = 128U;

    void set_max_memory_budget_bytes(std::size_t bytes) noexcept;
    [[nodiscard]] std::size_t max_memory_budget_bytes() const noexcept;
    [[nodiscard]] bool budget_exhausted() const noexcept;

    // F7.0-A3.1: controlled reset for long-lived REPL/server processes.
    void reset();

private:
    template <typename T, typename... Args>
    [[nodiscard]] ExprPtr make_uncached(Args&&... args) {
        T* node = static_cast<T*>(allocate(sizeof(T), alignof(T)));
        if (node == nullptr) {
            return ExprPtr{};
        }
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

    // F7.0-A3.5: memory budget bookkeeping (under alloc_mutex_).
    std::size_t bytes_allocated_{0U};
    std::size_t max_memory_budget_bytes_{0U};  // 0 = unlimited
    bool        budget_exhausted_{false};
    // F7.0-A3.7: Hash-DoS sticky flag (under alloc_mutex_).
    bool        hash_dos_detected_{false};

    // F1.3-NEW: bump allocator state protected by its own mutex.
    mutable std::mutex alloc_mutex_;
    std::vector<Block> blocks_;
    std::vector<std::vector<ExprNode*>> node_chunks_;
    std::size_t total_nodes_{0U};

    // F1.3-NEW / HPP-016: interning table sharded dynamically by hash % num_shards_.
    std::size_t num_shards_{16U};
    mutable std::unique_ptr<std::mutex[]> intern_shards_;
    std::unique_ptr<std::unordered_set<ExprPtr, ExprHash, ExprEqual>[]> intern_shard_tables_;

    // F1.3-NEW: hot-constant caches as atomics for lock-free fast-path reads.
    std::array<std::atomic<ExprPtr>, 5U> interned_constants_{};
    std::atomic<ExprPtr> interned_zero_{ExprPtr{}};
    std::atomic<ExprPtr> interned_one_{ExprPtr{}};
    std::atomic<ExprPtr> interned_neg_one_{ExprPtr{}};
};

// ---------------------------------------------------------------------------
// visit_expr<Visitor> — compile-time dispatch on ExprKind
// ---------------------------------------------------------------------------

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

}  // namespace cas
