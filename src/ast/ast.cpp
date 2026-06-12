/// @file ast.cpp
/// @brief AstArena lifecycle: constructor, destructor, allocator, reset.
///
/// F8.0 / Task 1.2 — Monolith extraction:
///   - structural_equal, expr_hash, expr_kind, expr_kind_name → ast_compare.cpp
///   - clone_into_arena                                        → ast_clone.cpp
///
/// This translation unit now focuses exclusively on AstArena memory management.
#include "cas/ast.hpp"

#include <cstddef>
#include <new>

namespace cas {

namespace {

void destroy_node(ExprNode* node) noexcept {
    switch (node->kind) {
    case ExprKind::IntegerLit:
        static_cast<IntegerLit*>(node)->~IntegerLit();
        break;
    case ExprKind::RationalLit:
        static_cast<RationalLit*>(node)->~RationalLit();
        break;
    case ExprKind::ComplexLit:
        static_cast<ComplexLit*>(node)->~ComplexLit();
        break;
    case ExprKind::DecimalLit:
        static_cast<DecimalLit*>(node)->~DecimalLit();
        break;
    case ExprKind::Symbol:
        static_cast<Symbol*>(node)->~Symbol();
        break;
    case ExprKind::Constant:
        static_cast<Constant*>(node)->~Constant();
        break;
    case ExprKind::Unary:
        static_cast<Unary*>(node)->~Unary();
        break;
    case ExprKind::Binary:
        static_cast<Binary*>(node)->~Binary();
        break;
    case ExprKind::FuncCall:
        static_cast<FuncCall*>(node)->~FuncCall();
        break;
    case ExprKind::Sum:
        static_cast<Sum*>(node)->~Sum();
        break;
    case ExprKind::Product:
        static_cast<Product*>(node)->~Product();
        break;
    case ExprKind::Integral:
        static_cast<Integral*>(node)->~Integral();
        break;
    case ExprKind::Derivative:
        static_cast<Derivative*>(node)->~Derivative();
        break;
    case ExprKind::Limit:
        static_cast<Limit*>(node)->~Limit();
        break;
    case ExprKind::RootOf:
        static_cast<RootOf*>(node)->~RootOf();
        break;
    case ExprKind::Matrix:
        static_cast<Matrix*>(node)->~Matrix();
        break;
    case ExprKind::SeriesExp:
        static_cast<SeriesExp*>(node)->~SeriesExp();
        break;
    case ExprKind::Quantity:
        static_cast<Quantity*>(node)->~Quantity();
        break;
    case ExprKind::Null:
        break;
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// AstArena move constructor / move-assign
// ---------------------------------------------------------------------------

// F1.3-NEW: move ctor — atomics and mutex arrays are not movable, so we
// transfer the data members that can be moved under a full lock of the source.
AstArena::AstArena(AstArena&& other) noexcept {
    std::lock_guard<std::mutex> alloc_lock(other.alloc_mutex_);
    for (auto& shard : other.intern_shards_) shard.lock();

    blocks_ = std::move(other.blocks_);
    node_chunks_ = std::move(other.node_chunks_);
    total_nodes_ = other.total_nodes_;
    other.total_nodes_ = 0;
    for (std::size_t i = 0; i < N_INTERN_SHARDS; ++i) {
        intern_shard_tables_[i] = std::move(other.intern_shard_tables_[i]);
    }

    for (std::size_t i = 0; i < interned_constants_.size(); ++i) {
        interned_constants_[i].store(
            other.interned_constants_[i].load(std::memory_order_relaxed),
            std::memory_order_relaxed);
        other.interned_constants_[i].store(ExprPtr{}, std::memory_order_relaxed);
    }
    interned_zero_.store(other.interned_zero_.load(std::memory_order_relaxed), std::memory_order_relaxed);
    interned_one_.store(other.interned_one_.load(std::memory_order_relaxed), std::memory_order_relaxed);
    interned_neg_one_.store(other.interned_neg_one_.load(std::memory_order_relaxed), std::memory_order_relaxed);
    other.interned_zero_.store(ExprPtr{}, std::memory_order_relaxed);
    other.interned_one_.store(ExprPtr{}, std::memory_order_relaxed);
    other.interned_neg_one_.store(ExprPtr{}, std::memory_order_relaxed);

    for (auto& shard : other.intern_shards_) shard.unlock();
}

AstArena& AstArena::operator=(AstArena&& other) noexcept {
    if (this != &other) {
        for (auto& chunk : node_chunks_)
            for (ExprNode* node : chunk)
                if (node) destroy_node(node);

        std::lock_guard<std::mutex> alloc_lock_other(other.alloc_mutex_);
        for (auto& shard : other.intern_shards_) shard.lock();
        std::lock_guard<std::mutex> alloc_lock_self(alloc_mutex_);

        blocks_ = std::move(other.blocks_);
        node_chunks_ = std::move(other.node_chunks_);
        total_nodes_ = other.total_nodes_;
        other.total_nodes_ = 0;
        for (std::size_t i = 0; i < N_INTERN_SHARDS; ++i) {
            intern_shard_tables_[i] = std::move(other.intern_shard_tables_[i]);
        }

        for (std::size_t i = 0; i < interned_constants_.size(); ++i) {
            interned_constants_[i].store(
                other.interned_constants_[i].load(std::memory_order_relaxed),
                std::memory_order_relaxed);
            other.interned_constants_[i].store(ExprPtr{}, std::memory_order_relaxed);
        }
        interned_zero_.store(other.interned_zero_.load(std::memory_order_relaxed), std::memory_order_relaxed);
        interned_one_.store(other.interned_one_.load(std::memory_order_relaxed), std::memory_order_relaxed);
        interned_neg_one_.store(other.interned_neg_one_.load(std::memory_order_relaxed), std::memory_order_relaxed);
        other.interned_zero_.store(ExprPtr{}, std::memory_order_relaxed);
        other.interned_one_.store(ExprPtr{}, std::memory_order_relaxed);
        other.interned_neg_one_.store(ExprPtr{}, std::memory_order_relaxed);

        for (auto& shard : other.intern_shards_) shard.unlock();
    }
    return *this;
}

// ---------------------------------------------------------------------------
// AstArena destructor / memory management
// ---------------------------------------------------------------------------

AstArena::~AstArena() {
    for (auto& chunk : node_chunks_) {
        for (ExprNode* node : chunk) {
            if (node != nullptr) {
                destroy_node(node);
            }
        }
    }
}

void* AstArena::allocate(std::size_t size, std::size_t alignment) {
    if (budget_exhausted_) return nullptr;  // F7.0-A3.5
    if (blocks_.empty()) {
        append_block(size + alignment);
        if (budget_exhausted_) return nullptr;
    }

    auto align_up = [](std::size_t value, std::size_t align) noexcept {
        const std::size_t remainder = value % align;
        return remainder == 0U ? value : value + (align - remainder);
    };

    Block* block = &blocks_.back();
    std::size_t offset = align_up(block->used, alignment);
    if (offset + size > block->capacity) {
        append_block(size + alignment);
        if (budget_exhausted_) return nullptr;
        block = &blocks_.back();
        offset = align_up(block->used, alignment);
    }

    std::byte* memory = block->data.get() + offset;
    block->used = offset + size;
    return memory;
}

void AstArena::append_block(std::size_t minimum_bytes) {
    const std::size_t capacity = minimum_bytes > DEFAULT_BLOCK_BYTES ? minimum_bytes : DEFAULT_BLOCK_BYTES;
    if (max_memory_budget_bytes_ > 0
        && bytes_allocated_ + capacity > max_memory_budget_bytes_) {
        budget_exhausted_ = true;
        return;
    }
    blocks_.push_back(Block{
        .data = std::make_unique<std::byte[]>(capacity),
        .capacity = capacity,
        .used = 0U,
    });
    bytes_allocated_ += capacity;
}

std::size_t AstArena::bytes_allocated() const noexcept {
    std::lock_guard<std::mutex> lock(alloc_mutex_);
    return bytes_allocated_;
}

void AstArena::set_max_memory_budget_bytes(std::size_t bytes) noexcept {
    std::lock_guard<std::mutex> lock(alloc_mutex_);
    max_memory_budget_bytes_ = bytes;
    if (bytes == 0U || bytes_allocated_ <= bytes) {
        budget_exhausted_ = false;
    }
}

std::size_t AstArena::max_memory_budget_bytes() const noexcept {
    std::lock_guard<std::mutex> lock(alloc_mutex_);
    return max_memory_budget_bytes_;
}

bool AstArena::budget_exhausted() const noexcept {
    std::lock_guard<std::mutex> lock(alloc_mutex_);
    return budget_exhausted_;
}

bool AstArena::hash_dos_detected() const noexcept {
    std::lock_guard<std::mutex> lock(alloc_mutex_);
    return hash_dos_detected_;
}

void AstArena::clear_hash_dos_flag() noexcept {
    std::lock_guard<std::mutex> lock(alloc_mutex_);
    hash_dos_detected_ = false;
}

std::size_t AstArena::size() const noexcept {
    std::lock_guard<std::mutex> lock(alloc_mutex_);
    return total_nodes_;
}

void AstArena::reset() {
    // F7.0-A3.1: lock everything in canonical order (shard → alloc).
    for (auto& shard : intern_shards_) shard.lock();
    std::lock_guard<std::mutex> alloc_lock(alloc_mutex_);

    for (auto& chunk : node_chunks_) {
        for (ExprNode* node : chunk) {
            if (node != nullptr) destroy_node(node);
        }
    }
    node_chunks_.clear();
    blocks_.clear();

    for (auto& table : intern_shard_tables_) table.clear();

    for (auto& cell : interned_constants_) {
        cell.store(ExprPtr{}, std::memory_order_release);
    }
    interned_zero_.store(ExprPtr{}, std::memory_order_release);
    interned_one_.store(ExprPtr{}, std::memory_order_release);
    interned_neg_one_.store(ExprPtr{}, std::memory_order_release);

    total_nodes_ = 0U;
    bytes_allocated_ = 0U;
    budget_exhausted_ = false;
    hash_dos_detected_ = false;

    for (auto& shard : intern_shards_) shard.unlock();
}

}  // namespace cas
