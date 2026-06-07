#include "cas/ast.hpp"

#include <cstddef>
#include <new>

namespace cas {

namespace {

bool optional_expr_equal(const std::optional<ExprPtr>& lhs, const std::optional<ExprPtr>& rhs) noexcept {
    if (lhs.has_value() != rhs.has_value()) {
        return false;
    }

    if (!lhs.has_value()) {
        return true;
    }

    return structural_equal(*lhs, *rhs);
}

template <typename Container>
bool expr_ptr_range_equal(const Container& lhs, const Container& rhs) noexcept {
    if (lhs.size() != rhs.size()) {
        return false;
    }

    for (std::size_t index = 0; index < lhs.size(); ++index) {
        if (!structural_equal(lhs[index], rhs[index])) {
            return false;
        }
    }

    return true;
}

bool big_int_equal(const BigInt& lhs, const BigInt& rhs) noexcept {
    return lhs == rhs;
}

void destroy_node(ExprNode* node) noexcept {
    switch (node->kind) {
    case ExprKind::IntegerLit:
        static_cast<IntegerLit*>(node)->~IntegerLit();
        break;
    case ExprKind::RationalLit:
    case ExprKind::ComplexLit:
        static_cast<ComplexLit*>(node)->~ComplexLit();
        break;

        static_cast<RationalLit*>(node)->~RationalLit();
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

// F1.3-NEW: move ctor — atomics and mutex arrays are not movable, so we
// transfer the data members that can be moved under a full lock of the source.
AstArena::AstArena(AstArena&& other) noexcept {
    // Lock all shards + alloc mutex of `other` to ensure no concurrent access.
    std::lock_guard<std::mutex> alloc_lock(other.alloc_mutex_);
    for (auto& shard : other.intern_shards_) shard.lock();

    blocks_ = std::move(other.blocks_);
    node_chunks_ = std::move(other.node_chunks_);
    total_nodes_ = other.total_nodes_;
    other.total_nodes_ = 0;
    // Move all per-shard interning tables (replaces single interning_table_).
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
        // Destroy current nodes first.
        for (auto& chunk : node_chunks_)
            for (ExprNode* node : chunk)
                if (node) destroy_node(node);

        // Lock both arenas to prevent races.
        std::lock_guard<std::mutex> alloc_lock_other(other.alloc_mutex_);
        for (auto& shard : other.intern_shards_) shard.lock();
        std::lock_guard<std::mutex> alloc_lock_self(alloc_mutex_);

        blocks_ = std::move(other.blocks_);
        node_chunks_ = std::move(other.node_chunks_);
        total_nodes_ = other.total_nodes_;
        other.total_nodes_ = 0;
        // Move all per-shard interning tables (replaces single interning_table_).
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
    if (blocks_.empty()) {
        append_block(size + alignment);
    }

    auto align_up = [](std::size_t value, std::size_t align) noexcept {
        const std::size_t remainder = value % align;
        return remainder == 0U ? value : value + (align - remainder);
    };

    Block* block = &blocks_.back();
    std::size_t offset = align_up(block->used, alignment);
    if (offset + size > block->capacity) {
        append_block(size + alignment);
        block = &blocks_.back();
        offset = align_up(block->used, alignment);
    }

    std::byte* memory = block->data.get() + offset;
    block->used = offset + size;
    return memory;
}

void AstArena::append_block(std::size_t minimum_bytes) {
    const std::size_t capacity = minimum_bytes > DEFAULT_BLOCK_BYTES ? minimum_bytes : DEFAULT_BLOCK_BYTES;
    blocks_.push_back(Block{
        .data = std::make_unique<std::byte[]>(capacity),
        .capacity = capacity,
        .used = 0U,
    });
}

std::size_t AstArena::size() const noexcept {
    // F1.3-NEW: use alloc_mutex_ (total_nodes_ is modified under it).
    std::lock_guard<std::mutex> lock(alloc_mutex_);
    return total_nodes_;
}

ExprKind expr_kind(ExprPtr expr) noexcept {
    return expr ? expr->kind : ExprKind::Null;
}

bool structural_equal(ExprPtr lhs, ExprPtr rhs) noexcept {
    if (lhs == rhs) {
        return true;
    }

    if (!lhs || !rhs) {
        return false;
    }

    if (lhs->kind != rhs->kind) {
        return false;
    }

    switch (lhs->kind) {
    case ExprKind::IntegerLit: {
        const auto& lhs_value = expr_ref<IntegerLit>(lhs);
        const auto& rhs_value = expr_ref<IntegerLit>(rhs);
        return big_int_equal(lhs_value.value, rhs_value.value);
    }
    case ExprKind::RationalLit: {
        const auto& lhs_value = expr_ref<RationalLit>(lhs);
        const auto& rhs_value = expr_ref<RationalLit>(rhs);
        return big_int_equal(lhs_value.numerator, rhs_value.numerator) &&
               big_int_equal(lhs_value.denominator, rhs_value.denominator);
    }

    case ExprKind::ComplexLit: {
        const auto& l = expr_ref<ComplexLit>(lhs);
        const auto& r = expr_ref<ComplexLit>(rhs);
        return big_int_equal(l.re_num, r.re_num) && big_int_equal(l.re_den, r.re_den) &&
               big_int_equal(l.im_num, r.im_num) && big_int_equal(l.im_den, r.im_den);
    }

    case ExprKind::DecimalLit:
        return expr_ref<DecimalLit>(lhs).text == expr_ref<DecimalLit>(rhs).text;
    case ExprKind::Symbol:
        return expr_ref<Symbol>(lhs).name == expr_ref<Symbol>(rhs).name;
    case ExprKind::Constant:
        return expr_ref<Constant>(lhs).value == expr_ref<Constant>(rhs).value;
    case ExprKind::Unary: {
        const auto& lhs_value = expr_ref<Unary>(lhs);
        const auto& rhs_value = expr_ref<Unary>(rhs);
        return lhs_value.op == rhs_value.op && structural_equal(lhs_value.operand, rhs_value.operand);
    }
    case ExprKind::Binary: {
        const auto& lhs_value = expr_ref<Binary>(lhs);
        const auto& rhs_value = expr_ref<Binary>(rhs);
        return lhs_value.op == rhs_value.op &&
               structural_equal(lhs_value.left, rhs_value.left) &&
               structural_equal(lhs_value.right, rhs_value.right);
    }
    case ExprKind::FuncCall: {
        const auto& lhs_value = expr_ref<FuncCall>(lhs);
        const auto& rhs_value = expr_ref<FuncCall>(rhs);
        if (lhs_value.func_id != rhs_value.func_id) return false;
        if (lhs_value.func_id == BuiltinOp::Unknown && lhs_value.name != rhs_value.name) return false;
        return expr_ptr_range_equal(lhs_value.args, rhs_value.args);
    }
    case ExprKind::Sum:
        return expr_ptr_range_equal(expr_ref<Sum>(lhs).terms, expr_ref<Sum>(rhs).terms);
    case ExprKind::Product:
        return expr_ptr_range_equal(expr_ref<Product>(lhs).factors, expr_ref<Product>(rhs).factors);
    case ExprKind::Integral: {
        const auto& lhs_value = expr_ref<Integral>(lhs);
        const auto& rhs_value = expr_ref<Integral>(rhs);
        return structural_equal(lhs_value.integrand, rhs_value.integrand) &&
               lhs_value.variable.name == rhs_value.variable.name &&
               optional_expr_equal(lhs_value.lower, rhs_value.lower) &&
               optional_expr_equal(lhs_value.upper, rhs_value.upper);
    }
    case ExprKind::Derivative: {
        const auto& lhs_value = expr_ref<Derivative>(lhs);
        const auto& rhs_value = expr_ref<Derivative>(rhs);
        return structural_equal(lhs_value.expression, rhs_value.expression) &&
               lhs_value.variable.name == rhs_value.variable.name &&
               lhs_value.order == rhs_value.order;
    }
    case ExprKind::Limit: {
        const auto& lhs_value = expr_ref<Limit>(lhs);
        const auto& rhs_value = expr_ref<Limit>(rhs);
        return structural_equal(lhs_value.expression, rhs_value.expression) &&
               lhs_value.variable.name == rhs_value.variable.name &&
               structural_equal(lhs_value.point, rhs_value.point) &&
               lhs_value.direction == rhs_value.direction;
    }
    case ExprKind::RootOf: {
        const auto& lhs_value = expr_ref<RootOf>(lhs);
        const auto& rhs_value = expr_ref<RootOf>(rhs);
        return structural_equal(lhs_value.polynomial, rhs_value.polynomial) &&
               lhs_value.variable.name == rhs_value.variable.name &&
               lhs_value.root_index == rhs_value.root_index;
    }
    case ExprKind::Matrix: {
        const auto& lhs_value = expr_ref<Matrix>(lhs);
        const auto& rhs_value = expr_ref<Matrix>(rhs);
        return lhs_value.rows == rhs_value.rows &&
               lhs_value.cols == rhs_value.cols &&
               expr_ptr_range_equal(lhs_value.elements, rhs_value.elements);
    }
    case ExprKind::SeriesExp: {
        const auto& l = expr_ref<SeriesExp>(lhs);
        const auto& r = expr_ref<SeriesExp>(rhs);
        if (l.var.name != r.var.name || l.order != r.order) return false;
        if (!structural_equal(l.point, r.point)) return false;
        if (l.terms.size() != r.terms.size()) return false;
        for (std::size_t i = 0; i < l.terms.size(); ++i) {
            if (l.terms[i].first != r.terms[i].first) return false;
            if (!structural_equal(l.terms[i].second, r.terms[i].second)) return false;
        }
        return true;
    }
    case ExprKind::Quantity: {
        const auto& l = expr_ref<Quantity>(lhs);
        const auto& r = expr_ref<Quantity>(rhs);
        return l.dimensions == r.dimensions && structural_equal(l.value, r.value);
    }
    case ExprKind::Null:
        return false;
    }

    return false;
}

std::size_t expr_hash(ExprPtr expr) noexcept {
    if (!expr) return 0;
    
    std::size_t cached = expr->hash_cache_.load(std::memory_order_relaxed);
    if (cached != 0) return cached;

    auto hash_combine = [](std::size_t& seed, std::size_t value) {
        seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };

    std::size_t seed = static_cast<std::size_t>(expr->kind);

    switch (expr->kind) {
    case ExprKind::IntegerLit:
    case ExprKind::ComplexLit: {
        const auto& node = expr_ref<ComplexLit>(expr);
        hash_combine(seed, node.re_num.hash());
        hash_combine(seed, node.re_den.hash());
        hash_combine(seed, node.im_num.hash());
        hash_combine(seed, node.im_den.hash());
        break;
    }

        hash_combine(seed, expr_ref<IntegerLit>(expr).value.hash());
        break;
    case ExprKind::RationalLit: {
        const auto& node = expr_ref<RationalLit>(expr);
        hash_combine(seed, node.numerator.hash());
        hash_combine(seed, node.denominator.hash());
        break;
    }
    case ExprKind::DecimalLit:
        hash_combine(seed, std::hash<std::string>{}(expr_ref<DecimalLit>(expr).text));
        break;
    case ExprKind::Symbol:
        hash_combine(seed, std::hash<std::string>{}(expr_ref<Symbol>(expr).name));
        break;
    case ExprKind::Constant:
        hash_combine(seed, static_cast<std::size_t>(expr_ref<Constant>(expr).value));
        break;
    case ExprKind::Unary: {
        const auto& node = expr_ref<Unary>(expr);
        hash_combine(seed, static_cast<std::size_t>(node.op));
        hash_combine(seed, expr_hash(node.operand));
        break;
    }
    case ExprKind::Binary: {
        const auto& node = expr_ref<Binary>(expr);
        hash_combine(seed, static_cast<std::size_t>(node.op));
        hash_combine(seed, expr_hash(node.left));
        hash_combine(seed, expr_hash(node.right));
        break;
    }
    case ExprKind::FuncCall: {
        const auto& node = expr_ref<FuncCall>(expr);
        hash_combine(seed, static_cast<std::size_t>(node.func_id));
        if (node.func_id == BuiltinOp::Unknown) {
            hash_combine(seed, std::hash<std::string>{}(node.name));
        }
        for (auto arg : node.args) hash_combine(seed, expr_hash(arg));
        break;
    }
    case ExprKind::Sum:
        for (auto term : expr_ref<Sum>(expr).terms) hash_combine(seed, expr_hash(term));
        break;
    case ExprKind::Product:
        for (auto factor : expr_ref<Product>(expr).factors) hash_combine(seed, expr_hash(factor));
        break;
    case ExprKind::Integral: {
        const auto& node = expr_ref<Integral>(expr);
        hash_combine(seed, expr_hash(node.integrand));
        hash_combine(seed, std::hash<std::string>{}(node.variable.name));
        if (node.lower) hash_combine(seed, expr_hash(*node.lower));
        if (node.upper) hash_combine(seed, expr_hash(*node.upper));
        break;
    }
    case ExprKind::Derivative: {
        const auto& node = expr_ref<Derivative>(expr);
        hash_combine(seed, expr_hash(node.expression));
        hash_combine(seed, std::hash<std::string>{}(node.variable.name));
        hash_combine(seed, node.order);
        break;
    }
    case ExprKind::Limit: {
        const auto& node = expr_ref<Limit>(expr);
        hash_combine(seed, expr_hash(node.expression));
        hash_combine(seed, std::hash<std::string>{}(node.variable.name));
        hash_combine(seed, expr_hash(node.point));
        hash_combine(seed, static_cast<std::size_t>(node.direction));
        break;
    }
    case ExprKind::RootOf: {
        const auto& node = expr_ref<RootOf>(expr);
        hash_combine(seed, expr_hash(node.polynomial));
        hash_combine(seed, std::hash<std::string>{}(node.variable.name));
        if (node.root_index) hash_combine(seed, *node.root_index);
        break;
    }
    case ExprKind::Matrix: {
        const auto& node = expr_ref<Matrix>(expr);
        hash_combine(seed, node.rows);
        hash_combine(seed, node.cols);
        for (auto elem : node.elements) hash_combine(seed, expr_hash(elem));
        break;
    }
    case ExprKind::SeriesExp: {
        const auto& node = expr_ref<SeriesExp>(expr);
        hash_combine(seed, std::hash<std::string>{}(node.var.name));
        hash_combine(seed, expr_hash(node.point));
        hash_combine(seed, static_cast<std::size_t>(node.order));
        for (const auto& [exp, coeff] : node.terms) {
            hash_combine(seed, static_cast<std::size_t>(exp));
            hash_combine(seed, expr_hash(coeff));
        }
        break;
    }
    case ExprKind::Quantity: {
        const auto& node = expr_ref<Quantity>(expr);
        hash_combine(seed, expr_hash(node.value));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.m));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.kg));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.s));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.A));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.K));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.mol));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.cd));
        break;
    }
    case ExprKind::Null:
        break;
    }

    std::size_t final_hash = (seed == 0) ? 1 : seed;
    expr->hash_cache_.store(final_hash, std::memory_order_relaxed);
    return final_hash;
}

std::string_view expr_kind_name(ExprKind kind) noexcept {
    switch (kind) {
    case ExprKind::Null:
        return "Null";
    case ExprKind::IntegerLit:
        return "IntegerLit";
    case ExprKind::RationalLit:
        return "RationalLit";
    case ExprKind::ComplexLit:
        return "ComplexLit";

    case ExprKind::DecimalLit:
        return "DecimalLit";
    case ExprKind::Symbol:
        return "Symbol";
    case ExprKind::Constant:
        return "Constant";
    case ExprKind::Unary:
        return "Unary";
    case ExprKind::Binary:
        return "Binary";
    case ExprKind::FuncCall:
        return "FuncCall";
    case ExprKind::Sum:
        return "Sum";
    case ExprKind::Product:
        return "Product";
    case ExprKind::Integral:
        return "Integral";
    case ExprKind::Derivative:
        return "Derivative";
    case ExprKind::Limit:
        return "Limit";
    case ExprKind::RootOf:
        return "RootOf";
    case ExprKind::Matrix:
        return "Matrix";
    case ExprKind::SeriesExp:
        return "SeriesExp";
    case ExprKind::Quantity:
        return "Quantity";
    }

    return "Unknown";
}

ExprPtr clone_into_arena(ExprPtr expr, AstArena& target, std::unordered_map<ExprPtr, ExprPtr>& cache) {
    if (!expr) {
        return expr;
    }

    if (auto it = cache.find(expr); it != cache.end()) {
        return it->second;
    }

    ExprPtr cloned;
    switch (expr->kind) {
    case ExprKind::IntegerLit:
        cloned = target.make<IntegerLit>(expr_ref<IntegerLit>(expr).value);
        break;
    case ExprKind::RationalLit: {
        const auto& node = expr_ref<RationalLit>(expr);
        cloned = target.make<RationalLit>(node.numerator, node.denominator);
        break;
    }
    case ExprKind::ComplexLit: {
        const auto& node = expr_ref<ComplexLit>(expr);
        cloned = target.make<ComplexLit>(node.re_num, node.re_den, node.im_num, node.im_den);
        break;
    }

    case ExprKind::DecimalLit:
        cloned = target.make<DecimalLit>(expr_ref<DecimalLit>(expr).text);
        break;
    case ExprKind::Symbol:
        cloned = target.make<Symbol>(expr_ref<Symbol>(expr).name);
        break;
    case ExprKind::Constant:
        cloned = target.make<Constant>(expr_ref<Constant>(expr).value);
        break;
    case ExprKind::Unary: {
        const auto& node = expr_ref<Unary>(expr);
        cloned = target.make<Unary>(node.op, clone_into_arena(node.operand, target, cache));
        break;
    }
    case ExprKind::Binary: {
        const auto& node = expr_ref<Binary>(expr);
        auto left = clone_into_arena(node.left, target, cache);
        auto right = clone_into_arena(node.right, target, cache);
        cloned = target.make<Binary>(node.op, left, right);
        break;
    }
    case ExprKind::FuncCall: {
        const auto& node = expr_ref<FuncCall>(expr);
        std::vector<ExprPtr> args;
        args.reserve(node.args.size());
        for (auto arg : node.args) {
            args.push_back(clone_into_arena(arg, target, cache));
        }
        cloned = target.make<FuncCall>(node.name, std::move(args));
        break;
    }
    case ExprKind::Sum: {
        const auto& node = expr_ref<Sum>(expr);
        std::vector<ExprPtr> terms;
        terms.reserve(node.terms.size());
        for (auto term : node.terms) {
            terms.push_back(clone_into_arena(term, target, cache));
        }
        cloned = target.make<Sum>(std::move(terms));
        break;
    }
    case ExprKind::Product: {
        const auto& node = expr_ref<Product>(expr);
        std::vector<ExprPtr> factors;
        factors.reserve(node.factors.size());
        for (auto factor : node.factors) {
            factors.push_back(clone_into_arena(factor, target, cache));
        }
        cloned = target.make<Product>(std::move(factors));
        break;
    }
    case ExprKind::Integral: {
        const auto& node = expr_ref<Integral>(expr);
        auto integrand = clone_into_arena(node.integrand, target, cache);
        auto lower = node.lower ? std::optional<ExprPtr>(clone_into_arena(*node.lower, target, cache)) : std::nullopt;
        auto upper = node.upper ? std::optional<ExprPtr>(clone_into_arena(*node.upper, target, cache)) : std::nullopt;
        // Node.variable is a Symbol (embedded, not ExprPtr), but it's an ExprNode in fact?
        // Wait, Symbol is a struct Symbol : ExprNode. Integral::variable is Symbol (not ExprPtr).
        // That's slightly inconsistent with how we handle other nodes.
        // Let's re-check Integral struct.
        cloned = target.make<Integral>(integrand, Symbol(node.variable.name), lower, upper);
        break;
    }
    case ExprKind::Derivative: {
        const auto& node = expr_ref<Derivative>(expr);
        auto expression = clone_into_arena(node.expression, target, cache);
        cloned = target.make<Derivative>(expression, Symbol(node.variable.name), node.order);
        break;
    }
    case ExprKind::Limit: {
        const auto& node = expr_ref<Limit>(expr);
        auto expression = clone_into_arena(node.expression, target, cache);
        auto point = clone_into_arena(node.point, target, cache);
        cloned = target.make<Limit>(expression, Symbol(node.variable.name), point, node.direction);
        break;
    }
    case ExprKind::RootOf: {
        const auto& node = expr_ref<RootOf>(expr);
        auto polynomial = clone_into_arena(node.polynomial, target, cache);
        cloned = target.make<RootOf>(polynomial, Symbol(node.variable.name), node.root_index);
        break;
    }
    case ExprKind::Matrix: {
        const auto& node = expr_ref<Matrix>(expr);
        std::vector<ExprPtr> cloned_elements;
        cloned_elements.reserve(node.elements.size());
        for (auto elem : node.elements) {
            cloned_elements.push_back(clone_into_arena(elem, target, cache));
        }
        cloned = target.make<Matrix>(node.rows, node.cols, std::move(cloned_elements));
        break;
    }
    case ExprKind::SeriesExp: {
        const auto& node = expr_ref<SeriesExp>(expr);
        std::vector<std::pair<long long, ExprPtr>> cloned_terms;
        cloned_terms.reserve(node.terms.size());
        for (const auto& [exp, coeff] : node.terms) {
            cloned_terms.push_back({exp, clone_into_arena(coeff, target, cache)});
        }
        cloned = target.make<SeriesExp>(Symbol(node.var.name), clone_into_arena(node.point, target, cache), std::move(cloned_terms), node.order);
        break;
    }
    case ExprKind::Quantity: {
        const auto& node = expr_ref<Quantity>(expr);
        cloned = target.make<Quantity>(clone_into_arena(node.value, target, cache), node.dimensions);
        break;
    }
    case ExprKind::Null:
        std::abort();
    }

    cache[expr] = cloned;
    return cloned;
}

}  // namespace cas
