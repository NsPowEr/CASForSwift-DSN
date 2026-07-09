#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <functional>
#include <utility>
#include <vector>

namespace cas::symbolic {

CASContext::CASContext() : rewrite_provider_(&default_rewrite_provider()) {
    intern_shards_ = arena_.num_shards();
    ::cas::algebra::register_algebraic_simplify_hook(*this);
}

void CASContext::define(const Symbol& symbol, ExprPtr value) {
    variables_[symbol.name] = value;
}

void CASContext::clear_variables() noexcept {
    variables_.clear();
}

std::optional<ExprPtr> CASContext::lookup(const Symbol& symbol) const {
    const auto found = variables_.find(symbol.name);
    if (found == variables_.end()) {
        return std::nullopt;
    }
    return found->second;
}

Assumptions& CASContext::assumptions() noexcept {
    return assumptions_;
}

const Assumptions& CASContext::assumptions() const noexcept {
    return assumptions_;
}

AstArena& CASContext::arena() noexcept {
    return arena_;
}

const AstArena& CASContext::arena() const noexcept {
    return arena_;
}

void CASContext::set_rewrite_provider(const RewriteProvider* provider) noexcept {
    rewrite_provider_ = provider;
}

const RewriteProvider* CASContext::rewrite_provider() const noexcept {
    return rewrite_provider_;
}

void CASContext::enable_trace(bool enabled) noexcept {
    trace_enabled_ = enabled;
    if (!enabled) {
        trace_.clear();
    }
}

const ComputationTrace& CASContext::get_trace() const noexcept {
    return trace_enabled_ ? trace_ : empty_trace();
}

void CASContext::set_timeout(std::chrono::milliseconds timeout) noexcept {
    timeout_ = timeout;
    // A30 contract: an explicit wall-clock budget means the caller owns the
    // operation budget, so the default deterministic ops gate steps aside
    // (long-running factorisation/Galois workloads legitimately exceed it).
    // An explicit set_max_operation_ops call always wins over this default.
    if (!max_operation_ops_explicit_) {
        max_operation_ops_ = 0U;
    }
}

Symbol CASContext::make_fresh_symbol(const std::string& prefix) {
    std::string base = prefix.empty() ? std::string("_g") : prefix;
    while (true) {
        ++fresh_symbol_counter_;
        std::string candidate = base + "_" + std::to_string(fresh_symbol_counter_);
        if (variables_.find(candidate) == variables_.end()) {
            return Symbol(candidate);
        }
    }
}

void CASContext::clear_caches() noexcept {
    simplify_cache_.clear();
    diff_cache_.clear();
    integrate_cache_.clear();
}

void CASContext::set_caching_enabled(bool enabled) noexcept {
    caching_enabled_ = enabled;
    if (!enabled) {
        clear_caches();
    }
}

void CASContext::set_cache_limit(std::size_t limit) noexcept {
    simplify_cache_.set_max_size(limit);
    diff_cache_.set_max_size(limit);
    integrate_cache_.set_max_size(limit);
}

std::size_t CASContext::get_cache_limit() const noexcept {
    return simplify_cache_.max_size();
}

CacheMetrics CASContext::get_simplify_metrics() const noexcept {
    return simplify_cache_.metrics();
}

CacheMetrics CASContext::get_diff_metrics() const noexcept {
    return diff_cache_.metrics();
}

CacheMetrics CASContext::get_integrate_metrics() const noexcept {
    return integrate_cache_.metrics();
}

Result<ExprPtr> CASContext::simplify(ExprPtr expr) {
    if (!expr) return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot simplify null expression"));

    {
        const std::uint64_t cur_rev = assumptions_.revision();
        if (cur_rev != last_assumptions_revision_) {
            clear_caches();
            last_assumptions_revision_ = cur_rev;
        }
    }

    if (caching_enabled_ && !trace_enabled_) {
        if (auto cached = simplify_cache_.get(expr)) {
            return ok(*cached);
        }
    }

    const bool owns_operation = !operation_active_;
    if (owns_operation) {
        operation_active_ = true;
        trace_capture_active_ = trace_enabled_;
        trace_.clear();
        ops_count_ = 0;
        operation_started_at_ = std::chrono::steady_clock::now();
    }
    auto result = symbolic::simplify(expr, *this);
    if (owns_operation) {
        operation_active_ = false;
        trace_capture_active_ = false;
        ops_count_ = 0;

        if (result.is_ok() && post_simplify_hook_) {
            auto hook_res = post_simplify_hook_(result.value(), *this);
            if (hook_res.is_ok() && hook_res.value() != result.value()) {
                auto resimplified = symbolic::simplify(hook_res.value(), *this);
                result = resimplified.is_ok() ? resimplified : hook_res;
            }
        }
    }

    if (caching_enabled_ && !trace_enabled_ && result.is_ok()) {
        simplify_cache_.put(expr, result.value());
    }

#ifndef NDEBUG
    if (owns_operation && result.is_ok()
        && !is_strictly_canonical(result.value())) {
        std::fprintf(stderr,
            "[F7.0-A4.2 WARN] simplify result is not strictly canonical\n");
    }
#endif

    return result;
}

void CASContext::collect_garbage(const std::vector<ExprPtr*>& external_roots) {
    AstArena new_arena;
    std::unordered_map<ExprPtr, ExprPtr> cache;

    for (auto& [name, expr] : variables_) {
        expr = clone_into_arena(expr, new_arena, cache);
    }

    assumptions_.update_roots(new_arena, cache);

    if (trace_enabled_) {
        for (auto& step : trace_) {
            step.target_before = clone_into_arena(step.target_before, new_arena, cache);
            step.target_after = clone_into_arena(step.target_after, new_arena, cache);
            step.root_after = clone_into_arena(step.root_after, new_arena, cache);
        }
    }

    for (auto* root_ptr : external_roots) {
        if (root_ptr && *root_ptr) {
            *root_ptr = clone_into_arena(*root_ptr, new_arena, cache);
        }
    }

    if (!simplify_cache_.empty()) {
        CacheContainer<ExprPtr, ExprPtr, ExprHash, ExprEqual> new_simplify_cache(simplify_cache_.max_size());
        for (auto& it : simplify_cache_) {
            new_simplify_cache.put(clone_into_arena(it.first, new_arena, cache), clone_into_arena(it.second.first, new_arena, cache));
        }
        new_simplify_cache.metrics() = simplify_cache_.metrics();
        simplify_cache_ = std::move(new_simplify_cache);
    }

    if (!diff_cache_.empty()) {
        CacheContainer<DiffKey, ExprPtr, DiffHash> new_diff_cache(diff_cache_.max_size());
        for (auto& it : diff_cache_) {
            const auto& key = it.first;
            DiffKey new_key = {
                .expr = clone_into_arena(key.expr, new_arena, cache),
                .var_name = key.var_name,
                .order = key.order
            };
            new_diff_cache.put(new_key, clone_into_arena(it.second.first, new_arena, cache));
        }
        new_diff_cache.metrics() = diff_cache_.metrics();
        diff_cache_ = std::move(new_diff_cache);
    }

    if (!integrate_cache_.empty()) {
        CacheContainer<IntegrateKey, ExprPtr, IntegrateHash> new_integrate_cache(integrate_cache_.max_size());
        for (auto& it : integrate_cache_) {
            const auto& key = it.first;
            IntegrateKey new_key = {
                .expr = clone_into_arena(key.expr, new_arena, cache),
                .var_name = key.var_name
            };
            new_integrate_cache.put(new_key, clone_into_arena(it.second.first, new_arena, cache));
        }
        new_integrate_cache.metrics() = integrate_cache_.metrics();
        integrate_cache_ = std::move(new_integrate_cache);
    }

    arena_ = std::move(new_arena);
}

} // namespace cas::symbolic
