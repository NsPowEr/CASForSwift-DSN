#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
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

void CASContext::clear_assumptions() noexcept {
    assumptions_ = Assumptions{};
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

// A51 — invariante: CHIUNQUE apra un'operazione top-level (cioe' porti
// `operation_active_` da false a true) DEVE inizializzarne il budget qui, e
// azzerarlo alla chiusura con `end_operation_budget`.
//
// Senza l'inizializzazione, `Simplifier::check_timeout` misura `elapsed` a
// partire dall'ULTIMA operazione top-level del contesto — o dall'epoch se non
// ce n'e' mai stata una, e allora `elapsed` vale l'uptime della macchina e il
// timeout scatta al primo controllo wall-clock. L'esito di un confronto finiva
// cosi' per dipendere dalla storia del contesto e dal carico invece che dai
// suoi operandi: e' esattamente il non-determinismo misurato in A51, dove una
// entry golden oscillava SKIP<->FAIL fra misure con lo STESSO binario.
//
// Il pattern era duplicato in quattro punti (CASContext::simplify,
// CASContext::substitute, algebra::polynomial_gcd, mathematically_equal) e uno
// solo — `mathematically_equal` — lo sbagliava. Incapsularlo rende
// l'invariante impossibile da violare per omissione.
void CASContext::begin_operation_budget(bool capture_trace) noexcept {
    trace_capture_active_ = capture_trace;
    trace_.clear();
    ops_count_ = 0;
    operation_started_at_ = std::chrono::steady_clock::now();
}

void CASContext::end_operation_budget() noexcept {
    trace_capture_active_ = false;
    // A51 — high-water mark del budget effettivamente consumato. Serve a
    // tarare il gate ops su MISURA invece che a intuito: senza questo dato
    // l'unico modo di sapere se il gate deterministico morde prima del
    // wall-clock e' dedurlo, e la deduzione era sbagliata (il commento di A30
    // dava 2'000'000 ops come "coperti da 10 s"; la misura dice ~150 s).
    if (ops_count_ > ops_high_water_) ops_high_water_ = ops_count_;
    ops_count_ = 0;
}

void CASContext::reset_ops_high_water() noexcept { ops_high_water_ = 0; }

CASContext::OperationScope::OperationScope(CASContext& ctx, bool capture_trace) noexcept
    : ctx_(ctx), owns_(!ctx.operation_active_) {
    if (owns_) {
        ctx_.operation_active_ = true;
        ctx_.begin_operation_budget(capture_trace);
    }
}

// A53 — variante con tetto proprio. Il budget effettivo e' il `min` fra quello
// del contesto e `ops_cap`: un motore interno puo' essere piu' parsimonioso del
// contesto, mai piu' generoso (sarebbe un aggiramento del gate del chiamante).
// Gate spento (0) = il chiamante possiede il budget per contratto A30, e nessun
// tetto interno glielo toglie: e' cio' che rende possibile un run di sola
// misura, dove imporre il tetto falserebbe proprio il dato da raccogliere.
CASContext::OperationScope::OperationScope(CASContext& ctx, bool capture_trace,
                                           std::uint64_t ops_cap) noexcept
    : ctx_(ctx), owns_(!ctx.operation_active_) {
    if (!owns_) return;
    const std::uint64_t current = ctx_.max_operation_ops_;
    if (current != 0U && ops_cap != 0U && ops_cap < current) {
        restore_max_ops_ = current;
        restore_max_ops_explicit_ = ctx_.max_operation_ops_explicit_;
        capped_ = true;
        ctx_.max_operation_ops_ = ops_cap;
        ctx_.max_operation_ops_explicit_ = true;
    }
    ctx_.operation_active_ = true;
    ctx_.begin_operation_budget(capture_trace);
}

CASContext::OperationScope::~OperationScope() noexcept {
    if (owns_) {
        ctx_.operation_active_ = false;
        ctx_.end_operation_budget();
        if (capped_) {
            ctx_.max_operation_ops_ = restore_max_ops_;
            ctx_.max_operation_ops_explicit_ = restore_max_ops_explicit_;
        }
    }
}

// A53 — vedi symbolic.hpp per la semantica. Il mark e' una COPIA del set, non
// la sua dimensione: `SideConditionSet::add` puo' RIMUOVERE una condizione piu'
// debole quando ne entra una che la subsume (Positive subsume NonZero, §3.4),
// quindi troncare alla dimensione di partenza non ricostruirebbe lo stato.
CASContext::SideConditionRollback::SideConditionRollback(CASContext& ctx)
    : ctx_(ctx), mark_(ctx.side_conditions_) {}

CASContext::SideConditionRollback::~SideConditionRollback() noexcept {
    if (!committed_) {
        ctx_.side_conditions_ = mark_;
    }
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
            // A31 fase 1 (spec §4.2): a cache hit must re-emit the
            // conditions that were recorded when this entry was first
            // computed, or a second simplify(same expr) would silently
            // under-report last_side_conditions() vs the first call.
            // A31 fase 2 (spec §10.1): a TOP-LEVEL hit must also reset the
            // accumulator first — without this, a hit following an unrelated
            // top-level call merged into that call's leftover set and
            // over-reported (a/a, b/b, a/a-hit -> {a,b} instead of {a}).
            // An inner hit (operation_active_) keeps merging without reset.
            if (!operation_active_) side_conditions_.clear();
            side_conditions_.merge(cached->conditions);
            return ok(cached->result);
        }
    }

    const bool owns_operation = !operation_active_;
    if (owns_operation) {
        operation_active_ = true;
        begin_operation_budget(trace_enabled_);
        side_conditions_.clear();
    }
    // A31 fase 1: snapshot before processing `expr` so the PUT below can
    // store exactly the conditions this specific call contributed, not the
    // full accumulated set of a multi-expression outer operation.
    const SideConditionSet conditions_mark = side_conditions_;
    auto result = symbolic::simplify(expr, *this);
    if (owns_operation) {
        operation_active_ = false;
        end_operation_budget();

        if (result.is_ok() && post_simplify_hook_) {
            auto hook_res = post_simplify_hook_(result.value(), *this);
            if (hook_res.is_ok() && hook_res.value() != result.value()) {
                auto resimplified = symbolic::simplify(hook_res.value(), *this);
                result = resimplified.is_ok() ? resimplified : hook_res;
            }
        }
    }

    if (caching_enabled_ && !trace_enabled_ && result.is_ok()) {
        simplify_cache_.put(expr, SimplifyCacheEntry{
            result.value(), side_conditions_.since(conditions_mark)});
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
        CacheContainer<ExprPtr, SimplifyCacheEntry, ExprHash, ExprEqual> new_simplify_cache(simplify_cache_.max_size());
        for (auto& it : simplify_cache_) {
            const SimplifyCacheEntry& entry = it.second.first;
            // A31 fase 1: each condition's `subject` is an ExprPtr rooted in
            // the OLD arena — must be re-interned exactly like every other
            // ExprPtr root, or last_side_conditions()/cache hits after a GC
            // would dereference stale memory.
            SideConditionSet new_conditions;
            new_conditions.set_max_size(entry.conditions.max_size());
            for (const auto& c : entry.conditions.items()) {
                new_conditions.add(DomainCondition{
                    c.kind, clone_into_arena(c.subject, new_arena, cache)});
            }
            new_simplify_cache.put(
                clone_into_arena(it.first, new_arena, cache),
                SimplifyCacheEntry{
                    clone_into_arena(entry.result, new_arena, cache),
                    std::move(new_conditions)});
        }
        new_simplify_cache.metrics() = simplify_cache_.metrics();
        simplify_cache_ = std::move(new_simplify_cache);
    }

    {
        // A31 fase 1: re-intern the live side_conditions_ accumulator too —
        // it holds ExprPtr subjects from the OLD arena until the next
        // top-level simplify() call clears it.
        SideConditionSet new_side_conditions;
        new_side_conditions.set_max_size(side_conditions_.max_size());
        for (const auto& c : side_conditions_.items()) {
            new_side_conditions.add(DomainCondition{
                c.kind, clone_into_arena(c.subject, new_arena, cache)});
        }
        side_conditions_ = std::move(new_side_conditions);
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
