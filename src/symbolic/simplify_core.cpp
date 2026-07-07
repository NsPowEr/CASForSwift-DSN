#include "simplify_impl.hpp"
#include "cas/error_helpers.hpp"


namespace cas::symbolic {
namespace detail {

Simplifier::Simplifier(
    AstArena& arena,
    const Assumptions* assumptions,
    const RewriteProvider* rewrite_provider,
    ComputationTrace* trace,
    bool trace_enabled,
    std::chrono::milliseconds timeout,
    std::chrono::steady_clock::time_point* operation_started_at,
    std::uint64_t* ops_count,
    CASContext* context)
    : arena_(arena),
      assumptions_(assumptions),
      rewrite_provider_(rewrite_provider),
      trace_(trace),
      trace_enabled_(trace_enabled),
      timeout_(timeout),
      operation_started_at_(operation_started_at),
      ops_count_(ops_count),
      context_(context) {}

Result<ExprPtr> Simplifier::simplify_expr(ExprPtr expr) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot simplify null expression"));
    }

    // L0-12 contract: the simplifier recursion budget stays governed by
    // max_simplification_depth (set_max_simplification_depth); the A20
    // max_recursion_depth parameter governs the numeric evaluators.
    const int max_depth = (context_ != nullptr) ? context_->max_simplification_depth() : MAX_SIMPLIFICATION_DEPTH;
    DepthGuard guard(max_depth);
    if (guard.exceeded()) {
        return make_unimplemented<ExprPtr>(
            "symbolic", "Simplifier::simplify_expr",
            "expression recursion depth",
            error::reason_codes::RECURSION_DEPTH_EXCEEDED,
            "Increase max_simplification_depth in CASContext",
            "A20",
            "Simplification recursion depth limit exceeded");
    }

    auto timeout = check_timeout();
    if (timeout.is_error()) {
        return fail<ExprPtr>(timeout.error());
    }

    CycleGuard cycle_guard(expr);
    if (cycle_guard.cycle_detected()) {
        return make_unimplemented<ExprPtr>(
            "symbolic", "Simplifier::simplify_expr",
            "expression cycle",
            error::reason_codes::CYCLE_DETECTED,
            "Ensure expression AST has no cyclic references",
            "A20",
            "Cyclic evaluation detected during simplification");
    }

    auto simplified = visit_expr(
        expr,
        [this, expr](const auto& node) -> Result<ExprPtr> {
            return simplify_node(expr, node);
        });
    
    return simplified;
}

Simplifier::ScopedFrame::ScopedFrame(Simplifier& owner, std::function<ExprPtr(ExprPtr)> builder)
    : owner_(owner), active_(owner.trace_enabled_) {
    if (active_) {
        owner_.root_frames_.push_back(std::move(builder));
    }
}

Simplifier::ScopedFrame::~ScopedFrame() {
    if (active_) {
        owner_.root_frames_.pop_back();
    }
}

Result<void> Simplifier::check_timeout() {
    // F7.0-A3.3: poll cancellation flag FIRST, before allocation-based timeout.
    // The interrupted_ atomic can be set asynchronously by an external
    // controller (server, REPL signal handler, watchdog thread). Reported as
    // Timeout with explicit message to avoid breaking ABI with a new
    // CASErrorKind enumerator.
    if (context_ != nullptr && context_->is_interrupted()) {
        return Result<void>(make_error(
            CASErrorKind::Timeout,
            "Operation cancelled by interrupt request"));
    }
    if (operation_started_at_ == nullptr || ops_count_ == nullptr) {
        return ok();
    }
    ++(*ops_count_);
    const std::uint64_t interval = (context_ != nullptr) ? context_->timeout_check_interval() : 1024U;
    if ((*ops_count_ % interval) != 0U) {
        return ok();
    }
    auto elapsed = std::chrono::steady_clock::now() - *operation_started_at_;
    if (elapsed >= timeout_) {
        return Result<void>(make_error(CASErrorKind::Timeout, "Symbolic operation timed out"));
    }
    return ok();
}

ExprPtr Simplifier::build_root_after(ExprPtr target_after) {
    ExprPtr root = target_after;
    for (auto it = root_frames_.rbegin(); it != root_frames_.rend(); ++it) {
        root = (*it)(root);
    }
    return root;
}

void Simplifier::append_trace(RuleId rule_id, ExprPtr before, ExprPtr after, bool allow_identity) {
    if (!trace_enabled_ || trace_ == nullptr) {
        return;
    }
    if (!allow_identity && before == after) {
        return;
    }

    trace_->push_back(TraceStep{
        .rule_id = rule_id,
        .depth = static_cast<std::uint8_t>(std::min<std::size_t>(root_frames_.size(), 255U)),
        .target_before = before,
        .target_after = after,
        .root_after = build_root_after(after),
    });
}

Result<ExprPtr> Simplifier::traced_result(RuleId rule_id, ExprPtr before, ExprPtr after) {
    append_trace(rule_id, before, after);
    return ok(after);
}

void Simplifier::append_assumption(ExprPtr target) {
    append_trace(RuleId::AssumptionApplied, target, target, true);
}

ExprPtr Simplifier::make_sum_target(const std::vector<ExprPtr>& terms) {
    if (terms.size() == 1U) {
        return terms.front();
    }
    return arena_.make<Sum>(terms);
}

ExprPtr Simplifier::make_product_target(const std::vector<ExprPtr>& factors) {
    if (factors.size() == 1U) {
        return factors.front();
    }
    return arena_.make<Product>(factors);
}

} // namespace detail

Result<ExprPtr> simplify(ExprPtr expr, AstArena& arena) {
    return detail::Simplifier(arena, nullptr, &default_rewrite_provider()).simplify_expr(expr);
}

Result<ExprPtr> simplify(ExprPtr expr, CASContext& context) {
    return detail::Simplifier(
        context.arena(),
        &context.assumptions(),
        context.rewrite_provider(),
        &context.trace_,
        context.trace_capture_active_,
        context.timeout_,
        &context.operation_started_at_,
        &context.ops_count_,
        &context).simplify_expr(expr);
}

} // namespace cas::symbolic
