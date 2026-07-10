// limit_branch_cut.cpp — HC-F8-PENDING-20 / Branch_Cut_Propagation.md §3.2:
// direction-limit table at branch-cut edges.
//
// For a one-sided limit of f(w(t)) as t → t₀ where w(t₀) lands on the branch
// cut of f (negative real axis for sqrt and ln) and Im(w(t)) genuinely
// depends on t, the value of the limit is the CUT-EDGE value determined by
// the side from which Im(w) approaches 0:
//
//   f      Im → 0⁺ (top edge)        Im → 0⁻ (bottom edge)
//   sqrt   +i·sqrt(−w₀)              −i·sqrt(−w₀)
//   ln     ln(−w₀) + iπ              ln(−w₀) − iπ
//
// The top edge coincides with the principal branch (arg = π included), which
// is what plain substitution produces; before this table the direction was
// silently ignored and the bottom edge got the top-edge value.
//
// The approach side is the sign of the first non-vanishing derivative of
// Im(w) at t₀ (valuation scan, same idiom as the signed-pole scan in
// limit_rules.cpp): for t → t₀⁺ the sign of Im is sign(c_k); for t → t₀⁻ it
// flips iff k is odd.
//
// Soundness scope (REGOLA ZERO):
//   - The table applies only when the expression IS the cut node f(w):
//     replacing a subterm by its limit value inside an arbitrary enclosing
//     expression is invalid in indeterminate forms.
//   - A composite expression containing an analyzed cut node returns a
//     structured Unimplemented instead of the silently edge-blind direct
//     substitution (was silent-wrong for bottom-edge approaches).
//   - Every undecidable step (non-literal real part, undecidable sign,
//     imaginary part not in a + b·I polynomial shape) bails to std::nullopt
//     → the previous engine path, changing nothing.
//
// Reference: Kahan (1987); Corless-Davenport-Jeffrey (2000) — see
// Branch_Cut_Propagation.md.

#include "cas/calculus.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "cas/unimplemented_info.hpp"
#include "calculus_internal.hpp"
#include "limit_internal.hpp"

#include <functional>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] bool is_imaginary_unit(ExprPtr e) {
    const auto* c = expr_cast<Constant>(e);
    return c != nullptr && c->value == MathConstant::I;
}

// Structural scan: does expr contain the imaginary unit anywhere?
[[nodiscard]] bool mentions_i(ExprPtr e) {
    if (!e) return false;
    if (is_imaginary_unit(e)) return true;
    if (expr_is<ComplexLit>(e)) return true;
    bool found = false;
    visit_expr(e, [&](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, Unary>) {
            found = mentions_i(node.operand);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            found = mentions_i(node.left) || mentions_i(node.right);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            for (ExprPtr a : node.args) if (mentions_i(a)) { found = true; break; }
        } else if constexpr (std::is_same_v<Node, Sum>) {
            for (ExprPtr t : node.terms) if (mentions_i(t)) { found = true; break; }
        } else if constexpr (std::is_same_v<Node, Product>) {
            for (ExprPtr f : node.factors) if (mentions_i(f)) { found = true; break; }
        }
    });
    return found;
}

// Cut-bearing function: real-axis cut on (−∞, 0). Arg is included because
// the simplifier's complex-log decomposition ln(z) = ln|z| + i·arg(z)
// relocates the ln cut into arg before the limit engine runs.
[[nodiscard]] bool has_negative_axis_cut(BuiltinOp op) {
    return op == BuiltinOp::Sqrt || op == BuiltinOp::Ln || op == BuiltinOp::Arg;
}

// Collect distinct FuncCall nodes f(w) with a negative-axis cut whose
// argument depends on `var` AND mentions I (fast pre-filter: without an
// imaginary component the path runs on the cut itself → principal branch,
// existing behaviour is correct).
void collect_cut_nodes(ExprPtr e, const Symbol& var, std::vector<ExprPtr>& out) {
    if (!e) return;
    if (const auto* fc = expr_cast<FuncCall>(e);
        fc != nullptr && fc->args.size() == 1U && has_negative_axis_cut(fc->func_id)
        && depends_on(fc->args[0], var) && mentions_i(fc->args[0])) {
        bool seen = false;
        for (ExprPtr x : out) if (x == e) { seen = true; break; }
        if (!seen) out.push_back(e);
        // Nested cuts inside the argument are collected too.
    }
    visit_expr(e, [&](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, Unary>) {
            collect_cut_nodes(node.operand, var, out);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            collect_cut_nodes(node.left, var, out);
            collect_cut_nodes(node.right, var, out);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            for (ExprPtr a : node.args) collect_cut_nodes(a, var, out);
        } else if constexpr (std::is_same_v<Node, Sum>) {
            for (ExprPtr t : node.terms) collect_cut_nodes(t, var, out);
        } else if constexpr (std::is_same_v<Node, Product>) {
            for (ExprPtr f : node.factors) collect_cut_nodes(f, var, out);
        }
    });
}

// Exact sign of a fully-simplified rational literal (incl. −literal);
// 0 when the value is zero or not exactly decidable (same contract as the
// scan in limit_rules.cpp).
[[nodiscard]] int exact_literal_sign(ExprPtr e) {
    if (const auto* il = expr_cast<IntegerLit>(e))
        return il->value.is_zero() ? 0 : (il->value.is_negative() ? -1 : 1);
    if (const auto* rl = expr_cast<RationalLit>(e))
        return rl->numerator.is_zero() ? 0 : (rl->numerator.is_negative() ? -1 : 1);
    if (const auto* u = expr_cast<Unary>(e); u != nullptr && u->op == UnaryOp::Neg) {
        const int inner = exact_literal_sign(u->operand);
        return -inner;
    }
    return 0;
}

// Split a SIMPLIFIED expression into (re, im) with w = re + im·I, purely
// structurally: Sum terms carrying an I factor go to im (accumulated — a
// multi-term imaginary part is kept whole), everything else to re. Any I
// found deeper than a top-level product factor → nullopt (undecidable).
[[nodiscard]] std::optional<std::pair<ExprPtr, ExprPtr>> split_re_im(
    ExprPtr w, AstArena& arena) {
    ExprPtr zero = arena.make<IntegerLit>(BigInt(0));

    // Declared as std::function for the single Neg-unwrapping recursion step.
    std::function<ExprPtr(ExprPtr)> imag_factor_of = [&](ExprPtr term) -> ExprPtr {
        if (is_imaginary_unit(term)) return arena.make<IntegerLit>(BigInt(1));
        if (const auto* u = expr_cast<Unary>(term);
            u != nullptr && u->op == UnaryOp::Neg) {
            if (is_imaginary_unit(u->operand)) return arena.make<IntegerLit>(BigInt(-1));
            ExprPtr inner = imag_factor_of(u->operand);
            if (inner == nullptr) return nullptr;
            return arena.make<Unary>(UnaryOp::Neg, inner);
        }
        if (const auto* prod = expr_cast<Product>(term)) {
            std::vector<ExprPtr> rest;
            bool found = false;
            for (ExprPtr f : prod->factors) {
                if (!found && is_imaginary_unit(f)) { found = true; continue; }
                rest.push_back(f);
            }
            if (!found) return nullptr;
            for (ExprPtr f : rest) {
                if (mentions_i(f)) return nullptr;  // I twice / nested → bail
            }
            if (rest.empty()) return arena.make<IntegerLit>(BigInt(1));
            if (rest.size() == 1U) return rest[0];
            return arena.make<Product>(std::move(rest));
        }
        return nullptr;
    };

    std::vector<ExprPtr> re_terms;
    std::vector<ExprPtr> im_terms;
    auto classify = [&](ExprPtr term) -> bool {
        if (ExprPtr b = imag_factor_of(term); b != nullptr) {
            im_terms.push_back(b);
            return true;
        }
        if (mentions_i(term)) return false;  // I hidden deeper → undecidable
        re_terms.push_back(term);
        return true;
    };

    if (const auto* sum = expr_cast<Sum>(w)) {
        for (ExprPtr term : sum->terms) {
            if (!classify(term)) return std::nullopt;
        }
    } else {
        if (!classify(w)) return std::nullopt;
    }

    auto join = [&](std::vector<ExprPtr>& v) -> ExprPtr {
        if (v.empty()) return zero;
        if (v.size() == 1U) return v[0];
        return arena.make<Sum>(std::move(v));
    };
    return std::make_pair(join(re_terms), join(im_terms));
}

}  // namespace

// Analyze one cut node argument. nullopt → undecidable, caller falls back.
std::optional<CutEdgeAnalysis> LimitEngine::analyze_cut_edge(
    ExprPtr arg, const Symbol& var, ExprPtr point) {
    auto w_simp = context_.simplify(arg);
    if (w_simp.is_error()) return std::nullopt;

    auto parts = split_re_im(w_simp.value(), arena_);
    if (!parts.has_value()) return std::nullopt;
    auto [re_w, im_w] = parts.value();

    // Real part at the point must be an exactly negative literal (on the cut).
    auto re0 = substitute_and_simplify(re_w, var, point);
    if (re0.is_error() || exact_literal_sign(re0.value()) != -1) return std::nullopt;

    // Im(w)(point) must vanish (we are AT the cut, not near it).
    auto im0 = substitute_and_simplify(im_w, var, point);
    if (im0.is_error() || !limit_is_zero(im0.value())) return std::nullopt;

    // Valuation scan: first k with (d^k/dt^k Im)(point) ≠ 0 decides the side.
    ExprPtr im_k = im_w;
    for (unsigned int k = 1U; k <= max_depth_budget_; ++k) {
        auto dk = diff(im_k, var, 1U, context_);
        if (dk.is_error()) return std::nullopt;
        im_k = dk.value();
        auto ck = substitute_and_simplify(im_k, var, point);
        if (ck.is_error()) return std::nullopt;
        if (limit_is_zero(ck.value())) continue;
        const int cs = exact_literal_sign(ck.value());
        if (cs == 0) return std::nullopt;  // undecidable sign → no guess

        CutEdgeAnalysis res;
        res.side_right = cs;
        res.side_left = ((k % 2U) == 1U) ? -cs : cs;
        ExprPtr neg = arena_.make<Unary>(UnaryOp::Neg, re0.value());
        auto neg_s = context_.simplify(neg);
        res.neg_re0 = neg_s.is_ok() ? neg_s.value() : neg;
        return res;
    }
    return std::nullopt;  // Im ≡ 0 up to budget → treat as on-cut real path
}

// Edge value of f on the cut at −neg_re0 < 0, from the given side (+1 top).
Result<ExprPtr> LimitEngine::cut_edge_value(
    BuiltinOp op, ExprPtr neg_re0, int side) {
    ExprPtr i_unit = arena_.make<Constant>(MathConstant::I);
    ExprPtr value = nullptr;
    if (op == BuiltinOp::Arg) {
        ExprPtr pi = arena_.make<Constant>(MathConstant::Pi);
        value = (side > 0)
            ? pi
            : static_cast<ExprPtr>(arena_.make<Unary>(UnaryOp::Neg, pi));
        return ok(value);
    }
    if (op == BuiltinOp::Sqrt) {
        ExprPtr root = arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{neg_re0});
        ExprPtr unsigned_val = arena_.make<Product>(std::vector<ExprPtr>{i_unit, root});
        value = (side > 0)
            ? unsigned_val
            : static_cast<ExprPtr>(arena_.make<Unary>(UnaryOp::Neg, unsigned_val));
    } else {  // Ln
        ExprPtr ln_abs = arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{neg_re0});
        ExprPtr pi = arena_.make<Constant>(MathConstant::Pi);
        ExprPtr ipi = arena_.make<Product>(std::vector<ExprPtr>{i_unit, pi});
        ExprPtr corr = (side > 0)
            ? ipi
            : static_cast<ExprPtr>(arena_.make<Unary>(UnaryOp::Neg, ipi));
        value = arena_.make<Sum>(std::vector<ExprPtr>{ln_abs, corr});
    }
    // No simplify here: inside the limit operation the simplifier runs with
    // operation_active_ and rewrites ln(4) → ln(16)/2, a form the equality
    // helpers cannot fold back. The caller canonicalizes on exit.
    return ok(value);
}

std::optional<Result<ExprPtr>> LimitEngine::try_branch_cut_directional(
    ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir,
    unsigned int depth) {
    // Fast pre-filters: finite point, expression mentioning both var and I.
    if (!expr || limit_is_infinity(point)) return std::nullopt;
    if (!mentions_i(expr)) return std::nullopt;

    std::vector<ExprPtr> cut_nodes;
    collect_cut_nodes(expr, var, cut_nodes);
    if (cut_nodes.empty()) return std::nullopt;

    // Case 1: the expression IS the cut node → direction-limit table.
    if (cut_nodes.size() == 1U && cut_nodes[0] == expr) {
        const auto& fc = expr_ref<FuncCall>(expr);
        auto edge = analyze_cut_edge(fc.args[0], var, point);
        if (!edge.has_value()) return std::nullopt;

        if (dir == LimitDirection::Both) {
            if (edge->side_right != edge->side_left) {
                return std::optional<Result<ExprPtr>>(fail<ExprPtr>(CASError{
                    .kind = CASErrorKind::Undefined,
                    .message = "Limite bilaterale attraverso il branch cut: i due "
                               "lati approcciano bordi opposti del taglio"}));
            }
            return std::optional<Result<ExprPtr>>(
                cut_edge_value(fc.func_id, edge->neg_re0, edge->side_right));
        }
        const int side = (dir == LimitDirection::Right) ? edge->side_right
                                                        : edge->side_left;
        return std::optional<Result<ExprPtr>>(
            cut_edge_value(fc.func_id, edge->neg_re0, side));
    }

    // Case 2: composite expression containing an ANALYZABLE cut subterm.
    // Plain substitution would silently return the top-edge (principal)
    // value regardless of direction, so it must not be trusted. Sound
    // decomposition by the algebra of limits: lim Σ = Σ lim and lim Π = Π lim
    // whenever every member limit exists FINITE (no indeterminate form is
    // possible among finite values). Members recurse through
    // compute_recursive, which re-enters this table on isolated cut nodes.
    bool any_analyzable = false;
    for (ExprPtr node : cut_nodes) {
        const auto& fc = expr_ref<FuncCall>(node);
        if (analyze_cut_edge(fc.args[0], var, point).has_value()) {
            any_analyzable = true;
            break;
        }
    }
    if (!any_analyzable) return std::nullopt;  // legacy path unchanged

    auto is_finite_value = [](ExprPtr v) -> bool {
        if (limit_is_infinity(v)) return false;
        if (const auto* c = expr_cast<Constant>(v)) {
            if (c->value == MathConstant::Indeterminate ||
                c->value == MathConstant::NaN) return false;
        }
        return true;
    };
    auto member_limits = [&](const std::vector<ExprPtr>& members)
        -> std::optional<std::vector<ExprPtr>> {
        std::vector<ExprPtr> out;
        out.reserve(members.size());
        for (ExprPtr m : members) {
            auto lm = compute_recursive(m, var, point, dir, depth + 1U);
            if (lm.is_error() || !is_finite_value(lm.value())) return std::nullopt;
            out.push_back(lm.value());
        }
        return out;
    };

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (auto vals = member_limits(sum->terms); vals.has_value()) {
            auto combined = context_.simplify(arena_.make<Sum>(std::move(*vals)));
            return std::optional<Result<ExprPtr>>(std::move(combined));
        }
    } else if (const auto* prod = expr_cast<Product>(expr)) {
        if (auto vals = member_limits(prod->factors); vals.has_value()) {
            auto combined = context_.simplify(arena_.make<Product>(std::move(*vals)));
            return std::optional<Result<ExprPtr>>(std::move(combined));
        }
    }

    // Not decomposable soundly (indeterminate member, ∞, or a non-Sum/Product
    // shape like a quotient through the cut): refuse explicitly (REGOLA ZERO:
    // incompleteness over silent-wrong).
    return std::optional<Result<ExprPtr>>(
        make_unimplemented<ExprPtr>(UnimplementedInfo{
            .module      = "calculus",
            .function    = "LimitEngine::try_branch_cut_directional",
            .input_shape = "indeterminate composite through a branch cut",
            .reason      = cas::error::reason_codes::SERIES_GENERAL,
            .suggestion  = "Decompose the limit so the cut function is the "
                           "outermost node, or resolve the indeterminate part",
            .ticket      = "HC-F8-PENDING-20"
        }));
}

}  // namespace cas::calculus
