// limit_mrv_compare.cpp — Gruntz §3.5 asymptotic growth comparison.
//
// Provides:
//   - poly_degree_wrt    (file-local helper)
//   - compare_growth     (external linkage, declared in calculus_internal.hpp)

#include "limit_mrv_internal.hpp"

#include <optional>
#include <vector>

namespace cas::calculus {

namespace {

// Returns polynomial degree of e w.r.t. var, or nullopt if not a polynomial.
std::optional<int> poly_degree_wrt(ExprPtr e, const Symbol& var) {
    if (!depends_on(e, var)) return 0;
    if (const auto* sym = expr_cast<Symbol>(e)) {
        return (sym->name == var.name) ? 1 : 0;
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        int deg = 0;
        for (auto t : sum->terms) {
            auto d = poly_degree_wrt(t, var);
            if (!d) return std::nullopt;
            deg = std::max(deg, *d);
        }
        return deg;
    }
    if (const auto* product = expr_cast<Product>(e)) {
        int deg = 0;
        for (auto f : product->factors) {
            auto d = poly_degree_wrt(f, var);
            if (!d) return std::nullopt;
            deg += *d;
        }
        return deg;
    }
    if (const auto* binary = expr_cast<Binary>(e)) {
        if (binary->op == BinaryOp::Pow) {
            auto base_deg = poly_degree_wrt(binary->left, var);
            if (!base_deg) return std::nullopt;
            if (*base_deg == 0) return 0;
            // base depends on var — need non-negative integer exponent
            if (const auto* exp_lit = expr_cast<IntegerLit>(binary->right)) {
                if (!exp_lit->value.is_negative()) {
                    auto eu = exp_lit->value.to_u64();
                    if (eu > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) return std::nullopt;
                    return *base_deg * static_cast<int>(eu);
                }
            }
            return std::nullopt;
        }
        if (binary->op == BinaryOp::Add) {
            auto l = poly_degree_wrt(binary->left, var);
            auto r = poly_degree_wrt(binary->right, var);
            if (!l || !r) return std::nullopt;
            return std::max(*l, *r);
        }
        if (binary->op == BinaryOp::Mul) {
            auto l = poly_degree_wrt(binary->left, var);
            auto r = poly_degree_wrt(binary->right, var);
            if (!l || !r) return std::nullopt;
            return *l + *r;
        }
    }
    return std::nullopt;  // FuncCall or unrecognized — not polynomial
}

} // namespace

// Growth comparison for x -> +infinity.
// Returns: +1 if a grows faster, -1 if b grows faster, 0 if same rate or undecidable.
// This is deliberately structural and recursive: exp(exp(x)) must dominate exp(x^n)
// because exp(x) dominates every polynomial x^n.
int compare_growth(ExprPtr a, ExprPtr b, const Symbol& var, symbolic::CASContext& ctx) {
    if (structural_equal(a, b)) return 0;

    auto positive_growth_part = [](ExprPtr e) -> ExprPtr {
        if (const auto* unary = expr_cast<Unary>(e)) {
            if (unary->op == UnaryOp::Neg) return unary->operand;
        }
        return e;
    };

    a = positive_growth_part(a);
    b = positive_growth_part(b);
    if (structural_equal(a, b)) return 0;

    // Dynamic rank matching limit_infinite.cpp::get_growth_rank (Cat 10):
    //   exp(arg) -> rank(arg) + 1  (so nested towers grow strictly)
    //   log(arg) -> max(rank(arg) - 1, 0)
    //
    // F2.1.a: depth bound derived from AST nesting. The recursion descends
    // strictly one AST level per call; an expression of N nodes can have
    // recursion depth ≤ N. Cap prevents stack overflow on
    // pathological or self-referential trees (return 0 = "unknown growth
    // class" which is conservative and forces the caller to bail).
    const int growth_rank_max_depth = ctx.mrv_growth_rank_max_depth();
    auto get_growth_rank_impl = [&](ExprPtr e, const auto& self, int depth) -> int {
        if (depth >= growth_rank_max_depth) return 0;
        if (!depends_on(e, var)) return 0;
        if (const auto* unary = expr_cast<Unary>(e)) {
            if (unary->op == UnaryOp::Neg) return self(unary->operand, self, depth + 1);
        }
        if (const auto* call = expr_cast<FuncCall>(e); call != nullptr && call->args.size() == 1U) {
            if (call->func_id == BuiltinOp::Exp) {
                auto arg_limit = try_infinite_limit(
                    call->args.front(),
                    var,
                    ctx.arena().make<Constant>(MathConstant::Infinity),
                    ctx);
                // exp of an arg that goes to -infinity decays to 0: rank 0.
                if (arg_limit.is_ok() && limit_is_infinity(arg_limit.value())
                    && expr_is<Unary>(arg_limit.value())) {
                    return 0;
                }
                int inner = self(call->args.front(), self, depth + 1);
                return inner + 1;
            }
            if (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) {
                int inner = self(call->args.front(), self, depth + 1);
                return inner > 1 ? inner - 1 : 0;
            }
        }
        if (const auto* product = expr_cast<Product>(e)) {
            int rank = 0;
            for (ExprPtr factor : product->factors) {
                rank = std::max(rank, self(factor, self, depth + 1));
            }
            return rank;
        }
        if (const auto* sum = expr_cast<Sum>(e)) {
            int rank = 0;
            for (ExprPtr term : sum->terms) {
                rank = std::max(rank, self(term, self, depth + 1));
            }
            return rank;
        }
        if (const auto* binary = expr_cast<Binary>(e)) {
            if (binary->op == BinaryOp::Mul || binary->op == BinaryOp::Div) {
                return std::max(self(binary->left, self, depth + 1),
                                self(binary->right, self, depth + 1));
            }
            if (binary->op == BinaryOp::Add || binary->op == BinaryOp::Sub) {
                return std::max(self(binary->left, self, depth + 1),
                                self(binary->right, self, depth + 1));
            }
            if (binary->op == BinaryOp::Pow) {
                const bool exp_dep = depends_on(binary->right, var);
                const bool base_dep = depends_on(binary->left, var);
                if (!exp_dep) return self(binary->left, self, depth + 1);
                if (!base_dep) return self(binary->right, self, depth + 1) + 1;
                int g = self(binary->right, self, depth + 1);
                int ln_f = self(binary->left, self, depth + 1);
                return g + (ln_f > 1 ? ln_f - 1 : 0) + 1;
            }
        }
        return 2;  // polynomial or mixed expression depending on var
    };
    auto get_growth_rank = [&](ExprPtr e) -> int {
        return get_growth_rank_impl(e, get_growth_rank_impl, 0);
    };

    int ra = get_growth_rank(a);
    int rb = get_growth_rank(b);

    if (ra != rb) return ra > rb ? 1 : -1;

    // Same coarse rank: use polynomial degree as tiebreak when both are polynomial.
    if (ra == 2) {
        auto da = poly_degree_wrt(a, var);
        auto db = poly_degree_wrt(b, var);
        if (da && db) {
            if (*da != *db) return *da > *db ? 1 : -1;
            return 0;
        }
    }

    // exp(f) vs exp(g) is decided by f vs g, recursively.  Applies at any
    // depth in the exponential tower (ra == rb >= 3 with both exp).
    if (ra >= 3) {
        const auto* ca = expr_cast<FuncCall>(a);
        const auto* cb = expr_cast<FuncCall>(b);
        if (ca && cb && ca->func_id == BuiltinOp::Exp && cb->func_id == BuiltinOp::Exp
            && !ca->args.empty() && !cb->args.empty()) {
            return compare_growth(ca->args.front(), cb->args.front(), var, ctx);
        }
    }

    return 0;
}

} // namespace cas::calculus
