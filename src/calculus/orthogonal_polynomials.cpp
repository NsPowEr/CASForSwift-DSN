// L3-04 closure: orthogonality patterns for the classical orthogonal
// polynomial families.
//
//   Legendre P_n on [-1, 1]:        ∫ P_m(x) P_n(x) dx = 2/(2n+1) · δ_{mn}
//   Hermite physicist H_n on R:     ∫ H_m(x) H_n(x) e^{-x²} dx = 2ⁿ·n!·√π · δ_{mn}
//   Hermite probabilist He_n on R:  ∫ He_m(x) He_n(x) e^{-x²/2} dx = n!·√(2π) · δ_{mn}
//
// Each pattern matches the canonical integrand shape (product of two
// polynomial calls, plus the appropriate weight for the Hermite cases),
// extracts the integer indices m and n, and returns the closed-form
// value.  Mismatches simply propagate as nullopt so the next pattern in
// the registry can try.

#include "integrate_definite_patterns.hpp"

#include "cas/builtin_functions.hpp"
#include "cas/error.hpp"

#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

// ─── shared helpers ────────────────────────────────────────────────────────

void flatten_mul_factors(ExprPtr expr, std::vector<ExprPtr>& out) {
    if (!expr) return;
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr f : product->factors) flatten_mul_factors(f, out);
        return;
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Mul) {
            flatten_mul_factors(binary->left, out);
            flatten_mul_factors(binary->right, out);
            return;
        }
    }
    out.push_back(expr);
}

[[nodiscard]] bool is_pos_infinity(ExprPtr expr) {
    const auto* c = expr_cast<Constant>(expr);
    return c != nullptr && c->value == MathConstant::Infinity;
}

[[nodiscard]] bool is_neg_infinity(ExprPtr expr) {
    const auto* u = expr_cast<Unary>(expr);
    return u != nullptr && u->op == UnaryOp::Neg && is_pos_infinity(u->operand);
}

[[nodiscard]] bool depends_on_var(ExprPtr expr, const Symbol& var) {
    if (!expr) return false;
    if (const auto* sym = expr_cast<Symbol>(expr)) return sym->name == var.name;
    bool found = false;
    visit_expr(expr, [&](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, Unary>) {
            if (depends_on_var(node.operand, var)) found = true;
        } else if constexpr (std::is_same_v<Node, Binary>) {
            if (depends_on_var(node.left, var)) found = true;
            if (depends_on_var(node.right, var)) found = true;
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            for (ExprPtr arg : node.args)
                if (depends_on_var(arg, var)) { found = true; break; }
        } else if constexpr (std::is_same_v<Node, Sum>) {
            for (ExprPtr term : node.terms)
                if (depends_on_var(term, var)) { found = true; break; }
        } else if constexpr (std::is_same_v<Node, Product>) {
            for (ExprPtr fac : node.factors)
                if (depends_on_var(fac, var)) { found = true; break; }
        }
    });
    return found;
}

// Match a FuncCall(op, [idx_lit, var_symbol]) and return the integer index.
// nullopt if the call shape does not fit.
[[nodiscard]] std::optional<BigInt> match_poly_call(
    ExprPtr expr, BuiltinOp op, const Symbol& var) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (!call || call->func_id != op || call->args.size() != 2U) return std::nullopt;
    const auto* idx_lit = expr_cast<IntegerLit>(call->args[0]);
    if (!idx_lit) return std::nullopt;
    if (idx_lit->value.is_negative()) return std::nullopt;
    const auto* arg_sym = expr_cast<Symbol>(call->args[1]);
    if (!arg_sym || arg_sym->name != var.name) return std::nullopt;
    return idx_lit->value;
}

[[nodiscard]] bool is_literal_rational(ExprPtr expr, long long num, long long den) {
    if (den == 1) {
        const auto* il = expr_cast<IntegerLit>(expr);
        return il != nullptr && il->value == BigInt(num);
    }
    const auto* rl = expr_cast<RationalLit>(expr);
    return rl != nullptr && rl->numerator == BigInt(num) && rl->denominator == BigInt(den);
}

[[nodiscard]] bool is_literal_neg_one(ExprPtr expr) {
    if (const auto* il = expr_cast<IntegerLit>(expr); il && il->value == BigInt(-1)) return true;
    if (const auto* un = expr_cast<Unary>(expr); un && un->op == UnaryOp::Neg) {
        const auto* il = expr_cast<IntegerLit>(un->operand);
        if (il && il->value == BigInt(1)) return true;
    }
    return false;
}

// Match exp(−x²) where x is the integration variable.  Used by HermiteH.
[[nodiscard]] bool match_exp_neg_xsq(ExprPtr expr, const Symbol& var) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (!call || call->func_id != BuiltinOp::Exp || call->args.size() != 1U) return false;
    ExprPtr arg = call->args.front();
    const auto* un = expr_cast<Unary>(arg);
    if (!un || un->op != UnaryOp::Neg) return false;
    const auto* pw = expr_cast<Binary>(un->operand);
    if (!pw || pw->op != BinaryOp::Pow) return false;
    const auto* sym = expr_cast<Symbol>(pw->left);
    if (!sym || sym->name != var.name) return false;
    const auto* exp2 = expr_cast<IntegerLit>(pw->right);
    return exp2 != nullptr && exp2->value == BigInt(2);
}

// Helper: does `expr` equal x² where x is `var`?
[[nodiscard]] bool is_var_squared(ExprPtr expr, const Symbol& var) {
    const auto* pw = expr_cast<Binary>(expr);
    if (!pw || pw->op != BinaryOp::Pow) return false;
    const auto* sym = expr_cast<Symbol>(pw->left);
    if (!sym || sym->name != var.name) return false;
    const auto* e2 = expr_cast<IntegerLit>(pw->right);
    return e2 != nullptr && e2->value == BigInt(2);
}

// Match exp(−x²/2).  Used by HermiteHe.
[[nodiscard]] bool match_exp_neg_xsq_over_2(ExprPtr expr, const Symbol& var) {
    const auto* call = expr_cast<FuncCall>(expr);
    if (!call || call->func_id != BuiltinOp::Exp || call->args.size() != 1U) return false;
    ExprPtr arg = call->args.front();
    // Accept Binary(Div, Unary(Neg, x²), 2)  ← parser output for "exp(-(x^2)/2)"
    if (const auto* div = expr_cast<Binary>(arg); div && div->op == BinaryOp::Div) {
        const auto* two = expr_cast<IntegerLit>(div->right);
        if (two && two->value == BigInt(2)) {
            if (const auto* neg = expr_cast<Unary>(div->left); neg && neg->op == UnaryOp::Neg) {
                if (is_var_squared(neg->operand, var)) return true;
            }
        }
    }
    // Accept Unary(Neg, Binary(Div, x², 2))
    if (const auto* un = expr_cast<Unary>(arg); un && un->op == UnaryOp::Neg) {
        if (const auto* div = expr_cast<Binary>(un->operand); div && div->op == BinaryOp::Div) {
            const auto* pw = expr_cast<Binary>(div->left);
            if (!pw || pw->op != BinaryOp::Pow) return false;
            const auto* sym = expr_cast<Symbol>(pw->left);
            if (!sym || sym->name != var.name) return false;
            const auto* exp2 = expr_cast<IntegerLit>(pw->right);
            if (!exp2 || exp2->value != BigInt(2)) return false;
            const auto* two = expr_cast<IntegerLit>(div->right);
            return two != nullptr && two->value == BigInt(2);
        }
    }
    if (const auto* mul = expr_cast<Binary>(arg); mul && mul->op == BinaryOp::Mul) {
        auto try_match = [&](ExprPtr coeff, ExprPtr pow_expr) -> bool {
            if (!is_literal_rational(coeff, -1, 2)) return false;
            const auto* pw = expr_cast<Binary>(pow_expr);
            if (!pw || pw->op != BinaryOp::Pow) return false;
            const auto* sym = expr_cast<Symbol>(pw->left);
            if (!sym || sym->name != var.name) return false;
            const auto* exp2 = expr_cast<IntegerLit>(pw->right);
            return exp2 != nullptr && exp2->value == BigInt(2);
        };
        return try_match(mul->left, mul->right) || try_match(mul->right, mul->left);
    }
    return false;
}

[[nodiscard]] BigInt factorial_bigint(const BigInt& n) {
    BigInt acc(1);
    BigInt i(2);
    while (i <= n) { acc = acc * i; i = i + BigInt(1); }
    return acc;
}

[[nodiscard]] BigInt pow2_bigint(const BigInt& n) {
    BigInt acc(1);
    BigInt i(0);
    while (i < n) { acc = acc * BigInt(2); i = i + BigInt(1); }
    return acc;
}

// Build the diagonal Legendre normalization 2/(2n+1).
[[nodiscard]] ExprPtr build_legendre_norm(const BigInt& n, AstArena& arena) {
    BigInt denom = BigInt(2) * n + BigInt(1);
    return arena.make<RationalLit>(BigInt(2), denom);
}

// Build the diagonal Hermite H normalization 2^n · n! · √π.
[[nodiscard]] ExprPtr build_hermite_h_norm(const BigInt& n, AstArena& arena) {
    const BigInt pow_two = pow2_bigint(n);
    const BigInt fact = factorial_bigint(n);
    const BigInt scalar = pow_two * fact;
    ExprPtr scalar_expr = arena.make<IntegerLit>(scalar);
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr sqrt_pi = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{pi});
    return arena.make<Binary>(BinaryOp::Mul, scalar_expr, sqrt_pi);
}

// Build the diagonal Hermite He normalization n! · √(2π).
[[nodiscard]] ExprPtr build_hermite_he_norm(const BigInt& n, AstArena& arena) {
    const BigInt fact = factorial_bigint(n);
    ExprPtr fact_expr = arena.make<IntegerLit>(fact);
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr two_pi = arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(2)), pi);
    ExprPtr sqrt_two_pi = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{two_pi});
    return arena.make<Binary>(BinaryOp::Mul, fact_expr, sqrt_two_pi);
}

// Detect a product-of-two-poly-calls pattern.  Returns (m, n) extracted
// indices and the residual non-poly factors (which must all be the
// weight if any, plus possibly numeric constants).
struct PolyProductMatch {
    BigInt m;
    BigInt n;
    std::vector<ExprPtr> other_factors;       // factors that are not the two poly calls
};

[[nodiscard]] std::optional<PolyProductMatch> match_two_poly_product(
    ExprPtr integrand_normalized,
    BuiltinOp op,
    const Symbol& var) {
    std::vector<ExprPtr> factors;
    flatten_mul_factors(integrand_normalized, factors);
    std::optional<BigInt> first_idx;
    std::optional<BigInt> second_idx;
    PolyProductMatch out;
    for (ExprPtr f : factors) {
        if (!first_idx.has_value()) {
            if (auto idx = match_poly_call(f, op, var)) { first_idx = idx; continue; }
        } else if (!second_idx.has_value()) {
            if (auto idx = match_poly_call(f, op, var)) { second_idx = idx; continue; }
        }
        // Pow(P_n(x), 2) → both indices identical, base contributes twice.
        if (const auto* binary = expr_cast<Binary>(f); binary && binary->op == BinaryOp::Pow) {
            if (auto idx = match_poly_call(binary->left, op, var)) {
                const auto* exp_lit = expr_cast<IntegerLit>(binary->right);
                if (exp_lit && exp_lit->value == BigInt(2) && !first_idx.has_value()) {
                    first_idx = idx;
                    second_idx = idx;
                    continue;
                }
            }
        }
        out.other_factors.push_back(f);
    }
    if (!first_idx.has_value() || !second_idx.has_value()) return std::nullopt;
    out.m = first_idx.value();
    out.n = second_idx.value();
    return out;
}

// Multiply `base` by every constant non-var factor in `others`.  Var-dependent
// factors disqualify the match (caller must check).
[[nodiscard]] Result<std::optional<ExprPtr>> apply_constant_factors(
    ExprPtr base,
    const std::vector<ExprPtr>& others,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr acc = base;
    for (ExprPtr f : others) {
        if (depends_on_var(f, var)) return ok(std::optional<ExprPtr>{});
        acc = arena.make<Binary>(BinaryOp::Mul, acc, f);
    }
    auto simp = ctx.simplify(acc);
    if (simp.is_error()) return fail<std::optional<ExprPtr>>(simp.error());
    return ok(std::optional<ExprPtr>(simp.value()));
}

}  // namespace

[[nodiscard]] Result<std::optional<ExprPtr>> pattern_legendre_orthogonality(const DefiniteContext& dc) {
    // Limits must be exactly [-1, 1].
    if (!is_literal_neg_one(dc.lower)) return ok(std::optional<ExprPtr>{});
    if (!is_literal_rational(dc.upper, 1, 1)) return ok(std::optional<ExprPtr>{});

    auto match = match_two_poly_product(dc.integrand, BuiltinOp::LegendreP, dc.var);
    if (!match.has_value()) return ok(std::optional<ExprPtr>{});

    AstArena& arena = dc.ctx.arena();
    ExprPtr base;
    if (match->m == match->n) {
        base = build_legendre_norm(match->n, arena);
    } else {
        base = arena.make<IntegerLit>(BigInt(0));
    }
    return apply_constant_factors(base, match->other_factors, dc.var, dc.ctx);
}

[[nodiscard]] Result<std::optional<ExprPtr>> pattern_hermite_h_orthogonality(const DefiniteContext& dc) {
    if (!is_neg_infinity(dc.lower) || !is_pos_infinity(dc.upper)) return ok(std::optional<ExprPtr>{});

    // Hermite H product needs an explicit exp(-x²) weight factor.
    std::vector<ExprPtr> factors;
    flatten_mul_factors(dc.integrand, factors);
    std::vector<ExprPtr> non_weight;
    bool weight_found = false;
    for (ExprPtr f : factors) {
        if (!weight_found && match_exp_neg_xsq(f, dc.var)) {
            weight_found = true;
            continue;
        }
        non_weight.push_back(f);
    }
    if (!weight_found) return ok(std::optional<ExprPtr>{});

    // Run the poly-product matcher directly on the flattened non-weight
    // factors — DO NOT call ctx.simplify(), which would expand HermiteH
    // into its polynomial form and destroy the pattern.
    AstArena& arena = dc.ctx.arena();
    ExprPtr residual = arena.make<IntegerLit>(BigInt(1));
    if (!non_weight.empty()) {
        if (non_weight.size() == 1U) residual = non_weight.front();
        else residual = arena.make<Product>(std::move(non_weight));
    }

    auto match = match_two_poly_product(residual, BuiltinOp::HermiteH, dc.var);
    if (!match.has_value()) return ok(std::optional<ExprPtr>{});

    ExprPtr base;
    if (match->m == match->n) {
        base = build_hermite_h_norm(match->n, arena);
    } else {
        base = arena.make<IntegerLit>(BigInt(0));
    }
    return apply_constant_factors(base, match->other_factors, dc.var, dc.ctx);
}

[[nodiscard]] Result<std::optional<ExprPtr>> pattern_hermite_he_orthogonality(const DefiniteContext& dc) {
    if (!is_neg_infinity(dc.lower) || !is_pos_infinity(dc.upper)) return ok(std::optional<ExprPtr>{});

    std::vector<ExprPtr> factors;
    flatten_mul_factors(dc.integrand, factors);
    std::vector<ExprPtr> non_weight;
    bool weight_found = false;
    for (ExprPtr f : factors) {
        if (!weight_found && match_exp_neg_xsq_over_2(f, dc.var)) {
            weight_found = true;
            continue;
        }
        non_weight.push_back(f);
    }
    if (!weight_found) return ok(std::optional<ExprPtr>{});

    AstArena& arena = dc.ctx.arena();
    ExprPtr residual = arena.make<IntegerLit>(BigInt(1));
    if (!non_weight.empty()) {
        if (non_weight.size() == 1U) residual = non_weight.front();
        else residual = arena.make<Product>(std::move(non_weight));
    }
    auto match = match_two_poly_product(residual, BuiltinOp::HermiteHe, dc.var);
    if (!match.has_value()) return ok(std::optional<ExprPtr>{});

    ExprPtr base;
    if (match->m == match->n) {
        base = build_hermite_he_norm(match->n, arena);
    } else {
        base = arena.make<IntegerLit>(BigInt(0));
    }
    return apply_constant_factors(base, match->other_factors, dc.var, dc.ctx);
}

}  // namespace cas::calculus
