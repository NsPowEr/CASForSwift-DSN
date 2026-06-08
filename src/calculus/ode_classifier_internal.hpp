#pragma once

#include "cas/ode.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <functional>

// Internal helpers shared between ode_classifier*.cpp translation units.
// NOT part of the public CAS API — do not include from include/cas/.

namespace cas::calculus {

[[nodiscard]] inline CASError ode_make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] inline bool is_zero_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    auto res = ctx.simplify(expr);
    if (res.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(res.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(res.value())) return rl->numerator.is_zero();
    return false;
}

[[nodiscard]] inline Result<ExprPtr> substitute_any(
    ExprPtr expr, ExprPtr target, ExprPtr replacement, symbolic::CASContext& ctx) {
    if (structural_equal(expr, target)) return ok(replacement);

    return visit_expr(expr, [&](const auto& node) -> Result<ExprPtr> {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, Unary>) {
            auto op = substitute_any(node.operand, target, replacement, ctx);
            if (op.is_error()) return op;
            return ok(ctx.arena().make<Unary>(node.op, op.value()));
        } else if constexpr (std::is_same_v<T, Binary>) {
            auto l = substitute_any(node.left, target, replacement, ctx);
            if (l.is_error()) return l;
            auto r = substitute_any(node.right, target, replacement, ctx);
            if (r.is_error()) return r;
            return ok(ctx.arena().make<Binary>(node.op, l.value(), r.value()));
        } else if constexpr (std::is_same_v<T, Sum>) {
            std::vector<ExprPtr> terms;
            for (auto t : node.terms) {
                auto r = substitute_any(t, target, replacement, ctx);
                if (r.is_error()) return r;
                terms.push_back(r.value());
            }
            return ok(ctx.arena().make<Sum>(std::move(terms)));
        } else if constexpr (std::is_same_v<T, Product>) {
            std::vector<ExprPtr> factors;
            for (auto f : node.factors) {
                auto r = substitute_any(f, target, replacement, ctx);
                if (r.is_error()) return r;
                factors.push_back(r.value());
            }
            return ok(ctx.arena().make<Product>(std::move(factors)));
        } else if constexpr (std::is_same_v<T, FuncCall>) {
            std::vector<ExprPtr> args;
            for (auto a : node.args) {
                auto r = substitute_any(a, target, replacement, ctx);
                if (r.is_error()) return r;
                args.push_back(r.value());
            }
            return ok(ctx.arena().make<FuncCall>(node.name, std::move(args)));
        } else if constexpr (std::is_same_v<T, Derivative>) {
            auto e = substitute_any(node.expression, target, replacement, ctx);
            if (e.is_error()) return e;
            return ok(ctx.arena().make<Derivative>(e.value(), Symbol(node.variable.name), node.order));
        } else {
            return ok(expr);
        }
    });
}

// Substitute `target` with `replacement` inside `expr`, but treat any
// Derivative subtree as an opaque atom.  Required when extracting
// polynomial-in-y coefficients from an ODE: a naïve recursion would also
// replace the inner Symbol(y) under D(y,x) and collapse the y' coefficient.
[[nodiscard]] inline Result<ExprPtr> substitute_y_shielded(
    ExprPtr expr, ExprPtr target, ExprPtr replacement, AstArena& arena) {
    if (structural_equal(expr, target)) return ok(replacement);
    return visit_expr(expr, [&](const auto& node) -> Result<ExprPtr> {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, Derivative>) {
            return ok(expr);  // opaque
        } else if constexpr (std::is_same_v<T, Unary>) {
            auto op = substitute_y_shielded(node.operand, target, replacement, arena);
            if (op.is_error()) return op;
            return ok(arena.make<Unary>(node.op, op.value()));
        } else if constexpr (std::is_same_v<T, Binary>) {
            auto l = substitute_y_shielded(node.left, target, replacement, arena);
            if (l.is_error()) return l;
            auto r = substitute_y_shielded(node.right, target, replacement, arena);
            if (r.is_error()) return r;
            return ok(arena.make<Binary>(node.op, l.value(), r.value()));
        } else if constexpr (std::is_same_v<T, Sum>) {
            std::vector<ExprPtr> terms;
            for (auto t : node.terms) {
                auto r = substitute_y_shielded(t, target, replacement, arena);
                if (r.is_error()) return r;
                terms.push_back(r.value());
            }
            return ok(arena.make<Sum>(std::move(terms)));
        } else if constexpr (std::is_same_v<T, Product>) {
            std::vector<ExprPtr> factors;
            for (auto f : node.factors) {
                auto r = substitute_y_shielded(f, target, replacement, arena);
                if (r.is_error()) return r;
                factors.push_back(r.value());
            }
            return ok(arena.make<Product>(std::move(factors)));
        } else if constexpr (std::is_same_v<T, FuncCall>) {
            std::vector<ExprPtr> args;
            for (auto a : node.args) {
                auto r = substitute_y_shielded(a, target, replacement, arena);
                if (r.is_error()) return r;
                args.push_back(r.value());
            }
            return ok(arena.make<FuncCall>(node.name, std::move(args)));
        } else {
            return ok(expr);
        }
    });
}

[[nodiscard]] inline std::optional<uint32_t> find_max_order(
    ExprPtr expr, const Symbol& y, const Symbol& x) {
    uint32_t max_order = 0;
    bool found = false;

    std::function<void(ExprPtr)> visit = [&](ExprPtr e) {
        if (const auto* deriv = expr_cast<Derivative>(e)) {
            if (const auto* sym = expr_cast<Symbol>(deriv->expression)) {
                if (sym->name == y.name && deriv->variable.name == x.name) {
                    if (deriv->order > max_order) max_order = deriv->order;
                    found = true;
                }
            }
        } else if (const auto* sym = expr_cast<Symbol>(e)) {
            if (sym->name == y.name) {
                found = true;
            }
        }

        if (const auto* un = expr_cast<Unary>(e)) visit(un->operand);
        else if (const auto* bin = expr_cast<Binary>(e)) { visit(bin->left); visit(bin->right); }
        else if (const auto* sum = expr_cast<Sum>(e)) { for (auto t : sum->terms) visit(t); }
        else if (const auto* prod = expr_cast<Product>(e)) { for (auto f : prod->factors) visit(f); }
        else if (const auto* func = expr_cast<FuncCall>(e)) { for (auto a : func->args) visit(a); }
        else if (const auto* intg = expr_cast<Integral>(e)) { visit(intg->integrand); }
        else if (const auto* der = expr_cast<Derivative>(e)) { visit(der->expression); }
    };
    visit(expr);

    if (found) return max_order;
    return std::nullopt;
}

[[nodiscard]] inline Result<ExprPtr> substitute_y_derivatives(
    ExprPtr E,
    const std::vector<ExprPtr>& y_ders,
    const std::vector<ExprPtr>& vals,
    symbolic::CASContext& ctx) {
    ExprPtr res = E;
    // Sostituzione in ordine inverso: partiamo dalle derivate più alte
    for (int i = static_cast<int>(y_ders.size()) - 1; i >= 0; --i) {
        auto sub_res = substitute_any(res, y_ders[i], vals[i], ctx);
        if (sub_res.is_error()) return sub_res;
        res = sub_res.value();
    }
    return ctx.simplify(res);
}

struct FirstOrderLinearInYPrime {
    ExprPtr A; // Coefficient of y'
    ExprPtr B; // Constant term (w.r.t y')
};

inline std::optional<FirstOrderLinearInYPrime> extract_first_order_yprime(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto eval_yp = [&](long long k) -> Result<ExprPtr> {
        ExprPtr val = arena.make<IntegerLit>(BigInt(k));
        auto s = substitute_any(E, y_prime, val, ctx);
        if (s.is_error()) return s;
        return ctx.simplify(s.value());
    };

    auto e0 = eval_yp(0); if (e0.is_error()) return std::nullopt;
    auto e1 = eval_yp(1); if (e1.is_error()) return std::nullopt;
    auto e2 = eval_yp(2); if (e2.is_error()) return std::nullopt;

    ExprPtr B = e0.value();
    auto A_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, e1.value(), B));
    if (A_raw.is_error()) return std::nullopt;
    ExprPtr A = A_raw.value();

    auto two = arena.make<IntegerLit>(BigInt(2));
    auto expected_e2 = ctx.simplify(arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, two, A), B));
    if (expected_e2.is_error()) return std::nullopt;
    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, e2.value(), expected_e2.value()), ctx);
    if (diff.is_error() || !is_zero_expr(diff.value(), ctx)) return std::nullopt;

    auto depends_on_yp = [&](ExprPtr e) {
        auto mo = find_max_order(e, y, x);
        return mo.has_value() && *mo >= 1U;
    };
    if (depends_on_yp(A) || depends_on_yp(B)) return std::nullopt;

    return FirstOrderLinearInYPrime{A, B};
}

// Helper used in try_bernoulli: extract the y-power from a monomial term.
inline std::optional<ExprPtr> get_y_power(ExprPtr term, const Symbol& y, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (!depends_on(term, y)) return arena.make<IntegerLit>(BigInt(0));

    if (const auto* un = expr_cast<Unary>(term)) {
        if (un->op == UnaryOp::Neg) {
            return get_y_power(un->operand, y, ctx);
        }
    }

    if (const auto* prod = expr_cast<Product>(term)) {
        ExprPtr y_pow = nullptr;
        for (auto f : prod->factors) {
            if (depends_on(f, y)) {
                if (y_pow) return std::nullopt;
                if (const auto* sym = expr_cast<Symbol>(f)) {
                    if (sym->name == y.name) y_pow = arena.make<IntegerLit>(BigInt(1));
                    else return std::nullopt;
                } else if (const auto* bin = expr_cast<Binary>(f)) {
                    if (bin->op == BinaryOp::Pow) {
                        if (const auto* sym = expr_cast<Symbol>(bin->left); sym && sym->name == y.name) {
                            if (!depends_on(bin->right, y)) y_pow = bin->right;
                            else return std::nullopt;
                        } else return std::nullopt;
                    } else return std::nullopt;
                } else return std::nullopt;
            }
        }
        return y_pow;
    } else {
        if (const auto* sym = expr_cast<Symbol>(term)) {
            if (sym->name == y.name) return arena.make<IntegerLit>(BigInt(1));
        } else if (const auto* bin = expr_cast<Binary>(term)) {
            if (bin->op == BinaryOp::Pow) {
                if (const auto* sym = expr_cast<Symbol>(bin->left); sym && sym->name == y.name) {
                    if (!depends_on(bin->right, y)) return bin->right;
                }
            }
        }
    }
    return std::nullopt;
}

// ---------------------------------------------------------------------------
// Forward declarations of ODE type detectors (defined in _first_order.cpp
// and _higher_order.cpp).  Called from classify_ode() in ode_classifier.cpp.
// ---------------------------------------------------------------------------

[[nodiscard]] std::optional<OdeClassification> try_separable(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime,
    symbolic::CASContext& ctx);

[[nodiscard]] std::optional<OdeClassification> try_exact(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime,
    symbolic::CASContext& ctx);

[[nodiscard]] std::optional<OdeClassification> try_bernoulli(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime,
    symbolic::CASContext& ctx);

[[nodiscard]] std::optional<OdeClassification> try_homogeneous(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime,
    symbolic::CASContext& ctx);

[[nodiscard]] std::optional<OdeClassification> try_riccati(
    ExprPtr E, const Symbol& y, const Symbol& x,
    ExprPtr y_symbol, ExprPtr y_prime,
    symbolic::CASContext& ctx);

[[nodiscard]] std::optional<OdeClassification> try_lagrange(
    ExprPtr E, const Symbol& y, const Symbol& x,
    ExprPtr y_symbol, ExprPtr y_prime,
    symbolic::CASContext& ctx);

} // namespace cas::calculus
