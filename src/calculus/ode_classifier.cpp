#include "cas/ode.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <functional>

namespace cas::calculus {

[[nodiscard]] static CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] static bool is_zero_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    auto res = ctx.simplify(expr);
    if (res.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(res.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(res.value())) return rl->numerator.is_zero();
    return false;
}

[[nodiscard]] static Result<ExprPtr> substitute_any(ExprPtr expr, ExprPtr target, ExprPtr replacement, symbolic::CASContext& ctx) {
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
[[nodiscard]] static Result<ExprPtr> substitute_y_shielded(
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

[[nodiscard]] static std::optional<uint32_t> find_max_order(ExprPtr expr, const Symbol& y, const Symbol& x) {
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

[[nodiscard]] static Result<ExprPtr> substitute_y_derivatives(ExprPtr E, const std::vector<ExprPtr>& y_ders, const std::vector<ExprPtr>& vals, symbolic::CASContext& ctx) {
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

static std::optional<FirstOrderLinearInYPrime> extract_first_order_yprime(
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

[[nodiscard]] static std::optional<OdeClassification> try_separable(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    auto lin = extract_first_order_yprime(E, y, x, y_prime, ctx);
    if (!lin) return std::nullopt;
    ExprPtr A = lin->A;
    ExprPtr B = lin->B;
    
    if (is_zero_expr(A, ctx)) return std::nullopt;
    
    AstArena& arena = ctx.arena();
    
    if (is_zero_expr(B, ctx)) {
        OdeClassification res(OdeType::Separable, E, y, x);
        res.components.push_back(arena.make<IntegerLit>(BigInt(1)));
        res.components.push_back(arena.make<IntegerLit>(BigInt(0)));
        return res;
    }

    if (!depends_on(A, x) && !depends_on(B, y)) {
        OdeClassification res(OdeType::Separable, E, y, x);
        res.components.push_back(A);
        res.components.push_back(ctx.simplify(arena.make<Unary>(UnaryOp::Neg, B)).value());
        return res;
    }

    std::pair<long long, long long> points[] = {{0,0}, {1,1}, {0,1}, {1,0}, {-1,-1}, {2,2}};
    for (auto [x0, y0] : points) {
        auto eval_xy = [&](ExprPtr expr, long long xv, long long yv) -> Result<ExprPtr> {
            auto s1 = substitute_any(expr, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(xv)), ctx);
            if (s1.is_error()) return s1;
            auto s2 = substitute_any(s1.value(), arena.make<Symbol>(y.name), arena.make<IntegerLit>(BigInt(yv)), ctx);
            if (s2.is_error()) return s2;
            return ctx.simplify(s2.value());
        };
        
        auto A0_res = eval_xy(A, x0, y0);
        auto B0_res = eval_xy(B, x0, y0);
        if (A0_res.is_error() || B0_res.is_error()) continue;
        ExprPtr A0 = A0_res.value();
        ExprPtr B0 = B0_res.value();
        
        if (is_zero_expr(A0, ctx) || is_zero_expr(B0, ctx)) continue;
        
        auto Ax_y0_res = substitute_any(A, arena.make<Symbol>(y.name), arena.make<IntegerLit>(BigInt(y0)), ctx);
        auto Bx_y0_res = substitute_any(B, arena.make<Symbol>(y.name), arena.make<IntegerLit>(BigInt(y0)), ctx);
        if (Ax_y0_res.is_error() || Bx_y0_res.is_error()) continue;
        ExprPtr Ax_y0 = ctx.simplify(Ax_y0_res.value()).value();
        ExprPtr Bx_y0 = ctx.simplify(Bx_y0_res.value()).value();
        
        auto Ax0_y_res = substitute_any(A, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(x0)), ctx);
        auto Bx0_y_res = substitute_any(B, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(x0)), ctx);
        if (Ax0_y_res.is_error() || Bx0_y_res.is_error()) continue;
        ExprPtr Ax0_y = ctx.simplify(Ax0_y_res.value()).value();
        ExprPtr Bx0_y = ctx.simplify(Bx0_y_res.value()).value();
        
        auto term1 = arena.make<Binary>(BinaryOp::Mul, A, 
            arena.make<Binary>(BinaryOp::Mul, Bx_y0, 
                arena.make<Binary>(BinaryOp::Mul, Bx0_y, A0)));
        auto term2 = arena.make<Binary>(BinaryOp::Mul, B, 
            arena.make<Binary>(BinaryOp::Mul, Ax0_y, 
                arena.make<Binary>(BinaryOp::Mul, B0, Ax_y0)));
        
        auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, term1, term2), ctx);
        if (diff.is_ok() && is_zero_expr(diff.value(), ctx)) {
            auto M_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div, 
                arena.make<Unary>(UnaryOp::Neg, Bx_y0), Ax_y0));
            auto N_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div,
                arena.make<Binary>(BinaryOp::Mul, Ax0_y, B0),
                arena.make<Binary>(BinaryOp::Mul, Bx0_y, A0)));
                
            if (M_raw.is_ok() && N_raw.is_ok()) {
                OdeClassification res(OdeType::Separable, E, y, x);
                res.components.push_back(N_raw.value());
                res.components.push_back(M_raw.value());
                return res;
            }
        }
    }
    
    return std::nullopt;
}

[[nodiscard]] static std::optional<OdeClassification> try_exact(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    auto lin = extract_first_order_yprime(E, y, x, y_prime, ctx);
    if (!lin) return std::nullopt;
    ExprPtr N = lin->A;
    ExprPtr M = lin->B;
    

    
    AstArena& arena = ctx.arena();
    
    auto dM_dy = cas::calculus::diff(M, Symbol(y.name), 1, ctx);
    auto dN_dx = cas::calculus::diff(N, Symbol(x.name), 1, ctx);
    if (dM_dy.is_error() || dN_dx.is_error()) return std::nullopt;
    
    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, dM_dy.value(), dN_dx.value()), ctx);
    if (diff.is_ok()) {
        if (is_zero_expr(diff.value(), ctx)) {
            OdeClassification res(OdeType::Exact, E, y, x);
            res.components.push_back(M);
            res.components.push_back(N);
            return res;
        }
    } else {
    }
    
    return std::nullopt;
}

static std::optional<ExprPtr> get_y_power(ExprPtr term, const Symbol& y, symbolic::CASContext& ctx) {
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

[[nodiscard]] static std::optional<OdeClassification> try_bernoulli(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    auto lin = extract_first_order_yprime(E, y, x, y_prime, ctx);
    if (!lin) return std::nullopt;
    ExprPtr A = lin->A;
    ExprPtr B = lin->B;
    
    AstArena& arena = ctx.arena();
    if (depends_on(A, y)) return std::nullopt;
    
    std::vector<ExprPtr> terms;
    std::function<void(ExprPtr, bool)> collect_terms = [&](ExprPtr e, bool neg) {
        if (const auto* sum = expr_cast<Sum>(e)) {
            for (auto t : sum->terms) collect_terms(t, neg);
        } else if (const auto* bin = expr_cast<Binary>(e)) {
            if (bin->op == BinaryOp::Add) {
                collect_terms(bin->left, neg);
                collect_terms(bin->right, neg);
            } else if (bin->op == BinaryOp::Sub) {
                collect_terms(bin->left, neg);
                collect_terms(bin->right, !neg);
            } else {
                terms.push_back(neg ? arena.make<Unary>(UnaryOp::Neg, e) : e);
            }
        } else if (const auto* un = expr_cast<Unary>(e)) {
            if (un->op == UnaryOp::Neg) collect_terms(un->operand, !neg);
            else terms.push_back(neg ? arena.make<Unary>(UnaryOp::Neg, e) : e);
        } else {
            terms.push_back(neg ? arena.make<Unary>(UnaryOp::Neg, e) : e);
        }
    };
    collect_terms(B, false);
    
    struct Group {
        ExprPtr power;
        std::vector<ExprPtr> coeffs;
    };
    std::vector<Group> groups;
    
    for (auto t : terms) {
        auto p_opt = get_y_power(t, y, ctx);
        if (!p_opt) return std::nullopt;
        ExprPtr p = *p_opt;
        
        ExprPtr yp;
        if (is_zero_expr(p, ctx)) {
            yp = arena.make<IntegerLit>(BigInt(1));
        } else if (const auto* il = expr_cast<IntegerLit>(p); il && il->value == BigInt(1)) {
            yp = arena.make<Symbol>(y.name);
        } else {
            yp = arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(y.name), p);
        }
        auto C_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div, t, yp));
        if (C_raw.is_error()) return std::nullopt;
        ExprPtr C = C_raw.value();
        
        bool found = false;
        for (auto& g : groups) {
            auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, g.power, p), ctx);
            if (diff.is_ok() && is_zero_expr(diff.value(), ctx)) {
                g.coeffs.push_back(C);
                found = true;
                break;
            }
        }
        if (!found) {
            groups.push_back({p, {C}});
        }
    }
    
    if (groups.size() != 2) return std::nullopt;
    
    int idx_1 = -1;
    for (size_t i = 0; i < groups.size(); ++i) {
        auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, groups[i].power, arena.make<IntegerLit>(BigInt(1))), ctx);
        if (diff.is_ok() && is_zero_expr(diff.value(), ctx)) {
            idx_1 = static_cast<int>(i);
            break;
        }
    }
    if (idx_1 == -1) return std::nullopt;
    
    int idx_n = 1 - idx_1;
    ExprPtr n_expr = groups[idx_n].power;
    
    ExprPtr P_raw = ctx.simplify(arena.make<Sum>(groups[idx_1].coeffs)).value();
    ExprPtr Q_raw = ctx.simplify(arena.make<Unary>(UnaryOp::Neg, arena.make<Sum>(groups[idx_n].coeffs))).value();
    
    auto P = ctx.simplify(arena.make<Binary>(BinaryOp::Div, P_raw, A));
    auto Q = ctx.simplify(arena.make<Binary>(BinaryOp::Div, Q_raw, A));
    if (P.is_error() || Q.is_error()) return std::nullopt;
    
    OdeClassification res(OdeType::Bernoulli, E, y, x);
    res.components.push_back(P.value());
    res.components.push_back(Q.value());
    res.components.push_back(n_expr);
    return res;
}

[[nodiscard]] static std::optional<OdeClassification> try_homogeneous(
    ExprPtr E, const Symbol& y, const Symbol& x, ExprPtr y_prime, symbolic::CASContext& ctx) {
    auto lin = extract_first_order_yprime(E, y, x, y_prime, ctx);
    if (!lin) return std::nullopt;
    ExprPtr A = lin->A;
    ExprPtr B = lin->B;
    
    AstArena& arena = ctx.arena();
    
    Symbol v_sym = ctx.make_fresh_symbol("v");
    ExprPtr v_expr = arena.make<Symbol>(v_sym.name);
    ExprPtr vx = arena.make<Binary>(BinaryOp::Mul, v_expr, arena.make<Symbol>(x.name));
    
    auto A1_res = substitute_any(A, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(1)), ctx);
    if (A1_res.is_error()) return std::nullopt;
    auto A1_v_res = substitute_any(A1_res.value(), arena.make<Symbol>(y.name), v_expr, ctx);
    if (A1_v_res.is_error()) return std::nullopt;
    ExprPtr A1 = ctx.simplify(A1_v_res.value()).value();
    
    auto B1_res = substitute_any(B, arena.make<Symbol>(x.name), arena.make<IntegerLit>(BigInt(1)), ctx);
    if (B1_res.is_error()) return std::nullopt;
    auto B1_v_res = substitute_any(B1_res.value(), arena.make<Symbol>(y.name), v_expr, ctx);
    if (B1_v_res.is_error()) return std::nullopt;
    ExprPtr B1 = ctx.simplify(B1_v_res.value()).value();
    
    auto Ax_res = substitute_any(A, arena.make<Symbol>(y.name), vx, ctx);
    if (Ax_res.is_error()) return std::nullopt;
    ExprPtr Ax = ctx.simplify(Ax_res.value()).value();
    
    auto Bx_res = substitute_any(B, arena.make<Symbol>(y.name), vx, ctx);
    if (Bx_res.is_error()) return std::nullopt;
    ExprPtr Bx = ctx.simplify(Bx_res.value()).value();
    
    auto term1 = arena.make<Binary>(BinaryOp::Mul, Ax, B1);
    auto term2 = arena.make<Binary>(BinaryOp::Mul, Bx, A1);
    
    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, term1, term2), ctx);
    if (diff.is_ok() && is_zero_expr(diff.value(), ctx)) {
        auto F_v_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Div,
            arena.make<Unary>(UnaryOp::Neg, B1), A1));
        if (F_v_raw.is_ok()) {
            OdeClassification res(OdeType::Homogeneous, E, y, x);
            res.parameter = v_sym;
            res.components.push_back(F_v_raw.value());
            return res;
        }
    }
    
    return std::nullopt;
}


// Detect the Riccati family  y' = q_0(x) + q_1(x)·y + q_2(x)·y²  with q_2 ≢ 0.
// Returns std::nullopt if `E` (the equation written as E = 0) does not fit the
// Riccati pattern; on success, fills `components` with [q_0, q_1, q_2] in this
// order.
//
// Algorithm (no hardcode, no shortcut):
//   1. View E as a polynomial in y of expected degree ≤ 2 with coefficients in
//      Q(x)[y'].  Recover the coefficients c_0(x,y'), c_1(x,y'), c_2(x,y') by
//      evaluating E at four y-values (Lagrange three-point fit + degree witness).
//   2. Reject if any of c_1, c_2 depends on y' (would make the equation
//      non-first-order or non-polynomial in y).
//   3. Reject if c_2 ≡ 0 (would be the already-handled linear branch).
//   4. Split c_0 as α(x)·y' + β(x).  Reject if α ≡ 0 (no first derivative) or
//      α depends on y' (E quadratic in y').
//   5. Normalise q_0 = -β/α, q_1 = -c_1/α, q_2 = -c_2/α.  Reject if any q_i
//      depends on y (genuine non-polynomial coefficient).
[[nodiscard]] static std::optional<OdeClassification> try_riccati(
    ExprPtr E,
    const Symbol& y,
    const Symbol& x,
    ExprPtr y_symbol,
    ExprPtr y_prime,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    auto eval_y = [&](long long k) -> Result<ExprPtr> {
        ExprPtr val = arena.make<IntegerLit>(BigInt(k));
        auto s = substitute_y_shielded(E, y_symbol, val, arena);
        if (s.is_error()) return s;
        return ctx.simplify(s.value());
    };

    auto e0  = eval_y(0);  if (e0.is_error())  return std::nullopt;
    auto e1  = eval_y(1);  if (e1.is_error())  return std::nullopt;
    auto em1 = eval_y(-1); if (em1.is_error()) return std::nullopt;
    auto e2  = eval_y(2);  if (e2.is_error())  return std::nullopt;

    auto two   = arena.make<IntegerLit>(BigInt(2));
    auto four  = arena.make<IntegerLit>(BigInt(4));
    auto half  = arena.make<RationalLit>(BigInt(1), BigInt(2));

    auto c0 = e0.value();
    auto c1_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Mul,
        half,
        arena.make<Binary>(BinaryOp::Sub, e1.value(), em1.value())));
    if (c1_raw.is_error()) return std::nullopt;
    auto c1 = c1_raw.value();

    auto c2_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Mul, half,
            arena.make<Binary>(BinaryOp::Add, e1.value(), em1.value())),
        c0));
    if (c2_raw.is_error()) return std::nullopt;
    auto c2 = c2_raw.value();

    // Degree-witness check: E|_{y=2} must equal 4·c_2 + 2·c_1 + c_0.
    auto expected_e2 = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<Binary>(BinaryOp::Mul, four, c2),
        arena.make<Binary>(BinaryOp::Mul, two,  c1),
        c0}));
    if (expected_e2.is_error()) return std::nullopt;
    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub, e2.value(),
        expected_e2.value()), ctx);
    if (diff.is_error() || !is_zero_expr(diff.value(), ctx)) return std::nullopt;

    // Reject c_2 ≡ 0 (linear-in-y case handled elsewhere).
    if (is_zero_expr(c2, ctx)) return std::nullopt;

    // c_1, c_2 must NOT depend on y' (otherwise structurally non-Riccati).
    auto depends_on_yprime = [&](ExprPtr e) {
        auto mo = find_max_order(e, y, x);
        return mo.has_value() && *mo >= 1U;
    };
    if (depends_on_yprime(c1) || depends_on_yprime(c2)) return std::nullopt;

    // Split c_0 into α(x)·y' + β(x) via two y'-evaluations.
    auto c0_at_yp = [&](long long k) -> Result<ExprPtr> {
        ExprPtr val = arena.make<IntegerLit>(BigInt(k));
        auto s = substitute_any(c0, y_prime, val, ctx);
        if (s.is_error()) return s;
        return ctx.simplify(s.value());
    };
    auto beta_res  = c0_at_yp(0);    if (beta_res.is_error())  return std::nullopt;
    auto c0_at_1   = c0_at_yp(1);    if (c0_at_1.is_error())   return std::nullopt;
    auto c0_at_2   = c0_at_yp(2);    if (c0_at_2.is_error())   return std::nullopt;

    auto alpha_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        c0_at_1.value(), beta_res.value()));
    if (alpha_raw.is_error()) return std::nullopt;
    auto alpha = alpha_raw.value();

    // Affine witness: c_0|_{y'=2} must equal 2·α + β.
    auto expected_c0_at_2 = ctx.simplify(arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, two, alpha),
        beta_res.value()));
    if (expected_c0_at_2.is_error()) return std::nullopt;
    auto c0_diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub,
        c0_at_2.value(), expected_c0_at_2.value()), ctx);
    if (c0_diff.is_error() || !is_zero_expr(c0_diff.value(), ctx)) return std::nullopt;

    if (is_zero_expr(alpha, ctx)) return std::nullopt;          // no y' term
    if (depends_on_yprime(alpha) || depends_on_yprime(beta_res.value()))
        return std::nullopt;                                    // c_0 not affine in y'

    // q_0 = -β/α, q_1 = -c_1/α, q_2 = -c_2/α.
    auto neg = [&](ExprPtr e) { return arena.make<Unary>(UnaryOp::Neg, e); };
    auto divα = [&](ExprPtr num) {
        return ctx.simplify(arena.make<Binary>(BinaryOp::Div, num, alpha));
    };
    auto q0 = divα(neg(beta_res.value()));
    auto q1 = divα(neg(c1));
    auto q2 = divα(neg(c2));
    if (q0.is_error() || q1.is_error() || q2.is_error()) return std::nullopt;

    // Final sanity: q_i must not depend on y (otherwise spurious split).
    auto deps_y = [&](ExprPtr e) { return depends_on(e, y); };
    if (deps_y(q0.value()) || deps_y(q1.value()) || deps_y(q2.value()))
        return std::nullopt;

    OdeClassification res(OdeType::Riccati, /*equation*/E, y, x);
    res.components.push_back(q0.value());
    res.components.push_back(q1.value());
    res.components.push_back(q2.value());
    return res;
}

// Detect Lagrange family:  E ≡ 0  ⇔  y = x·F(p) + G(p),  p ≡ y'.
// Sub-cases:
//   F(p) ≡ p  →  Clairaut, components = [G(p)], parameter = p.
//   F(p) ≢ p  →  d'Alembert, components = [F(p), G(p)], parameter = p.
//
// Algorithm (no hardcode, no shortcut):
//   1. View E as polynomial in y of degree ≤ 1: c_1·y + c_0(x, y').
//      c_1 = E|_{y=1} - E|_{y=0}.  Degree witness E|_{y=2} = 2c_1 + c_0.
//   2. Reject if c_1 depends on x, y, or y' (Lagrange has constant y-coeff).
//   3. Reject if c_1 ≡ 0 (no y term).
//   4. View c_0 as affine in x: c_0 = α(y')·x + β(y').
//      α = c_0|_{x=1} - c_0|_{x=0}.  Witness c_0|_{x=2} = 2α + β.
//   5. Reject if α or β depends on x or y.
//   6. F(p) = -α/c_1, G(p) = -β/c_1, with y' substituted by a fresh p.
//   7. Reject if F or G depends on x or y, or if F ≡ 0 (degenerate).
//   8. Classify: F ≡ p → Clairaut, else → d'Alembert.
[[nodiscard]] static std::optional<OdeClassification> try_lagrange(
    ExprPtr E,
    const Symbol& y,
    const Symbol& x,
    ExprPtr y_symbol,
    ExprPtr y_prime,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    auto eval_y = [&](long long k) -> Result<ExprPtr> {
        ExprPtr val = arena.make<IntegerLit>(BigInt(k));
        auto s = substitute_y_shielded(E, y_symbol, val, arena);
        if (s.is_error()) return s;
        return ctx.simplify(s.value());
    };
    auto e0 = eval_y(0); if (e0.is_error()) return std::nullopt;
    auto e1 = eval_y(1); if (e1.is_error()) return std::nullopt;
    auto e2 = eval_y(2); if (e2.is_error()) return std::nullopt;

    ExprPtr c0 = e0.value();
    auto c1_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        e1.value(), e0.value()));
    if (c1_raw.is_error()) return std::nullopt;
    ExprPtr c1 = c1_raw.value();

    auto two = arena.make<IntegerLit>(BigInt(2));
    auto expected_e2 = ctx.simplify(arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, two, c1), c0));
    if (expected_e2.is_error()) return std::nullopt;
    auto diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub,
        e2.value(), expected_e2.value()), ctx);
    if (diff.is_error() || !is_zero_expr(diff.value(), ctx)) return std::nullopt;

    if (is_zero_expr(c1, ctx)) return std::nullopt;

    auto depends_on_yprime = [&](ExprPtr e) {
        auto mo = find_max_order(e, y, x);
        return mo.has_value() && *mo >= 1U;
    };
    if (depends_on(c1, x) || depends_on(c1, y) || depends_on_yprime(c1))
        return std::nullopt;

    // c_0 = α(y')·x + β(y') affine-in-x extraction.
    auto c0_at_x = [&](long long k) -> Result<ExprPtr> {
        ExprPtr val = arena.make<IntegerLit>(BigInt(k));
        auto s = substitute_any(c0, arena.make<Symbol>(x.name), val, ctx);
        if (s.is_error()) return s;
        return ctx.simplify(s.value());
    };
    auto beta_raw = c0_at_x(0); if (beta_raw.is_error()) return std::nullopt;
    auto c0_at_1  = c0_at_x(1); if (c0_at_1.is_error())  return std::nullopt;
    auto c0_at_2  = c0_at_x(2); if (c0_at_2.is_error())  return std::nullopt;
    ExprPtr beta = beta_raw.value();

    auto alpha_raw = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        c0_at_1.value(), beta));
    if (alpha_raw.is_error()) return std::nullopt;
    ExprPtr alpha = alpha_raw.value();

    auto expected_c0_at_2 = ctx.simplify(arena.make<Binary>(BinaryOp::Add,
        arena.make<Binary>(BinaryOp::Mul, two, alpha), beta));
    if (expected_c0_at_2.is_error()) return std::nullopt;
    auto c0_diff = algebra::expand(arena.make<Binary>(BinaryOp::Sub,
        c0_at_2.value(), expected_c0_at_2.value()), ctx);
    if (c0_diff.is_error() || !is_zero_expr(c0_diff.value(), ctx))
        return std::nullopt;

    if (depends_on(alpha, x) || depends_on(alpha, y)) return std::nullopt;
    if (depends_on(beta,  x) || depends_on(beta,  y)) return std::nullopt;

    // F(p) = -α/c_1, G(p) = -β/c_1 with y' replaced by fresh p.
    Symbol p_sym = ctx.make_fresh_symbol("p");
    ExprPtr p_expr = arena.make<Symbol>(p_sym.name);

    auto neg_div = [&](ExprPtr num) -> Result<ExprPtr> {
        ExprPtr q = arena.make<Binary>(BinaryOp::Div,
            arena.make<Unary>(UnaryOp::Neg, num), c1);
        return ctx.simplify(q);
    };
    auto F_in_yp = neg_div(alpha); if (F_in_yp.is_error()) return std::nullopt;
    auto G_in_yp = neg_div(beta);  if (G_in_yp.is_error()) return std::nullopt;

    auto F_in_p = substitute_any(F_in_yp.value(), y_prime, p_expr, ctx);
    if (F_in_p.is_error()) return std::nullopt;
    auto G_in_p = substitute_any(G_in_yp.value(), y_prime, p_expr, ctx);
    if (G_in_p.is_error()) return std::nullopt;
    auto F_p = ctx.simplify(F_in_p.value());
    auto G_p = ctx.simplify(G_in_p.value());
    if (F_p.is_error() || G_p.is_error()) return std::nullopt;

    if (depends_on(F_p.value(), x) || depends_on(F_p.value(), y))
        return std::nullopt;
    if (depends_on(G_p.value(), x) || depends_on(G_p.value(), y))
        return std::nullopt;
    if (is_zero_expr(F_p.value(), ctx)) return std::nullopt;  // degenerate

    auto f_minus_p = algebra::expand(arena.make<Binary>(BinaryOp::Sub,
        F_p.value(), p_expr), ctx);
    bool is_clairaut = f_minus_p.is_ok() && is_zero_expr(f_minus_p.value(), ctx);

    if (is_clairaut) {
        OdeClassification res(OdeType::Clairaut, /*equation*/E, y, x);
        res.parameter = p_sym;
        res.components.push_back(G_p.value());
        return res;
    }
    OdeClassification res(OdeType::DAlembert, /*equation*/E, y, x);
    res.parameter = p_sym;
    res.components.push_back(F_p.value());
    res.components.push_back(G_p.value());
    return res;
}

[[nodiscard]] Result<OdeClassification> classify_ode(
    ExprPtr equation,
    const Symbol& y,
    const Symbol& x,
    symbolic::CASContext& ctx) {
    
    AstArena& arena = ctx.arena();
    ExprPtr eq_lhs = equation;
    if (const auto* bin = expr_cast<Binary>(equation)) {
        if (bin->op == BinaryOp::Equal) {
            auto sub = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, bin->left, bin->right));
            if (sub.is_error()) return fail<OdeClassification>(sub.error());
            eq_lhs = sub.value();
        }
    }

    auto E_res = algebra::expand(eq_lhs, ctx);
    if (E_res.is_error()) return fail<OdeClassification>(E_res.error());
    ExprPtr E = E_res.value();

    auto max_order_opt = find_max_order(E, y, x);
    if (!max_order_opt) return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    
    uint32_t n = *max_order_opt;
    
    std::vector<ExprPtr> y_ders;
    y_ders.push_back(arena.make<Symbol>(y.name));
    for (uint32_t i = 1; i <= n; ++i) {
        y_ders.push_back(arena.make<Derivative>(y_ders[0], Symbol(x.name), i));
    }
    
    auto get_val = [&](const std::vector<long long>& vals) -> Result<ExprPtr> {
        std::vector<ExprPtr> sub_vals;
        for (long long v : vals) {
            sub_vals.push_back(algebra::poly_make_integer(arena, v));
        }
        return substitute_y_derivatives(E, y_ders, sub_vals, ctx);
    };
    
    std::vector<long long> zero_vals(n + 1, 0);
    auto v0_res = get_val(zero_vals);
    if (v0_res.is_error()) return fail<OdeClassification>(v0_res.error());
    ExprPtr v0 = v0_res.value();
    
    std::vector<ExprPtr> coeffs(n + 1);
    for (uint32_t i = 0; i <= n; ++i) {
        std::vector<long long> one_vals(n + 1, 0);
        one_vals[i] = 1;
        auto vi_res = get_val(one_vals);
        if (vi_res.is_error()) return fail<OdeClassification>(vi_res.error());
        
        auto a_i = ctx.simplify(arena.make<Binary>(BinaryOp::Sub, vi_res.value(), v0));
        if (a_i.is_error()) return fail<OdeClassification>(a_i.error());
        coeffs[i] = a_i.value();
    }
    
    std::vector<ExprPtr> L_terms;
    for (uint32_t i = 0; i <= n; ++i) {
        L_terms.push_back(arena.make<Binary>(BinaryOp::Mul, coeffs[i], y_ders[i]));
    }
    L_terms.push_back(v0);
    ExprPtr L = arena.make<Sum>(std::move(L_terms));
    
    auto diff_res = algebra::expand(arena.make<Binary>(BinaryOp::Sub, E, L), ctx);
    if (diff_res.is_error() || !is_zero_expr(diff_res.value(), ctx)) {
        // The linear template did not match.  For order-1 ODEs we still have
        // a chance to recognise nonlinear families (F5.3): Riccati first,
        // then Lagrange (Clairaut / d'Alembert).
        if (n == 1U) {
            auto separable = try_separable(E, y, x, y_ders[1], ctx);
            if (separable.has_value()) return ok(std::move(*separable));

            auto exact = try_exact(E, y, x, y_ders[1], ctx);
            if (exact.has_value()) return ok(std::move(*exact));

            auto bernoulli = try_bernoulli(E, y, x, y_ders[1], ctx);
            if (bernoulli.has_value()) return ok(std::move(*bernoulli));

            auto homogeneous = try_homogeneous(E, y, x, y_ders[1], ctx);
            if (homogeneous.has_value()) return ok(std::move(*homogeneous));

            auto riccati = try_riccati(E, y, x, y_ders[0], y_ders[1], ctx);
            if (riccati.has_value()) return ok(std::move(*riccati));
            auto lagrange = try_lagrange(E, y, x, y_ders[0], y_ders[1], ctx);
            if (lagrange.has_value()) return ok(std::move(*lagrange));
        }
        return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    }
    
    auto deps_y = [&](ExprPtr e) { return depends_on(e, y); };
    if (deps_y(v0)) return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    for (const auto& a_i : coeffs) {
        if (deps_y(a_i)) return ok(OdeClassification(OdeType::Unknown, equation, y, x));
    }
    
    bool constant_coeffs = true;
    for (const auto& a_i : coeffs) {
        if (depends_on(a_i, x)) {
            constant_coeffs = false;
            break;
        }
    }
    
    if (n >= 1 && !is_zero_expr(coeffs[n], ctx)) {
        if (constant_coeffs) {
            OdeType t = OdeType::LinearNthOrderConstantCoeff;
            if (n == 1) t = OdeType::Linear1stOrder;
            else if (n == 2) t = OdeType::Linear2ndOrderConstantCoeff;
            
            OdeClassification res(t, equation, y, x);
            // Salva come a_n, a_{n-1}, ..., a_0, e infine f(x) = -v0
            for (int i = static_cast<int>(n); i >= 0; --i) {
                res.components.push_back(coeffs[i]);
            }
            res.components.push_back(ctx.simplify(arena.make<Unary>(UnaryOp::Neg, v0)).value());
            return ok(res);
        } else {
             if (n == 2) {
                 OdeClassification res(OdeType::Linear2ndOrderRationalCoeff, equation, y, x);
                 for (int i = 2; i >= 0; --i) {
                     res.components.push_back(coeffs[i]);
                 }
                 res.components.push_back(ctx.simplify(arena.make<Unary>(UnaryOp::Neg, v0)).value());
                 return ok(res);
             }
        }
    }
    
    return ok(OdeClassification(OdeType::Unknown, equation, y, x));
}

[[nodiscard]] Result<ExprPtr> solve_ode(ExprPtr equation, const Symbol& y, const Symbol& x, symbolic::CASContext& ctx) {
    auto class_res = classify_ode(equation, y, x, ctx);
    if (class_res.is_error()) return fail<ExprPtr>(class_res.error());
    
    const auto& classification = class_res.value();
    
    switch (classification.type) {
        case OdeType::Separable:
        case OdeType::Linear1stOrder:
        case OdeType::Bernoulli:
        case OdeType::Exact:
            return solve_ode_1st_order(classification, ctx);

        case OdeType::Riccati:
        case OdeType::Clairaut:
        case OdeType::DAlembert:
            return solve_ode_nonlinear(classification, ctx);

        case OdeType::Linear2ndOrderConstantCoeff:
        case OdeType::Linear2ndOrderRationalCoeff:
        case OdeType::LinearNthOrderConstantCoeff:
            return solve_ode_advanced(classification, ctx);

        default:
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Tipo di ODE non riconosciuto o non supportato analiticamente."));
    }
}

} // namespace cas::calculus
