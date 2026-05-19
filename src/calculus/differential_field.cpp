#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"
#include "calculus_internal.hpp"
#include <algorithm>
#include <string>

namespace cas::calculus {

namespace {

// F2.1.b: recursion depth bound. AST descent is strictly one level per
// call; an expression of N nodes has nesting depth ≤ N. Cap at 4096
// matches the simplifier's MAX_SIMPLIFICATION_DEPTH (300) with margin
// for AST trees that have not been simplified yet. Beyond this depth,
// return Unimplemented so the caller can bail diagnostically rather
// than blowing the stack.
constexpr unsigned int kVisitRecursiveMaxDepth = 4096U;

template <typename F>
Result<void> visit_recursive_impl(ExprPtr expr, F&& f, unsigned int depth) {
    if (depth >= kVisitRecursiveMaxDepth) {
        return fail<void>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Differential field visit recursion budget exceeded",
            .hint = std::nullopt,
        });
    }
    auto res = f(expr);
    if (res.is_error()) return res;

    return visit_expr(expr, [&](const auto& node) -> Result<void> {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, Unary>) {
            return visit_recursive_impl(node.operand, f, depth + 1U);
        } else if constexpr (std::is_same_v<T, Binary>) {
            auto r1 = visit_recursive_impl(node.left, f, depth + 1U);
            if (r1.is_error()) return r1;
            return visit_recursive_impl(node.right, f, depth + 1U);
        } else if constexpr (std::is_same_v<T, FuncCall>) {
            for (ExprPtr arg : node.args) {
                auto r = visit_recursive_impl(arg, f, depth + 1U);
                if (r.is_error()) return r;
            }
        } else if constexpr (std::is_same_v<T, Sum>) {
            for (ExprPtr term : node.terms) {
                auto r = visit_recursive_impl(term, f, depth + 1U);
                if (r.is_error()) return r;
            }
        } else if constexpr (std::is_same_v<T, Product>) {
            for (ExprPtr factor : node.factors) {
                auto r = visit_recursive_impl(factor, f, depth + 1U);
                if (r.is_error()) return r;
            }
        }
        return ok();
    });
}

template <typename F>
Result<void> visit_recursive(ExprPtr expr, F&& f) {
    return visit_recursive_impl(expr, std::forward<F>(f), 0U);
}

ExprPtr substitute_pattern(ExprPtr expr, ExprPtr pattern, ExprPtr replacement, AstArena& arena) {
    if (structural_equal(expr, pattern)) return replacement;

    return visit_expr(expr, [&](const auto& node) -> ExprPtr {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, Unary>) {
            return arena.make<Unary>(node.op, substitute_pattern(node.operand, pattern, replacement, arena));
        } else if constexpr (std::is_same_v<T, Binary>) {
            return arena.make<Binary>(node.op, 
                substitute_pattern(node.left, pattern, replacement, arena),
                substitute_pattern(node.right, pattern, replacement, arena));
        } else if constexpr (std::is_same_v<T, FuncCall>) {
            std::vector<ExprPtr> args;
            for (ExprPtr arg : node.args) args.push_back(substitute_pattern(arg, pattern, replacement, arena));
            return arena.make<FuncCall>(node.name, std::move(args));
        } else if constexpr (std::is_same_v<T, Sum>) {
            std::vector<ExprPtr> terms;
            for (ExprPtr term : node.terms) terms.push_back(substitute_pattern(term, pattern, replacement, arena));
            return arena.make<Sum>(std::move(terms));
        } else if constexpr (std::is_same_v<T, Product>) {
            std::vector<ExprPtr> factors;
            for (ExprPtr factor : node.factors) factors.push_back(substitute_pattern(factor, pattern, replacement, arena));
            return arena.make<Product>(std::move(factors));
        }
        return expr;
    });
}

} // namespace

Result<DifferentialField> DifferentialField::build(ExprPtr expr, const Symbol& x, symbolic::CASContext& ctx) {
    DifferentialField field(x);
    auto res = field.add_extension(expr, ctx);
    if (res.is_error()) return fail<DifferentialField>(res.error());
    return ok(field);
}

Result<void> DifferentialField::add_extension(ExprPtr expr, [[maybe_unused]] symbolic::CASContext& ctx) {
    return visit_recursive(expr, [&](ExprPtr node) -> Result<void> {
        if (const auto* call = expr_cast<FuncCall>(node)) {
            ExtensionType type;
            if (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) type = ExtensionType::Logarithmic;
            else if (call->func_id == BuiltinOp::Exp) type = ExtensionType::Exponential;
            else return ok();

            bool exists = false;
            for (const auto& ext : extensions_) {
                if (structural_equal(ext.argument, call->args[0]) && ext.type == type) {
                    exists = true;
                    break;
                }
            }

            if (!exists) {
                Symbol t("t_" + std::to_string(extensions_.size()));
                extensions_.push_back({type, call->args[0], t});
            }
        }
        return ok();
    });
}

Result<ExprPtr> DifferentialField::to_field_generators(ExprPtr expr, symbolic::CASContext& ctx) const {
    ExprPtr current = expr;
    AstArena& arena = ctx.arena();

    for (const auto& ext : extensions_) {
        std::string func_name = (ext.type == ExtensionType::Logarithmic) ? "ln" : "exp";
        ExprPtr pattern = arena.make<FuncCall>(func_name, std::vector<ExprPtr>{ext.argument});
        ExprPtr replacement = arena.make<Symbol>(ext.t_var.name);
        
        current = substitute_pattern(current, pattern, replacement, arena);
    }
    return ok(current);
}

Result<ExprPtr> DifferentialField::from_field_generators(ExprPtr expr, symbolic::CASContext& ctx) const {
    ExprPtr current = expr;
    AstArena& arena = ctx.arena();

    for (auto it = extensions_.rbegin(); it != extensions_.rend(); ++it) {
        std::string func_name = (it->type == ExtensionType::Logarithmic) ? "ln" : "exp";
        ExprPtr pattern = arena.make<Symbol>(it->t_var.name);
        ExprPtr replacement = arena.make<FuncCall>(func_name, std::vector<ExprPtr>{it->argument});
        
        current = substitute_pattern(current, pattern, replacement, arena);
    }
    return ok(current);
}

Result<ExprPtr> DifferentialField::derive(ExprPtr expr, symbolic::CASContext& ctx) const {
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        if (sym->name == base_var_.name) return ok(ctx.arena().make<IntegerLit>(BigInt(1)));
        
        for (const auto& ext : extensions_) {
            if (sym->name == ext.t_var.name) {
                auto du = derive(ext.argument, ctx);
                if (du.is_error()) return du;
                
                if (ext.type == ExtensionType::Logarithmic) {
                    return ok(ctx.arena().make<Binary>(BinaryOp::Div, du.value(), ext.argument));
                } else {
                    return ok(ctx.arena().make<Product>(std::vector<ExprPtr>{
                        ctx.arena().make<Symbol>(ext.t_var.name), du.value()
                    }));
                }
            }
        }
    }
    
    return diff(expr, base_var_, 1U, ctx);
}

Result<ExprPtr> compute_resultant(ExprPtr A, ExprPtr B, const Symbol& x, symbolic::CASContext& ctx) {
    return algebra::polynomial_resultant(A, B, x, ctx);
}

// Horowitz-Ostrogradsky one-step reduction for ∫ P/V^n (n >= 2, V squarefree).
// Finds A of deg < deg(V) such that (n-1)*A*V' ≡ -P (mod V) via extended GCD.
// Returns updated P and adds -A/((n-1)*V^(n-1)) to rational_part.
static bool ho_reduce_step(
    algebra::PolyExpr& P_poly,     // in/out: numerator poly (updated)
    const algebra::PolyExpr& V_poly,
    unsigned int n,
    ExprPtr V,
    const Symbol& x,
    ExprPtr& rational_part,        // accumulates rational contributions
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    // V' as poly
    auto V_expr_res = algebra::polynomial_to_expr(V_poly, x, ctx);
    if (V_expr_res.is_error()) return false;
    auto Vd_res = diff(V_expr_res.value(), x, 1U, ctx);
    if (Vd_res.is_error()) return false;
    auto Vd_poly_res = algebra::parse_polynomial(Vd_res.value(), x, ctx);
    if (Vd_poly_res.is_error()) return false;
    const auto& Vd_poly = Vd_poly_res.value();

    // (n-1)*V' as poly
    algebra::PolyExpr nmo_poly;
    nmo_poly.push_back(algebra::poly_make_integer(arena, static_cast<long long>(n - 1)));
    auto nmo_Vd_res = algebra::poly_multiply(nmo_poly, Vd_poly, ctx);
    if (nmo_Vd_res.is_error()) return false;

    // Extended GCD: s*(n-1)*V' + t*V = g  (g = constant since gcd(V,V')=1 for squarefree V)
    auto xgcd_res = algebra::poly_extended_gcd(nmo_Vd_res.value(), V_poly, ctx);
    if (xgcd_res.is_error()) return false;
    const auto& xgcd = xgcd_res.value();

    // P mod V
    auto P_rem_res = algebra::divide_poly_with_remainder(P_poly, V_poly, ctx);
    if (P_rem_res.is_error()) return false;
    const auto& P_rem = P_rem_res.value().remainder;

    // A_raw = -P_rem * s
    auto neg_P_res = algebra::poly_negate(P_rem, ctx);
    if (neg_P_res.is_error()) return false;
    auto A_raw_res = algebra::poly_multiply(neg_P_res.value(), xgcd.s, ctx);
    if (A_raw_res.is_error()) return false;

    // A_poly = A_raw mod V  (deg A < deg V)
    auto A_rem_res = algebra::divide_poly_with_remainder(A_raw_res.value(), V_poly, ctx);
    if (A_rem_res.is_error()) return false;
    algebra::PolyExpr A_poly = A_rem_res.value().remainder;

    // Divide by gcd (should be a nonzero constant)
    if (!algebra::is_zero_poly(xgcd.gcd) && algebra::poly_degree(xgcd.gcd) == 0) {
        ExprPtr g_val = algebra::leading_coefficient(xgcd.gcd);
        auto A_scaled = algebra::poly_divide_by_scalar(A_poly, g_val, ctx);
        if (A_scaled.is_error()) return false;
        A_poly = A_scaled.value();
    }

    // rational_part += A / V^(n-1)
    // From H-O derivation: d/dx[A/V^(n-1)] = (A'V - (n-1)AV')/V^n, so P/V^n = d/dx[A/V^(n-1)] + B/V^(n-1)
    auto A_expr_res = algebra::polynomial_to_expr(A_poly, x, ctx);
    if (A_expr_res.is_error()) return false;
    ExprPtr A_expr = A_expr_res.value();
    ExprPtr Vn1 = arena.make<Binary>(BinaryOp::Pow, V, arena.make<IntegerLit>(BigInt(n - 1)));
    ExprPtr contrib = arena.make<Binary>(BinaryOp::Div, A_expr, Vn1);
    rational_part = arena.make<Binary>(BinaryOp::Add, rational_part, contrib);

    // Compute A' (polynomial derivative of A_poly)
    algebra::PolyExpr A_prime_poly;
    if (A_poly.size() > 1) {
        A_prime_poly.resize(A_poly.size() - 1, nullptr);
        for (std::size_t k = 1; k < A_poly.size(); ++k) {
            ExprPtr coeff = A_poly[k];
            if (!coeff) continue;
            ExprPtr k_expr = algebra::poly_make_integer(arena, static_cast<long long>(k));
            auto term = algebra::poly_simplify_expr(arena.make<Binary>(BinaryOp::Mul, k_expr, coeff), ctx);
            if (term.is_ok()) A_prime_poly[k - 1] = term.value();
        }
        algebra::normalize_poly(A_prime_poly);
    }

    // New P = (P - A'*V + (n-1)*A*V') / V  (H-O standard formula)
    auto A_Vd_res = algebra::poly_multiply(A_poly, Vd_poly, ctx);
    if (A_Vd_res.is_error()) return false;
    auto nmo_A_Vd_res = algebra::poly_multiply(nmo_poly, A_Vd_res.value(), ctx);
    if (nmo_A_Vd_res.is_error()) return false;
    auto new_P_sum_res = algebra::poly_add(P_poly, nmo_A_Vd_res.value(), ctx);
    if (new_P_sum_res.is_error()) return false;
    algebra::PolyExpr new_P_full = new_P_sum_res.value();
    if (!algebra::is_zero_poly(A_prime_poly)) {
        auto A_prime_V_res = algebra::poly_multiply(A_prime_poly, V_poly, ctx);
        if (A_prime_V_res.is_ok()) {
            auto sub_res = algebra::poly_subtract(new_P_full, A_prime_V_res.value(), ctx);
            if (sub_res.is_ok()) new_P_full = sub_res.value();
        }
    }
    auto new_P_div_res = algebra::divide_poly_with_remainder(new_P_full, V_poly, ctx);
    if (new_P_div_res.is_error()) return false;
    P_poly = new_P_div_res.value().quotient;

    return true;
}

/// Performs Hermite reduction on P/Q with respect to generator t_var
Result<HermiteReduction> hermite_reduce(
    ExprPtr P, ExprPtr Q, const Symbol& t_var, const DifferentialField& field, symbolic::CASContext& ctx) {

    // If we are in the base field (no extensions), use the exact RatPoly-based reduction.
    // This is more robust for rational functions over Q(parameters).
    if (field.extensions().empty() && t_var.name == field.base_var().name) {
        return hermite_reduction_exact(P, Q, t_var, ctx);
    }

    AstArena& arena = ctx.arena();
    const Symbol& x = field.base_var();

    // Square-free factorization of Q
    auto sqf_res = algebra::square_free_factorization(Q, x, ctx);
    if (sqf_res.is_error()) {
        return ok(HermiteReduction{
            .rational_part = arena.make<IntegerLit>(BigInt(0)),
            .remaining_P = P,
            .remaining_Q = Q
        });
    }
    const auto& sqf = sqf_res.value();

    bool has_repeated = false;
    for (const auto& f : sqf.factors) {
        if (f.multiplicity >= 2) { has_repeated = true; break; }
    }
    if (!has_repeated) {
        return ok(HermiteReduction{
            .rational_part = arena.make<IntegerLit>(BigInt(0)),
            .remaining_P = P,
            .remaining_Q = Q
        });
    }

    // Parse current numerator as polynomial
    auto P_poly_res = algebra::parse_polynomial(P, x, ctx);
    if (P_poly_res.is_error()) {
        return ok(HermiteReduction{
            .rational_part = arena.make<IntegerLit>(BigInt(0)),
            .remaining_P = P,
            .remaining_Q = Q
        });
    }
    algebra::PolyExpr P_poly = P_poly_res.value();
    ExprPtr rational_part = arena.make<IntegerLit>(BigInt(0));
    ExprPtr current_Q = Q;

    for (const auto& factor : sqf.factors) {
        if (factor.multiplicity < 2) continue;
        ExprPtr V = factor.factor;
        auto V_poly_res = algebra::parse_polynomial(V, x, ctx);
        if (V_poly_res.is_error()) continue;
        const auto& V_poly = V_poly_res.value();

        for (unsigned int n = factor.multiplicity; n >= 2; --n) {
            if (!ho_reduce_step(P_poly, V_poly, n, V, x, rational_part, ctx)) break;

            // Remove one power of V from the denominator
            auto Q_poly_res = algebra::parse_polynomial(current_Q, x, ctx);
            if (Q_poly_res.is_error()) break;
            auto Q_div = algebra::divide_poly_with_remainder(Q_poly_res.value(), V_poly, ctx);
            if (Q_div.is_error()) break;
            auto new_Q_res = algebra::polynomial_to_expr(Q_div.value().quotient, x, ctx);
            if (new_Q_res.is_error()) break;
            current_Q = new_Q_res.value();
        }
    }

    auto P_expr_res = algebra::polynomial_to_expr(P_poly, x, ctx);
    ExprPtr remaining_P = P_expr_res.is_ok() ? P_expr_res.value() : P;

    {
        auto s = ctx.simplify(rational_part);
        if (s.is_ok()) rational_part = s.value();
    }
    {
        auto s = ctx.simplify(remaining_P);
        if (s.is_ok()) remaining_P = s.value();
    }
    {
        auto s = ctx.simplify(current_Q);
        if (s.is_ok()) current_Q = s.value();
    }

    return ok(HermiteReduction{
        .rational_part = rational_part,
        .remaining_P = remaining_P,
        .remaining_Q = current_Q
    });
}

Result<ExprPtr> integrate_rothstein_trager(
    ExprPtr P, ExprPtr Q, [[maybe_unused]] const Symbol& t_var, const DifferentialField& field, symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();
    const Symbol& x = field.base_var();

    auto dQ_res = diff(Q, x, 1U, ctx);
    if (dQ_res.is_error()) return fail<ExprPtr>(dQ_res.error());
    ExprPtr dQ = dQ_res.value();

    // t must be a fresh variable not appearing in P or Q.
    const Symbol t_fresh = ctx.make_fresh_symbol("rt_t");
    ExprPtr t_sym = arena.make<Symbol>(t_fresh.name);
    ExprPtr t_dQ = arena.make<Product>(std::vector<ExprPtr>{t_sym, dQ});
    ExprPtr A = arena.make<Binary>(BinaryOp::Sub, P, t_dQ);

    auto R_res = compute_resultant(A, Q, x, ctx);
    if (R_res.is_error()) return fail<ExprPtr>(R_res.error());
    ExprPtr R = R_res.value();

    // Check P == 0 early (integral of 0 is 0)
    {
        bool P_zero = false;
        if (const auto* il = expr_cast<IntegerLit>(P)) P_zero = il->value.is_zero();
        if (P_zero) return ok(arena.make<IntegerLit>(BigInt(0)));
    }

    auto roots_res = algebra::solve_polynomial(R, t_fresh, ctx);
    if (roots_res.is_error()) return fail<ExprPtr>(roots_res.error());

    // Deduplicate roots (structural equality)
    std::vector<ExprPtr> unique_roots;
    for (ExprPtr r : roots_res.value()) {
        bool dup = false;
        for (ExprPtr u : unique_roots) { if (structural_equal(r, u)) { dup = true; break; } }
        if (!dup) unique_roots.push_back(r);
    }

    if (unique_roots.empty()) {
        // No rational roots found but P != 0 → integral has non-rational parts (arctan etc.)
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Rothstein-Trager: resultant has no rational roots; integral requires arctan or algebraic extensions",
        });
    }

    std::vector<ExprPtr> integral_terms;
    for (ExprPtr root : unique_roots) {
        ExprPtr root_dQ = arena.make<Product>(std::vector<ExprPtr>{root, dQ});
        ExprPtr A_root = arena.make<Binary>(BinaryOp::Sub, P, root_dQ);
        { auto s = ctx.simplify(A_root); if (s.is_ok()) A_root = s.value(); }

        auto v_res = algebra::polynomial_gcd(A_root, Q, x, ctx);
        if (v_res.is_error()) continue;
        ExprPtr v = v_res.value();

        // Use ln(|v|) for real-domain correctness
        ExprPtr abs_v = arena.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{v});
        ExprPtr ln_v = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{abs_v});
        ExprPtr term = arena.make<Product>(std::vector<ExprPtr>{root, ln_v});
        integral_terms.push_back(term);
    }

    if (integral_terms.empty()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Rothstein-Trager: no valid GCD terms computed (roots may be algebraic/RootOf)",
        });
    } else if (integral_terms.size() == 1) {
        return ok(integral_terms[0]);
    } else {
        return ok(arena.make<Sum>(std::move(integral_terms)));
    }
}

} // namespace cas::calculus
