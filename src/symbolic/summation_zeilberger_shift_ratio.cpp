// summation_zeilberger_shift_ratio.cpp
// Pochhammer-based F(sym+delta)/F(sym) shift-ratio, split out of
// summation_zeilberger_helpers.cpp (T-048 anti-monolith).

#include "summation_zeilberger_helpers.hpp"
#include "summation_zeilberger_internal.hpp"
#include "cas/algebra.hpp"  // for algebra::together
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/symbolic.hpp"
#include <cstdio>
#include <optional>
#include <utility>
#include <vector>

namespace cas::symbolic::zeilberger_detail {

namespace {

// Computes the linear shift Δ = arg(sym+1) - arg(sym) when arg is linear in
// sym (returns std::nullopt for non-linear args).  Linear forms of arg cover
// the binomial / hypergeometric case Γ(a·sym + b) with a ∈ Z.
[[nodiscard]] std::optional<long long> integer_linear_shift(
    ExprPtr arg, const Symbol& sym, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto probe = [&](long long c) -> std::optional<ExprPtr> {
        ExprPtr cv = arena.make<IntegerLit>(BigInt(c));
        ExprPtr cv1 = arena.make<IntegerLit>(BigInt(c + 1));
        auto a1 = ctx.substitute(arg, sym, cv1);
        auto a0 = ctx.substitute(arg, sym, cv);
        if (a1.is_error() || a0.is_error()) return std::nullopt;
        
        ExprPtr diff = arena.make<Binary>(BinaryOp::Sub, a1.value(), a0.value());
        auto tog = algebra::together(diff, ctx);
        if (tog.is_error()) {
            return std::nullopt;
        }
        auto exp = algebra::expand(tog.value(), ctx);
        if (exp.is_error()) return std::nullopt;
        auto simp = ctx.simplify(exp.value());
        if (simp.is_error()) return std::nullopt;
        return simp.value();
    };
    auto to_int = [](ExprPtr e) -> std::optional<long long> {
        if (const auto* il = expr_cast<IntegerLit>(e)) {
            if (il->value.bit_length() > 31U) return std::nullopt;
            long long mag = static_cast<long long>(il->value.abs().to_u64());
            return il->value.is_negative() ? -mag : mag;
        }
        if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
            if (const auto* il = expr_cast<IntegerLit>(un->operand);
                il && il->value.bit_length() <= 31U) {
                return -static_cast<long long>(il->value.abs().to_u64());
            }
        }
        return std::nullopt;
    };
    // Two probes confirm linearity (linear arg has constant first-difference).
    auto d1_opt = probe(7);  if (!d1_opt) return std::nullopt;
    auto d2_opt = probe(31); if (!d2_opt) return std::nullopt;
    auto d1 = to_int(*d1_opt);  if (!d1) return std::nullopt;
    auto d2 = to_int(*d2_opt);  if (!d2) return std::nullopt;
    if (*d1 != *d2) return std::nullopt;  // arg is not linear in sym
    return *d1;
}

// For a Gamma leaf Γ(arg) with integer linear shift Δ, compute the explicit
// rational expression Γ(arg+Δ)/Γ(arg) = ∏_{i=0..Δ-1} (arg+i) for Δ>0 (the
// Pochhammer falling/rising factorial).  Returns 1 for Δ=0, reciprocal form
// for Δ<0.  Leaves arg in symbolic form — no simplification of arg internals.
[[nodiscard]] ExprPtr pochhammer_from_shift(
    ExprPtr arg, long long delta, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (delta == 0) return arena.make<IntegerLit>(BigInt(1));
    long long mag = std::abs(delta);
    std::vector<ExprPtr> factors;
    factors.reserve(static_cast<std::size_t>(mag));
    if (delta > 0) {
        // Γ(arg+Δ)/Γ(arg) = arg·(arg+1)···(arg+Δ-1)
        for (long long i = 0; i < mag; ++i) {
            ExprPtr fi = (i == 0)
                ? arg
                : arena.make<Binary>(BinaryOp::Add, arg,
                    arena.make<IntegerLit>(BigInt(i)));
            factors.push_back(fi);
        }
    } else {
        // Γ(arg-|Δ|)/Γ(arg) = 1/((arg-1)·(arg-2)···(arg-|Δ|))
        for (long long i = 1; i <= mag; ++i) {
            ExprPtr fi = arena.make<Binary>(BinaryOp::Sub, arg,
                arena.make<IntegerLit>(BigInt(i)));
            factors.push_back(fi);
        }
    }
    ExprPtr prod;
    if (factors.size() == 1U) prod = factors.front();
    else prod = arena.make<Product>(std::move(factors));
    if (delta > 0) return prod;
    return arena.make<Binary>(BinaryOp::Div,
        arena.make<IntegerLit>(BigInt(1)), prod);
}

// Leaf-level Pochhammer decomposition of F(sym+1)/F(sym).  Decomposes F into
// multiplicative leaves, then per leaf:
//   - Γ(arg) with linear-in-sym arg of integer slope Δ →  Pochhammer form
//     (Γ(arg+Δ)/Γ(arg) explicit rational, never goes through simplifier).
//   - non-Γ leaf: substitute sym → sym+1 and form raw ratio.
//
// Avoids the simplifier's known limitation of not cancelling gamma factors
// inside Sum aggregates post-`together` (HARDCODE_LEDGER F5.7-ZEIL-GAMMA-RATIO).
// Returns std::nullopt if any Γ leaf has non-integer-linear shift in sym
// (caller falls back to the general expand+cancel path).
[[nodiscard]] std::optional<std::pair<ExprPtr, ExprPtr>> compute_shift_ratio_via_pochhammer(
    ExprPtr F, const Symbol& sym, symbolic::CASContext& ctx, long long total_delta) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> num_leaves, den_leaves;
    collect_product_factors(F, num_leaves, den_leaves);

    ExprPtr sym_delta = arena.make<Binary>(BinaryOp::Add,
        arena.make<Symbol>(sym), arena.make<IntegerLit>(BigInt(total_delta)));

    std::vector<ExprPtr> ratio_num, ratio_den;
    auto process_leaf = [&](ExprPtr L, bool in_num) -> bool {
        if (const auto* fc = expr_cast<FuncCall>(L);
            fc && fc->func_id == BuiltinOp::Gamma && fc->args.size() == 1U) {
            auto slope = integer_linear_shift(fc->args[0], sym, ctx);
            if (!slope) return false;
            long long delta = (*slope) * total_delta;
            ExprPtr poch = pochhammer_from_shift(fc->args[0], delta, ctx);
            // poch is Γ(arg+Δ)/Γ(arg). If Δ < 0, it's 1/Prod.
            if (delta >= 0) {
                if (in_num) ratio_num.push_back(poch);
                else        ratio_den.push_back(poch);
            } else {
                // poch = 1/Prod.
                if (const auto* bin = expr_cast<Binary>(poch)) {
                    if (in_num) ratio_den.push_back(bin->right);
                    else        ratio_num.push_back(bin->right);
                } else {
                    if (in_num) ratio_num.push_back(poch);
                    else        ratio_den.push_back(poch);
                }
            }
            return true;
        }
        // Non-Γ leaf: substitute sym→sym+delta, ratio_factor = L'/L.
        auto Ls_res = ctx.substitute(L, sym, sym_delta);
        if (Ls_res.is_error()) return false;
        if (in_num) {
            ratio_num.push_back(Ls_res.value());
            ratio_den.push_back(L);
        } else {
            ratio_num.push_back(L);
            ratio_den.push_back(Ls_res.value());
        }
        return true;
    };
    for (ExprPtr L : num_leaves) if (!process_leaf(L, true))  return std::nullopt;
    for (ExprPtr L : den_leaves) if (!process_leaf(L, false)) return std::nullopt;

    auto build = [&](std::vector<ExprPtr>& fs) -> ExprPtr {
        if (fs.empty()) return arena.make<IntegerLit>(BigInt(1));
        if (fs.size() == 1U) return fs[0];
        return arena.make<Product>(std::move(fs));
    };
    
    ExprPtr N = build(ratio_num);
    ExprPtr D = build(ratio_den);
    // Structural cancellation.
    ExprPtr ratio = arena.make<Binary>(BinaryOp::Div, N, D);
    ExprPtr cancelled = cancel_common_factors_in_ratio(ratio, ctx);
    
    // Extract numerator and denominator structurally to avoid simplifier distributing Div over Sum.
    std::vector<ExprPtr> nf, df;
    std::function<void(ExprPtr, bool)> collect = [&](ExprPtr f, bool in_den) {
        if (!f) return;
        if (const auto* prod = expr_cast<Product>(f)) {
            for (ExprPtr child : prod->factors) collect(child, in_den);
        } else if (const auto* bin = expr_cast<Binary>(f)) {
            if (bin->op == BinaryOp::Mul) {
                collect(bin->left, in_den);
                collect(bin->right, in_den);
            } else if (bin->op == BinaryOp::Div) {
                collect(bin->left, in_den);
                collect(bin->right, !in_den);
            } else if (bin->op == BinaryOp::Pow) {
                if (const auto* il = expr_cast<IntegerLit>(bin->right)) {
                    long long p = static_cast<long long>(il->value.abs().to_u64());
                    bool neg = il->value.is_negative();
                    for (long long i = 0; i < p; ++i) collect(bin->left, in_den ^ neg);
                } else {
                    (in_den ? df : nf).push_back(f);
                }
            } else {
                (in_den ? df : nf).push_back(f);
            }
        } else if (const auto* un = expr_cast<Unary>(f)) {
            if (un->op == UnaryOp::Neg) {
                nf.push_back(arena.make<IntegerLit>(BigInt(-1)));
                collect(un->operand, in_den);
            } else {
                (in_den ? df : nf).push_back(f);
            }
        } else {
            (in_den ? df : nf).push_back(f);
        }
    };
    collect(cancelled, false);
    
    ExprPtr final_N = build(nf);
    ExprPtr final_D = build(df);
    auto s_num = ctx.simplify(final_N);
    auto s_den = ctx.simplify(final_D);
    return std::make_pair(s_num.is_ok() ? s_num.value() : final_N,
                          s_den.is_ok() ? s_den.value() : final_D);
}

}  // namespace

std::optional<std::pair<ExprPtr, ExprPtr>> compute_shift_ratio(
    ExprPtr F, const Symbol& sym, symbolic::CASContext& ctx, long long delta) {
    if (delta == 0) return std::make_pair(
        ctx.arena().make<IntegerLit>(BigInt(1)),
        ctx.arena().make<IntegerLit>(BigInt(1))
    );
    // Fast-path: leaf-level Pochhammer decomposition.
    if (auto poch = compute_shift_ratio_via_pochhammer(F, sym, ctx, delta);
        poch.has_value()) {
        return poch;
    }
    AstArena& arena = ctx.arena();
    ExprPtr sym_delta = arena.make<Binary>(BinaryOp::Add,
        arena.make<Symbol>(sym),
        arena.make<IntegerLit>(BigInt(delta)));
    auto F_shifted = ctx.substitute(F, sym, sym_delta);
    if (F_shifted.is_error()) return std::nullopt;
    
    // Fallback path.
    ExprPtr F_exp  = expand_gamma_int_shifts(F, ctx);
    ExprPtr Fs_exp = expand_gamma_int_shifts(F_shifted.value(), ctx);
    ExprPtr ratio_raw = arena.make<Binary>(BinaryOp::Div, Fs_exp, F_exp);
    ExprPtr cancelled = cancel_common_factors_in_ratio(ratio_raw, ctx);
    
    std::vector<ExprPtr> nf, df;
    std::function<void(ExprPtr, bool)> collect = [&](ExprPtr f, bool in_den) {
        if (!f) return;
        if (const auto* prod = expr_cast<Product>(f)) {
            for (ExprPtr child : prod->factors) collect(child, in_den);
        } else if (const auto* bin = expr_cast<Binary>(f)) {
            if (bin->op == BinaryOp::Mul) {
                collect(bin->left, in_den);
                collect(bin->right, in_den);
            } else if (bin->op == BinaryOp::Div) {
                collect(bin->left, in_den);
                collect(bin->right, !in_den);
            } else if (bin->op == BinaryOp::Pow) {
                if (const auto* il = expr_cast<IntegerLit>(bin->right)) {
                    long long p = static_cast<long long>(il->value.abs().to_u64());
                    bool neg = il->value.is_negative();
                    for (long long i = 0; i < p; ++i) collect(bin->left, in_den ^ neg);
                } else {
                    (in_den ? df : nf).push_back(f);
                }
            } else {
                (in_den ? df : nf).push_back(f);
            }
        } else if (const auto* un = expr_cast<Unary>(f)) {
            if (un->op == UnaryOp::Neg) {
                nf.push_back(arena.make<IntegerLit>(BigInt(-1)));
                collect(un->operand, in_den);
            } else {
                (in_den ? df : nf).push_back(f);
            }
        } else {
            (in_den ? df : nf).push_back(f);
        }
    };
    collect(cancelled, false);
    
    auto build = [&](std::vector<ExprPtr>& fs) -> ExprPtr {
        if (fs.empty()) return arena.make<IntegerLit>(BigInt(1));
        if (fs.size() == 1U) return fs[0];
        return arena.make<Product>(std::move(fs));
    };
    ExprPtr final_N = build(nf);
    ExprPtr final_D = build(df);
    auto s_num = ctx.simplify(final_N);
    auto s_den = ctx.simplify(final_D);
    return std::make_pair(s_num.is_ok() ? s_num.value() : final_N,
                          s_den.is_ok() ? s_den.value() : final_D);
}

}  // namespace cas::symbolic::zeilberger_detail
