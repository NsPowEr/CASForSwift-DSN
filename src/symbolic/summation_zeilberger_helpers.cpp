// F5.7 — Zeilberger helpers: gamma-shift expansion and structural cancellation.
//
// The CAS simplifier intentionally defers Γ(z+n) reduction inside Product
// contexts so reflection identities (Γ(z)·Γ(1−z) = π/sin(πz)) can fire.  For
// Zeilberger's algorithm we need the opposite behaviour: ratios like
// Γ(k+2)/Γ(k+1) must reduce to (k+1) so the hypergeometric ratio r(k) emerges
// as a plain rational function.  This file implements that reduction
// explicitly without touching the global simplifier policy.

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

// Extract (base, n) from arg = base + n with n ∈ ℤ \ {0}.
std::optional<std::pair<ExprPtr, long long>> extract_integer_shift(
    ExprPtr arg, AstArena& arena) {
    auto try_get_int = [](ExprPtr e, long long& out) -> bool {
        if (const auto* il = expr_cast<IntegerLit>(e)) {
            if (il->value.bit_length() > 31U) return false;
            long long mag = static_cast<long long>(il->value.abs().to_u64());
            out = il->value.is_negative() ? -mag : mag;
            return true;
        }
        if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
            if (const auto* il2 = expr_cast<IntegerLit>(un->operand)) {
                if (il2->value.bit_length() > 31U) return false;
                out = -static_cast<long long>(il2->value.abs().to_u64());
                return true;
            }
        }
        return false;
    };
    auto build_residual = [&arena](std::vector<ExprPtr> rest) -> ExprPtr {
        if (rest.empty()) return arena.make<IntegerLit>(BigInt(0));
        if (rest.size() == 1U) return rest[0];
        return arena.make<Sum>(std::move(rest));
    };
    if (const auto* sum = expr_cast<Sum>(arg)) {
        long long n_total = 0;
        std::vector<ExprPtr> rest;
        for (ExprPtr t : sum->terms) {
            long long v;
            if (try_get_int(t, v)) n_total += v;
            else rest.push_back(t);
        }
        if (n_total != 0)
            return std::make_pair(build_residual(std::move(rest)), n_total);
    }
    if (const auto* bin = expr_cast<Binary>(arg)) {
        if (bin->op == BinaryOp::Add) {
            long long v;
            if (try_get_int(bin->right, v)) return std::make_pair(bin->left, v);
            if (try_get_int(bin->left,  v)) return std::make_pair(bin->right, v);
        }
        if (bin->op == BinaryOp::Sub) {
            long long v;
            if (try_get_int(bin->right, v)) return std::make_pair(bin->left, -v);
        }
    }
    return std::nullopt;
}


}  // namespace

bool struct_equal(ExprPtr a, ExprPtr b) {
    if (a == b) return true;
    if (!a || !b) return false;
    if (const auto* sa = expr_cast<Symbol>(a))
        if (const auto* sb = expr_cast<Symbol>(b)) return sa->name == sb->name;
    if (const auto* ia = expr_cast<IntegerLit>(a))
        if (const auto* ib = expr_cast<IntegerLit>(b)) return ia->value == ib->value;
    if (const auto* ra = expr_cast<RationalLit>(a))
        if (const auto* rb = expr_cast<RationalLit>(b))
            return ra->numerator == rb->numerator &&
                   ra->denominator == rb->denominator;
    if (const auto* ca = expr_cast<Constant>(a))
        if (const auto* cb = expr_cast<Constant>(b)) return ca->value == cb->value;
    if (const auto* fa = expr_cast<FuncCall>(a))
        if (const auto* fb = expr_cast<FuncCall>(b)) {
            if (fa->func_id != fb->func_id || fa->args.size() != fb->args.size())
                return false;
            for (std::size_t i = 0; i < fa->args.size(); ++i)
                if (!struct_equal(fa->args[i], fb->args[i])) return false;
            return true;
        }
    if (const auto* ba = expr_cast<Binary>(a))
        if (const auto* bb = expr_cast<Binary>(b))
            return ba->op == bb->op &&
                   struct_equal(ba->left,  bb->left) &&
                   struct_equal(ba->right, bb->right);
    if (const auto* ua = expr_cast<Unary>(a))
        if (const auto* ub = expr_cast<Unary>(b))
            return ua->op == ub->op && struct_equal(ua->operand, ub->operand);
    if (const auto* su_a = expr_cast<Sum>(a))
        if (const auto* su_b = expr_cast<Sum>(b)) {
            if (su_a->terms.size() != su_b->terms.size()) return false;
            for (std::size_t i = 0; i < su_a->terms.size(); ++i)
                if (!struct_equal(su_a->terms[i], su_b->terms[i])) return false;
            return true;
        }
    if (const auto* pa = expr_cast<Product>(a))
        if (const auto* pb = expr_cast<Product>(b)) {
            if (pa->factors.size() != pb->factors.size()) return false;
            for (std::size_t i = 0; i < pa->factors.size(); ++i)
                if (!struct_equal(pa->factors[i], pb->factors[i])) return false;
            return true;
        }
    return false;
}

ExprPtr expand_gamma_int_shifts(ExprPtr e, symbolic::CASContext& ctx) {
    if (!e) return e;
    AstArena& arena = ctx.arena();
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        std::vector<ExprPtr> new_args;
        new_args.reserve(fc->args.size());
        for (ExprPtr a : fc->args) new_args.push_back(expand_gamma_int_shifts(a, ctx));
        if (fc->func_id == BuiltinOp::Gamma && new_args.size() == 1U) {
            // Normalise the gamma argument so nested Adds (e.g. (k+1)+1)
            // collapse to a canonical k+m form before integer-shift extraction.
            auto arg_simp = ctx.simplify(new_args[0]);
            if (arg_simp.is_ok()) new_args[0] = arg_simp.value();
            auto shift = extract_integer_shift(new_args[0], arena);
            if (shift) {
                ExprPtr base = shift->first;
                long long n = shift->second;
                // Canonicalise base so gamma bases match structurally across
                // F and F_shifted (otherwise Binary(Sub,n,k) vs Sum([n,-k])
                // mismatch breaks the rationalize-gammas pass).
                auto base_simp = ctx.simplify(base);
                if (base_simp.is_ok()) base = base_simp.value();
                std::vector<ExprPtr> gb_args{base};
                ExprPtr gamma_base =
                    arena.make<FuncCall>(BuiltinOp::Gamma, std::move(gb_args));
                if (n > 0) {
                    ExprPtr result = gamma_base;
                    for (long long i = 0; i < n; ++i) {
                        ExprPtr factor = (i == 0)
                            ? base
                            : arena.make<Binary>(BinaryOp::Add, base,
                                arena.make<IntegerLit>(BigInt(i)));
                        result = arena.make<Binary>(BinaryOp::Mul, result, factor);
                    }
                    return result;
                }
                long long m = -n;
                ExprPtr den = arena.make<IntegerLit>(BigInt(1));
                for (long long i = 1; i <= m; ++i) {
                    ExprPtr factor = arena.make<Binary>(BinaryOp::Sub, base,
                        arena.make<IntegerLit>(BigInt(i)));
                    den = arena.make<Binary>(BinaryOp::Mul, den, factor);
                }
                return arena.make<Binary>(BinaryOp::Div, gamma_base, den);
            }
        }
        return arena.make<FuncCall>(fc->func_id, std::move(new_args));
    }
    if (const auto* bin = expr_cast<Binary>(e))
        return arena.make<Binary>(bin->op,
            expand_gamma_int_shifts(bin->left, ctx),
            expand_gamma_int_shifts(bin->right, ctx));
    if (const auto* un = expr_cast<Unary>(e))
        return arena.make<Unary>(un->op,
            expand_gamma_int_shifts(un->operand, ctx));
    if (const auto* sum = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> nt;
        for (ExprPtr t : sum->terms) nt.push_back(expand_gamma_int_shifts(t, ctx));
        return arena.make<Sum>(std::move(nt));
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        std::vector<ExprPtr> nf;
        for (ExprPtr f : prod->factors) nf.push_back(expand_gamma_int_shifts(f, ctx));
        return arena.make<Product>(std::move(nf));
    }
    return e;
}

ExprPtr cancel_common_factors_in_ratio(
    ExprPtr ratio, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> num_factors, den_factors;
    collect_product_factors(ratio, num_factors, den_factors);

    // Cancel matching factors structurally.
    std::vector<bool> n_used(num_factors.size(), false);
    std::vector<bool> d_used(den_factors.size(), false);
    for (std::size_t i = 0; i < num_factors.size(); ++i) {
        if (n_used[i]) continue;
        for (std::size_t j = 0; j < den_factors.size(); ++j) {
            if (d_used[j]) continue;
            if (struct_equal(num_factors[i], den_factors[j])) {
                n_used[i] = true; d_used[j] = true; break;
            }
        }
    }
    std::vector<ExprPtr> new_n, new_d;
    for (std::size_t i = 0; i < num_factors.size(); ++i)
        if (!n_used[i]) new_n.push_back(num_factors[i]);
    for (std::size_t i = 0; i < den_factors.size(); ++i)
        if (!d_used[i]) new_d.push_back(den_factors[i]);

    auto build = [&arena](std::vector<ExprPtr> fs) -> ExprPtr {
        if (fs.empty()) return arena.make<IntegerLit>(BigInt(1));
        if (fs.size() == 1U) return fs[0];
        return arena.make<Product>(std::move(fs));
    };
    return arena.make<Binary>(BinaryOp::Div,
        build(std::move(new_n)), build(std::move(new_d)));
}

}  // namespace cas::symbolic::zeilberger_detail
