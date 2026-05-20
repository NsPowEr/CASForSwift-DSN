#include "polynomial_groebner_f4.hpp"

#include "algebra_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include <cstdint>
#include <limits>
#include <optional>

namespace cas::algebra {
namespace {

[[nodiscard]] Result<std::uint64_t> exact_nonnegative_exponent(const BigInt& value) {
    if (value.is_negative()) {
        return fail<std::uint64_t>(make_error(
            CASErrorKind::Unimplemented,
            "F4 polynomial conversion requires non-negative integer exponents"));
    }
    if (value.bit_length() > 63U) {
        return fail<std::uint64_t>(make_error(
            CASErrorKind::Overflow,
            "F4 polynomial exponent exceeds supported explicit expansion bound"));
    }
    return ok(value.to_u64());
}

[[nodiscard]] PolyF4 add_poly(PolyF4 lhs, const PolyF4& rhs) {
    for (const auto& [mon, coeff] : rhs.terms) {
        lhs.terms[mon] = lhs.terms[mon] + coeff;
        if (lhs.terms[mon].numerator().is_zero()) lhs.terms.erase(mon);
    }
    return lhs;
}

[[nodiscard]] PolyF4 subtract_poly(PolyF4 lhs, const PolyF4& rhs) {
    for (const auto& [mon, coeff] : rhs.terms) {
        lhs.terms[mon] = lhs.terms[mon] - coeff;
        if (lhs.terms[mon].numerator().is_zero()) lhs.terms.erase(mon);
    }
    return lhs;
}

[[nodiscard]] Result<PolyF4> multiply_poly(const PolyF4& lhs, const PolyF4& rhs, std::size_t n_vars) {
    PolyF4 res;
    for (const auto& [mon_l, coeff_l] : lhs.terms) {
        for (const auto& [mon_r, coeff_r] : rhs.terms) {
            Monomial mon(n_vars);
            for (std::size_t k = 0; k < n_vars; ++k) {
                const std::uint64_t exp_sum =
                    static_cast<std::uint64_t>(mon_l[k]) + static_cast<std::uint64_t>(mon_r[k]);
                if (exp_sum > std::numeric_limits<unsigned int>::max()) {
                    return fail<PolyF4>(make_error(
                        CASErrorKind::Overflow,
                        "F4 monomial exponent overflow during polynomial multiplication"));
                }
                mon[k] = static_cast<unsigned int>(exp_sum);
            }
            res.terms[mon] = res.terms[mon] + coeff_l * coeff_r;
            if (res.terms[mon].numerator().is_zero()) res.terms.erase(mon);
        }
    }
    return ok(std::move(res));
}

[[nodiscard]] Result<PolyF4> pow_poly(PolyF4 base, std::uint64_t exponent, std::size_t n_vars) {
    PolyF4 res;
    res.terms[Monomial(n_vars, 0)] = Rational(1);
    while (exponent > 0U) {
        if ((exponent % 2U) == 1U) {
            auto next = multiply_poly(res, base, n_vars);
            if (next.is_error()) return next;
            res = std::move(next.value());
        }
        exponent /= 2U;
        if (exponent > 0U) {
            auto next_base = multiply_poly(base, base, n_vars);
            if (next_base.is_error()) return next_base;
            base = std::move(next_base.value());
        }
    }
    return ok(std::move(res));
}

} // namespace

Result<PolyF4> expr_to_f4(ExprPtr expr, const std::vector<Symbol>& vars, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<PolyF4>(make_error(CASErrorKind::InvalidArgument, "Null expression in expr_to_f4"));
    }

    const std::size_t n = vars.size();
    auto get_var_idx = [&](const std::string& name) -> std::optional<std::size_t> {
        for (std::size_t i = 0; i < n; ++i) {
            if (vars[i].name == name) return i;
        }
        return std::nullopt;
    };

    if (const auto* il = expr_cast<IntegerLit>(expr)) {
        PolyF4 p;
        if (!il->value.is_zero()) p.terms[Monomial(n, 0)] = Rational(il->value);
        return ok(std::move(p));
    }
    if (const auto* rl = expr_cast<RationalLit>(expr)) {
        PolyF4 p;
        if (!rl->numerator.is_zero()) p.terms[Monomial(n, 0)] = Rational(rl->numerator, rl->denominator);
        return ok(std::move(p));
    }
    if (expr_cast<DecimalLit>(expr)) {
        return fail<PolyF4>(make_error(
            CASErrorKind::Unimplemented,
            "DecimalLit is not accepted in exact F4 polynomial conversion"));
    }
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        auto idx = get_var_idx(sym->name);
        if (idx) {
            PolyF4 p;
            Monomial mon(n, 0);
            mon[*idx] = 1;
            p.terms[mon] = Rational(1);
            return ok(std::move(p));
        }
        return fail<PolyF4>(make_error(
            CASErrorKind::Unimplemented,
            "Parametric constants not supported in expr_to_f4: " + sym->name));
    }
    if (const auto* un = expr_cast<Unary>(expr)) {
        if (un->op == UnaryOp::Neg) {
            auto r = expr_to_f4(un->operand, vars, ctx);
            if (r.is_error()) return r;
            for (auto& [mon, coeff] : r.value().terms) coeff = -coeff;
            return r;
        }
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        auto lhs = expr_to_f4(bin->left, vars, ctx);
        if (lhs.is_error()) return lhs;

        if (bin->op == BinaryOp::Pow) {
            const auto* il = expr_cast<IntegerLit>(bin->right);
            if (!il) {
                return fail<PolyF4>(make_error(
                    CASErrorKind::Unimplemented,
                    "F4 polynomial conversion requires integer exponents"));
            }
            auto exponent = exact_nonnegative_exponent(il->value);
            if (exponent.is_error()) return fail<PolyF4>(exponent.error());
            return pow_poly(lhs.value(), exponent.value(), n);
        }

        auto rhs = expr_to_f4(bin->right, vars, ctx);
        if (rhs.is_error()) return rhs;

        if (bin->op == BinaryOp::Add) return ok(add_poly(lhs.value(), rhs.value()));
        if (bin->op == BinaryOp::Sub) return ok(subtract_poly(lhs.value(), rhs.value()));
        if (bin->op == BinaryOp::Mul) return multiply_poly(lhs.value(), rhs.value(), n);
    }
    if (const auto* s = expr_cast<Sum>(expr)) {
        PolyF4 res;
        for (auto term : s->terms) {
            auto r = expr_to_f4(term, vars, ctx);
            if (r.is_error()) return r;
            res = add_poly(std::move(res), r.value());
        }
        return ok(std::move(res));
    }
    if (const auto* p = expr_cast<Product>(expr)) {
        PolyF4 res;
        res.terms[Monomial(n, 0)] = Rational(1);
        for (auto factor : p->factors) {
            auto r = expr_to_f4(factor, vars, ctx);
            if (r.is_error()) return r;
            auto next = multiply_poly(res, r.value(), n);
            if (next.is_error()) return next;
            res = std::move(next.value());
        }
        return ok(std::move(res));
    }

    return fail<PolyF4>(make_error(CASErrorKind::Unimplemented, "expr_to_f4 unsupported: " + debug_print(expr)));
}

[[nodiscard]] Result<ExprPtr> f4_to_expr(const PolyF4& p, const std::vector<Symbol>& vars, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> sum_terms;
    AstArena& arena = ctx.arena();
    for (const auto& [mon, coeff] : p.terms) {
        ExprPtr term_expr = make_rational_expr(arena, coeff);
        for (std::size_t i = 0; i < vars.size(); ++i) {
            if (mon[i] > 0U) {
                ExprPtr var_expr = arena.make<Symbol>(vars[i].name);
                if (mon[i] > 1U) {
                    var_expr = arena.make<Binary>(BinaryOp::Pow, var_expr, arena.make<IntegerLit>(BigInt(mon[i])));
                }
                term_expr = arena.make<Binary>(BinaryOp::Mul, term_expr, var_expr);
            }
        }
        sum_terms.push_back(term_expr);
    }
    if (sum_terms.empty()) return ok(arena.make<IntegerLit>(BigInt(0)));
    if (sum_terms.size() == 1U) return ok(sum_terms[0]);
    return ok(arena.make<Sum>(sum_terms));
}

} // namespace cas::algebra
