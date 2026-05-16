#include "calculus_internal.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"
#include "../symbolic/simplify_impl.hpp"

#include "cas/algebra.hpp"
#include <iostream>
#include <algorithm>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cas::calculus {

bool MRVCompare::operator()(ExprPtr lhs, ExprPtr rhs) const noexcept {
    return symbolic::canonical_compare(lhs, rhs) < 0;
}

namespace {

void collect_mrv_candidates(ExprPtr e, const Symbol& var, MRVSet& candidates, symbolic::CASContext& ctx) {
    if (!depends_on(e, var)) return;

    if (const auto* sym = expr_cast<Symbol>(e)) {
        if (sym->name == var.name) candidates.insert(e);
        return;
    }

    if (const auto* call = expr_cast<FuncCall>(e)) {
        if (call->func_id == BuiltinOp::Exp) {
            auto lim_u = limit(call->args.front(), var, ctx.arena().make<Constant>(MathConstant::Infinity), LimitDirection::Both, ctx);
            if (lim_u.is_ok() && limit_is_infinity(lim_u.value())) {
                if (expr_is<Unary>(lim_u.value())) {
                    // exp(-inf) -> tends to 0, but its reciprocal exp(inf) is in MRV
                    candidates.insert(ctx.arena().make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{ctx.arena().make<Unary>(UnaryOp::Neg, call->args.front())}));
                } else {
                    candidates.insert(e);
                }
            }
        } else if (call->func_id == BuiltinOp::Ln) {
            candidates.insert(e);
        }
        for (auto arg : call->args) collect_mrv_candidates(arg, var, candidates, ctx);
    } else if (const auto* unary = expr_cast<Unary>(e)) {
        collect_mrv_candidates(unary->operand, var, candidates, ctx);
    } else if (const auto* binary = expr_cast<Binary>(e)) {
        collect_mrv_candidates(binary->left, var, candidates, ctx);
        collect_mrv_candidates(binary->right, var, candidates, ctx);
    } else if (const auto* sum = expr_cast<Sum>(e)) {
        for (auto term : sum->terms) collect_mrv_candidates(term, var, candidates, ctx);
    } else if (const auto* product = expr_cast<Product>(e)) {
        for (auto factor : product->factors) collect_mrv_candidates(factor, var, candidates, ctx);
    }
}

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
    auto get_growth_rank_impl = [&](ExprPtr e, const auto& self) -> int {
        if (!depends_on(e, var)) return 0;
        if (const auto* unary = expr_cast<Unary>(e)) {
            if (unary->op == UnaryOp::Neg) return self(unary->operand, self);
        }
        if (const auto* call = expr_cast<FuncCall>(e); call != nullptr && call->args.size() == 1U) {
            if (call->func_id == BuiltinOp::Exp) {
                auto arg_limit = try_infinite_limit(
                    call->args.front(),
                    var,
                    ctx.arena().make<Constant>(MathConstant::Infinity),
                    ctx.arena());
                // exp of an arg that goes to -infinity decays to 0: rank 0.
                if (arg_limit.is_ok() && limit_is_infinity(arg_limit.value())
                    && expr_is<Unary>(arg_limit.value())) {
                    return 0;
                }
                int inner = self(call->args.front(), self);
                return inner + 1;
            }
            if (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) {
                int inner = self(call->args.front(), self);
                return inner > 1 ? inner - 1 : 0;
            }
        }
        if (const auto* product = expr_cast<Product>(e)) {
            int rank = 0;
            for (ExprPtr factor : product->factors) {
                rank = std::max(rank, self(factor, self));
            }
            return rank;
        }
        if (const auto* sum = expr_cast<Sum>(e)) {
            int rank = 0;
            for (ExprPtr term : sum->terms) {
                rank = std::max(rank, self(term, self));
            }
            return rank;
        }
        if (const auto* binary = expr_cast<Binary>(e)) {
            if (binary->op == BinaryOp::Mul || binary->op == BinaryOp::Div) {
                return std::max(self(binary->left, self), self(binary->right, self));
            }
            if (binary->op == BinaryOp::Add || binary->op == BinaryOp::Sub) {
                return std::max(self(binary->left, self), self(binary->right, self));
            }
            if (binary->op == BinaryOp::Pow) {
                const bool exp_dep = depends_on(binary->right, var);
                const bool base_dep = depends_on(binary->left, var);
                if (!exp_dep) return self(binary->left, self);
                if (!base_dep) return self(binary->right, self) + 1;
                int g = self(binary->right, self);
                int ln_f = self(binary->left, self);
                return g + (ln_f > 1 ? ln_f - 1 : 0) + 1;
            }
        }
        return 2;  // polynomial or mixed expression depending on var
    };
    auto get_growth_rank = [&](ExprPtr e) -> int {
        return get_growth_rank_impl(e, get_growth_rank_impl);
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

struct LeadingPower {
    long long power{};
    ExprPtr coefficient{};
};

[[nodiscard]] bool is_exact_zero(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator.is_zero();
    }
    symbolic::detail::LiteralRational literal;
    auto exact = symbolic::detail::try_get_exact_rational(expr, literal);
    return exact.is_ok() && exact.value() && literal.value.numerator().is_zero();
}

[[nodiscard]] bool is_exact_nonzero(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return !integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return !rational->numerator.is_zero();
    }
    symbolic::detail::LiteralRational literal;
    auto exact = symbolic::detail::try_get_exact_rational(expr, literal);
    return exact.is_ok() && exact.value() && !literal.value.numerator().is_zero();
}

[[nodiscard]] bool is_negative_infinity(ExprPtr point) {
    const auto* unary = expr_cast<Unary>(point);
    return unary != nullptr && unary->op == UnaryOp::Neg && limit_is_infinity(unary->operand);
}

[[nodiscard]] std::optional<long long> integer_value(ExprPtr expr) {
    constexpr std::uint64_t kMaxI64 =
        static_cast<std::uint64_t>(std::numeric_limits<long long>::max());
    constexpr std::uint64_t kMinI64Magnitude = kMaxI64 + 1ULL;

    auto bigint_to_i64 = [&](const BigInt& value) -> std::optional<long long> {
        const std::uint64_t magnitude = value.abs().to_u64();
        if (!value.is_negative()) {
            if (magnitude > kMaxI64) return std::nullopt;
            return static_cast<long long>(magnitude);
        }
        if (magnitude > kMinI64Magnitude) return std::nullopt;
        if (magnitude == kMinI64Magnitude) return std::numeric_limits<long long>::min();
        return -static_cast<long long>(magnitude);
    };

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return bigint_to_i64(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        if (rational->denominator == BigInt(1)) {
            return bigint_to_i64(rational->numerator);
        }
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<ExprPtr> simplify_binary(
    BinaryOp op,
    ExprPtr lhs,
    ExprPtr rhs,
    symbolic::CASContext& ctx) {
    auto expr = ctx.arena().make<Binary>(op, lhs, rhs);
    auto simplified = ctx.simplify(expr);
    if (simplified.is_error()) {
        std::cerr << "simplify_binary error: " << simplified.error().message << " for op " << (int)op << "\n";
        return std::nullopt;
    }
    return simplified.value();
}

[[nodiscard]] std::optional<ExprPtr> simplify_unary_neg(ExprPtr operand, symbolic::CASContext& ctx) {
    auto simplified = ctx.simplify(ctx.arena().make<Unary>(UnaryOp::Neg, operand));
    if (simplified.is_error()) return std::nullopt;
    return simplified.value();
}

[[nodiscard]] std::optional<LeadingPower> leading_power_w(ExprPtr expr, const Symbol& w_var, symbolic::CASContext& ctx);

[[nodiscard]] std::optional<LeadingPower> combine_sum_terms(
    const LeadingPower& lhs,
    const LeadingPower& rhs,
    bool is_sub,
    symbolic::CASContext& ctx) {
    if (lhs.power < rhs.power) return lhs;
    if (rhs.power < lhs.power) {
        if (!is_sub) return rhs;
        auto neg_coeff = simplify_unary_neg(rhs.coefficient, ctx);
        if (!neg_coeff.has_value()) return std::nullopt;
        return LeadingPower{.power = rhs.power, .coefficient = *neg_coeff};
    }

    ExprPtr rhs_coefficient = rhs.coefficient;
    if (is_sub) {
        auto neg_coeff = simplify_unary_neg(rhs.coefficient, ctx);
        if (!neg_coeff.has_value()) return std::nullopt;
        rhs_coefficient = *neg_coeff;
    }
    auto combined = simplify_binary(BinaryOp::Add, lhs.coefficient, rhs_coefficient, ctx);
    if (!combined.has_value()) return std::nullopt;
    if (is_exact_zero(*combined)) {
        // A cancelled leading term requires the next series term; do not guess.
        return std::nullopt;
    }
    if (!is_exact_nonzero(*combined)) {
        return std::nullopt;
    }
    return LeadingPower{.power = lhs.power, .coefficient = *combined};
}

[[nodiscard]] std::optional<LeadingPower> leading_power_w(
    ExprPtr expr,
    const Symbol& w_var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (!depends_on(expr, w_var)) {
        return LeadingPower{.power = 0LL, .coefficient = expr};
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        if (symbol->name == w_var.name) {
            return LeadingPower{
                .power = 1LL,
                .coefficient = limit_make_integer(arena, 1),
            };
        }
        return std::nullopt;
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) return std::nullopt;
        auto inner = leading_power_w(unary->operand, w_var, ctx);
        if (!inner.has_value()) return std::nullopt;
        auto coeff = simplify_unary_neg(inner->coefficient, ctx);
        if (!coeff.has_value()) return std::nullopt;
        return LeadingPower{.power = inner->power, .coefficient = *coeff};
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Add || binary->op == BinaryOp::Sub) {
            auto left = leading_power_w(binary->left, w_var, ctx);
            auto right = leading_power_w(binary->right, w_var, ctx);
            if (!left.has_value() || !right.has_value()) return std::nullopt;
            return combine_sum_terms(*left, *right, binary->op == BinaryOp::Sub, ctx);
        }
        if (binary->op == BinaryOp::Mul || binary->op == BinaryOp::Div) {
            auto left = leading_power_w(binary->left, w_var, ctx);
            auto right = leading_power_w(binary->right, w_var, ctx);
            if (!left.has_value() || !right.has_value()) return std::nullopt;
            auto coeff = simplify_binary(binary->op, left->coefficient, right->coefficient, ctx);
            if (!coeff.has_value()) return std::nullopt;
            const long long power = binary->op == BinaryOp::Mul
                ? left->power + right->power
                : left->power - right->power;
            return LeadingPower{.power = power, .coefficient = *coeff};
        }
        if (binary->op == BinaryOp::Pow) {
            auto exponent = integer_value(binary->right);
            if (!exponent.has_value()) return std::nullopt;
            auto base = leading_power_w(binary->left, w_var, ctx);
            if (!base.has_value()) return std::nullopt;
            auto coeff = simplify_binary(
                BinaryOp::Pow,
                base->coefficient,
                limit_make_integer(arena, *exponent),
                ctx);
            if (!coeff.has_value()) return std::nullopt;
            return LeadingPower{
                .power = base->power * (*exponent),
                .coefficient = *coeff,
            };
        }
        return std::nullopt;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        LeadingPower result{
            .power = 0LL,
            .coefficient = limit_make_integer(arena, 1),
        };
        for (ExprPtr factor : product->factors) {
            auto term = leading_power_w(factor, w_var, ctx);
            if (!term.has_value()) return std::nullopt;
            result.power += term->power;
            auto coeff = simplify_binary(BinaryOp::Mul, result.coefficient, term->coefficient, ctx);
            if (!coeff.has_value()) return std::nullopt;
            result.coefficient = *coeff;
        }
        return result;
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.empty()) {
            return LeadingPower{
                .power = 0LL,
                .coefficient = limit_make_integer(arena, 0),
            };
        }
        auto current = leading_power_w(sum->terms.front(), w_var, ctx);
        if (!current.has_value()) return std::nullopt;
        for (std::size_t index = 1; index < sum->terms.size(); ++index) {
            auto next = leading_power_w(sum->terms[index], w_var, ctx);
            if (!next.has_value()) return std::nullopt;
            auto merged = combine_sum_terms(*current, *next, false, ctx);
            if (!merged.has_value()) {
                // A cancelled dominant term needs the next Laurent term. Factor out
                // the current valuation first, then use Taylor on the regular part.
                const long long scale_power = -current->power;
                ExprPtr scaled_expr = expr;
                if (scale_power != 0LL) {
                    ExprPtr w_expr = arena.make<Symbol>(w_var.name);
                    ExprPtr scale = scale_power == 1LL
                        ? w_expr
                        : arena.make<Binary>(
                            BinaryOp::Pow,
                            w_expr,
                            limit_make_integer(arena, scale_power));
                    auto scaled = ctx.simplify(arena.make<Binary>(BinaryOp::Mul, expr, scale));
                    if (scaled.is_error()) return std::nullopt;
                    scaled_expr = scaled.value();
                }
                auto series = taylor_series(scaled_expr, w_var, limit_make_integer(arena, 0), 6U, ctx);
                if (series.is_error()) return std::nullopt;
                if (structural_equal(series.value().polynomial, scaled_expr)) return std::nullopt;
                auto regular = leading_power_w(series.value().polynomial, w_var, ctx);
                if (!regular.has_value()) return std::nullopt;
                return LeadingPower{
                    .power = current->power + regular->power,
                    .coefficient = regular->coefficient,
                };
            }
            current = merged;
        }
        return current;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        std::vector<ExprPtr> args_coeffs;
        for (ExprPtr arg : call->args) {
            auto arg_leading = leading_power_w(arg, w_var, ctx);
            if (!arg_leading.has_value()) return std::nullopt;
            if (arg_leading->power < 0) return std::nullopt;
            if (arg_leading->power > 0) {
                args_coeffs.push_back(limit_make_integer(arena, 0));
            } else {
                args_coeffs.push_back(arg_leading->coefficient);
            }
        }
        auto res_coeff = ctx.simplify(arena.make<FuncCall>(call->func_id, std::move(args_coeffs)));
        if (res_coeff.is_error()) return std::nullopt;
        return LeadingPower{.power = 0LL, .coefficient = res_coeff.value()};
    }

    return std::nullopt;
}

[[nodiscard]] std::optional<int> exact_sign(ExprPtr expr) {
    symbolic::detail::LiteralRational lr;
    auto exact = symbolic::detail::try_get_exact_rational(expr, lr);
    if (exact.is_error() || !exact.value()) return std::nullopt;
    if (lr.value.numerator().is_zero()) return 0;
    return lr.value.numerator().is_negative() ? -1 : 1;
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_quotient_valuation_limit(
    ExprPtr expr,
    const Symbol& w_var,
    AstArena& arena,
    symbolic::CASContext& ctx) {
    auto quotient = extract_quotient_view(expr, arena);
    if (!quotient.has_value()) return std::nullopt;

    auto numerator = leading_power_w(quotient->numerator, w_var, ctx);
    auto denominator = leading_power_w(quotient->denominator, w_var, ctx);
    if (!numerator.has_value() || !denominator.has_value()) return std::nullopt;

    if (numerator->power > denominator->power) {
        return ok(limit_make_integer(arena, 0));
    }

    if (numerator->power == denominator->power) {
        auto ratio = ctx.simplify(arena.make<Binary>(
            BinaryOp::Div,
            numerator->coefficient,
            denominator->coefficient));
        if (ratio.is_error()) return ratio;
        if (depends_on(ratio.value(), w_var)) return std::nullopt;
        return ratio;
    }

    auto ratio = ctx.simplify(arena.make<Binary>(
        BinaryOp::Div,
        numerator->coefficient,
        denominator->coefficient));
    if (ratio.is_error()) return ratio;
    auto sign = exact_sign(ratio.value());
    if (!sign.has_value() || *sign == 0) return std::nullopt;
    if (*sign > 0) return ok(arena.make<Constant>(MathConstant::Infinity));
    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_leading_power_limit(
    ExprPtr expr,
    const Symbol& w_var,
    AstArena& arena,
    symbolic::CASContext& ctx) {
    auto leading = leading_power_w(expr, w_var, ctx);
    if (!leading.has_value()) return std::nullopt;

    if (leading->power > 0) {
        return ok(limit_make_integer(arena, 0));
    }
    if (leading->power == 0) {
        auto coeff = ctx.simplify(leading->coefficient);
        if (coeff.is_ok() && !depends_on(coeff.value(), w_var)) {
            return coeff;
        }
        return std::nullopt;
    }

    auto coeff = ctx.simplify(leading->coefficient);
    if (coeff.is_error()) return coeff;
    auto sign = exact_sign(coeff.value());
    if (!sign.has_value() || *sign == 0) return std::nullopt;
    if (*sign > 0) {
        return ok(arena.make<Constant>(MathConstant::Infinity));
    }
    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
}

struct ExponentialTerm {
    ExprPtr exponent{};
    long long coefficient{};
};

[[nodiscard]] bool safe_mul_i64(long long lhs, long long rhs, long long& out) {
    if (lhs == 0 || rhs == 0) {
        out = 0;
        return true;
    }
    if (lhs == std::numeric_limits<long long>::min() || rhs == std::numeric_limits<long long>::min()) {
        return false;
    }
    const auto abs_lhs = lhs < 0 ? -lhs : lhs;
    const auto abs_rhs = rhs < 0 ? -rhs : rhs;
    if (abs_lhs > std::numeric_limits<long long>::max() / abs_rhs) return false;
    out = lhs * rhs;
    return true;
}

[[nodiscard]] bool append_scaled_exponential_term(
    std::vector<ExponentialTerm>& terms,
    ExprPtr exponent,
    long long coefficient,
    AstArena& arena) {
    if (coefficient == 0) return true;
    if (const auto* unary = expr_cast<Unary>(exponent)) {
        if (unary->op == UnaryOp::Neg) {
            if (coefficient == std::numeric_limits<long long>::min()) return false;
            return append_scaled_exponential_term(terms, unary->operand, -coefficient, arena);
        }
    }

    if (const auto* binary = expr_cast<Binary>(exponent)) {
        if (binary->op == BinaryOp::Mul) {
            if (auto lhs_coeff = integer_value(binary->left)) {
                long long scaled{};
                if (!safe_mul_i64(coefficient, *lhs_coeff, scaled)) return false;
                return append_scaled_exponential_term(terms, binary->right, scaled, arena);
            }
            if (auto rhs_coeff = integer_value(binary->right)) {
                long long scaled{};
                if (!safe_mul_i64(coefficient, *rhs_coeff, scaled)) return false;
                return append_scaled_exponential_term(terms, binary->left, scaled, arena);
            }
        }
    }

    if (const auto* product = expr_cast<Product>(exponent)) {
        long long scaled = coefficient;
        std::vector<ExprPtr> symbolic_factors;
        for (ExprPtr factor : product->factors) {
            if (auto factor_coeff = integer_value(factor)) {
                if (!safe_mul_i64(scaled, *factor_coeff, scaled)) return false;
            } else {
                symbolic_factors.push_back(factor);
            }
        }
        if (symbolic_factors.empty()) return true;
        ExprPtr symbolic_exponent = symbolic_factors.size() == 1U
            ? symbolic_factors.front()
            : arena.make<Product>(std::move(symbolic_factors));
        return append_scaled_exponential_term(terms, symbolic_exponent, scaled, arena);
    }

    terms.push_back(ExponentialTerm{.exponent = exponent, .coefficient = coefficient});
    return true;
}

[[nodiscard]] bool collect_exponential_product_terms(
    ExprPtr expr,
    const Symbol& var,
    long long multiplier,
    std::vector<ExponentialTerm>& terms,
    AstArena& arena) {
    if (multiplier == 0) return true;
    if (!expr) return false;

    if (!depends_on(expr, var)) {
        return integer_value(expr).has_value();
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::Exp && call->args.size() == 1U) {
            return append_scaled_exponential_term(terms, call->args.front(), multiplier, arena);
        }
        return false;
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) return false;
        return collect_exponential_product_terms(unary->operand, var, multiplier, terms, arena);
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Mul) {
            return collect_exponential_product_terms(binary->left, var, multiplier, terms, arena) &&
                   collect_exponential_product_terms(binary->right, var, multiplier, terms, arena);
        }
        if (binary->op == BinaryOp::Div) {
            long long denominator_multiplier{};
            if (!safe_mul_i64(multiplier, -1LL, denominator_multiplier)) return false;
            return collect_exponential_product_terms(binary->left, var, multiplier, terms, arena) &&
                   collect_exponential_product_terms(binary->right, var, denominator_multiplier, terms, arena);
        }
        if (binary->op == BinaryOp::Pow) {
            const auto exponent = integer_value(binary->right);
            if (!exponent.has_value()) return false;
            long long scaled{};
            if (!safe_mul_i64(multiplier, *exponent, scaled)) return false;
            return collect_exponential_product_terms(binary->left, var, scaled, terms, arena);
        }
        return false;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (!collect_exponential_product_terms(factor, var, multiplier, terms, arena)) return false;
        }
        return true;
    }

    return integer_value(expr).has_value();
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_exponential_product_limit(
    ExprPtr expr,
    const Symbol& var,
    AstArena& arena,
    symbolic::CASContext& ctx) {
    std::vector<ExponentialTerm> raw_terms;
    if (!collect_exponential_product_terms(expr, var, 1LL, raw_terms, arena) || raw_terms.empty()) {
        return std::nullopt;
    }

    std::vector<ExponentialTerm> terms;
    for (const auto& raw : raw_terms) {
        if (!depends_on(raw.exponent, var)) continue;
        bool merged = false;
        for (auto& term : terms) {
            if (structural_equal(term.exponent, raw.exponent)) {
                term.coefficient += raw.coefficient;
                merged = true;
                break;
            }
        }
        if (!merged) terms.push_back(raw);
    }

    terms.erase(
        std::remove_if(terms.begin(), terms.end(), [](const ExponentialTerm& term) {
            return term.coefficient == 0;
        }),
        terms.end());
    if (terms.empty()) return std::nullopt;

    while (!terms.empty()) {
        std::size_t dominant = 0U;
        bool undecidable_tie = false;
        for (std::size_t i = 1U; i < terms.size(); ++i) {
            const int cmp = compare_growth(terms[i].exponent, terms[dominant].exponent, var, ctx);
            if (cmp > 0) {
                dominant = i;
                undecidable_tie = false;
            } else if (cmp == 0 && !structural_equal(terms[i].exponent, terms[dominant].exponent)) {
                undecidable_tie = true;
            }
        }

        long long dominant_coefficient = terms[dominant].coefficient;
        for (std::size_t i = 0U; i < terms.size(); ++i) {
            if (i == dominant) continue;
            if (structural_equal(terms[i].exponent, terms[dominant].exponent)) {
                dominant_coefficient += terms[i].coefficient;
            }
        }

        if (undecidable_tie) return std::nullopt;
        if (dominant_coefficient > 0) return ok(arena.make<Constant>(MathConstant::Infinity));
        if (dominant_coefficient < 0) return ok(limit_make_integer(arena, 0));

        ExprPtr dominant_expr = terms[dominant].exponent;
        terms.erase(
            std::remove_if(terms.begin(), terms.end(), [&](const ExponentialTerm& term) {
                return structural_equal(term.exponent, dominant_expr);
            }),
            terms.end());
    }

    return std::nullopt;
}

} // namespace

MRVSet mrv_set(ExprPtr e, const Symbol& var, symbolic::CASContext& ctx) {
    MRVSet candidates;
    collect_mrv_candidates(e, var, candidates, ctx);
    if (candidates.empty()) return {};

    MRVSet mrv;
    for (auto c : candidates) {
        if (mrv.empty()) {
            mrv.insert(c);
        } else {
            int cmp = compare_growth(c, *mrv.begin(), var, ctx);
            if (cmp > 0) {
                mrv.clear();
                mrv.insert(c);
            } else if (cmp == 0) {
                mrv.insert(c);
            }
        }
    }
    return mrv;
}

Result<ExprPtr> rewrite_mrv(ExprPtr e, const MRVSet& mrv, ExprPtr w, const Symbol& var, symbolic::CASContext& ctx) {
    for (ExprPtr m : mrv) {
        if (structural_equal(e, m)) return ok(w);
    }

    if (expr_is<Symbol>(e) || expr_is<IntegerLit>(e) || expr_is<RationalLit>(e) || 
        expr_is<Constant>(e) || expr_is<DecimalLit>(e)) {
        return ok(e);
    }

    if (const auto* call = expr_cast<FuncCall>(e)) {
        if (call->func_id == BuiltinOp::Exp) {
            ExprPtr neg_arg = ctx.arena().make<Unary>(UnaryOp::Neg, call->args.front());
            auto s_neg = ctx.simplify(neg_arg);
            if (s_neg.is_ok()) {
                ExprPtr exp_neg = ctx.arena().make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{s_neg.value()});
                for (ExprPtr m : mrv) {
                    if (structural_equal(exp_neg, m)) {
                        // If the reciprocal MRV element exp(-u) maps to 1/w,
                        // then exp(u) maps to w.
                        if (const auto* bin = expr_cast<Binary>(w)) {
                            if (bin->op == BinaryOp::Div) return ok(bin->right);
                        }
                    }
                }
            }
        }
        std::vector<ExprPtr> new_args;
        for (auto arg : call->args) {
            auto res = rewrite_mrv(arg, mrv, w, var, ctx);
            if (res.is_error()) return res;
            new_args.push_back(res.value());
        }
        return ok(ctx.arena().make<FuncCall>(call->name, std::move(new_args)));
    }
    if (const auto* unary = expr_cast<Unary>(e)) {
        auto res = rewrite_mrv(unary->operand, mrv, w, var, ctx);
        if (res.is_error()) return res;
        return ok(ctx.arena().make<Unary>(unary->op, res.value()));
    }
    if (const auto* binary = expr_cast<Binary>(e)) {
        auto l = rewrite_mrv(binary->left, mrv, w, var, ctx);
        if (l.is_error()) return l;
        auto r = rewrite_mrv(binary->right, mrv, w, var, ctx);
        if (r.is_error()) return r;
        return ok(ctx.arena().make<Binary>(binary->op, l.value(), r.value()));
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> new_terms;
        for (auto t : sum->terms) {
            auto res = rewrite_mrv(t, mrv, w, var, ctx);
            if (res.is_error()) return res;
            new_terms.push_back(res.value());
        }
        return ok(ctx.arena().make<Sum>(std::move(new_terms)));
    }
    if (const auto* product = expr_cast<Product>(e)) {
        std::vector<ExprPtr> new_factors;
        for (auto f : product->factors) {
            auto res = rewrite_mrv(f, mrv, w, var, ctx);
            if (res.is_error()) return res;
            new_factors.push_back(res.value());
        }
        return ok(ctx.arena().make<Product>(std::move(new_factors)));
    }

    return ok(e);
}

Result<ExprPtr> compute_limit_mrv(ExprPtr expr, const Symbol& var, ExprPtr point, symbolic::CASContext& ctx) {
    if (!limit_is_infinity(point)) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV is currently defined only for infinite limits",
            .hint = std::nullopt
        });
    }

    if (is_negative_infinity(point)) {
        AstArena& arena = ctx.arena();
        ExprPtr neg_var = arena.make<Unary>(UnaryOp::Neg, arena.make<Symbol>(var.name));
        auto transformed = ctx.substitute(expr, var, neg_var);
        if (transformed.is_error()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "MRV normalization for -infinity failed",
                .hint = std::nullopt
            });
        }
        auto transformed_simplified = ctx.simplify(transformed.value());
        if (transformed_simplified.is_error()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "MRV normalization simplify failed",
                .hint = std::nullopt
            });
        }
        return compute_limit_mrv(
            transformed_simplified.value(),
            var,
            arena.make<Constant>(MathConstant::Infinity),
            ctx);
    }

    if (auto exp_product_limit = try_exponential_product_limit(expr, var, ctx.arena(), ctx);
        exp_product_limit.has_value()) {
        return *exp_product_limit;
    }

    MRVSet mrv = mrv_set(expr, var, ctx);
    if (mrv.empty()) return ok(expr);

    AstArena& arena = ctx.arena();
    // Fresh MRV variable.  Iterates make_fresh_symbol until the candidate
    // is structurally absent from the input expression — guards against
    // capture even when the user pre-populated the context with
    // similarly-prefixed symbols.
    Symbol w_var = ctx.make_fresh_symbol("mrv_w");
    while (depends_on(expr, w_var)) {
        w_var = ctx.make_fresh_symbol("mrv_w");
    }
    std::string w_name = w_var.name;
    ExprPtr w_sym = arena.make<Symbol>(w_name);

    ExprPtr replacement = arena.make<Binary>(BinaryOp::Div, limit_make_integer(arena, 1), w_sym);
    auto rewritten = rewrite_mrv(expr, mrv, replacement, var, ctx);
    if (rewritten.is_error()) {
        std::cerr << "rewrite_mrv error\n";
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV rewrite failed",
            .hint = std::nullopt
        });
    }

    auto simplified = ctx.simplify(rewritten.value());
    if (simplified.is_error()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV simplified form unavailable",
            .hint = std::nullopt
        });
    }

    if (auto quotient_limit = try_quotient_valuation_limit(simplified.value(), w_var, arena, ctx);
        quotient_limit.has_value()) {
        return *quotient_limit;
    }

    if (auto leading_limit = try_leading_power_limit(simplified.value(), w_var, arena, ctx);
        leading_limit.has_value()) {
        return *leading_limit;
    }

    auto combined_res = algebra::together(simplified.value(), ctx);
    ExprPtr final_rewritten = combined_res.is_ok() ? combined_res.value() : simplified.value();

    if (auto quotient_limit = try_quotient_valuation_limit(final_rewritten, w_var, arena, ctx);
        quotient_limit.has_value()) {
        return *quotient_limit;
    }

    auto direct = ctx.substitute(final_rewritten, w_var, limit_make_integer(arena, 0));
    if (direct.is_ok()) {
        auto direct_simplified = ctx.simplify(direct.value());
        if (direct_simplified.is_ok() && !depends_on(direct_simplified.value(), w_var)) {
            return direct_simplified;
        }
    }

    auto leading = leading_power_w(final_rewritten, w_var, ctx);
    if (!leading.has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV leading power is not decidable",
            .hint = std::nullopt
        });
    }
    if (leading->power > 0) {
        return ok(limit_make_integer(arena, 0));
    }
    if (leading->power == 0) {
        auto coeff = ctx.simplify(leading->coefficient);
        if (coeff.is_ok() && !depends_on(coeff.value(), w_var)) {
            return coeff;
        }
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV leading coefficient is not decidable at w -> 0+",
            .hint = std::nullopt
        });
    }

    auto coeff = ctx.simplify(leading->coefficient);
    if (coeff.is_error()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV leading coefficient simplify failed",
            .hint = std::nullopt
        });
    }
    auto sign = exact_sign(coeff.value());
    if (!sign.has_value() || *sign == 0) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "MRV pole sign is not exactly decidable",
            .hint = std::nullopt
        });
    }
    if (*sign > 0) {
        return ok(arena.make<Constant>(MathConstant::Infinity));
    }
    return ok(arena.make<Unary>(UnaryOp::Neg, arena.make<Constant>(MathConstant::Infinity)));
}

} // namespace cas::calculus
