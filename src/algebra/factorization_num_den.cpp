// factorization_num_den.cpp — split_num_den / apart_num_den / together helpers.
// Extracted from factorization_polynomials.cpp (A1 anti-monolith split, F2 Block A).
//
// HC-F8-PENDING-20-RESIDUE (2026-06-13): `together()` now performs polynomial
// GCD content reduction on the (N, D) pair returned by `apart_num_den` so that
// equivalent denominator forms (e.g. y^4+x²·y² vs (x²+y²)·y²) are unified.
// See .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Together_Polynomial_GCD_Reduction.md
// for the formal contract. Reduction is gated by CASContext params
// `together_gcd_enabled` / `together_gcd_max_degree` / `together_gcd_max_symbols`
// and falls back silently to the unreduced rational on any error path —
// `together` remains total.

#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/ast_nodes.hpp"
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace cas::algebra {

// A47 — combinazione STRUTTURALE dei due lati di una frazione.
//
// `multiply_exprs`/`add_exprs`/`subtract_exprs` (algebra_core.cpp) semplificano
// ognuna il proprio risultato, ma qui l'espressione costruita finisce SEMPRE in
// `normalize_rational_parts`, che semplifica numeratore e denominatore prima di
// qualunque test di forma (`is_zero_expr`). Quelle semplificazioni intermedie
// sono quindi lavoro ripetuto sull'accumulato che cresce: sei `simplify` per
// combinazione invece di due, con l'aggregazione a coppie che le paga a ogni
// passo. Misurato su `∫dx/cos²`: `apart_num_den` = 16.5 s su 17.3 s totali di
// test, 42 sole `add_parts` (≈394 ms l'una).
//
// L'invariante resta quello di prima — i predicati di forma girano sull'output
// di `normalize_rational_parts`, che non cambia. Un tentativo precedente (A47)
// aveva provato a togliere anche QUELLE, e li' la logica si rompe davvero:
// `is_zero_expr`/`is_one_expr` si appoggiano alla forma semplificata.
//
// A56 — `add_parts`/`subtract_parts` applicano inoltre la scorciatoia esatta
// `N₁/D + N₂/D = (N₁+N₂)/D` quando i due denominatori sono strutturalmente
// uguali, evitando il cross-moltiplica su D² per il caso più comune (denominatori
// condivisi). Era stata rimossa in A47 per un falso negativo del verificatore
// IBP, chiuso da A54 (vedi commento al sito di `add_parts`).
[[nodiscard]] static ExprPtr combine_mul_raw(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return ctx.arena().make<Product>(std::vector<ExprPtr>{lhs, rhs});
}

[[nodiscard]] static ExprPtr combine_add_raw(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return ctx.arena().make<Sum>(std::vector<ExprPtr>{lhs, rhs});
}

[[nodiscard]] static ExprPtr combine_sub_raw(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs);
}

[[nodiscard]] static Result<RationalParts> normalize_rational_parts(RationalParts parts, symbolic::CASContext& ctx) {
    auto numerator = simplify_expr(parts.numerator, ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }

    auto denominator = simplify_expr(parts.denominator, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }
    if (is_zero_expr(denominator.value())) {
        return fail<RationalParts>(make_error(CASErrorKind::Undefined, "Denominatore nullo in apart_num_den"));
    }

    if (is_zero_expr(numerator.value())) {
        return ok(RationalParts{
            .numerator = make_integer(ctx.arena(), 0),
            .denominator = make_integer(ctx.arena(), 1),
        });
    }

    return ok(RationalParts{
        .numerator = numerator.value(),
        .denominator = denominator.value(),
    });
}

[[nodiscard]] static Result<RationalParts> make_atomic_parts(ExprPtr expr, symbolic::CASContext& ctx) {
    auto cloned = clone_into_context(expr, ctx);
    if (cloned.is_error()) {
        return fail<RationalParts>(cloned.error());
    }
    return normalize_rational_parts(
        RationalParts{
            .numerator = cloned.value(),
            .denominator = make_integer(ctx.arena(), 1),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> multiply_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    return normalize_rational_parts(
        RationalParts{
            .numerator = combine_mul_raw(lhs.numerator, rhs.numerator, ctx),
            .denominator = combine_mul_raw(lhs.denominator, rhs.denominator, ctx),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> divide_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    if (is_zero_expr(rhs.numerator)) {
        return fail<RationalParts>(make_error(CASErrorKind::Undefined, "Divisione per zero in apart_num_den"));
    }

    auto numerator = multiply_exprs(lhs.numerator, rhs.denominator, ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = multiply_exprs(lhs.denominator, rhs.numerator, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }
    return normalize_rational_parts(
        RationalParts{
            .numerator = numerator.value(),
            .denominator = denominator.value(),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> add_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    // A56 — denominatori strutturalmente uguali (già in forma normalizzata da
    // normalize_rational_parts): N₁/D + N₂/D = (N₁+N₂)/D evita il cross-moltiplica
    // su D² e la GCD-reduction a valle che lo disfa. Rientrata dopo A54 (che ha
    // chiuso il falso negativo del verificatore IBP che l'aveva fatta rimuovere
    // in A47): ~23% su together per la famiglia trigonometrica (∫dx/cos²).
    if (structural_equal(lhs.denominator, rhs.denominator)) {
        return normalize_rational_parts(
            RationalParts{
                .numerator = combine_add_raw(lhs.numerator, rhs.numerator, ctx),
                .denominator = lhs.denominator,
            },
            ctx);
    }
    ExprPtr lhs_scaled = combine_mul_raw(lhs.numerator, rhs.denominator, ctx);
    ExprPtr rhs_scaled = combine_mul_raw(rhs.numerator, lhs.denominator, ctx);
    return normalize_rational_parts(
        RationalParts{
            .numerator = combine_add_raw(lhs_scaled, rhs_scaled, ctx),
            .denominator = combine_mul_raw(lhs.denominator, rhs.denominator, ctx),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> subtract_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    // A56 — stesso argomento di add_parts; qui il trigger empirico della
    // famiglia (misurato in A56: il grosso del guadagno arriva da subtract, non da add).
    if (structural_equal(lhs.denominator, rhs.denominator)) {
        return normalize_rational_parts(
            RationalParts{
                .numerator = combine_sub_raw(lhs.numerator, rhs.numerator, ctx),
                .denominator = lhs.denominator,
            },
            ctx);
    }
    ExprPtr lhs_scaled = combine_mul_raw(lhs.numerator, rhs.denominator, ctx);
    ExprPtr rhs_scaled = combine_mul_raw(rhs.numerator, lhs.denominator, ctx);
    return normalize_rational_parts(
        RationalParts{
            .numerator = combine_sub_raw(lhs_scaled, rhs_scaled, ctx),
            .denominator = combine_mul_raw(lhs.denominator, rhs.denominator, ctx),
        },
        ctx);
}

[[nodiscard]] static Result<RationalParts> pow_parts(const RationalParts& parts, const IntegerExponent& exponent, symbolic::CASContext& ctx) {
    if (exponent.negative && is_zero_expr(parts.numerator)) {
        return fail<RationalParts>(make_error(
            CASErrorKind::Undefined,
            "Una potenza negativa richiede una base razionale non nulla"));
    }

    auto numerator = pow_expr(parts.numerator, exponent.magnitude, ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = pow_expr(parts.denominator, exponent.magnitude, ctx);
    if (denominator.is_error()) {
        return fail<RationalParts>(denominator.error());
    }

    if (!exponent.negative) {
        return normalize_rational_parts(
            RationalParts{
                .numerator = numerator.value(),
                .denominator = denominator.value(),
            },
            ctx);
    }

    return normalize_rational_parts(
        RationalParts{
            .numerator = denominator.value(),
            .denominator = numerator.value(),
        },
        ctx);
}

Result<RationalParts> split_num_den(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<RationalParts>(make_error(CASErrorKind::InvalidArgument, "Espressione nulla in apart_num_den"));
    }
    if (contains_decimal_literal(expr)) {
        return fail<RationalParts>(make_error(
            CASErrorKind::Unimplemented,
            "I literal decimali non sono supportati in apart_num_den"));
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            auto operand = split_num_den(unary->operand, ctx);
            if (operand.is_error()) {
                return fail<RationalParts>(operand.error());
            }
            auto negated = negate_expr(operand.value().numerator, ctx);
            if (negated.is_error()) {
                return fail<RationalParts>(negated.error());
            }
            return normalize_rational_parts(
                RationalParts{
                    .numerator = negated.value(),
                    .denominator = operand.value().denominator,
                },
                ctx);
        }
        return make_atomic_parts(expr, ctx);
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        auto lhs = split_num_den(binary->left, ctx);
        if (lhs.is_error()) {
            return fail<RationalParts>(lhs.error());
        }

        switch (binary->op) {
        case BinaryOp::Add: {
            auto rhs = split_num_den(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<RationalParts>(rhs.error());
            }
            return add_parts(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Sub: {
            auto rhs = split_num_den(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<RationalParts>(rhs.error());
            }
            return subtract_parts(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Mul: {
            auto rhs = split_num_den(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<RationalParts>(rhs.error());
            }
            return multiply_parts(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Div: {
            auto rhs = split_num_den(binary->right, ctx);
            if (rhs.is_error()) {
                return fail<RationalParts>(rhs.error());
            }
            return divide_parts(lhs.value(), rhs.value(), ctx);
        }
        case BinaryOp::Pow: {
            auto exponent = parse_integer_exponent(binary->right);
            if (exponent.is_error()) {
                return fail<RationalParts>(exponent.error());
            }
            return pow_parts(lhs.value(), exponent.value(), ctx);
        }
        case BinaryOp::Mod:
            return fail<RationalParts>(make_error(
                CASErrorKind::Unimplemented,
                "Il modulo non e' supportato in apart_num_den"));
        case BinaryOp::Equal:
        case BinaryOp::Less:
        case BinaryOp::Greater:
        case BinaryOp::LessEqual:
        case BinaryOp::GreaterEqual:
            return fail<RationalParts>(make_error(
                CASErrorKind::InvalidArgument,
                "Comparison operators not supported in apart_num_den"));
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        auto result = normalize_rational_parts(
            RationalParts{
                .numerator = make_integer(ctx.arena(), 0),
                .denominator = make_integer(ctx.arena(), 1),
            },
            ctx);
        if (result.is_error()) {
            return fail<RationalParts>(result.error());
        }
        for (ExprPtr term : sum->terms) {
            auto current = split_num_den(term, ctx);
            if (current.is_error()) {
                return fail<RationalParts>(current.error());
            }
            result = add_parts(result.value(), current.value(), ctx);
            if (result.is_error()) {
                return fail<RationalParts>(result.error());
            }
        }
        return result;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        auto result = normalize_rational_parts(
            RationalParts{
                .numerator = make_integer(ctx.arena(), 1),
                .denominator = make_integer(ctx.arena(), 1),
            },
            ctx);
        if (result.is_error()) {
            return fail<RationalParts>(result.error());
        }
        for (ExprPtr factor : product->factors) {
            auto current = split_num_den(factor, ctx);
            if (current.is_error()) {
                return fail<RationalParts>(current.error());
            }
            result = multiply_parts(result.value(), current.value(), ctx);
            if (result.is_error()) {
                return fail<RationalParts>(result.error());
            }
        }
        return result;
    }

    return make_atomic_parts(expr, ctx);
}

// HC-F8-PENDING-20-RESIDUE: walk an ExprPtr and collect every Symbol name into
// `out`. Used to drive `polynomial_exact_divide` (univariate-in-var) over the
// shared indeterminate set of (N, D).
static void collect_symbol_names(ExprPtr expr, std::set<std::string>& out) {
    if (!expr) return;
    if (const auto* s = expr_cast<Symbol>(expr)) {
        out.insert(s->name);
        return;
    }
    if (const auto* u = expr_cast<Unary>(expr)) {
        collect_symbol_names(u->operand, out);
        return;
    }
    if (const auto* b = expr_cast<Binary>(expr)) {
        collect_symbol_names(b->left, out);
        collect_symbol_names(b->right, out);
        return;
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr t : sum->terms) collect_symbol_names(t, out);
        return;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) collect_symbol_names(f, out);
        return;
    }
    if (const auto* fc = expr_cast<FuncCall>(expr)) {
        for (ExprPtr a : fc->args) collect_symbol_names(a, out);
        return;
    }
}

// HC-F8-PENDING-20-RESIDUE: reduce (N, D) by polynomial content GCD if the
// engine can compute it. Returns the *unchanged* pair on any soft failure
// (GCD unimplemented, exact-divide remainder non-zero, guard fired). Hard
// errors are surfaced.
struct ReducedRational {
    ExprPtr numerator;
    ExprPtr denominator;
};
static Result<ReducedRational> reduce_rational_by_gcd(
    ExprPtr N, ExprPtr D, symbolic::CASContext& ctx)
{
    ReducedRational identity{N, D};
    if (!ctx.together_gcd_enabled()) return ok(identity);
    if (!N || !D) return ok(identity);
    if (is_zero_expr(N) || is_one_expr(D) || is_one_expr(N)) return ok(identity);

    // Guard: cap symbol count to avoid catastrophic GCD on wide expressions.
    std::set<std::string> shared;
    collect_symbol_names(N, shared);
    std::set<std::string> dvars;
    collect_symbol_names(D, dvars);
    // Intersect → symbols common to both (potential cancellation indeterminates).
    std::vector<std::string> common;
    std::set_intersection(shared.begin(), shared.end(), dvars.begin(), dvars.end(),
                          std::back_inserter(common));
    if (common.empty()) return ok(identity);
    if (common.size() > ctx.together_gcd_max_symbols()) return ok(identity);

    // Compute the multivariate GCD; soft-fail to identity on any error.
    auto g_res = polynomial_gcd_multivariate(N, D, ctx);
    if (g_res.is_error()) return ok(identity);
    ExprPtr g = g_res.value();
    if (is_zero_expr(g) || is_one_expr(g)) return ok(identity);

    // Pick the lexicographically-first shared symbol as the main variable for
    // univariate exact-divide. Any var present in g works because g divides
    // both N and D exactly; we deliberately choose deterministically.
    Symbol main_var(common.front());

    // Guard: skip reduction if degree in main_var blows up past the cap.
    auto deg_g = polynomial_degree(g, main_var, ctx);
    if (deg_g.is_error()) return ok(identity);
    if (deg_g.value() > ctx.together_gcd_max_degree()) return ok(identity);

    auto N_q = polynomial_exact_divide(N, g, main_var, ctx);
    if (N_q.is_error()) return ok(identity);  // non-exact → fall back unreduced
    auto D_q = polynomial_exact_divide(D, g, main_var, ctx);
    if (D_q.is_error()) return ok(identity);

    return ok(ReducedRational{N_q.value(), D_q.value()});
}

Result<ExprPtr> together(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "together richiede un'espressione non nulla"));
    }

    auto parts = apart_num_den(expr, ctx);
    if (parts.is_error()) {
        return fail<ExprPtr>(parts.error());
    }

    ExprPtr N = parts.value().numerator;
    ExprPtr D = parts.value().denominator;

    // HC-F8-PENDING-20-RESIDUE: polynomial GCD content reduction on (N, D).
    auto reduced = reduce_rational_by_gcd(N, D, ctx);
    if (reduced.is_ok()) {
        N = reduced.value().numerator;
        D = reduced.value().denominator;
    }

    if (is_one_expr(D)) {
        return ok(N);
    }
    return divide_exprs(N, D, ctx);
}

Result<RationalParts> apart_num_den(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<RationalParts>(make_error(
            CASErrorKind::InvalidArgument,
            "apart_num_den richiede un'espressione non nulla"));
    }
    return split_num_den(expr, ctx);
}

} // namespace cas::algebra
