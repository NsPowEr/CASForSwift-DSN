#include "cas/algebra.hpp"
#include "cas/ast_debug.hpp"
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <vector>

namespace cas::algebra {

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
    auto numerator = multiply_exprs(lhs.numerator, rhs.numerator, ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = multiply_exprs(lhs.denominator, rhs.denominator, ctx);
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
    auto lhs_scaled = multiply_exprs(lhs.numerator, rhs.denominator, ctx);
    if (lhs_scaled.is_error()) {
        return fail<RationalParts>(lhs_scaled.error());
    }
    auto rhs_scaled = multiply_exprs(rhs.numerator, lhs.denominator, ctx);
    if (rhs_scaled.is_error()) {
        return fail<RationalParts>(rhs_scaled.error());
    }
    auto numerator = add_exprs(lhs_scaled.value(), rhs_scaled.value(), ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = multiply_exprs(lhs.denominator, rhs.denominator, ctx);
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

[[nodiscard]] static Result<RationalParts> subtract_parts(const RationalParts& lhs, const RationalParts& rhs, symbolic::CASContext& ctx) {
    auto lhs_scaled = multiply_exprs(lhs.numerator, rhs.denominator, ctx);
    if (lhs_scaled.is_error()) {
        return fail<RationalParts>(lhs_scaled.error());
    }
    auto rhs_scaled = multiply_exprs(rhs.numerator, lhs.denominator, ctx);
    if (rhs_scaled.is_error()) {
        return fail<RationalParts>(rhs_scaled.error());
    }
    auto numerator = subtract_exprs(lhs_scaled.value(), rhs_scaled.value(), ctx);
    if (numerator.is_error()) {
        return fail<RationalParts>(numerator.error());
    }
    auto denominator = multiply_exprs(lhs.denominator, rhs.denominator, ctx);
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

Result<ExprPtr> together(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "together richiede un'espressione non nulla"));
    }

    auto parts = apart_num_den(expr, ctx);
    if (parts.is_error()) {
        return fail<ExprPtr>(parts.error());
    }
    if (is_one_expr(parts.value().denominator)) {
        return ok(parts.value().numerator);
    }
    return divide_exprs(parts.value().numerator, parts.value().denominator, ctx);
}

Result<RationalParts> apart_num_den(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return fail<RationalParts>(make_error(
            CASErrorKind::InvalidArgument,
            "apart_num_den richiede un'espressione non nulla"));
    }
    return split_num_den(expr, ctx);
}

Result<ExprPtr> polynomial_gcd(ExprPtr p, ExprPtr q, const Symbol& var, symbolic::CASContext& ctx) {
    const bool owns_operation = !ctx.operation_active_;
    if (owns_operation) {
        ctx.operation_active_ = true;
        ctx.trace_capture_active_ = ctx.trace_enabled_;
        ctx.trace_.clear();
        ctx.ops_count_ = 0;
        ctx.operation_started_at_ = std::chrono::steady_clock::now();
    }

    auto record_trace = [&](symbolic::RuleId rule_id, ExprPtr result) {
        if (!ctx.trace_enabled_) {
            return;
        }
        ctx.trace_.push_back(symbolic::TraceStep{
            .rule_id = rule_id,
            .depth = 0U,
            .target_before = p,
            .target_after = q,
            .root_after = result,
        });
    };

    auto finalize = [&]() {
        if (owns_operation) {
            ctx.operation_active_ = false;
            ctx.trace_capture_active_ = false;
            ctx.ops_count_ = 0;
        }
    };

    auto result = [&]() -> Result<ExprPtr> {
    if (!p || !q) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd richiede due espressioni polinomiali non nulle"));
    }

    // Try multivariate GCD first if there are multiple variables
    auto p_multi = parse_multivariate_polynomial(p, ctx);
    auto q_multi = parse_multivariate_polynomial(q, ctx);
    
    if (p_multi.is_ok() && q_multi.is_ok()) {
        auto vars = p_multi.value().variables();
        auto vars_q = q_multi.value().variables();
        vars.insert(vars.end(), vars_q.begin(), vars_q.end());
        std::sort(vars.begin(), vars.end(), [](const auto& a, const auto& b){ return a.name < b.name; });
        vars.erase(std::unique(vars.begin(), vars.end(), [](const auto& a, const auto& b){ return a.name == b.name; }), vars.end());

        if (vars.size() > 1 || (vars.size() == 1 && vars[0].name != var.name)) {
            // Use Heuristic GCD
            auto gcd_res = gcd_heuristic(p_multi.value(), q_multi.value());
            if (gcd_res.is_error()) {
                // Fallback to Modular GCD
                gcd_res = gcd_modular(p_multi.value(), q_multi.value());
            }

            if (gcd_res.is_ok()) {
                return multivariate_to_expr(gcd_res.value(), ctx);
            }
        }
    }

    auto left = parse_polynomial(p, var, ctx);
    if (left.is_error()) {
        return fail<ExprPtr>(left.error());
    }
    auto right = parse_polynomial(q, var, ctx);
    if (right.is_error()) {
        return fail<ExprPtr>(right.error());
    }

    auto left_integer = poly_to_integer_poly(left.value());
    auto right_integer = poly_to_integer_poly(right.value());
    if (left_integer.is_ok() && right_integer.is_ok()) {
        IntegerGcdResult gcd_result =
            gcd_integer_poly_with_subresultant(left_integer.value(), right_integer.value());
        IntPoly gcd_poly = std::move(gcd_result.gcd);
        symbolic::RuleId path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallback;
        switch (gcd_result.path) {
        case IntegerGcdPath::Subresultant:
            path_rule = symbolic::RuleId::PolynomialGcdSubresultant;
            break;
        case IntegerGcdPath::PrimitiveFallbackPsi:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallbackPsi;
            break;
        case IntegerGcdPath::PrimitiveFallbackBeta:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallbackBeta;
            break;
        case IntegerGcdPath::PrimitiveFallback:
            path_rule = symbolic::RuleId::PolynomialGcdPrimitiveFallback;
            break;
        }

        if (is_zero_integer_poly(gcd_poly)) {
            ExprPtr zero = make_integer(ctx.arena(), 0);
            record_trace(path_rule, zero);
            return ok(zero);
        }

        auto gcd_expr = integer_coefficients_to_expr(gcd_poly, var, ctx);
        if (gcd_expr.is_error()) {
            return fail<ExprPtr>(gcd_expr.error());
        }

        auto traced_result = ctx.simplify(gcd_expr.value());
        if (traced_result.is_ok()) {
            record_trace(path_rule, traced_result.value());
        }
        return traced_result;
    }

    if (is_zero_poly(left.value()) && is_zero_poly(right.value())) {
        ExprPtr zero = make_integer(ctx.arena(), 0);
        record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, zero);
        return ok(zero);
    }
    if (is_zero_poly(left.value())) {
        auto normalized = normalize_poly_monic(right.value(), ctx);
        if (normalized.is_error()) {
            return fail<ExprPtr>(normalized.error());
        }
        auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
        if (traced_result.is_ok()) {
            record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
        }
        return traced_result;
    }
    if (is_zero_poly(right.value())) {
        auto normalized = normalize_poly_monic(left.value(), ctx);
        if (normalized.is_error()) {
            return fail<ExprPtr>(normalized.error());
        }
        auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
        if (traced_result.is_ok()) {
            record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
        }
        return traced_result;
    }

    PolyExpr a = std::move(left.value());
    PolyExpr b = std::move(right.value());

    while (!is_zero_poly(b)) {
        auto division = divide_poly_with_remainder(a, b, ctx);
        if (division.is_error()) {
            return fail<ExprPtr>(division.error());
        }
        a = std::move(b);
        b = std::move(division.value().remainder);
    }

    auto normalized = normalize_poly_monic(a, ctx);
    if (normalized.is_error()) {
        return fail<ExprPtr>(normalized.error());
    }
    auto traced_result = polynomial_to_expr(normalized.value(), var, ctx);
    if (traced_result.is_ok()) {
        record_trace(symbolic::RuleId::PolynomialGcdSymbolicEuclidean, traced_result.value());
    }
    return traced_result;
    }();

    finalize();
    return result;
}

namespace {

struct ExtensionInfo {
    Symbol var;
    RatPoly min_poly;
    ExprPtr original_expr;
};

[[nodiscard]] ExprPtr replace_subexpression(ExprPtr expr, ExprPtr target, ExprPtr replacement, AstArena& arena) {
    if (!expr) {
        return ExprPtr{};
    }
    if (structural_equal(expr, target)) {
        return replacement;
    }

    return visit_expr(
        expr,
        [&](const auto& node) -> ExprPtr {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, Symbol> ||
                std::is_same_v<Node, Constant>) {
                return arena.make<Node>(node);
            } else if constexpr (std::is_same_v<Node, Unary>) {
                return arena.make<Unary>(node.op, replace_subexpression(node.operand, target, replacement, arena));
            } else if constexpr (std::is_same_v<Node, Binary>) {
                return arena.make<Binary>(
                    node.op,
                    replace_subexpression(node.left, target, replacement, arena),
                    replace_subexpression(node.right, target, replacement, arena));
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> args;
                args.reserve(node.args.size());
                for (ExprPtr arg : node.args) {
                    args.push_back(replace_subexpression(arg, target, replacement, arena));
                }
                return arena.make<FuncCall>(node.name, std::move(args));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> terms;
                terms.reserve(node.terms.size());
                for (ExprPtr term : node.terms) {
                    terms.push_back(replace_subexpression(term, target, replacement, arena));
                }
                return arena.make<Sum>(std::move(terms));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> factors;
                factors.reserve(node.factors.size());
                for (ExprPtr factor : node.factors) {
                    factors.push_back(replace_subexpression(factor, target, replacement, arena));
                }
                return arena.make<Product>(std::move(factors));
            } else if constexpr (std::is_same_v<Node, Integral>) {
                return arena.make<Integral>(
                    replace_subexpression(node.integrand, target, replacement, arena),
                    node.variable,
                    node.lower.has_value()
                        ? std::optional<ExprPtr>(replace_subexpression(*node.lower, target, replacement, arena))
                        : std::nullopt,
                    node.upper.has_value()
                        ? std::optional<ExprPtr>(replace_subexpression(*node.upper, target, replacement, arena))
                        : std::nullopt);
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                return arena.make<Derivative>(
                    replace_subexpression(node.expression, target, replacement, arena),
                    node.variable,
                    node.order);
            } else if constexpr (std::is_same_v<Node, Limit>) {
                return arena.make<Limit>(
                    replace_subexpression(node.expression, target, replacement, arena),
                    node.variable,
                    replace_subexpression(node.point, target, replacement, arena),
                    node.direction);
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                return arena.make<RootOf>(
                    replace_subexpression(node.polynomial, target, replacement, arena),
                    node.variable,
                    node.root_index);
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::vector<ExprPtr> elements;
                elements.reserve(node.elements.size());
                for (ExprPtr element : node.elements) {
                    elements.push_back(replace_subexpression(element, target, replacement, arena));
                }
                return arena.make<Matrix>(node.rows, node.cols, std::move(elements));
            } else {
                return ExprPtr{};
            }
        });
}

[[nodiscard]] Result<ExtensionInfo> get_extension_info(ExprPtr extension, symbolic::CASContext& ctx) {
    // Caso sqrt(n) o Pow(n, 1/2)
    auto try_parse_radical = [&](ExprPtr base, ExprPtr exp) -> std::optional<ExtensionInfo> {
        auto parsed_exp = parse_polynomial(exp, Symbol{"_unused"}, ctx);
        auto exp_rational = poly_to_rational_poly(parsed_exp.is_ok() ? parsed_exp.value() : PolyExpr{});
        if (exp_rational.is_ok() && exp_rational.value().size() == 1) {
            Rational r = exp_rational.value().constant_term();
            if (r.numerator() == BigInt(1) && r.denominator() == BigInt(2)) {
                auto base_poly = parse_polynomial(base, Symbol{"_unused"}, ctx);
                if (base_poly.is_ok()) {
                    auto base_rat = poly_to_rational_poly(base_poly.value());
                    if (base_rat.is_ok() && base_rat.value().size() == 1) {
                        Rational b = base_rat.value().constant_term();
                        RatPoly mp;
                        mp.push_back(-b);
                        mp.push_back(Rational(0));
                        mp.push_back(Rational(1));
                        return ExtensionInfo{Symbol{"__alpha"}, mp, extension};
                    }
                }
            }
        }
        return std::nullopt;
    };

    if (const auto* binary = expr_cast<Binary>(extension)) {
        if (binary->op == BinaryOp::Pow) {
            if (auto info = try_parse_radical(binary->left, binary->right)) return ok(*info);
        }
    }
    if (const auto* func = expr_cast<FuncCall>(extension)) {
        if (func->func_id == BuiltinOp::Sqrt && func->args.size() == 1) {
            if (auto info = try_parse_radical(func->args[0], ctx.arena().make<RationalLit>(BigInt(1), BigInt(2)))) return ok(*info);
        }
    }

    // Se l'estensione è già un simbolo, forse il polinomio minimo è fornito altrove?
    // Per ora supportiamo solo radicali espliciti.
    return fail<ExtensionInfo>(make_error(CASErrorKind::Unimplemented, "Estensione algebrica non supportata (solo sqrt(n) per ora)"));
}

[[nodiscard]] BigInt lcm_positive(BigInt lhs, BigInt rhs) {
    lhs = lhs.abs();
    rhs = rhs.abs();
    if (lhs.is_zero() || rhs.is_zero()) {
        return BigInt(0);
    }
    return (lhs / gcd(lhs, rhs)) * rhs;
}

[[nodiscard]] Result<Factorization> factor_over_rationals(ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx) {
    auto normalized = expand(poly, ctx);
    if (normalized.is_error()) {
        return fail<Factorization>(normalized.error());
    }
    auto parsed = parse_polynomial(normalized.value(), var, ctx);
    if (parsed.is_error()) {
        return fail<Factorization>(parsed.error());
    }
    auto rational_poly = poly_to_rational_poly(parsed.value());
    if (rational_poly.is_error()) {
        return fail<Factorization>(make_error(
            rational_poly.error().kind,
            rational_poly.error().message + ": " + debug_print(normalized.value())));
    }

    BigInt scale(1);
    for (const Rational& coefficient : rational_poly.value().coefficients()) {
        scale = lcm_positive(scale, coefficient.denominator());
    }

    PolyExpr integer_coefficients;
    integer_coefficients.reserve(rational_poly.value().size());
    for (const Rational& coefficient : rational_poly.value().coefficients()) {
        const Rational scaled = coefficient * Rational(scale);
        if (scaled.denominator() != BigInt(1)) {
            return fail<Factorization>(make_error(
                CASErrorKind::InternalError,
                "Clearing rational denominators did not produce integer coefficients"));
        }
        integer_coefficients.push_back(
            scaled.numerator().is_zero() ? ExprPtr{} : ctx.arena().make<IntegerLit>(scaled.numerator()));
    }
    normalize_poly(integer_coefficients);
    if (integer_coefficients.empty()) {
        return fail<Factorization>(make_error(CASErrorKind::InvalidArgument, "Zero polynomial has no canonical factorization"));
    }

    auto integer_expr = polynomial_to_expr(integer_coefficients, var, ctx);
    if (integer_expr.is_error()) {
        return fail<Factorization>(integer_expr.error());
    }
    return factor_over_integers(integer_expr.value(), var, ctx);
}

[[nodiscard]] Result<ExprPtr> simplify_factor(ExprPtr factor, symbolic::CASContext& ctx) {
    auto simplified = ctx.simplify(factor);
    if (simplified.is_error()) {
        return fail<ExprPtr>(simplified.error());
    }
    return simplified;
}

[[nodiscard]] Result<std::optional<Factorization>> factor_radical_special_cases(
    ExprPtr poly,
    const Symbol& var,
    const ExtensionInfo& ext_info,
    symbolic::CASContext& ctx) {
    if (ext_info.min_poly.size() != 3U ||
        ext_info.min_poly[1] != Rational(BigInt(0)) ||
        ext_info.min_poly[2] != Rational(BigInt(1))) {
        return ok(std::optional<Factorization>{});
    }
    const Rational radicand = -ext_info.min_poly[0];

    auto parsed = parse_polynomial(poly, var, ctx);
    if (parsed.is_error()) {
        return fail<std::optional<Factorization>>(parsed.error());
    }
    auto rational_poly = poly_to_rational_poly(parsed.value());
    if (rational_poly.is_error()) {
        return ok(std::optional<Factorization>{});
    }

    auto x = ctx.arena().make<Symbol>(var.name);
    auto alpha = ext_info.original_expr;
    Factorization result{.content = make_integer(ctx.arena(), 1), .factors = {}};

    if (rational_poly.value().size() == 3U &&
        rational_poly.value()[0] == -radicand &&
        rational_poly.value()[1] == Rational(BigInt(0)) &&
        rational_poly.value()[2] == Rational(BigInt(1))) {
        auto left = simplify_factor(ctx.arena().make<Binary>(BinaryOp::Sub, x, alpha), ctx);
        if (left.is_error()) return fail<std::optional<Factorization>>(left.error());
        auto right = simplify_factor(ctx.arena().make<Binary>(BinaryOp::Add, x, alpha), ctx);
        if (right.is_error()) return fail<std::optional<Factorization>>(right.error());
        result.factors.push_back({left.value(), 1});
        result.factors.push_back({right.value(), 1});
        return ok(std::optional<Factorization>(std::move(result)));
    }

    if (rational_poly.value().size() == 5U &&
        rational_poly.value()[0] == Rational(BigInt(1)) &&
        rational_poly.value()[1] == Rational(BigInt(0)) &&
        rational_poly.value()[2] == Rational(BigInt(2)) - radicand &&
        rational_poly.value()[3] == Rational(BigInt(0)) &&
        rational_poly.value()[4] == Rational(BigInt(1))) {
        auto x_squared = ctx.arena().make<Binary>(BinaryOp::Pow, x, make_integer(ctx.arena(), 2));
        auto alpha_x = ctx.arena().make<Binary>(BinaryOp::Mul, alpha, x);
        auto one = make_integer(ctx.arena(), 1);
        auto first = simplify_factor(
            ctx.arena().make<Sum>(std::vector<ExprPtr>{
                x_squared,
                ctx.arena().make<Unary>(UnaryOp::Neg, alpha_x),
                one,
            }),
            ctx);
        if (first.is_error()) return fail<std::optional<Factorization>>(first.error());
        auto second = simplify_factor(
            ctx.arena().make<Sum>(std::vector<ExprPtr>{x_squared, alpha_x, one}),
            ctx);
        if (second.is_error()) return fail<std::optional<Factorization>>(second.error());
        result.factors.push_back({first.value(), 1});
        result.factors.push_back({second.value(), 1});
        return ok(std::optional<Factorization>(std::move(result)));
    }

    return ok(std::optional<Factorization>{});
}

[[nodiscard]] Result<PolyExpr> reduce_poly_mod_m(PolyExpr p, const Symbol& y, const RatPoly& m_y, symbolic::CASContext& ctx) {
    for (auto& coeff : p.coefficients()) {
        if (!coeff) continue;
        auto poly_y = parse_polynomial(coeff, y, ctx);
        if (poly_y.is_ok()) {
            auto rat_y = poly_to_rational_poly(poly_y.value());
            if (rat_y.is_ok()) {
                auto [q, rem] = div_rem_rational_poly(rat_y.value(), m_y);
                std::vector<ExprPtr> rem_expr_coeffs;
                for (const auto& c : rem.coefficients()) {
                    rem_expr_coeffs.push_back(ctx.arena().make<RationalLit>(c.numerator(), c.denominator()));
                }
                auto int_poly_res = poly_to_integer_poly(PolyExpr(rem_expr_coeffs));
                auto rem_expr = polynomial_to_expr(PolyExpr(integer_coefficients_to_poly(int_poly_res.is_ok() ? int_poly_res.value() : IntPoly{}, ctx.arena())), y, ctx); // Simplification bypass
                // Actually easier:
                std::vector<ExprPtr> rem_coeffs;
                for (const auto& c : rem.coefficients()) {
                    rem_coeffs.push_back(ctx.arena().make<RationalLit>(c.numerator(), c.denominator()));
                }
                auto rem_poly_expr = polynomial_to_expr(PolyExpr(rem_coeffs), y, ctx);
                if (rem_poly_expr.is_ok()) coeff = rem_poly_expr.value();
            }
        }
    }
    normalize_poly(p);
    return ok(p);
}

[[nodiscard]] Result<PolyExpr> gcd_over_extension_internal(
    PolyExpr A,
    PolyExpr B,
    [[maybe_unused]] const Symbol& x,
    const Symbol& y,
    const RatPoly& m_y,
    symbolic::CASContext& ctx) {
    
    auto reduce = [&](PolyExpr p) { return reduce_poly_mod_m(std::move(p), y, m_y, ctx); };

    A = reduce(std::move(A)).value();
    B = reduce(std::move(B)).value();

    while (!is_zero_poly(B)) {
        ExprPtr lc_b = leading_coefficient(B);
        auto poly_lc = parse_polynomial(lc_b, y, ctx);
        if (poly_lc.is_error()) return fail<PolyExpr>(poly_lc.error());
        auto rat_lc = poly_to_rational_poly(poly_lc.value());
        if (rat_lc.is_error()) return fail<PolyExpr>(rat_lc.error());

        auto [g, s, t] = extended_gcd_rational_poly(rat_lc.value(), m_y);
        if (g.degree() > 0 || g.is_zero()) {
            return fail<PolyExpr>(make_error(CASErrorKind::InvalidArgument, "Coefficiente non invertibile nell'estensione (m(y) non irriducibile?)"));
        }
        
        Rational inv_g = Rational(BigInt(1)) / g.constant_term();
        std::vector<ExprPtr> inv_coeffs;
        for (const auto& c : s.coefficients()) {
            inv_coeffs.push_back(ctx.arena().make<RationalLit>((c * inv_g).numerator(), (c * inv_g).denominator()));
        }
        auto inv_lc_expr = polynomial_to_expr(PolyExpr(inv_coeffs), y, ctx).value();

        // B_monic = B * inv_lc
        auto b_monic_res = poly_multiply(B, PolyExpr({inv_lc_expr}), ctx);
        if (b_monic_res.is_error()) return fail<PolyExpr>(b_monic_res.error());
        PolyExpr B_monic = reduce(std::move(b_monic_res.value())).value();

        auto div_res = divide_poly_with_remainder(A, B_monic, ctx);
        if (div_res.is_error()) return fail<PolyExpr>(div_res.error());
        
        A = std::move(B);
        B = reduce(std::move(div_res.value().remainder)).value();
    }

    return normalize_poly_monic(A, ctx);
}

} // namespace

Result<Factorization> factor_polynomial(
    ExprPtr poly,
    const Symbol& var,
    symbolic::CASContext& ctx,
    std::optional<ExprPtr> extension) {
    
    if (!extension.has_value()) {
        return factor_over_integers(poly, var, ctx);
    }

    auto ext_info_res = get_extension_info(extension.value(), ctx);
    if (ext_info_res.is_error()) return fail<Factorization>(ext_info_res.error());
    auto ext_info = ext_info_res.value();

    auto special = factor_radical_special_cases(poly, var, ext_info, ctx);
    if (special.is_error()) return fail<Factorization>(special.error());
    if (special.value().has_value()) return ok(std::move(*special.value()));

    // 1. Preparazione f(x, y)
    auto alpha = ctx.arena().make<Symbol>(ext_info.var.name);
    auto f_x_y_raw = replace_subexpression(poly, ext_info.original_expr, alpha, ctx.arena());
    auto f_x_y_expr = ctx.simplify(f_x_y_raw);
    if (f_x_y_expr.is_error()) return fail<Factorization>(f_x_y_expr.error());
    
    auto f_poly_res = parse_polynomial(f_x_y_expr.value(), var, ctx);
    if (f_poly_res.is_error()) return fail<Factorization>(f_poly_res.error());
    PolyExpr f_poly = f_poly_res.value();

    // 2. Norm(f) = Res_y(m(y), f(x, y))
    std::vector<ExprPtr> m_coeffs;
    for (const auto& c : ext_info.min_poly.coefficients()) {
        m_coeffs.push_back(ctx.arena().make<RationalLit>(c.numerator(), c.denominator()));
    }
    auto m_y_expr = polynomial_to_expr(PolyExpr(m_coeffs), ext_info.var, ctx).value();
    
    int s = 0;
    while (s < 10) {
        // g(x) = Norm(f(x - sy, y))
        ExprPtr f_shifted_expr;
        if (s == 0) {
            f_shifted_expr = f_x_y_expr.value();
        } else {
            auto shift = ctx.arena().make<Binary>(BinaryOp::Mul, make_integer(ctx.arena(), s), ctx.arena().make<Symbol>(ext_info.var.name));
            auto x_minus_sy = ctx.arena().make<Binary>(BinaryOp::Sub, ctx.arena().make<Symbol>(var.name), shift);
            f_shifted_expr = ctx.substitute(f_x_y_expr.value(), var, x_minus_sy).value();
        }
        auto expanded_shift = expand(f_shifted_expr, ctx);
        if (expanded_shift.is_error()) return fail<Factorization>(expanded_shift.error());
        f_shifted_expr = expanded_shift.value();

        auto norm_res = polynomial_resultant(m_y_expr, f_shifted_expr, ext_info.var, ctx);
        if (norm_res.is_error()) return fail<Factorization>(norm_res.error());
        ExprPtr g_x = norm_res.value();

        // 3. Fattorizza g(x) in Q[x]
        auto g_factors_res = factor_over_rationals(g_x, var, ctx);
        if (g_factors_res.is_error()) return fail<Factorization>(g_factors_res.error());
        auto g_factors = g_factors_res.value();

        // Verifica se Norm è square-free (per semplicità, controlliamo se ci sono fattori con molteplicità > 1)
        bool square_free = true;
        for (const auto& fact : g_factors.factors) {
            if (fact.multiplicity > 1) {
                square_free = false;
                break;
            }
        }

        if (square_free || s > 5) {
            // 4. GCD: f_i(x, y) = gcd(f(x, y), g_i(x + sy, y))
            Factorization result;
            result.content = make_integer(ctx.arena(), 1); // Trager usually returns monic factors
            
            for (const auto& gi_fact : g_factors.factors) {
                ExprPtr gi_shifted_expr = gi_fact.factor;
                if (s != 0) {
                    auto shift = ctx.arena().make<Binary>(BinaryOp::Mul, make_integer(ctx.arena(), s), ctx.arena().make<Symbol>(ext_info.var.name));
                    auto x_plus_sy = ctx.arena().make<Binary>(BinaryOp::Add, ctx.arena().make<Symbol>(var.name), shift);
                    gi_shifted_expr = ctx.substitute(gi_fact.factor, var, x_plus_sy).value();
                }

                auto gi_poly_res = parse_polynomial(gi_shifted_expr, var, ctx);
                if (gi_poly_res.is_error()) continue;

                auto gcd_res = gcd_over_extension_internal(f_poly, gi_poly_res.value(), var, ext_info.var, ext_info.min_poly, ctx);
                if (gcd_res.is_ok() && poly_degree(gcd_res.value()) > 0) {
                    auto factor_expr = polynomial_to_expr(gcd_res.value(), var, ctx).value();
                    // Replace __alpha back with original expression
                    factor_expr = ctx.substitute(factor_expr, ext_info.var, ext_info.original_expr).value();
                    result.factors.push_back({factor_expr, 1});
                }
            }
            return ok(result);
        }
        s++;
    }

    return fail<Factorization>(make_error(CASErrorKind::InternalError, "Impossibile trovare uno shift per la fattorizzazione di Trager"));
}


} // namespace cas::algebra
