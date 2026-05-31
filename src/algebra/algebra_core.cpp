#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/error_helpers.hpp"
#include "algebra_internal.hpp"
#include <algorithm>
#include <chrono>
#include <map>
#include <utility>
#include <vector>

namespace cas::algebra {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message, std::optional<std::string> hint) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::move(hint),
    };
}

[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator.is_zero();
    }
    return false;
}

[[nodiscard]] bool is_one_expr(ExprPtr expr) {
    static const BigInt one(1);
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == one;
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == one && rational->denominator == one;
    }
    return false;
}

[[nodiscard]] bool contains_decimal_literal(ExprPtr expr) {
    if (!expr) {
        return false;
    }
    if (expr_is<DecimalLit>(expr)) {
        return true;
    }

    return visit_expr(expr, [](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (
            std::is_same_v<Node, IntegerLit> ||
            std::is_same_v<Node, RationalLit> ||
            std::is_same_v<Node, Symbol> ||
            std::is_same_v<Node, Constant>) {
            return false;
        } else if constexpr (std::is_same_v<Node, Unary>) {
            return contains_decimal_literal(node.operand);
        } else if constexpr (std::is_same_v<Node, Binary>) {
            return contains_decimal_literal(node.left) || contains_decimal_literal(node.right);
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            return std::any_of(node.args.begin(), node.args.end(), [](ExprPtr arg) { return contains_decimal_literal(arg); });
        } else if constexpr (std::is_same_v<Node, Sum>) {
            return std::any_of(node.terms.begin(), node.terms.end(), [](ExprPtr term) { return contains_decimal_literal(term); });
        } else if constexpr (std::is_same_v<Node, Product>) {
            return std::any_of(node.factors.begin(), node.factors.end(), [](ExprPtr factor) { return contains_decimal_literal(factor); });
        } else if constexpr (std::is_same_v<Node, Integral>) {
            return contains_decimal_literal(node.integrand) ||
                   (node.lower.has_value() && contains_decimal_literal(*node.lower)) ||
                   (node.upper.has_value() && contains_decimal_literal(*node.upper));
        } else if constexpr (std::is_same_v<Node, Derivative>) {
            return contains_decimal_literal(node.expression);
        } else if constexpr (std::is_same_v<Node, Limit>) {
            return contains_decimal_literal(node.expression) || contains_decimal_literal(node.point);
        } else if constexpr (std::is_same_v<Node, RootOf>) {
            return contains_decimal_literal(node.polynomial);
        } else if constexpr (std::is_same_v<Node, Matrix>) {
            return std::any_of(node.elements.begin(), node.elements.end(), [](ExprPtr element) {
                return contains_decimal_literal(element);
            });
        } else {
            return false;
        }
    });
}

[[nodiscard]] Result<ExprPtr> simplify_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    auto simplified = ctx.simplify(expr);
    if (simplified.is_error()) {
        return fail<ExprPtr>(simplified.error());
    }
    return simplified;
}

[[nodiscard]] Result<ExprPtr> clone_into_context(ExprPtr expr, symbolic::CASContext& ctx) {
    auto cloned = symbolic::materialize_expr(expr, ctx.arena());
    if (cloned.is_error()) {
        return fail<ExprPtr>(cloned.error());
    }
    return simplify_expr(cloned.value(), ctx);
}

[[nodiscard]] Result<ExprPtr> add_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return simplify_expr(ctx.arena().make<Sum>(std::vector<ExprPtr>{lhs, rhs}), ctx);
}

[[nodiscard]] Result<ExprPtr> negate_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    return simplify_expr(ctx.arena().make<Unary>(UnaryOp::Neg, expr), ctx);
}

[[nodiscard]] Result<ExprPtr> multiply_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return simplify_expr(ctx.arena().make<Product>(std::vector<ExprPtr>{lhs, rhs}), ctx);
}

[[nodiscard]] Result<ExprPtr> divide_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    return simplify_expr(ctx.arena().make<Binary>(BinaryOp::Div, lhs, rhs), ctx);
}

[[nodiscard]] Result<ExprPtr> pow_expr(ExprPtr base, std::size_t exponent, symbolic::CASContext& ctx) {
    if (exponent == 0U) {
        return ok(make_integer(ctx.arena(), 1));
    }
    if (exponent == 1U) {
        return ok(base);
    }
    return simplify_expr(
        ctx.arena().make<Binary>(BinaryOp::Pow, base, make_integer(ctx.arena(), static_cast<long long>(exponent))),
        ctx);
}

[[nodiscard]] Result<ExprPtr> subtract_exprs(ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    auto negated = negate_expr(rhs, ctx);
    if (negated.is_error()) {
        return fail<ExprPtr>(negated.error());
    }
    return add_exprs(lhs, negated.value(), ctx);
}

[[nodiscard]] BigInt pow_bigint_nonnegative(BigInt base, unsigned int exponent) {
    BigInt result(1);
    while (exponent > 0U) {
        if ((exponent & 1U) != 0U) {
            result *= base;
        }
        exponent >>= 1U;
        if (exponent > 0U) {
            base *= base;
        }
    }
    return result;
}

using MultivariateExponentKey = std::vector<std::pair<std::string, unsigned int>>;

[[nodiscard]] MultivariateExponentKey make_multivariate_key(
    const std::vector<std::pair<Symbol, unsigned int>>& factors) {
    MultivariateExponentKey key;
    key.reserve(factors.size());
    for (const auto& [symbol, exponent] : factors) {
        if (exponent == 0U) {
            continue;
        }
        key.emplace_back(symbol.name, exponent);
    }

    std::sort(key.begin(), key.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });

    MultivariateExponentKey canonical;
    canonical.reserve(key.size());
    for (const auto& [name, exponent] : key) {
        if (!canonical.empty() && canonical.back().first == name) {
            canonical.back().second += exponent;
        } else {
            canonical.emplace_back(name, exponent);
        }
    }
    return canonical;
}

[[nodiscard]] std::vector<std::pair<Symbol, unsigned int>> make_multivariate_exponents(
    const MultivariateExponentKey& key) {
    std::vector<std::pair<Symbol, unsigned int>> exponents;
    exponents.reserve(key.size());
    for (const auto& [name, exponent] : key) {
        exponents.emplace_back(Symbol(name), exponent);
    }
    return exponents;
}

[[nodiscard]] Result<BigInt> expr_to_integer_value_for_multivariate(ExprPtr value) {
    if (const auto* integer = expr_cast<IntegerLit>(value)) {
        return ok(integer->value);
    }
    return fail<BigInt>(make_error(CASErrorKind::InvalidArgument, "Valore non intero nella valutazione multivariata", std::nullopt));
}

[[nodiscard]] Result<ExprPtr> build_multivariate_monomial_expr(const MultivariateTerm& term, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> products;
    if (term.coefficient != BigInt(1)) {
        products.push_back(ctx.arena().make<IntegerLit>(term.coefficient));
    } else if (term.factors.empty()) {
        products.push_back(ctx.arena().make<IntegerLit>(term.coefficient));
    }
    
    for (const auto& [symbol, exponent] : term.factors) {
        if (exponent == 0) continue;
        ExprPtr var_expr = ctx.arena().make<Symbol>(symbol.name);
        if (exponent == 1) {
            products.push_back(var_expr);
        } else {
            products.push_back(ctx.arena().make<Binary>(
                BinaryOp::Pow,
                var_expr,
                ctx.arena().make<IntegerLit>(BigInt(static_cast<long long>(exponent)))
            ));
        }
    }
    
    if (products.empty()) {
        return ok(make_integer(ctx.arena(), 1));
    }
    if (products.size() == 1) {
        return ok(products[0]);
    }
    return ok(ctx.arena().make<Product>(std::move(products)));
}

[[nodiscard]] std::vector<MultivariateTerm> normalize_multivariate_terms(
    const std::vector<MultivariateTerm>& terms) {
    std::map<MultivariateExponentKey, BigInt> accumulated;
    for (const MultivariateTerm& term : terms) {
        if (term.coefficient.is_zero()) {
            continue;
        }
        accumulated[make_multivariate_key(term.factors)] += term.coefficient;
    }

    std::vector<MultivariateTerm> normalized;
    normalized.reserve(accumulated.size());
    for (const auto& [key, coefficient] : accumulated) {
        if (coefficient.is_zero()) {
            continue;
        }
        normalized.push_back(MultivariateTerm{
            .coefficient = coefficient,
            .factors = make_multivariate_exponents(key),
        });
    }
    return normalized;
}

MultivariatePolynomial::MultivariatePolynomial(std::vector<MultivariateTerm> terms)
    : terms_(normalize_multivariate_terms(terms)) {}

bool MultivariatePolynomial::is_zero() const noexcept {
    return terms_.empty();
}

std::size_t MultivariatePolynomial::total_degree() const noexcept {
    std::size_t degree = 0U;
    for (const MultivariateTerm& term : terms_) {
        std::size_t term_degree = 0U;
        for (const auto& [symbol, exponent] : term.factors) {
            static_cast<void>(symbol);
            term_degree += exponent;
        }
        degree = std::max(degree, term_degree);
    }
    return degree;
}

std::vector<Symbol> MultivariatePolynomial::variables() const {
    std::vector<std::string> names;
    for (const MultivariateTerm& term : terms_) {
        for (const auto& [symbol, exponent] : term.factors) {
            if (exponent == 0U) {
                continue;
            }
            names.push_back(symbol.name);
        }
    }

    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());

    std::vector<Symbol> result;
    result.reserve(names.size());
    for (const std::string& name : names) {
        result.emplace_back(name);
    }
    return result;
}

const std::vector<MultivariateTerm>& MultivariatePolynomial::terms() const noexcept {
    return terms_;
}

Result<MultivariatePolynomial> MultivariatePolynomial::evaluate_at(const Symbol& var, ExprPtr value) const {
    std::vector<MultivariateTerm> evaluated_terms;
    evaluated_terms.reserve(terms_.size());

    // Se value non è un intero, per ora falliamo come da specifiche attuali del test (integers only)
    const auto* integer_val = expr_cast<IntegerLit>(value);
    if (!integer_val) {
        // F0.8-MIGRATED
        return make_unimplemented<MultivariatePolynomial>(
            "algebra", "MultivariatePolynomial::evaluate_at",
            "non-IntegerLit value",
            cas::error::reason_codes::ALGEBRA_MULTIVAR_NON_INTEGER,
            "Extend evaluate_at to accept RationalLit and RootOf values",
            "F0.8");
    }

    for (const auto& term : terms_) {
        BigInt new_coeff = term.coefficient;
        std::vector<std::pair<Symbol, unsigned int>> new_factors;
        
        for (const auto& factor : term.factors) {
            const Symbol& symbol = factor.first;
            const unsigned int exponent = factor.second;
            if (symbol.name == var.name) {
                // Valuta: coeff *= value^exponent
                for (std::size_t i = 0; i < static_cast<std::size_t>(exponent); ++i) {
                    new_coeff = new_coeff * integer_val->value;
                }
            } else {
                new_factors.push_back({symbol, exponent});
            }
        }
        
        evaluated_terms.push_back(MultivariateTerm{
            .coefficient = new_coeff,
            .factors = std::move(new_factors)
        });
    }

    return ok(MultivariatePolynomial(std::move(evaluated_terms)));
}

Result<ExprPtr> MultivariatePolynomial::evaluate_at_rational(
    const Symbol& var, const Rational& value, AstArena& arena) const {
    Rational result(BigInt(0));
    for (const auto& term : terms_) {
        bool has_other_vars = false;
        for (const auto& factor : term.factors) {
            if (factor.first.name != var.name) { has_other_vars = true; break; }
        }
        if (has_other_vars) {
            // F0.8-MIGRATED
            return make_unimplemented<ExprPtr>(
                "algebra", "MultivariatePolynomial::evaluate_at_rational",
                "term with remaining free variables after substitution",
                cas::error::reason_codes::ALGEBRA_MULTIVAR_REMAINING_VARS,
                "Add recursive substitution for all free variables before calling evaluate_at_rational",
                "F0.8");
        }
        unsigned int exponent = 0;
        for (const auto& factor : term.factors) {
            if (factor.first.name == var.name) { exponent = factor.second; break; }
        }
        Rational contrib(term.coefficient);
        for (unsigned int i = 0; i < exponent; ++i) contrib = contrib * value;
        result = result + contrib;
    }
    if (result.denominator() == BigInt(1)) {
        return ok(static_cast<ExprPtr>(arena.make<IntegerLit>(result.numerator())));
    }
    return ok(static_cast<ExprPtr>(arena.make<RationalLit>(result.numerator(), result.denominator())));
}

Result<std::vector<ExprPtr>> MultivariatePolynomial::to_univariate_coefficients(const Symbol& var, symbolic::CASContext& ctx) const {
    // Raggruppa i termini per grado della variabile target
    std::map<std::size_t, std::vector<MultivariateTerm>> groups;
    std::size_t max_deg = 0;

    for (const auto& term : terms_) {
        std::size_t deg = 0;
        for (const auto& v : term.factors) {
            if (v.first.name == var.name) {
                deg = v.second;
                break;
            }
        }
        groups[deg].push_back(term);
        if (deg > max_deg) max_deg = deg;
    }

    std::vector<ExprPtr> result(max_deg + 1);

    for (const auto& [deg, terms] : groups) {
        // Costruisci il coefficiente come somma dei termini (senza la variabile var)
        std::vector<ExprPtr> coeff_parts;
        for (const auto& term : terms) {
            MultivariateTerm reduced = term;
            // Rimuovi var dai fattori
            for (auto it = reduced.factors.begin(); it != reduced.factors.end(); ++it) {
                if (it->first.name == var.name) {
                    reduced.factors.erase(it);
                    break;
                }
            }
            auto part = build_multivariate_monomial_expr(reduced, ctx);
            if (part.is_error()) return fail<std::vector<ExprPtr>>(part.error());
            coeff_parts.push_back(part.value());
        }
        
        if (coeff_parts.empty()) {
            result[deg] = make_integer(ctx.arena(), 0);
        } else if (coeff_parts.size() == 1) {
            result[deg] = coeff_parts[0];
        } else {
            result[deg] = ctx.arena().make<Sum>(std::move(coeff_parts));
        }
    }

    // Assicura che gli slot vuoti siano inizializzati a zero
    for (auto& coeff : result) {
        if (!coeff) {
            coeff = make_integer(ctx.arena(), 0);
        }
    }

    return ok(std::move(result));
}

MultivariatePolynomial MultivariatePolynomial::operator+(const MultivariatePolynomial& other) const {
    std::vector<MultivariateTerm> combined;
    combined.reserve(terms_.size() + other.terms_.size());
    combined.insert(combined.end(), terms_.begin(), terms_.end());
    combined.insert(combined.end(), other.terms_.begin(), other.terms_.end());
    return MultivariatePolynomial(std::move(combined));
}

MultivariatePolynomial MultivariatePolynomial::operator*(const MultivariatePolynomial& other) const {
    if (is_zero() || other.is_zero()) {
        return MultivariatePolynomial{};
    }

    std::vector<MultivariateTerm> products;
    products.reserve(terms_.size() * other.terms_.size());
    for (const MultivariateTerm& lhs : terms_) {
        for (const MultivariateTerm& rhs : other.terms_) {
            std::vector<std::pair<Symbol, unsigned int>> factors = lhs.factors;
            factors.insert(factors.end(), rhs.factors.begin(), rhs.factors.end());
            products.push_back(MultivariateTerm{
                .coefficient = lhs.coefficient * rhs.coefficient,
                .factors = std::move(factors),
            });
        }
    }
    return MultivariatePolynomial(std::move(products));
}

} // namespace cas::algebra
