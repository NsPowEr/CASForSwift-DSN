#include "integrate_definite_patterns.hpp"

#include "cas/builtin_functions.hpp"
#include "cas/error.hpp"

#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

// Flatten a multiplicative product into a list of factors.  Handles Product,
// Binary(Mul), and treats any other node as a single factor.
void flatten_mul_factors(ExprPtr expr, std::vector<ExprPtr>& out) {
    if (!expr) return;
    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr f : product->factors) flatten_mul_factors(f, out);
        return;
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Mul) {
            flatten_mul_factors(binary->left, out);
            flatten_mul_factors(binary->right, out);
            return;
        }
    }
    out.push_back(expr);
}

[[nodiscard]] bool is_zero_lit(ExprPtr expr) {
    const auto* lit = expr_cast<IntegerLit>(expr);
    return lit && lit->value.is_zero();
}

[[nodiscard]] bool is_var_symbol(ExprPtr expr, const Symbol& var) {
    const auto* sym = expr_cast<Symbol>(expr);
    return sym && sym->name == var.name;
}

[[nodiscard]] bool is_bessel_j(ExprPtr expr) {
    const auto* call = expr_cast<FuncCall>(expr);
    return call && call->func_id == BuiltinOp::BesselJ && call->args.size() == 2U;
}

[[nodiscard]] bool is_bessel_zero(ExprPtr expr) {
    const auto* call = expr_cast<FuncCall>(expr);
    return call && call->func_id == BuiltinOp::BesselZero && call->args.size() == 2U;
}

// Test whether `expr` depends on `var` (structural search for the symbol).
[[nodiscard]] bool depends_on_var(ExprPtr expr, const Symbol& var) {
    if (!expr) return false;
    if (const auto* sym = expr_cast<Symbol>(expr)) return sym->name == var.name;
    bool found = false;
    visit_expr(expr, [&](const auto& node) {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, Unary>) {
            if (depends_on_var(node.operand, var)) found = true;
        } else if constexpr (std::is_same_v<Node, Binary>) {
            if (depends_on_var(node.left, var)) found = true;
            if (depends_on_var(node.right, var)) found = true;
        } else if constexpr (std::is_same_v<Node, FuncCall>) {
            for (ExprPtr arg : node.args)
                if (depends_on_var(arg, var)) { found = true; break; }
        } else if constexpr (std::is_same_v<Node, Sum>) {
            for (ExprPtr term : node.terms)
                if (depends_on_var(term, var)) { found = true; break; }
        } else if constexpr (std::is_same_v<Node, Product>) {
            for (ExprPtr fac : node.factors)
                if (depends_on_var(fac, var)) { found = true; break; }
        }
    });
    return found;
}

// Try to extract α from `arg = α·var` where α is independent of var.
// Returns nullopt if not of that form (e.g. var-quadratic, var-free, x+const).
[[nodiscard]] std::optional<ExprPtr> extract_linear_coefficient(
    ExprPtr arg,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    if (!arg) return std::nullopt;
    if (is_var_symbol(arg, var)) {
        return ok(ExprPtr{ctx.arena().make<IntegerLit>(BigInt(1))}).value();
    }
    // arg / var must be independent of var, and arg must be linear in var
    // (no constant term).  Use diff and substitute at 0 cheaply:
    //   if arg is α·var (linear, zero constant term), then α = arg / var
    //   and arg|_{var=0} = 0.
    // Detect this via structural Mul.
    std::vector<ExprPtr> factors;
    flatten_mul_factors(arg, factors);
    std::vector<ExprPtr> coef_factors;
    bool found_var = false;
    for (ExprPtr f : factors) {
        if (is_var_symbol(f, var) && !found_var) {
            found_var = true;
            continue;
        }
        if (depends_on_var(f, var)) return std::nullopt;
        coef_factors.push_back(f);
    }
    if (!found_var) return std::nullopt;
    AstArena& arena = ctx.arena();
    if (coef_factors.empty()) return ExprPtr{arena.make<IntegerLit>(BigInt(1))};
    if (coef_factors.size() == 1U) return coef_factors.front();
    return ExprPtr{arena.make<Product>(std::move(coef_factors))};
}

// Decide whether two expressions are mathematically equal: difference simplifies to zero.
[[nodiscard]] Result<bool> mathematically_equal_simple(
    ExprPtr lhs,
    ExprPtr rhs,
    symbolic::CASContext& ctx) {
    if (structural_equal(lhs, rhs)) return ok(true);
    ExprPtr diff = ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs);
    auto simplified = ctx.simplify(diff);
    if (simplified.is_error()) return fail<bool>(simplified.error());
    const auto* lit = expr_cast<IntegerLit>(simplified.value());
    if (lit && lit->value.is_zero()) return ok(true);
    return ok(false);
}

// Build BesselJ(order+1, BesselZero(order, m)).
[[nodiscard]] ExprPtr build_j_nu_plus_one_at_zero(
    ExprPtr order,
    ExprPtr m_index,
    AstArena& arena) {
    ExprPtr one = arena.make<IntegerLit>(BigInt(1));
    ExprPtr order_plus_one = arena.make<Binary>(BinaryOp::Add, order, one);
    ExprPtr zero_node = arena.make<FuncCall>(BuiltinOp::BesselZero, std::vector<ExprPtr>{order, m_index});
    return arena.make<FuncCall>(BuiltinOp::BesselJ, std::vector<ExprPtr>{order_plus_one, zero_node});
}

}  // namespace

[[nodiscard]] Result<std::optional<ExprPtr>> pattern_bessel_orthogonality(const DefiniteContext& dc) {
    // Limits: lower must be exactly 0 (literal), upper is the scale parameter a.
    if (!is_zero_lit(dc.lower)) return ok(std::optional<ExprPtr>{});
    ExprPtr a_param = dc.upper;
    if (!a_param) return ok(std::optional<ExprPtr>{});
    if (depends_on_var(a_param, dc.var)) return ok(std::optional<ExprPtr>{});

    // Flatten factors of normalized integrand: must be exactly { x, BesselJ(ν, α·x), BesselJ(ν, β·x) }
    // possibly with extra var-independent constant factors that all multiply through.
    std::vector<ExprPtr> factors;
    flatten_mul_factors(dc.integrand_normalized, factors);

    bool found_x = false;
    std::vector<ExprPtr> bessel_terms;       // BesselJ FuncCall nodes
    std::vector<ExprPtr> extra_const_factors; // constant in var
    for (ExprPtr f : factors) {
        if (is_var_symbol(f, dc.var) && !found_x) {
            found_x = true;
            continue;
        }
        if (is_bessel_j(f)) {
            bessel_terms.push_back(f);
            continue;
        }
        // Handle BesselJ(...)^n: replicate the base node n times in bessel_terms.
        if (const auto* binary = expr_cast<Binary>(f)) {
            if (binary->op == BinaryOp::Pow && is_bessel_j(binary->left)) {
                const auto* exp_lit = expr_cast<IntegerLit>(binary->right);
                if (exp_lit && !exp_lit->value.is_negative() && !exp_lit->value.is_zero()
                    && exp_lit->value.bit_length() <= 8U) {
                    const std::size_t count = static_cast<std::size_t>(exp_lit->value.to_u64());
                    for (std::size_t i = 0; i < count; ++i) bessel_terms.push_back(binary->left);
                    continue;
                }
            }
        }
        if (!depends_on_var(f, dc.var)) {
            extra_const_factors.push_back(f);
            continue;
        }
        // Some other var-dependent factor: not orthogonality pattern.
        return ok(std::optional<ExprPtr>{});
    }
    std::fprintf(stderr, "[BO] found_x=%d bessel_terms=%zu extras=%zu total_factors=%zu\n",
                 (int)found_x, bessel_terms.size(), extra_const_factors.size(), factors.size());
    if (!found_x) return ok(std::optional<ExprPtr>{});
    if (bessel_terms.size() != 2U) return ok(std::optional<ExprPtr>{});

    // Extract ν, α from each BesselJ(ν, arg).
    const auto& call0 = expr_ref<FuncCall>(bessel_terms[0]);
    const auto& call1 = expr_ref<FuncCall>(bessel_terms[1]);
    ExprPtr nu_0 = call0.args[0];
    ExprPtr nu_1 = call1.args[0];

    // ν must match between the two factors.
    auto nu_equal = mathematically_equal_simple(nu_0, nu_1, dc.ctx);
    if (nu_equal.is_error()) return fail<std::optional<ExprPtr>>(nu_equal.error());
    if (!nu_equal.value()) return ok(std::optional<ExprPtr>{});

    if (depends_on_var(nu_0, dc.var)) return ok(std::optional<ExprPtr>{});

    auto alpha_opt = extract_linear_coefficient(call0.args[1], dc.var, dc.ctx);
    if (!alpha_opt.has_value()) return ok(std::optional<ExprPtr>{});
    auto beta_opt = extract_linear_coefficient(call1.args[1], dc.var, dc.ctx);
    if (!beta_opt.has_value()) return ok(std::optional<ExprPtr>{});

    AstArena& arena = dc.ctx.arena();

    // Each α_i must equal BesselZero(ν, m_i) / a.  Equivalently, α_i · a must
    // simplify to a BesselZero(ν, m_i) node.
    auto resolve_index = [&](ExprPtr alpha) -> Result<std::optional<ExprPtr>> {
        ExprPtr product = arena.make<Binary>(BinaryOp::Mul, alpha, a_param);
        auto simplified = dc.ctx.simplify(product);
        if (simplified.is_error()) return fail<std::optional<ExprPtr>>(simplified.error());
        if (!is_bessel_zero(simplified.value())) return ok(std::optional<ExprPtr>{});
        const auto& call = expr_ref<FuncCall>(simplified.value());
        auto same_nu = mathematically_equal_simple(call.args[0], nu_0, dc.ctx);
        if (same_nu.is_error()) return fail<std::optional<ExprPtr>>(same_nu.error());
        if (!same_nu.value()) return ok(std::optional<ExprPtr>{});
        return ok(std::optional<ExprPtr>(call.args[1]));
    };
    auto m_opt = resolve_index(alpha_opt.value());
    if (m_opt.is_error()) return fail<std::optional<ExprPtr>>(m_opt.error());
    if (!m_opt.value().has_value()) return ok(std::optional<ExprPtr>{});
    auto k_opt = resolve_index(beta_opt.value());
    if (k_opt.is_error()) return fail<std::optional<ExprPtr>>(k_opt.error());
    if (!k_opt.value().has_value()) return ok(std::optional<ExprPtr>{});

    // Decide same root index: integers compare structurally; symbolic compare via simplify.
    auto same_index = mathematically_equal_simple(m_opt.value().value(), k_opt.value().value(), dc.ctx);
    if (same_index.is_error()) return fail<std::optional<ExprPtr>>(same_index.error());

    // Assemble the constant factor multiplier (from extra_const_factors).
    auto mul_extras = [&](ExprPtr base) -> ExprPtr {
        for (ExprPtr c : extra_const_factors) {
            base = arena.make<Binary>(BinaryOp::Mul, base, c);
        }
        return base;
    };

    if (!same_index.value()) {
        ExprPtr zero = arena.make<IntegerLit>(BigInt(0));
        auto simplified = dc.ctx.simplify(zero);
        if (simplified.is_error()) return fail<std::optional<ExprPtr>>(simplified.error());
        return ok(std::optional<ExprPtr>(simplified.value()));
    }

    // m == k: result = (a^2 / 2) · J_{ν+1}(BesselZero(ν, m))^2 · extras.
    ExprPtr a_sq = arena.make<Binary>(BinaryOp::Pow, a_param, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr half_a_sq = arena.make<Binary>(BinaryOp::Div, a_sq, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr j_value = build_j_nu_plus_one_at_zero(nu_0, m_opt.value().value(), arena);
    ExprPtr j_sq = arena.make<Binary>(BinaryOp::Pow, j_value, arena.make<IntegerLit>(BigInt(2)));
    ExprPtr result = arena.make<Binary>(BinaryOp::Mul, half_a_sq, j_sq);
    result = mul_extras(result);
    auto simplified = dc.ctx.simplify(result);
    if (simplified.is_error()) return fail<std::optional<ExprPtr>>(simplified.error());
    return ok(std::optional<ExprPtr>(simplified.value()));
}

}  // namespace cas::calculus
