// differential_field_generators.cpp — mapping between the original
// log/exp expression form and the tower generator symbols t_1..t_n.
// Split out of differential_field.cpp (anti-monolith, A27 cleanup): keeps
// build()/add_extension()/derive() together in the base file and the
// generator-representation helpers (used by the Risch integration pipeline
// to keep K = C(x, t_1..t_n) differentially closed) here.

#include "cas/differential_algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/error_helpers.hpp"

#include <string>
#include <type_traits>
#include <vector>

namespace cas::calculus {

namespace {

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

Result<ExprPtr> DifferentialField::to_field_generators(ExprPtr expr, symbolic::CASContext& ctx) const {
    ExprPtr current = expr;
    AstArena& arena = ctx.arena();

    // Substitute OUTERMOST generators first (reverse tower order).  A higher
    // generator's pattern is a function whose argument still contains the lower
    // generators' original forms (e.g. t_1 = exp(exp(x)) has pattern exp(exp(x)),
    // which embeds t_0 = exp(x)).  Replacing t_0 first would rewrite the inner
    // exp(x) → t_0 *inside* that pattern, turning exp(exp(x)) into exp(t_0) so the
    // t_1 pattern no longer matches — leaving the outer transcendental unmapped
    // (silent-wrong: it is then mis-integrated as a constant coefficient).
    // from_field_generators already unwinds in this same reverse order.
    for (auto it = extensions_.rbegin(); it != extensions_.rend(); ++it) {
        std::string func_name = (it->type == ExtensionType::Logarithmic) ? "ln" : "exp";
        ExprPtr pattern = arena.make<FuncCall>(func_name, std::vector<ExprPtr>{it->argument});
        ExprPtr replacement = arena.make<Symbol>(it->t_var.name);

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

Result<ExprPtr> DifferentialField::derive_in_generators(ExprPtr expr, symbolic::CASContext& ctx) const {
    auto d = derive(expr, ctx);
    if (d.is_error()) return d;
    auto gen = to_field_generators(d.value(), ctx);
    if (gen.is_error()) return gen;
    if (auto s = ctx.simplify(gen.value()); s.is_ok()) return s;
    return gen;
}

} // namespace cas::calculus
