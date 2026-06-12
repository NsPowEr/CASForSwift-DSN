/// @file ast_clone.cpp
/// @brief Deep-copy of AST trees across AstArena instances.
///
/// Extracted from src/ast/ast.cpp (F8.0 / Task 1.2).
/// Contains: clone_into_arena.
#include "cas/ast.hpp"

#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace cas {

ExprPtr clone_into_arena(ExprPtr expr, AstArena& target, std::unordered_map<ExprPtr, ExprPtr>& cache) {
    if (!expr) {
        return expr;
    }

    if (auto it = cache.find(expr); it != cache.end()) {
        return it->second;
    }

    ExprPtr cloned;
    switch (expr->kind) {
    case ExprKind::IntegerLit:
        cloned = target.make<IntegerLit>(expr_ref<IntegerLit>(expr).value);
        break;
    case ExprKind::RationalLit: {
        const auto& node = expr_ref<RationalLit>(expr);
        cloned = target.make<RationalLit>(node.numerator, node.denominator);
        break;
    }
    case ExprKind::ComplexLit: {
        const auto& node = expr_ref<ComplexLit>(expr);
        cloned = target.make<ComplexLit>(node.re_num, node.re_den, node.im_num, node.im_den);
        break;
    }
    case ExprKind::DecimalLit:
        cloned = target.make<DecimalLit>(expr_ref<DecimalLit>(expr).text);
        break;
    case ExprKind::Symbol:
        cloned = target.make<Symbol>(expr_ref<Symbol>(expr).name);
        break;
    case ExprKind::Constant:
        cloned = target.make<Constant>(expr_ref<Constant>(expr).value);
        break;
    case ExprKind::Unary: {
        const auto& node = expr_ref<Unary>(expr);
        cloned = target.make<Unary>(node.op, clone_into_arena(node.operand, target, cache));
        break;
    }
    case ExprKind::Binary: {
        const auto& node = expr_ref<Binary>(expr);
        auto left  = clone_into_arena(node.left,  target, cache);
        auto right = clone_into_arena(node.right, target, cache);
        cloned = target.make<Binary>(node.op, left, right);
        break;
    }
    case ExprKind::FuncCall: {
        const auto& node = expr_ref<FuncCall>(expr);
        std::vector<ExprPtr> args;
        args.reserve(node.args.size());
        for (auto arg : node.args) {
            args.push_back(clone_into_arena(arg, target, cache));
        }
        cloned = target.make<FuncCall>(node.name, std::move(args));
        break;
    }
    case ExprKind::Sum: {
        const auto& node = expr_ref<Sum>(expr);
        std::vector<ExprPtr> terms;
        terms.reserve(node.terms.size());
        for (auto term : node.terms) {
            terms.push_back(clone_into_arena(term, target, cache));
        }
        cloned = target.make<Sum>(std::move(terms));
        break;
    }
    case ExprKind::Product: {
        const auto& node = expr_ref<Product>(expr);
        std::vector<ExprPtr> factors;
        factors.reserve(node.factors.size());
        for (auto factor : node.factors) {
            factors.push_back(clone_into_arena(factor, target, cache));
        }
        cloned = target.make<Product>(std::move(factors));
        break;
    }
    case ExprKind::Integral: {
        const auto& node = expr_ref<Integral>(expr);
        auto integrand = clone_into_arena(node.integrand, target, cache);
        auto lower = node.lower ? std::optional<ExprPtr>(clone_into_arena(*node.lower, target, cache)) : std::nullopt;
        auto upper = node.upper ? std::optional<ExprPtr>(clone_into_arena(*node.upper, target, cache)) : std::nullopt;
        cloned = target.make<Integral>(integrand, Symbol(node.variable.name), lower, upper);
        break;
    }
    case ExprKind::Derivative: {
        const auto& node = expr_ref<Derivative>(expr);
        auto expression = clone_into_arena(node.expression, target, cache);
        cloned = target.make<Derivative>(expression, Symbol(node.variable.name), node.order);
        break;
    }
    case ExprKind::Limit: {
        const auto& node = expr_ref<Limit>(expr);
        auto expression = clone_into_arena(node.expression, target, cache);
        auto point = clone_into_arena(node.point, target, cache);
        cloned = target.make<Limit>(expression, Symbol(node.variable.name), point, node.direction);
        break;
    }
    case ExprKind::RootOf: {
        const auto& node = expr_ref<RootOf>(expr);
        auto polynomial = clone_into_arena(node.polynomial, target, cache);
        if (node.isolating_bound) {
            cloned = target.make<RootOf>(
                polynomial, Symbol(node.variable.name),
                *node.isolating_bound, node.root_index);
        } else {
            cloned = target.make<RootOf>(
                polynomial, Symbol(node.variable.name), node.root_index);
        }
        break;
    }
    case ExprKind::Matrix: {
        const auto& node = expr_ref<Matrix>(expr);
        std::vector<ExprPtr> elems;
        elems.reserve(node.elements.size());
        for (auto elem : node.elements) {
            elems.push_back(clone_into_arena(elem, target, cache));
        }
        cloned = target.make<Matrix>(node.rows, node.cols, std::move(elems));
        break;
    }
    case ExprKind::SeriesExp: {
        const auto& node = expr_ref<SeriesExp>(expr);
        std::vector<std::pair<long long, ExprPtr>> cloned_terms;
        cloned_terms.reserve(node.terms.size());
        for (const auto& [exp, coeff] : node.terms) {
            cloned_terms.push_back({exp, clone_into_arena(coeff, target, cache)});
        }
        cloned = target.make<SeriesExp>(Symbol(node.var.name), clone_into_arena(node.point, target, cache), std::move(cloned_terms), node.order);
        break;
    }
    case ExprKind::Quantity: {
        const auto& node = expr_ref<Quantity>(expr);
        cloned = target.make<Quantity>(clone_into_arena(node.value, target, cache), node.dimensions);
        break;
    }
    case ExprKind::Null:
        std::abort();
    }

    cache[expr] = cloned;
    return cloned;
}

}  // namespace cas
