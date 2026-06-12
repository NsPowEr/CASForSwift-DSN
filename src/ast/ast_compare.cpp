/// @file ast_compare.cpp
/// @brief Structural equality, hashing, and ExprKind naming for AST nodes.
///
/// Extracted from src/ast/ast.cpp (F8.0 / Task 1.2).
/// Contains: structural_equal, expr_hash, expr_kind, expr_kind_name.
#include "cas/ast.hpp"

#include <cstddef>
#include <string>

namespace cas {

namespace {

bool optional_expr_equal(const std::optional<ExprPtr>& lhs, const std::optional<ExprPtr>& rhs) noexcept {
    if (lhs.has_value() != rhs.has_value()) {
        return false;
    }
    if (!lhs.has_value()) {
        return true;
    }
    return structural_equal(*lhs, *rhs);
}

template <typename Container>
bool expr_ptr_range_equal(const Container& lhs, const Container& rhs) noexcept {
    if (lhs.size() != rhs.size()) {
        return false;
    }
    for (std::size_t index = 0; index < lhs.size(); ++index) {
        if (!structural_equal(lhs[index], rhs[index])) {
            return false;
        }
    }
    return true;
}

bool big_int_equal(const BigInt& lhs, const BigInt& rhs) noexcept {
    return lhs == rhs;
}

}  // namespace

// ---------------------------------------------------------------------------
// expr_kind
// ---------------------------------------------------------------------------

ExprKind expr_kind(ExprPtr expr) noexcept {
    return expr ? expr->kind : ExprKind::Null;
}

// ---------------------------------------------------------------------------
// structural_equal
// ---------------------------------------------------------------------------

bool structural_equal(ExprPtr lhs, ExprPtr rhs) noexcept {
    if (lhs == rhs) {
        return true;
    }
    if (!lhs || !rhs) {
        return false;
    }
    if (lhs->kind != rhs->kind) {
        return false;
    }

    switch (lhs->kind) {
    case ExprKind::IntegerLit: {
        const auto& lv = expr_ref<IntegerLit>(lhs);
        const auto& rv = expr_ref<IntegerLit>(rhs);
        return big_int_equal(lv.value, rv.value);
    }
    case ExprKind::RationalLit: {
        const auto& lv = expr_ref<RationalLit>(lhs);
        const auto& rv = expr_ref<RationalLit>(rhs);
        return big_int_equal(lv.numerator, rv.numerator) &&
               big_int_equal(lv.denominator, rv.denominator);
    }
    case ExprKind::ComplexLit: {
        const auto& l = expr_ref<ComplexLit>(lhs);
        const auto& r = expr_ref<ComplexLit>(rhs);
        return big_int_equal(l.re_num, r.re_num) && big_int_equal(l.re_den, r.re_den) &&
               big_int_equal(l.im_num, r.im_num) && big_int_equal(l.im_den, r.im_den);
    }
    case ExprKind::DecimalLit:
        return expr_ref<DecimalLit>(lhs).text == expr_ref<DecimalLit>(rhs).text;
    case ExprKind::Symbol:
        return expr_ref<Symbol>(lhs).name == expr_ref<Symbol>(rhs).name;
    case ExprKind::Constant:
        return expr_ref<Constant>(lhs).value == expr_ref<Constant>(rhs).value;
    case ExprKind::Unary: {
        const auto& lv = expr_ref<Unary>(lhs);
        const auto& rv = expr_ref<Unary>(rhs);
        return lv.op == rv.op && structural_equal(lv.operand, rv.operand);
    }
    case ExprKind::Binary: {
        const auto& lv = expr_ref<Binary>(lhs);
        const auto& rv = expr_ref<Binary>(rhs);
        return lv.op == rv.op &&
               structural_equal(lv.left, rv.left) &&
               structural_equal(lv.right, rv.right);
    }
    case ExprKind::FuncCall: {
        const auto& lv = expr_ref<FuncCall>(lhs);
        const auto& rv = expr_ref<FuncCall>(rhs);
        if (lv.func_id != rv.func_id) return false;
        if (lv.func_id == BuiltinOp::Unknown && lv.name != rv.name) return false;
        return expr_ptr_range_equal(lv.args, rv.args);
    }
    case ExprKind::Sum:
        return expr_ptr_range_equal(expr_ref<Sum>(lhs).terms, expr_ref<Sum>(rhs).terms);
    case ExprKind::Product:
        return expr_ptr_range_equal(expr_ref<Product>(lhs).factors, expr_ref<Product>(rhs).factors);
    case ExprKind::Integral: {
        const auto& lv = expr_ref<Integral>(lhs);
        const auto& rv = expr_ref<Integral>(rhs);
        return structural_equal(lv.integrand, rv.integrand) &&
               lv.variable.name == rv.variable.name &&
               optional_expr_equal(lv.lower, rv.lower) &&
               optional_expr_equal(lv.upper, rv.upper);
    }
    case ExprKind::Derivative: {
        const auto& lv = expr_ref<Derivative>(lhs);
        const auto& rv = expr_ref<Derivative>(rhs);
        return structural_equal(lv.expression, rv.expression) &&
               lv.variable.name == rv.variable.name &&
               lv.order == rv.order;
    }
    case ExprKind::Limit: {
        const auto& lv = expr_ref<Limit>(lhs);
        const auto& rv = expr_ref<Limit>(rhs);
        return structural_equal(lv.expression, rv.expression) &&
               lv.variable.name == rv.variable.name &&
               structural_equal(lv.point, rv.point) &&
               lv.direction == rv.direction;
    }
    case ExprKind::RootOf: {
        const auto& lv = expr_ref<RootOf>(lhs);
        const auto& rv = expr_ref<RootOf>(rhs);
        if (!structural_equal(lv.polynomial, rv.polynomial)) return false;
        if (lv.variable.name != rv.variable.name) return false;
        if (lv.root_index != rv.root_index) return false;
        if (lv.isolating_bound.has_value() != rv.isolating_bound.has_value())
            return false;
        if (lv.isolating_bound.has_value()) {
            const auto& lb = lv.isolating_bound.value();
            const auto& rb = rv.isolating_bound.value();
            if (lb.low_num != rb.low_num || lb.low_den != rb.low_den) return false;
            if (lb.high_num != rb.high_num || lb.high_den != rb.high_den) return false;
        }
        return true;
    }
    case ExprKind::Matrix: {
        const auto& lv = expr_ref<Matrix>(lhs);
        const auto& rv = expr_ref<Matrix>(rhs);
        return lv.rows == rv.rows &&
               lv.cols == rv.cols &&
               expr_ptr_range_equal(lv.elements, rv.elements);
    }
    case ExprKind::SeriesExp: {
        const auto& l = expr_ref<SeriesExp>(lhs);
        const auto& r = expr_ref<SeriesExp>(rhs);
        if (l.var.name != r.var.name || l.order != r.order) return false;
        if (!structural_equal(l.point, r.point)) return false;
        if (l.terms.size() != r.terms.size()) return false;
        for (std::size_t i = 0; i < l.terms.size(); ++i) {
            if (l.terms[i].first != r.terms[i].first) return false;
            if (!structural_equal(l.terms[i].second, r.terms[i].second)) return false;
        }
        return true;
    }
    case ExprKind::Quantity: {
        const auto& l = expr_ref<Quantity>(lhs);
        const auto& r = expr_ref<Quantity>(rhs);
        return l.dimensions == r.dimensions && structural_equal(l.value, r.value);
    }
    case ExprKind::Null:
        return false;
    }

    return false;
}

// ---------------------------------------------------------------------------
// expr_hash
// ---------------------------------------------------------------------------

std::size_t expr_hash(ExprPtr expr) noexcept {
    if (!expr) return 0;

    std::size_t cached = expr->hash_cache_.load(std::memory_order_relaxed);
    if (cached != 0) return cached;

    auto hash_combine = [](std::size_t& seed, std::size_t value) {
        seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };

    std::size_t seed = static_cast<std::size_t>(expr->kind);

    switch (expr->kind) {
    case ExprKind::IntegerLit:
        hash_combine(seed, expr_ref<IntegerLit>(expr).value.hash());
        break;
    case ExprKind::ComplexLit: {
        const auto& node = expr_ref<ComplexLit>(expr);
        hash_combine(seed, node.re_num.hash());
        hash_combine(seed, node.re_den.hash());
        hash_combine(seed, node.im_num.hash());
        hash_combine(seed, node.im_den.hash());
        break;
    }
    case ExprKind::RationalLit: {
        const auto& node = expr_ref<RationalLit>(expr);
        hash_combine(seed, node.numerator.hash());
        hash_combine(seed, node.denominator.hash());
        break;
    }
    case ExprKind::DecimalLit:
        hash_combine(seed, std::hash<std::string>{}(expr_ref<DecimalLit>(expr).text));
        break;
    case ExprKind::Symbol:
        hash_combine(seed, std::hash<std::string>{}(expr_ref<Symbol>(expr).name));
        break;
    case ExprKind::Constant:
        hash_combine(seed, static_cast<std::size_t>(expr_ref<Constant>(expr).value));
        break;
    case ExprKind::Unary: {
        const auto& node = expr_ref<Unary>(expr);
        hash_combine(seed, static_cast<std::size_t>(node.op));
        hash_combine(seed, expr_hash(node.operand));
        break;
    }
    case ExprKind::Binary: {
        const auto& node = expr_ref<Binary>(expr);
        hash_combine(seed, static_cast<std::size_t>(node.op));
        hash_combine(seed, expr_hash(node.left));
        hash_combine(seed, expr_hash(node.right));
        break;
    }
    case ExprKind::FuncCall: {
        const auto& node = expr_ref<FuncCall>(expr);
        hash_combine(seed, static_cast<std::size_t>(node.func_id));
        if (node.func_id == BuiltinOp::Unknown) {
            hash_combine(seed, std::hash<std::string>{}(node.name));
        }
        for (auto arg : node.args) hash_combine(seed, expr_hash(arg));
        break;
    }
    case ExprKind::Sum:
        for (auto term : expr_ref<Sum>(expr).terms) hash_combine(seed, expr_hash(term));
        break;
    case ExprKind::Product:
        for (auto factor : expr_ref<Product>(expr).factors) hash_combine(seed, expr_hash(factor));
        break;
    case ExprKind::Integral: {
        const auto& node = expr_ref<Integral>(expr);
        hash_combine(seed, expr_hash(node.integrand));
        hash_combine(seed, std::hash<std::string>{}(node.variable.name));
        if (node.lower) hash_combine(seed, expr_hash(*node.lower));
        if (node.upper) hash_combine(seed, expr_hash(*node.upper));
        break;
    }
    case ExprKind::Derivative: {
        const auto& node = expr_ref<Derivative>(expr);
        hash_combine(seed, expr_hash(node.expression));
        hash_combine(seed, std::hash<std::string>{}(node.variable.name));
        hash_combine(seed, node.order);
        break;
    }
    case ExprKind::Limit: {
        const auto& node = expr_ref<Limit>(expr);
        hash_combine(seed, expr_hash(node.expression));
        hash_combine(seed, std::hash<std::string>{}(node.variable.name));
        hash_combine(seed, expr_hash(node.point));
        hash_combine(seed, static_cast<std::size_t>(node.direction));
        break;
    }
    case ExprKind::RootOf: {
        const auto& node = expr_ref<RootOf>(expr);
        hash_combine(seed, expr_hash(node.polynomial));
        hash_combine(seed, std::hash<std::string>{}(node.variable.name));
        if (node.root_index) hash_combine(seed, *node.root_index);
        if (node.isolating_bound) {
            const auto& b = *node.isolating_bound;
            hash_combine(seed, b.low_num.hash());
            hash_combine(seed, b.low_den.hash());
            hash_combine(seed, b.high_num.hash());
            hash_combine(seed, b.high_den.hash());
        }
        break;
    }
    case ExprKind::Matrix: {
        const auto& node = expr_ref<Matrix>(expr);
        hash_combine(seed, node.rows);
        hash_combine(seed, node.cols);
        for (auto elem : node.elements) hash_combine(seed, expr_hash(elem));
        break;
    }
    case ExprKind::SeriesExp: {
        const auto& node = expr_ref<SeriesExp>(expr);
        hash_combine(seed, std::hash<std::string>{}(node.var.name));
        hash_combine(seed, expr_hash(node.point));
        hash_combine(seed, static_cast<std::size_t>(node.order));
        for (const auto& [exp, coeff] : node.terms) {
            hash_combine(seed, static_cast<std::size_t>(exp));
            hash_combine(seed, expr_hash(coeff));
        }
        break;
    }
    case ExprKind::Quantity: {
        const auto& node = expr_ref<Quantity>(expr);
        hash_combine(seed, expr_hash(node.value));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.m));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.kg));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.s));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.A));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.K));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.mol));
        hash_combine(seed, static_cast<std::size_t>(node.dimensions.cd));
        break;
    }
    case ExprKind::Null:
        break;
    }

    std::size_t final_hash = (seed == 0) ? 1 : seed;
    expr->hash_cache_.store(final_hash, std::memory_order_relaxed);
    return final_hash;
}

// ---------------------------------------------------------------------------
// expr_kind_name
// ---------------------------------------------------------------------------

std::string_view expr_kind_name(ExprKind kind) noexcept {
    switch (kind) {
    case ExprKind::Null:        return "Null";
    case ExprKind::IntegerLit:  return "IntegerLit";
    case ExprKind::RationalLit: return "RationalLit";
    case ExprKind::ComplexLit:  return "ComplexLit";
    case ExprKind::DecimalLit:  return "DecimalLit";
    case ExprKind::Symbol:      return "Symbol";
    case ExprKind::Constant:    return "Constant";
    case ExprKind::Unary:       return "Unary";
    case ExprKind::Binary:      return "Binary";
    case ExprKind::FuncCall:    return "FuncCall";
    case ExprKind::Sum:         return "Sum";
    case ExprKind::Product:     return "Product";
    case ExprKind::Integral:    return "Integral";
    case ExprKind::Derivative:  return "Derivative";
    case ExprKind::Limit:       return "Limit";
    case ExprKind::RootOf:      return "RootOf";
    case ExprKind::Matrix:      return "Matrix";
    case ExprKind::SeriesExp:   return "SeriesExp";
    case ExprKind::Quantity:    return "Quantity";
    }
    return "Unknown";
}

}  // namespace cas
