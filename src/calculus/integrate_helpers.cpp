#include "integrate_engine.hpp"
#include "cas/error.hpp"
#include "cas/algebra.hpp"
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value) {
    return arena.make<IntegerLit>(BigInt(value));
}

[[nodiscard]] ExprPtr make_rational(AstArena& arena, long long numerator, long long denominator) {
    return arena.make<RationalLit>(BigInt(numerator), BigInt(denominator));
}

[[nodiscard]] ExprPtr make_rational(AstArena& arena, const Rational& value) {
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] ExprPtr make_unary(AstArena& arena, UnaryOp op, ExprPtr operand) {
    return arena.make<Unary>(op, operand);
}

[[nodiscard]] ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs) {
    return arena.make<Binary>(op, lhs, rhs);
}

[[nodiscard]] ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms) {
    return arena.make<Sum>(std::move(terms));
}

[[nodiscard]] ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors) {
    return arena.make<Product>(std::move(factors));
}

[[nodiscard]] ExprPtr make_function(AstArena& arena, std::string name, std::vector<ExprPtr> args) {
    return arena.make<FuncCall>(std::move(name), std::move(args));
}

[[nodiscard]] std::string canonical_function_name(const std::string& name) {
    return name;
}

[[nodiscard]] bool depends_on(ExprPtr expr, const Symbol& var) {
    if (!expr) return false;
    if (const auto* s = expr_cast<Symbol>(expr)) return s->name == var.name;
    bool dep = false;
    visit_expr(expr, [&](const auto& node) {
        using NodeT = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<NodeT, Symbol>) {
            if (node.name == var.name) dep = true;
        } else if constexpr (std::is_same_v<NodeT, Unary>) {
            if (depends_on(node.operand, var)) dep = true;
        } else if constexpr (std::is_same_v<NodeT, Binary>) {
            if (depends_on(node.left, var) || depends_on(node.right, var)) dep = true;
        } else if constexpr (std::is_same_v<NodeT, FuncCall>) {
            for (auto a : node.args) if (depends_on(a, var)) dep = true;
        } else if constexpr (std::is_same_v<NodeT, Sum>) {
            for (auto t : node.terms) if (depends_on(t, var)) dep = true;
        } else if constexpr (std::is_same_v<NodeT, Product>) {
            for (auto f : node.factors) if (depends_on(f, var)) dep = true;
        }
    });
    return dep;
}

[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var) {
    if (const auto* s = expr_cast<Symbol>(expr)) return s->name == var.name;
    return false;
}

[[nodiscard]] bool is_one(ExprPtr expr) {
    if (const auto* i = expr_cast<IntegerLit>(expr)) return i->value == BigInt(1);
    if (const auto* r = expr_cast<RationalLit>(expr)) return r->numerator == BigInt(1) && r->denominator == BigInt(1);
    return false;
}

[[nodiscard]] bool is_negative_one(ExprPtr expr) {
    if (const auto* i = expr_cast<IntegerLit>(expr)) return i->value == BigInt(-1);
    return false;
}

[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr) {
    if (const auto* i = expr_cast<IntegerLit>(expr)) return Rational(i->value);
    if (const auto* r = expr_cast<RationalLit>(expr)) return Rational(r->numerator, r->denominator);
    if (const auto* u = expr_cast<Unary>(expr); u != nullptr && u->op == UnaryOp::Neg) {
        auto value = exact_scalar_from_expr(u->operand);
        if (value.has_value()) return -value.value();
    }
    if (const auto* b = expr_cast<Binary>(expr); b != nullptr) {
        auto left = exact_scalar_from_expr(b->left);
        auto right = exact_scalar_from_expr(b->right);
        if (b->op == BinaryOp::Mul && left.has_value() && right.has_value()) return left.value() * right.value();
        if (b->op == BinaryOp::Div && left.has_value() && right.has_value() && !right->numerator().is_zero()) return left.value() / right.value();
    }
    if (const auto* p = expr_cast<Product>(expr); p != nullptr) {
        Rational result(1);
        for (ExprPtr factor : p->factors) {
            auto value = exact_scalar_from_expr(factor);
            if (!value.has_value()) return std::nullopt;
            result *= value.value();
        }
        return result;
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<AffineArgument> extract_affine_argument(ExprPtr expr, const Symbol& var) {
    if (is_same_symbol(expr, var)) return AffineArgument{Rational(1), Rational(0)};
    if (!depends_on(expr, var)) {
        auto scalar = exact_scalar_from_expr(expr);
        if (scalar.has_value()) return AffineArgument{Rational(0), scalar.value()};
        return std::nullopt;
    }
    if (const auto* u = expr_cast<Unary>(expr); u != nullptr && u->op == UnaryOp::Neg) {
        auto inner = extract_affine_argument(u->operand, var);
        if (inner.has_value()) return AffineArgument{-inner->coefficient, -inner->constant};
    }
    if (const auto* b = expr_cast<Binary>(expr)) {
        if (b->op == BinaryOp::Add) {
            auto l = extract_affine_argument(b->left, var);
            auto r = extract_affine_argument(b->right, var);
            if (l && r) return AffineArgument{l->coefficient + r->coefficient, l->constant + r->constant};
        }
        if (b->op == BinaryOp::Sub) {
            auto l = extract_affine_argument(b->left, var);
            auto r = extract_affine_argument(b->right, var);
            if (l && r) return AffineArgument{l->coefficient - r->coefficient, l->constant - r->constant};
        }
        if (b->op == BinaryOp::Mul) {
            auto left_scalar = exact_scalar_from_expr(b->left);
            auto right_scalar = exact_scalar_from_expr(b->right);
            if (left_scalar.has_value()) {
                auto right = extract_affine_argument(b->right, var);
                if (right.has_value()) return AffineArgument{left_scalar.value() * right->coefficient, left_scalar.value() * right->constant};
            }
            if (right_scalar.has_value()) {
                auto left = extract_affine_argument(b->left, var);
                if (left.has_value()) return AffineArgument{right_scalar.value() * left->coefficient, right_scalar.value() * left->constant};
            }
        }
        if (b->op == BinaryOp::Div) {
            auto denominator = exact_scalar_from_expr(b->right);
            if (denominator.has_value() && !denominator->numerator().is_zero()) {
                auto numerator = extract_affine_argument(b->left, var);
                if (numerator.has_value()) return AffineArgument{numerator->coefficient / denominator.value(), numerator->constant / denominator.value()};
            }
        }
    }
    if (const auto* p = expr_cast<Product>(expr); p != nullptr) {
        Rational scalar(1);
        ExprPtr affine_factor{};
        for (ExprPtr factor : p->factors) {
            auto value = exact_scalar_from_expr(factor);
            if (value.has_value()) {
                scalar *= value.value();
                continue;
            }
            if (affine_factor) return std::nullopt;
            affine_factor = factor;
        }
        if (affine_factor) {
            auto affine = extract_affine_argument(affine_factor, var);
            if (affine.has_value()) return AffineArgument{scalar * affine->coefficient, scalar * affine->constant};
        }
    }
    if (const auto* s = expr_cast<Sum>(expr); s != nullptr) {
        Rational coeff(0);
        Rational constant(0);
        for (ExprPtr term : s->terms) {
            auto affine = extract_affine_argument(term, var);
            if (!affine.has_value()) return std::nullopt;
            coeff += affine->coefficient;
            constant += affine->constant;
        }
        return AffineArgument{coeff, constant};
    }
    return std::nullopt;
}

[[nodiscard]] bool matches_square_of_variable(ExprPtr expr, const Symbol& var) {
    const auto* b = expr_cast<Binary>(expr);
    if (b && b->op == BinaryOp::Pow && is_same_symbol(b->left, var) && is_rational_value(b->right, 2, 1)) {
        return true;
    }
    if (const auto* p = expr_cast<Product>(expr); p != nullptr && p->factors.size() == 2U) {
        return is_same_symbol(p->factors[0], var) && is_same_symbol(p->factors[1], var);
    }
    return false;
}

[[nodiscard]] bool is_rational_value(ExprPtr expr, long long numerator, long long denominator) {
    if (const auto* i = expr_cast<IntegerLit>(expr)) return denominator == 1 && i->value == BigInt(numerator);
    if (const auto* r = expr_cast<RationalLit>(expr)) return r->numerator == BigInt(numerator) && r->denominator == BigInt(denominator);
    return false;
}

[[nodiscard]] bool matches_one_plus_square(ExprPtr expr, const Symbol& var) {
    const auto* b = expr_cast<Binary>(expr);
    if (b && b->op == BinaryOp::Add) {
        return (is_one(b->left) && matches_square_of_variable(b->right, var)) ||
               (is_one(b->right) && matches_square_of_variable(b->left, var));
    }
    if (const auto* sum = expr_cast<Sum>(expr); sum != nullptr && sum->terms.size() == 2U) {
        return (is_one(sum->terms[0]) && matches_square_of_variable(sum->terms[1], var)) ||
               (is_one(sum->terms[1]) && matches_square_of_variable(sum->terms[0], var));
    }
    return false;
}

[[nodiscard]] bool matches_one_minus_square(ExprPtr expr, const Symbol& var) {
    const auto* b = expr_cast<Binary>(expr);
    if (!b || b->op != BinaryOp::Sub) return false;
    return is_one(b->left) && matches_square_of_variable(b->right, var);
}

// Try to extract C from C^2 (the constant factor in x^2 ± C^2 patterns).
// Returns true and sets base=C if the expression is C^2 with C independent of var.
[[nodiscard]] static bool extract_constant_squared(ExprPtr expr, const Symbol& var, ExprPtr& base) {
    if (depends_on(expr, var)) return false;
    if (const auto* bp = expr_cast<Binary>(expr); bp && bp->op == BinaryOp::Pow) {
        if (is_rational_value(bp->right, 2, 1) && !depends_on(bp->left, var)) {
            base = bp->left;
            return true;
        }
    }
    // Accept any constant expression C — note: the caller uses base as-is in 1/(2*base),
    // so only accept when we know base itself (not base^2). For plain symbols/integers, use sqrt.
    // Conservative: only accept the Binary(Pow, C, 2) form.
    return false;
}

[[nodiscard]] bool matches_square_plus_constant_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base) {
    // Matches x^2 + C^2 in Sum form (after simplification)
    if (const auto* sum = expr_cast<Sum>(expr); sum && sum->terms.size() == 2U) {
        if (matches_square_of_variable(sum->terms[0], var) && extract_constant_squared(sum->terms[1], var, constant_base)) return true;
        if (matches_square_of_variable(sum->terms[1], var) && extract_constant_squared(sum->terms[0], var, constant_base)) return true;
    }
    // Binary Add form
    if (const auto* b = expr_cast<Binary>(expr); b && b->op == BinaryOp::Add) {
        if (matches_square_of_variable(b->left, var) && extract_constant_squared(b->right, var, constant_base)) return true;
        if (matches_square_of_variable(b->right, var) && extract_constant_squared(b->left, var, constant_base)) return true;
    }
    return false;
}

[[nodiscard]] bool matches_square_minus_constant_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base) {
    // Matches x^2 - C^2: Binary(Sub) form or Sum([x^2, Unary(Neg, C^2)]) form
    if (const auto* b = expr_cast<Binary>(expr); b && b->op == BinaryOp::Sub) {
        if (matches_square_of_variable(b->left, var) && extract_constant_squared(b->right, var, constant_base)) return true;
    }
    if (const auto* sum = expr_cast<Sum>(expr); sum && sum->terms.size() == 2U) {
        for (int i = 0; i < 2; ++i) {
            ExprPtr t_sq = sum->terms[i];
            ExprPtr t_neg = sum->terms[1-i];
            if (!matches_square_of_variable(t_sq, var)) continue;
            // t_neg should be Unary(Neg, C^2) or Product([-1, C^2])
            ExprPtr inner{};
            if (const auto* u = expr_cast<Unary>(t_neg); u && u->op == UnaryOp::Neg) inner = u->operand;
            else if (const auto* p = expr_cast<Product>(t_neg); p && p->factors.size() == 2U) {
                if (is_negative_one(p->factors[0])) inner = p->factors[1];
                else if (is_negative_one(p->factors[1])) inner = p->factors[0];
            }
            if (inner && extract_constant_squared(inner, var, constant_base)) return true;
        }
    }
    return false;
}

[[nodiscard]] bool matches_constant_square_minus_variable_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base) {
    // Matches C^2 - x^2: Binary(Sub) form or Sum([C^2, Unary(Neg, x^2)])
    if (const auto* b = expr_cast<Binary>(expr); b && b->op == BinaryOp::Sub) {
        if (extract_constant_squared(b->left, var, constant_base) && matches_square_of_variable(b->right, var)) return true;
    }
    if (const auto* sum = expr_cast<Sum>(expr); sum && sum->terms.size() == 2U) {
        for (int i = 0; i < 2; ++i) {
            ExprPtr t_csq = sum->terms[i];
            ExprPtr t_neg = sum->terms[1-i];
            ExprPtr inner{};
            if (const auto* u = expr_cast<Unary>(t_neg); u && u->op == UnaryOp::Neg) inner = u->operand;
            else if (const auto* p = expr_cast<Product>(t_neg); p && p->factors.size() == 2U) {
                if (is_negative_one(p->factors[0])) inner = p->factors[1];
                else if (is_negative_one(p->factors[1])) inner = p->factors[0];
            }
            if (inner && extract_constant_squared(t_csq, var, constant_base) && matches_square_of_variable(inner, var)) return true;
        }
    }
    return false;
}

[[nodiscard]] bool matches_reciprocal_sqrt_one_minus_square(ExprPtr expr, const Symbol& var) {
    const auto* b = expr_cast<Binary>(expr);
    if (!b || b->op != BinaryOp::Pow || !is_negative_one(b->right)) return false;
    const auto* f = expr_cast<FuncCall>(b->left);
    return f && f->func_id == BuiltinOp::Sqrt && f->args.size() == 1 && matches_one_minus_square(f->args[0], var);
}

[[nodiscard]] Result<ExprPtr> integrate_improper_residues(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    auto parts = algebra::apart_num_den(expr, ctx);
    if (parts.is_error()) return fail<ExprPtr>(parts.error());
    auto roots = algebra::solve_polynomial(parts.value().denominator, var, ctx);
    if (roots.is_error()) return fail<ExprPtr>(roots.error());

    ExprPtr sum_residues = make_integer(ctx.arena(), 0);
    ExprPtr i_const = ctx.arena().make<Constant>(MathConstant::I);

    for (auto root : roots.value()) {
        auto img_part_res = ctx.simplify(ctx.arena().make<FuncCall>("im", std::vector<ExprPtr>{root}));
        double img_val = 0.0;
        if (img_part_res.is_ok()) {
            if (const auto* dec = expr_cast<DecimalLit>(img_part_res.value())) img_val = dec->to_double();
            else img_val = 1.0; 
        }

        if (img_val > 0) {
            auto q_prime = diff(parts.value().denominator, var, 1U, ctx);
            if (q_prime.is_error()) continue;
            auto p_val = ctx.substitute(parts.value().numerator, var, root);
            auto q_p_val = ctx.substitute(q_prime.value(), var, root);
            if (p_val.is_error() || q_p_val.is_error()) continue;
            auto res = ctx.simplify(make_binary(ctx.arena(), BinaryOp::Div, p_val.value(), q_p_val.value()));
            if (res.is_ok()) sum_residues = make_sum(ctx.arena(), {sum_residues, res.value()});
        }
    }
    ExprPtr pi_const = ctx.arena().make<Constant>(MathConstant::Pi);
    ExprPtr factor = make_product(ctx.arena(), {make_integer(ctx.arena(), 2), pi_const, i_const, sum_residues});
    return ctx.simplify(factor);
}

} // namespace cas::calculus::integrate_detail
