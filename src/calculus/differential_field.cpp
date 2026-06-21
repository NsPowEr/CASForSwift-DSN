#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "cas/error_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"
#include "calculus_internal.hpp"
#include <algorithm>
#include <string>

namespace cas::calculus {

namespace {

template <typename F>
Result<void> visit_recursive_impl(ExprPtr expr, F&& f, unsigned int max_depth, unsigned int depth) {
    if (depth >= max_depth) {
        return fail<void>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Differential field visit recursion budget exceeded",
            .hint = std::nullopt,
        });
    }
    auto child_res = visit_expr(expr, [&](const auto& node) -> Result<void> {
        using T = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<T, Unary>) {
            return visit_recursive_impl(node.operand, f, max_depth, depth + 1U);
        } else if constexpr (std::is_same_v<T, Binary>) {
            auto r1 = visit_recursive_impl(node.left, f, max_depth, depth + 1U);
            if (r1.is_error()) return r1;
            return visit_recursive_impl(node.right, f, max_depth, depth + 1U);
        } else if constexpr (std::is_same_v<T, FuncCall>) {
            for (ExprPtr arg : node.args) {
                auto r = visit_recursive_impl(arg, f, max_depth, depth + 1U);
                if (r.is_error()) return r;
            }
        } else if constexpr (std::is_same_v<T, Sum>) {
            for (ExprPtr term : node.terms) {
                auto r = visit_recursive_impl(term, f, max_depth, depth + 1U);
                if (r.is_error()) return r;
            }
        } else if constexpr (std::is_same_v<T, Product>) {
            for (ExprPtr factor : node.factors) {
                auto r = visit_recursive_impl(factor, f, max_depth, depth + 1U);
                if (r.is_error()) return r;
            }
        }
        return ok();
    });
    if (child_res.is_error()) return child_res;

    return f(expr);
}

template <typename F>
Result<void> visit_recursive(ExprPtr expr, F&& f, unsigned int max_depth) {
    return visit_recursive_impl(expr, std::forward<F>(f), max_depth, 0U);
}

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

Result<DifferentialField> DifferentialField::build(ExprPtr expr, const Symbol& x, symbolic::CASContext& ctx) {
    DifferentialField field(x);
    auto res = field.add_extension(expr, ctx);
    if (res.is_error()) return fail<DifferentialField>(res.error());
    return ok(field);
}

Result<void> DifferentialField::add_extension(ExprPtr expr, symbolic::CASContext& ctx) {
    return visit_recursive(expr, [&](ExprPtr node) -> Result<void> {
        if (const auto* call = expr_cast<FuncCall>(node)) {
            ExtensionType type;
            if (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log) type = ExtensionType::Logarithmic;
            else if (call->func_id == BuiltinOp::Exp) type = ExtensionType::Exponential;
            else return ok();

            bool exists = false;
            for (const auto& ext : extensions_) {
                if (structural_equal(ext.argument, call->args[0]) && ext.type == type) {
                    exists = true;
                    break;
                }
            }

            if (!exists) {
                ExprPtr u = call->args[0];

                if (type == ExtensionType::Logarithmic) {
                    auto du_res = this->derive(u, ctx);
                    if (du_res.is_ok()) {
                        ExprPtr integrand = ctx.arena().make<Binary>(BinaryOp::Div, du_res.value(), u);
                        
                        bool is_zero = false;
                        if (auto s = ctx.simplify(integrand); s.is_ok()) {
                            integrand = s.value();
                            if (auto* il = expr_cast<IntegerLit>(integrand)) is_zero = il->value.is_zero();
                            else if (auto* rl = expr_cast<RationalLit>(integrand)) is_zero = rl->numerator.is_zero();
                        }
                        
                        if (!is_zero) {
                            bool is_dependent = false;
                            auto w_res = integrate(integrand, base_var_, ctx);
                            if (w_res.is_ok()) {
                                auto gen_res = this->to_field_generators(w_res.value(), ctx);
                                if (gen_res.is_ok()) {
                                    bool has_func_call = false;
                                    auto check_tree = [&](auto& self, ExprPtr e) -> void {
                                        if (!e) return;
                                        visit_expr(e, [&](const auto& n) -> ExprPtr {
                                            using T = std::decay_t<decltype(n)>;
                                            if constexpr (std::is_same_v<T, FuncCall>) {
                                                has_func_call = true;
                                            } else if constexpr (std::is_same_v<T, Unary>) {
                                                self(self, n.operand);
                                            } else if constexpr (std::is_same_v<T, Binary>) {
                                                self(self, n.left);
                                                self(self, n.right);
                                            } else if constexpr (std::is_same_v<T, Sum>) {
                                                for (const auto& a : n.terms) self(self, a);
                                            } else if constexpr (std::is_same_v<T, Product>) {
                                                for (const auto& a : n.factors) self(self, a);
                                            }
                                            return ExprPtr{};
                                        });
                                    };
                                    check_tree(check_tree, gen_res.value());
                                    if (!has_func_call) is_dependent = true;
                                }
                            }
                            if (is_dependent) {
                                return make_unimplemented<void>(
                                    "calculus", "add_extension",
                                    "Algebraically dependent logarithmic extension detected",
                                    cas::error::reason_codes::RISCH_LOG_EXTENSION_GENERAL,
                                    "Implement exact dependent extension replacement",
                                    "F5.1");
                            }
                        } else {
                            return ok();
                        }
                    }
                } else if (type == ExtensionType::Exponential) {
                    auto du_res = this->derive(u, ctx);
                    if (du_res.is_ok()) {
                        ExprPtr du = du_res.value();
                        
                        bool is_zero = false;
                        if (auto s = ctx.simplify(du); s.is_ok()) {
                            du = s.value();
                            if (auto* il = expr_cast<IntegerLit>(du)) is_zero = il->value.is_zero();
                            else if (auto* rl = expr_cast<RationalLit>(du)) is_zero = rl->numerator.is_zero();
                        }
                        
                        if (!is_zero) {
                            bool is_dependent = false;
                            auto int_res = integrate(du, base_var_, ctx);
                            if (int_res.is_ok()) {
                                auto gen_res = this->to_field_generators(int_res.value(), ctx);
                                if (gen_res.is_ok()) {
                                    bool has_non_log = false;
                                    auto check_tree = [&](auto& self, ExprPtr e) -> void {
                                        if (!e) return;
                                        visit_expr(e, [&](const auto& n) -> ExprPtr {
                                            using T = std::decay_t<decltype(n)>;
                                            if constexpr (std::is_same_v<T, Symbol>) {
                                                if (n.name == base_var_.name) has_non_log = true;
                                                else {
                                                    for (const auto& ext : extensions_) {
                                                        if (ext.t_var.name == n.name && ext.type == ExtensionType::Exponential) {
                                                            has_non_log = true;
                                                        }
                                                    }
                                                }
                                            } else if constexpr (std::is_same_v<T, FuncCall>) {
                                                has_non_log = true;
                                            } else if constexpr (std::is_same_v<T, Unary>) {
                                                self(self, n.operand);
                                            } else if constexpr (std::is_same_v<T, Binary>) {
                                                self(self, n.left);
                                                self(self, n.right);
                                            } else if constexpr (std::is_same_v<T, Sum>) {
                                                for (const auto& a : n.terms) self(self, a);
                                            } else if constexpr (std::is_same_v<T, Product>) {
                                                for (const auto& a : n.factors) self(self, a);
                                            }
                                            return ExprPtr{};
                                        });
                                    };
                                    check_tree(check_tree, gen_res.value());
                                    if (!has_non_log) is_dependent = true;
                                }
                            }
                            if (is_dependent) {
                                return make_unimplemented<void>(
                                    "calculus", "add_extension",
                                    "Algebraically dependent exponential extension detected",
                                    cas::error::reason_codes::RISCH_EXPONENTIAL_DE,
                                    "Implement exact dependent extension replacement",
                                    "F5.1");
                            }
                        } else {
                            return ok();
                        }
                    }
                }

                Symbol t("t_" + std::to_string(extensions_.size()));
                extensions_.push_back({type, call->args[0], t});
            }
        }
        return ok();
    }, ctx.diff_field_max_visit_depth());
}

Result<ExprPtr> DifferentialField::to_field_generators(ExprPtr expr, symbolic::CASContext& ctx) const {
    ExprPtr current = expr;
    AstArena& arena = ctx.arena();

    for (const auto& ext : extensions_) {
        std::string func_name = (ext.type == ExtensionType::Logarithmic) ? "ln" : "exp";
        ExprPtr pattern = arena.make<FuncCall>(func_name, std::vector<ExprPtr>{ext.argument});
        ExprPtr replacement = arena.make<Symbol>(ext.t_var.name);
        
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

Result<ExprPtr> DifferentialField::derive(ExprPtr expr, symbolic::CASContext& ctx) const {
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        if (sym->name == base_var_.name) return ok(ctx.arena().make<IntegerLit>(BigInt(1)));
        
        for (const auto& ext : extensions_) {
            if (sym->name == ext.t_var.name) {
                auto du = derive(ext.argument, ctx);
                if (du.is_error()) return du;
                
                if (ext.type == ExtensionType::Logarithmic) {
                    return ok(ctx.arena().make<Binary>(BinaryOp::Div, du.value(), ext.argument));
                } else {
                    return ok(ctx.arena().make<Product>(std::vector<ExprPtr>{
                        ctx.arena().make<Symbol>(ext.t_var.name), du.value()
                    }));
                }
            }
        }
    }
    
    return diff(expr, base_var_, 1U, ctx);
}

} // namespace cas::calculus
