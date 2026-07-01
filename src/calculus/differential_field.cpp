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

// True iff any tower generator symbol t_i occurs anywhere in `expr`.  Drives the
// split in DifferentialField::derive: a generator-free subexpression is a pure
// function of the base variable (and field constants), so the engine's
// differentiator is authoritative; only when a generator is present must we apply
// the tower-aware derivation so the known D(t_i) propagate through the operators.
bool contains_generator(ExprPtr expr, const std::vector<DifferentialExtension>& exts) {
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        for (const auto& ext : exts)
            if (sym->name == ext.t_var.name) return true;
        return false;
    }
    if (const auto* u = expr_cast<Unary>(expr))
        return contains_generator(u->operand, exts);
    if (const auto* b = expr_cast<Binary>(expr))
        return contains_generator(b->left, exts) || contains_generator(b->right, exts);
    if (const auto* fc = expr_cast<FuncCall>(expr)) {
        for (ExprPtr a : fc->args) if (contains_generator(a, exts)) return true;
        return false;
    }
    if (const auto* s = expr_cast<Sum>(expr)) {
        for (ExprPtr t : s->terms) if (contains_generator(t, exts)) return true;
        return false;
    }
    if (const auto* p = expr_cast<Product>(expr)) {
        for (ExprPtr f : p->factors) if (contains_generator(f, exts)) return true;
        return false;
    }
    return false;
}

// The derivation D of the differential field K(t_1,…,t_n), applied structurally so
// the known D(t_i) propagate through the field operators (Bronstein §3.4).  A
// generator-free subexpression is a pure function of the base variable, where the
// engine's differentiator is authoritative; only generator-bearing nodes need the
// tower-aware Leibniz rules.  `depth` is bounded by ctx.max_recursion_depth()
// (CASContext budget, hardcode-catalog category 1) — an adversarially deep AST
// yields a diagnostic Unimplemented, never a stack overflow.
[[maybe_unused]] Result<ExprPtr> derive_in_field(const DifferentialField& field, ExprPtr expr,
                                symbolic::CASContext& ctx, std::size_t depth) {
    AstArena& arena = ctx.arena();
    const auto& exts = field.extensions();

    if (depth >= ctx.max_recursion_depth()) {
        return make_unimplemented<ExprPtr>(
            "calculus", "DifferentialField::derive",
            "AST nesting exceeds the configured recursion budget",
            "DIFFERENTIAL_FIELD_DERIVE_DEPTH",
            "Raise CASContext::set_max_recursion_depth for genuinely deep towers");
    }

    // Generator-free subexpression: a pure function of the base variable (and
    // field constants).  The engine's differentiator is authoritative here — it
    // treats every non-base symbol as a constant, exactly the field semantics for
    // parameters (bare base_var → 1, parameters → 0, sin(x) → cos(x)…).  This also
    // preserves the previous behaviour for all such inputs.
    if (!contains_generator(expr, exts)) {
        return diff(expr, field.base_var(), 1U, ctx);
    }

    // Bare generator symbol: D(t_i) = u'/u (log) or t_i·u' (exp).
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        for (const auto& ext : exts) {
            if (sym->name == ext.t_var.name) {
                auto du = derive_in_field(field, ext.argument, ctx, depth + 1U);
                if (du.is_error()) return du;
                if (ext.type == ExtensionType::Logarithmic) {
                    return ok(arena.make<Binary>(BinaryOp::Div, du.value(), ext.argument));
                }
                return ok(arena.make<Product>(std::vector<ExprPtr>{
                    arena.make<Symbol>(ext.t_var.name), du.value()}));
            }
        }
        // contains_generator was true but this symbol is not a generator:
        // unreachable, but stay total.
        return ok(arena.make<IntegerLit>(BigInt(0)));
    }

    // Unary: D(−a) = −D(a).
    if (const auto* un = expr_cast<Unary>(expr)) {
        if (un->op == UnaryOp::Neg) {
            auto d = derive_in_field(field, un->operand, ctx, depth + 1U);
            if (d.is_error()) return d;
            return ok(arena.make<Unary>(UnaryOp::Neg, d.value()));
        }
        return make_unimplemented<ExprPtr>(
            "calculus", "DifferentialField::derive",
            "unary operator other than negation applied to a field generator",
            "DIFFERENTIAL_FIELD_DERIVE_NONFIELD",
            "Only elementary field operators (+,−,·,/,^const) admit a closed-form "
            "derivation inside K(t_1,…,t_n)");
    }

    // Sum: D(Σ a_i) = Σ D(a_i).
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms) {
            auto d = derive_in_field(field, t, ctx, depth + 1U);
            if (d.is_error()) return d;
            terms.push_back(d.value());
        }
        return ok(arena.make<Sum>(std::move(terms)));
    }

    // Product: D(∏ a_i) = Σ_i D(a_i)·∏_{j≠i} a_j.
    if (const auto* prod = expr_cast<Product>(expr)) {
        const auto& fs = prod->factors;
        std::vector<ExprPtr> terms;
        terms.reserve(fs.size());
        for (std::size_t i = 0; i < fs.size(); ++i) {
            auto d = derive_in_field(field, fs[i], ctx, depth + 1U);
            if (d.is_error()) return d;
            std::vector<ExprPtr> factors;
            factors.reserve(fs.size());
            factors.push_back(d.value());
            for (std::size_t j = 0; j < fs.size(); ++j)
                if (j != i) factors.push_back(fs[j]);
            terms.push_back(arena.make<Product>(std::move(factors)));
        }
        return ok(arena.make<Sum>(std::move(terms)));
    }

    // Binary.
    if (const auto* bin = expr_cast<Binary>(expr)) {
        switch (bin->op) {
            case BinaryOp::Add:
            case BinaryOp::Sub: {
                auto dl = derive_in_field(field, bin->left, ctx, depth + 1U);
                if (dl.is_error()) return dl;
                auto dr = derive_in_field(field, bin->right, ctx, depth + 1U);
                if (dr.is_error()) return dr;
                return ok(arena.make<Binary>(bin->op, dl.value(), dr.value()));
            }
            case BinaryOp::Mul: {
                // (u·v)' = u'·v + u·v'
                auto du = derive_in_field(field, bin->left, ctx, depth + 1U);
                if (du.is_error()) return du;
                auto dv = derive_in_field(field, bin->right, ctx, depth + 1U);
                if (dv.is_error()) return dv;
                ExprPtr t1 = arena.make<Binary>(BinaryOp::Mul, du.value(), bin->right);
                ExprPtr t2 = arena.make<Binary>(BinaryOp::Mul, bin->left, dv.value());
                return ok(arena.make<Binary>(BinaryOp::Add, t1, t2));
            }
            case BinaryOp::Div: {
                // (u/v)' = (u'·v − u·v')/v²
                auto du = derive_in_field(field, bin->left, ctx, depth + 1U);
                if (du.is_error()) return du;
                auto dv = derive_in_field(field, bin->right, ctx, depth + 1U);
                if (dv.is_error()) return dv;
                ExprPtr num1 = arena.make<Binary>(BinaryOp::Mul, du.value(), bin->right);
                ExprPtr num2 = arena.make<Binary>(BinaryOp::Mul, bin->left, dv.value());
                ExprPtr num = arena.make<Binary>(BinaryOp::Sub, num1, num2);
                ExprPtr den = arena.make<Binary>(BinaryOp::Mul, bin->right, bin->right);
                return ok(arena.make<Binary>(BinaryOp::Div, num, den));
            }
            case BinaryOp::Pow: {
                // D(u^n) = n·u^{n−1}·u' — valid only for an exponent that is
                // constant in the field (generator-free AND base-independent).  A
                // base_var- or generator-dependent exponent needs log(u), outside
                // K → diagnostic Unimplemented (REGOLA ZERO: no silent wrong
                // result).  Any rational constant exponent is admitted (not just
                // integers): the check is D(exponent) ≡ 0, not a type test.
                ExprPtr base = bin->left;
                ExprPtr exp = bin->right;
                auto nonfield = [&]() {
                    return make_unimplemented<ExprPtr>(
                        "calculus", "DifferentialField::derive",
                        "power with a non-constant exponent over a field generator "
                        "(would require log of the base, outside K(t_1,…,t_n))",
                        "DIFFERENTIAL_FIELD_DERIVE_NONFIELD",
                        "Represent base^exponent as an explicit exp/log tower "
                        "generator before differentiating");
                };
                if (contains_generator(exp, exts)) return nonfield();
                auto dexp = diff(exp, field.base_var(), 1U, ctx);
                if (dexp.is_error()) return dexp;
                ExprPtr dexp_v = dexp.value();
                if (auto s = ctx.simplify(dexp_v); s.is_ok()) dexp_v = s.value();
                bool exp_const = false;
                if (const auto* il = expr_cast<IntegerLit>(dexp_v)) exp_const = il->value.is_zero();
                else if (const auto* rl = expr_cast<RationalLit>(dexp_v)) exp_const = rl->numerator.is_zero();
                if (!exp_const) return nonfield();

                auto du = derive_in_field(field, base, ctx, depth + 1U);
                if (du.is_error()) return du;
                ExprPtr n_minus_1 = arena.make<Binary>(BinaryOp::Sub, exp,
                                        arena.make<IntegerLit>(BigInt(1)));
                ExprPtr u_pow = arena.make<Binary>(BinaryOp::Pow, base, n_minus_1);
                ExprPtr coeff = arena.make<Binary>(BinaryOp::Mul, exp, u_pow);
                return ok(arena.make<Binary>(BinaryOp::Mul, coeff, du.value()));
            }
            default:
                return make_unimplemented<ExprPtr>(
                    "calculus", "DifferentialField::derive",
                    "binary operator without a field derivation rule applied to a "
                    "generator",
                    "DIFFERENTIAL_FIELD_DERIVE_NONFIELD",
                    "Only +,−,·,/,^const admit a closed-form derivation in "
                    "K(t_1,…,t_n)");
        }
    }

    // FuncCall wrapping a generator (e.g. sin(t)) or any other node: outside the
    // polynomial/rational field representation the Risch pipeline constructs.
    return make_unimplemented<ExprPtr>(
        "calculus", "DifferentialField::derive",
        "expression node with a field generator that is not a polynomial/rational "
        "combination of the generators",
        "DIFFERENTIAL_FIELD_DERIVE_NONFIELD",
        "The Risch differential-field solver only derives rational functions of the "
        "tower generators");
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

Result<ExprPtr> DifferentialField::derive(ExprPtr expr, symbolic::CASContext& ctx) const {
    return derive_in_field(*this, expr, ctx, 0U);
}

} // namespace cas::calculus
