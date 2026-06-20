#include "simplify_impl.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/numeric.hpp"
#include "cas/algebra.hpp"
#include "../algebra/polynomial_internal.hpp"

namespace cas::symbolic::detail {

// ── simplify_node(FuncCall) — pure dispatcher ─────────────────────────────────

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const FuncCall& node) {
    std::vector<ExprPtr> args;
    args.reserve(node.args.size());
    for (std::size_t i = 0; i < node.args.size(); ++i) {
        auto simplify_arg = [&]() -> Result<ExprPtr> {
            if (trace_enabled_) {
                std::vector<ExprPtr> current_args = node.args;
                for (std::size_t j = 0; j < args.size(); ++j) current_args[j] = args[j];
                ScopedFrame frame(*this, [this, current_args = std::move(current_args), i,
                                          name = node.name](ExprPtr value) mutable {
                    current_args[i] = value;
                    return arena_.make<FuncCall>(name, std::move(current_args));
                });
                return simplify_expr(node.args[i]);
            }
            return simplify_expr(node.args[i]);
        };
        auto res = simplify_arg();
        if (res.is_error()) return res;
        args.push_back(res.value());
    }

    const ExprPtr target_before = expr_ptr_sequence_identical(args, node.args)
        ? original
        : (trace_enabled_ ? arena_.make<FuncCall>(node.func_id, args) : ExprPtr{});

    if (rewrite_provider_ != nullptr && may_rewrite_function_call(node.func_id, args)) {
        ExprPtr rewrite_target = expr_ptr_sequence_identical(args, node.args)
            ? original
            : arena_.make<FuncCall>(node.func_id, args);
        auto rewritten = rewrite_provider_->try_rewrite(
            rewrite_target, arena_, assumptions_, context_);
        if (rewritten.is_ok() && rewritten.value() != rewrite_target) {
            append_trace(RuleId::RewriteProviderApplied, rewrite_target, rewritten.value());
            return simplify_expr(rewritten.value());
        }
    }

    const BuiltinOp op = node.func_id;

    switch (op) {
    case BuiltinOp::Sin:
    case BuiltinOp::Cos:
    case BuiltinOp::Tan:
    case BuiltinOp::Asin:
    case BuiltinOp::Acos:
    case BuiltinOp::Atan:
        return simplify_funcall_trig(original, op, std::move(args), target_before);

    case BuiltinOp::Exp:
    case BuiltinOp::Ln:
    case BuiltinOp::Log:
    case BuiltinOp::Sqrt:
        return simplify_funcall_exp_log_sqrt(original, op, std::move(args), target_before);

    case BuiltinOp::Gamma:
    case BuiltinOp::Digamma:
    case BuiltinOp::Polygamma:
    case BuiltinOp::Beta:
    case BuiltinOp::Pochhammer:
    case BuiltinOp::Erf:
    case BuiltinOp::Zeta:
    case BuiltinOp::Hypergeometric0F1:
    case BuiltinOp::Hypergeometric1F1:
    case BuiltinOp::Hypergeometric2F1:
    case BuiltinOp::EllipticK:
    case BuiltinOp::EllipticE:
    case BuiltinOp::EllipticPi:
    case BuiltinOp::EllipticF:
        return simplify_funcall_special(original, op, std::move(args), target_before);

    case BuiltinOp::BesselJ:
    case BuiltinOp::BesselY:
    case BuiltinOp::BesselI:
    case BuiltinOp::BesselK:
    case BuiltinOp::BesselZero:
    case BuiltinOp::ChebyshevT:
    case BuiltinOp::ChebyshevU:
    case BuiltinOp::HermiteH:
    case BuiltinOp::HermiteHe:
    case BuiltinOp::JacobiP:
    case BuiltinOp::LambertW:
    case BuiltinOp::LaguerreL:
    case BuiltinOp::LegendreP:
        return simplify_funcall_bessel_orthogonal(original, op, std::move(args), target_before);

    case BuiltinOp::Sinh:
    case BuiltinOp::Cosh:
    case BuiltinOp::Tanh:
    case BuiltinOp::Coth:
        return simplify_funcall_hyperbolic(original, op, std::move(args), target_before);

    case BuiltinOp::Factorial:
    case BuiltinOp::Binomial:
    case BuiltinOp::Erfc:
        return simplify_funcall_combinatorial(original, op, std::move(args), target_before);

    case BuiltinOp::Re:
    case BuiltinOp::Im:
    case BuiltinOp::Conj:
    case BuiltinOp::Abs:
    case BuiltinOp::Arg:
    case BuiltinOp::Sign:
        return simplify_funcall_complex(original, op, std::move(args), target_before);

    default:
        break;
    }

    // Matrix operations — Det/Rank/Trace/Inv/Transpose
    if (context_ != nullptr && args.size() == 1U && expr_is<Matrix>(args.front())) {
        const auto& m_node = expr_ref<Matrix>(args.front());
        cas::linalg::MatrixExpr m_expr(m_node.rows, m_node.cols, m_node.elements);
        if (op == BuiltinOp::Det) {
            auto res = cas::linalg::determinant(m_expr, *context_);
            if (res.is_ok()) return ok(res.value());
        }
        if (op == BuiltinOp::Rank) {
            auto res = cas::linalg::rank(m_expr, *context_);
            if (res.is_ok()) return ok(make_integer(arena_,
                static_cast<long long>(res.value())));
        }
        if (op == BuiltinOp::Trace) {
            auto res = cas::linalg::trace(m_expr, *context_);
            if (res.is_ok()) return ok(res.value());
        }
        if (op == BuiltinOp::Inv) {
            auto res = cas::linalg::inverse(m_expr, *context_);
            if (res.is_ok()) return ok(arena_.make<Matrix>(
                res.value().rows(), res.value().cols(), res.value().elements()));
        }
        if (op == BuiltinOp::Transpose) {
            auto res = cas::linalg::transpose(m_expr);
            if (res.is_ok()) return ok(arena_.make<Matrix>(
                res.value().rows(), res.value().cols(), res.value().elements()));
        }
    }

    // N — numeric evaluation
    if (op == BuiltinOp::N && args.size() == 2U) {
        // N(expr, digits) — MPFR arbitrary-precision evaluation
        LiteralRational drat;
        if (auto ex = try_get_exact_rational(args[1], drat);
            ex.is_ok() && ex.value()
            && !drat.value.numerator().is_negative()
            && !drat.value.numerator().is_zero()) {
            const unsigned int digits = static_cast<unsigned int>(
                drat.value.numerator().to_u64()
                / (drat.value.denominator().is_zero()
                    ? 1ULL : drat.value.denominator().to_u64()));
            auto mpfr_res = cas::numeric::eval_mpfr(args[0], digits);
            if (mpfr_res.is_ok())
                return ok(arena_.make<DecimalLit>(mpfr_res.value()));
        }
    }
    if (op == BuiltinOp::N && args.size() == 1U) {
        if (const auto* matrix = expr_cast<Matrix>(args.front())) {
            std::vector<ExprPtr> numeric_elements;
            for (auto elem : matrix->elements) {
                auto val = cas::numeric::eval(elem);
                if (val.is_ok())
                    numeric_elements.push_back(arena_.make<DecimalLit>(val.value()));
                else
                    numeric_elements.push_back(elem);
            }
            return ok(arena_.make<Matrix>(
                matrix->rows, matrix->cols, std::move(numeric_elements)));
        }
        auto val = cas::numeric::eval(args.front());
        if (val.is_ok()) return ok(arena_.make<DecimalLit>(val.value()));
    }

    // L2-26 Piecewise: piecewise(cond_1, expr_1, ..., cond_k, expr_k, default).
    // args.size() must be odd: even index = condition, odd index = branch
    // expression, last arg = default fallback. Walk conditions; if exact
    // True (IntegerLit 1) return that branch; if exact False (IntegerLit 0)
    // drop pair; otherwise keep undecided pair. If no decided branch and
    // all dropped, return default.
    if (op == BuiltinOp::Piecewise && args.size() >= 1U && args.size() % 2U == 1U) {
        // Only default branch → return it directly.
        if (args.size() == 1U) return ok(args[0]);
        std::vector<ExprPtr> remaining;
        ExprPtr default_branch = args.back();
        bool any_dropped_or_decided = false;
        for (std::size_t i = 0; i + 1 < args.size(); i += 2) {
            ExprPtr cond = args[i];
            ExprPtr branch = args[i + 1];
            if (const auto* il = expr_cast<IntegerLit>(cond)) {
                any_dropped_or_decided = true;
                if (!il->value.is_zero()) {
                    // Definitely true → branch selected. Pre-pending remaining
                    // pairs ahead of selected branch is impossible (earlier
                    // ones must be evaluated first); they're still undecided.
                    if (remaining.empty()) return ok(branch);
                    // Otherwise: append branch as new default; further pairs
                    // ignored (sequential semantics).
                    remaining.push_back(arena_.make<IntegerLit>(BigInt(1)));
                    remaining.push_back(branch);
                    default_branch = branch;
                    break;
                }
                // false → drop pair
                continue;
            }
            remaining.push_back(cond);
            remaining.push_back(branch);
        }
        if (any_dropped_or_decided) {
            if (remaining.empty()) return ok(default_branch);
            remaining.push_back(default_branch);
            return ok(arena_.make<FuncCall>(BuiltinOp::Piecewise, std::move(remaining)));
        }
    }

    if (expr_ptr_sequence_identical(args, node.args)) return ok(original);
    return ok(arena_.make<FuncCall>(node.func_id, std::move(args)));
}

// ── Other node types ──────────────────────────────────────────────────────────

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Integral& node) {
    return simplify_passthrough(original, node);
}
Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Derivative& node) {
    return simplify_passthrough(original, node);
}
Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Limit& node) {
    return simplify_passthrough(original, node);
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const RootOf& node) {
    auto simplified_poly = simplify_expr(node.polynomial);
    if (simplified_poly.is_error()) return simplified_poly;

    if (context_) {
        auto poly_res = cas::algebra::parse_polynomial(
            simplified_poly.value(), node.variable, *context_);
        if (poly_res.is_ok()
            && cas::algebra::poly_degree(poly_res.value())
               <= context_->max_rootof_explicit_degree()) {
            auto roots_res = cas::algebra::solve_polynomial(
                simplified_poly.value(), node.variable, *context_);
            if (roots_res.is_ok()) {
                const auto& roots = roots_res.value();
                if (roots.size() == 1) {
                    if (roots[0]->kind != ExprKind::RootOf) return ok(roots[0]);
                } else if (node.root_index.has_value()
                           && *node.root_index < roots.size()) {
                    if (roots[*node.root_index]->kind != ExprKind::RootOf)
                        return ok(roots[*node.root_index]);
                }
            }
        }
    }

    if (simplified_poly.value() == node.polynomial) return ok(original);
    return ok(arena_.make<RootOf>(
        simplified_poly.value(), node.variable, node.root_index));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Matrix& node) {
    std::vector<ExprPtr> elements;
    elements.reserve(node.elements.size());
    for (std::size_t i = 0; i < node.elements.size(); ++i) {
        auto simplify_el = [&]() -> Result<ExprPtr> {
            if (trace_enabled_) {
                std::vector<ExprPtr> current = node.elements;
                for (std::size_t j = 0; j < elements.size(); ++j) current[j] = elements[j];
                ScopedFrame frame(*this,
                    [this, current = std::move(current), i,
                     r = node.rows, c = node.cols](ExprPtr value) mutable {
                        current[i] = value;
                        return arena_.make<Matrix>(r, c, std::move(current));
                    });
                return simplify_expr(node.elements[i]);
            }
            return simplify_expr(node.elements[i]);
        };
        auto res = simplify_el();
        if (res.is_error()) return res;
        elements.push_back(res.value());
    }
    if (expr_ptr_sequence_identical(elements, node.elements)) return ok(original);
    return ok(arena_.make<Matrix>(node.rows, node.cols, std::move(elements)));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const SeriesExp& node) {
    std::vector<std::pair<long long, ExprPtr>> simplified_terms;
    simplified_terms.reserve(node.terms.size());
    bool changed = false;

    auto simplified_point = [&]() -> Result<ExprPtr> {
        if (trace_enabled_) {
            ScopedFrame frame(*this, [this, &node](ExprPtr value) {
                return arena_.make<SeriesExp>(node.var, value, node.terms, node.order);
            });
            return simplify_expr(node.point);
        }
        return simplify_expr(node.point);
    }();
    if (simplified_point.is_error()) return simplified_point;
    if (simplified_point.value() != node.point) changed = true;

    for (std::size_t i = 0; i < node.terms.size(); ++i) {
        auto simplified_coeff = [&]() -> Result<ExprPtr> {
            if (trace_enabled_) {
                ScopedFrame frame(*this, [this, &node, i](ExprPtr value) {
                    auto terms = node.terms;
                    terms[i].second = value;
                    return arena_.make<SeriesExp>(
                        node.var, node.point, std::move(terms), node.order);
                });
                return simplify_expr(node.terms[i].second);
            }
            return simplify_expr(node.terms[i].second);
        }();
        if (simplified_coeff.is_error()) return simplified_coeff;
        if (simplified_coeff.value() != node.terms[i].second) changed = true;
        simplified_terms.push_back({node.terms[i].first, simplified_coeff.value()});
    }

    if (!changed) return ok(original);
    return ok(arena_.make<SeriesExp>(
        node.var, simplified_point.value(), std::move(simplified_terms), node.order));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Quantity& node) {
    auto res = simplify_expr(node.value);
    if (res.is_error()) return res;
    if (res.value() == node.value) return ok(original);
    return ok(arena_.make<Quantity>(res.value(), node.dimensions));
}

template <typename Node>
Result<ExprPtr> Simplifier::simplify_passthrough(ExprPtr original, const Node&) {
    return ok(original);
}

// ── Rewrite guards ────────────────────────────────────────────────────────────

bool Simplifier::may_rewrite_function_call(BuiltinOp op,
    const std::vector<ExprPtr>& args) const {
    if (args.size() != 1U) return false;
    if (op == BuiltinOp::Tan) return true;
    if (is_parity_rewrite_function(op)) {
        if (is_zero_expr(args.front()))
            return op == BuiltinOp::Sin || op == BuiltinOp::Cos || op == BuiltinOp::Tan
                || op == BuiltinOp::Sinh || op == BuiltinOp::Tanh;
        const auto* unary = expr_cast<Unary>(args.front());
        return unary != nullptr && unary->op == UnaryOp::Neg;
    }
    if (op == BuiltinOp::Exp)
        return is_zero_expr(args.front()) || is_one_expr(args.front());
    if (op == BuiltinOp::Ln) {
        if (is_one_expr(args.front()) || is_constant_expr(args.front(), MathConstant::E))
            return true;
        if (const auto* quot = expr_cast<Binary>(args.front());
            quot != nullptr && quot->op == BinaryOp::Div)
            return is_known_positive(quot->left) && is_known_positive(quot->right);
        if (const auto* pow = expr_cast<Binary>(args.front());
            pow != nullptr && pow->op == BinaryOp::Pow)
            return is_constant_expr(pow->left, MathConstant::E)
                || is_known_positive(pow->left);
        if (const auto* prod = expr_cast<Product>(args.front()))
            return std::all_of(prod->factors.begin(), prod->factors.end(),
                [this](ExprPtr f) { return is_known_positive(f); });
        const auto* sqrt = expr_cast<FuncCall>(args.front());
        return sqrt != nullptr && sqrt->func_id == BuiltinOp::Sqrt
            && sqrt->args.size() == 1U && is_known_positive(sqrt->args.front());
    }
    if (op != BuiltinOp::Sqrt) {
        if (op == BuiltinOp::Det || op == BuiltinOp::Rank || op == BuiltinOp::Trace
            || op == BuiltinOp::Inv || op == BuiltinOp::Transpose || op == BuiltinOp::N)
            return !args.empty() && (expr_is<Matrix>(args.front()) || op == BuiltinOp::N);
        return false;
    }
    if (const auto* quot = expr_cast<Binary>(args.front());
        quot != nullptr && quot->op == BinaryOp::Div)
        return is_known_nonnegative(quot->left) && is_known_positive(quot->right);
    const auto* prod = expr_cast<Product>(args.front());
    return prod != nullptr && std::all_of(prod->factors.begin(), prod->factors.end(),
        [this](ExprPtr f) { return is_known_nonnegative(f); });
}

bool Simplifier::may_rewrite_sum_terms(const std::vector<ExprPtr>& terms) const {
    ExprPtr s_arg, c_arg;
    for (ExprPtr term : terms) {
        const auto* pow = expr_cast<Binary>(term);
        if (pow == nullptr || pow->op != BinaryOp::Pow) continue;
        if (auto exp = try_get_integer_exponent(pow->right);
            !exp.has_value() || *exp != BigInt(2)) continue;
        const auto* call = expr_cast<FuncCall>(pow->left);
        if (call == nullptr || call->args.size() != 1U) continue;
        if (call->func_id == BuiltinOp::Sin) {
            if (c_arg && call->args.front() == c_arg) return true;
            if (!s_arg) s_arg = call->args.front();
        } else if (call->func_id == BuiltinOp::Cos) {
            if (s_arg && call->args.front() == s_arg) return true;
            if (!c_arg) c_arg = call->args.front();
        }
    }
    return false;
}

bool Simplifier::may_rewrite_power(ExprPtr base, ExprPtr exponent) const {
    if (!is_constant_expr(base, MathConstant::E)) return false;
    const auto* ln = expr_cast<FuncCall>(exponent);
    return ln != nullptr && ln->func_id == BuiltinOp::Ln
        && ln->args.size() == 1U && is_known_positive(ln->args.front());
}

// ── Sign predicates ───────────────────────────────────────────────────────────
// is_known_positive / _nonnegative / _negative / is_assumed_nonzero are defined
// in simplify_sign_predicates.cpp (anti-monolith split 2026-06-19).

} // namespace cas::symbolic::detail
