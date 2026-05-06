#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <chrono>
#include <functional>
#include <vector>

namespace cas::symbolic {

class Substituter {
public:
    Substituter(CASContext& context, const Symbol& variable, ExprPtr value)
        : context_(context),
          variable_(variable),
          value_(value) {}

    [[nodiscard]] Result<ExprPtr> substitute_expr(ExprPtr expr) {
        if (!expr) {
            return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot substitute inside null expression"));
        }

        ++context_.ops_count_;
        if ((context_.ops_count_ % context_.timeout_check_interval()) == 0U &&
            std::chrono::steady_clock::now() - context_.operation_started_at_ >= context_.timeout_) {
            return fail<ExprPtr>(make_error(CASErrorKind::Timeout, "Symbolic operation timed out"));
        }

        return visit_expr(
            expr,
            [this, expr](const auto& node) -> Result<ExprPtr> {
                return substitute_node(expr, node);
            });
    }

private:
    class ScopedFrame {
    public:
        ScopedFrame(Substituter& owner, std::function<ExprPtr(ExprPtr)> builder)
            : owner_(owner),
              active_(owner.context_.trace_capture_active_) {
            if (active_) {
                owner_.root_frames_.push_back(std::move(builder));
            }
        }

        ~ScopedFrame() {
            if (active_) {
                owner_.root_frames_.pop_back();
            }
        }

    private:
        Substituter& owner_;
        bool active_{false};
    };

    [[nodiscard]] ExprPtr build_root_after(ExprPtr target_after) {
        ExprPtr root = target_after;
        for (auto it = root_frames_.rbegin(); it != root_frames_.rend(); ++it) {
            root = (*it)(root);
        }
        return root;
    }

    void append_trace(RuleId rule_id, ExprPtr before, ExprPtr after) {
        if (!context_.trace_capture_active_) {
            return;
        }

        context_.trace_.push_back(TraceStep{
            .rule_id = rule_id,
            .depth = static_cast<std::uint8_t>(std::min<std::size_t>(root_frames_.size(), 255U)),
            .target_before = before,
            .target_after = after,
            .root_after = build_root_after(after),
        });
    }

    [[nodiscard]] Result<ExprPtr> clone_value() {
        return materialize_expr(value_, context_.arena());
    }

    template <typename Node>
    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Node&) {
        return ok(original);
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Symbol& node) {
        if (node.name != variable_.name) {
            return ok(original);
        }

        auto cloned = clone_value();
        if (cloned.is_error()) {
            return cloned;
        }
        append_trace(RuleId::SubstituteSymbol, original, cloned.value());
        return cloned;
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const IntegerLit&) {
        return ok(original);
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const RationalLit&) {
        return ok(original);
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const DecimalLit&) {
        return ok(original);
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Constant&) {
        return ok(original);
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr, const ExprNode&) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot substitute inside null expression node"));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Unary& node) {
        ScopedFrame frame(*this, [this, op = node.op](ExprPtr operand) {
            return context_.arena().make<Unary>(op, operand);
        });
        auto operand = substitute_expr(node.operand);
        if (operand.is_error()) {
            return operand;
        }
        if (operand.value() == node.operand) {
            return ok(original);
        }
        return ok(context_.arena().make<Unary>(node.op, operand.value()));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Binary& node) {
        Result<ExprPtr> left = fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
        {
            ScopedFrame frame(*this, [this, op = node.op, right = node.right](ExprPtr left_value) {
                return context_.arena().make<Binary>(op, left_value, right);
            });
            left = substitute_expr(node.left);
            if (left.is_error()) {
                return left;
            }
        }
        Result<ExprPtr> right = fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
        {
            ScopedFrame frame(*this, [this, op = node.op, left_value = left.value()](ExprPtr right_value) {
                return context_.arena().make<Binary>(op, left_value, right_value);
            });
            right = substitute_expr(node.right);
            if (right.is_error()) {
                return right;
            }
        }

        if (left.value() == node.left && right.value() == node.right) {
            return ok(original);
        }
        return ok(context_.arena().make<Binary>(node.op, left.value(), right.value()));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const FuncCall& node) {
        std::vector<ExprPtr> args;
        args.reserve(node.args.size());
        for (std::size_t index = 0; index < node.args.size(); ++index) {
            std::vector<ExprPtr> current_args = node.args;
            for (std::size_t rewritten = 0; rewritten < args.size(); ++rewritten) {
                current_args[rewritten] = args[rewritten];
            }
            ScopedFrame frame(*this, [this, current_args = std::move(current_args), index, name = node.name](ExprPtr value) mutable {
                current_args[index] = value;
                return context_.arena().make<FuncCall>(name, std::move(current_args));
            });
            auto substituted = substitute_expr(node.args[index]);
            if (substituted.is_error()) {
                return substituted;
            }
            args.push_back(substituted.value());
        }

        if (std::equal(args.begin(), args.end(), node.args.begin())) {
            return ok(original);
        }
        return ok(context_.arena().make<FuncCall>(node.name, std::move(args)));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Sum& node) {
        std::vector<ExprPtr> terms;
        terms.reserve(node.terms.size());
        for (std::size_t index = 0; index < node.terms.size(); ++index) {
            std::vector<ExprPtr> current_terms = node.terms;
            for (std::size_t rewritten = 0; rewritten < terms.size(); ++rewritten) {
                current_terms[rewritten] = terms[rewritten];
            }
            ScopedFrame frame(*this, [this, current_terms = std::move(current_terms), index](ExprPtr value) mutable {
                current_terms[index] = value;
                return context_.arena().make<Sum>(std::move(current_terms));
            });
            auto substituted = substitute_expr(node.terms[index]);
            if (substituted.is_error()) {
                return substituted;
            }
            terms.push_back(substituted.value());
        }

        if (std::equal(terms.begin(), terms.end(), node.terms.begin())) {
            return ok(original);
        }
        return ok(context_.arena().make<Sum>(std::move(terms)));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Product& node) {
        std::vector<ExprPtr> factors;
        factors.reserve(node.factors.size());
        for (std::size_t index = 0; index < node.factors.size(); ++index) {
            std::vector<ExprPtr> current_factors = node.factors;
            for (std::size_t rewritten = 0; rewritten < factors.size(); ++rewritten) {
                current_factors[rewritten] = factors[rewritten];
            }
            ScopedFrame frame(*this, [this, current_factors = std::move(current_factors), index](ExprPtr value) mutable {
                current_factors[index] = value;
                return context_.arena().make<Product>(std::move(current_factors));
            });
            auto substituted = substitute_expr(node.factors[index]);
            if (substituted.is_error()) {
                return substituted;
            }
            factors.push_back(substituted.value());
        }

        if (std::equal(factors.begin(), factors.end(), node.factors.begin())) {
            return ok(original);
        }
        return ok(context_.arena().make<Product>(std::move(factors)));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Integral& node) {
        const bool shadows_variable = node.variable.name == variable_.name;
        ExprPtr integrand = ExprPtr{};
        if (shadows_variable) {
            integrand = node.integrand;
        } else {
            ScopedFrame frame(*this, [this, variable = node.variable, lower = node.lower, upper = node.upper](ExprPtr value) {
                return context_.arena().make<Integral>(value, variable, lower, upper);
            });
            auto substituted = substitute_expr(node.integrand);
            if (substituted.is_error()) {
                return substituted;
            }
            integrand = substituted.value();
        }

        std::optional<ExprPtr> lower;
        if (node.lower.has_value()) {
            ScopedFrame frame(*this, [this, integrand, variable = node.variable, upper = node.upper](ExprPtr value) {
                return context_.arena().make<Integral>(integrand, variable, value, upper);
            });
            auto substituted = substitute_expr(*node.lower);
            if (substituted.is_error()) {
                return substituted;
            }
            lower = substituted.value();
        }

        std::optional<ExprPtr> upper;
        if (node.upper.has_value()) {
            ScopedFrame frame(*this, [this, integrand, variable = node.variable, lower](ExprPtr value) {
                return context_.arena().make<Integral>(integrand, variable, lower, value);
            });
            auto substituted = substitute_expr(*node.upper);
            if (substituted.is_error()) {
                return substituted;
            }
            upper = substituted.value();
        }

        if (integrand == node.integrand && lower == node.lower && upper == node.upper) {
            return ok(original);
        }
        return ok(context_.arena().make<Integral>(integrand, node.variable, lower, upper));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Derivative& node) {
        const bool shadows_variable = node.variable.name == variable_.name;
        if (shadows_variable) {
            return ok(original);
        }

        ScopedFrame frame(*this, [this, variable = node.variable, order = node.order](ExprPtr value) {
            return context_.arena().make<Derivative>(value, variable, order);
        });
        auto substituted = substitute_expr(node.expression);
        if (substituted.is_error()) {
            return substituted;
        }

        if (substituted.value() == node.expression) {
            return ok(original);
        }
        return ok(context_.arena().make<Derivative>(substituted.value(), node.variable, node.order));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Limit& node) {
        const bool shadows_variable = node.variable.name == variable_.name;
        ExprPtr expression = ExprPtr{};
        if (shadows_variable) {
            expression = node.expression;
        } else {
            ScopedFrame frame(*this, [this, variable = node.variable, point = node.point, direction = node.direction](ExprPtr value) {
                return context_.arena().make<Limit>(value, variable, point, direction);
            });
            auto substituted = substitute_expr(node.expression);
            if (substituted.is_error()) {
                return substituted;
            }
            expression = substituted.value();
        }

        ScopedFrame frame(*this, [this, expression, variable = node.variable, direction = node.direction](ExprPtr value) {
            return context_.arena().make<Limit>(expression, variable, value, direction);
        });
        auto substituted_point = substitute_expr(node.point);
        if (substituted_point.is_error()) {
            return substituted_point;
        }

        if (expression == node.expression && substituted_point.value() == node.point) {
            return ok(original);
        }
        return ok(context_.arena().make<Limit>(expression, node.variable, substituted_point.value(), node.direction));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const RootOf& node) {
        const bool shadows_variable = node.variable.name == variable_.name;
        if (shadows_variable) {
            return ok(original);
        }

        ScopedFrame frame(*this, [this, variable = node.variable, index = node.root_index](ExprPtr value) {
            return context_.arena().make<RootOf>(value, variable, index);
        });
        auto substituted = substitute_expr(node.polynomial);
        if (substituted.is_error()) {
            return substituted;
        }

        if (substituted.value() == node.polynomial) {
            return ok(original);
        }
        return ok(context_.arena().make<RootOf>(substituted.value(), node.variable, node.root_index));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const Matrix& node) {
        std::vector<ExprPtr> elements;
        elements.reserve(node.elements.size());
        for (std::size_t index = 0; index < node.elements.size(); ++index) {
            std::vector<ExprPtr> current_elements = node.elements;
            for (std::size_t rewritten = 0; rewritten < elements.size(); ++rewritten) {
                current_elements[rewritten] = elements[rewritten];
            }
            ScopedFrame frame(*this, [this, current_elements = std::move(current_elements), index, rows = node.rows, cols = node.cols](ExprPtr value) mutable {
                current_elements[index] = value;
                return context_.arena().make<Matrix>(rows, cols, std::move(current_elements));
            });
            auto substituted = substitute_expr(node.elements[index]);
            if (substituted.is_error()) {
                return substituted;
            }
            elements.push_back(substituted.value());
        }

        if (std::equal(elements.begin(), elements.end(), node.elements.begin())) {
            return ok(original);
        }
        return ok(context_.arena().make<Matrix>(node.rows, node.cols, std::move(elements)));
    }

    [[nodiscard]] Result<ExprPtr> substitute_node(ExprPtr original, const SeriesExp& node) {
        const bool shadows_variable = node.var.name == variable_.name;
        
        // Point is ALWAYS substituted
        auto substituted_point = substitute_expr(node.point);
        if (substituted_point.is_error()) return substituted_point;

        if (shadows_variable) {
            if (substituted_point.value() == node.point) return ok(original);
            return ok(context_.arena().make<SeriesExp>(node.var, substituted_point.value(), node.terms, node.order));
        }

        // Coefficients are substituted
        std::vector<std::pair<long long, ExprPtr>> terms;
        terms.reserve(node.terms.size());
        bool changed = (substituted_point.value() != node.point);
        for (std::size_t i = 0; i < node.terms.size(); ++i) {
            auto substituted_coeff = substitute_expr(node.terms[i].second);
            if (substituted_coeff.is_error()) return substituted_coeff;
            if (substituted_coeff.value() != node.terms[i].second) changed = true;
            terms.push_back({node.terms[i].first, substituted_coeff.value()});
        }

        if (!changed) return ok(original);
        return ok(context_.arena().make<SeriesExp>(node.var, substituted_point.value(), std::move(terms), node.order));
    }

    CASContext& context_;
    const Symbol& variable_;
    ExprPtr value_;
    std::vector<std::function<ExprPtr(ExprPtr)>> root_frames_;
};

Result<ExprPtr> CASContext::substitute(ExprPtr expr, const Symbol& variable, ExprPtr value) {
    const bool owns_operation = !operation_active_;
    if (owns_operation) {
        operation_active_ = true;
        trace_capture_active_ = trace_enabled_;
        trace_.clear();
        ops_count_ = 0;
        operation_started_at_ = std::chrono::steady_clock::now();
    }
    auto result = symbolic::substitute(expr, variable, value, *this);
    if (owns_operation) {
        operation_active_ = false;
        trace_capture_active_ = false;
        ops_count_ = 0;
    }
    return result;
}

Result<ExprPtr> substitute(ExprPtr expr, const Symbol& variable, ExprPtr value, CASContext& context) {
    Substituter substituter(context, variable, value);
    auto substituted = substituter.substitute_expr(expr);
    if (substituted.is_error()) {
        return substituted;
    }
    return context.simplify(substituted.value());
}

} // namespace cas::symbolic
