#include "CasGuiSession.hpp"

#include "cas/builtin_functions.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/numeric/sampler.hpp"
#include "cas/parser.hpp"
#include "cas/trace.hpp"

#include <algorithm>
#include <unordered_set>
#include <sstream>
#include <type_traits>

namespace cas::gui {
namespace {

[[nodiscard]] std::string format_latex_or_empty(formatter::LaTeXFormatter& formatter, ExprPtr expr) {
    return expr ? formatter.format(expr) : std::string{};
}

[[nodiscard]] std::string format_double(double value) {
    std::ostringstream stream;
    stream.precision(17);
    stream << value;
    return stream.str();
}

[[nodiscard]] CASError gui_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] Result<ExprPtr> resolve_definitions_impl(
    ExprPtr expr,
    symbolic::CASContext& context,
    std::unordered_set<std::string> blocked_names) {
    if (!expr) {
        return fail<ExprPtr>(gui_error(CASErrorKind::InvalidArgument, "Cannot resolve null expression"));
    }

    return visit_expr(
        expr,
        [&](const auto& node) -> Result<ExprPtr> {
            using Node = std::decay_t<decltype(node)>;

            if constexpr (std::is_same_v<Node, Symbol>) {
                if (blocked_names.contains(node.name)) {
                    return ok(expr);
                }
                const auto resolved = context.lookup(node);
                if (!resolved.has_value()) {
                    return ok(expr);
                }
                auto nested_blocked = blocked_names;
                nested_blocked.insert(node.name);
                return resolve_definitions_impl(*resolved, context, std::move(nested_blocked));
            } else if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, Constant>) {
                return ok(expr);
            } else if constexpr (std::is_same_v<Node, Unary>) {
                auto operand = resolve_definitions_impl(node.operand, context, blocked_names);
                if (operand.is_error()) {
                    return operand;
                }
                if (operand.value() == node.operand) {
                    return ok(expr);
                }
                return ok(context.arena().make<Unary>(node.op, operand.value()));
            } else if constexpr (std::is_same_v<Node, Binary>) {
                auto left = resolve_definitions_impl(node.left, context, blocked_names);
                if (left.is_error()) {
                    return left;
                }
                auto right = resolve_definitions_impl(node.right, context, blocked_names);
                if (right.is_error()) {
                    return right;
                }
                if (left.value() == node.left && right.value() == node.right) {
                    return ok(expr);
                }
                return ok(context.arena().make<Binary>(node.op, left.value(), right.value()));
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> args;
                bool changed = false;
                args.reserve(node.args.size());
                for (ExprPtr arg : node.args) {
                    auto rewritten = resolve_definitions_impl(arg, context, blocked_names);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || (rewritten.value() != arg);
                    args.push_back(rewritten.value());
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(context.arena().make<FuncCall>(node.name, std::move(args)));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> terms;
                bool changed = false;
                terms.reserve(node.terms.size());
                for (ExprPtr term : node.terms) {
                    auto rewritten = resolve_definitions_impl(term, context, blocked_names);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || (rewritten.value() != term);
                    terms.push_back(rewritten.value());
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(context.arena().make<Sum>(std::move(terms)));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> factors;
                bool changed = false;
                factors.reserve(node.factors.size());
                for (ExprPtr factor : node.factors) {
                    auto rewritten = resolve_definitions_impl(factor, context, blocked_names);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || (rewritten.value() != factor);
                    factors.push_back(rewritten.value());
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(context.arena().make<Product>(std::move(factors)));
            } else if constexpr (std::is_same_v<Node, Integral>) {
                auto inner_blocked = blocked_names;
                inner_blocked.insert(node.variable.name);
                auto integrand = resolve_definitions_impl(node.integrand, context, std::move(inner_blocked));
                if (integrand.is_error()) {
                    return integrand;
                }
                std::optional<ExprPtr> lower = std::nullopt;
                std::optional<ExprPtr> upper = std::nullopt;
                bool changed = integrand.value() != node.integrand;
                if (node.lower.has_value()) {
                    auto rewritten = resolve_definitions_impl(*node.lower, context, blocked_names);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    lower = rewritten.value();
                    changed = changed || (rewritten.value() != *node.lower);
                }
                if (node.upper.has_value()) {
                    auto rewritten = resolve_definitions_impl(*node.upper, context, blocked_names);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    upper = rewritten.value();
                    changed = changed || (rewritten.value() != *node.upper);
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(context.arena().make<Integral>(integrand.value(), node.variable, lower, upper));
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                auto inner_blocked = blocked_names;
                inner_blocked.insert(node.variable.name);
                auto expression = resolve_definitions_impl(node.expression, context, std::move(inner_blocked));
                if (expression.is_error()) {
                    return expression;
                }
                if (expression.value() == node.expression) {
                    return ok(expr);
                }
                return ok(context.arena().make<Derivative>(expression.value(), node.variable, node.order));
            } else if constexpr (std::is_same_v<Node, Limit>) {
                auto inner_blocked = blocked_names;
                inner_blocked.insert(node.variable.name);
                auto expression = resolve_definitions_impl(node.expression, context, std::move(inner_blocked));
                if (expression.is_error()) {
                    return expression;
                }
                auto point = resolve_definitions_impl(node.point, context, blocked_names);
                if (point.is_error()) {
                    return point;
                }
                if (expression.value() == node.expression && point.value() == node.point) {
                    return ok(expr);
                }
                return ok(context.arena().make<Limit>(expression.value(), node.variable, point.value(), node.direction));
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                auto inner_blocked = blocked_names;
                inner_blocked.insert(node.variable.name);
                auto polynomial = resolve_definitions_impl(node.polynomial, context, std::move(inner_blocked));
                if (polynomial.is_error()) {
                    return polynomial;
                }
                if (polynomial.value() == node.polynomial) {
                    return ok(expr);
                }
                return ok(context.arena().make<RootOf>(polynomial.value(), node.variable, node.root_index));
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::vector<ExprPtr> elements;
                bool changed = false;
                elements.reserve(node.elements.size());
                for (ExprPtr element : node.elements) {
                    auto rewritten = resolve_definitions_impl(element, context, blocked_names);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || (rewritten.value() != element);
                    elements.push_back(rewritten.value());
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(context.arena().make<Matrix>(node.rows, node.cols, std::move(elements)));
            } else if constexpr (std::is_same_v<Node, SeriesExp>) {
                auto point = resolve_definitions_impl(node.point, context, blocked_names);
                if (point.is_error()) {
                    return point;
                }
                auto inner_blocked = blocked_names;
                inner_blocked.insert(node.var.name);
                bool changed = point.value() != node.point;
                std::vector<std::pair<long long, ExprPtr>> terms;
                terms.reserve(node.terms.size());
                for (const auto& [power, coeff] : node.terms) {
                    auto rewritten = resolve_definitions_impl(coeff, context, inner_blocked);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || (rewritten.value() != coeff);
                    terms.push_back({power, rewritten.value()});
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(context.arena().make<SeriesExp>(node.var, point.value(), std::move(terms), node.order));
            } else if constexpr (std::is_same_v<Node, Quantity>) {
                auto value = resolve_definitions_impl(node.value, context, blocked_names);
                if (value.is_error()) {
                    return value;
                }
                if (value.value() == node.value) {
                    return ok(expr);
                }
                return ok(context.arena().make<Quantity>(value.value(), node.dimensions));
            } else {
                return ok(expr);
            }
        });
}

} // namespace

Result<ExprPtr> CasGuiSession::parse(std::string_view input) {
    const std::string owned_input(input);
    Lexer lexer(owned_input);
    auto tokens = lexer.tokenize();
    if (tokens.is_error()) {
        return fail<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), context_.arena());
    return parser.parse();
}

Result<ExprPtr> CasGuiSession::resolve_definitions(ExprPtr expr) {
    return resolve_definitions_impl(expr, context_, {});
}

ComputeResult CasGuiSession::simplify(std::string_view input) {
    ComputeResult output;

    auto parsed = parse(input);
    if (parsed.is_error()) {
        output.error = parsed.error().message;
        return output;
    }

    ExprPtr expression = parsed.value();
    if (const auto* equality = expr_cast<Binary>(expression);
        equality != nullptr && equality->op == BinaryOp::Equal) {
        if (const auto* symbol = expr_cast<Symbol>(equality->left)) {
            auto resolved_rhs = resolve_definitions_impl(
                equality->right,
                context_,
                std::unordered_set<std::string>{symbol->name});
            if (resolved_rhs.is_error()) {
                output.error = resolved_rhs.error().message;
                return output;
            }
            context_.enable_trace(true);
            auto simplified_rhs = context_.simplify(resolved_rhs.value());
            const auto trace = context_.get_trace();
            context_.enable_trace(false);
            if (simplified_rhs.is_error()) {
                output.error = simplified_rhs.error().message;
                return output;
            }
            context_.define(*symbol, simplified_rhs.value());
            expression = context_.arena().make<Binary>(BinaryOp::Equal, equality->left, simplified_rhs.value());

            formatter::TextFormatter text_formatter;
            formatter::LaTeXFormatter latex_formatter;
            formatter::Ascii2DFormatter ascii_formatter;
            output.ok = true;
            output.text = text_formatter.format(expression);
            output.latex = latex_formatter.format(expression);
            output.ascii = ascii_formatter.format(expression);
            output.representations.push_back({"text", "Text", output.text});
            output.representations.push_back({"ascii", "ASCII", output.ascii});
            output.representations.push_back({"binding", "Binding", symbol->name});
            output.steps.reserve(trace.size());
            for (const auto& step : trace) {
                output.steps.push_back(ComputeResult::Step{
                    static_cast<std::uint16_t>(step.rule_id),
                    step.depth,
                    format_latex_or_empty(latex_formatter, step.target_before),
                    format_latex_or_empty(latex_formatter, step.target_after),
                    format_latex_or_empty(latex_formatter, step.root_after),
                });
            }
            return output;
        }
    }

    auto resolved = resolve_definitions(expression);
    if (resolved.is_error()) {
        output.error = resolved.error().message;
        return output;
    }

    context_.enable_trace(true);
    auto simplified = context_.simplify(resolved.value());
    const auto trace = context_.get_trace();
    context_.enable_trace(false);
    if (simplified.is_error()) {
        output.error = simplified.error().message;
        return output;
    }

    formatter::TextFormatter text_formatter;
    formatter::LaTeXFormatter latex_formatter;
    formatter::Ascii2DFormatter ascii_formatter;

    const ExprPtr expr = simplified.value();
    output.ok = true;
    output.text = text_formatter.format(expr);
    output.latex = latex_formatter.format(expr);
    output.ascii = ascii_formatter.format(expr);
    output.representations.push_back({"text", "Text", output.text});
    output.representations.push_back({"ascii", "ASCII", output.ascii});

    auto numeric = numeric::eval(expr);
    if (numeric.is_ok()) {
        output.numeric_value = numeric.value();
        output.representations.push_back(
            {"numeric", "Numeric", format_double(output.numeric_value.value())});
    }

    output.steps.reserve(trace.size());
    for (const auto& step : trace) {
        output.steps.push_back(ComputeResult::Step{
            static_cast<std::uint16_t>(step.rule_id),
            step.depth,
            format_latex_or_empty(latex_formatter, step.target_before),
            format_latex_or_empty(latex_formatter, step.target_after),
            format_latex_or_empty(latex_formatter, step.root_after),
        });
    }

    return output;
}

Result<std::vector<PlotSample>> CasGuiSession::sample_2d(
    std::string_view input,
    std::string variable,
    double x_min,
    double x_max) {
    auto parsed = parse(input);
    if (parsed.is_error()) {
        return fail<std::vector<PlotSample>>(parsed.error());
    }

    auto resolved = resolve_definitions_impl(
        parsed.value(),
        context_,
        std::unordered_set<std::string>{variable});
    if (resolved.is_error()) {
        return fail<std::vector<PlotSample>>(resolved.error());
    }

    numeric::AdaptiveSampler sampler(resolved.value(), std::move(variable));
    auto sampled = sampler.sample(x_min, x_max);
    if (sampled.is_error()) {
        return fail<std::vector<PlotSample>>(sampled.error());
    }

    std::vector<PlotSample> points;
    points.reserve(sampled.value().size());
    for (const auto& point : sampled.value()) {
        points.push_back({point.x, point.y});
    }

    return ok(std::move(points));
}

std::vector<std::string> CasGuiSession::list_functions() const {
    std::vector<std::string> funcs;
    // Builtin operations from the core enum
    for (int i = 1; i <= static_cast<int>(BuiltinOp::HermiteHe); ++i) {
        auto name = builtin_op_name(static_cast<BuiltinOp>(i));
        if (name != "unknown") {
            funcs.emplace_back(name);
        }
    }

    // Top-level CAS commands not necessarily in BuiltinOp enum but supported by parser/calculus
    const std::vector<std::string> extra = {
        "integrate", "integral", "diff", "derivative", "partial_diff",
        "limit", "lim", "sum", "summation", "taylor", "laurent", "series",
        "solve", "csolve", "factor", "expand", "simplify", "rewrite",
        "collect", "together", "apart", "gcd", "lcm", "resultant", "discriminant",
        "jacobian", "hessian", "gradient", "implicit_diff", "residue",
        "asymptotes", "groebner", "clear", "eval", "N"
    };
    funcs.insert(funcs.end(), extra.begin(), extra.end());

    std::sort(funcs.begin(), funcs.end());
    funcs.erase(std::unique(funcs.begin(), funcs.end()), funcs.end());
    return funcs;
}

std::vector<std::pair<std::string, std::string>> CasGuiSession::list_variables() const {
    std::vector<std::pair<std::string, std::string>> vars;
    formatter::LaTeXFormatter latex_formatter;
    for (const auto& [name, expr] : context_.variables()) {
        vars.push_back({name, latex_formatter.format(expr)});
    }
    std::sort(vars.begin(), vars.end());
    return vars;
}

std::vector<CasGuiSession::StoredDefinition> CasGuiSession::snapshot_definitions() const {
    std::vector<StoredDefinition> definitions;
    definitions.reserve(context_.variables().size());
    for (const auto& [name, expr] : context_.variables()) {
        auto round_trip = to_round_trip_text(expr);
        if (round_trip.is_ok()) {
            definitions.push_back({name, round_trip.value()});
        }
    }
    std::sort(
        definitions.begin(),
        definitions.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.name < rhs.name; });
    return definitions;
}

Result<void> CasGuiSession::restore_definitions(const std::vector<StoredDefinition>& definitions) {
    clear_definitions();
    for (const auto& definition : definitions) {
        auto parsed = parse(definition.value_text);
        if (parsed.is_error()) {
            return fail<void>(parsed.error());
        }
        auto resolved = resolve_definitions_impl(
            parsed.value(),
            context_,
            std::unordered_set<std::string>{definition.name});
        if (resolved.is_error()) {
            return fail<void>(resolved.error());
        }
        auto simplified = context_.simplify(resolved.value());
        if (simplified.is_error()) {
            return fail<void>(simplified.error());
        }
        context_.define(Symbol(definition.name), simplified.value());
    }
    return ok();
}

void CasGuiSession::clear_definitions() {
    context_.clear_variables();
}

} // namespace cas::gui
