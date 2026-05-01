#include "CasGuiSession.hpp"

#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/numeric/sampler.hpp"
#include "cas/parser.hpp"

namespace cas::gui {

Result<ExprPtr> CasGuiSession::parse(std::string_view input) {
    auto tokens = Lexer(std::string(input)).tokenize();
    if (tokens.is_error()) {
        return fail<ExprPtr>(tokens.error());
    }

    Parser parser(tokens.value(), context_.arena());
    return parser.parse();
}

ComputeResult CasGuiSession::simplify(std::string_view input) {
    ComputeResult output;

    auto parsed = parse(input);
    if (parsed.is_error()) {
        output.error = parsed.error().message;
        return output;
    }

    auto simplified = context_.simplify(parsed.value());
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

    auto numeric = numeric::eval(expr);
    if (numeric.is_ok()) {
        output.numeric_value = numeric.value();
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

    numeric::AdaptiveSampler sampler(parsed.value(), std::move(variable));
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

} // namespace cas::gui
