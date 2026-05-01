#pragma once

#include "cas/expr.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cas::gui {

struct PlotSample {
    double x{};
    double y{};
};

struct ComputeResult {
    bool ok{false};
    std::string error;
    std::string text;
    std::string latex;
    std::string ascii;
    std::optional<double> numeric_value;
};

class CasGuiSession {
public:
    [[nodiscard]] ComputeResult simplify(std::string_view input);

    [[nodiscard]] Result<std::vector<PlotSample>> sample_2d(
        std::string_view input,
        std::string variable,
        double x_min,
        double x_max);

private:
    [[nodiscard]] Result<ExprPtr> parse(std::string_view input);

    symbolic::CASContext context_;
};

} // namespace cas::gui
