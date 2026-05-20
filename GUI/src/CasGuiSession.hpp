#pragma once

#include "cas/expr.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::gui {

struct PlotSample {
    double x{};
    double y{};
};

struct ComputeResult {
    bool ok{false};
    bool interrupted{false};
    std::string error;
    std::string text;
    std::string latex;
    std::string ascii;
    std::optional<double> numeric_value;
    struct Representation {
        std::string id;
        std::string label;
        std::string value;
    };
    struct Step {
        std::uint16_t rule_id{0};
        std::uint8_t depth{0};
        std::string before_latex;
        std::string after_latex;
        std::string root_latex;
    };
    std::vector<Representation> representations;
    std::vector<Step> steps;
};

class CasGuiSession {
public:
    struct StoredDefinition {
        std::string name;
        std::string value_text;
    };

    [[nodiscard]] ComputeResult simplify(std::string_view input);

    [[nodiscard]] Result<std::vector<PlotSample>> sample_2d(
        std::string_view input,
        std::string variable,
        double x_min,
        double x_max);

    [[nodiscard]] std::vector<std::string> list_functions() const;
    [[nodiscard]] std::vector<std::pair<std::string, std::string>> list_variables() const;
    [[nodiscard]] std::vector<StoredDefinition> snapshot_definitions() const;
    [[nodiscard]] Result<void> restore_definitions(const std::vector<StoredDefinition>& definitions);
    void clear_definitions();
    void interrupt() { context_.interrupt(); }

private:
    [[nodiscard]] Result<ExprPtr> parse(std::string_view input);
    [[nodiscard]] Result<ExprPtr> resolve_definitions(ExprPtr expr);

    symbolic::CASContext context_;
};

} // namespace cas::gui
