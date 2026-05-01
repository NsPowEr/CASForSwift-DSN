#pragma once

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/result.hpp"
#include <vector>
#include <memory>

namespace cas::calculus {

/// Types of elementary transcendental extensions
enum class ExtensionType {
    Base,       // The base variable, typically x, where D(x) = 1
    Logarithmic,// t' = u'/u for some u in the base field
    Exponential // t' = u' * t for some u in the base field
};

/// A single step in the differential extension tower
struct DifferentialExtension {
    ExtensionType type;
    ExprPtr argument; // The inner expression (u)
    Symbol t_var;     // The new variable (t) representing log(u) or exp(u)
};

/// Represents a differential field K(t_1, ..., t_n)
class DifferentialField {
public:
    explicit DifferentialField(const Symbol& base_var) : base_var_(base_var) {}

    /// Construct a tower from an expression, finding all logs and exps
    static Result<DifferentialField> build(ExprPtr expr, const Symbol& x, symbolic::CASContext& ctx);

    const Symbol& base_var() const { return base_var_; }
    const std::vector<DifferentialExtension>& extensions() const { return extensions_; }
    
    /// Compute the formal derivative D(expr) in the field
    Result<ExprPtr> derive(ExprPtr expr, symbolic::CASContext& ctx) const;

    /// Replace original functions (log/exp) with generator variables (t_i)
    Result<ExprPtr> to_field_generators(ExprPtr expr, symbolic::CASContext& ctx) const;

    /// Substitute generator variables back to original expressions
    Result<ExprPtr> from_field_generators(ExprPtr expr, symbolic::CASContext& ctx) const;

private:
    Symbol base_var_;
    std::vector<DifferentialExtension> extensions_;

    Result<void> add_extension(ExprPtr expr, symbolic::CASContext& ctx);
};

struct HermiteReduction {
    ExprPtr rational_part; // A/B
    ExprPtr remaining_P;   // C
    ExprPtr remaining_Q;   // D (square-free)
};

/// Performs Hermite reduction on P/Q with respect to generator t_var
Result<HermiteReduction> hermite_reduce(
    ExprPtr P, ExprPtr Q, const Symbol& t_var, const DifferentialField& field, symbolic::CASContext& ctx);

/// Evaluates the logarithmic part of the integral using Rothstein-Trager
Result<ExprPtr> integrate_rothstein_trager(
    ExprPtr P, ExprPtr Q, const Symbol& t_var, const DifferentialField& field, symbolic::CASContext& ctx);

} // namespace cas::calculus
