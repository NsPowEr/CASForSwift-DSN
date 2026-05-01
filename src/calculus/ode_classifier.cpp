#include "cas/ode.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"

namespace cas::calculus {

[[nodiscard]] static CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] Result<OdeClassification> classify_ode(
    ExprPtr equation,
    const Symbol& y,
    const Symbol& x,
    [[maybe_unused]] symbolic::CASContext& ctx) {
    
    return ok(OdeClassification(OdeType::Unknown, equation, y, x));
}

[[nodiscard]] Result<ExprPtr> solve_ode(ExprPtr equation, const Symbol& y, const Symbol& x, symbolic::CASContext& ctx) {
    auto class_res = classify_ode(equation, y, x, ctx);
    if (class_res.is_error()) return fail<ExprPtr>(class_res.error());
    
    const auto& classification = class_res.value();
    
    switch (classification.type) {
        case OdeType::Separable:
        case OdeType::Linear1stOrder:
        case OdeType::Bernoulli:
        case OdeType::Exact:
            return solve_ode_1st_order(classification, ctx);
            
        case OdeType::Linear2ndOrderConstantCoeff:
        case OdeType::Linear2ndOrderRationalCoeff:
            return solve_ode_advanced(classification, ctx);
            
        default:
            return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Tipo di ODE non riconosciuto o non supportato analiticamente."));
    }
}

} // namespace cas::calculus
