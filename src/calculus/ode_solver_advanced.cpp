#include "cas/ode.hpp"
#include "calculus_internal.hpp"
#include "cas/algebra.hpp"

namespace cas::calculus {

[[nodiscard]] static CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] Result<ExprPtr> solve_ode_advanced(const OdeClassification& classification, [[maybe_unused]] symbolic::CASContext& ctx) {
    if (classification.type == OdeType::Linear2ndOrderConstantCoeff) {
        // y'' + ay' + by = f(x)
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "ODE lineare 2° ordine a coefficienti costanti non ancora implementata."));
    }
    
    if (classification.type == OdeType::Linear2ndOrderRationalCoeff) {
        // Algoritmo di Kovacic
        if (classification.components.size() < 2) {
             return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Componenti p(x) o q(x) mancanti per Kovacic."));
        }
        
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Algoritmo di Kovacic in fase di sviluppo. Nessuna soluzione Liouvilliana trovata per ora."));
    }
    
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Metodo avanzato ODE non implementato."));
}

} // namespace cas::calculus
