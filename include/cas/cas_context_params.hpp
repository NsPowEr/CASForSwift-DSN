#pragma once

// cas_context_params.hpp — Configurable algorithm parameter storage for CASContext.
//
// Defines CASContextParams: a base struct that carries every tuneable knob
// (thresholds, budgets, quality parameters) exposed by CASContext.
//
// Refactored to inherit from:
// - CASContextSimplifierParams (simplifier parameters)
// - CASContextCalculusParams (calculus & integration parameters)
// - CASContextAlgebraParams (algebra & arithmetic parameters)
// keeping each file under 500 lines per Mandato di Guardia Architetturale.

#include "cas_context_simplifier_params.hpp"
#include "cas_context_calculus_params.hpp"
#include "cas_context_algebra_params.hpp"

namespace cas::symbolic {

struct CASContextParams :
    public CASContextSimplifierParams,
    public CASContextCalculusParams,
    public CASContextAlgebraParams
{
};

}  // namespace cas::symbolic
