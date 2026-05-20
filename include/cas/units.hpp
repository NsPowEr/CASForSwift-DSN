// CAS-L3-08 — SI unit conversion factor registry.
//
// Each registered unit name maps to a pair (scale, SIDimensions):
//   1·cm = 0.01·m       → scale = 1/100, dim = {m=1}
//   1·ft = 0.3048·m     → scale = 3048/10000, dim = {m=1}
//   1·Hz = 1·s^-1       → scale = 1, dim = {s=-1}
//
// All scales stored as Rational (no float losses). Public API:
//
//   make_quantity_from_unit(value, "cm", ctx) → Quantity{value/100, {m=1}}
//   convert_quantity(qty, "ft", ctx)          → Quantity{qty.value·scale, target_dim}
//                                                fails if dims mismatch.

#pragma once

#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"

#include <string>

namespace cas::symbolic { class CASContext; }

namespace cas::units {

// Look up a registered unit by name. Returns nullopt if unknown.
struct UnitInfo {
    Rational scale_to_si;
    SIDimensions dimensions;
};

[[nodiscard]] std::optional<UnitInfo> lookup_unit(const std::string& name) noexcept;

// Wrap a numeric value with a unit, converting scale to SI base.
//   make_quantity_from_unit(5, "cm", ctx) → Quantity(5/100, {m=1})
[[nodiscard]] Result<ExprPtr> make_quantity_from_unit(
    ExprPtr value, const std::string& unit_name, symbolic::CASContext& ctx);

// Convert a Quantity to expression in the requested unit (scale-back).
//   qty = Quantity(0.05, {m=1}); convert(qty, "cm") → Quantity(5, {m=1})
//                                                     same dim, value scaled.
// Fails (Unimplemented) if Quantity dimensions don't match unit dimensions.
[[nodiscard]] Result<ExprPtr> convert_quantity(
    ExprPtr quantity, const std::string& unit_name, symbolic::CASContext& ctx);

}  // namespace cas::units
