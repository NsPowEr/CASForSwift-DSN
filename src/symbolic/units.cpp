// CAS-L3-08 — SI unit conversion implementation.

#include "cas/units.hpp"

#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"

#include <unordered_map>

namespace cas::units {

namespace {

// Build SIDimensions helper. SI exponents fit in int16_t.
[[nodiscard]] constexpr SIDimensions D(int m_, int kg_, int s_, int A_, int K_,
                                       int mol_, int cd_) noexcept {
    return SIDimensions{
        static_cast<int16_t>(m_),
        static_cast<int16_t>(kg_),
        static_cast<int16_t>(s_),
        static_cast<int16_t>(A_),
        static_cast<int16_t>(K_),
        static_cast<int16_t>(mol_),
        static_cast<int16_t>(cd_)};
}

// Static registry. Scale stored as Rational(num, den).
struct RegistryEntry {
    const char* name;
    long long scale_num;
    long long scale_den;
    SIDimensions dim;
};

// Base SI + common derived + non-SI conversions. All scales exact rationals.
constexpr RegistryEntry kUnitRegistry[] = {
    // Length (m=1)
    {"m",   1,    1,    D(1,0,0,0,0,0,0)},
    {"cm",  1,    100,  D(1,0,0,0,0,0,0)},
    {"mm",  1,    1000, D(1,0,0,0,0,0,0)},
    {"km",  1000, 1,    D(1,0,0,0,0,0,0)},
    {"ft",  3048, 10000, D(1,0,0,0,0,0,0)},   // 0.3048 m
    {"in",  254,  10000, D(1,0,0,0,0,0,0)},   // 0.0254 m
    {"mi",  1609344, 1000, D(1,0,0,0,0,0,0)}, // 1609.344 m

    // Mass (kg=1)
    {"kg",  1,    1,    D(0,1,0,0,0,0,0)},
    {"g",   1,    1000, D(0,1,0,0,0,0,0)},
    {"mg",  1,    1000000, D(0,1,0,0,0,0,0)},
    {"lb",  453592, 1000000, D(0,1,0,0,0,0,0)}, // 0.453592 kg

    // Time (s=1)
    {"s",   1,    1,    D(0,0,1,0,0,0,0)},
    {"ms",  1,    1000, D(0,0,1,0,0,0,0)},
    {"min", 60,   1,    D(0,0,1,0,0,0,0)},
    {"h",   3600, 1,    D(0,0,1,0,0,0,0)},

    // Current
    {"A",   1,    1,    D(0,0,0,1,0,0,0)},
    {"mA",  1,    1000, D(0,0,0,1,0,0,0)},

    // Temperature
    {"K",   1,    1,    D(0,0,0,0,1,0,0)},

    // Amount, Luminosity
    {"mol", 1,    1,    D(0,0,0,0,0,1,0)},
    {"cd",  1,    1,    D(0,0,0,0,0,0,1)},

    // Derived (exact)
    {"Hz",  1,    1,    D(0,0,-1,0,0,0,0)},
    {"N",   1,    1,    D(1,1,-2,0,0,0,0)},
    {"J",   1,    1,    D(2,1,-2,0,0,0,0)},
    {"W",   1,    1,    D(2,1,-3,0,0,0,0)},
    {"Pa",  1,    1,    D(-1,1,-2,0,0,0,0)},
    {"C",   1,    1,    D(0,0,1,1,0,0,0)},   // Coulomb (electric charge)
    {"V",   1,    1,    D(2,1,-3,-1,0,0,0)}, // Volt
    {"Ohm", 1,    1,    D(2,1,-3,-2,0,0,0)}, // Ohm

    // Non-SI energy
    {"cal", 4184, 1000, D(2,1,-2,0,0,0,0)},  // 4.184 J
    {"eV",  16021766340, 100000000000000000ULL, D(2,1,-2,0,0,0,0)}, // ≈ 1.602e-19 J
};

}  // namespace

std::optional<UnitInfo> lookup_unit(const std::string& name) noexcept {
    for (const auto& e : kUnitRegistry) {
        if (name == e.name) {
            return UnitInfo{
                .scale_to_si = Rational(BigInt(e.scale_num), BigInt(e.scale_den)),
                .dimensions = e.dim,
            };
        }
    }
    return std::nullopt;
}

Result<ExprPtr> make_quantity_from_unit(ExprPtr value,
                                        const std::string& unit_name,
                                        symbolic::CASContext& ctx) {
    auto info = lookup_unit(unit_name);
    if (!info) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::Unimplemented,
            "Unknown unit: " + unit_name, std::nullopt});
    }
    AstArena& arena = ctx.arena();
    // Multiply value by scale (Rational → expression).
    ExprPtr scale_expr = (info->scale_to_si.denominator() == BigInt(1))
        ? static_cast<ExprPtr>(arena.make<IntegerLit>(info->scale_to_si.numerator()))
        : static_cast<ExprPtr>(arena.make<RationalLit>(
            info->scale_to_si.numerator(), info->scale_to_si.denominator()));
    ExprPtr scaled = arena.make<Product>(std::vector<ExprPtr>{value, scale_expr});
    auto simp = ctx.simplify(scaled);
    if (simp.is_ok()) scaled = simp.value();
    return ok(static_cast<ExprPtr>(arena.make<Quantity>(scaled, info->dimensions)));
}

Result<ExprPtr> convert_quantity(ExprPtr quantity,
                                 const std::string& unit_name,
                                 symbolic::CASContext& ctx) {
    const auto* qty = expr_cast<Quantity>(quantity);
    if (!qty) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::InvalidArgument,
            "convert_quantity: expression is not a Quantity", std::nullopt});
    }
    auto info = lookup_unit(unit_name);
    if (!info) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::Unimplemented,
            "Unknown target unit: " + unit_name, std::nullopt});
    }
    if (!(qty->dimensions == info->dimensions)) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::Undefined,
            "convert_quantity: dimensional mismatch", std::nullopt});
    }
    // value / scale_to_si → value in requested unit
    AstArena& arena = ctx.arena();
    // value_in_unit = qty.value · (den/num) where scale = num/den.
    // Equivalently: qty.value / Rational(num, den).
    ExprPtr inv_scale = (info->scale_to_si.numerator() == BigInt(1))
        ? static_cast<ExprPtr>(arena.make<IntegerLit>(info->scale_to_si.denominator()))
        : static_cast<ExprPtr>(arena.make<RationalLit>(
            info->scale_to_si.denominator(), info->scale_to_si.numerator()));
    ExprPtr scaled_back = arena.make<Product>(
        std::vector<ExprPtr>{qty->value, inv_scale});
    auto simp = ctx.simplify(scaled_back);
    if (simp.is_ok()) scaled_back = simp.value();
    return ok(static_cast<ExprPtr>(arena.make<Quantity>(scaled_back, qty->dimensions)));
}

}  // namespace cas::units
