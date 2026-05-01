#pragma once

#include "cas/ast.hpp"

#include <cstdint>
#include <vector>

namespace cas::symbolic {

enum class RuleId : std::uint16_t {
    Unknown = 0,
    RewriteProviderApplied,
    AssumptionApplied,
    SubstituteSymbol,
    SimplifyAddZero,
    SimplifyMultiplyByOne,
    SimplifyMultiplyByZero,
    SimplifyDivideByOne,
    SimplifySelfDivisionAssumedNonzero,
    SimplifyPowerZero,
    SimplifyPowerOne,
    SimplifyZeroPowerPositive,
    SimplifyNegativePower,
    SimplifyMergePowers,
    SimplifyFlattenNestedPowers,
    SimplifyCollectLikeTerms,
    SimplifyDistributeProductOverSum,
    SimplifySinZero,
    SimplifySinPi,
    SimplifyCosZero,
    SimplifyCosPi,
    SimplifyExpZero,
    SimplifyExpOne,
    SimplifyExpSum,
    SimplifyLnOne,
    SimplifyLnE,
    SimplifyLnExp,
    SimplifyExpLnPositive,
    SimplifyLnPowerPositiveBase,
    SimplifyLnPositiveProduct,
    SimplifyLnSqrtPositive,
    SimplifySqrtPositiveProduct,
    SimplifySqrtSquare,
    SimplifyTrigPythagorean,
    PolynomialGcdSubresultant,
    PolynomialGcdPrimitiveFallback,
    PolynomialGcdPrimitiveFallbackPsi,
    PolynomialGcdPrimitiveFallbackBeta,
    PolynomialGcdSymbolicEuclidean,
};

struct TraceStep {
    static constexpr std::uint32_t API_VERSION = 1;

    RuleId rule_id{RuleId::Unknown};
    std::uint8_t depth{0};
    ExprPtr target_before{};
    ExprPtr target_after{};
    ExprPtr root_after{};
};

using ComputationTrace = std::vector<TraceStep>;

}  // namespace cas::symbolic
