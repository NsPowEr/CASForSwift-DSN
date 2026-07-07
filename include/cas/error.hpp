#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include "cas/unimplemented_info.hpp"

namespace cas {

enum class CASErrorKind : std::uint8_t {
    Unsupported,
    Unimplemented,
    Undefined,
    Overflow,
    Timeout,
    ParseError,
    AssumptionConflict,
    InvalidArgument,
    InternalError,
    DivisionByZero,

};

// ---------------------------------------------------------------------------
// F0.8 — Structured diagnostic payload for Unimplemented errors
// ---------------------------------------------------------------------------

/// Reason codes for Unimplemented errors. Use the constants defined in
/// namespace cas::error::reason_codes rather than raw strings.
namespace error::reason_codes {
    static constexpr const char* RISCH_LOG_EXTENSION_GENERAL   = "RISCH_LOG_EXT_GENERAL";
    static constexpr const char* RISCH_EXPONENTIAL_DE          = "RISCH_EXP_DE";
    static constexpr const char* RISCH_NO_POLYNOMIAL_SOLUTION  = "RISCH_NO_POLY_SOL";
    static constexpr const char* RISCH_SINGULAR_SYSTEM         = "RISCH_SINGULAR_SYS";
    static constexpr const char* RISCH_NO_MATCH                = "RISCH_NO_MATCH";
    static constexpr const char* INTEGRATE_DECIMAL_INPUT       = "INTEGRATE_DECIMAL_INPUT";
    static constexpr const char* INTEGRATE_FACTORIAL           = "INTEGRATE_FACTORIAL";
    static constexpr const char* INTEGRATE_DEPTH_EXCEEDED      = "INTEGRATE_DEPTH_EXCEEDED";
    static constexpr const char* INTEGRATE_UNKNOWN_EXPR        = "INTEGRATE_UNKNOWN_EXPR";
    static constexpr const char* INTEGRATE_MODULO              = "INTEGRATE_MODULO";
    static constexpr const char* INTEGRATE_COMPARISON          = "INTEGRATE_COMPARISON";
    static constexpr const char* INTEGRATE_NO_STRATEGY         = "INTEGRATE_NO_STRATEGY";
    static constexpr const char* INTEGRATE_TRIG_GENERAL        = "INTEGRATE_TRIG_GENERAL";
    static constexpr const char* DIFF_UNKNOWN_FUNCTION         = "DIFF_UNKNOWN_FUNCTION";
    static constexpr const char* DIFF_UNARY_OP                 = "DIFF_UNARY_OP";
    static constexpr const char* ALGEBRA_MULTIVAR_NON_INTEGER  = "ALGEBRA_MULTIVAR_NON_INT";
    static constexpr const char* ALGEBRA_MULTIVAR_REMAINING_VARS = "ALGEBRA_MULTIVAR_REMAIN_VARS";
    static constexpr const char* ALGEBRA_EXPAND_DECIMAL        = "ALGEBRA_EXPAND_DECIMAL";
    static constexpr const char* POLY_GCD_MULTIVAR             = "POLY_GCD_MULTIVAR";
    static constexpr const char* POLY_FACTOR_EXTENSION         = "POLY_FACTOR_EXTENSION";
    static constexpr const char* LAPLACE_UNKNOWN_FORM          = "LAPLACE_UNKNOWN_FORM";
    static constexpr const char* ODE_UNSUPPORTED_TYPE          = "ODE_UNSUPPORTED_TYPE";
    static constexpr const char* SUMMATION_GENERAL             = "SUMMATION_GENERAL";
    static constexpr const char* SERIES_GENERAL                = "SERIES_GENERAL";
    static constexpr const char* LINALG_GENERAL                = "LINALG_GENERAL";
    static constexpr const char* LINALG_ZERO_PIVOT             = "LINALG_ZERO_PIVOT";
    static constexpr const char* LINALG_SINGULAR_MATRIX        = "LINALG_SINGULAR_MATRIX";
    static constexpr const char* LINALG_LINEAR_DEPENDENT       = "LINALG_LINEAR_DEPENDENT";
    static constexpr const char* LINALG_SMITH_NON_INTEGER      = "LINALG_SMITH_NON_INTEGER";
    static constexpr const char* LINALG_RREF_UNDECIDABLE_PIVOT = "LINALG_RREF_UNDECIDABLE_PIVOT";
    static constexpr const char* NUMERIC_DECIMAL_OP            = "NUMERIC_DECIMAL_OP";
    static constexpr const char* NUMERIC_UNSUPPORTED_CONSTANT  = "NUMERIC_UNSUPPORTED_CONSTANT";
    static constexpr const char* NUMERIC_UNSUPPORTED_UNARY_OP  = "NUMERIC_UNSUPPORTED_UNARY_OP";
    static constexpr const char* NUMERIC_UNSUPPORTED_BINARY_OP = "NUMERIC_UNSUPPORTED_BINARY_OP";
    static constexpr const char* NUMERIC_UNSUPPORTED_FUNCTION  = "NUMERIC_UNSUPPORTED_FUNCTION";
    static constexpr const char* NUMERIC_UNSUPPORTED_NODE_TYPE = "NUMERIC_UNSUPPORTED_NODE_TYPE";
    static constexpr const char* SYMBOLIC_DEGREE_TOO_LARGE     = "SYMBOLIC_DEGREE_TOO_LARGE";
    static constexpr const char* SYMBOLIC_SERIES_INVERSION     = "SYMBOLIC_SERIES_INVERSION";
    static constexpr const char* SYMBOLIC_NORMAL_FORM_DEPTH    = "SYMBOLIC_NORMAL_FORM_DEPTH";
    static constexpr const char* SYMBOLIC_UNITS_UNKNOWN        = "SYMBOLIC_UNITS_UNKNOWN";
    static constexpr const char* SYMBOLIC_ZETA_OVERFLOW        = "SYMBOLIC_ZETA_OVERFLOW";
    static constexpr const char* RECURSION_DEPTH_EXCEEDED     = "RECURSION_DEPTH_EXCEEDED";
    static constexpr const char* CYCLE_DETECTED                = "CYCLE_DETECTED";
    static constexpr const char* GENERIC                       = "GENERIC";
}  // namespace error::reason_codes

// ---------------------------------------------------------------------------

struct CASError {
    static constexpr std::uint32_t API_VERSION = 2;

    CASErrorKind kind{CASErrorKind::InternalError};
    std::string message;
    std::optional<std::string> hint;

    /// Structured diagnostic payload (present only for Unimplemented errors
    /// created via make_unimplemented()).
    /// Default-initialised to nullopt so existing designated-init aggregate
    /// constructions {.kind=..., .message=..., .hint=...} remain valid.
    std::optional<UnimplementedInfo> payload{std::nullopt};

    /// Produce a multi-line user-facing diagnostic string.
    /// Format:
    ///   [Unimplemented] module=<m> function=<f>
    ///     Input shape: <i>
    ///     Reason: <r>
    ///     Suggestion: <s>
    ///     Ticket: <t>
    /// Falls back to the plain message when no payload is present.
    [[nodiscard]] std::string format_user_message() const;
};

// ---------------------------------------------------------------------------
// Helper: build a structured Unimplemented CASError
// ---------------------------------------------------------------------------

[[nodiscard]] inline CASError make_unimplemented_error(
    UnimplementedInfo info,
    std::string brief_message = {})
{
    std::string msg = brief_message.empty()
        ? (info.reason + ": " + info.input_shape)
        : std::move(brief_message);
    return CASError{
        .kind    = CASErrorKind::Unimplemented,
        .message = std::move(msg),
        .hint    = std::nullopt,
        .payload = std::move(info),
    };
}

[[nodiscard]] inline CASError make_unimplemented_error(
    std::string module,
    std::string function,
    std::string input_shape,
    std::string reason,
    std::string suggestion,
    std::string ticket = {},
    std::string brief_message = {})
{
    return make_unimplemented_error(
        UnimplementedInfo{
            .module      = std::move(module),
            .function    = std::move(function),
            .input_shape = std::move(input_shape),
            .reason      = std::move(reason),
            .suggestion  = std::move(suggestion),
            .ticket      = std::move(ticket),
        },
        std::move(brief_message));
}

// ---------------------------------------------------------------------------
// Typed helper that returns a Result<T> error directly.
// Usage:  return make_unimplemented<ExprPtr>(...);
// ---------------------------------------------------------------------------

}  // namespace cas

// For the typed Result<T> helper, include "cas/error_helpers.hpp"
// (which itself includes "cas/result.hpp"). Do NOT include it here to
// avoid a circular dependency: result.hpp → error.hpp → result.hpp.
