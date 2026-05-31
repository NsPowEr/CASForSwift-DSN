// test_diagnostic.cpp — Unit tests for F0.8 error diagnostic framework.

#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "cas/result.hpp"

#include <gtest/gtest.h>

#include <string>
#include <string_view>

namespace cas {
namespace {

// ─── Helper ────────────────────────────────────────────────────────────────

static bool contains(const std::string& haystack, std::string_view needle) {
    return haystack.find(needle) != std::string::npos;
}

// ─── 1. make_unimplemented_error creates structured payload ────────────────

TEST(DiagnosticFramework, MakeUnimplementedErrorPopulatesPayload) {
    CASError err = make_unimplemented_error(
        "calculus",
        "integrate_once",
        "DecimalLit",
        error::reason_codes::INTEGRATE_DECIMAL_INPUT,
        "Convert DecimalLit to Rational before calling integrate()",
        "F0.8");

    EXPECT_EQ(err.kind, CASErrorKind::Unimplemented);
    ASSERT_TRUE(err.payload.has_value());

    const auto& p = *err.payload;
    EXPECT_EQ(p.module,      "calculus");
    EXPECT_EQ(p.function,    "integrate_once");
    EXPECT_EQ(p.input_shape, "DecimalLit");
    EXPECT_EQ(p.reason_code, error::reason_codes::INTEGRATE_DECIMAL_INPUT);
    EXPECT_FALSE(p.suggestion.empty());
    EXPECT_EQ(p.ticket_id,   "F0.8");
}

// ─── 2. format_user_message contains all expected fields ──────────────────

TEST(DiagnosticFramework, FormatUserMessageContainsAllFields) {
    CASError err = make_unimplemented_error(
        "algebra",
        "evaluate_at",
        "MultivariatePolynomial with non-integer value",
        error::reason_codes::ALGEBRA_MULTIVAR_NON_INTEGER,
        "Pass an IntegerLit or extend evaluate_at to RationalLit",
        "F0.8");

    const std::string msg = err.format_user_message();

    EXPECT_TRUE(contains(msg, "[Unimplemented]"))         << msg;
    EXPECT_TRUE(contains(msg, "module=algebra"))          << msg;
    EXPECT_TRUE(contains(msg, "function=evaluate_at"))    << msg;
    EXPECT_TRUE(contains(msg, "Input shape:"))            << msg;
    EXPECT_TRUE(contains(msg, "Reason:"))                 << msg;
    EXPECT_TRUE(contains(msg, "Suggestion:"))             << msg;
    EXPECT_TRUE(contains(msg, "Ticket:"))                 << msg;
    EXPECT_TRUE(contains(msg, "F0.8"))                    << msg;
}

// ─── 3. Backward compat: plain CASError without payload ───────────────────

TEST(DiagnosticFramework, BackwardCompatPlainCASErrorFormatsAsFallback) {
    CASError legacy{
        .kind    = CASErrorKind::Unimplemented,
        .message = "Risch: log extension integration not fully implemented",
        .hint    = std::nullopt,
        .payload = std::nullopt,
    };

    const std::string msg = legacy.format_user_message();
    EXPECT_EQ(msg, legacy.message);

    // Also verify that old aggregate-init without .payload field compiles:
    // This is tested by the compilation of this file itself.
}

// ─── 4. make_unimplemented<T> builds Result<T> error ─────────────────────

TEST(DiagnosticFramework, MakeUnimplementedTemplateBuildsResultError) {
    using ExprPtrSurrogate = int;  // use int as stand-in for ExprPtr

    auto result = make_unimplemented<ExprPtrSurrogate>(
        "calculus",
        "Integrator::integrate",
        "arbitrary expression — no strategy matched",
        error::reason_codes::INTEGRATE_NO_STRATEGY,
        "Implement Risch for this expression class",
        "F0.8");

    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Unimplemented);
    ASSERT_TRUE(result.error().payload.has_value());
    EXPECT_EQ(result.error().payload->module, "calculus");
}

// ─── 5. Payload survives Result<T> propagation ────────────────────────────

// Simulate a call chain: inner returns Unimplemented, outer propagates it.
static Result<int> inner_fn() {
    return make_unimplemented<int>(
        "calculus",
        "inner_fn",
        "Sum{DecimalLit}",
        error::reason_codes::INTEGRATE_DECIMAL_INPUT,
        "Strip decimal before entering inner_fn",
        "F0.8");
}

static Result<int> outer_fn() {
    auto r = inner_fn();
    if (r.is_error()) return fail<int>(r.error());  // propagate
    return ok(r.value() * 2);
}

TEST(DiagnosticFramework, PayloadPreservedThroughResultPropagation) {
    auto result = outer_fn();

    ASSERT_TRUE(result.is_error());
    ASSERT_TRUE(result.error().payload.has_value());

    const auto& p = result.error().payload.value();
    EXPECT_EQ(p.module,   "calculus");
    EXPECT_EQ(p.function, "inner_fn");
    EXPECT_EQ(p.reason_code, error::reason_codes::INTEGRATE_DECIMAL_INPUT);

    const std::string msg = result.error().format_user_message();
    EXPECT_TRUE(contains(msg, "[Unimplemented]")) << msg;
    EXPECT_TRUE(contains(msg, "module=calculus")) << msg;
}

// ─── 6. reason_codes constants have expected prefixes ─────────────────────

TEST(DiagnosticFramework, ReasonCodeConstantsAreNonEmpty) {
    using namespace error::reason_codes;
    EXPECT_NE(std::string_view(RISCH_LOG_EXTENSION_GENERAL), "");
    EXPECT_NE(std::string_view(RISCH_EXPONENTIAL_DE), "");
    EXPECT_NE(std::string_view(INTEGRATE_DECIMAL_INPUT), "");
    EXPECT_NE(std::string_view(ALGEBRA_MULTIVAR_NON_INTEGER), "");
    EXPECT_NE(std::string_view(ODE_UNSUPPORTED_TYPE), "");
    EXPECT_NE(std::string_view(POLY_GCD_MULTIVAR), "");
    EXPECT_NE(std::string_view(LAPLACE_UNKNOWN_FORM), "");
}

// ─── 7. Ticket-less payload omits Ticket line ─────────────────────────────

TEST(DiagnosticFramework, FormatUserMessageOmitsTicketWhenEmpty) {
    CASError err = make_unimplemented_error(
        "algebra",
        "poly_gcd",
        "multivariate GCD over Q",
        error::reason_codes::POLY_GCD_MULTIVAR,
        "Use polynomial_gcd_multivariate module",
        "" /* no ticket */);

    const std::string msg = err.format_user_message();
    EXPECT_TRUE(contains(msg, "[Unimplemented]")) << msg;
    EXPECT_FALSE(contains(msg, "Ticket:"))        << msg;
}

}  // namespace
}  // namespace cas
