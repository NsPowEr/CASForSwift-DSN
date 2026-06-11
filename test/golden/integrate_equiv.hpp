#pragma once
// F7.5 follow-up — antiderivative equivalence check.
//
// Two antiderivatives of the same integrand differ at most by a
// constant. The corpus integrate area frequently shows this:
//
//   integrate(1/x, x):
//     CAS    = ln(abs(x))
//     Maxima = log(x)
//
// Both have derivative 1/x, so they are equal as antiderivatives
// (up to a constant and on the principal real branch). Direct
// `mathematically_equal` rejects them because `ln(abs(x))` and
// `log(x)` are not algebraically identical (one has the abs guard
// the other lacks). This helper differentiates both sides w.r.t.
// the integration variable, simplifies, and re-runs
// `mathematically_equal`. If derivatives match, the antiderivative
// pair is accepted as a corpus PASS.

#include "corpus_runner.hpp"

#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

#include <string>

namespace cas::golden {

[[nodiscard]] inline Result<bool> antiderivative_equivalent(
    const std::string& input_str,
    cas::ExprPtr cas_value,
    cas::ExprPtr maxima_value,
    cas::symbolic::CASContext& ctx)
{
    auto cmd = parse_command(input_str);
    if (cmd.arg_strs.size() < 2U) {
        return ok(false);
    }
    // Integration variable is the second positional argument.
    std::string var_name = cmd.arg_strs[1];
    // Trim surrounding whitespace defensively.
    while (!var_name.empty() && (var_name.front() == ' ' || var_name.front() == '\t'))
        var_name.erase(var_name.begin());
    while (!var_name.empty() && (var_name.back() == ' ' || var_name.back() == '\t'))
        var_name.pop_back();
    if (var_name.empty()) return ok(false);

    cas::Symbol var{var_name};
    auto d_cas = cas::calculus::diff(cas_value, var, 1U, ctx);
    if (!d_cas.is_ok()) return ok(false);
    auto d_max = cas::calculus::diff(maxima_value, var, 1U, ctx);
    if (!d_max.is_ok()) return ok(false);
    auto cs = ctx.simplify(d_cas.value());
    if (!cs.is_ok()) return ok(false);
    auto ms = ctx.simplify(d_max.value());
    if (!ms.is_ok()) return ok(false);
    auto eq = cas::symbolic::mathematically_equal(cs.value(), ms.value(), ctx);
    if (!eq.is_ok()) return ok(false);
    return ok(eq.value());
}

}  // namespace cas::golden
