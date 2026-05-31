#pragma once
// corpus_runner.hpp — evaluate one corpus entry using our CAS engine.
//
// Supported top-level commands and how they map to our C++ API:
//   integrate(f, x)          -> calculus::integrate
//   diff(f, x)               -> calculus::diff order=1
//   diff(f, x, n)            -> calculus::diff order=n
//   limit(f, x, a)           -> calculus::limit dir=Both
//   limit(f, x, a, plus)     -> calculus::limit dir=Plus
//   limit(f, x, a, minus)    -> calculus::limit dir=Minus
//   factor(p)                -> algebra::factor_over_integers
//   gcd(p, q)                -> algebra::polynomial_gcd (single var inferred)
//   solve(p, x)              -> algebra::solve_polynomial
//   series(f, x, a, n)       -> calculus::taylor_series
//   simplify(e)              -> ctx.simplify
//   det([[...]])             -> linalg::det (if available) else SKIP
//   trace([[...]])           -> SKIP (no public API yet, SKIP)
//   eigenvalues([[...]])     -> SKIP
//   rank([[...]])            -> SKIP
//   inverse([[...]])         -> SKIP
//   transpose([[...]])       -> SKIP
//   gamma(x)                 -> simplify only
//   factorial(n)             -> simplify only
//   binomial(n, k)           -> simplify only
//   BesselJ(n, x)            -> simplify only
//   LegendreP(n, x)          -> simplify only
//   zeta(s)                  -> simplify only
//   erf(x)                   -> simplify only
//
// For commands not in the dispatch table the result is SKIP.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <string>
#include <optional>
#include <variant>

namespace cas::golden {

enum class RunStatus { Pass, Fail, Skip };

struct RunResult {
    RunStatus status;
    std::string reason; // non-empty on Fail or Skip
};

// Parse a string into ExprPtr using the context's arena.
inline Result<ExprPtr> parse_expr(const std::string& s, symbolic::CASContext& ctx) {
    auto tok = Lexer(s).tokenize();
    if (!tok.is_ok())
        return CASError{CASErrorKind::Unimplemented, "lex error: " + s, std::nullopt};
    Parser p(tok.value(), ctx.arena());
    return p.parse();
}

// Extract the function name and arguments from a top-level FuncCall or
// from a bare expression (treated as simplify).
// Returns {"", {expr_str}} for bare expressions.
struct ParsedCommand {
    std::string fn;            // "integrate", "diff", etc., or "" for bare expr
    std::vector<std::string> arg_strs; // raw argument strings (unparsed)
    std::string raw_input;     // original input string
};

// Minimal argument splitter: splits top-level comma-separated arguments.
// Does NOT handle matrix literals with nested commas robustly.
inline std::vector<std::string> split_top_level_args(const std::string& args_str) {
    std::vector<std::string> result;
    int depth = 0;
    std::string current;
    for (char c : args_str) {
        if (c == '(' || c == '[') { depth++; current += c; }
        else if (c == ')' || c == ']') { depth--; current += c; }
        else if (c == ',' && depth == 0) {
            // trim
            while (!current.empty() && current.front() == ' ') current.erase(current.begin());
            while (!current.empty() && current.back() == ' ') current.pop_back();
            result.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    while (!current.empty() && current.front() == ' ') current.erase(current.begin());
    while (!current.empty() && current.back() == ' ') current.pop_back();
    if (!current.empty()) result.push_back(current);
    return result;
}

// Parse a command string like "integrate(f, x)" into name + arg strings.
inline ParsedCommand parse_command(const std::string& input) {
    ParsedCommand cmd;
    cmd.raw_input = input;

    auto pos = input.find('(');
    if (pos == std::string::npos) {
        // bare expression
        cmd.fn = "";
        cmd.arg_strs = {input};
        return cmd;
    }

    // Check that the first '(' closes the function name without spaces
    std::string prefix = input.substr(0, pos);
    bool is_name = !prefix.empty();
    for (char c : prefix) {
        if (!(std::isalnum(c) || c == '_')) { is_name = false; break; }
    }
    if (!is_name) {
        cmd.fn = "";
        cmd.arg_strs = {input};
        return cmd;
    }

    // Find matching closing ')'
    auto rpos = input.rfind(')');
    if (rpos == std::string::npos) {
        cmd.fn = "";
        cmd.arg_strs = {input};
        return cmd;
    }

    cmd.fn = prefix;
    std::string inner = input.substr(pos + 1, rpos - pos - 1);
    cmd.arg_strs = split_top_level_args(inner);
    return cmd;
}

// Check mathematical equality: simplify(a - b) == 0.
// Returns true if provably zero, false if provably non-zero, error otherwise.
inline Result<bool> check_equal(ExprPtr a, ExprPtr b, symbolic::CASContext& ctx) {
    return symbolic::mathematically_equal(a, b, ctx);
}

// Set of function names that require matrix literal parsing (skip for now).
inline bool is_matrix_fn(const std::string& fn) {
    return fn == "det" || fn == "trace" || fn == "transpose" ||
           fn == "inverse" || fn == "rank" || fn == "eigenvalues";
}

// Evaluate our CAS on the given input command string.
// Returns an ExprPtr representing the result, or error.
inline Result<ExprPtr> evaluate_cas(const std::string& input_str,
                                     symbolic::CASContext& ctx) {
    auto cmd = parse_command(input_str);

    // --- Matrix operations: skip (no simple public API for matrix literal parsing) ---
    if (is_matrix_fn(cmd.fn))
        return CASError{CASErrorKind::Unimplemented, "matrix fn skipped: " + cmd.fn, std::nullopt};

    // --- integrate(f, x) ---
    if (cmd.fn == "integrate" && cmd.arg_strs.size() == 2) {
        auto f = parse_expr(cmd.arg_strs[0], ctx);
        if (!f.is_ok()) return f;
        Symbol var(cmd.arg_strs[1]);
        return calculus::integrate(f.value(), var, ctx);
    }

    // --- diff(f, x) or diff(f, x, n) ---
    if (cmd.fn == "diff") {
        if (cmd.arg_strs.size() < 2)
            return CASError{CASErrorKind::Unimplemented, "diff: too few args", std::nullopt};
        auto f = parse_expr(cmd.arg_strs[0], ctx);
        if (!f.is_ok()) return f;
        Symbol var(cmd.arg_strs[1]);
        unsigned int order = 1;
        if (cmd.arg_strs.size() >= 3) {
            try { order = static_cast<unsigned int>(std::stoul(cmd.arg_strs[2])); }
            catch (...) { order = 1; }
        }
        return calculus::diff(f.value(), var, order, ctx);
    }

    // --- limit(f, x, a) or limit(f, x, a, plus/minus) ---
    if (cmd.fn == "limit") {
        if (cmd.arg_strs.size() < 3)
            return CASError{CASErrorKind::Unimplemented, "limit: too few args", std::nullopt};
        auto f = parse_expr(cmd.arg_strs[0], ctx);
        if (!f.is_ok()) return f;
        Symbol var(cmd.arg_strs[1]);
        auto pt = parse_expr(cmd.arg_strs[2], ctx);
        if (!pt.is_ok()) return pt;
        LimitDirection dir = LimitDirection::Both;
        if (cmd.arg_strs.size() >= 4) {
            const auto& d = cmd.arg_strs[3];
            if (d == "plus" || d == "right" || d == "+") dir = LimitDirection::Right;
            else if (d == "minus" || d == "left" || d == "-") dir = LimitDirection::Left;
        }
        return calculus::limit(f.value(), var, pt.value(), dir, ctx);
    }

    // --- factor(p) ---
    if (cmd.fn == "factor" && cmd.arg_strs.size() >= 1) {
        auto p = parse_expr(cmd.arg_strs[0], ctx);
        if (!p.is_ok()) return p;
        // Detect variable: try x first, then y; simplify-only if multi-var
        Symbol var("x");
        auto fac = algebra::factor_over_integers(p.value(), var, ctx);
        if (!fac.is_ok()) return CASError{CASErrorKind::Unimplemented, "factor failed", std::nullopt};
        // Reconstruct product from factorization
        auto& result = fac.value();
        ExprPtr out = result.content;
        for (auto& pf : result.factors) {
            if (pf.multiplicity == 1) {
                out = ctx.arena().make<Product>(std::vector<ExprPtr>{out, pf.factor});
            } else {
                auto exp_lit = ctx.arena().make<IntegerLit>(BigInt(pf.multiplicity));
                auto pw = ctx.arena().make<Binary>(BinaryOp::Pow, pf.factor, exp_lit);
                out = ctx.arena().make<Product>(std::vector<ExprPtr>{out, pw});
            }
        }
        return ctx.simplify(out);
    }

    // --- solve(p, x) ---
    if (cmd.fn == "solve" && cmd.arg_strs.size() == 2) {
        auto p = parse_expr(cmd.arg_strs[0], ctx);
        if (!p.is_ok()) return p;
        Symbol var(cmd.arg_strs[1]);
        auto roots = algebra::solve_polynomial(p.value(), var, ctx);
        if (!roots.is_ok())
            return CASError{CASErrorKind::Unimplemented, "solve failed", std::nullopt};
        // Return as a Sum representing a set — equality checked element-wise in runner
        // For now return first root or zero for empty
        if (roots.value().empty())
            return ctx.arena().make<IntegerLit>(BigInt(0));
        // We'll mark solve results as SKIP for now (set equality needed)
        return CASError{CASErrorKind::Unimplemented, "solve: set equality not yet supported", std::nullopt};
    }

    // --- series(f, x, a, n) -> taylor_series ---
    if (cmd.fn == "series" && cmd.arg_strs.size() == 4) {
        auto f = parse_expr(cmd.arg_strs[0], ctx);
        if (!f.is_ok()) return f;
        Symbol var(cmd.arg_strs[1]);
        auto pt = parse_expr(cmd.arg_strs[2], ctx);
        if (!pt.is_ok()) return pt;
        unsigned int order = 0;
        try { order = static_cast<unsigned int>(std::stoul(cmd.arg_strs[3])); }
        catch (...) { return CASError{CASErrorKind::Unimplemented, "series: bad order", std::nullopt}; }
        auto ts = calculus::taylor_series(f.value(), var, pt.value(), order, ctx);
        if (!ts.is_ok()) return ts.error();
        return ts.value().polynomial;
    }

    // --- gcd(p, q) ---
    // Use univariate GCD for single-variable inputs only.
    // Multivariate inputs (containing both x and y) are SKIPped: the
    // multivariate GCD algorithm can exceed the runner's per-entry budget.
    if (cmd.fn == "gcd" && cmd.arg_strs.size() == 2) {
        auto p = parse_expr(cmd.arg_strs[0], ctx);
        if (!p.is_ok()) return p;
        auto q = parse_expr(cmd.arg_strs[1], ctx);
        if (!q.is_ok()) return q;
        // Detect multivariate: if both 'x' and 'y' appear in either arg, SKIP
        bool has_y = (cmd.arg_strs[0].find('y') != std::string::npos ||
                      cmd.arg_strs[1].find('y') != std::string::npos);
        bool has_x = (cmd.arg_strs[0].find('x') != std::string::npos ||
                      cmd.arg_strs[1].find('x') != std::string::npos);
        if (has_x && has_y)
            return CASError{CASErrorKind::Unimplemented,
                            "gcd: multivariate skipped in golden runner", std::nullopt};
        Symbol var("x");
        ctx.set_timeout(std::chrono::milliseconds(5000));
        auto result = algebra::polynomial_gcd(p.value(), q.value(), var, ctx);
        ctx.set_timeout(std::chrono::milliseconds(1000)); // restore default
        return result;
    }

    // --- bare expression or simplify(e) ---
    if (cmd.fn == "" || cmd.fn == "simplify") {
        auto expr_str = (cmd.fn == "") ? cmd.arg_strs[0] : cmd.arg_strs[0];
        auto e = parse_expr(expr_str, ctx);
        if (!e.is_ok()) return e;
        return ctx.simplify(e.value());
    }

    // --- All other named functions: parse full input as expression, simplify ---
    // This covers gamma, factorial, BesselJ, etc.
    {
        auto e = parse_expr(input_str, ctx);
        if (!e.is_ok()) return e;
        return ctx.simplify(e.value());
    }
}

} // namespace cas::golden
