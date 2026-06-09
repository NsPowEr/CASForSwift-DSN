#pragma once
// maxima_parser.hpp — translate Maxima text output into our CAS input syntax,
// then parse it through our own Lexer+Parser.
//
// Documented differences between Maxima output syntax and our parser:
//
// | Maxima output | Our parser equivalent | Notes                            |
// |---------------|-----------------------|----------------------------------|
// | %pi           | pi                    | Maxima uses % prefix for constants|
// | %e            | e                     | Euler's number                   |
// | %i            | i                     | Imaginary unit                   |
// | %gamma        | EulerGamma            | Euler-Mascheroni (if supported)  |
// | log(x)        | log(x)                | Natural log — same               |
// | sqrt(x)       | sqrt(x)               | Same                             |
// | x^2           | x^2                   | Same                             |
// | x*y           | x*y                   | Same                             |
// | -1/2          | -1/2                  | Same                             |
// | (%o1) expr    | expr                  | Strip output label               |
// | display2d:f   | (suppressed)          | Batch mode removes 2D layout     |
// | \n continued  | single line           | Maxima may break long lines      |
// | inf           | inf                   | Same positive infinity           |
// | minf          | -inf                  | Maxima uses minf for -∞          |
// | false         | (SKIP)                | Maxima returns false on failure  |
// | done          | (SKIP)                | Assignment side-effect result    |

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <sstream>
#include <string>
#include <regex>
#include <vector>

namespace cas::golden {

// Normalise a raw Maxima output string to our parser syntax.
// Returns empty string if the output should be SKIPped.
inline std::string normalize_maxima_output(const std::string& raw) {
    std::string s = raw;

    // Remove output labels like (%o1) or (%o12)
    {
        std::regex label_re(R"(\(%o\d+\)\s*)");
        s = std::regex_replace(s, label_re, "");
    }

    // Strip leading/trailing whitespace and newlines
    while (!s.empty() && (s.front() == ' ' || s.front() == '\n' || s.front() == '\r' || s.front() == '\t'))
        s.erase(s.begin());
    while (!s.empty() && (s.back() == ' ' || s.back() == '\n' || s.back() == '\r' || s.back() == '\t' || s.back() == '$' || s.back() == ';'))
        s.pop_back();

    // Skip sentinel values
    if (s.empty() || s == "false" || s == "done" || s == "TIMEOUT" || s.rfind("ERROR:", 0) == 0)
        return "";

    // %pi -> pi
    {
        std::regex pi_re(R"(%pi)");
        s = std::regex_replace(s, pi_re, "pi");
    }
    // %e -> e  (but only standalone, not in words like "else")
    {
        std::regex e_re(R"(%e\b)");
        s = std::regex_replace(s, e_re, "e");
    }
    // %i -> i
    {
        std::regex i_re(R"(%i)");
        s = std::regex_replace(s, i_re, "i");
    }
    // %gamma -> EulerGamma (if we ever support it)
    {
        std::regex g_re(R"(%gamma)");
        s = std::regex_replace(s, g_re, "EulerGamma");
    }
    // minf -> -inf
    {
        std::regex minf_re(R"(\bminf\b)");
        s = std::regex_replace(s, minf_re, "-inf");
    }

    // Maxima outputs e^x (with bare 'e') instead of exp(x).
    // Convert: e^(...) -> exp(...) when 'e' is a standalone base.
    // Pattern: matches 'e^' followed by a simple token or parenthesised expr.
    // We do a simple iterative replacement: e^(expr) -> exp(expr)
    // and e^x (bare symbol/number) -> exp(x).
    {
        // First handle e^(grouped): e^(expr) -> exp(expr)
        std::regex e_paren_re(R"(\be\^(\([^)]*\)))");
        s = std::regex_replace(s, e_paren_re, "exp$1");
        // Then handle e^token (no parens): e^token -> exp(token)
        // token = word chars, digits, or ^-chain like x^2
        std::regex e_token_re(R"(\be\^([A-Za-z0-9_]+(?:\^[A-Za-z0-9_]+)*))");
        s = std::regex_replace(s, e_token_re, "exp($1)");
    }

    // Maxima sec(x)^2 — keep as-is; our parser handles sec.
    // Maxima log(abs(x)) — normalise to log(abs(x)); our parser handles.
    // Maxima abs(x) — same.

    // Maxima series() output (unevaluated) — SKIP
    if (s.find("series(") != std::string::npos ||
        s.find("taylor(") != std::string::npos)
        return "";

    // Maxima sometimes outputs "unknown(expr, x, 0, n)" for unrecognised series
    if (s.find("unknown(") != std::string::npos)
        return "";

    // Maxima sometimes outputs "gamma_incomplete(...)" — SKIP (no mapping yet)
    if (s.find("gamma_incomplete") != std::string::npos)
        return "";

    // Maxima sometimes outputs "bessel_j(...)" — SKIP (different arg order)
    if (s.find("bessel_j(") != std::string::npos)
        return "";

    // Maxima sometimes outputs "erf(%i*x)" — SKIP (complex erf)
    if (s.find("erf(i*") != std::string::npos || s.find("erf(-i") != std::string::npos)
        return "";

    // Maxima sometimes outputs "sec(x)^2" — our parser accepts sec, keep as-is.

    // Maxima sometimes outputs "asin(x)" — same as our atan, keep.

    // Maxima sometimes outputs "asinh(x)", "acosh(x)" — keep.

    // Maxima sometimes outputs "log(..." — keep (natural log).

    // Maxima sometimes outputs expressions with "where { v1 = ... }" — strip
    {
        auto where_pos = s.find(" where ");
        if (where_pos != std::string::npos)
            s = s.substr(0, where_pos);
    }

    // Maxima sometimes outputs matrix(row,row,...) — map back to [[...]]
    // We leave this for now; matrix inputs are already complex to parse back.
    // SKIP matrix outputs from Maxima (they require special round-trip).
    if (s.rfind("matrix(", 0) == 0)
        return "";

    // Maxima list output like [x=2, x=-2] from solve — skip, use algebraic equiv
    if (!s.empty() && s.front() == '[')
        return "";

    // Maxima integral output may contain 'integrate(...)' (unsolved) — skip
    if (s.find("integrate(") != std::string::npos)
        return "";

    return s;
}

// Parse a normalised Maxima output string into an ExprPtr using our CAS parser.
// Returns error if parse fails.
inline Result<ExprPtr> parse_maxima_expr(const std::string& normalised, AstArena& arena) {
    if (normalised.empty())
        return CASError{CASErrorKind::Unimplemented, "empty/skip Maxima output", std::nullopt};

    auto tok_result = Lexer(normalised).tokenize();
    if (!tok_result.is_ok())
        return CASError{CASErrorKind::Unimplemented, "Maxima output lex error: " + normalised, std::nullopt};

    Parser p(tok_result.value(), arena);
    return p.parse();
}

// Extract the meaningful last expression line from a Maxima .out file content.
// Maxima with display2d:false emits lines like:
//   (%o1) expr
// We want the last non-blank, non-label-only line.
inline std::string extract_maxima_result_line(const std::string& file_content) {
    // Two-pass scan:
    //  1) collect candidate lines (skip blank/comment-style headers and
    //     the input echo `display2d:false`);
    //  2) join continuation lines belonging to the final expression.
    //
    // Maxima with display2d:false can wrap long list outputs across
    // multiple lines without an explicit continuation marker, e.g.
    //
    //   solve([x^6-1],[x])
    //   [x = (sqrt(3)*%i+1)/2,x = (sqrt(3)*%i-1)/2,x = -1,
    //    x = -((sqrt(3)*%i+1)/2),x = -((sqrt(3)*%i-1)/2),x = 1]
    //
    // We need both lines as a single string.
    std::vector<std::string> lines;
    {
        std::istringstream ss(file_content);
        std::string line;
        while (std::getline(ss, line)) {
            while (!line.empty() && (line.front() == ' ' || line.front() == '\t'))
                line.erase(line.begin());
            while (!line.empty() && (line.back() == ' ' || line.back() == '\r'))
                line.pop_back();
            if (line.empty()) continue;
            // Skip the display2d:false setting echo.
            if (line == "display2d:false" || line == "display2d:false;" ||
                line == "display2d:false$")
                continue;
            // Skip lines that are pure prompt labels like "(%i1)" without
            // body, but keep "(%o1) expr".
            lines.push_back(line);
        }
    }
    if (lines.empty()) return "";

    // Find the last line that starts an expression: it either starts with
    // `(` (output label) or any other non-paren char, but we must NOT
    // count input-echo lines (e.g. `solve([x^6-1],[x])` is the input
    // echo, not the answer). Heuristic: when a multi-line list `[…]`
    // straddles lines, the answer starts with `[` somewhere in the tail.
    //
    // Strategy:
    //   - Walk from the end backwards.
    //   - If the last line ends with `]` and contains the matching `[`
    //     anchor on an earlier line, join from that anchor onward.
    //   - Else: return the last single line (existing behaviour).
    auto strip_label = [](std::string s) {
        // Strip leading "(%oN)" prefix.
        if (s.rfind("(%o", 0) == 0) {
            auto rp = s.find(')');
            if (rp != std::string::npos) {
                s = s.substr(rp + 1);
                while (!s.empty() && (s.front() == ' ' || s.front() == '\t'))
                    s.erase(s.begin());
            }
        }
        return s;
    };

    const std::string& last = lines.back();
    if (!last.empty() && last.back() == ']') {
        // Walk backwards to find the line starting the list.
        std::ptrdiff_t i = static_cast<std::ptrdiff_t>(lines.size()) - 1;
        for (; i >= 0; --i) {
            std::string stripped = strip_label(lines[static_cast<std::size_t>(i)]);
            if (!stripped.empty() && stripped.front() == '[') break;
        }
        if (i >= 0) {
            std::string joined;
            for (std::size_t k = static_cast<std::size_t>(i); k < lines.size(); ++k) {
                joined += strip_label(lines[k]);
            }
            return joined;
        }
    }

    return strip_label(lines.back());
}

} // namespace cas::golden
