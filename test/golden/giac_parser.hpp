#pragma once
// giac_parser.hpp — parse Giac 2.0.0 (icas) raw output text into our AST.
//
// A35 — Giac's surface syntax for the results actually produced by the
// golden corpus areas is very close to ours (same function names: ln, atan,
// exp, sqrt, abs, gamma, Ei... — verified against build-golden/giac_out/*
// samples written by run_golden_giac.sh). The only observed deltas:
//   - solve() results: `list[r1,r2,...]` instead of a bare set — no `x = r`
//     equation wrapper (giac returns bare root values, unlike Maxima).
//   - matrix-construction results carry a `matrix[[..],[..]]` tag;
//     det/trace/transpose/elementwise ops return plain `[[..],[..]]`.
//   - series() results carry giac's own truncation marker `order_size(x)`,
//     which has no equivalent node in our AST — left deliberately
//     unparseable (SKIP is sound: same choice maxima_parser.hpp makes for
//     ITS OWN truncation markers, never silently coerced).
//
// GPL-3.0-or-later compliance (CLAUDE.md §6): this is a plain TEXT parser
// over giac's fork/exec output (.giac.out files written by
// run_golden_giac.sh); no giac source is consulted, modified, or embedded.

#include "solve_set_equal.hpp" // split_maxima_list: dialect-agnostic bracket splitting

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/result.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace cas::golden {

// ---------------------------------------------------------------------------
// A35 — read a <idx>.giac.out file written by run_golden_giac.sh: first
// line is the verdict tag (ANSWERED|UNEVALUATED|TIMEOUT|ERROR), the rest
// (when ANSWERED) is the raw result text (possibly multi-line — join with
// a single space, giac's tail -5 output is short).
// ---------------------------------------------------------------------------
struct GiacFileContent {
    std::string tag;    // "" if the file does not exist
    std::string result; // only meaningful when tag == "ANSWERED"
};

inline GiacFileContent read_giac_file(const std::string& giac_dir, int idx) {
    if (giac_dir.empty()) return {};
    std::ifstream f(giac_dir + "/" + std::to_string(idx) + ".giac.out");
    if (!f.is_open()) return {};
    GiacFileContent out;
    std::getline(f, out.tag);
    std::string line;
    bool first = true;
    while (std::getline(f, line)) {
        if (!first) out.result += " ";
        out.result += line;
        first = false;
    }
    return out;
}

// Strip a `list[...]` / `matrix[...]` container tag, leaving the bare
// bracketed literal our own parser (and split_maxima_list) already handle.
inline std::string strip_giac_container_tag(const std::string& s) {
    for (const char* tag : {"list", "matrix"}) {
        const std::size_t len = std::string(tag).size();
        if (s.size() > len && s.compare(0, len, tag) == 0 && s[len] == '[')
            return s.substr(len);
    }
    return s;
}

// Parse a scalar giac result into an ExprPtr via our own lexer/parser.
inline Result<ExprPtr> parse_giac_expr(const std::string& raw, AstArena& arena) {
    std::string s = strip_giac_container_tag(raw);
    if (s.find("order_size(") != std::string::npos) {
        return CASError{CASErrorKind::Unimplemented,
                        "giac series truncation marker has no CAS equivalent: " + s,
                        std::nullopt};
    }
    if (s.empty())
        return CASError{CASErrorKind::Unimplemented, "empty giac output", std::nullopt};
    auto tok = Lexer(s).tokenize();
    if (!tok.is_ok())
        return CASError{CASErrorKind::Unimplemented, "giac output lex error: " + s, std::nullopt};
    Parser p(tok.value(), arena);
    return p.parse();
}

// Parse a giac `list[r1, r2, ...]` solve result into a vector of root
// ExprPtrs (bare values, no `x = r` wrapper). `list[]` is a valid empty set.
inline Result<std::vector<ExprPtr>> parse_giac_solve_list(
    const std::string& raw, AstArena& arena) {
    std::string s = strip_giac_container_tag(raw); // "list[...]" -> "[...]"
    auto elements = split_maxima_list(s);           // top-level comma split
    std::vector<ExprPtr> out;
    out.reserve(elements.size());
    for (const auto& el : elements) {
        auto parsed = parse_giac_expr(el, arena);
        if (!parsed.is_ok()) return parsed.error();
        out.push_back(parsed.value());
    }
    return ok(std::move(out));
}

}  // namespace cas::golden
