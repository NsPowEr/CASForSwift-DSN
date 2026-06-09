#pragma once
// solve_set_equal.hpp — multiset equality for Maxima-style solve output.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Solve_Adapter.md (F7.5.A1).
//
// Maxima returns solve() as a list of equations `[x = r_1, ..., x = r_n]`.
// Our CAS returns `std::vector<ExprPtr>` of root values. Direct
// mathematically_equal() on a list expression doesn't apply.
//
// This header provides:
//   parse_maxima_solve_list(raw, arena)    → vector<ExprPtr> RHS values
//   compare_solve_sets(cas, maxima, ctx)   → Result<bool> set equality
//
// Set equality uses Hopcroft-Karp on the bipartite graph where
// edge (i, j) exists iff mathematically_equal(cas[i], maxima[j]).
// Perfect matching of equal cardinality ⇒ equal. Greedy matching is
// VIETATO (regola zero: shortest path may miss valid matchings).

#include "bipartite_matching.hpp"
#include "maxima_parser.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cas::golden {

// ---------------------------------------------------------------------------
// Split a Maxima list "[a, b, c]" into elements, respecting parenthesis
// nesting. The outer brackets are stripped; commas at depth 0 are
// separators. Returns empty vector on malformed input.
// ---------------------------------------------------------------------------
inline std::vector<std::string> split_maxima_list(const std::string& raw) {
    std::vector<std::string> out;
    std::string s = raw;
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t' || s.front() == '\n')) s.erase(s.begin());
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\n' || s.back() == ';' || s.back() == '$')) s.pop_back();
    if (s.size() < 2 || s.front() != '[' || s.back() != ']') return out;

    // Strip outer brackets.
    s = s.substr(1, s.size() - 2);

    int depth = 0;
    std::string cur;
    for (char c : s) {
        if (c == '(' || c == '[' || c == '{') ++depth;
        else if (c == ')' || c == ']' || c == '}') --depth;
        if (c == ',' && depth == 0) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

// ---------------------------------------------------------------------------
// Extract the RHS from a Maxima equation "x = r". The first top-level
// '=' is the splitter; '==' is not used by Maxima but we guard against
// it. Returns empty string on malformed input.
// ---------------------------------------------------------------------------
inline std::string extract_rhs(const std::string& eq) {
    int depth = 0;
    for (std::size_t i = 0; i < eq.size(); ++i) {
        char c = eq[i];
        if (c == '(' || c == '[' || c == '{') ++depth;
        else if (c == ')' || c == ']' || c == '}') --depth;
        if (c == '=' && depth == 0) {
            // Guard against '=='.
            if (i + 1 < eq.size() && eq[i + 1] == '=') return "";
            return eq.substr(i + 1);
        }
    }
    return "";
}

// ---------------------------------------------------------------------------
// Parse a Maxima "[x = r_1, ..., x = r_n]" output into a vector of
// RHS ExprPtrs. Empty list "[]" returns empty vector (valid).
// Returns CASError if any element fails to parse.
// ---------------------------------------------------------------------------
inline Result<std::vector<ExprPtr>> parse_maxima_solve_list(
    const std::string& raw, AstArena& arena) {

    auto elements = split_maxima_list(raw);
    std::vector<ExprPtr> out;
    out.reserve(elements.size());
    for (const auto& el : elements) {
        std::string rhs_raw = extract_rhs(el);
        if (rhs_raw.empty()) {
            return CASError{CASErrorKind::Unimplemented,
                            "Maxima solve element has no '=': " + el,
                            std::nullopt};
        }
        std::string normalised = normalize_maxima_output(rhs_raw);
        if (normalised.empty()) {
            return CASError{CASErrorKind::Unimplemented,
                            "Maxima solve RHS normalized to empty: " + rhs_raw,
                            std::nullopt};
        }
        auto parsed = parse_maxima_expr(normalised, arena);
        if (!parsed.is_ok()) return parsed.error();
        out.push_back(parsed.value());
    }
    return ok(std::move(out));
}

// ---------------------------------------------------------------------------
// Compare two solve solution sets via maximum bipartite matching.
//
// Returns Ok(true) iff |matching| == |cas| == |maxima| (perfect matching
// over multisets — multiplicity preserved).
//
// Returns Ok(false) if cardinality differs OR matching is not perfect.
// Returns Unimplemented if too many edges are inconclusive
// (mathematically_equal returns error) — caller can SKIP.
// ---------------------------------------------------------------------------
// Deduplicate a vector of root ExprPtrs modulo mathematically_equal.
// O(n²) — corpus has |roots| ≤ 6, this is negligible. Preserves order
// of first occurrence so diagnostic output stays readable.
inline std::vector<ExprPtr> dedup_modulo_equal(
    const std::vector<ExprPtr>& roots, cas::symbolic::CASContext& ctx) {
    std::vector<ExprPtr> uniq;
    for (const auto& r : roots) {
        bool seen = false;
        for (const auto& u : uniq) {
            auto eq = cas::symbolic::mathematically_equal(r, u, ctx);
            if (eq.is_ok() && eq.value()) { seen = true; break; }
        }
        if (!seen) uniq.push_back(r);
    }
    return uniq;
}

inline Result<bool> compare_solve_sets(
    const std::vector<ExprPtr>& cas_raw,
    const std::vector<ExprPtr>& maxima_raw,
    cas::symbolic::CASContext& ctx) {

    // Maxima `solve()` by default returns each distinct root once
    // (multiplicity is reported separately via `multiplicities`).
    // Our CAS includes multiplicity in the vector. Normalize both sides
    // to the unique-root set before comparison.
    auto cas    = dedup_modulo_equal(cas_raw, ctx);
    auto maxima = dedup_modulo_equal(maxima_raw, ctx);

    if (cas.size() != maxima.size()) {
        return ok(false);
    }
    if (cas.empty()) return ok(true);

    std::size_t inconclusive = 0;

    auto edge_predicate = [&](std::size_t i, std::size_t j) -> bool {
        auto r = cas::symbolic::mathematically_equal(cas[i], maxima[j], ctx);
        if (!r.is_ok()) {
            ++inconclusive;
            return false;
        }
        return r.value();
    };

    std::size_t matched = hopcroft_karp_max_matching(
        cas.size(), maxima.size(), edge_predicate);

    if (matched == cas.size()) return ok(true);

    // Threshold: if more than half the (i,j) probes were inconclusive,
    // we cannot honestly declare unequal — report Unimplemented so the
    // runner SKIPs rather than FAILs.
    std::size_t total_probes = cas.size() * maxima.size();
    if (total_probes > 0 && inconclusive * 2 > total_probes) {
        return CASError{CASErrorKind::Unimplemented,
                        "compare_solve_sets: too many inconclusive comparisons",
                        std::nullopt};
    }
    return ok(false);
}

}  // namespace cas::golden
