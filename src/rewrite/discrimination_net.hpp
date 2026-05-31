#pragma once
// F1.5 — Discrimination Net for fast rewrite rule lookup.
//
// A discrimination net (Forgy 1982, "RETE: A Fast Algorithm for the Many
// Pattern / Many Object Pattern Match Problem", Artificial Intelligence 19)
// is a trie-like structure that avoids testing each rule independently.
// Here we implement a two-level net:
//
//   Level 0: ExprKind of pattern root.
//   Level 1 (FuncCall only): function name (string key).
//
// Lookup: given an expression, retrieve only rules whose pattern root kind
// matches the expression's kind (and function name for FuncCall).
//
// Complexity:
//   Insert: O(1) amortised.
//   Lookup: O(|candidates|) where |candidates| << |all_rules| for
//           large rule sets with diverse pattern kinds.
//
// Pattern type filters supported (via is_wildcard_name convention):
//   x_         — any expression (plain wildcard)
//   x_Integer  — matches IntegerLit only
//   x_Positive — matches expressions known to be positive literals
//   x_Symbol   — matches Symbol only
//
// AC matching: pre-sort of commutative children is handled by the existing
// match_sequence_ac_internal (Hoffmann-O'Donnell style backtracking with
// expr_weight ordering).  The net does not duplicate that logic.
//
// Usage:
//   DiscriminationNet net;
//   net.insert(rule);                  // add RewriteRule
//   auto cands = net.lookup(expr);     // returns const ref to candidates
//   for (const RewriteRule* r : cands) { /* try match */ }

#include "cas/symbolic.hpp"
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace cas::rewrite {

// ── Typed wildcard filter ─────────────────────────────────────────────────────

// Constraint carried by a typed wildcard pattern symbol name.
enum class WildcardType : std::uint8_t {
    Any,      // x_          — unconstrained
    Integer,  // x_Integer   — IntegerLit only
    Positive, // x_Positive  — literal positive
    Symbol,   // x_Symbol    — Symbol only
};

// Parse the typed wildcard constraint from a symbol name such as "x_Integer".
// Returns WildcardType::Any for plain wildcards ending in '_' with no suffix.
[[nodiscard]] WildcardType wildcard_type_from_name(const std::string& name) noexcept;

// Test whether an expression satisfies a WildcardType constraint.
[[nodiscard]] bool expr_satisfies_wildcard_type(ExprPtr expr, WildcardType wt) noexcept;

// ── DiscriminationNet ─────────────────────────────────────────────────────────

class DiscriminationNet {
public:
    DiscriminationNet() = default;

    // Insert a rule into the net.
    // The rule is indexed by (pattern root ExprKind) and, for FuncCall patterns,
    // additionally by the function name.
    // Rules whose pattern is a plain wildcard (Symbol ending in '_') are placed
    // in the universal bucket (matched against every expression).
    void insert(const symbolic::RewriteRule& rule);

    // Retrieve candidate rules for the given expression.
    // Returns a reference to an internal vector; valid until next insert().
    // The caller must still perform full pattern matching on each candidate.
    [[nodiscard]] const std::vector<const symbolic::RewriteRule*>&
    lookup(ExprPtr expr) const noexcept;

    // Number of distinct rules stored.
    [[nodiscard]] std::size_t size() const noexcept { return total_rules_; }

    // Remove all rules.
    void clear() noexcept;

private:
    // Per-kind buckets.  kNumKinds covers all ExprKind values (0 .. 17).
    static constexpr std::size_t kNumKinds =
        static_cast<std::size_t>(ExprKind::Quantity) + 1U;

    std::vector<const symbolic::RewriteRule*> kind_buckets_[kNumKinds];

    // FuncCall patterns additionally keyed by function name.
    std::unordered_map<std::string,
                       std::vector<const symbolic::RewriteRule*>> funcall_buckets_;

    // Wildcard (Symbol "x_") patterns matched against every expression.
    std::vector<const symbolic::RewriteRule*> universal_bucket_;

    // Empty vector returned when no candidates exist.
    std::vector<const symbolic::RewriteRule*> empty_;

    // Scratch buffer for lookup() — avoids heap allocation on hot path.
    mutable std::vector<const symbolic::RewriteRule*> scratch_;

    std::size_t total_rules_{0U};

    // True if pattern root is a plain wildcard Symbol (ends in '_').
    [[nodiscard]] static bool is_universal_wildcard(ExprPtr pattern) noexcept;
};

} // namespace cas::rewrite
