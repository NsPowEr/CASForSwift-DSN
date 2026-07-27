#pragma once
// symbolic_api.hpp — funzioni libere di `cas::symbolic` (rewriting, matching,
// simplify/substitute/equal).
//
// NON includere direttamente: e' incluso in fondo a `cas/symbolic.hpp`, dopo
// le definizioni di `CASContext`, `RewriteRule` e `MatchMap` da cui queste
// firme dipendono. Estratto da symbolic.hpp per il limite anti-monolito
// (500 righe): il file era a 486 e ogni aggiunta lo sforava.

namespace cas::symbolic {

[[nodiscard]] int canonical_compare(ExprPtr lhs, ExprPtr rhs) noexcept;

// F7.0-A4.2: post-simplify canonical-form invariant check.
// Returns true if `expr` and every reachable sub-expression respect the
// invariants that simplify() is supposed to maintain:
//   - Sum::terms sorted by polynomial degree descending then canonical_compare,
//     no nested Sum, no exact-zero IntegerLit / RationalLit summand.
//   - Product::factors sorted by canonical_compare, no nested Product, no
//     exact-one IntegerLit / RationalLit factor.
//   - All Sum/Product nodes have ≥ 2 operands (singletons collapsed).
//
// Used in DEBUG builds via an assert at the end of CASContext::simplify()
// to catch invariant violations close to their source. In release builds
// the function is still available for explicit verification but the assert
// is compiled out.
[[nodiscard]] bool is_strictly_canonical(ExprPtr expr) noexcept;
[[nodiscard]] TermOrderRelation compare_rewrite_terms(ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] bool rewrite_rule_is_oriented(const RewriteRule& rule);
[[nodiscard]] bool is_strongly_normalizing(const std::vector<RewriteRule>& rules);
[[nodiscard]] bool match_pattern(ExprPtr expr, ExprPtr pattern, MatchMap& out_matches);
[[nodiscard]] bool match_ac_pattern(ExprPtr expr, ExprPtr pattern, MatchMap& out_matches);
[[nodiscard]] Result<ExprPtr> apply_rule(ExprPtr expr, const RewriteRule& rule, TraversalStrategy strategy, AstArena& arena);
[[nodiscard]] Result<ExprPtr> apply_rule_set(ExprPtr expr, const std::vector<RewriteRule>& rules, AstArena& arena);
[[nodiscard]] Result<ExprPtr> materialize_expr(ExprPtr expr, AstArena& arena);
[[nodiscard]] const RewriteProvider& default_rewrite_provider();
[[nodiscard]] Result<ExprPtr> simplify(ExprPtr expr, AstArena& arena);
[[nodiscard]] Result<ExprPtr> simplify(ExprPtr expr, CASContext& context);
[[nodiscard]] Result<ExprPtr> substitute(ExprPtr expr, const Symbol& variable, ExprPtr value, CASContext& context);
[[nodiscard]] Result<bool> mathematically_equal(ExprPtr lhs, ExprPtr rhs, CASContext& context);

// L2-19: decidable subset of transcendental equivalence via Risch-style
// log/exp/trig normalisation.  Applies opt-in expansion rules to lhs/rhs
// (log(x*y) -> log(x)+log(y) under x>0, y>0; exp(x+y) -> exp(x)*exp(y);
// exp(n*ln(x)) -> x^n under x>0; sin^2+cos^2 collapse) before delegating
// to mathematically_equal.  Returns false (not Unimplemented) for cases
// outside the decidable subset — Richardson's theorem precludes a total
// decision procedure.
[[nodiscard]] Result<bool> mathematically_equal_subset_risch(
    ExprPtr lhs, ExprPtr rhs, CASContext& context);

}  // namespace cas::symbolic
