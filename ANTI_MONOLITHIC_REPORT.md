# Anti-Monolithic Report — CAS Engine C++ — 2026-06-11

> Generato post F7.5.H2 audit T3-Opus (vedi `AUDIT_CAS_F7.5_2026-06-11.md`
> §4) per soddisfare condition **C1 — Anti-monolith split or formal
> waiver**. Sostituisce il gate informale per il limite 500 LOC con un
> piano di chiusura formale tracciato per Fase 8.

## Status

28 file `.cpp`/`.hpp` superano il limite 500 LOC dichiarato in
CLAUDE.md "STANDARD TECNICI E ANTI-DEBITO". Nessuno è stato splittato
durante F7.5 perché:

1. Tutti i file >500 LOC contengono algoritmi coesi (es. `simplify_arithmetic.cpp`
   = simplifier core arithmetic + Pow + Sum normalization) il cui split
   richiede analisi semantica approfondita per evitare circular include
   o duplicazione di dichiarazioni interne.
2. F7.5 ha priorità per chiusura aggregato corpus (94.5% raggiunto) +
   F7.5.F1 Extended-Real Phase 2 (HC-F70-A43 chiuso).
3. Split prematuro durante F7.5.B/C/D in corso introdurrebbe merge
   conflict massicci con le PR research-grade pending (F7.5.B2 Hermite,
   F7.5.B3 Risch transcendental).

## Waiver formale: HC-F8-MONOLITH-WAIVER

Iscritto in `HARDCODE_LEDGER.md` come `HC-F8-MONOLITH-WAIVER`:
- 14 file `>600` LOC = **tier-1**: split obbligatorio prima di Fase 8 entry
  per qualsiasi PR research-grade (Risch, Galois ≥6, CAD).
- 14 file `500-600` LOC = **tier-2**: split tollerato durante Fase 8
  migration, da chiudere entro fine Fase 8.

Nessun nuovo file `>500` LOC permesso post 2026-06-11 senza ledger
entry esplicita.

## Tier-1 (>600 LOC) — split obbligatorio inizio Fase 8

| LOC | File | Cohesion-based split plan |
|---|---|---|
| 833 | `src/symbolic/simplify_arithmetic.cpp` | → `simplify_arithmetic.cpp` (Binary dispatch + helpers, ~350) + `simplify_arithmetic_pow.cpp` (Pow expansion + Chebyshev/DeMoivre, ~300) + `simplify_arithmetic_sum_collapse.cpp` (Sum collapse + canonicalization, ~180) |
| 753 | `src/ast/ast.cpp` | → `ast.cpp` (constructor + arena allocation, ~300) + `ast_compare.cpp` (canonical_compare + structural_equal, ~250) + `ast_hash.cpp` (expr_hash + interning helpers, ~200) |
| 713 | `src/symbolic/assumptions.cpp` | → `assumptions.cpp` (storage + query API, ~350) + `assumptions_inference.cpp` (inference rules positivity/nonzero/interval, ~360) |
| 701 | `include/cas/ast.hpp` | → header split via include chain: `ast_kinds.hpp` (enum + forward decl, ~150) + `ast_nodes.hpp` (struct defs, ~300) + `ast.hpp` (umbrella include, ~80) |
| 693 | `src/calculus/differentiate.cpp` | → `differentiate.cpp` (Symbol + Binary dispatch, ~300) + `differentiate_func.cpp` (FuncCall dispatch tabella derivate elementari, ~250) + `differentiate_special.cpp` (Bessel/Gamma/Erf/HyperGeometric, ~150) |
| 661 | `src/symbolic/simplify_exp_log.cpp` | → `simplify_exp_log.cpp` (exp/ln base rules, ~340) + `simplify_exp_log_denest.cpp` (Landau denesting + Cardano, ~320) |
| 658 | `src/symbolic/context_core.cpp` | → `context_core.cpp` (CASContext lifecycle + arena, ~320) + `context_params.cpp` (configurable params accessors, ~340) |
| 658 | `src/calculus/residue_theorem.cpp` | → `residue_theorem.cpp` (entry + Jordan lemma dispatch, ~330) + `residue_theorem_poles.cpp` (pole-order arbitrary handling, ~330) |
| 644 | `src/calculus/ode_solver_frobenius.cpp` | → `ode_solver_frobenius.cpp` (indicial + recurrence, ~340) + `ode_solver_frobenius_log.cpp` (log-term case differing-by-integer, ~300) |
| 642 | `src/algebra/algebraic_tower_primitive_internal.hpp` | → `algebraic_tower_primitive_internal.hpp` (Trager shift + composite norm, ~350) + `algebraic_tower_primitive_lift.hpp` (factor lift in tower arithmetic, ~290) |
| 620 | `src/algebra/algebraic_tower_bridge.cpp` | → `algebraic_tower_bridge.cpp` (RootOf↔Q(α) conversion, ~310) + `algebraic_tower_bridge_polynomial.cpp` (polynomial-in-α simplify, ~310) |
| 615 | `src/rewrite/builtin_rewrite.cpp` | → `builtin_rewrite.cpp` (dispatcher + caching, ~310) + `builtin_rewrite_rules.cpp` (rule table per BuiltinOp, ~305) |
| 609 | `src/foundation/bigint.cpp` | → `bigint.cpp` (representation + I/O, ~300) + `bigint_div.cpp` (Knuth D long division + helpers, ~309) — moltiplicazione già splittata su `bigint_mul_*.cpp` |
| 604 | `src/symbolic/term_order.cpp` | → `term_order.cpp` (canonical_compare dispatcher, ~310) + `term_order_func.cpp` (FuncCall ordering rules, ~294) |

## Tier-2 (500-600 LOC) — split tollerato durante Fase 8

| LOC | File | Note split (post-tier-1) |
|---|---|---|
| 596 | `src/algebra/polynomial_gcd_brown_modular.cpp` | Brown core + CRT lifting separabili |
| 584 | `src/algebra/algebraic_number_bridge.cpp` | API + conversion separabili |
| 577 | `src/calculus/differential_field.cpp` | Visitor + field ops separabili |
| 573 | `src/algebra/solve_polynomial.cpp` | Cardano/Ferrari + quadratic separabili |
| 572 | `src/calculus/limit.cpp` | Engine + dispatch separabili |
| 572 | `src/algebra/polynomial_groebner_fglm.cpp` | BFS loop + shape-lemma separabili |
| 566 | `src/algebra/polynomial_gcd_fp_recursive.cpp` | Recursive GCD + base case separabili |
| 562 | `src/calculus/summation_abramov.cpp` | Algoritmo Abramov + helpers |
| 551 | `src/linalg/matrix_smith.cpp` | Smith Z + Smith Q[x] separabili |
| 551 | `include/cas/cas_context_params.hpp` | Params header — split per dominio |
| 545 | `src/symbolic/summation_zeilberger_helpers.cpp` | Telescoping + rational summation |
| 541 | `src/algebra/polynomial_groebner_f5.cpp` | F5 + signatures |
| 532 | `src/calculus/laplace_transform.cpp` | Forward + inverse |
| 521 | `src/calculus/orthogonal_polynomials.cpp` | Chebyshev + Hermite + Legendre |

## Regola operativa enforcement

1. CI gate `scripts/check_file_size.sh` (PLAN_HP_PRIME_PARITY.md F0.6)
   resta attivo. Solo i file in waiver tier-1+tier-2 sono esentati,
   tutti gli altri devono passare ≤500 LOC.
2. Whitelist documentata in `scripts/check_file_size.sh` con
   commento riferimento a questo report.
3. Ogni nuovo split deve includere update di questa lista (rimuovere
   il file splittato) + verifica build + suite quick green.

## Ownership

Tier-1 split: F8.0 prerequisite block — T1-Sonnet (split meccanico),
~3-5 giorni. Tier-2 split: durante Fase 8 mainstream, T1-Sonnet
opportunistico, ~1-2 settimane spalmate.

## Sign-off audit C1

Waiver formale documentata + ledger entry creata + piano split
dettagliato per ciascun file. Condition C1 audit `AUDIT_CAS_F7.5_2026-06-11.md`
**risolta come waiver** (path B della raccomandazione, riga 261-263).
