# SESSION HANDOFF — F5 Closure Multi-Session

> **Data ultimo aggiornamento**: 2026-06-03
> **Sessione**: B1 + B2 (a,b,c) + B3 completate; pre-B2b cleanup chiuso (Regola1 simplify_exp_log + restore QR 8x8 verification + RootOf dispatch eigen/jordan); B4-B9 da eseguire
> **Orchestrator**: Opus 4.7
> **Modalità utente**: caveman ultra (interno) + no shortcut + no hardcode

---

## OBIETTIVO MACRO

Chiudere **definitivamente F5 (Strato L2 Calculus)** senza hardcode, senza scorciatoie, senza debiti nascosti. Multi-session OK. Ogni blocco chiude il proprio sotto-item completamente.

---

## STATO BLOCCHI

### ✅ B1 — Quick wins / ledger cleanup (CHIUSO 2026-06-02)

| HPP | File | Categoria | Fix |
|---|---|---|---|
| HPP-005 | `src/calculus/integrate.cpp:38-87, :244-275` | Regola 1 / Cat 4 | `approx_bound` ora `BigFloat` (256-bit MPFR), `cos_zero_in_range` tutto BigFloat. CMake: MPFR include path early-find aggiunto a cas_calculus |
| HPP-008 | `src/calculus/ode_solver_advanced.cpp:228-230` | Cat 7 | `ctx.make_fresh_symbol("C")` |
| HPP-009 | `src/calculus/ode_solver_1st_order.cpp:37-38, :59-60` | Cat 7 | Già pre-fixato, ledger aggiornato |
| HPP-010 | `src/linalg/matrix_solve.cpp:137-143` | Cat 7 | `ctx.make_fresh_symbol("c")` |
| HPP-011 | `src/linalg/matrix_eigenvalues.cpp:191-193` | Cat 7 | `ctx.make_fresh_symbol("lambda")` |
| HPP-012 | `src/linalg/matrix_jordan.cpp:96-98` | Cat 7 | `ctx.make_fresh_symbol("lambda")` |
| HPP-015 | `src/symbolic/simplify_special_fn.cpp` (5 siti) | Cat 1+4 | Aggiunti 2 param in `cas_context_params.hpp`: `max_special_fn_integer_arg_bits_` (def 16) + `max_bernoulli_index_bits_` (def 30). 3 test in `test_special_functions.cpp` |

**Verifica B1**: 74/74 SpecialFn + 65/65 linalg + 14/14 ODE + 69/69 Integrate. Zero regressioni.

Nota: `AcidTest.Test5_ExpansionStress` (986ms vs 500ms soglia) era già rosso pre-fix (verificato via `git stash` + rerun). Non blocking.

---

### ✅ B2 — F5.3 ODE nonlinear 1st-order + Frobenius log-term (CHIUSO 2026-06-03)

**B2a Riccati** (commit `12e6982`):
- Classifier via Lagrange three-point fit + shielded substitution.
- Solver: Bernoulli reduction (q_0=0) + closed-form for constant coefficients
  (Δ>0/=0/<0 → tanh / rational / tan branches).
- Variable Riccati senza particolare → Unimplemented diagnostico (no HC).

**B2b Clairaut / d'Alembert** (commit `b535bab`):
- `try_lagrange` in `ode_classifier.cpp`: estrae c_1 (y-coeff) e c_0 = α(y')·x + β(y'),
  normalizza F(p) = -α/c_1, G(p) = -β/c_1 con p fresh symbol.
- Clairaut (F ≡ p): `Equal(y, C·x + G(C))` + parametrica singolare
  `x = -G'(p), y = G(p) - p·G'(p)` → output `FuncCall("GeneralAndSingular", [...])`.
- d'Alembert (F ≢ p): riduzione a `dx/dp + P·x = Q`, solve via `solve_ode_1st_order`,
  output `FuncCall("ParametricSolution", [Equal(x,X(p,C)), Equal(y,X·F+G)])`.
- Branch convention: parameter p marcato positivo in assumptions per ramo principale.

**B2c Frobenius log-term** (commit `d9fa86b`):
- `integer_gap(r_1, r_2, ctx)`: detect r_1 - r_2 ∈ Z_{>0}.
- `build_log_branch`: recurrence completa con risoluzione del coefficiente c
  alla resonance step (n = N). h_m precomputato per coupling c·y_1.
- y_2 = c·ln(x)·y_1 + x^{r_2}·Σ b_n x^n, con c = -S_N / h_0 e h_0 = I'(r_1) = N.
- Skip se num_terms < N (resonance non raggiunta nella troncatura).
- Secondary resonance al n > N → Unimplemented diagnostico.

**Output test**:
- 4/4 OdeTest Riccati (B2a) + 2/2 OdeTest Clairaut/d'Alembert (B2b)
- 4/4 FrobeniusTest (Euler, rational, ±1, **BesselOrder1_ResonantLogBranch**)
- 1950/1950 non-stress + 24/24 AcidTest + 3/3 SupremeStressTest + QR 8x8 6.6s

**Riferimenti**:
- Ince, *Ordinary Differential Equations* (1926) §§4.4, 4.6, 5.3
- Coddington-Levinson 1955 §4.8 (Frobenius II)
- Bronstein-Schiermayer 1998 ODE classification

---

### 🟡 B3 — F5.2 Gruntz dynamic growth rank (DA INIZIARE)

**Scope**:
- Sostituire `GrowthRank` statico in `limit_infinite.cpp` con confronto asintotico dinamico via algoritmo Gruntz §3.5.
- Refactor `compare_growth()` in `limit_mrv.cpp:115-160` per usare leading coefficient e ricorsione MRV set.
- Cancellation tower con livelli arbitrari.

**Effort stimato**: ~500 LOC. Opus required (research-grade).

**Riferimenti**:
- Gruntz, *On Computing Limits in a Symbolic Manipulation System* (PhD ETH 1996) §3.5

---

### 🟡 B4 — F5.6 Residue theorem deg≥5 generico (DA INIZIARE)

**Sub-block B4a — Complex root isolator (Aberth/Pan)**:
- Nuovo modulo `src/numeric/complex_root_isolator.cpp`.
- Aberth iteration con BigFloat (MPFR-C++ NTL-style complex).
- Bound errore tramite Gerschgorin + Smale α-theory.

**Sub-block B4b — Residue driver deg≥5**:
- Estendere `integrate_rational_full_real_line` in `residue_theorem.cpp:352`.
- Per ogni fattore irreducibile `m(x)` ∈ Q[x] grado ≥3 senza radici reali:
  - Isola complex roots via B4a
  - Filtra Im(α) > 0
  - Per ogni upper-half root: compute Res(N/D, α) ∈ Q(α) via `residue()`
  - Somma `2πi · ΣRes`
  - Verifica imaginary cancellation simbolicamente

**Effort stimato**: ~1500 LOC. Multi-session. Opus required.

**Riferimenti**:
- Pan, *Univariate Polynomials: Nearly Optimal Algorithms for Numerical Factorization* (1997)
- Bronstein, *Symbolic Integration I* §5.6

---

### 🟡 B5 — F5.5 Puiseux series + Padé non-Q (DA INIZIARE)

**Sub-block B5a — Puiseux series via Newton polygon**:
- Algoritmo Newton polygon per branch point algebrici.
- Iterative Hensel su Q[[t^{1/n}]].

**Sub-block B5b — Padé non-Q**:
- Estendere `pade.cpp` per supportare Taylor con coeff algebriche/transcendentali (π, e, sqrt) via Q(α)[[t]].

**Effort stimato**: ~1500 LOC totali. Multi-session.

**Riferimenti**:
- Walker, *Algebraic Curves* (Princeton 1950) Ch.V
- Brent-Kung, *Fast algorithms for manipulating formal power series* (1978)

---

### 🟡 B6 — F5.7 Zeilberger/WZ + Abramov (DA INIZIARE)

**Scope**:
- Zeilberger creative telescoping per somme ipergeometriche.
- WZ pair method (Petkovšek-Wilf-Zeilberger 1996).
- Abramov rational summation per chiusura closed-form.

**Effort stimato**: ~2000 LOC. Multi-session. Opus required.

**Riferimenti**:
- Petkovšek-Wilf-Zeilberger, *A = B* (1996), Ch.4-6
- SymPy `sympy.concrete.summations` (citazione algoritmica, non copia)

---

### 🟡 B7 — F5.8 Fourier/Mellin/Z transforms + Inv Laplace via Bronstein (DA INIZIARE)

**Sub-block B7a — Fourier transform**: tabella + contour per rationale.
**Sub-block B7b — Mellin transform**: tabella + contour residue.
**Sub-block B7c — Z-transform**: tabella + algoritmico per sequenze rationale.
**Sub-block B7d — Inverse Laplace**: Bronstein residue method per `F(s)` rationale.

**Effort stimato**: ~3000 LOC totali. Multi-session.

**Riferimenti**:
- Abramowitz-Stegun, *Handbook of Mathematical Functions* §29 (Laplace)
- Bronstein 1992, *Integration of Elementary Functions* (PhD)

---

### 🟡 B8 — F5.9 Hypergeometric pFq + Special functions completi (DA INIZIARE)

**Sub-block B8a — pFq(a;b;z) first-class**: nuovo AST node + contiguous relations + Gauss 2F1(a,b;c;1) + Saalschütz.
**Sub-block B8b — Bessel/Erf/Zeta complete**: closed-form recognition rules.
**Sub-block B8c — Jacobi P_n^(α,β) + Elliptic K/E/Π/F**: builtin + simplification.

**Effort stimato**: ~3500 LOC. Multi-session.

**Riferimenti**:
- Andrews-Askey-Roy, *Special Functions* (CUP 1999)
- Olver, *NIST Handbook of Mathematical Functions* §§13-19

---

### 🟡 B9 — F5.1 Risch Bronstein cap.7-9 + HPP-007 (DA INIZIARE)

**Sub-block B9a — cap.7 parametric problems**: ParametricLogarithmicDerivative + PolynomialReduce.
**Sub-block B9b — cap.8 Risch DE general**: Y' + f·Y = g in differential field.
**Sub-block B9c — cap.9 structure theorem**: detect algebraic dependencies between generators.
**Sub-block B9d — HPP-007 closure**: dopo cap.9, trial constants `{±1,±1/2,±2}` in `integrate_risch.cpp:623` rimpiazzati da `c = integrand/D(F)` via residue field equation.

**Effort stimato**: ~5000 LOC. Multi-session (research grade, top priority HP Prime parity blocker).

**Riferimenti**:
- Bronstein, *Symbolic Integration I — Transcendental Functions* (Springer 2005) Ch.7-9
- Davenport, *On the Risch Differential Equation Problem* (SIGSAM 1981)

---

## RIPARTENZA — COMANDO

1. `/model opus`
2. Apri questo file
3. Comando utente: `"riprendi handoff B<n>"` (es. `"riprendi handoff B2"`)
4. Orchestrator:
   - Legge sezione B<n> qui sopra
   - Legge file `src/calculus/<modulo>` di riferimento
   - Apre task list (se presente) o crea
   - Procede in modalità no-hardcode con commit-per-fix granulare

---

## REGOLE PERMANENTI (estratto CLAUDE.md)

- **REGOLA 0**: no shortcut, no hardcode-of-passage non documentato in `HARDCODE_LEDGER.md`
- **Regola 1**: zero `int64_t`/`double` nel core simbolico (use `BigInt`, `Rational`, `BigFloat`)
- **Anti-loop**: 3 tentativi falliti → STOP + report stallo (no rewrite cieco)
- **Anti-monolito**: ≤500 LOC per file sorgente
- **Test**: solo confronto strutturale o equivalenza matematica
- **Trace**: `D(integrate(f)) ≡ f` per ogni integral closure

---

## FILE TOCCATI IN SESSIONE B1 (2026-06-02)

- `include/cas/cas_context_params.hpp` — 2 nuovi param (HPP-015)
- `src/symbolic/simplify_special_fn.cpp` — 5 bail-out configurabili
- `src/calculus/ode_solver_advanced.cpp` — fresh symbol Ci
- `src/linalg/matrix_solve.cpp` — fresh symbol c
- `src/linalg/matrix_eigenvalues.cpp` — fresh symbol lambda
- `src/linalg/matrix_jordan.cpp` — fresh symbol lambda
- `src/calculus/integrate.cpp` — approx_bound BigFloat + cos_zero_in_range BigFloat
- `CMakeLists.txt` — MPFR_INCLUDE_DIR_EARLY per cas_calculus
- `test/unit/symbolic/test_special_functions.cpp` — 3 nuovi test HPP-015
- `HARDCODE_LEDGER.md` — 7 voci chiuse (HPP-005, 008, 009, 010, 011, 012, 015)
- `SESSION_HANDOFF.md` — questo file riscritto

Nessun commit creato in B1 (l'utente deciderà al resume).

---

## COMMIT SUGGERITO B1 (se utente lo richiede)

```
feat(F5/B1): close 7 HPP entries — fresh symbols + Bernoulli/special-fn budget + BigFloat singularity check

- HPP-005: integrate.cpp approx_bound + cos_zero_in_range now BigFloat (256-bit MPFR);
  eliminates double round-trips and to_u64() truncation that violated CLAUDE.md
  Regola 1 in the core symbolic path.  CMake: MPFR include path early-find added
  for cas_calculus.
- HPP-008: ode_solver_advanced.cpp Ci constant uses ctx.make_fresh_symbol("C").
- HPP-009: ode_solver_1st_order.cpp verified pre-fixed; ledger updated.
- HPP-010: matrix_solve.cpp free-variable symbol uses ctx.make_fresh_symbol("c").
- HPP-011: matrix_eigenvalues.cpp lambda uses ctx.make_fresh_symbol("lambda").
- HPP-012: matrix_jordan.cpp lambda uses ctx.make_fresh_symbol("lambda").
- HPP-015: simplify_special_fn.cpp bit_length budgets configurable via two new
  CASContextParams (max_special_fn_integer_arg_bits, max_bernoulli_index_bits).
  Five bail-out sites updated.  Three new tests in test_special_functions.cpp.

Verification: 74/74 SpecialFn + 65/65 linalg + 14/14 ODE + 69/69 Integrate PASS.
AcidTest.Test5_ExpansionStress regression is pre-existing (986ms vs 500ms threshold);
verified via git-stash bisection.

Prepares F5/B2 (ODE Riccati/Clairaut/d'Alembert + Frobenius log-term).
```
