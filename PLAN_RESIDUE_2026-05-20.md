# Piano Implementativo — Task Aperte Residue + Chiusura Debiti
> Data: 2026-05-20 (sessione post STEP 34)
> Vincoli: CLAUDE.md REGOLA ZERO + ANTI-HARDCODE + L0→L1→L2→L3
> Stato test: 21/21 LaplaceTest + OdeLaplaceTest, -Werror pulita

---

## A. Inventario Debiti Rilevati (19 voci)

### A.1 Test DISABLED (5)
| # | Test | Task | Causa |
|---|---|---|---|
| D1 | `ComplexLogBranchTest.LnOfOnePlusI` | L2-17 | `exp(ln(z)) = z` su z=1+i: simplifier non canonicalizza abs(1+i)→√2 + arg(1+i)→π/4 in forma roundtripable |
| D2 | `OdeLaplaceTest.FirstOrderForced` | L3-07/10 | Y(s)=1/(s(s+1)): PFD ritorna termini ma inverse Laplace bridge fragile |
| D3 | `EquivalenceSubsetRischTest.ExpOfLogSum` | L2-19 | Branch-cut strict pre-esistente |
| D4-5 | `FactorizationTowerTest.SplitsX4Minus...` | L3-06 | Tower deg-4 perf ASan |

### A.2 GTEST_SKIP (5)
| # | Test | Task | Causa |
|---|---|---|---|
| S1 | `test_residue_theorem.cpp:138` | L2-22 | Double-pole Q(α) reduction |
| S2 | `test_property_based.cpp:168` | meta | factor random failure tolerance |
| S3-6 | `test_factorization_trager.cpp` x4 | L3-06 | Galois extension edge cases |

### A.3 Code TODO / Follow-up (4)
| # | File | Task | Causa |
|---|---|---|---|
| T1 | `laplace_transform.cpp:13` | L3-07 | Inverse Laplace PFD bridge marginale |
| T2 | `matrix_lu.cpp:85` | L3-17 | Zero pivot → permuted LU mancante |
| T3 | `galois.cpp:10` | L3-18 | Deg 4 V4/S4 conservative — A4/D4/C4 richiede resolvent cubic |
| T4 | gcd_probabilistic | L3-15 | Wrapper thin, no Brown's algorithm reale |

### A.4 Architettura / Wiring (5)
| # | Voce | Task | Causa |
|---|---|---|---|
| A1 | `strict_branch_cuts_` flag exposed | L2-21 | Non wired in simplify rules |
| A2 | Interval sin/cos critical-point | L3-13 | Conservative width>π → [-1,1] |
| A3 | Quantity^non-int exponent | L3-08 | Solo integer power gestito |
| A4 | Quantity·non-Quantity Mul | L3-08 | Coefficienti scalari non-Quantity richiedono regola |
| A5 | `max_denesting_depth` CASContext | L1-12 | Cascade tramite simplify_expr — no explicit field |

---

## B. Piano Task Aperte + Integrazione Debiti

### Ordinamento ottimale (logica dipendenze + ROI)

```
Fase 1 (Chiusura debiti — abilita test future):
  STEP 35: L3-17 permuted LU (T2)          ~3h
  STEP 36: L3-18 Galois deg 4 advanced (T3) ~5h
  STEP 37: L3-15 Brown modular GCD MVP (T4) ~6h

Fase 2 (Estensioni matematiche):
  STEP 38: L3-04 JacobiP polynomial         ~4h
  STEP 39: L2-16 Cambio variabile auto      ~6h
  STEP 40: L3-14 Hensel multivariate MVP    ~12h
  STEP 41: L3-17 QR via Gram-Schmidt        ~6h
  STEP 42: L3-09 FGLM zero-dim MVP          ~10h

Fase 3 (Research / Heavy):
  STEP 43: L3-05 Zeilberger Gosper          ~12h
  STEP 44: L2-18 Sparse interp Ben-Or-Tiwari ~8h
  STEP 45: L3-10 Frobenius advanced         ~10h
  STEP 46: L3-02 CAD Collins MVP            ~30h (research-grade)

Cleanup (fix DISABLED ove possibile):
  STEP 47: L3-07 PFD-inverse-Laplace bridge harden (D2, T1) ~3h
  STEP 48: L2-21 strict mode wiring (A1)    ~4h
  STEP 49: L3-13 sin/cos critical point (A2) ~3h
  STEP 50: L3-08 Quantity coverage gaps (A3, A4) ~3h
```

**Effort totale stimato**: ~130h. Suddivisibile in micro-step da committare singolarmente.

---

## C. Dettaglio Implementativo per Step

### STEP 35 — L3-17 Permuted LU (T2)
**Diagnosi**: `lu_decompose` fallisce su matrici con A[k][k]=0 senza pivoting. Esempio: `[[0,1],[1,0]]` errore.

**Algoritmo**: Aggiungere `lu_decompose_with_pivoting` che ritorna `(P, L, U)` con `P·A = L·U`. Partial pivoting: ad ogni step k, scegliere riga `argmax_i |U[i][k]|` (`PivotScore` symbolic già esistente). Permutare righe.

**Files**:
- `include/cas/linalg/Matrix.hpp`: aggiungi struct `PLUDecomposition{P_perm, L, U}` + decl `lu_decompose_pivoted(M)`.
- `src/linalg/matrix_lu.cpp`: implementazione (~60 LOC).

**Test**:
- `[[0,1],[1,0]]` produce P swap + L·U = swapped.
- `[[2,3],[4,7]]` no swap (no zero pivot) → identico a lu_decompose.
- Anti-hardcode: 3 random matrici con zero entries, verify P·A = L·U.

### STEP 36 — L3-18 Galois deg 4 advanced (T3)
**Algoritmo (Lazard 1990 / Buhler-Reverter)**:
1. Depress quartic: monic `x⁴+ax³+bx²+cx+d` → `y⁴+py²+qy+r` via `y=x-a/4`.
2. Resolvent cubic: `R(z) = z³ - 2p·z² + (p²-4r)·z + q²`.
3. Factor R(z) over Q via solve_polynomial.
4. Dispatch:
   - R reducibile + 3 rational roots + disc(quartic) square → V4
   - R reducibile + 1 rational root + ... → D4
   - R reducibile + 1 rational root + condition C → C4
   - R irreducibile + disc square → A4
   - R irreducibile + disc non-square → S4

**Files**: `src/algebra/galois.cpp` extend (~80 LOC).

**Test** (5 nuovi):
- x⁴+1 (Φ₈) → V4
- x⁴-2 → D4
- x⁴+x²+1 → D4
- x⁴-10x²+1 (sqrt 2 + sqrt 3) → V4
- x⁴+x+1 → S4

### STEP 37 — L3-15 Brown's modular GCD MVP (T4)
**Algoritmo (Brown 1971)**:
1. Sel primi p_1, ..., p_k coprime al lc.
2. Per ogni p: compute gcd_{Z_p}[x] via Euclidean.
3. CRT lift back a Z[x] usando bound Mignotte.
4. Verify certified divisibility.

**Files**: `src/algebra/polynomial_gcd_modular.cpp` riscrivi (~200 LOC).

**Test**: 5 esempi univariati noti.

### STEP 38 — L3-04 JacobiP polynomial
**Algoritmo (Bonnet recurrence)**:
- `P_0^{(α,β)}(x) = 1`
- `P_1^{(α,β)}(x) = (α-β)/2 + (α+β+2)·x/2`
- `2(n+1)(n+α+β+1)(2n+α+β) · P_{n+1}^{(α,β)}(x) = ((2n+α+β+1)·(α²-β²) + (2n+α+β)·(2n+α+β+1)·(2n+α+β+2)·x) · P_n^{(α,β)}(x) - 2(n+α)(n+β)(2n+α+β+2) · P_{n-1}^{(α,β)}(x)`

**Files**:
- `include/cas/builtin_functions.hpp`: aggiungi `JacobiP` BuiltinOp.
- `src/symbolic/simplify_special_fn.cpp`: rule JacobiP(n,α,β,x) per n integer ≥ 0.

**Test**: P_0/P_1/P_2/P_3 + orthogonality cross-check.

### STEP 39 — L2-16 Cambio variabile auto
**Algoritmo**:
- Estendi `integrate_once` per detect pattern `f(g(x)) · g'(x)`:
  1. Per ogni sub-expression `g(x)` con symbolic substructure, calcola `g'(x) = diff(g, x)`.
  2. Test se integrand = `outer(g(x)) · g'(x)`: split factors, check structural match.
  3. Se match: sostituisci `g(x) = u`, integra `outer(u) du`, back-sub `u = g(x)`.

**Files**: `src/calculus/integrate_substitution.cpp` extend (~100 LOC).

**Test**: `∫ 2x·exp(x²) dx = exp(x²)`, `∫ sin(x)·cos(x) dx = sin²(x)/2`, etc.

### STEP 40 — L3-14 Hensel multivariate MVP
**Algoritmo (Wang 1976)**:
1. Univariate factorization mod p^k esistente.
2. Multivariate: scegli punto valutazione (x_2,...,x_n)=(α_2,...,α_n).
3. Specializza poly → univariato in x_1.
4. Fattorizza univariato.
5. Hensel lift su Z_{p^k}[x_1, x_2, ..., x_n] iterativo per dimensione.
6. Reconstruct factors via leading coefficient correction.

**Files**: `src/algebra/polynomial_hensel.cpp` extend o nuovo `polynomial_hensel_multivariate.cpp` (~300 LOC).

**Test**: 3 polinomi bi-variati noti (`(x+y+1)(x-y+1)`, etc.).

### STEP 41 — L3-17 QR via Gram-Schmidt
**Algoritmo (Gram-Schmidt symbolic)**:
- Q = orthonormal basis di colonne.
- R = upper triangular (Q^T · A).
- Q[i] = (a[i] - Σ proj_{Q[j]}(a[i])) / ||...||.

**Note**: norma sqrt(sum sq) produce sqrt symbolic — Q simbolico ma complesso. Limitarsi a 2×2/3×3 MVP.

**Files**: `src/linalg/matrix_qr.cpp` (~80 LOC).

**Test**: A 2×2 orthogonal verify Q·R = A, Q·Q^T = I.

### STEP 42 — L3-09 FGLM MVP
**Algoritmo (Faugère-Gianni-Lazard-Mora 1993)**:
- Input: Groebner basis G_grevlex per ideal zero-dim I.
- Output: Groebner basis G_lex.
- Approach: linear-algebra over Q[x_1,...,x_n]/I, multiply by x_i, build matrix, find linear dependencies.

**Files**: `src/algebra/polynomial_groebner_fglm.cpp` (~250 LOC).

**Test**: cyclic-3 grevlex → lex.

### STEP 43 — L3-05 Zeilberger Gosper
**Algoritmo (Gosper 1978 + Zeilberger 1990)**:
- Gosper: indefinite hypergeometric summation.
- Zeilberger: extended Gosper for definite sums.

**Files**: `src/calculus/zeilberger.cpp` (~400 LOC research).

### STEP 44 — L2-18 Sparse interpolation Ben-Or-Tiwari
**Algoritmo (Ben-Or-Tiwari 1988)**:
- Input: black-box f(x) presumed sparse t terms.
- Eval f(1), f(ω), ..., f(ω^{2t-1}) where ω = primitive root.
- Berlekamp-Massey su sequenza → recover exponents.
- Vandermonde solve → coefficients.

**Files**: `src/algebra/sparse_interpolation.cpp` (~200 LOC).

### STEP 45 — L3-10 Frobenius advanced
**Algoritmo**: extend `solve_ode_frobenius_at_zero` per gestire roots-differ-by-integer case (logarithmic term).

**Files**: `src/calculus/ode_solver_frobenius.cpp` extend.

### STEP 46 — L3-02 CAD MVP
**Algoritmo (Collins 1975)**: Cylindrical Algebraic Decomposition. Massive scope — defer single-var partial.

### STEP 47-50 — Cleanup debiti

#### STEP 47: PFD inverse Laplace bridge (D2, T1)
Verifica `partial_fractions` output struttura ed adatta inverse_laplace_transform parser.

#### STEP 48: L2-21 strict mode wiring (A1)
In `simplify_arithmetic.cpp` e `simplify_exp_log.cpp`: prima di applicare `exp(ln(x)) → x` o `ln(x·y) → ln(x)+ln(y)`, checkare `ctx.strict_branch_cuts()`. Se true, refuse senza positivity esplicita.

#### STEP 49: L3-13 Interval sin/cos critical points (A2)
Sostituire conservative `[-1, 1]` con detection di π/2 + 2kπ in [lo, hi]. Modular arithmetic via MPFR floor.

#### STEP 50: L3-08 Quantity gaps (A3, A4)
- Quantity^non-int: error onesto se exp non integer.
- Quantity·non-Quantity: già handled by L3-08 STEP 14 group_by_dim — verify.

---

## D. Vincoli CLAUDE.md mantenuti

- L0→L1→L2→L3: ✓ Fase 0 chiusa, lavoro su L2/L3 autorizzato.
- Anti-hardcode: ✓ ogni step include anti-hardcode test.
- BigInt-only nei core: ✓ BigFloat solo in src/numeric/.
- Structural sharing: ✓ return identity quando children invariati.
- 500-line limit: ✓ nuovi file split se necessario.
- -Werror: ✓ verifica build pulita ogni step.
- Result<T> monadico: ✓ no throw/catch.
- Anti-furbizia: ✓ Unimplemented onesto vs silent wrong.

---

## E. Esecuzione

**Modalità**: micro-step, ogni step → 1 commit dedicato + push origin/main.

**Verification per ogni step**:
1. Build `-Werror` pulita.
2. Regression broad: Acid + AcidComplex + suite specifica del step verde.
3. Test nuovi anti-hardcode verdi (incl. negative cases).
4. CAS_TASKS.md sync con stato + evidenza commit.

**Prossimo step suggerito**: STEP 35 (Permuted LU) — basso effort, alto ROI, abilita matrici singolari.
