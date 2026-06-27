> # ⚠️ SUPERSEDED — file storico (non aggiornare)
> Tutti i task di questo file sono stati **consolidati e ri-verificati a codice** in
> **[`TASKLIST_MASTER.md`](TASKLIST_MASTER.md)** (audit 2026-06-26, single source of truth).
> Conservato solo come storico/contesto. Le voci di stato qui dentro sono **obsolete**.

---

# Handoff — CAS Engine session (continua in nuova chat)

Branch: `chore/claude-scaling-tooling`. Tutto buildato/testato sotto **ASan+UBSan**
(dir `build`, Debug). Gate suite: `bash scripts/test_quick.sh` (cap 600s, foreground).
Mai `run_in_background` per la suite intera. Build: `cmake --build build --target cas_foundation_tests -j8`.

## Stato: commit di questa sessione (tutti full-suite-verified, 0 fail)
- `9555373` refactor anti-monolito: split trig linearization da simplify_arithmetic_power.cpp
- `d03ff6d` T-004 cos(7π/16) flaky → co-function reduction (kill Chebyshev T₇ blow-up)
- `89eae9e` T-017 ln(a+bi) esatto su tutta la famiglia 30°/45°/60°, tutti i quadranti
- `4597332` T-025 mathematically_equal Φ_n ciclotomico composito ↔ exp(2πi/n)
- `67be548` test: copertura integrali definiti FTC + gap documentati
- `6c5ae48` feat: ∫1/√(Ax²+Bx+C) generale (completamento del quadrato → arcsin/ln)

Ultima suite completa: **2440 passed / 0 failed / 2 skipped** (pre-esistenti).
Tree pulito tranne `?? TASKLIST_MASTER.md` (planning doc pre-esistente, NON committare).

## Task corrente: T-055 — ∫xᵏ/√(1−x²) per k≥1 (IN CORSO, analisi fatta, codice NON scritto)

### Mappa coperture (probe già eseguito, 2026-06-20)
- `x/√(1-x²)` (k=1) → **OK** (esistente, x/√(a²-x²) in integrate_core)
- `x²/√(1-x²)` (k=2) → **OK** (`integrate_xsq_over_sqrt_quadratic`)
- `x³,x⁴,x⁵/√(1-x²)` (k≥3) → **NO_STRATEGY** ← QUESTO È IL GAP
- `asin(x)`, `x·asin(x)` → OK
- `x²·asin(x)`, `x³·asin(x)` → **NO_STRATEGY** (bloccati dal gap k≥3: la IBP riduce
  ∫xⁿ·asin → xⁿ⁺¹/(n+1)·asin − 1/(n+1)·∫xⁿ⁺¹/√(1−x²); implementare k≥3 li sblocca)

### Matematica (derivata e verificata a mano)
Radicando generale **R = c − x²** (c>0 razionale; A=−1, B=0, C=c). Niente hardcode di c=1.
Formula di riduzione:

    I_k = ∫ xᵏ/√(c−x²) dx = [ (k−1)·c·I_{k-2} − xᵏ⁻¹·√(c−x²) ] / k        (k≥2)

Derivazione: d/dx[xᵏ⁻¹√(c−x²)] = (k−1)c·xᵏ⁻²/√(c−x²) − k·xᵏ/√(c−x²)  ⟹ formula.
Basi (costruirle ESPLICITAMENTE per non dipendere dalla copertura engine):
- I₀ = ∫1/√(c−x²) = arcsin(x/√c)   [già coperto da `integrate_inverse_sqrt_quadratic`, commit 6c5ae48]
- I₁ = ∫x/√(c−x²) = −√(c−x²)

Per radicando generale **c − d·x²** (d>0): I_k = [(k−1)c·I_{k-2} − xᵏ⁻¹√(c−dx²)]/(k·d),
I₀=(1/√d)arcsin(√d·x/√c), I₁=−(1/d)√(c−dx²). Scegliere se generalizzare a d≠1 o
scopare a d=1 (c−x², che è ciò che la catena asin produce). Raccomandato: gestire c−x²
(A=−1,B=0,C>0) pulito con Unimplemented esplicito fuori range (Cat 4 OK, mai silenzioso).

### Punto di innesto preciso
`src/calculus/integrate_product_power.cpp:210-264` — blocco `HC-IBP-VDU` che riconosce
la shape `c · x² · (1/√R)` (accetta sia `Pow(sqrt(R),-1)` sia `Div(1,sqrt(R))`) e chiama
`integrate_xsq_over_sqrt_quadratic`. **Generalizzare**: estrarre la potenza intera k di x
(non solo ==2) e, per k≥3 (o k≥2 unificato), chiamare un nuovo helper ricorsivo
`integrate_monomial_over_sqrt_quadratic(k, R, var)` che applica la riduzione sopra,
bottoming su k=0 (arcsin via `integrate_inverse_sqrt_quadratic`) e k=1 (−√(c−x²) costruito).

Dichiarare il metodo in `src/calculus/integrate_engine.hpp` (classe `Integrator`, ~riga 62,
vicino a `integrate_xsq_over_sqrt_quadratic`). Implementarlo in
`src/calculus/integrate_trig_substitution.cpp` (dove vivono i fratelli quadratic-radical e
dove sta già `integrate_inverse_sqrt_quadratic`). Helper disponibili lì:
`make_sum/make_product/make_binary/make_function/make_rational/make_integer/make_unary`,
`exact_scalar_from_expr(ExprPtr)→optional<Rational>`, `algebra::univariate_coefficients`.

### Verifica (OBBLIGATORIA, indipendente dal simplifier)
Pattern già usato in `test/unit/calculus/test_integrate_ibp_coverage.cpp`:
`numeric::eval(diff(integrate(f)), env)` vs `numeric::eval(f, env)` su più sample point
del dominio reale (|x|<1 per √(1−x²)), tolleranza 1e-7. cas_core linka cas_numeric →
`numeric::eval` disponibile nei test (`#include "cas/numeric.hpp"`, env `unordered_map<string,double>`).
Aggiungere casi: x³,x⁴,x⁵/√(1−x²) e x²·asin(x), x³·asin(x) (questi ultimi sbloccati
automaticamente). NON usare round-trip simbolico (il simplifier non riconcilia sempre √-power forms).

### Regole da rispettare (CLAUDE.md)
- REGOLA ZERO: algoritmo generale, NO pattern-table chiusa, NO hardcode di c=1.
- Anti-monolito: integrate_trig_substitution.cpp è a ~337 LOC, integrate_product_power.cpp
  a 374, integrate_core.cpp ≤500, integrate_power.cpp 143 — c'è spazio ma ricontrolla con
  `bash scripts/check_file_size.sh` (hard block >500, Gate 1).
- Gate pre-commit: `bash scripts/debt_gate.sh` (4 gate). Commit solo dopo suite completa verde.
- Niente `git reset --hard`/`git restore --source`; solo `git stash push` per backup.
- Commit msg termina con:
  `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`
  `Claude-Session: https://claude.ai/code/session_01DXrcfQM4bmS2EMi7UtTVy3`
- Dopo commit: marcare T-055 (TaskUpdate id 55) — ma resta "parziale" finché anche acos
  non è coperto (acos = π/2 − asin, dovrebbe seguire). Aggiornare ledger se c'è entry.

## Pattern di lavoro che ha funzionato (mantienilo)
"Probe-first": scrivi test che FALLISCE sulla feature dichiarata-coperta, trova il gap reale,
fixa alla radice. Il task list + HARDCODE_LEDGER **sovrastimano** il lavoro aperto e marcano
"RISOLTO/CHIUSO" cose con gap residui reali (vedi memoria `tasklist-ledger-stale-2026-06-19`).
Verifica SEMPRE prima di implementare. `python3 scripts/ledger_index.py` status flags inaffidabili.

## Note tecniche apprese
- `1/√(1-x²)` canonicalizza a `Pow(Sqrt(1-x²), -1)`, non Div.
- builtin disponibili: Asin, Acos, Ln, Sqrt. **NO Asinh/Acosh** con regole diff/numeric →
  per A>0 usa forma ln(x+√…); per A<0 usa arcsin. (asinh/acosh sono emessi altrove ma NON
  hanno diff/numeric support — non usarli in codice verificato numericamente.)
- `is_known_positive` ora copre √x; `is_known_negative` copre base^n (odd); simplifier ha
  regola `(−x)^n` intero → ±x^n (commit 89eae9e).
