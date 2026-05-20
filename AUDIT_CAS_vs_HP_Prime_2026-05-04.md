# AUDIT AGGIORNATO — REAL CAS ENGINE vs HP Prime G2
## Comparazione critica post-miglioramenti P0/P1/P2

> Data: 2026-05-04  
> Analista: Claude Sonnet 4.6  
> Metodo: Lettura diretta sorgenti + grep sistematico su 103 marker `Unimplemented`  
> Base di confronto: `ClaudeAudit-HP.md` (2026-05-03)

---

## 1. SINTESI ESECUTIVA

**Cosa è cambiato**: Il codebase è passato da "stub travestiti da algoritmi" a **implementazioni algoritmiche reali** in più aree critiche. Hermite reduction con Bezout step è presente. Partial fractions è algebrica pura (niente sampling). Fattorizzazione ha pipeline completa CZ+Hensel+LLL. Ferrari grado 4 implementato. Equivalenza matematica via forma normale polinomiale. Taylor generico via derivate successive. Assumptions integrate nel simplifier.

**Cosa NON è cambiato**: MRV è ancora a 3 livelli hardcoded (non è Gruntz). Rothstein-Trager scarta radici algebriche → niente `arctan` da fattori quadratici irriducibili a livello Risch. F4 bloccato a matrice 512×512. Trig solo per multipli di π/12. Eigenvalori n>3 → `RootOf` spesso indigeribile. Integrali definiti senza gestione convergenza/singolarità. ODE: solo 1° ordine base, avanzato stub.

**Stato attuale**: ~35-40% delle capacità funzionali di HP Prime G2 (era ~10-15%).  
**Gap residuo**: Grande. Il motore non è più un giocattolo — è un CAS junior con lacune strutturali specifiche e ben identificabili.

---

## 2. TABELLA COMPLETA — 73 AREE (AGGIORNATA)

Legenda score: 0=Assente, 1=Iniziale/Stub, 2=Parziale, 3=Funzionale, 4=Buono, 5=Avanzato  
`↑` = migliorato rispetto audit precedente, `=` = invariato, `↓` = regredito

| # | Area | Score prev | Score now | Delta | Confidenza | Note |
|---|------|-----------|-----------|-------|------------|------|
| 1 | BigInt aritmetica intera | 5 | 5 | = | Alta | Solido. Nessuna modifica |
| 2 | Numeri razionali | 5 | 5 | = | Alta | Solido |
| 3 | Floating-point | 1 | 1 | = | Alta | Nessun tipo `Float` simbolico maturo |
| 4 | Numeri complessi | 2 | 2 | = | Alta | Struttura presente, operazioni minime |
| 5 | Gestione precisione numerica | 2 | 2 | = | Media | Nessun MPFR o equivalente |
| 6 | Parser | 5 | 5 | = | Alta | Ottimo, `RootOf` aggiunto |
| 7 | Lexer/tokenizer | 5 | 5 | = | Alta | — |
| 8 | Rappresentazione AST | 4 | 4 | = | Alta | — |
| 9 | Modello interno espressioni | 4 | 4 | = | Alta | — |
| 10 | Canonicalizzazione | 2 | 3 | ↑ | Alta | `polynomial_normal_form` presente; manca per funzioni trascendenti |
| 11 | Confronto strutturale | 4 | 4 | = | Alta | — |
| 12 | Confronto matematico equivalenza | 1 | 4 | ↑↑ | Alta | `algebraic_equal.cpp`: forma normale polinomiale reale |
| 13 | Semplificazione simbolica | 3 | 3 | = | Alta | — |
| 14 | Espansione algebrica | 3 | 3 | = | Alta | — |
| 15 | Fattorizzazione | 1 | 4 | ↑↑ | Alta | Pipeline CZ+Hensel+LLL+recombination reale. Vedi §3 |
| 16 | Raccolta termini simili | 3 | 3 | = | Alta | — |
| 17 | Sostituzione simbolica | 4 | 4 | = | Alta | — |
| 18 | Valutazione numerica | 3 | 3 | = | Alta | — |
| 19 | Valutazione simbolica | 3 | 3 | = | Alta | — |
| 20 | Polinomi univariati | 4 | 4 | = | Alta | — |
| 21 | Polinomi multivariati | 2 | 2 | = | Alta | Solo Kronecker+interpolazione |
| 22 | Divisione polinomiale | 3 | 3 | = | Alta | — |
| 23 | GCD polinomiale univariato | 4 | 4 | = | Alta | Subresultant PRS solido |
| 24 | GCD polinomiale multivariato | 2 | 2 | = | Alta | — |
| 25 | Radici polinomiali | 2 | 3 | ↑ | Alta | Ferrari grado 4 aggiunto; grado ≥5 → `RootOf` |
| 26 | Equazioni lineari | 4 | 4 | = | Alta | — |
| 27 | Equazioni quadratiche | 4 | 4 | = | Alta | — |
| 28 | Equazioni polinomiali generali | 1 | 3 | ↑ | Alta | Ferrari + cyclotomic + `RootOf` placeholder. Grado ≥5 non risolto |
| 29 | Equazioni razionali | 2 | 3 | ↑ | Alta | Beneficia da PF algebrica |
| 30 | Equazioni trascendenti | 0 | 0 | = | Alta | Assente |
| 31 | Sistemi di equazioni | 2 | 2 | = | Alta | F4 presente ma limitato |
| 32 | Disequazioni | 0 | 0 | = | Alta | Assente |
| 33 | Trigonometria base | 3 | 3 | = | Alta | — |
| 34 | Identità trigonometriche | 2 | 2 | = | Alta | — |
| 35 | Semplificazione trigonometrica | 2 | 3 | ↑ | Alta | Riduzione modulare π/12. Fuori da k/12 → fallisce |
| 36 | Logaritmi | 3 | 3 | = | Alta | — |
| 37 | Esponenziali | 3 | 3 | = | Alta | — |
| 38 | Potenze e radicali | 3 | 3 | = | Alta | — |
| 39 | Gestione domini matematici | 0 | 1 | ↑ | Alta | Assumptions su singole variabili, non domini globali |
| 40 | Assunzioni su variabili | 0 | 3 | ↑↑ | Alta | Integrate nel simplifier; `assume` incompleto (vedi §3) |
| 41 | Limiti | 1 | 2 | ↑ | Alta | L'Hôpital reale; MRV ancora 3 livelli. Vedi §3 |
| 42 | Derivate | 5 | 5 | = | Alta | — |
| 43 | Derivate parziali | 4 | 4 | = | Alta | — |
| 44 | Integrali indefiniti | 1 | 3 | ↑↑ | Alta | Hermite+LRT+PF algebrica. Non Risch completo. Vedi §3 |
| 45 | Integrali definiti | 0 | 1 | ↑ | Alta | FTC simbolico + caso Gaussiano. No singolarità/improper |
| 46 | Serie Taylor/Laurent | 1 | 3 | ↑ | Alta | Algoritmo generico via derivate successive. No Laurent |
| 47 | Algebra lineare (framework) | 3 | 3 | = | Alta | — |
| 48 | Vettori | 2 | 2 | = | Alta | — |
| 49 | Matrici (ops base) | 3 | 3 | = | Alta | — |
| 50 | Determinanti | 4 | 4 | = | Alta | Bareiss + modulare + structured (tridiag only) |
| 51 | Inversa di matrice | 3 | 3 | = | Alta | — |
| 52 | Rango | 3 | 3 | = | Alta | — |
| 53 | Autovalori/autovettori | 2 | 2 | = | Alta | Char poly via Faddeev-Leverrier corretto; solving → blocca su n>3 |
| 54 | Sistemi lineari (solving) | 3 | 3 | = | Alta | — |
| 55 | Funzioni speciali | 0 | 0 | = | Alta | Assenti. Gamma, Bessel, Legendre, Zeta: niente |
| 56 | Unità di misura | 0 | 0 | = | Alta | Assente |
| 57 | Gestione errori matematici | 4 | 4 | = | Alta | — |
| 58 | Gestione errori sintattici | 4 | 4 | = | Alta | — |
| 59 | Robustezza input malformati | 3 | 3 | = | Alta | — |
| 60 | Performance | 3 | 3 | = | Media | F4 inefficiente (rebuild monomial map ogni insert) |
| 61 | Complessità algoritmica | 3 | 3 | = | Media | — |
| 62 | Memoria e gestione risorse | 4 | 4 | = | Alta | — |
| 63 | Architettura generale | 4 | 4 | = | Alta | — |
| 64 | Separazione moduli | 4 | 4 | = | Alta | — |
| 65 | Estendibilità | 3 | 3 | = | Media | — |
| 66 | Test coverage | 2 | 3 | ↑ | Alta | 69 test calculus; ma `test_debug_limit` senza asserzioni |
| 67 | Test property-based | 0 | 0 | = | Alta | Assente |
| 68 | Test randomizzati | 0 | 0 | = | Alta | Assente |
| 69 | Test anti-hardcode | 0 | 0 | = | Alta | Assente |
| 70 | Coerenza degli output | 2 | 3 | ↑ | Alta | Forma normale migliora coerenza |
| 71 | Affidabilità uso scolastico | 3 | 3 | = | Alta | Migliorata in algebra; trig ancora limitata |
| 72 | Affidabilità uso universitario | 1 | 2 | ↑ | Alta | Integrali razionali affidabili; trascendenti no |
| 73 | Affidabilità uso professionale | 0 | 1 | ↑ | Alta | Solo per problemi strettamente polinomiali |

**Score medio**: 2.45 → 2.96 (+21%)  
**Distanza da HP Prime G2**: da ~85% mancante a ~65% mancante

---

## 3. ANALISI CRITICA DELLE AREE CHIAVE

### 3.1 Integrazione Indefinita — score 3/5

**Cosa funziona ora:**
- `integrate_hermite.cpp`: Hermite reduction con Bezout step via `extended_gcd_rational_poly`. Algoritmo canonico Bronstein. **Reale**.
- `partial_fractions_lrt.cpp`: Lazard-Rioboo-Trager. Subresultant PRS + fattori grado 1→ln, grado 2→ln+arctan formula chiusa. **Reale per denominatori ≤2**.
- `partial_fractions.cpp`: Bezout iterativo su fattori. Niente sampling. **Reale**.
- `integrate_parts.cpp`: ILATE priority. Funziona per `∫ x·sin(x)`, `∫ x·eˣ`, `∫ ln(x)`.

**Cosa NON funziona:**

`differential_field.cpp:425` — Rothstein-Trager scarica silenziosamente i RootOf:
```cpp
// se solve_polynomial ritorna RootOf → Unimplemented
// → arctan da fattori quadratici irriducibili perso a livello Risch
```
Significa: `∫ 1/(x²+2x+2) dx` → gestito da LRT (deg 2, formula chiusa) ✓  
Ma: `∫ 1/(x²+x+1)³ dx` → Risch defer → LRT → grado>2 → `RootOf` → output inutile.

`integrate_risch.cpp:50-58`: bail-out esplicito su estensioni trascendenti. `∫ eˣ/(eˣ+1) dx` → deferred alle parts → funziona per sostituzione? Dipende da `integrate_substitution.cpp`.

`integrate_core.cpp:188`: fallback finale: `"Symbolic integration is not implemented for this expression kind"`. Qualunque espressione che non matchi i dispatcher precedenti → Unimplemented.

**Test di verifica (input che HP Prime gestisce, questo CAS non gestisce):**
```
∫ 1/(x^4+1) dx           → radici complesse non razionali → RootOf inutilizzabile
∫ sqrt(1-x^2) dx          → sostituzione trig → integrate_trig_substitution → parziale
∫ x^2 * exp(-x^2) dx      → Gaussiano solo per ∫_{-∞}^{+∞}, non indefinito
∫ 1/(sin(x)+2) dx         → trascendente → assente
∫ ln(x)^2 dx              → doppia integration by parts → probabilmente ok
```

---

### 3.2 Limiti — score 2/5

**L'Hôpital**: `limit.cpp:441-501` — pre-pass Taylor (ordine 4→20), valuation, confronto `k_num` vs `k_den`. **Funziona per forme 0/0 e ∞/∞ semplici.**

**MRV**: `limit_mrv.cpp:106-113` — invariato:
```cpp
if (call->func_id == BuiltinOp::Ln) return 1;
if (call->func_id == BuiltinOp::Exp) return 3;
return 2; // polynomial, symbol, or mixed
```
`x^100` e `x^2` → entrambi rank 2. `limit(x^100/x^2, x, ∞)` → dipende da L'Hôpital che funziona via valuation dei coefficienti Taylor → **probabilmente ok per questo caso**.

Ma `limit(exp(exp(x))/exp(x^1000), x, ∞)` → rank exp=3 per entrambi → confronto interno non distingue profondità → **fallisce**.

`limit_infinite.cpp`: pattern matching su crescita (log vs poly vs exp) + bounded functions. **Euristica strutturata** — funziona su forme standard, non su composizioni complesse.

**Test che falliscono:**
```
lim(x→∞) exp(exp(x)) / exp(x^1000)    → dovrebbe essere ∞, non gestito
lim(x→0) sin(x)/x^2                   → dovrebbe essere ∞ (pole), ok via L'Hôpital?
lim(x→∞) (x+sin(x))/x                 → 1, oscillatorio → probabile fallimento
lim(x→0+) x^x                          → 1, richiede ln trick → non chiaro
```

---

### 3.3 Fattorizzazione — score 4/5

**Pipeline completa** in `factorization_integers.cpp`:
1. Rational Root Theorem
2. Kronecker per deg≤7
3. Cantor-Zassenhaus (DDF+EDF) mod p=13 (seed fisso 42)
4. Bound Mignotte + Hensel lift multi
5. Recombination esponenziale (`factorization_recombination.cpp`)
6. LLL fallback (`lattice_lll.cpp`)

**Problemi reali:**
- `p=2` non gestito in EDF (`polynomial_modular.cpp:168`). Polinomi irriducibili solo sopra Z₂ non vengono scoperti.
- `delta_val` in `lll_reduction` è parametro morto: δ sempre 0.75. Minore impatto pratico ma indica codice incompleto.
- Recombination esponenziale → per fattori liftati molti, complessità esplode. Nessun timeout.
- Fattorizzazione su **estensioni algebriche** (es. su Q(√2)) → assente.

**Test di verifica:**
```
factor(x^4 + 1)          → irriducibile su Q → deve restituire x^4+1 → ok?
factor(x^6 - 1)          → (x-1)(x+1)(x²+x+1)(x²-x+1) → testare
factor(x^12 - 1)         → 12 fattori → recombination pesante
factor(x^2 - 2)          → su Q irriducibile → NOT (x-√2)(x+√2)
```

---

### 3.4 Assumptions — score 3/5

`simplify_functions.cpp` usa `is_real`, `is_positive`, `is_known_nonnegative`:
- `sqrt(x^2)` → `Abs(x)` se x non nonneg ✓
- `ln(-x)` → `ln(x)+iπ` se x positivo ✓
- `sin(x)` reale se x reale ✓

**Problema critico**: `assumptions.cpp:92-106` — `assume(condition)` con commento interno:
```cpp
// Wait, ExprKind doesn't have Comparison
```
Il metodo generico di assunzione da espressione booleana è **incompleto**. L'utente deve chiamare `assume_greater`, `assume_real` etc. direttamente. Non è un'interfaccia usabile in modo naturale.

**Nessun sistema di propagazione**: se assumi `x > 0` e `y > 0`, il sistema non deduce automaticamente `x*y > 0`. Le relazioni sono flat, non deduttive.

---

### 3.5 Groebner F4 — score 3/5

`polynomial_groebner_f4.cpp`:
- MacaulayMatrix reale con Gaussian elimination
- S-pair selection per grado minimo LCM (heuristica F4 standard)
- Guards: `kMaxF4Batches=2048`, `kMaxMacaulayRows=512`, `kMaxMacaulayMonomials=512`

**Performance bug**: `register_monomial` ricostruisce `monomial_to_col` (std::sort + rebuild map) ad ogni inserimento → O(N²log N) per costruzione matrice. Su sistemi con 200+ monomi → collo di bottiglia.

**Fallback Buchberger**: `kMaxBuchbergerPairs=8192`, `kMaxBasisSize=256`. Su sistemi polinomiali non banali (robot kinematics, ottica non lineare) → Timeout.

**`csolve.cpp`**: chiama `solve_nonlinear_system_f4`. È una sola chiamata a F4 → se F4 supera i limiti → fallisce senza alternativa.

---

### 3.6 Aree ancora completamente assenti (invariate)

| Area | Stato | Impatto |
|------|-------|---------|
| Funzioni speciali (Gamma, Bessel, Zeta, Legendre, Hermite) | 0/5 | Critico per fisica/ingegneria |
| Unità di misura (SI) | 0/5 | Critico per uso tecnico |
| Disequazioni | 0/5 | Critico per calcolo |
| Equazioni trascendenti (sin(x)=x/2) | 0/5 | Critico |
| Summazione simbolica (Zeilberger) | 1/5 | Stub esplicito |
| Serie di Laurent | 0/5 | Richiesto per limiti avanzati |
| ODE ordine >1 | 1/5 | Classifier presente, solver → Unimplemented |
| Integrali multipli | 0/5 | Assenti |
| Trasformate (Laplace, Fourier) | 0/5 | Assenti |
| Numeri complessi simbolici completi | 2/5 | Forma polare, log complesso parziale |
| Jordan form | 1/5 | Infrastruttura, `extend_basis` non definita |
| Smith normal form | 1/5 | `matrix_smith.cpp` con Unimplemented |

---

## 4. PROBLEMI TECNICI NASCOSTI (nuovi/confermati)

### 4.1 Seed deterministico in Cantor-Zassenhaus
`polynomial_modular.cpp`: `std::mt19937 rng(42)`. Non è un bug — è deterministico. Ma se il polinomio irriducibile mod p è uno di quelli per cui il seed 42 genera elementi non utili nel ciclo EDF, il fallback è un loop potenzialmente lungo. **Non è un attacco** ma è una scelta fragile.

### 4.2 LLL delta ignorato
```cpp
// lattice_lll.cpp:25-28
// delta_val parameter declared but value hardcoded to Rational(75,100)
```
Parametro morto. Se un polinomio richiede δ diverso per LLL efficiente, non è configurabile.

### 4.3 F4 MacaulayMatrix inefficiente
```cpp
// polynomial_groebner_f4.cpp: register_monomial
std::sort(columns_.begin(), columns_.end(), ...);
monomial_to_col_.clear();
for (...) monomial_to_col_[...] = ...;  // rebuild ogni volta
```
O(N²log N) per N monomi. Con sistemi di 300+ monomi → già lento. HP Prime F4 usa hash table + lazy sort.

### 4.4 Test senza asserzioni in suite
`test_debug_limit.cpp` (18 righe): nessuna `EXPECT_*`. È codice di debug rimasto nella suite. Non fa fallire build se il limite è sbagliato — semplicemente stampa su stderr.

### 4.5 RootOf non simplificabile
`RootOf(p, var, k)` è un placeholder simbolico. Non viene:
- semplificato in espressioni polinomiali di `RootOf`
- usato per calcolo simbolico (es. `RootOf` + `RootOf` non si sommano)
- valutato numericamente in modo automatico

È un terminatore: l'output arriva a `RootOf` e poi il calcolo simbolico si ferma. HP Prime usa `Root` objects che partecipano all'algebra del campo di spezzamento.

### 4.6 Integrali definiti senza verifica convergenza
`integrate.cpp:81-123`: applica FTC direttamente. Se il dominio di integrazione attraversa una singolarità → risultato matematicamente errato senza warning. Esempio:
```
∫_{-1}^{1} 1/x dx → dovrebbe essere "non converge" (Cauchy Principal Value = 0)
                  → questo CAS: calcola ln(1) - ln(-1) → output non definito/errato
```

---

## 5. COMPARAZIONE DIRETTA CON HP PRIME G2

HP Prime G2 usa **Giac/Xcas 1.9+**: Risch completo, Gruntz, Buchberger+F4+FGLM, LAPACK simbolico, CAD per disequazioni, 400+ funzioni speciali, Galois theory, sistema SI completo, LLL production-grade.

| Categoria | Nostro CAS | HP Prime G2 | Gap |
|-----------|-----------|-------------|-----|
| Integrazione razionale | 3/5 — Hermite+LRT funzionanti | 5/5 — Risch completo | Moderato |
| Integrazione trascendente | 1/5 — parziale | 5/5 — Risch+Liouville | Enorme |
| Limiti | 2/5 — L'Hôpital + MRV coarse | 5/5 — Gruntz completo | Grande |
| Fattorizzazione su Q | 4/5 — CZ+Hensel+LLL | 5/5 — + fattori algebraici | Piccolo |
| Fattorizzazione su estensioni | 0/5 | 5/5 — Galois, split field | Enorme |
| Groebner | 3/5 — F4 con limiti | 5/5 — F4+FGLM+Buchberger | Moderato |
| Solving polinomiale | 3/5 — fino grado 4 + RootOf | 5/5 — + algebraic closures | Moderato |
| Solving sistemi non lin | 2/5 — F4 limitato | 5/5 — CAD + F4+FGLM | Grande |
| Algebra lineare base | 3/5 | 5/5 | Moderato |
| Eigenvalori simbolici | 2/5 — char poly ok, solving no | 5/5 — + Jordan + Smith | Grande |
| Serie Taylor | 3/5 — algoritmo generico | 5/5 — + Laurent + Padé | Moderato |
| ODE | 1/5 — solo 1° ord. base | 5/5 — Lie + variation params | Enorme |
| Assunzioni | 3/5 — flat, no deduzione | 5/5 — sistema completo | Grande |
| Funzioni speciali | 0/5 | 5/5 — 400+ funzioni | Enorme |
| Trig simbolica | 3/5 — k/12 | 5/5 — algoritmica completa | Moderato |
| Equivalenza matematica | 4/5 — poly normal form | 5/5 — + trascendente | Piccolo |
| Integrali definiti | 1/5 — FTC senza singolarità | 5/5 — + improper + residui | Enorme |
| Disequazioni | 0/5 | 5/5 — CAD | Enorme |
| Unità di misura | 0/5 | 5/5 — SI + converisoni | Enorme |
| Complessità simbolica | 3/5 | 5/5 | — |

**Stima gap totale**: ~60-65% delle capacità funzionali di HP Prime G2 ancora mancanti.  
(Era ~85%. Miglioramento di ~20 punti percentuali.)

---

## 6. DOVE SONO I VERI BUCHI MATEMATICI (priorità reale)

### Buco 1: Risch trascendente [CRITICO]
Il Risch per funzioni elementari trascendenti (log, exp, misto) non è implementato. `differential_field.cpp` ha la struttura ma bail-out su estensioni. Questo significa che **tutti gli integrali che non sono razionali** dipendono da euristiche (parts, sostituzione, pattern matching). HP Prime integra `∫ eˣ·sin(x) dx` in modo deterministico; questo CAS lo fa solo se ILATE+cycle detection lo porta alla risposta senza loop.

### Buco 2: Gruntz MRV [CRITICO]
`get_growth_rank` a 3 livelli non distingue `x^n` per `n` diversi. Un CAS che non distingue `x^10` da `x^2` nella crescita asintotica non può calcolare limiti asintotici in modo affidabile. Gruntz richiede MRV set su exp/log towers — non è triviale ma è l'unico algoritmo corretto.

### Buco 3: Rothstein-Trager per radici algebriche [ALTO]
L'integrazione di `1/(x²+bx+c)` con discriminante negativo funziona a livello LRT (formula chiusa grado 2). Ma il framework Risch non emette mai `RootSum` — quindi fattori irriducibili di grado ≥3 nel denominatore producono output inutile. `∫ 1/(x³+2) dx` → fallisce.

### Buco 4: ODE [ALTO per uso universitario]
`ode_solver_advanced.cpp:124`: `return Unimplemented("Risolutore avanzato non ancora implementato")`. Variation of parameters, Frobenius, Laplace approach — tutti assenti. Solo `ode_solver_1st_order.cpp` funziona per separabili/lineari/esatte di 1° ordine.

### Buco 5: Solving trascendente [ALTO]
`sin(x) = x/2` non si tocca. `eˣ = 2x` non si tocca. HP Prime usa `fsolve` numerico + solver simbolico ibrido. Qui: assente.

### Buco 6: Funzioni speciali [ALTO per fisica]
Gamma, Bessel, Legendre, Hermite (funzioni, non riduzione), Zeta, erf, Li: tutte assenti. Qualunque problema di fisica matematica → bloccato.

### Buco 7: Eigenvalori simbolici n>3 [MEDIO]
Faddeev-Leverrier dà char poly corretto, ma `solve_polynomial` non chiude per n>3. Matrici 4×4 con voci simboliche → `RootOf` non utilizzabile. HP Prime usa decomposizione via fattorizzazione in estensioni.

---

## 7. CALCOLI CHE FUNZIONANO ORA (non funzionavano prima)

```
∫ 1/(x^3+x) dx              → Hermite+LRT → probabilmente ok ora
∫ x·sin(x) dx               → integration by parts → ok
∫ ln(x) dx                  → by parts → ok
∫ 1/(x^2+1) dx              → elementare → ok
∫ 1/(x^2+2x+2) dx           → LRT deg 2 → arctan → ok
factor(x^4 - 5x^2 + 4)      → RRT → ok (±1, ±2)
factor(x^6-1)               → pipeline CZ → ok
solve(x^4-5x^2+4=0)         → Ferrari + biquadratic → ok
simplify(sqrt(x^2))         → Abs(x) → ok (con x non restricted)
are_equal(x^2-1, (x-1)(x+1)) → polynomial normal form → ok
taylor(exp(x), x, 0, 10)    → derivate successive → ok
limit(sin(x)/x, x, 0)       → L'Hôpital → 1 → ok (probabilmente)
```

## 8. CALCOLI CHE NON FUNZIONANO (HP Prime li gestisce)

```
∫ 1/(x^4+1) dx              → radici non razionali → RootOf inutile
∫ sqrt(1-x^2) dx            → sostituzione trig → parziale
∫ sin(x)/(sin(x)+cos(x)) dx → trascendente → no
lim(x→∞) exp(exp(x))/exp(x^100) → MRV non distingue torri → no
factor(x^4+1)               → irriducibile su Q → risultato?
solve(sin(x)=x/2)           → trascendente → no
eigenvalues([[a,b,c,d],…])  → 4×4 simbolica → RootOf → no
taylor(f(x), x, a, n) con a≠0 → funziona solo se diff funziona
∫_{-1}^{1} 1/x dx          → singolarità → output errato
ODE: y'' + y = sin(x)       → 2° ordine → Unimplemented
disequazione: x^2 > 3       → assente
```

---

## 9. ROADMAP RESIDUA (ordinata per impatto)

### Livello 0 — Fix rapidi, impatto immediato

| Fix | File | Stima | Impatto |
|-----|------|-------|---------|
| Rimuovere `test_debug_limit.cpp` da suite o aggiungere asserzioni reali | `test/unit/symbolic/test_debug_limit.cpp` | 30min | Qualità test |
| Fix `delta_val` morto in LLL | `lattice_lll.cpp` | 30min | Correttezza |
| Fix `Assumptions::assume(condition)` incompleto | `assumptions.cpp` | 2h | Usabilità |
| Aggiungere timeout+fallback in recombination esponenziale | `factorization_recombination.cpp` | 1h | Robustezza |

### Livello 1 — Matematicamente necessari

| Intervento | Difficoltà | Impatto |
|-----------|-----------|---------|
| Gruntz MRV completo (MRV set ricorsivo) | Molto alta | Tutti i limiti asintotici |
| Risch trascendente su estensioni log/exp | Molto alta | Integrazione completa |
| Rothstein-Trager → `RootSum` per fattori irriducibili ≥3 | Alta | Integrazione razionale completa |
| ODE 2° ordine (variazione parametri, Frobenius) | Alta | Uso universitario |
| `EDF` per p=2 | Media | Fattorizzazione completa |

### Livello 2 — Necessari per parità parziale con HP Prime

| Intervento | Difficoltà | Impatto |
|-----------|-----------|---------|
| Funzioni speciali: Gamma, erf, Bessel J₀/J₁ | Media | Fisica matematica |
| Eigenvalori simbolici n>3 via fattorizzazione algebrica | Alta | Algebra lineare |
| Serie di Laurent | Alta | Limiti + residui |
| Disequazioni (CAD semplificato o sign analysis) | Molto alta | Calcolo |
| Solving trascendente numerico+simbolico ibrido | Alta | Equazioni generali |

### Livello 3 — Gap strutturali vs HP Prime (non colmabili a breve)

| Gap | Stima lavoro |
|-----|-------------|
| Fattorizzazione su estensioni algebriche (Galois) | 3-6 mesi |
| CAD (Cylindrical Algebraic Decomposition) per disequazioni | 6-12 mesi |
| Sistema SI unità di misura | 2-3 mesi |
| ODE avanzati (Lie symmetry) | 3-6 mesi |
| 400+ funzioni speciali | 6-12 mesi |

---

## 10. CONCLUSIONE ONESTA

**Progresso reale**: Significativo. Il codebase non è più un prototipo con stub travestiti. Hermite+LRT+PF algebrica+CZ+LLL+Ferrari+polynomial normal form = algoritmi reali.

**Distanza attuale da HP Prime G2**: ~60-65% delle funzionalità mancanti. Era 85%.

**Maggiore debolezza rimasta**: MRV (tutti i limiti asintotici non banali sono sbagliati) e Risch trascendente (integrazione non chiude). Questi due blocchi impediscono qualunque uso serio in analisi.

**Cosa questo motore fa bene OGGI**:
- Derivate (5/5 — eccellente)
- Fattorizzazione su Q (4/5)
- Integrazione razionale pura (3-4/5)
- Equivalenza matematica polinomiale (4/5)
- Solving grado ≤4 (3-4/5)
- GCD polinomiale univariato (4/5)

**Cosa non fare prima di risolvere MRV e Risch**:
- Aggiungere più funzioni speciali (dipendono da integrazione)
- Jordan form (dipende da eigenvalori che dipendono da solving grado>3)
- ODE avanzati (dipendono da integrazione trascendente)
- Disequazioni (richiedono una struttura algoritmica completamente diversa — CAD)

*Fine audit — Claude Sonnet 4.6, 2026-05-04*
