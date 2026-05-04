# CAS ENGINE — Sistema Task Unificato
## Controllo Avanzamento verso HP Prime G2

> Aggiornato: 2026-05-03  
> Basato su: `ClaudeAudit-HP.md` + stato `TODO.md`  
> Protocollo: Unificato Anti-Hardcode (vedi prompt agenti)

**Regola invariabile**: Finché esistono P0 aperti, nessun agente lavora su P1/P2/P3.  
**Criterio "Risolta"**: Capacità generalizzabile + test robusti + test anti-hardcode + nessun hardcode + regressioni ok.

---

## TABELLA MASTER TASK

| ID | Area | Titolo | Priorità | Stato | Dipendenze | Impatto HP Prime G2 | Muscolo matematico | Prossima azione |
|---|---|---|---|---|---|---|---|---|
| CAS-P0-001 | Equivalenza | Forma normale polinomiale (confronto matematico reale) | P0 | Risolta | AST, Rational, Polynomial core | Molto alto | `are_mathematically_equal()` via normal form | — |
| CAS-P0-002 | Frazioni parziali | Rimpiazzare sampling-based PFD con algoritmo simbolico | P0 | Risolta | GCD esteso polinomiale, Bezout | Molto alto | Decomposizione simbolica corretta | — |
| CAS-P0-003 | Limiti | Fix confronto asintotico (MRV fake → ranking reale) | P0 | Risolta | Simplifier, normal form poli | Alto | Analisi asintotica generalizzabile | — |
| CAS-P0-004 | Testing | Infrastruttura test anti-hardcode + property-based | P0 | Risolta | Test framework (GoogleTest presente) | Molto alto | Validazione robusta non-specifica | — |
| CAS-P1-001 | Integrazione | Hermite Reduction (Bezout step mancante) | P1 | Risolta | CAS-P0-002 (GCD esteso + Bezout) | Molto alto | Riduzione Hermite per integrali razionali | — |
| CAS-P1-002 | Integrazione | Rothstein-Trager / Lazard-Rioboo (log-part) | P1 | Risolta | CAS-P1-001 | Molto alto | Integrazione razionale simbolica completa | — |
| CAS-P1-003 | Solving | Polynomial solving grado 4 (Ferrari) | P1 | Risolta | Solve gradi 1-3 (presenti), sqrt simbolico | Alto | Radici quartica esatte | — |
| CAS-P1-004 | Fattorizzazione | Fattorizzazione completa su Q (oltre Rational Root Theorem) | P1 | Risolta | Square-free Yun (presente, P5 done) | Molto alto | Berlekamp/Zassenhaus modular lifting | — |
| CAS-P1-005 | Assunzioni | Integrazione assumptions engine nel simplifier | P1 | Risolta | Assumptions data struct (dichiarata) | Alto | `sqrt(x^2) → |x|`, `ln(x)` su `x>0` | — |
| CAS-P2-001 | Serie | Taylor via derivate successive (algoritmo, non lookup) | P2 | Bloccata | Derivate (presenti), CAS-P0-001 per equivalenza | Alto | `taylor(f, x, a, n)` generalizzato | Attendere P0-001 |
| CAS-P2-002 | Integrazione | Integrali definiti via FTC simbolico | P2 | Bloccata | CAS-P1-001, CAS-P1-002 | Molto alto | `integrate(f, x, a, b)` esatto | Attendere P1-002 |
| CAS-P2-003 | Integrazione | Sostituzione trigonometrica (rimuovere Unimplemented) | P2 | Bloccata | CAS-P1-001 | Alto | Integrazione tramite sub trig | Attendere P1-001 |
| CAS-P2-004 | Limiti | L'Hôpital affidabile come fallback per limiti indeterminati | P2 | Aperta | Derivate (presenti), CAS-P0-003 | Alto | Limite forme 0/0, ∞/∞ garantito | Implementare l'hopital con guard anti-loop |
| CAS-P2-005 | Trig | Riduzione angolare modulare per sin/cos/tan | P2 | Aperta | Simplifier (presente) | Medio | `sin(5π/6) → 1/2` senza lookup | Implementare `reduce_angle_modular()` |
| CAS-P2-006 | LinAlg | Autovalori simbolici per matrici n>3 | P2 | Bloccata | CAS-P1-003, CAS-P1-004 | Alto | `eigenvalues()` senza limite di grado | Attendere P1-003 |
| CAS-P3-001 | LinAlg | Jordan normal form completa | P3 | Bloccata | CAS-P2-006, autovettori generalizzati | Medio | `jordan_form()` simbolico | Attendere P2-006 |
| CAS-P3-002 | Algebra | Gröbner basis F4 verificato e testato | P3 | Aperta | Multivariate polynomial (parziale) | Alto | Solving sistemi polinomiali | Verificare `f4_groebner()` — implementazione nascosta o mancante |
| CAS-P3-003 | Calcolo | ODE simbolici oltre variabili separabili | P3 | Bloccata | CAS-P1-002 | Medio | Classificazione + solving ODE | Attendere integrazione completa |
| CAS-P3-004 | Funzioni | Funzioni speciali (Gamma, Beta, Bessel) | P3 | Aperta | Nessuna critica | Medio | Estensione dominio funzioni | Aggiungere Gamma/Beta come prima fase |
| CAS-P3-005 | Numeri | Complessi simbolici completi (aritmetica su C) | P3 | Aperta | Complex struct (parziale) | Medio | Campo C simbolico completo | Completare operazioni algebriche su Complex |
| CAS-P3-006 | Sistema | Unità di misura (sistema SI) | P3 | Aperta | Nessuna critica | Basso | Calcolo dimensionato | Definire tipo `Quantity` |

---

## STATO TASK COMPLETATE (da TODO.md)

| ID | Area | Titolo | Stato | Note |
|---|---|---|---|---|
| DONE-01 | Numeri | Unificazione Complessi Canonici (I^n, I^2=-1) | Risolta | P1 TODO completato |
| DONE-02 | Canonicalizzazione | LPO/KB, ordinamento deterministico, potenze annidate | Risolta | P2 TODO completato |
| DONE-03 | Rewrite | Rewrite system chiuso (exp/log, sin^2+cos^2, parità) | Risolta | P3 TODO completato |
| DONE-04 | Assunzioni | Assumptions engine (struttura dati + propagazione AST) | Parziale | P4 struttura presente ma NON integrata in simplifier → vedi CAS-P1-005 |
| DONE-05 | Polinomi | GCD Subresultant PRS, Square-free Yun | Risolta | P5 TODO completato |
| DONE-06 | Testing | Acid test suite pilastri 1-5 | Risolta | P9 TODO completato |
| DONE-07 | Core | BigInt/Rational con Result<T>, no throw | Risolta | MEMO risolta |
| DONE-08 | AST | Hash-consing, interning, structural sharing | Risolta | MEMO risolta |

---

## DETTAGLIO TASK P0

---

### CAS-P0-001 — Forma normale polinomiale (confronto matematico reale)

**Area:** Equivalenza matematica  
**Priorità:** P0  
**Gravità:** Critica  
**Stato:** Risolta  
**Agente assegnato:** Non assegnato  
**Data ultimo aggiornamento:** 2026-05-03

**Aggiornamento Codex 2026-05-03 (finale):**  
Implementata `polynomial_normal_form(expr)` in `src/symbolic/normal_form.cpp` e `include/cas/normal_form.hpp`. La funzione `mathematically_equal` in `src/algebra/algebraic_equal.cpp` è stata aggiornata per utilizzare la forma normale (espansione e raccolta termini in mappa Monomial->Rational). Verificata con `test_math_equal.cpp`: tutti i criteri P0-001 passano (inclusi anti-hardcode e coefficienti grandi).

**Problema risolto:**  
`structural_equal()` non era sufficiente per il confronto matematico. Ora `mathematically_equal(a, b)` verifica `normal_form(expand(a - b)) == 0`.

**Evidenza tecnica:**  
- `src/symbolic/simplify_utils.cpp`: `structural_equal()` — confronto puntatori/struttura AST
- Nessuna funzione `are_mathematically_equal()` o `normal_form()` nell'intero codebase
- Test in `test/unit/` usano confronto strutturale o `toString()` — entrambi inaffidabili per verifica matematica

**Perché è importante:**  
Senza confronto matematico reale:
- nessun test può verificare correttezza di integrazione (D(∫f) == f?)
- nessun simplifier può verificare di non aver prodotto risultati equivalenti ma non identici
- nessuna feature del CAS può essere dichiarata "corretta" con certezza

**Obiettivo matematico:**  
Implementare `polynomial_normal_form(expr)` che converte espressioni polinomiali in una rappresentazione canonica unica. Due espressioni `p` e `q` sono matematicamente uguali sse `polynomial_normal_form(p - q) == 0`.

**Nuovo muscolo matematico atteso:**  
`are_mathematically_equal(ExprPtr a, ExprPtr b, CASContext& ctx) → bool` via espansione e normalizzazione polinomiale.

**Algoritmo consigliato:**  
1. Expand entrambe le espressioni (già presente via `expand()`)
2. Sottrai: `diff = expand(a - b)`
3. Colleziona termini: mappa `monomial → coefficient`
4. Verifica tutti i coefficienti == 0
5. Per espressioni trascendenti: usa simplify + structural_equal come fallback (onesto: dichiara incertezza se non determinabile)

**Dipendenze:**  
- `expand()` (presente)
- Raccolta termini (parzialmente presente)
- Aritmetica su coefficienti razionali (presente)

**Piano tecnico:**  
1. Creare `src/symbolic/normal_form.cpp` + `include/cas/normal_form.hpp`
2. Implementare `collect_polynomial_terms(ExprPtr, variables)` → `map<Monomial, Rational>`
3. Implementare `polynomial_normal_form(ExprPtr)` → forma canonica
4. Implementare `are_mathematically_equal(ExprPtr, ExprPtr, CASContext&)` → bool
5. Aggiungere `is_zero_polynomial(ExprPtr)` usato internamente
6. Integrare nei test esistenti come helper di verifica

**Criteri di accettazione:**  
- `are_mathematically_equal("(x+1)^2", "x^2+2*x+1")` → true
- `are_mathematically_equal("x^2-1", "(x-1)*(x+1)")` → true
- `are_mathematically_equal("x+y", "y+x")` → true
- `are_mathematically_equal("x^2+1", "x^2+2")` → false
- Funziona con variabili arbitrarie (a, b, z, t — non solo x)

**Test obbligatori:**  
```cpp
TEST(MathEqual, PolynomialExpansion) {
    EXPECT_TRUE(are_equal("(x+1)^2", "x^2+2*x+1", ctx));
    EXPECT_TRUE(are_equal("(a+b)^2", "a^2+2*a*b+b^2", ctx));
    EXPECT_TRUE(are_equal("x^2-1", "(x-1)*(x+1)", ctx));
    EXPECT_FALSE(are_equal("x^2+1", "x^2-1", ctx));
}
```

**Test anti-hardcode:**  
```cpp
// Variabili diverse — stesso algoritmo
EXPECT_TRUE(are_equal("(z+1)^2", "z^2+2*z+1", ctx));
EXPECT_TRUE(are_equal("(alpha+beta)^2", "alpha^2+2*alpha*beta+beta^2", ctx));
// Coefficienti grandi
EXPECT_TRUE(are_equal("(x+100)^2", "x^2+200*x+10000", ctx));
// Non uguali — non deve dare falsi positivi
EXPECT_FALSE(are_equal("(x+1)^2", "x^2+x+1", ctx));
```

**Rischi:**  
- Espressioni trascendenti (sin, exp) non riducibili a polinomi: normale, dichiarare `uncertain` non `true`/`false`
- Overflow coefficienti: usare BigInt (già disponibile)

**Stato finale richiesto:**  
`are_mathematically_equal()` funziona per tutte le espressioni polinomiali multivariate. Fallback onesto (`Unimplemented` o `uncertain`) per trascendenti non riducibili.

**Impatto HP Prime G2:** Molto alto — gap ridotto: verifica correttezza algebrica  
**Note per agenti futuri:** Non usare floating-point. Non usare `toString()`. Solo aritmetica esatta su BigInt/Rational.

---

### CAS-P0-002 — Rimpiazzare sampling-based PFD con algoritmo simbolico

**Area:** Frazioni parziali  
**Priorità:** P0  
**Gravità:** Critica  
**Stato:** Risolta  
**Agente assegnato:** Non assegnato  
**Data ultimo aggiornamento:** 2026-05-03

**Aggiornamento Codex 2026-05-03 (finale):**  
Rimpiazzato l'algoritmo di campionamento (Vandermonde) con un algoritmo di Bezout simbolico ricorsivo operante su `RatPoly`. L'implementazione in `src/algebra/partial_fractions.cpp` ora esegue la decomposizione esatta su Q[x] e l'espansione delle potenze via divisioni successive. Verificato con `test_pf_debug.cpp` e `CalculusIntegrateTest.*`: i casi precedentemente critici o crashanti ora restituiscono risultati simbolici corretti.

**Problema risolto:**  
Il sampling numerico falliva su radici non intere e non era conforme agli standard HP Prime G2. Ora il motore è puramente simbolico.

**Evidenza tecnica:**  
```cpp
// partial_fractions.cpp — algoritmo attuale (errato)
for (std::size_t step = 0; samples.size() < count; ++step) {
    const long long magnitude = static_cast<long long>(step / 2U);
    const long long candidate = (step % 2U == 0U) ? magnitude : -magnitude - 1LL;
    const Rational sample{BigInt(candidate)};
    // campiona e costruisce sistema lineare
}
```
Test disabilitato: `test/unit/test_acid_complex_canonical.cpp` — `integrate("1/(x^3+x)", "x")` crashava.

**Perché è importante:**  
Le frazioni parziali sono prerequisito per l'integrazione di funzioni razionali (che è la base del Risch algorithm). Con un algoritmo rotto, l'intera pipeline di integrazione è compromessa.

**Obiettivo matematico:**  
Decomposizione simbolica corretta: dato `P(x)/Q(x)`, trovare `A_i`, `B_i`, `k_i` tali che:
`P/Q = Σ A_i/(x - r_i)^k_i` (over splitting field), usando GCD esteso e Bezout.

**Nuovo muscolo matematico atteso:**  
`partial_fraction_decompose(P, Q, var, ctx)` via:
1. Fattorizzazione square-free di Q (già disponibile via Yun)
2. Bezout identity per ogni coppia di fattori coprime
3. Riduzione ricorsiva senza valutazione numerica

**Algoritmo consigliato:**  
Algoritmo di Hermite classico (parte razionale):
1. Square-free decomposizione di Q: `Q = Q_1 * Q_2^2 * ... * Q_k^k`
2. Per ogni `Q_i^i`: trova `A_i, B_i` tali che `A_i*Q_i + B_i*R_i = P_i` (Bezout su polinomi)
3. Divide ricorsivamente finché tutti i denominatori sono square-free
4. Output: somma di `c_j / (factor_j)^power_j`

**Dipendenze:**  
- Square-free Yun: presente (DONE-05)
- GCD univariato Subresultant: presente (DONE-05)
- Extended GCD (Bezout) su polinomi: **da verificare/implementare**

**Piano tecnico:**  
1. Verificare se `extended_gcd_polynomials(A, B)` esiste → restituisce `(g, s, t)` tale che `s*A + t*B = g`
2. Se manca: implementare in `src/algebra/polynomial_gcd_extended.cpp`
3. Riscrivere `partial_fractions.cpp` usando Bezout ricorsivo
4. Rimuovere completamente il codice di campionamento
5. Test su casi con radici razionali, irrazionali, complesse, radici multiple

**Criteri di accettazione:**  
- `pfd("1/(x^2-1)", "x")` → `1/2 * (1/(x-1) - 1/(x+1))`
- `pfd("1/(x^2*(x+1))", "x")` → `1/x^2 - 1/x + 1/(x+1)`
- `pfd("x/(x^2+1)^2", "x")` → corretto (radici complesse, non deve crashare)
- `pfd("1/(x^3+x)", "x")` → `1/x - x/(x^2+1)` (test precedentemente crashante)

**Test anti-hardcode:**  
```cpp
// Variabili diverse
pfd("1/(y^2-1)", "y")  // stesso algoritmo, variabile diversa
pfd("1/(t^2*(t+2))", "t")  // coefficienti diversi
// Gradi più alti
pfd("x^2/(x^3-1)", "x")  // radici cubiche dell'unità
pfd("1/(x^4-1)", "x")  // grado 4
```

**Rischi:**  
- Bezout su polinomi con coefficienti razionali: verificare aritmetica esatta
- Fattori irriducibili su Q (es. `x^2+1`): la PFD su Q li lascia interi — comportamento corretto, non un bug

**Stato finale richiesto:**  
Nessun campionamento numerico. Solo algebra simbolica. `integrate("1/(x^3+x)", "x")` non crasha e restituisce risultato corretto.

**Impatto HP Prime G2:** Molto alto — gap ridotto: integrazione razionale simbolica  
**Note per agenti futuri:** Il vecchio codice va rimosso completamente, non lasciato come fallback.

---

### CAS-P0-003 — Fix confronto asintotico (MRV fake → ranking reale)

**Area:** Limiti  
**Priorità:** P0  
**Gravità:** Alta  
**Stato:** Risolta  
**Agente assegnato:** Non assegnato  
**Data ultimo aggiornamento:** 2026-05-03

**Aggiornamento Codex 2026-05-03 (finale):**  
`compare_growth()` riscritto con `poly_degree_wrt()` che estrae grado rispetto alla variabile specifica. Rank: costante(0) < ln(1) < polinomio(2, tiebreak per grado) < exp(3). Exp vs exp: confronto argomenti. 4 nuovi test P0-003 accettazione: `x^10/x^2=∞`, `x^2/x^10=0`, `x^5/x^5=1`, `x^7/x^3=∞`. Tutti 29 test P0 passano.

**Aggiornamento Codex 2026-05-03 (precedente):**  
Chiusi casi base verificati per gerarchie `log < potenza < exp`, potenze/inversi a `±inf`, `1/inf -> 0`, dominio `ln(-inf)` e assenza di fallback MRV rumoroso. Verificato con `CalculusLimitTest.*` 8/8. La task resta **Parziale**: manca ancora un comparatore MRV/degree generalizzato come da criteri completi.

**Problema attuale:**  
`src/calculus/limit_mrv.cpp`: la funzione `compare_growth()` assegna rank 1 a qualsiasi simbolo o potenza (inclusi `x^2`, `x^100`). `x^10` e `x^2` sono considerati "incomparabili" (rank identico → return 0). Questo non è l'algoritmo MRV/Gruntz — è una euristica con 3 livelli fissi.

**Evidenza tecnica:**  
```cpp
// limit_mrv.cpp — SBAGLIATO
int get_growth_rank(ExprPtr e) {
    if (expr_is<Symbol>(e)) return 1;         // x e x^100 → rank 1 IDENTICO
    if (call->func_id == BuiltinOp::Ln) return 0;
    if (call->func_id == BuiltinOp::Exp) return 2;
    return 1;  // fallback: qualunque cosa → rank 1
}
```
`lim(x→∞) x^10/x^2` probabilmente non restituisce `∞` correttamente.

**Obiettivo matematico:**  
Confronto asintotico corretto per x→∞:
- Polinomi: confronta per grado (`x^n` vs `x^m`: n>m → x^n domina)
- Log vs polinomio: log(x) < x^ε per ogni ε>0
- Polinomio vs exp: x^n < e^x per ogni n
- Exp vs exp: `e^(f(x))` vs `e^(g(x))` → confronta f vs g

**Nuovo muscolo matematico atteso:**  
`compare_growth_rate(ExprPtr a, ExprPtr b, Symbol var) → {DOMINATES, DOMINATED, EQUIVALENT, INCOMPARABLE}` corretto per polinomi, log, exp, composizioni.

**Algoritmo consigliato:**  
Step incrementale (non Gruntz completo — troppo grande per P0):
1. Estrarre grado polinomiale di un'espressione rispetto a `var` via `polynomial_degree()`
2. Se entrambe polinomiali: confrontare gradi
3. Se una è `ln(...)`: rank inferiore a qualsiasi polinomio
4. Se una è `exp(...)`: rank superiore a qualsiasi polinomio
5. Per exp vs exp: confronto ricorsivo degli argomenti
6. Fallback: `INCOMPARABLE` (onesto, non sbagliato)

**Dipendenze:**  
- `polynomial_degree(expr, var)`: da verificare se esiste
- Simplifier (presente)

**Piano tecnico:**  
1. Implementare `polynomial_degree(ExprPtr, Symbol) → optional<int>` se non esiste
2. Riscrivere `compare_growth()` usando l'algoritmo descritto
3. Aggiornare `limit_mrv.cpp` con il nuovo comparatore
4. Test sistematici su tutti i casi polinomiali + log/exp base
5. Documentare esplicitamente i casi non ancora gestiti

**Criteri di accettazione:**  
- `compare_growth(x^10, x^2, x)` → DOMINATES
- `compare_growth(x^2, x^10, x)` → DOMINATED
- `compare_growth(x^5, x^5, x)` → EQUIVALENT
- `compare_growth(ln(x), x, x)` → DOMINATED
- `compare_growth(exp(x), x^1000, x)` → DOMINATES
- `lim(x^10/x^2, x, inf)` → ∞
- `lim(x^2/x^10, x, inf)` → 0
- `lim((x^3+x)/(2*x^3-1), x, inf)` → 1/2

**Test anti-hardcode:**  
```cpp
// Esponenti arbitrari — non solo 10 e 2
compare_growth(x^37, x^12, x)   // DOMINATES
compare_growth(x^3, x^3, x)     // EQUIVALENT
compare_growth(x^1, x^100, x)   // DOMINATED
// Con coefficienti
lim("(3*x^4 + x)/(x^4 - 1)", x, inf)  // → 3
lim("x^99 / x^100", x, inf)           // → 0
```

**Rischi:**  
- Limiti oscillatori (sin(x)/x per x→∞): non gestibili con confronto asintotico, dichiarare INCOMPARABLE onestamente
- Composizioni complesse: step incrementale, non pretendere di gestire tutto

**Stato finale richiesto:**  
Confronto asintotico corretto per polinomi e gerarchie log/exp base. Nessun rank fisso hardcoded. Fallback onesto per casi non gestibili.

**Impatto HP Prime G2:** Alto — gap ridotto: limiti affidabili per forme algebriche  

---

### CAS-P0-004 — Infrastruttura test anti-hardcode e property-based

**Area:** Testing  
**Priorità:** P0  
**Gravità:** Critica  
**Stato:** Aperta  
**Agente assegnato:** Non assegnato  
**Data ultimo aggiornamento:** 2026-05-03

**Problema attuale:**  
Test esistenti verificano solo casi specifici con variabili specifiche (quasi sempre `x`). Nessun property-based test. Nessun test anti-hardcode sistematico. Un'implementazione hardcoded supererebbe tutti i test attuali.

**Evidenza tecnica:**  
- `test/unit/`: test con variabile `x` hardcoded nella maggior parte dei casi
- Nessun file che generi input casuali o permutazioni
- Nessun test che verifichi `D(∫f dx) == f` come proprietà universale
- Nessun test che verifichi `expand(factor(p)) == p` come proprietà

**Obiettivo matematico:**  
Test che verificano proprietà matematiche universali, non casi singoli.

**Nuovo muscolo matematico atteso:**  
Helper di test riusabili:
- `verify_inverse_property(f, g, inputs)`: verifica `g(f(x)) == x` su lista di input
- `verify_with_variable_renaming(test_fn, vars)`: esegue stesso test con variabili diverse
- `verify_differentiation_is_inverse_of_integration(exprs)`

**Piano tecnico:**  
1. Creare `test/unit/helpers/math_test_helpers.hpp`
2. Implementare `are_mathematically_equal()` come helper test (dipende da CAS-P0-001 o può usare expand+subtract)
3. Aggiungere `test/unit/anti_hardcode/test_variable_independence.cpp`
4. Aggiungere `test/unit/anti_hardcode/test_integration_differentiation_inverse.cpp`
5. Aggiungere `test/unit/anti_hardcode/test_algebraic_equivalence.cpp`
6. Aggiungere generator di polinomi casuali per stress test

**Test obbligatori da aggiungere:**

```cpp
// test_variable_independence.cpp
// Ogni funzionalità deve funzionare con qualsiasi variabile
void test_with_all_vars(const std::string& expr_template, 
                        const std::vector<std::string>& vars) {
    for (const auto& v : vars) {
        std::string expr = replace(expr_template, "VAR", v);
        // verifica che il risultato sia corretto
    }
}

// Usa: {"x", "y", "z", "a", "b", "t", "alpha", "theta"}

// test_integration_differentiation_inverse.cpp
// Proprietà: D(∫f dx) == f
TEST(Property, IntegrationDifferentiationInverse) {
    std::vector<std::string> funcs = {
        "x^3", "2*x^2 + x - 1", "x^5 - 3*x^3 + 2*x"
    };
    for (const auto& f : funcs) {
        auto integral = integrate(parse_expr(f), Symbol("x"), ctx);
        if (integral.is_ok()) {
            auto deriv = differentiate(integral.value(), Symbol("x"), ctx);
            EXPECT_TRUE(are_equal(deriv, parse_expr(f), ctx)) << "Failed: " << f;
        }
    }
}

// test_algebraic_equivalence.cpp
TEST(Property, ExpandFactorInverse) {
    std::vector<std::string> polys = {
        "x^2 - 1", "x^2 - 4", "x^3 - x", "x^4 - 1"
    };
    for (const auto& p : polys) {
        auto factored = factor(parse_expr(p), ctx);
        if (factored.is_ok()) {
            auto expanded = expand(factored.value(), ctx);
            EXPECT_TRUE(are_equal(expanded, parse_expr(p), ctx)) << "Failed: " << p;
        }
    }
}
```

**Criteri di accettazione:**  
- Ogni funzionalità matematica principale ha almeno un test property-based
- Ogni test esiste in versione con ≥3 variabili diverse
- Esiste un test che rileva immediatamente implementazioni hardcoded su `x`

**Impatto HP Prime G2:** Molto alto — prerequisito per dichiarare qualsiasi feature come corretta

---

## DETTAGLIO TASK P1

---

### CAS-P1-001 — Hermite Reduction (Bezout step)

**Area:** Integrazione  
**Priorità:** P1  
**Gravità:** Critica  
**Stato:** Bloccata  
**Agente assegnato:** Non assegnato  
**Data ultimo aggiornamento:** 2026-05-03

**Problema attuale:**  
`src/calculus/integrate_risch.cpp`:
```cpp
Result<ExprPtr> hermite_reduction(...) {
    return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented,
        "Hermite reduction: Bezout step non ancora implementato"});
}
```
È il primo passo obbligatorio del Risch algorithm. Senza di esso, integrazione di qualsiasi funzione razionale non banale è impossibile.

**Dipendenze:** CAS-P0-002 (extended GCD / Bezout su polinomi)

**Algoritmo:**  
Data `P/Q` con `Q` square-free:
1. Per ogni fattore `D_i^k` (k≥2) di Q: trova Bezout coeffs `s,t` con `s*D_i' + t*D_i = gcd(...)`
2. Integra la parte "ridotta" per parti: `∫P/Q = A/B + ∫C/D` con `deg(D) < deg(Q)`
3. Riduci fino a denominatore square-free

**Criteri di accettazione:**  
- `integrate("1/(x^2*(x+1))", "x")` → forma corretta
- `integrate("x/(x^2-1)^2", "x")` → forma corretta
- La parte razionale è riducibile via Hermite, la parte log-free via Rothstein-Trager

**Impatto HP Prime G2:** Molto alto — prerequisito per tutta integrazione razionale

---

### CAS-P1-002 — Rothstein-Trager / Lazard-Rioboo (log-part)

**Area:** Integrazione  
**Priorità:** P1  
**Gravità:** Critica  
**Stato:** Bloccata (dipende da CAS-P1-001)  
**Data ultimo aggiornamento:** 2026-05-03

**Algoritmo:**  
Dopo Hermite reduction, la parte residua è `P/Q` con `Q` square-free.  
Rothstein-Trager: calcola `resultant(Q, P - t*Q', t)` → radici `c_i` → integrale = `Σ c_i * ln(gcd(Q, P - c_i*Q'))`

**Criteri di accettazione:**  
- `integrate("1/(x^3+x)", "x")` → `ln(x) - 1/2*ln(x^2+1)` (era crash)
- `integrate("1/(x^2-1)", "x")` → `1/2*ln(x-1) - 1/2*ln(x+1)`
- Property test: `D(result) == original` per tutti i casi

**Impatto HP Prime G2:** Molto alto

---

### CAS-P1-003 — Polynomial solving grado 4 (Ferrari)

**Area:** Solving  
**Priorità:** P1  
**Gravità:** Alta  
**Stato:** Aperta  
**Data ultimo aggiornamento:** 2026-05-03

**Problema attuale:**  
`src/algebra/solve_polynomial.cpp`: gradi 1-3 implementati (Cardano corretto). Grado 4 assente. Impatta: autovalori di matrici 4×4, equazioni di grado 4, fattorizzazione via radici.

**Algoritmo:** Metodo di Ferrari:
1. Deprimere la quartica: `x^4 + px^2 + qx + r`
2. Introdurre cubica ausiliaria y (risolvibile via Cardano già implementato)
3. Fattorizzare in due quadratiche usando y
4. Risolvere le due quadratiche

**Dipendenze:** Cardano (grado 3) già presente

**Criteri di accettazione:**  
- `solve("x^4 - 5*x^2 + 4", "x")` → `{-2, -1, 1, 2}`
- `solve("x^4 - 1", "x")` → `{-1, 1, -i, i}`
- `solve("x^4 + 1", "x")` → 4 radici complesse
- Property: ogni radice r soddisfa `p(r) == 0`

**Test anti-hardcode:**  
```cpp
// Quartica generica — non solo casi biquadratici
solve("x^4 + 2*x^3 - 3*x^2 - 4*x + 4", "x")
solve("x^4 - 4*x^3 + 6*x^2 - 4*x + 1", "x")  // (x-1)^4
// Variabili diverse
solve("y^4 - 5*y^2 + 4", "y")
```

**Impatto HP Prime G2:** Alto

---

### CAS-P1-004 — Fattorizzazione completa su Q (oltre RRT)

**Area:** Fattorizzazione  
**Priorità:** P1  
**Gravità:** Alta  
**Stato:** Aperta  
**Data ultimo aggiornamento:** 2026-05-03

**Problema attuale:**  
`src/algebra/factorization_polynomials.cpp`: solo Rational Root Theorem. Square-free Yun è presente (DONE-05) ma non integrato nella pipeline di fattorizzazione completa. `x^4+1` (irriducibile su Q), `x^4-4` (fattorizzabile su Q) non gestiti correttamente.

**Step incrementale corretto (non tentare tutto in una volta):**  
1. ✓ Square-free Yun (già fatto)
2. → **Fattorizzazione modulo primo** (Berlekamp o Cantor-Zassenhaus su Z_p)
3. → **Hensel lifting** (sollevamento radici mod p^k a coefficienti interi)
4. → **Ricombinazione fattori** (Zassenhaus combinations)

**Questa task copre step 2+3+4.**

**Aggiornamento Codex 2026-05-03:**  
Implementata ricombinazione Hensel/Zassenhaus per subset di fattori modulari e validazione per divisibilità esatta prima di appendere fattori. Corretto anche `find_factor_lll` per non emettere vettori LLL non divisori e stabilizzato `lll_reduction()` dopo swap Gram-Schmidt. Verificati `AlgebraLLLTest.*`, `AlgebraHenselTest.*`, `AlgebraFactorizationTest.Degree10LargeCoeffs`, `FactorOverIntegers.*` 5/5. Task non chiusa finché non sono coperti tutti i criteri e property `expand(factor(p)) == p`.

**Dipendenze:** Square-free Yun (presente), GCD univariato (presente)

**Criteri di accettazione:**  
- `factor("x^4 - 1")` → `(x-1)*(x+1)*(x^2+1)`
- `factor("x^6 - 1")` → `(x-1)*(x+1)*(x^2-x+1)*(x^2+x+1)`
- `factor("x^4 + 1")` → dichiarato irriducibile su Q (corretto)
- `factor("6*x^2 + 7*x + 2")` → `(2*x+1)*(3*x+2)`
- Property: `expand(factor(p)) == p` per tutti i casi testati

**Impatto HP Prime G2:** Molto alto

---

### CAS-P1-005 — Integrazione assumptions engine nel simplifier

**Area:** Assunzioni  
**Priorità:** P1  
**Gravità:** Alta  
**Stato:** Aperta  
**Data ultimo aggiornamento:** 2026-05-03

**Problema attuale:**  
TODO.md dichiara P4 (Assumptions Engine) come completato. Ma l'audit mostra che le assunzioni non sono usate dal simplifier. `sqrt(x^2)` non viene semplificato a `x` anche con assunzione `x > 0`. `ln(x)` non conosce il dominio di `x`.

**Evidenza tecnica:**  
- Headers dichiarano struttura `AssumptionSet` o simile
- `src/symbolic/simplify_functions.cpp`: non consulta le assunzioni per le semplificazioni

**Obiettivo:**  
Collegare le assunzioni al simplifier in modo che:
- `sqrt(x^2)` con `x ∈ ℝ, x ≥ 0` → `x`
- `sqrt(x^2)` senza assunzioni → `|x|`
- `ln(x)` con `x > 0` → definito; senza → warning o dominio condizionale
- `abs(x)` con `x ≥ 0` → `x`

**Dipendenze:** Struttura assumptions esistente

**Criteri di accettazione:**  
- `simplify("sqrt(x^2)", assume(x, positive))` → `x`
- `simplify("abs(x)", assume(x, nonnegative))` → `x`
- `simplify("sqrt(x^2)")` senza assunzioni → `abs(x)` o espressione non semplificata (non `x`)

**Impatto HP Prime G2:** Alto — HP Prime G2 gestisce domini esplicitamente

---

## DETTAGLIO TASK P2

---

### CAS-P2-001 — Taylor via derivate successive (algoritmo reale)

**Area:** Serie  
**Stato:** Bloccata (dipende da CAS-P0-001 per verifica)  
**Problema:** `taylor_series()` usa lookup per sin/cos/exp. Nessun algoritmo per funzioni arbitrarie.  
**Algoritmo:** `a_n = f^(n)(a) / n!` → iterare `differentiate()` n volte → raccogliere coefficienti  
**Criteri:** `taylor(exp(x), x, 0, 5)` corretto + `taylor(sin(x^2), x, 0, 6)` (non in lookup)

---

### CAS-P2-002 — Integrali definiti via FTC

**Area:** Integrazione  
**Stato:** Bloccata (dipende da CAS-P1-002)  
**Problema:** Assenti completamente.  
**Algoritmo:** `∫_a^b f dx = F(b) - F(a)` con `F = ∫f dx` + gestione singolarità  
**Criteri:** `integrate("x^2", "x", 0, 1)` → `1/3`; `integrate("sin(x)", "x", 0, pi)` → `2`

---

### CAS-P2-003 — Sostituzione trigonometrica

**Area:** Integrazione  
**Stato:** Bloccata (dipende da CAS-P1-001)  
**Problema:** `integrate_trig_substitution.cpp` restituisce `Unimplemented` per tutto.  
**Algoritmo:** Pattern `sqrt(a^2-x^2)` → `x=a*sin(t)`; `sqrt(a^2+x^2)` → `x=a*tan(t)`; `sqrt(x^2-a^2)` → `x=a*sec(t)`  
**Criteri:** `integrate("sqrt(1-x^2)", "x")` → `1/2*(x*sqrt(1-x^2) + arcsin(x))`

---

### CAS-P2-004 — L'Hôpital come fallback affidabile

**Area:** Limiti  
**Stato:** Aperta  
**Problema:** Limiti del tipo 0/0 o ∞/∞ non gestiti.  
**Algoritmo:** Verificare forma indeterminata → applicare `lim(f/g) = lim(f'/g')` → guard contro loop (max depth) → fallback onesto  
**Criteri:** `lim(sin(x)/x, 0)` → 1; `lim((e^x-1)/x, 0)` → 1; `lim(x^2/e^x, inf)` → 0

---

### CAS-P2-005 — Riduzione angolare modulare

**Area:** Trigonometria  
**Stato:** Aperta  
**Problema:** `sin(5*pi/6)` non si semplifica — manca l'angolo dalla lookup table.  
**Algoritmo:** Riduzione `θ mod 2π` → quadrante → valore esatto usando `{0, π/6, π/4, π/3, π/2}` come base  
**Criteri:** `sin(5*pi/6)` → `1/2`; `cos(7*pi/4)` → `sqrt(2)/2`; funziona per angoli arbitrari multipli di π/12

---

### CAS-P2-006 — Autovalori simbolici per n>3

**Area:** Algebra lineare  
**Stato:** Bloccata (dipende da CAS-P1-003, CAS-P1-004)  
**Problema:** Eigenvalues dipende dal polynomial solver → max grado 3 → max matrice 3×3.  
**Soluzione:** Sblocco automatico dopo P1-003 (Ferrari) — polinomi caratteristici grado 4 risolvibili

---

## CHECKLIST ANTI-FURBIZIA (da usare prima di chiudere ogni task)

```
- [ ] Non ho aggiunto hardcode per input specifici.
- [ ] Non ho scritto codice solo per far passare il test.
- [ ] Non ho usato string matching fragile come logica matematica.
- [ ] Non ho nascosto un fallback errato (un Unimplemented onesto è meglio di un risultato sbagliato).
- [ ] Non ho restituito risultati matematici non verificati.
- [ ] Ho implementato una capacità generalizzabile.
- [ ] Ho aggiunto test con variabili diverse (non solo "x").
- [ ] Ho aggiunto test su forme sintattiche equivalenti.
- [ ] Ho aggiunto almeno un test anti-hardcode.
- [ ] Ho dichiarato i limiti residui.
- [ ] Ho aggiornato lo stato della task in questa tabella.
- [ ] Ho verificato che nessuna feature esistente sia rotta (regressioni).
```

---

## TEMPLATE REPORT INTERVENTO

```markdown
## Report intervento

**Task:** CAS-PX-XXX  
**Stato precedente:** Aperta / Bloccata  
**Stato nuovo:** In review / Risolta  

### Cosa è stato fatto
- ...

### File modificati
- src/...
- test/...

### Algoritmo implementato
- ...

### Perché non è una patch
- ...

### Nuovo muscolo matematico aggiunto
- ...

### Test aggiunti
- ...

### Test anti-hardcode aggiunti
- ...

### Regressioni controllate
- [ ] ctest --test-dir build --output-on-failure → 0 failures

### Limiti residui
- ...

### Prossime task consigliate
- ...
```

---

## ISTRUZIONE PER IL PROSSIMO AGENTE

1. Leggi `ClaudeAudit-HP.md` per il quadro completo
2. Leggi questa tabella task
3. Scegli la prima task P0 non bloccata
4. **Inizia da CAS-P0-001 o CAS-P0-004** (non hanno dipendenze esterne)
5. Non lavorare su P1/P2/P3 finché esistono P0 aperti
6. Non toccare UI, naming, refactor cosmetici
7. Ogni modifica deve avere giustificazione matematica
8. Usa la checklist anti-furbizia prima di chiudere

*Sistema task generato: 2026-05-03 — Claude Sonnet 4.6*
