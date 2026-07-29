# Executive Summary: Real-CAS-Engine C++ Architecture Audit (v2.0)

**Giudizio Complessivo:** CAS con fondazioni C++ solide (AST immutabile, BigInt, Assumptions base) fino al livello dell'algebra polinomiale su $\mathbb{Z}$. Soffre di un design "rigido" per i domini algebrici e presenta scorciatoie architetturali gravi (fallback a `double`) nei moduli di calcolo avanzato (Integration).
**Stato Parità vs HP Prime G2:** ~30-40% (Coerente con il piano F3).
**Nota di correzione:** Un precedente audit automatizzato ha erroneamente classificato come assente il sistema di Assumptions. Una verifica manuale e profonda del codice (`src/symbolic/assumptions.cpp` e `simplify_exp_log.cpp`) conferma che il motore deduce correttamente le limitazioni (es. `sqrt(x^2) -> abs(x)` gestito). I difetti rimanenti, tuttavia, sono reali e documentati di seguito riga per riga.

---

# Matrice Comparativa Reale e Scrupolosa

| Area | Stato nel progetto (Evidenza Codice) | Stato Giac / HP Prime | Gap Reale | Voto | Rischio mat. | Intervento Correttivo |
|---|---|---|---|---|---|---|
| **Aritmetica e AST** | Eccellente. `BigInt` custom limb-based. `AstArena` garantisce allocazione rapida e `ExprPtr` immutabili. | Multi-precisione ottimizzata in C. | Performance su FFT enormi. | A | Basso | Ottimizzazioni di basso livello (Fase post-parità). |
| **Assumptions** | Esiste `class Assumptions` con inferenza logica (`>`,`<`, `Real`, `Positive`). Gestisce `sqrt(x^2)`. | Deduzioni su domini misti, branch cut completi nel piano complesso. | Tracking di polidromia complessa. | B | Medio | Estendere a branch cuts complessi rigorosi. |
| **Algebra Polinomiale** | Implementata `MultivariatePolynomial` usando `std::vector<MultivariateTerm>` fissato su `BigInt`. | Sistema a domini generico (es. `FractionField(PolyRing(ZZ))`). | **Architetturale.** L'hardcoding su `BigInt` rende l'algebra in $\mathbb{Q}$ o $\mathbb{C}$ dipendente da hack (LRT) o da AST. | C | Medio | Rifattorizzare come `Polynomial<Ring>`. |
| **Gröbner Basis (F4)** | Implementato F4 (`polynomial_groebner_f4.cpp`). FGLM presente. | Solido e usato come motore di default per `solve`. | In `csolve.cpp:71`, se F4 fallisce su un sistema 2x2, il codice fa fallback su Resultant. Questo maschera bug dell'F4. | B- | Basso | Rimuovere il fallback in `csolve`; F4 deve poter chiudere i sistemi 2x2 nativamente. |
| **Solvers Polinomiali** | Risolve fino al grado 3. **Il grado 4 è commentato/disabilitato** in `solve_polynomial.cpp:478`. `RootOf` limitato. | Formule esatte per grado 4; `RootOf` completo su $\mathbb{Q}(\alpha)$. | Mancano la quartica e l'algebra dei campi di spezzamento per `RootOf`. | C | Basso | Ripristinare Ferrari/Lodovico per la quartica. Generalizzare `RootOf`. |
| **Integrazione** | Risch framework abbozzato. `integrate.cpp` usa **`double` e `M_PI`** per rilevare le singolarità trigonometriche (L254-L268). | Integrazione esatta di Liouville-Risch senza alcuna contaminazione floating-point. | **Critico.** La dipendenza dal floating-point (M_PI, stod) distrugge il contratto di un CAS esatto. | D | Critico | Sostituire il controllo `dlo/dhi` con la Sturm sequence o isolamento radici reali esatte. |
| **Limiti e Serie** | Gruntz algorithm implementato ma con limiti arbitrari (`limit_mrv.cpp:130` cap ricorsione a 1024). Bails out onesticamente su timeout. | Espansioni asintotiche formali di ordine variabile. | Dipendenza da limiti di profondità anziché da teoria di ordine asintotico. | C | Basso | Passare da euristiche di ricorsione a `FormalPowerSeries`. |

---

# Focus: Hardcode e Scorciatoie Tecniche (Verificate nel sorgente)

Questi sono i punti esatti in cui il progetto "bara" o adotta soluzioni fragili per far passare i test o per coprire mancanze algoritmiche.

### 1. Infiltrazione di `double` nell'Integrazione Simbolica (VIOLAZIONE REGOLA 1)
*   **File:** `src/calculus/integrate.cpp`
*   **Righe:** 247-268
*   **Codice:** 
    ```cpp
    double n = std::stod(r.numerator().decimal());
    double c = rat_to_double(linear->first);
    double x_base = (M_PI / 2.0 - d) / c;
    ```
*   **Problema:** Per capire se un integrale definito come $\int \tan(x) dx$ ha un polo nell'intervallo, il motore converte gli estremi in `double` e usa `M_PI`. In un CAS, se l'utente chiede limiti contenenti costanti simboliche o frazioni infinite, il `double` perde precisione e può dichiarare falso un polo reale (o inventarne uno).
*   **Fix:** Costruire un valutatore di intervalli rigoroso (Interval Arithmetic) o usare root isolation tramite sequenze di Sturm senza conversione a `double`.

### 2. Disattivazione Silenziosa della Quartica
*   **File:** `src/algebra/solve_polynomial.cpp`
*   **Righe:** 478
*   **Codice:** `if (deg == 4U) { /* disattivata */ }`
*   **Problema:** L'equazione di 4° grado ha una formula chiusa ben nota. Commentarla costringe l'engine a declassare la risposta a `RootOf` inerte o fallire, causando un gap immotivato per un CAS di livello L2/L3.
*   **Fix:** Implementare il risolutore quartico (metodo di Ferrari) con corretta scelta dei branch cuts per la risolvente cubica.

### 3. Fallback Gröbner 2x2 su Subresultant
*   **File:** `src/algebra/csolve.cpp`
*   **Righe:** 71-96
*   **Problema:** Il sistema di equazioni tenta l'F4 solver. Se questo restituisce un errore, e il sistema è 2x2, interviene un hack che usa l'eliminazione per risultanti. Questo è un "cerotto" (band-aid) su un bug non diagnosticato nel sistema Gröbner.
*   **Fix:** Capire perché F4 fallisce su casi test 2x2 e correggere F4/FGLM alla radice.

### 4. Gestione Budget "Onesta ma Limitante" (I numerosi `Unimplemented`)
*   **File:** `polynomial_gcd_zippel_prony.cpp`, `factor_multivariate_hensel.cpp`, etc.
*   **Problema:** Il codice è pieno di `return make_unimplemented(...)`. Dal punto di vista ingegneristico è corretto (non dà risultati sbagliati, `Result<T>` è usato bene). Dal punto di vista matematico, significa che l'algoritmo non è completo e che appena si esce dalla "comfort zone" dei test case semplici, l'algoritmo alza bandiera bianca.
*   **Fix:** Manca l'estensione multi-prime per Zippel e un lifting robusto in Hensel.

---

# Gap Strutturali Architetturali Reali

1.  **Mancanza di Astrazione sui Coefficienti (Domain System)**
    Attualmente `MultivariatePolynomial` mappa su `BigInt`. Per fare polinomi su $\mathbb{Q}$, il motore aggira il problema portando a denominatore comune (content) e usando l'algebra in $\mathbb{Z}$. Per estensioni algebriche $\mathbb{Q}(\alpha)$, deve ricorrere all'AST (`ExprPtr`). In Giac/HP Prime, l'algebra lavora su template/interfacce generiche `Ring`. Questo è il più grande freno all'espansione.

2.  **Integrazione Risch Limitata**
    L'infrastruttura Risch c'è (`integrate_risch.cpp`), ma molti sotto-algoritmi (es. estensioni logaritmiche generali) terminano con `Unimplemented`. L'integrazione di funzioni elementari si affida ancora a moduli di pattern matching (`integrate_definite_patterns.cpp`).

---

# Conclusione e Priorità Tecniche

**Valutazione Oggettiva:**
Il progetto **NON mi sta prendendo in giro**: l'architettura L0-L2 esiste. Hai costruito un motore solido sul lato C++ (memoria, AST, type-safety con Result monadico) e implementato algoritmi seri (GCD modulare, F4, Risch-Norman parziale, limit MRV, Assumptions).

Ma per avvicinarsi a **Giac**, devi estirpare le scorciatoie (hack) inserite per far passare i test.

**Le 3 azioni imperative:**
1.  **P0 - Epura i `double`:** Vai in `src/calculus/integrate.cpp` ed elimina `stod` e `M_PI`. Scrivi un root isolator o un boundary check interamente simbolico/razionale. È una questione di principio per un CAS.
2.  **P1 - F4 e Quartica:** Rimuovi il paracadute del 2x2 in `csolve.cpp` e fai funzionare F4. De-commenta e implementa il risolutore di grado 4.
3.  **P2 - Astrazione Polinomiale:** Crea una vera classe template `Polynomial<T>` in modo da non essere vincolato a `BigInt` per coefficienti, rendendo triviali i calcoli su $\mathbb{Q}$ o $\mathbb{Z}_p$.
