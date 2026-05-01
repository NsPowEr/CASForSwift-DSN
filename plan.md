Piano Implementazione CAS — Target HP Prime G2
[REV 3 — audit chirurgico 2026-04-29, 3 agenti in parallelo]

Stato ACID post-REV2: 18/21 pass. Falliscono: Test 9, 19, 28.
Fasi 0, 5 implementate correttamente. Fase 3 già risolta (binomial ok, Test 5 = 64ms).
Questo documento sostituisce completamente REV 2.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SEZIONE A — CAUSE RADICE ESATTE DEI 3 TEST CHE FALLISCONO
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TEST 28 (ln(ln(x+e))/ln(ln(x)) → 1)
  Causa root: limit.cpp:51-56 chiama compute_limit_mrv prima di qualsiasi altro
  path. MRV sostituisce x=1/w, ottiene ln(ln(1/w+e))/ln(ln(1/w)). taylor_series
  su questa espressione (che ha singolarità essenziale in w=0) restituisce
  ok(polynomial=0). Quindi compute_limit_mrv restituisce ok(0). Essendo ok, la
  funzione compute() ritorna immediatamente con 0 senza raggiungere la via
  try_log_log_limit (che sarebbe corretta).

  try_log_log_limit (limit.cpp:187-218) è matematicamente corretta: riceve
  ln(A)/ln(B), calcola lim(A/B) ricorsivamente, se finito non-zero ritorna 1.
  Il problema è che viene chiamata solo a limit.cpp:173 con guard
  !limit_is_infinity(point) — quindi MAI per x→∞.

  Fix minimo: 6 righe in limit.cpp, PRIMA della chiamata MRV a riga 51.

TEST 19 (sistema Gröbner)
  Causa root: polynomial_groebner.cpp restituisce vettore vuoto (stub, riga 24).
  f4_groebner è implementata correttamente in polynomial_groebner_f4.cpp:135-239
  ma è una free function NON dichiarata in nessun header → linkage impossibile
  senza forward declaration. solve_nonlinear_system_f4 (riga 278-327) chiama
  f4_groebner ma poi IGNORA il risultato e hardcoda la soluzione Test 8.
  MultivariatePolynomial rifiuta coefficienti razionali (riga 37).
  MonomialOrder enum: non esiste nel codebase (solo MonomialLexComparator hardcoded).

TEST 9 (radici ciclotomiche x^6-1)
  Causa root: solve_polynomial.cpp:488-494 emette RootOf generici per gradi ≥ 5
  senza tentare riconoscimento ciclotomico. Nessun file polynomial_cyclotomic.cpp
  esiste. x^6-1 = Φ1·Φ2·Φ3·Φ6 = (x-1)(x+1)(x²+x+1)(x²-x+1) — tutti di grado ≤ 4,
  quindi già risolubili. Il blocco è prima: solve_by_factoring non è chiamato
  correttamente per x^6-1 (va a RootOf prima di tentare fattorizzazione).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SEZIONE B — PANORAMICA FASI (AGGIORNATA REV 3)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Priorità │ Fase │ Problema              │ Effort │ ACID test sbloccato
─────────┼──────┼───────────────────────┼────────┼────────────────────
  1      │  2b  │ Test 28 (log limit)   │ 0.5 d  │ Test 28
  2      │   9  │ Test 9 (ciclotomico)  │ 1 d    │ Test 9
  3      │   1  │ Gröbner F4 wiring     │ 4 d    │ Test 19
  4      │   4  │ csolve N-D            │ 2 d    │ (dipende da Fase 1)
  5      │   6  │ ODE const-coeff       │ 4 d    │ ODE lineari 2° ordine
  6      │  2c  │ Serie di potenze full │ 4 d    │ ODE Frobenius, Risch fbk
  7      │   7  │ Gosper + Zeilberger   │ 5 d    │ somme ipergeometriche
  8      │  6b  │ Kovacic L=1,2         │ 5 d    │ ODE razionali 2° ordine

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASE 2b — FIX TEST 28: lim ln(ln(x+e))/ln(ln(x)) = 1
           Effort: 0.5 giorni
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

File: src/calculus/limit.cpp

Inserire QUESTE RIGHE tra la riga 49 (fine del blocco try_infinite_limit)
e la riga 51 (inizio if point_is_pos_inf → compute_limit_mrv):

    // ln(A)/ln(B) → 1 quando lim(A/B) finito e non-zero, anche per x→∞
    // DEVE precedere MRV: MRV fallisce su nested logs (taylor_series restituisce
    // ok(0) erroneamente per espressioni con singolarità essenziale in ω=0).
    if (point_is_pos_inf || point_is_neg_inf) {
        auto q_inf = extract_quotient_view(simplified_expr.value(), arena_);
        if (q_inf.has_value()) {
            LimitDirection inf_dir = point_is_pos_inf
                ? LimitDirection::Right : LimitDirection::Left;
            if (auto ll = try_log_log_limit(q_inf.value(), var,
                                             simplified_point.value(), inf_dir, 0U)) {
                return ll.value();
            }
        }
    }

Perché funziona:
  try_log_log_limit(ln(ln(x+e))/ln(ln(x)), x, ∞):
    a = ln(x+e), b = ln(x)
    Verifica a→∞ e b→∞: substitute x=∞ → ln(∞)=∞ ✓
    Calcola L = compute_recursive(ln(x+e)/ln(x), x, ∞, Right, 1):
      Questo è ancora ln/ln → try_log_log_limit con a2=x+e, b2=x (depth=2)
      L2 = compute_recursive((x+e)/x, x, ∞, Right, 2)
           = compute_recursive(1+e/x, x, ∞, Right, 2)
           → substitute 1+e/∞ = 1 ✓ → restituisce ok(1)
      L2=1, non zero, non ∞ → return ok(1) ✓
    L=1, non zero, non ∞ → return ok(1) ✓
  Test 28 passa. ✓

ZERO rischio regressione: try_log_log_limit ha guard interno depth>3U.
Se l'espressione non è ln(A)/ln(B) restituisce std::nullopt e la
funzione non viene chiamata.

Test di non-regressione da eseguire subito dopo:
  Test 1 (Gruntz lim x→0) — deve restare ✓
  Test 17 (Squeeze lim x→0) — deve restare ✓
  Test 26 (Schanuel limits) — deve restare ✓

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASE 9 — FIX TEST 9: solve(x^6-1, x) → radici esplicite
          Effort: 1 giorno
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Causa: solve_polynomial.cpp:488-494 emette RootOf per gradi ≥ 5 SENZA
tentare fattorizzazione. solve_by_factoring è chiamata solo per deg ≤ 4.
x^6-1 fattorizza in polinomi di grado ≤ 2, tutti risolubili.

AUDIT NECESSARIO prima di scrivere codice:
  Leggere solve_polynomial.cpp righe 520-548 (funzione solve_polynomial).
  Leggere solve_polynomial.cpp righe 419-495 (funzione solve_factor).
  Identificare esattamente dove è il branch che salta a RootOf per grado 6.

Sub-task A — Estendere solve_by_factoring a gradi arbitrari
  File: src/algebra/solve_polynomial.cpp
  Funzione: solve_polynomial() (riga ~497)

  Attuale logica (circa righe 524-548):
    if (degree <= 4) → dispatch a solve_degree_N o solve_by_factoring
    else → RootOf emission loop

  Fix: PRIMA dell'emissione RootOf, tentare solve_by_factoring:
    auto factored = solve_by_factoring(poly, var, ctx);
    if (factored.is_ok() && !factored.value().empty()) return factored;
    // solo se factoring fallisce → RootOf

  Verificare che solve_by_factoring (già esistente) chiami
  factorization_polynomials e poi ricorsivamente solve_factor su ogni fattore.

Sub-task B — Aggiungere riconoscimento ciclotomico come fallback elegante
  (opzionale ma migliora output)

  File nuovo: src/algebra/polynomial_cyclotomic.cpp (~ 100 LoC)

  Tabella Φ_m per m = 1..30 (questi coprono tutti i casi pratici):
    Φ_1 = x-1, Φ_2 = x+1, Φ_3 = x²+x+1, Φ_4 = x²+1,
    Φ_5 = x⁴+x³+x²+x+1, Φ_6 = x²-x+1, Φ_7 = x⁶+...+1,
    Φ_8 = x⁴+1, Φ_9 = x⁶+x³+1, Φ_10 = x⁴-x³+x²-x+1, ...

  Costruzione Φ_m via Möbius inversion:
    Φ_m(x) = ∏_{d|m} (x^d - 1)^{μ(m/d)}
    dove μ è la funzione di Möbius.
    Precomputa per m ≤ 30 al primo utilizzo (cache statica).

  Funzione: optional<int> is_cyclotomic(const IntPoly& p)
    Per ogni m tale che φ(m) == deg(p): confronta p con Φ_m(x)
    Se coincide: return m

  Funzione: vector<ExprPtr> cyclotomic_roots(int m, Symbol x, AstArena& a)
    Per k in [1, m-1] con gcd(k,m)==1:
      root = Binary(Pow, Constant(E),
                    Binary(Mul, Constant(I),
                           Binary(Mul, RationalLit(2k, m), Constant(Pi))))
    Return roots (φ(m) elementi)

  Integrazione: solve_factor() PRIMA dell'emissione RootOf:
    Se IntPoly ottenibile: test is_cyclotomic → emetti radici esplicite

  Nota: per Test 9 la soluzione corretta avviene già con Sub-task A (solve_by_factoring).
  Il riconoscimento ciclotomico è un miglioramento estetico dell'output (RootOf → e^(2πik/n)).

Test per Fase 9:
  solve(x^6-1, x) → 6 radici esplicite (non RootOf)
  solve(x^4-1, x) → {1,-1,i,-i} (già funziona, regression check)
  solve(x^2+x+1, x) → 2 radici complesse (Φ_3)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASE 1 — GRÖBNER F4: wiring completo
          Effort: 4 giorni
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

STATO CONFERMATO DAL AUDIT:
  f4_groebner:  IMPLEMENTATA (righe 135-239), libera, NON in nessun header.
  to_f4:        IMPLEMENTATA (riga 243, static), MultivariatePolynomial→PolyF4.
  f4_to_expr:   IMPLEMENTATA (riga 257, static [[maybe_unused]]).
  PolyF4:       struct definita solo in polynomial_groebner_f4.cpp, non esportata.
  MonomialOrder: NON ESISTE nel codebase.
  solve_nonlinear_system_f4: linee 294-327 HARDCODED per Test 8.
  polynomial_groebner.cpp: 27 righe, restituisce {}.
  MultivariatePolynomial:   Unimplemented per coefficienti razionali (riga 37).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Sotto-task 1.0 — Crea header polynomial_groebner_f4.hpp (0.5 d)
  File nuovo: src/algebra/polynomial_groebner_f4.hpp

  Sposta FUORI dal .cpp le definizioni di:
    using Monomial = std::vector<unsigned int>;
    struct MonomialLexComparator { ... };
    struct MonomialGRevLexComparator { ... }; ← NUOVO (vedi sotto)
    enum class MonomialOrder { Lex, GRevLex };
    struct PolyF4 { ... };
    struct Pair { ... };
    class MacaulayMatrix { ... };

  Forward-declare le funzioni libere:
    std::vector<PolyF4> f4_groebner(std::vector<PolyF4> G,
                                    MonomialOrder order = MonomialOrder::GRevLex);
    Result<ExprPtr> f4_to_expr(const PolyF4& p, const std::vector<Symbol>& vars,
                               symbolic::CASContext& ctx);
    Result<PolyF4> expr_to_f4(ExprPtr expr, const std::vector<Symbol>& vars,
                              symbolic::CASContext& ctx);  // NUOVA

  IMPORTANTE: polynomial_groebner_f4.cpp deve includere questo header.
  Non rimuovere le definizioni — solo aggiungere l'include nel .cpp.

  MonomialGRevLexComparator (graded reverse lex):
    bool operator()(Monomial a, Monomial b):
      uint da = sum(a), db = sum(b)
      if da != db: return da > db   // grado totale decrescente
      for k = n-1 downto 0:        // lex inversa sulle variabili
        if a[k] != b[k]: return a[k] < b[k]
      return false
  Nota: GRevLex è più efficiente di Lex per F4 perché produce catene di riduzione
  più corte. Obbligatorio per sistemi di grado elevato (Cyclic-3, ecc.).

Sotto-task 1.1 — Aggiorna f4_groebner per MonomialOrder parametrico (0.5 d)
  File: src/algebra/polynomial_groebner_f4.cpp

  Modifica firma:
    std::vector<PolyF4> f4_groebner(std::vector<PolyF4> G,
                                    MonomialOrder order = MonomialOrder::GRevLex)

  MacaulayMatrix deve usare il comparatore corretto:
    if (order == MonomialOrder::GRevLex):
      usa MonomialGRevLexComparator
    else:
      usa MonomialLexComparator

  Il comparatore va passato/propagato alla struttura interna della MacaulayMatrix.
  Opzione: template su Comparator, oppure std::function<bool(Monomial,Monomial)>.
  Scegliere std::function per semplicità (nessun template extra).

  Aggiungere Buchberger product criterion (riduce S-polinomi inutili):
    In f4_groebner, nella inizializzazione delle coppie (riga ~141-146):
    Prima di aggiungere la coppia (i,j), verificare:
      bool coprime = true;
      for (size_t k = 0; k < num_vars; k++)
        if (G[i].leading_monomial()[k] > 0 && G[j].leading_monomial()[k] > 0)
          { coprime = false; break; }
      if (coprime) continue; // S-poly riduce a 0 per criterio del prodotto

Sotto-task 1.2 — Implementa expr_to_f4 (0.5 d)
  Funzione nuova: Result<PolyF4> expr_to_f4(ExprPtr, const vector<Symbol>&, CASContext&)
  File: src/algebra/polynomial_groebner_f4.cpp + header

  Algoritmo ricorsivo su ExprPtr:
    IntegerLit(n): PolyF4 con termine {monomial_zero → Rational(n,1)}
    RationalLit(p,q): PolyF4 con termine {monomial_zero → Rational(p,q)}
    Symbol(s): se s in vars[k]: PolyF4 con termine {e_k → Rational(1,1)}
               altrimenti: tratta come costante parametrica (errore Unimplemented)
    Unary(Neg, e): negate PolyF4 coefficients
    Binary(Add, a, b): somma PolyF4 termine per termine
    Binary(Sub, a, b): differenza PolyF4
    Binary(Mul, a, b): moltiplicazione PolyF4 (convoluzione monomiali)
    Binary(Pow, base, IntegerLit(n)): potenza PolyF4 via binary exponentiation
    Sum(terms): fold con Add
    Product(factors): fold con Mul
    DecimalLit: decimal_to_rational poi RationalLit

  Moltiplicazione PolyF4 (necessaria):
    Result<PolyF4> poly_f4_multiply(const PolyF4& a, const PolyF4& b):
      Per ogni (ma, ca) in a.terms, (mb, cb) in b.terms:
        monomial_prod[k] = ma[k] + mb[k]
        coeff_prod = ca * cb
        Accumula in result.terms[monomial_prod] += coeff_prod

  NOTA: to_f4 esistente usa MultivariatePolynomial (solo interi). expr_to_f4
  bypassa MultivariatePolynomial e lavora direttamente su ExprPtr. NON rimuovere
  to_f4 (potrebbe essere usata da solve_nonlinear_system_f4).

Sotto-task 1.3 — Aggiungere inter-riduzione post-F4 (0.5 d)
  File: src/algebra/polynomial_groebner_f4.cpp

  Funzione: void inter_reduce(vector<PolyF4>& G, MonomialOrder order)
    Per ogni i:
      Riduci G[i] rispetto a {G[j] : j ≠ i} via riduzione standard:
        Finché esiste monomial m in G[i] divisibile per LM(G[j]):
          quotiente = coefficiente di m / LC(G[j])
          moltiplicatore = quotiente * G[j] * (m / LM(G[j]))
          G[i] -= moltiplicatore
      Rendi G[i] monico (G[i].make_monic())

  Chiama inter_reduce dopo f4_groebner nel wrapper polynomial_groebner.cpp.
  La base ridotta è unica → risultato deterministico.

  Funzione utility: bool f4_poly_divides(Monomial lm, const PolyF4& g)
    Verifica se LM(g) divide lm termine per termine.

  Funzione utility: PolyF4 f4_normal_form(PolyF4 f, const vector<PolyF4>& G, MonomialOrder)
    Riduzione di f modulo G (per ideal membership: f ∈ (G) ↔ NF(f,G) = 0).

Sotto-task 1.4 — Rewrite polynomial_groebner.cpp (0.5 d)
  File: src/algebra/polynomial_groebner.cpp

  Sostituire completamente il corpo:

    #include "polynomial_groebner_f4.hpp"

    Result<vector<ExprPtr>> polynomial_groebner(
        const vector<ExprPtr>& equations,
        const vector<Symbol>& variables,
        CASContext& ctx) {

      if (equations.empty()) return ok(vector<ExprPtr>{});

      // Converti a PolyF4 usando coefficienti razionali
      vector<PolyF4> F;
      F.reserve(equations.size());
      for (ExprPtr eq : equations) {
        auto r = expr_to_f4(eq, variables, ctx);
        if (r.is_error()) return fail<vector<ExprPtr>>(r.error());
        if (r.value().terms.empty()) continue; // polinomio zero
        F.push_back(r.value());
      }

      if (F.empty()) return ok(vector<ExprPtr>{});

      // Calcola base di Gröbner (GRevLex per efficienza)
      auto G = f4_groebner(F, MonomialOrder::GRevLex);

      // Inter-riduzione → base ridotta unica
      inter_reduce(G, MonomialOrder::GRevLex);

      // Converti tornando a ExprPtr
      vector<ExprPtr> result;
      result.reserve(G.size());
      for (const PolyF4& g : G) {
        auto e = f4_to_expr(g, variables, ctx);
        if (e.is_error()) return fail<vector<ExprPtr>>(e.error());
        result.push_back(e.value());
      }
      return ok(result);
    }

Sotto-task 1.5 — Riscrivere solve_nonlinear_system_f4 (1.5 d)
  File: src/algebra/polynomial_groebner_f4.cpp, righe 278-327

  CANCELLARE completamente le righe 294-327 (hardcode Test 8).

  Nuovo algoritmo — Shape Lemma back-substitution:
  (funziona quando l'ideale è 0-dimensionale e la GB è in lex order)

  Struttura della nuova funzione:
    1. Converti equazioni in PolyF4 via expr_to_f4
    2. Calcola GB in Lex order (non GRevLex, perché shape lemma richiede Lex)
       G_lex = f4_groebner(F, MonomialOrder::Lex)
       inter_reduce(G_lex, MonomialOrder::Lex)
    3. Verifica 0-dimensionalità: conta monomiali della "staircase"
       (monomiali non divisibili da alcun LM(g)) — deve essere finito
    4. Trova il "pure" polynomial: unico g in G_lex che dipende solo da var[n-1]
       Identificazione: tutti i monomiali di g hanno exp=0 per var[0..n-2]
    5. Solve pure polynomial: solve_polynomial(f4_to_expr(g_pure, vars, ctx), vars[n-1], ctx)
    6. Per ogni radice r_{n-1}:
       a. Sostituisci var[n-1] = r_{n-1} in tutti i restanti g ∈ G_lex
       b. Ottieni sistema triangolare in var[0..n-2]
       c. Ricorri su sistema ridotto
    7. Build solution tuples: [[x_0=v0, x_1=v1, ...]] per ogni soluzione

  NOTA: la back-substitution ricorsiva è O(n²) ed è sufficiente per n ≤ 10.
  Per n > 10 usare FGLM (Fase 4).

  Gestione fallimenti:
    Se il sistema non è 0-dimensionale: restituisci Unimplemented
      "Sistema non 0-dimensionale: infinità di soluzioni o ideale non primo"
    Se pure polynomial non trovato: restituisci Unimplemented
      "Base di Gröbner non in forma shape lemma — sistema probabilmente non 0-dim"

Sotto-task 1.6 — Fix MultivariatePolynomial per Q (0.5 d)
  File: src/algebra/polynomial_multivariate.cpp:27-37

  Sostituire il blocco RationalLit:
    ATTUALE (riga 27-38):
      if (const auto* rat = expr_cast<RationalLit>(expr)) {
        if (rat->denominator == BigInt(1)) { ... return ok con coeff intero ... }
        return fail(Unimplemented, "supporta solo coefficienti interi");
      }

    NUOVO:
      if (const auto* rat = expr_cast<RationalLit>(expr)) {
        // Rappresenta p/q come MultivariateTerm con coeff = numeratore,
        // poi alla fine della parse scala tutti i coefficienti.
        // ALTERNATIVA più semplice: rifiuta MultivariatePolynomial per Q e usa
        // expr_to_f4 invece (che già supporta Q). In questo caso aggiungi un
        // commento: "usare expr_to_f4 per sistemi con coefficienti razionali".
        // IMPLEMENTAZIONE RACCOMANDATA: converti RationalLit a intero scalato.
        // Accumula tutti i denominatori, alla fine moltiplica il numeratore.
        // Per ora: restituisci coeff = numeratore, denominatore separato.
        // MultivariateTerm deve avere un campo Rational coefficient.
      }

  NOTA: se MultivariatePolynomial usa BigInt per coefficienti e cambiare il tipo
  introduce troppi refactoring, NON cambiare MultivariatePolynomial. Invece:
  expr_to_f4 (Sotto-task 1.2) gestisce già coefficienti razionali. Il Gröbner
  wrapper usa expr_to_f4 direttamente, non MultivariatePolynomial.
  La riga 37 rimane Unimplemented per MultivariatePolynomial ma è irrilevante
  per il path Gröbner.

Test per Fase 1:
  {x²+y²-1, x-y} → GB = {2y²-1, x-y}   ← Test base
  {x+y-2, x-y} → {x-1, y-1}              ← Sistema lineare
  Test 19 ACID deve passare
  Cyclic-3: {x+y+z, xy+yz+xz, xyz-1} → termine entro 5 secondi

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASE 4 — csolve N-D via Gröbner
          Effort: 2 giorni (richiede Fase 1 completata)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

STATO CONFERMATO:
  csolve.cpp:113: Unimplemented per N ≠ 2.
  Path 2×2 via resultant polynomial_resultant() funziona correttamente.
  polynomial_resultant: IMPLEMENTATA con subresultant PRS.

Algoritmo:
  File: src/algebra/csolve.cpp — rimuovi il blocco di errore a riga 113
  File nuovo: src/algebra/csolve_groebner.cpp (~200 LoC)

  Funzione: Result<ExprPtr> csolve_nd(
      const vector<ExprPtr>& equations,
      const vector<Symbol>& variables,
      CASContext& ctx)

  Algoritmo:
    1. Se equations.size() == 2 e variables.size() == 2: delega al path resultante
       esistente (fast-path, invariato).
    2. Chiama polynomial_groebner(equations, variables, ctx) → G_lex in lex order
       [polynomial_groebner internamente chiama f4_groebner con GRevLex poi converte]
       ATTENZIONE: per shape lemma abbiamo bisogno di Lex. Aggiungere overload:
         polynomial_groebner_lex(equations, variables, ctx) che usa Lex direttamente.
    3. Esegui shape lemma back-substitution (stessa logica di solve_nonlinear_system_f4
       dopo il refactoring). Riutilizza la funzione estratta.
    4. Costruisci la matrice soluzione: ogni riga = una soluzione,
       ogni colonna = una variabile. Formato identico al path 2×2 esistente.

  Integrazione in csolve.cpp:
    Rimuovi: return fail(Unimplemented, "csolve supporta solo sistemi 2x2")
    Aggiungi: return csolve_nd(equations_vec, variables_vec, ctx)

  Test:
    sistema 3×3 lineare → soluzione unica razionale esatta
    sistema 2×2 non-lineare di regressione (uguale a prima) → stesso output
    {x²+y=2, x+y²=2} → 4 soluzioni

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASE 6 — ODE SIMBOLICHE: classificatore + coeff. costanti
          Effort: 4 giorni
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

STATO CONFERMATO:
  ode_classifier.cpp:17: restituisce SEMPRE OdeType::Unknown.
  ode_solver_advanced.cpp:14: Unimplemented per ConstantCoeff.
  ode_solver_advanced.cpp:23: Unimplemented per Kovacic.
  OdeClassification struct: {type, equation, y, x, components: vector<ExprPtr>}
  OdeType enum: Unknown, Separable, Linear1stOrder, Bernoulli, Exact,
                Linear2ndOrderConstantCoeff, Linear2ndOrderRationalCoeff

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Sub-task 6.0 — Implementa classify_ode per 2° ordine (1 giorno)
  File: src/calculus/ode_classifier.cpp (attuale: 42 righe, sempre Unknown)

  Algoritmo di classificazione 2° ordine costante:
  Dato l'equazione come espressione, cercare la forma:
    a₂(x)·y''(x) + a₁(x)·y'(x) + a₀(x)·y(x) = f(x)
  dove a₂, a₁, a₀ sono costanti (non dipendono da x) e y è la funzione incognita.

  Step 1 — Riordina equazione:
    Normalizza: sposta tutto a sinistra (equazione = 0)
    Usa `algebra::together(eq, ctx)` se necessario

  Step 2 — Identifica derivate:
    Cerca nella AST: FuncCall(Derivative, [y, x, IntegerLit(2)]) = y''
    Cerca: FuncCall(Derivative, [y, x, IntegerLit(1)]) = y'
    Cerca: Symbol(y.name) = y

  Step 3 — Estrai coefficienti:
    collect_coefficient(expr, target):
      Somma tutti i termini Product([c, target]) → restituisce c
      Se target appare senza coefficiente → aggiungi 1

    a₂ = collect_coefficient(eq, y'')
    a₁ = collect_coefficient(eq, y')
    a₀ = collect_coefficient(eq, y)
    f  = resto di eq (termini senza y, y', y'')

  Step 4 — Verifica costanza:
    Funzione: bool is_constant_wrt(ExprPtr expr, const Symbol& x)
      Calcola diff(expr, x, ctx). Se risultato = 0 → costante.
      ATTENZIONE: diff su costante deve restituire 0, non errore.

    Se a₂ ≠ 0 E is_constant_wrt(a₂, x) E is_constant_wrt(a₁, x) E
             is_constant_wrt(a₀, x):
      return OdeClassification{
        type: Linear2ndOrderConstantCoeff,
        equation: original,
        y: y, x: x,
        components: {a₂, a₁, a₀, f}  // ordine esatto
      }

  Step 5 — Fallback:
    Verifica Linear1stOrder (solo y' e y, costanti o razionali in x)
    Verifica Separable (f(y)·g(x) = 0 forma)
    Altrimenti: Unknown

  NOTA: ode_solver_1st_order.cpp:69 ha già un solver per 1° ordine.
  La classify_ode per 1° ordine era presumably già implementata ma non funzionante
  (dato che restituiva Unknown). Fixare anche per completezza.

Sub-task 6.1 — Implementa solve_ode_advanced per coeff. costanti (2 giorni)
  File: src/calculus/ode_solver_advanced.cpp, riga 12-14 (stub)

  Estrai da components: [a₂, a₁, a₀, f]

  PARTE A — Soluzione omogenea (f=0):

    Step 1: Polinomio caratteristico
      λ² * a₂ + λ * a₁ + a₀ = 0
      p(λ) = a₂λ² + a₁λ + a₀
      Costruisci come IntPoly: {a₀, a₁, a₂} (coeffs little-endian)
      Oppure usa la formula della radice: λ = (-a₁ ± √(a₁²-4a₂a₀)) / (2a₂)

    Step 2: Calcola discriminante D = a₁² - 4a₂a₀
      Valuta D come Rational (se a₂, a₁, a₀ sono Rational)
      Oppure chiama context_.simplify(a₁²-4a₂a₀) e verifica tipo

    Step 3: Classifica radici
      Case D > 0 (discriminante positivo, radici reali distinte):
        λ₁ = (-a₁ + sqrt(D)) / (2a₂)
        λ₂ = (-a₁ - sqrt(D)) / (2a₂)
        Semplifica: per D = n² perfetto → λ₁, λ₂ razionali
        y_h = C₁·e^(λ₁x) + C₂·e^(λ₂x)

      Case D = 0 (radice doppia):
        λ = -a₁ / (2a₂)
        y_h = (C₁ + C₂·x)·e^(λx)

      Case D < 0 (coppia coniugata complessa):
        α = -a₁ / (2a₂)
        β = sqrt(-D) / (2a₂)
        y_h = e^(αx) · (C₁·cos(βx) + C₂·sin(βx))

    Step 4: Costruisci y_h come ExprPtr
      C₁, C₂: Symbol("_C1"), Symbol("_C2") — costanti di integrazione
      Costruisci somma usando arena_.make<>()

  PARTE B — Soluzione particolare (f ≠ 0) via variazione dei parametri:

    Data y_h con base {y₁, y₂}:
      Wronskiano W = y₁·y₂' - y₂·y₁'
        W è calcolato via diff(y₁, x) e diff(y₂, x) con context_.diff()
      Soluzioni particolari u₁, u₂:
        u₁ = integrate(-y₂·f/a₂/W, x, ctx)
        u₂ = integrate(y₁·f/a₂/W, x, ctx)
      y_p = y₁·u₁ + y₂·u₂

    Se integrate fallisce: restituisci solo y_h (manca soluzione particolare)
    con WARN nel trace "variazione parametri: integrale non risolto"

    ATTENZIONE: f può essere 0. In tal caso y_p = 0, skip.

  PARTE C — Soluzione completa:
    y = y_h + y_p
    Semplifica con context_.simplify()

  Test:
    y'' - y = 0             → C₁·eˣ + C₂·e^{-x}
    y'' + 4y = 0            → C₁·cos(2x) + C₂·sin(2x)
    y'' - 3y' + 2y = 0      → C₁·eˣ + C₂·e^{2x}
    y'' + y = sin(x)        → y_h + y_p (risonanza: x·cos/sin)
    y'' + 2y' + y = 0       → (C₁ + C₂x)e^{-x} (radice doppia)

Sub-task 6.2 — Kovacic L=1 (3 giorni, separato in Fase 6b)
  Vedere sezione Fase 6b più avanti nel documento.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASE 7 — SOMMAZIONE: Gosper + Zeilberger
          Effort: 5 giorni
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

STATO CONFERMATO:
  summation.cpp:63: Unimplemented per tutto tranne Basel hardcode.
  Nessuna infrastruttura ipergeometrica esiste.
  Infrastruttura polinomiale disponibile: poly_multiply, poly_divide_by_scalar,
  divide_poly_with_remainder, poly_extended_gcd (ora implementato), polynomial_gcd.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Sub-task 7.0 — Riconoscimento termine ipergeometrico (0.5 d)
  File nuovo: src/calculus/sum_hypergeometric.cpp

  Funzione: Result<RatPoly, RatPoly> hypergeometric_ratio(
      ExprPtr t_k, const Symbol& k, CASContext& ctx)
    Calcola r(k) = t(k+1)/t(k) come funzione razionale in k.
    Steps:
      1. Sostituisci k → k+1 in t_k: t_kp1 = substitute(t_k, k, k+1)
      2. Forma t_kp1/t_k
      3. Chiama algebra::together() per razionalizzare
      4. Estrai numeratore e denominatore come RatPoly (polinomi in k)
      5. Verifica che il risultato non contenga k in modo trascendente
         (es. k^k o e^k non ipergeometrico)
    Restituisce (numeratore, denominatore) della funzione razionale r(k).

  Funzione: bool is_hypergeometric(ExprPtr t_k, const Symbol& k, CASContext& ctx)
    Chiama hypergeometric_ratio e verifica che non sia error.

  Casi ipergeometrici comuni:
    k!     → ratio = k+1          ✓
    C(n,k) → ratio = (n-k)/(k+1)  ✓ (n parametro)
    a^k    → ratio = a             ✓
    1/k!   → ratio = 1/(k+1)      ✓
    k      → ratio = (k+1)/k      ✓

Sub-task 7.1 — Algoritmo di Gosper (2 d)
  File nuovo: src/calculus/sum_gosper.cpp (~200 LoC)
  File nuovo: src/calculus/sum_gosper.hpp

  Funzione pubblica:
    Result<ExprPtr> gosper_antidifference(ExprPtr t_k, const Symbol& k, CASContext& ctx)
    Restituisce S tale che S(k+1) - S(k) = t_k, oppure error(Unimplemented).

  Algoritmo Gosper completo (Gosper 1978, Petkovšek–Wilf–Zeilberger §5):

  STEP 1 — Calcola r(k) = t(k+1)/t(k):
    auto [p, q] = hypergeometric_ratio(t_k, k, ctx)  // p,q RatPoly in k
    Dividi: p/q (già razionalizzato)

  STEP 2 — Gosper-Petkovšek factoring (trova a,b,c tale che r = a·c(k+1)/b/c(k)):
    Obiettivo: trovare polinomi a(k), b(k), c(k) tale che:
      r(k) = a(k)/b(k) · c(k+1)/c(k)
    con gcd(a(k), b(k+j)) = 1 per tutti j = 0, 1, 2, ...

    Algoritmo:
      a₀ = p, b₀ = q (estrai da r = p/q già ridotto)
      Cerca "shift content": trova il massimo c(k) polinomio tale che
        gcd(a₀(k), b₀(k+j)) = 1 per ogni j ≥ 0
      Metodo:
        Per ogni j = 0, 1, ..., deg(a₀)+deg(b₀):
          g_j = gcd(a₀(k), b₀(k+j)) [gcd come polinomi in k]
          Se g_j ≠ 1:
            a₀(k) = a₀(k) / g_j(k)       // cancella dal numeratore
            b₀(k) = b₀(k) / g_j(k-j)    // cancella dal denominatore (shift)
            c(k) *= g_j(k)                // accumula in c

        Risultato finale: a = a₀, b = b₀, c = c prodotto
        GCD polinomiale: usare polynomial_gcd (già implementato) convertendo
        RatPoly → IntPoly se coefficienti razionali o tenendo come PolyExpr.

  STEP 3 — Bound sul grado di x(k):
    deg_x = max(deg(a) - deg(b), 0)
    Se deg(a) == deg(b) e leading_coeff(a) == leading_coeff(b): deg_x = 0
    Se deg(a) < deg(b): NO soluzione polinomiale → return Unimplemented

  STEP 4 — Equazione funzionale a(k)·x(k+1) - b(k-1)·x(k) = c(k):
    Sostituisci x(k) = Σ_{i=0}^{deg_x} d_i · k^i  (coefficienti unknown d_i)
    Sostituisci nell'equazione, espandi, raggruppa per potenze di k
    Ottieni sistema lineare: M · [d_0,...,d_{deg_x}]^T = [rhs coefficienti]

    Soluzione sistema lineare: usa linalg esistente oppure
    Gaussian elimination su matrice di Rational coefficienti.

    Se sistema ha soluzione unica o parametrica (d_i ∈ Q): x(k) è trovato.
    Se sistema inconsistente: return error(Unimplemented, "Gosper: no antidifference")

  STEP 5 — Antidifferenza:
    S(k) = b(k-1) / c(k) · x(k) · t_k
    Semplifica e ritorna.

  Verifica: S(k+1) - S(k) deve essere = t_k. Da asserire in debug mode.

  Test di Gosper:
    Σ k → n(n+1)/2  [t_k = k, ratio = (k+1)/k]
    Σ 2^k → 2^{n+1}-1  [t_k = 2^k, ratio = 2]
    Σ k·2^k → (n-1)·2^{n+1}+2  [più complesso]
    Σ C(n,k) = 2^n  [Zeilberger, non Gosper puro]

Sub-task 7.2 — Algoritmo di Zeilberger (2 d)
  File nuovo: src/calculus/sum_zeilberger.cpp (~200 LoC)
  File nuovo: src/calculus/sum_zeilberger.hpp

  Funzione pubblica:
    Result<ExprPtr> zeilberger_definite_sum(
        ExprPtr F_nk, const Symbol& n, const Symbol& k,
        ExprPtr lower, ExprPtr upper, CASContext& ctx)
    Restituisce f(n) = Σ_k F(n,k) se trovata la ricorrenza, altrimenti Unimplemented.

  Algoritmo (Zeilberger 1990, Petkovšek–Wilf–Zeilberger §6):

  STEP 1 — Verifica che F(n,k) sia ipergeometrica propria in k:
    r_k(n,k) = F(n,k+1)/F(n,k) deve essere razionale in n,k.
    Chiamare hypergeometric_ratio con la variabile k.

  STEP 2 — Cerca ricorrenza di ordine J (J = 1, 2, 3):
    Per J = 1, 2, 3:
      Cerca a_0(n), ..., a_J(n) polinomi in n tale che:
        Σ_{j=0}^J a_j(n) · F(n+j, k) = G(n, k+1) - G(n, k)
      dove G(n,k) = R(n,k) · F(n,k) con R razionale in n,k.

      Metodo — Gosper multiparametrico:
        Definiamo T(n,k) = Σ_j a_j(n) · F(n+j,k) [combinazione lineare con a_j unknown]
        Calcoliamo il ratio r(n,k) = T(n,k+1)/T(n,k)
        Applichiamo Gosper a T per trovare G tale che T = G(k+1) - G(k)
        Le condizioni a_j si ottengono imponendo che il sistema sia consistente

        In pratica: si fissa il grado in n degli a_j, si parametrizza
        a_j(n) = Σ c_{j,m} n^m, si sostituisce in T, si applica Gosper,
        si risolve il sistema lineare in {c_{j,m}}.

  STEP 3 — Se trovata ricorrenza:
    Restituisce f(n) come soluzione della ricorrenza.
    Per ricorrenze del 1° ordine: f(n) = const · prodotto di fattori
    Per ricorrenze di ordine superiore: usa Petkovšek Hyper (Sub-task 7.3)

  STEP 4 — Somme definite:
    Se lower e upper sono concreti (es. 0 e n):
      Applica Teorema Fondamentale del Calcolo Discreto:
        Σ_{k=lower}^{upper} t_k = S(upper+1) - S(lower)
        dove S è l'antidifferenza di Gosper

  NOTA CRITICA: Zeilberger richiede la composizione di Gosper con parametri
  aggiuntivi. Il grado di difficoltà è ALTO. Si raccomanda di implementare
  prima casi speciali e poi generalizzare:

    FASE 7.2a — Somme polinomiali: Σ_k k^m = Bernoulli polynomial
      Se t_k = k^m: usa formula di Bernoulli B_m(n+1)/m+1
      Precomputa B_m per m ≤ 10 (tabella statica).

    FASE 7.2b — Somme geometriche: Σ_k a^k·P(k) via Gosper
    FASE 7.2c — Identità ipergeometriche: Zeilberger completo

Sub-task 7.3 — Integrazione in summation.cpp (0.5 d)
  File: src/calculus/summation.cpp

  Sostituire il blocco Unimplemented a riga 63 con:

    // 1. Gosper (somma indefinita / antidifferenza)
    if (is_indefinite_sum) {  // lower = k, upper = symbol (non numerico)
        auto res = gosper_antidifference(term, var, ctx);
        if (res.is_ok()) return res;
    }

    // 2. Gosper per somme definite via TFC
    auto anti = gosper_antidifference(term, var, ctx);
    if (anti.is_ok()) {
        auto at_upper = substitute(anti, var, upper+1)
        auto at_lower = substitute(anti, var, lower)
        return simplify(at_upper - at_lower)
    }

    // 3. Riconoscimento polinomiale Σk^m = Bernoulli
    if (is_polynomial_in_var(term, var)) {
        return sum_polynomial_formula(term, var, lower, upper, ctx);
    }

    // 4. Zeilberger per somme dipendenti da parametro esterno
    // (se term contiene simboli diversi da var)
    // return zeilberger_definite_sum(term, n_param, var, lower, upper, ctx);

    return fail(Unimplemented, "Somma non ipergeometrica o non polinomiale")

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASE 2c — RING SERIE TRONCATO K[[t]]/(t^N)
           Effort: 4 giorni
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

MOTIVO: limit_series.cpp usa differenziazione ripetuta che fallisce per
funzioni composte (es. exp(exp(x)-1)). Ring troncato è più robusto e veloce.
PREREQUISITO per: ODE Frobenius (Fase 6b), asintotici complessi.

Struttura dati:
  include/cas/series.hpp:

    struct TruncatedSeries {
        Symbol var;
        ExprPtr center;     // punto di espansione (es. IntegerLit(0))
        int order;          // N: troncamento a t^N
        int valuation;      // esponente minimo (Laurent: può essere negativo)
        int puiseux_denom;  // q per Puiseux t^(k/q), di default 1
        std::vector<ExprPtr> coeffs;  // coeffs[i] = coeff di (t-center)^(i+valuation)
        // Invariante: coeffs.size() <= order - valuation
    };

Operazioni del ring (src/calculus/series_core.cpp):
  Result<TruncatedSeries> series_add(a, b)
  Result<TruncatedSeries> series_mul(a, b)  // prodotto di Cauchy
  Result<TruncatedSeries> series_neg(a)
  Result<TruncatedSeries> series_inv(a)     // Newton: richiede a.coeffs[0] ≠ 0
  Result<TruncatedSeries> series_compose(f, g)  // Horner: g.valuation ≥ 1

  Newton inversion (a.coeffs[0] invertibile):
    b_0 = 1 / a_coeffs[0]
    b_{k} = -(1/a_0) · Σ_{j=1}^k a_j · b_{k-j}  (iterazione)

Funzioni trascendenti (src/calculus/series_transcendental.cpp):
  series_exp(a, N):  se a.coeffs[0] == 0 (necessario!)
    Usa ODE: B'=A'·B, B(0)=1 → B[k] = (1/k) · Σ_{j=1}^k j·A[j]·B[k-j]
  series_log(a, N):  a.coeffs[0] = 1 (log(1+...)  normalizzato)
    Integra A'/A: C[k] = (A[k] - Σ_{j=1}^{k-1} C[j]·A[k-j]) / A[0]
  series_pow(a, alpha, N): exp(alpha·log(a)) usando series_exp + series_log
  series_sin, series_cos via Taylor table o ODE

Expansion da AST (src/calculus/series_expand.cpp):
  Result<TruncatedSeries> series_expand(ExprPtr expr, Symbol t, ExprPtr center, int order, CASContext& ctx)

  Ricorsione:
    IntegerLit, RationalLit, Constant, Symbol ≠ t: serie costante
    Symbol == t: {valuation=1, coeffs=[1]}
    Binary(Add/Sub/Mul): operazioni ring
    Binary(Pow, e, IntegerLit(n)):
      se e = t: serie monomiale
      altrimenti: series_pow con alpha=n
    FuncCall(Exp, [a]): series_expand(a) poi series_exp
    FuncCall(Ln, [a]): series_expand(a) poi normalize poi series_log
    FuncCall(Sin/Cos/...): via composizione

  Integrazione con limit_series.cpp:
    Nuovo path in taylor_series(): se differenziazione ripetuta fallisce,
    prova series_expand come fallback.

  Test:
    series(1/(1-x), x, 0, 6) → 1+x+x²+x³+x⁴+x⁵
    series(tan(x), x, 0, 8) → x+x³/3+2x⁵/15+17x⁷/315
    series(exp(exp(x)-1), x, 0, 4) → 1+x+x²+x³·5/6+...

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASE 6b — KOVACIC L=1 (ODE 2° ordine razionale)
           Effort: 5 giorni (dopo Fase 6 e Fase 2c)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Problema: y'' = r(x)·y dove r ∈ Q(x).
Trasformazione: y = exp(∫ω) → ω' + ω² = r (equazione di Riccati).

File nuovo: src/calculus/ode_kovacic.cpp (~600 LoC)
Riferimento: Kovacic J.J., "An algorithm for solving second order linear
  homogeneous differential equations", J. Symbolic Computation 2 (1986), 3-43.

Caso L=1 (soluzione Liouvilliana: ω ∈ Q(x)):
  Step 1 — Analisi dei poli di r:
    Sia r = p/q, radici di q = poli di r.
    Per ogni polo c di r di ordine n, calcola la parte principale di ω in c:
      Se n=2: coeff leading di (1/2 ± √(1/4+a_2)) dove a_2 coeff quadratico
    Calcola le quantità E_c per ogni polo (finite set di candidate).

  Step 2 — Bound sul grado del componente polinomiale di ω:
    ω = p(x)/q(x) + parte "intera" di grado d
    d = max(0, max_polo - grado_polo) secondo tabella Kovacic

  Step 3 — Costruisci ω con coefficienti unknown:
    ω = A(x) + Σ (parti principali da poli)
    dove A(x) = polinomio di grado d con coefficienti {a_0,...,a_d} unknown

  Step 4 — Sostituisci in ω' + ω² = r:
    Espandi, raccogli per potenze di x.
    Sistema lineare in {a_i} → solve via linalg esistente.

  Step 5 — Se soluzione: y = exp(∫ω dx)
    integrate(ω, x, ctx) poi exp(...)

  PREREQUISITI per Fase 6b:
    integrate() funziona per funzioni razionali (Hermite ok dopo poly_xgcd fix ✓)
    factorization_polynomials per trovare radici di q(x) (esiste ✓)
    linalg (Matrix) per sistema lineare (esiste ✓)

  Test:
    y'' - y = 0         → exp(±x)   [r = 1, ω = ±1]
    y'' - x²y = 0       → Airy (no Liouvilliana)
    2x²y'' + 3xy' - y = 0 → Cauchy-Euler (ω razionale)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
INFRASTRUTTURA CROSS-CUTTING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CI e Test:
  Dopo ogni fase: esegui l'intera suite ACID e tutte le regressioni.
  Aggiungere test specifici per ogni fase nuova.

  Suite Gröbner: test/unit/algebra/test_groebner.cpp (da creare)
    GB di sistemi lineari, quadratici, Cyclic-3
  Suite ODE: test/unit/calculus/test_ode.cpp (da creare)
    ODE coeff costanti, omogenee e non
  Suite Summation: test/unit/calculus/test_summation.cpp (da creare)
    Gosper casi base, Zeilberger binomiale

Benchmark gate:
  Ogni fase: eseguire scripts/benchmark.sh e verificare nessuna regressione
  rispetto a baseline_release.txt.

Contatore Unimplemented:
  Prima di ogni fase: cat $(find src -name "*.cpp") | grep -c "Unimplemented"
  Deve diminuire dopo ogni fase completata.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
TIMELINE REALISTICA
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Giorno  0:     Fase 2b — Fix Test 28 (6 righe, 0.5d)
Giorno  0-1:   Fase 9  — Fix Test 9 (solve_by_factoring + ciclotomico)
Giorno  1-5:   Fase 1  — Gröbner F4 wiring (4d, 6 sub-task)
Giorno  5-7:   Fase 4  — csolve N-D (2d, richiede Fase 1)
Giorno  7-11:  Fase 6  — ODE classifier + const-coeff (4d)
Giorno  11-15: Fase 2c — Ring serie troncato (4d)
Giorno  15-20: Fase 7  — Gosper + Zeilberger (5d)
Giorno  20-25: Fase 6b — Kovacic L=1 (5d, richiede Fase 6 e 2c)

Al completamento Fase 2b+9: 20/21 ACID test.
Al completamento Fase 1+4:  ~55-60% parità HP Prime G2 (GB sbloccato).
Al completamento tutte fasi: ~80-85% parità su operazioni comuni.

Gap permanente fuori scope di questo piano:
  Laplace/Fourier transforms, plotting 2D/3D (ImGui, Fase 10),
  statistica, finanza, Gröbner per ideali positivo-dimensionali,
  soluzioni in forma chiusa per ODE non-lineari.
