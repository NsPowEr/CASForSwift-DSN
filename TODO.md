# REAL CAS ENGINE — Execution TODO (HP-like Consolidation)

## Pilastri Architetturali (Roadmap Consolidata)

- [x] **P1: Unificazione Complessi Canonici**
  - [x] Semplificazione $I^n$ e $I^2 = -1$
  - [x] Raggruppamento $a + bI$ in `Sum`
  - [x] Valutazione `sqrt(-x) -> I * sqrt(x)`
  - [x] Audit e rimozione residui `Symbol("i")`

- [x] **P2: Canonicalizzazione Forte**
  - [x] Unificazione sistemi di ordinamento (LPO/KB)
  - [x] Ordinamento deterministico indipendente dai puntatori
  - [x] Appiattimento potenze annidate $(a^b)^c$
  - [x] Valutazione quadrati perfetti in `sqrt`

- [x] **P3: Rewrite System Chiuso**
  - [x] Audit regole `exp`/`log` per terminazione
  - [x] Regole core (sin^2+cos^2, parità) orientate
  - [x] Integrazione `context` in `RewriteProvider`

- [x] **P4: Assumptions Engine**
  - [x] Propagazione proprietà AST-wide (positive, nonnegative, real)
  - [x] Integrazione in `simplify` (sqrt(x^2), abs, ecc.)

- [x] **P5: Polynomial Core**
  - [x] GCD subresultante PRS (Brown/Collins) robusto
  - [x] Fattorizzazione square-free (Yun)
  - [x] Integrazione GCD in `simplify` (via RewriteProvider)

- [ ] **P6: solve() e RootOf**
  - [ ] Gradi 1-4 analitici esatti
  - [ ] Fallback `RootOf` onesto

- [ ] **P7: Limit Pipeline**
  - [ ] Pipeline: Sostituzione -> Taylor -> Crescita
  - [ ] Gestione poli e singolarità logaritmiche

- [ ] **P8: Integration Core (Risch)**
  - [ ] Hermite Reduction
  - [ ] Logarithmic part (Lazard-Rioboo-Trager)

- [x] **P9: Acid Test Suite**
  - [x] Creazione `test/unit/test_acid_complex_canonical.cpp`
  - [x] Copertura Pilastri 1-5

- [x] **P10: Roadmap Sincronizzata**
  - [x] Sostituzione TODO.md con i 10 pilastri
