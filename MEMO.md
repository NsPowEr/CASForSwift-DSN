---                                                                         
  VERDETTO CRITICO — REAL_CAS_ENGINE_C++ @ 2026-04-24                                                                   
                                                                                                                        
  1. Falla ARCHITETTURALE #1 — INTEGRAZIONE E MODULARITÀ (RISOLTO STRUTTURALMENTE)
                                                                                                                        
  - [x] Split `integrate.cpp` e tutti i monoliti simbolici/algebrici (< 500 righe).
  - [x] Implementazione `u-substitution` reale.
  - [x] Infrastruttura per Hermite Reduction predisposta.
                                                                                                                        
  ---                                                                                                                   
  2. Falla ARCHITETTURALE #2 — ECCEZIONI NEL CORE (RISOLTO)
                                                                                                                        
  - [x] Riscritto `BigInt` e `Rational` con `Result<T>`.
  - [x] Rimossi tutti i `throw` dal core simbolico.
  - [x] Build pulita con -Wall -Wextra -Werror.
                                                                                                                        
  ---                                                                                                                   
  3. Falla ARCHITETTURALE #3 — ORDINAMENTO E TERMINAZIONE (IN CORSO)
                                                                                                                        
  - [x] Implementazione base LPO (Lexicographic Path Ordering).
  - [ ] Integrazione LPO in Knuth-Bendix per orientamento automatico regole.
                                                                                                                        
  ---                                                                                     
  7. Falla ARCHITETTURALE #7 — MANCA HASH-CONSING (RISOLTO)

  Rule #2 in CLAUDE.md: «Structural Sharing via Hash-consing».
  - [x] Implementata tabella interning in `AstArena`.
  - [x] Aggiunto supporto thread-safe tramite `std::mutex`.
  - [x] Ottimizzato hashing per `BigInt` (niente stringhe intermedie).
  - [x] Garantita identità puntatore per espressioni strutturalmente uguali.

  ---  4. Falla #4 — FATTORIZZAZIONE (PIANIFICATO)
  ... (resta uguale)
                                                                                                                        
  ---                                                                                                                   
  5. Falla #5 — BIGINT A LIMBS MA ALGORITMI SCHOOLBOOK (PIANIFICATO)
  ... (resta uguale)
