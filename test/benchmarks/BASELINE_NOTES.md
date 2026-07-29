# Note sul baseline benchmark

## 2026-07-16 — baseline rigenerato (metodologia mediana-di-5)

- **Perché**: il baseline 2026-06-01 (run singolo) era invalidato dalle
  condizioni di misura, NON dal codice. Certificato via A/B worktree: il
  codice al commit `5964be9` (quello del vecchio baseline), ricompilato e
  misurato 2026-07-16 nelle stesse condizioni di HEAD, dà
  `poly_gcd_subresultant` 1.862ms e `poly_expand_complex` 1.091ms — identici
  a HEAD (1.841/1.079). Zero regressione software su quelle metriche.
- **Regressione software residua reale** (A/B, condizioni omogenee):
  `simplify_basic` +21% (0.063→0.076) e `arena_alloc` +14% (0.794→0.907) —
  costo funzionale di A30 (ops-budget per-nodo) + A31 (side-conditions),
  deliberato, non bug. Se in futuro va ridotto: profilare quei due path.
- **Condizioni di questa misura**: mediana-di-5, ma con carico utente non
  nullo sulla macchina (load 1-min ~20-30 con 10 CPU; poco dopo la misura il
  load è salito a 33+ e le stesse metriche sono raddoppiate). Il baseline è
  quindi CONSERVATIVO (numeri gonfiati dal carico). **Da ririgenerare a
  macchina scarica** (`bash scripts/benchmark.sh --update-baseline`, load
  1-min < 5): produrrà numeri più bassi e un gate più severo.
- **Lezione di metodo**: prima di dichiarare una regressione dal gate,
  (1) guardare il warning sul load stampato dallo script, (2) in dubbio,
  A/B worktree — confronto relativo a condizioni omogenee, mai numeri
  assoluti attraverso settimane/toolchain.
