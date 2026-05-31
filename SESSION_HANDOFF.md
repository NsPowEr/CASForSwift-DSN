# SESSION HANDOFF — Resume Point per Prossima Sessione

> **Data**: 2026-05-21
> **Sessione precedente**: setup pre-F0 di `PLAN_HP_PRIME_PARITY.md` (parità HP Prime G2)
> **Orchestrator**: Opus 4.7

---

## STATO ATTUALE (cosa è già stato fatto)

### Infrastruttura golden oracle
- [x] Maxima **5.49.0** installato via Homebrew (`/opt/homebrew/bin/maxima`)
- [x] Smoke test passato: `integrate(sin(x)^2,x)` → output simbolico corretto
- [x] **CLAUDE.md Regola 6** aggiunta: divieto assoluto modifica sorgenti Maxima (GPL-2.0-only + integrità oracolo)
- [x] **Manifest SHA-256** generato: `scripts/maxima_5.49.0_manifest.sha256` (310 entries: 1 binario + 309 `.lisp` core)
- [x] **Script integrità**: `scripts/verify_maxima_integrity.sh` (verde su sistema corrente)
- [x] **Script rigenerazione manifest**: `scripts/regenerate_maxima_manifest.sh` (per re-bootstrap)
- [x] **PLAN pin aggiornato**: 5.47.0 → 5.49.0 (linee 30, 89)

### Plan analizzato
- [x] `PLAN_HP_PRIME_PARITY.md` letto e compreso (845 righe, 8 fasi F0-F7 + F8 research)
- [x] Inventory progetto fatto: HARDCODE_LEDGER.md sparse, CAS_TASKS.md 539 LOC, `.APROJECT_REFERENCES/algorithms/` mancante, scripts/benchmark.sh presente, test/{unit,regression,fuzz,bench} presenti

---

## PROSSIMA SESSIONE — RESUME POINT

**OBIETTIVO**: avviare **FASE 0** (sanitizzazione tracker + infrastruttura CI) come da PLAN sez. F0.1-F0.8.

**Tier routing**: F0 = lavoro meccanico → **T1-Sonnet** single subagent, BLOCCO completo (no spawn multipli).

### Checklist F0 da eseguire (1-2 giorni AI)

#### F0.1 — Retroclassificazione CAS_TASKS.md (12 voci ottimistiche → Parziale)
Voci da retro-classificare con evidenza file:riga:
- L1-01 Gruntz MRV (rank statico)
- L1-02 Risch (log-extension assente)
- L1-08 GCD multivariato (Brown stub)
- L1-17 Pivot Bareiss (magic 1000/500/400 in `matrix_ops.cpp:243`)
- L2-04 Smith (solo Z)
- L2-06 fsolve (tolerance hardcoded)
- L2-08 Polar/log (ln complesso parziale)
- L2-19 Equivalenza Risch
- L2-22 Residue (solo quad+biquad)
- L3-04 Funzioni speciali
- L3-06 Trager tower (solo 2 livelli)
- L3-18 Galois (deg ≤4)

#### F0.2 — HARDCODE_LEDGER.md: append 15 voci
Lista esatta in PLAN.md sez. F0.2.

#### F0.3 — Coverage baseline
- gcov/lcov in CMake (`-fprofile-arcs -ftest-coverage`)
- Run suite completa → `coverage/index.html`
- Documentare baseline per modulo

#### F0.4 — Property-based test framework
- `rapidcheck` come git submodule in `test/`
- 6 property minime:
  - `gcd(a,b)·lcm(a,b) ≡ a·b` su Z[x]
  - `D(integrate(f)) ≡ f` su 200 integrali
  - `f(roots(f)) ≡ 0`
  - `A·inv(A) ≡ I` su 5×5
  - `factor(p) → ∏ ≡ p`
  - `partial_fractions(f) → riassemblata ≡ f`

#### F0.5 — Golden test vs Maxima 5.49.0
- Corpus **2000 input** (~200 per area) — generare con script Python o bash
- Reference: `maxima --very-quiet --batch-string="..."`
- Confronto AST equivalenza (no toString)
- Pin: Bronstein "Symbolic Integration I" esercizi numerati come sub-corpus calculus
- **Integrare hook**: ogni golden run chiama `scripts/verify_maxima_integrity.sh` prima → fail se mismatch hash

#### F0.6 — Build infrastructure CI
- Sanitizer (ASan + UBSan + TSan) target `cas_tests_san`
- Benchmark gate: regression ≥10% blocca merge
- Anti-monolith: `scripts/check_file_size.sh` fail se `.cpp/.hpp` >500 righe
- Diff-coverage CI: gate ≥85% sui nuovi line per PR
- ccache + Ninja default
- Mutation testing setup (`mutest` o `mull`)

#### F0.7 — Doc-per-algorithm template
- `.APROJECT_REFERENCES/algorithms/_TEMPLATE.md` (sezioni: Algoritmo, Riferimento bibliografico, Dominio I/O, Complessità, Limiti, Esempi accettati/rifiutati, Test certificatore, Performance baseline)
- PR senza doc → blocco merge

#### F0.8 — Error diagnostic framework
- Estendere `CASErrorKind::Unimplemented` con payload `{module, function, input_shape_descriptor, reason_code, suggestion, ticket_id}`
- Helper `make_unimplemented(...)` obbligatorio
- Audit grep `return Unimplemented(` → assicurare payload strutturato

### Exit gate F0
- Tracker sanato (12 voci)
- Ledger completo (15 voci)
- Coverage baseline + dashboard
- ≥50 property tests verdi
- Golden suite 2000 input eseguibile, pass-rate baseline misurato
- CI sanitizer + bench + anti-monolith + diff-coverage attivi
- Template doc disponibile
- Diagnostic framework rolled out

---

## COMANDO PER RIPARTIRE (prossima sessione)

1. **`/model opus`** (se non già su Opus 4.7)
2. **Apri questo file**: `SESSION_HANDOFF.md`
3. **Comando utente**: `"esegui F0"` o `"riprendi handoff"`
4. **Orchestrator Opus**:
   - Legge `PLAN_HP_PRIME_PARITY.md` sez. F0
   - Spawna **1 subagent T1-Sonnet** con prompt completo per blocco F0
   - Attende report compatto
   - Audit T2-Sonnet-thinking pre-review interlock (cattura ~50% problemi)
   - Audit finale T3-Opus (gate Risolta/Parziale)
   - Presenta gate utente per approvazione F1

---

## VINCOLI PERMANENTI DA RICORDARE

### CLAUDE.md regole architetturali (NON NEGOZIABILI)
1. BigInt only (no `int64_t`/`double` nel core simbolico)
2. Structural sharing (ExprPtr identity O(1))
3. Memory Arena bump allocator
4. Moltiplicazione implicita parser
5. DecimalLit confine simbolico/numerico
6. **Maxima sorgenti immutabili** (GPL-2.0-only + integrità oracolo)

### REGOLA ZERO
Mai via facile. Mai shortcut. Hardcode-of-passage → ledger obbligatorio.

### Tier routing
- T1-Sonnet: meccanico (refactor, ledger, test boilerplate, algoritmi canonici well-documented)
- T2-Sonnet-thinking: invariant tracking (Brown GCD, Trager, ODE classifier, branch cuts)
- T3-Opus: research-grade + audit (Risch structure theorem, Gruntz §3.5, CAD, Galois ≥5, audit periodici)

### Anti-pattern subagent
- Mai 1-task-1-agent (overhead 5-10×). Batch 5-15 task per spawn.
- Max 3-5 subagent paralleli.
- Audit sequenziale dopo paralleli.

### Gates utente obbligatori
Dopo ogni fase F_n: report + audit + coverage delta + pass-rate + ledger updates + richiesta approvazione esplicita per F_{n+1}.

### Cross-phase regression policy
F_n scopre bug in F_{n-k} → STOP F_n + rollback + fix come PR separata + re-audit F_{n-k} + riprende F_n solo dopo re-certificazione.

### API stability freeze
A chiusura ogni fase: public API in `include/cas/*.hpp` congelata. Semver bump + deprecation per modifiche cross-fase.

---

## TIMELINE GLOBALE

| Fase | Effort AI | Mix routing |
|---|---|---|
| F0 sanitizzazione | 1-2 giorni | 100% T1-Sonnet |
| F1 L0 foundation | 2-3 settimane | 60% Sonnet / 40% Opus |
| F2 L1 poly univariati | 1-2 settimane | 70% / 30% |
| F3 L2 multivar + alg ext | 3-4 settimane | 40% / 60% |
| F4 L2 linalg | 5-7 giorni | 80% / 20% |
| F5 L2 calculus | 4-6 settimane | 30% / 70% |
| F6 L3 numerica+complex+units | 2-3 settimane | 50% / 50% |
| F7 plotting + acceptance | 1-2 settimane | 70% / 30% |

**Totale**: ~3-4 mesi seriale, ~2-3 mesi con parallelismo aggressivo.

**Budget token stimato**: ~6.5M total. Costo ~$300-500 con tier-routing.

---

## DOMANDE APERTE DA CHIUDERE A INIZIO SESSIONE

1. **rapidcheck**: git submodule o vendor in `third_party/`?
2. **Golden corpus 2000 input**: generato fresh da Maxima oppure parte da seed esistente?
3. **CI runner**: locale (Docker?) o GitHub Actions?
4. **Diff-coverage tool**: `diff-cover` (Python) o `codecov`?
5. **Mutation testing**: `mutest`, `mull`, o `mutate++`?

---

## FILE TOCCATI IN QUESTA SESSIONE

- `CLAUDE.md` — aggiunta Regola 6 (Maxima immutabile)
- `PLAN_HP_PRIME_PARITY.md` — pin 5.47.0 → 5.49.0 (linee 30, 89)
- `scripts/maxima_5.49.0_manifest.sha256` — **NEW** (310 hash)
- `scripts/verify_maxima_integrity.sh` — **NEW** (executable)
- `scripts/regenerate_maxima_manifest.sh` — **NEW** (executable)
- `SESSION_HANDOFF.md` — **NEW** (questo file)

Nessun file sorgente C++ modificato. Nessun commit creato (utente deciderà al resume).

---

## COMMIT SUGGERITO PRIMA DI CHIUDERE

Se l'utente vuole snapshottare lo stato pre-F0:

```
chore(infra): pin Maxima 5.49.0 as immutable oracle (CLAUDE.md Rule 6)

- Install Maxima 5.49.0 via Homebrew as primary golden test reference.
- Add CLAUDE.md Rule 6: Maxima sources are GPL-2.0-only and MUST NOT be
  modified (preserves oracle independence + avoids derivative-work copyleft).
- Generate SHA-256 manifest of Maxima binary + 309 core .lisp files
  (scripts/maxima_5.49.0_manifest.sha256).
- Add scripts/verify_maxima_integrity.sh (CI gate) and
  scripts/regenerate_maxima_manifest.sh (bootstrap utility).
- Update PLAN_HP_PRIME_PARITY.md pin 5.47.0 → 5.49.0.
- Add SESSION_HANDOFF.md as resume point for next session.

Prepares F0 (tracker sanitization + CI infrastructure).
```
