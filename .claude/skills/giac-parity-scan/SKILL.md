---
name: giac-parity-scan
description: Misura il gap CAS vs Giac 2.0.0 (target di parità) sul corpus golden. Genera scoreboard PARITY_GIAC.md per area (giac% vs CAS%) e candidati task. I gap trovati vanno SEMPRE elaborati come task in TASKLIST_MASTER.md (single source of truth) — lo scoreboard è solo misura. Giac fork/exec only (GPL-3.0-or-later, CLAUDE.md §6).
---

# giac-parity-scan

Misura oggettiva della distanza dal target di parità (Giac 2.0.0) sul corpus
golden condiviso. Produce `PARITY_GIAC.md` (scoreboard per area) e una lista di
candidati task.

## Regola non negoziabile — single source of truth

`PARITY_GIAC.md` è un **artefatto di misura rigenerabile**, NON un tracker.
Ogni gap che merita lavoro va **elaborato come task dentro `TASKLIST_MASTER.md`**
(formato `### A<N> · titolo — [E·C·S·R]`, con evidenza a codice), e lo stato si
aggiorna SOLO lì. Vietato creare tracker paralleli o todo-list nello scoreboard.

## Vincoli legali (CLAUDE.md §6)

- Giac GPL-3.0-or-later: solo fork/exec (`icas '<expr>'`), parsing testuale.
- VIETATO consultare/copiare sorgenti giac o derivarne algoritmi.
- Integrità binario prima di ogni scan: `bash scripts/giac_integrity.sh`.

## Protocollo

1. **Integrità oracle**:
   ```bash
   bash scripts/giac_integrity.sh
   ```
   Fail → STOP, indagare (binario cambiato = misure non confrontabili).

2. **Scan giac per area** (corpus condiviso con la golden suite Maxima):
   ```bash
   for area in diff integrate limit simplify factor solve series gcd special_fn bronstein; do
     corpus="test/golden/corpus/${area}/basic.jsonl"
     [ "$area" = bronstein ] && corpus="test/golden/corpus/bronstein/integrals.jsonl"
     bash scripts/run_golden_giac.sh "$corpus" "build-golden/giac_out/${area}"
   done
   ```
   Timeout per entry: `--per-entry-timeout N` (default 30s, env `GIAC_PER_ENTRY_TIMEOUT`).
   NON lanciare durante un gate bg (memoria: no gate concorrenti).

3. **Lato CAS**: riusa `build-golden/golden_<area>.json` esistenti; se stantii
   (codice cambiato dopo l'ultima misura) rigenera prima con
   `bash scripts/run_golden_measurement.sh --skip-maxima` (ricorda: il runner
   NON si rebuilda da solo — memoria `golden-runner-rebuild-gotcha`).

4. **Scoreboard**:
   ```bash
   python3 scripts/giac_parity_report.py
   ```
   Output: `PARITY_GIAC.md` + stampa dei candidati task (aree con Δ oltre soglia,
   default 10 pp, override `--delta-threshold-pct` / env `PARITY_DELTA_THRESHOLD_PCT`).

5. **Elaborazione task (obbligatoria se ci sono candidati)**:
   - Per ogni area candidata: verifica A CODICE la causa (module-locator /
     graphify; mai fidarsi del solo report — memoria `tasklist-ledger-stale`).
   - Scrivi la task in `TASKLIST_MASTER.md` sezione A, con: evidenza (esempi
     di input che giac risolve e noi no), causa sospetta file:riga, spec di
     riferimento in `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/` (se manca →
     prima `spec-researcher`, REGOLA 0.1).
   - Nessun fix inline durante lo scan: misura e task-writing soltanto.

## Interpretazione verdetti giac

- `ANSWERED`: forma chiusa prodotta. NON implica correttezza — per il verdetto
  matematico usare l'agente `giac-golden-diff` o la golden suite.
- `UNEVALUATED`: giac restituisce la chiamata non valutata → gap anche per giac
  (utile: se NOI rispondiamo dove giac non risponde, siamo avanti su quell'entry).
- `TIMEOUT`/`ERROR`: contare a parte; non gonfiano il gap.

## Fase 2 — cross-diff per-entry (A35, FATTA 2026-07-29)

Verdetto per-entry (non solo area-level), con equivalenza matematica VERA
(non solo ANSWERED vs pass) via `giac_parser.hpp` + confronto diretto nel
runner C++ (stessa `mathematically_equal`/`compare_solve_sets`/
`antiderivative_equivalent` già usata contro Maxima). Usa la fase 1 (scan
`.giac.out` già presenti) come dato d'ingresso — nessun nuovo scan giac.

```bash
# per ogni area d'interesse (richiede maxima_out/<area> e giac_out/<area> già presenti):
build/cas_golden_runner test/golden/corpus/<area>/basic.jsonl \
    build-golden/maxima_out/<area> --giac-dir build-golden/giac_out/<area> \
    --per-entry-json build-golden/per_entry_<area>.jsonl

# cross-diff su una o più aree:
python3 scripts/giac_per_entry_diff.py build-golden/per_entry_*.jsonl
```

Il report elenca: "giac risolve, noi no" (candidati diretti per task in
`TASKLIST_MASTER.md` — verificare a codice PRIMA), e "entrambi rispondono ma
divergono" (quasi sempre equivalenza vera non ancora provata dal motore, non
un errore — verificare con `numeric-certify` prima di aprire task). Il
verdetto giac NON tocca `AreaStats`/il ratchet: resta ancorato a Maxima.
`--giac-dir` è opt-in — senza, il runner si comporta esattamente come prima
di A35 (verificato bit-per-bit: 980/4/42 invariato).
