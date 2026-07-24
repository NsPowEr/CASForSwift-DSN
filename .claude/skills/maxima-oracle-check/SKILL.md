---
name: maxima-oracle-check
description: Verifica integrità Maxima 5.49.0 (SHA-256 manifest) ed esegue golden diff su un'espressione contro l'oracle indipendente. Maxima usato SOLO come verificatore esterno (fork/exec + parsing output testuale). Vietato copiare sorgenti Lisp/algoritmi Maxima nel codice CAS (GPL-2.0-only, CLAUDE.md §6).
---

# maxima-oracle-check

Confronto golden contro Maxima come **oracle indipendente**. Mai sorgente di implementazione.

## Vincoli legali e scientifici (CLAUDE.md §6)

- Maxima è GPL-2.0-only. Sorgenti `/opt/homebrew/Cellar/maxima/5.49.0/**` **immutabili**.
- **Vietato**: copiare codice Lisp, algoritmi, tabelle, costanti dai `*.lisp` Maxima nel codice CAS (deriverebbe copyleft + invalida indipendenza oracle).
- **Ammesso**: invocazione `maxima --very-quiet --batch-string="..."` come processo separato; parsing output testuale per confronto AST.
- Implementazioni CAS derivano da specifiche formali in `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/`, non da Maxima.

## Protocollo

1. **Integrità binario + manifest** (obbligatorio prima di ogni golden run):
   ```bash
   bash scripts/verify_maxima_integrity.sh
   ```
   Confronta SHA-256 di `maxima` + core `*.lisp` contro `scripts/maxima_5.49.0_manifest.sha256`. Fail → STOP, indagare contaminazione/aggiornamento non autorizzato.

2. **Golden diff su espressione singola** (interfaccia REALE: corpus jsonl — i flag `--expr/--op` NON esistono):
   ```bash
   TMP=$(mktemp -d); printf '%s\n' '{"input": "<EXPR>", "area": "<simplify|factor|integrate|diff|limit|solve|series|special_fn|matrix|bronstein|gcd>", "ref": "adhoc_1"}' > "$TMP/mini.jsonl"
   bash scripts/run_golden_maxima.sh "$TMP/mini.jsonl" "$TMP/maxrefs"
   build/cas_golden_runner "$TMP/mini.jsonl" "$TMP/maxrefs" --json "$TMP/report.json"
   ```
   Runner stale? `ninja -C build cas_golden_runner` prima di misurare. Mismatch → bug nel CAS, non nel test.

3. **Cross-check SymPy** (BSD, secondo oracle di conferma):
   ```bash
   python3 scripts/run_golden_sympy.py "$TMP/mini.jsonl" "$TMP/sympyrefs"
   ```
   Stessa policy: spunto, non copia sorgente.

4. **Log obbligatorio**:
   - `maxima --version` + hash binario;
   - mismatch hash → build fail immediato (CI).

## Uso tipico

- Validare un nuovo modulo (es. Risch sub-case, Zeilberger telescopia) prima del merge.
- Verificare che un fix di bug non introduca regressione matematica su un input noto.
- Sanity check pre-rilascio.

## Output atteso

- PASS → "AST equivalente a Maxima per <EXPR>" + log hash.
- FAIL → diff AST atteso/prodotto; mai modificare il test per aggirare il fail (CLAUDE.md §REGOLA 0.2). Fix nel motore.
