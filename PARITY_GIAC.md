# PARITY GIAC — Scoreboard CAS vs Giac 2.0.0

> Generato: 2026-07-19T15:54:44 · `python3 scripts/giac_parity_report.py`
> **Artefatto di MISURA, non tracker**: ogni gap da chiudere va elaborato
> come task in `TASKLIST_MASTER.md` (single source of truth).
> giac% = risposte in forma chiusa di giac sul corpus (copertura oracle).
> CAS% = pass del CAS vs Maxima sullo stesso corpus (correttezza+copertura).
> Δ = giac% − CAS%: positivo grande ⇒ area dove giac ci batte.

| Area | Entries | giac answered | giac uneval | giac timeout/err | giac% | CAS pass/fail/skip | CAS% | Δ pp |
|---|---|---|---|---|---|---|---|---|
| bronstein | 90 | 86 | 4 | 0 | 96% | 53/9/28 | 59% | +37 |
| diff | 80 | 80 | 0 | 0 | 100% | 80/0/0 | 100% | +0 |
| factor | 99 | 99 | 0 | 0 | 100% | 99/0/0 | 100% | +0 |
| gcd | 81 | 81 | 0 | 0 | 100% | 70/0/11 | 86% | +14 |
| integrate | 140 | 140 | 0 | 0 | 100% | 121/1/18 | 86% | +14 |
| limit | 99 | 98 | 0 | 1 | 99% | 88/0/11 | 89% | +10 |
| matrix | 79 | 79 | 0 | 0 | 100% | 79/0/0 | 100% | +0 |
| series | 81 | 81 | 0 | 0 | 100% | 72/0/9 | 89% | +11 |
| simplify | 116 | 104 | 12 | 0 | 90% | 116/0/0 | 100% | -10 |
| solve | 81 | 81 | 0 | 0 | 100% | 81/0/0 | 100% | +0 |
| special_fn | 80 | 34 | 46 | 0 | 42% | 75/0/5 | 94% | -51 |

⚠ **Dati CAS stantii** (mtime `golden_<area>.json` più vecchio di 7g rispetto al più recente — il Δ misura un motore vecchio, rigenerare con `run_golden_measurement.sh --area <a> --skip-maxima` dopo rebuild del runner):
- `bronstein`: 28 giorni indietro

## Aree oltre soglia (Δ > 10 pp) — candidati task

Da elaborare in `TASKLIST_MASTER.md` (formato `### A<N> · titolo — [E·C·S·R]`),
previa verifica a codice della causa (REGOLA memoria: mai fidarsi del solo report).
Caveat fase 1 (fino ad A35): niente attribuzione per-entry — il Δ è area-level;
parte degli skip CAS sono limiti del RUNNER golden (matrix literal, gcd multivariato),
non del motore: distinguere a codice prima di aprire task sul motore.

- **bronstein** (Δ +37 pp) — lato CAS: pass 53 / fail 9 / skip 28 → indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.
- **gcd** (Δ +14 pp) — lato CAS: pass 70 / fail 0 / skip 11 → indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.
- **integrate** (Δ +14 pp) — lato CAS: pass 121 / fail 1 / skip 18 → indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.
- **limit** (Δ +10 pp) — lato CAS: pass 88 / fail 0 / skip 11 → indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.
- **series** (Δ +11 pp) — lato CAS: pass 72 / fail 0 / skip 9 → indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.

