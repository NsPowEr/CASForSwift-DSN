# PARITY GIAC — Scoreboard CAS vs Giac 2.0.0

> Generato: 2026-07-24T03:21:23 · `python3 scripts/giac_parity_report.py`
> **Artefatto di MISURA, non tracker**: ogni gap da chiudere va elaborato
> come task in `TASKLIST_MASTER.md` (single source of truth).
> giac% = risposte in forma chiusa di giac sul corpus (copertura oracle).
> CAS% = pass del CAS vs Maxima sullo stesso corpus (correttezza+copertura).
> Δ = giac% − CAS%: positivo grande ⇒ area dove giac ci batte.

| Area | Entries | giac answered | giac uneval | giac timeout/err | giac% | CAS pass/fail/skip | CAS% | Δ pp |
|---|---|---|---|---|---|---|---|---|
| bronstein | 90 | 86 | 4 | 0 | 96% | 60/5/25 | 67% | +29 |
| diff | 80 | 80 | 0 | 0 | 100% | 80/0/0 | 100% | +0 |
| factor | 99 | 99 | 0 | 0 | 100% | 99/0/0 | 100% | +0 |
| gcd | 81 | 81 | 0 | 0 | 100% | 78/3/0 | 96% | +4 |
| integrate | 140 | 140 | 0 | 0 | 100% | 122/2/16 | 87% | +13 |
| limit | 99 | 98 | 0 | 1 | 99% | 88/0/11 | 89% | +10 |
| matrix | 79 | 79 | 0 | 0 | 100% | 79/0/0 | 100% | +0 |
| series | 81 | 81 | 0 | 0 | 100% | 79/0/2 | 98% | +2 |
| simplify | 116 | 104 | 12 | 0 | 90% | 116/0/0 | 100% | -10 |
| solve | 81 | 81 | 0 | 0 | 100% | 81/0/0 | 100% | +0 |
| special_fn | 80 | 34 | 46 | 0 | 42% | 79/0/1 | 99% | -56 |

## Aree oltre soglia (Δ > 10 pp) — candidati task

Da elaborare in `TASKLIST_MASTER.md` (formato `### A<N> · titolo — [E·C·S·R]`),
previa verifica a codice della causa (REGOLA memoria: mai fidarsi del solo report).
Caveat fase 1 (fino ad A35): niente attribuzione per-entry — il Δ è area-level;
parte degli skip CAS sono limiti del RUNNER golden (matrix literal, gcd multivariato),
non del motore: distinguere a codice prima di aprire task sul motore.

- **bronstein** (Δ +29 pp) — lato CAS: pass 60 / fail 5 / skip 25 → indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.
- **integrate** (Δ +13 pp) — lato CAS: pass 122 / fail 2 / skip 16 → indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.
- **limit** (Δ +10 pp) — lato CAS: pass 88 / fail 0 / skip 11 → indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.

