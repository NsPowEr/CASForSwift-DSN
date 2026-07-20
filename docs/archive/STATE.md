> # ⚠️ SUPERSEDED — file storico (non aggiornare)
> Tutti i task di questo file sono stati **consolidati e ri-verificati a codice** in
> **[`TASKLIST_MASTER.md`](TASKLIST_MASTER.md)** (audit 2026-06-26, single source of truth).
> Conservato solo come storico/contesto. Le voci di stato qui dentro sono **obsolete**.

---

# STATE — Stato vivo del progetto

> **Unico file di stato corrente.** I report/handoff/summary datati vivono in
> `docs/archive/` (storico, non caricare in context salvo necessità). Aggiorna
> QUESTO file a fine sessione invece di creare un nuovo `*_SESSION_SUMMARY_*.md`.
>
> Ultimo aggiornamento: 2026-06-18

## Snapshot

- **Fasi chiuse**: F0 / F1 / F2 / F3 (suite ~1872 test).
- **Debiti F3 aperti** (ledgered): F3.1 / F3.2 / F3.3 / F3.4 / F3.5 + HC-F36.
- **Corpus golden** (ultimo misurato, 2026-06-10): aggregato ~83.9% (753/898).
  Gap residui: limit, special_fn, integrate, bronstein. Target: aggregato ≥86%,
  bronstein ≥70%, special_fn ≥82%.
- **Hardcode ledger**: vedi `python3 scripts/ledger_index.py stats` (NON aprire
  il markdown da 183 KB per una singola voce).

## Dove guardare

| Cosa | File / comando |
|---|---|
| Regole architetturali | `CLAUDE.md` (indice) → `docs/rules/` (dettaglio) |
| Task aperti / id | `python3 scripts/ledger_index.py task <CAS-Lx-yy>` |
| Hardcode di passaggio | `python3 scripts/ledger_index.py hc <id>` / `open` |
| Ricerca trasversale | `python3 scripts/ledger_index.py search "<testo>"` |
| Piani attivi | `PLAN_HP_PRIME_PARITY.md`, `PLAN_F3_F8_GAP_CLOSURE.md`, `PLAN_TASKS_REMAINING.md` |
| Definition of Done | `DEFINITION_OF_DONE.md` |
| Navigazione codice | `graphify query "<domanda>"` (NON grep raw a 4x scala) |
| Storico sessioni | `docs/archive/` |

## Prossimi passi

Vedi `PLAN_NEXT_SESSIONS.md` e i `TODO.md` / `TODO_PH8.md` (pilastri P6–P8
ancora aperti: solve/RootOf, limit pipeline, integration core Risch).
