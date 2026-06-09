# Runner Robustness — per-entry timeout + output cap

> Task F7.5.A3 — stabilizza `cas_golden_runner` su corpus pesanti
> (bronstein hang, integrate truncate). Misurazione end-to-end di
> tutte le aree senza interruzioni.

## Problema

1. **Bronstein entry hang**: corpus `bronstein/integrals.jsonl` entry
   0 manda il runner in loop / ricorsione profonda. Senza timeout
   per-entry, nessuna entry successiva viene eseguita.
2. **Integrate output explosion**: `integrate(sin(x)^N, x)` per N ≥ 5
   produce espressione di centinaia di KB stampata a stdout. Buffer
   shell tronca output globale (vedi F7_GOLDEN_CORPUS_REPORT note
   integrate 116/140).
3. **Nessun watchdog globale**: una entry mal formata può consumare
   tutta la sessione.

## Specifica formale

### A. Timeout per-entry (`SIGALRM` + `setitimer`)

POSIX timer `ITIMER_REAL`:
- arm prima di `evaluate_cas` con `value = ctx_options.per_entry_timeout_seconds`
  (default 30 s, configurabile via CLI `--per-entry-timeout N`);
- handler `SIGALRM` setta flag `volatile sig_atomic_t timeout_flag = 1`
  + chiama `cas::symbolic::CASContext::interrupt()` (atomic) per
  arrestare loop simbolici;
- disarm dopo l'entry;
- se `timeout_flag` raised → entry contata come SKIP con messaggio
  `"TIMEOUT after N s"`, NON FAIL.

**Vincolo (CLAUDE.md A3.3)**: il cancellation token atomico è già
integrato nei poll-point principali del simplifier e substitute. Il
runner aggancia solo l'orchestrazione di alto livello.

**Vietato**: thread separati con `pthread_kill` (segfault risk),
`alarm()` plain (intero second-precision, non integrabile con i
poll-points), retry su timeout (ledger antiloop).

### B. Output cap per entry

`format_expr` può ritornare megabyte. Cap configurabile:

```cpp
constexpr std::size_t kMaxFormatBytes = 4096;  // 4 KB per entry default
std::string truncated_format(ExprPtr e) {
    std::string s = format_expr(e);
    if (s.size() > kMaxFormatBytes) {
        s.resize(kMaxFormatBytes - 32);
        s += " ... <truncated " + std::to_string(format_expr(e).size() -
              kMaxFormatBytes) + " bytes>";
    }
    return s;
}
```

Applicato solo a output diagnostico (FAIL CAS/Maxima dump). La
comparazione `mathematically_equal` resta sull'AST completo — la cap
non altera la correttezza, solo la verbosità.

### C. Configurazione CLI

```
cas_golden_runner <corpus.jsonl> <maxima_dir>
                  [--json out.json]
                  [--per-entry-timeout N]   (default 30)
                  [--format-cap-bytes N]    (default 4096)
                  [--verbose]
```

### D. Cleanup

Garantire `setitimer(ITIMER_REAL, &zero, NULL)` su ogni path di uscita
(normale, signal, exception → non applicabile, no throw).

## Acceptance criteria

- Corpus `bronstein/integrals.jsonl` (90 entry) eseguito end-to-end,
  almeno 1 entry passa (smoke test, no algoritmica). Aggregato finale
  riportato senza hang.
- Corpus `integrate/basic.jsonl` (140 entry): runner completa
  140/140 (oggi 116/140 troncato dal buffer).
- Nessun output > 4 KB per entry.
- Timeout default 30 s ragionevole (tutti i test correntemente passing
  completano entro questo budget; verificato benchmarking).

## File da modificare/creare

- `test/golden/main.cpp` — integrazione `setitimer` + cap.
- `test/golden/runner_timeout.hpp` (nuovo) — incapsulamento POSIX
  timer (RAII), forward-compat con altri runner.
- `test/golden/runner_format.hpp` (nuovo) — `truncated_format`.

## Vincoli (CLAUDE.md)

- 500 LOC max. main.cpp oggi 432 LOC; +60 stimati post-integrazione.
- No throw/catch. Signal handler usa `volatile sig_atomic_t` (POSIX).
- No hardcode-of-passage. Default 30 s + 4 KB documentati come
  parametri CLI.
- Cleanup garantito.
- Spec read first (REGOLA 0.1).

## Note implementative

Il signal handler **non** può chiamare funzioni non-async-signal-safe.
Si limita a:
1. Settare `timeout_flag = 1` (volatile sig_atomic_t).
2. Chiamare `ctx.interrupt()` SOLO se l'interrupt è dichiarato
   async-signal-safe (verifica: `std::atomic<bool>::store` con
   `memory_order_relaxed` su lock-free atomic è OK per POSIX).

Se il check `ctx.interrupted()` non viene letto dai poll-point
correntemente, il timer scade ma la entry continua a girare. In quel
caso il runner ha bisogno di un secondo livello: timeout hard tramite
processo figlio. Fuori scope F7.5.A3 — tracciato come ledger entry
F7.5-A3-HARD-TIMEOUT se il signal soft non basta.
