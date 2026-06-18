# Maxima Reference Oracle — Sorgente NON Modificabile (Licenza GPL-2.0-only)

> Dettaglio della Regola 6 di `CLAUDE.md`. Maxima 5.49.0 (Homebrew bottle,
> `/opt/homebrew/Cellar/maxima/5.49.0/`) è il **reference oracle** primario della
> golden test suite (F0.5).

## Divieto assoluto
- È **vietato** modificare per qualsiasi ragione i sorgenti, i binari, i file `.lisp`, `.mac`, `.fas`, `.dem` o qualunque altro artefatto contenuto in `/opt/homebrew/Cellar/maxima/`, `/opt/homebrew/share/maxima/`, o qualunque altra installazione di Maxima sul sistema.
- È **vietato** patchare, ricompilare con flag personalizzati, o re-distribuire Maxima alterato.
- È **vietato** ridistribuire output di Maxima embedded nei nostri binari senza rispetto della GPL-2.0-only.

## Motivazione
Maxima è rilasciato sotto **GPL-2.0-only**. Qualsiasi modifica al sorgente trasformerebbe il codice CAS Engine in un *derivative work* soggetto a copyleft, invalidando la nostra licenza proprietaria e invalidando l'uso di Maxima come oracolo indipendente nei test (oracolo modificato = non più indipendente, perde validità scientifica).

## Uso ammesso
- Invocazione via `maxima --very-quiet --batch-string="..."` come processo separato (fork/exec).
- Parsing dell'output testuale di Maxima per confronto AST.
- Pin esatto della versione (`5.49.0`) documentato in `PLAN_HP_PRIME_PARITY.md` e CI.
- Citazione di Maxima come reference nei doc, con link alla licenza GPL-2.0-only.

## Verifica integrità
- Lo script `scripts/verify_maxima_integrity.sh` deve calcolare lo SHA-256 dell'eseguibile Maxima e dei file `*.lisp` core, e fallire la build se non corrisponde al manifesto pinned `scripts/maxima_5.49.0_manifest.sha256`.
- Ogni golden run logga `maxima --version` + hash binario; mismatch → build fail.

## Conseguenze violazione
Rigetto immediato della PR + audit completo per identificare contaminazione GPL nel sorgente CAS Engine.
