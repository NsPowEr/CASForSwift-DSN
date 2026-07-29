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
- Pin esatto della versione (`5.49.0`) documentato in `docs/archive/PLAN_HP_PRIME_PARITY.md` (storico) e nel manifest `scripts/maxima_5.49.0_manifest.sha256`.
- Citazione di Maxima come reference nei doc, con link alla licenza GPL-2.0-only.

## Verifica integrità
- Lo script `scripts/verify_maxima_integrity.sh` deve calcolare lo SHA-256 dell'eseguibile Maxima e dei file `*.lisp` core, e fallire la build se non corrisponde al manifesto pinned `scripts/maxima_5.49.0_manifest.sha256`.
- Ogni golden run logga `maxima --version` + hash binario; mismatch → build fail.

## Conseguenze violazione
Rigetto immediato della PR + audit completo per identificare contaminazione GPL nel sorgente CAS Engine.

---

# Giac — Secondo Oracle e Target di Parità (GPL-3.0-or-later)

*Aggiunto 2026-07-19 (F7.5.G1 attivato). Stesse regole di Maxima, licenza diversa.*

## Installazione pinned
- Giac **2.0.0** arm64 nativo: `~/xcas-oracle/` (payload estratto dal `xcas_mac.dmg.gz` ufficiale di B. Parisse, build 2026-03-02), CLI `icas` symlinkato in `/opt/homebrew/bin/icas`.
- Shim ambiente (NON modifica di giac): `~/xcas-oracle/libs_arm64/libreadline.8.2.dylib` è un symlink alla readline Homebrew — il pacchetto upstream ometteva la dylib arm64.
- Integrità: `bash scripts/giac_integrity.sh` (manifest `scripts/giac_2.0.0_manifest.sha256`); `--pin` solo dopo upgrade autorizzato dall'utente.

## Divieto assoluto
- Vietato modificare/patchare/ricompilare binari o sorgenti giac, o consultarne il codice (C++/tabelle interne) per derivare implementazioni: GPL-3.0 → derivative work copyleft + oracolo non più indipendente.
- Doppio ruolo di giac (oracolo di verifica **e** target di parità) NON cambia la regola: la parità si raggiunge implementando gli algoritmi dalle spec di `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/`, mai leggendo come lo fa giac.

## Uso ammesso
- Fork/exec di `icas` con comandi su stdin (una espressione per riga, terminata da `;`) + parsing testuale dell'output.
- Formato output: `N>> <input>` eco, risultato sulla riga successiva, `// dclock1 <t>` = timing. Banner iniziale e warning locale/doc = rumore ignorabile.
- Agent dedicato: `.claude/agents/giac-golden-diff.md`.
