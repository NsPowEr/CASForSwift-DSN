#!/usr/bin/env python3
"""giac_per_entry_diff.py — A35 fase 2: cross-diff per-entry CAS vs Giac.

Consuma il dump --per-entry-json di cas_golden_runner (una riga JSON per
entry: {idx, area, ref, input, verdict, cas_output, giac_verdict,
giac_output}), prodotto passando ANCHE --giac-dir <area-dir> cosi' il
runner confronta il proprio risultato con quello di Giac (letto dai file
<idx>.giac.out gia' scritti da run_golden_giac.sh) usando la STESSA
equivalenza vera gia' in uso per Maxima (mathematically_equal /
compare_solve_sets / antiderivative_equivalent) — non solo ANSWERED vs pass.

A differenza della fase 1 (giac_parity_report.py, solo aggregati per area),
qui il verdetto e' PER ENTRY: si puo' produrre la lista esatta "giac risolve
(pass), noi no (skip/fail)" invece di un delta a livello di area.

Uso:
  # 1. genera il dump (per area; richiede maxima_out/ e giac_out/<area>/ gia' presenti):
  build/cas_golden_runner test/golden/corpus/<area>/basic.jsonl \\
      build-golden/maxima_out/<area> --giac-dir build-golden/giac_out/<area> \\
      --per-entry-json build-golden/per_entry_<area>.jsonl

  # 2. cross-diff su uno o piu' dump:
  python3 scripts/giac_per_entry_diff.py build-golden/per_entry_*.jsonl

giac GPL-3.0-or-later (CLAUDE.md §6): questo script legge solo l'OUTPUT
testuale gia' estratto dal runner C++ (fork/exec di icas, mai sorgenti giac).
"""

import argparse
import glob
import json
import sys


def load_entries(paths):
    entries = []
    for pattern in paths:
        for path in sorted(glob.glob(pattern)) or [pattern]:
            try:
                with open(path) as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        entries.append(json.loads(line))
            except FileNotFoundError:
                print(f"WARN: file non trovato, saltato: {path}", file=sys.stderr)
    return entries


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dumps", nargs="+",
                     help="uno o piu' file/pattern --per-entry-json (es. build-golden/per_entry_*.jsonl)")
    ap.add_argument("--out", default=None,
                     help="se dato, scrive il report anche su file (oltre a stdout)")
    args = ap.parse_args()

    entries = load_entries(args.dumps)
    if not entries:
        print("ERROR: nessuna entry letta — controlla i path", file=sys.stderr)
        return 1

    compared = [e for e in entries if e.get("giac_verdict") not in (None, "not_compared")]
    giac_wins = [e for e in compared
                 if e["giac_verdict"] == "pass" and e["verdict"] in ("skip", "fail")]
    both_fail_diverge = [e for e in compared
                          if e["giac_verdict"] == "fail" and e["verdict"] == "pass"]
    cas_wins = [e for e in compared
                if e["verdict"] == "pass" and e["giac_verdict"] in ("fail", "unevaluated", "timeout", "error")]

    lines = []
    lines.append(f"# Cross-diff per-entry CAS vs Giac (A35 fase 2)")
    lines.append("")
    lines.append(f"Entry totali lette: {len(entries)} · con confronto giac: {len(compared)} "
                 f"({len(entries) - len(compared)} senza --giac-dir o non ancora arrivate al confronto)")
    lines.append("")

    by_area = {}
    for e in compared:
        by_area.setdefault(e["area"], {"total": 0, "cas_pass": 0, "giac_pass": 0,
                                        "giac_wins": 0, "cas_wins": 0})
        a = by_area[e["area"]]
        a["total"] += 1
        if e["verdict"] == "pass": a["cas_pass"] += 1
        if e["giac_verdict"] == "pass": a["giac_pass"] += 1
    for e in giac_wins: by_area[e["area"]]["giac_wins"] += 1
    for e in cas_wins: by_area[e["area"]]["cas_wins"] += 1

    lines.append("| Area | Entry confrontate | CAS pass | Giac pass | Giac risolve/noi no | CAS risolve/giac no |")
    lines.append("|---|---|---|---|---|---|")
    for area in sorted(by_area):
        a = by_area[area]
        lines.append(f"| {area} | {a['total']} | {a['cas_pass']} | {a['giac_pass']} "
                     f"| {a['giac_wins']} | {a['cas_wins']} |")
    lines.append("")

    lines.append(f"## Giac risolve, noi no ({len(giac_wins)} entry)")
    lines.append("")
    lines.append("Candidati diretti per nuove task in `TASKLIST_MASTER.md` — verificare a")
    lines.append("codice PRIMA (REGOLA memoria: mai fidarsi del solo report) se il gap e' un")
    lines.append("limite del motore o del runner golden.")
    lines.append("")
    for e in giac_wins:
        lines.append(f"- **[{e['area']}#{e['idx']}]** `{e['input']}` — CAS: {e['verdict']} "
                     f"(`{e['cas_output'] or '(vuoto)'}`) · Giac: `{e['giac_output']}`")
    lines.append("")

    if both_fail_diverge:
        lines.append(f"## Entrambi rispondono ma divergono ({len(both_fail_diverge)} entry)")
        lines.append("")
        lines.append("CAS pass (vs Maxima) ma `mathematically_equal` non prova l'uguaglianza col")
        lines.append("valore di giac — quasi sempre una forma diversa non ancora riconosciuta")
        lines.append("(es. `cosh(x)` vs `(exp(x)+1/exp(x))/2`), non un errore: verificare col")
        lines.append("certificato numerico multi-punto (skill `numeric-certify`) prima di aprire task.")
        lines.append("")
        for e in both_fail_diverge:
            lines.append(f"- **[{e['area']}#{e['idx']}]** `{e['input']}` — CAS: `{e['cas_output']}` "
                         f"· Giac: `{e['giac_output']}`")
        lines.append("")

    report = "\n".join(lines) + "\n"
    print(report)
    if args.out:
        with open(args.out, "w") as f:
            f.write(report)
        print(f"Report scritto anche su: {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
