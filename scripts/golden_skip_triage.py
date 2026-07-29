#!/usr/bin/env python3
"""golden_skip_triage.py — classifica gli SKIP per-entry della golden suite.

Il JSON per area (`build-golden/golden_<area>.json`) porta solo aggregati
pass/fail/skip: uno skip puo' nascere da cause del tutto diverse (limite del
runner, gap del motore, oracolo Maxima muto, parser Maxima, equivalenza
inconcludente) e trattarle come un'unica voce porta a diagnosi sbagliate
(caso A37/gcd: 11 "gap motore" erano in realta' un limite del runner).

Questo script legge i log per-entry prodotti da run_golden_measurement.sh
(`build-golden/log_<area>.txt`, riga `  SKIP [ idx] input => motivo`) e li
classifica in bucket disgiunti:

  runner    limite del RUNNER golden (dispatch assente, arg-count, area skippata
            by design) — si chiude modificando test/golden/, non il motore
  engine    Unimplemented emesso dal MOTORE — gap reale, candidato task A<N>
  oracle    Maxima non ha prodotto output utilizzabile — ne' nostro ne' loro
  parser    output Maxima non parsato dal NOSTRO maxima_parser.hpp — gap nostro
            ma nel tooling di confronto, non nel motore
  equality  CAS e Maxima hanno entrambi risposto ma mathematically_equal e'
            inconcludente — gap simplify/equivalenza (la risposta CAS puo'
            essere giusta)

Uso:
  python3 scripts/golden_skip_triage.py                 # tabella tutte le aree
  python3 scripts/golden_skip_triage.py --area limit    # una sola area
  python3 scripts/golden_skip_triage.py --list          # + ogni entry skippata
  python3 scripts/golden_skip_triage.py --bucket engine # solo un bucket
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG_GLOB = os.path.join(ROOT, "build-golden", "log_*.txt")

SKIP_RE = re.compile(r"^\s*SKIP\s*\[\s*(\d+)\]\s(.*)$")

# (bucket, regex sul motivo). Primo match vince: l'ordine e' significativo,
# i pattern del runner sono i piu' specifici e vanno testati per primi.
RULES: list[tuple[str, re.Pattern[str]]] = [
    ("runner", re.compile(r"matrix fn skipped")),
    ("runner", re.compile(r"too few args|bad order|not a solve\(")),
    ("runner", re.compile(r"set equality not yet supported")),
    ("runner", re.compile(r"no dispatch kind")),
    ("oracle", re.compile(r"no Maxima output")),
    # Maxima itself errored out on the entry (needs an assume, unsupported
    # form, ...): its text lands in the "parse fail" branch but the gap is on
    # the oracle side, not in our parser.
    ("oracle", re.compile(r"an error\. To debug")),
    ("oracle", re.compile(r"Maxima (scalar|eig) (empty|inner empty)")),
    ("oracle", re.compile(r"\(Maxima: '")),
    ("parser", re.compile(r"Maxima (parse|matrix parse|eig elem parse)")),
    ("equality", re.compile(r"inconclusive")),
]


def classify(reason: str) -> str:
    for bucket, rx in RULES:
        if rx.search(reason):
            return bucket
    # Tutto cio' che resta e' un Unimplemented risalito dal motore attraverso
    # evaluate_cas: e' il bucket che genera task.
    return "engine"


def normalize(reason: str) -> str:
    """Chiave di raggruppamento: motivo senza payload variabile."""
    r = re.sub(r"\d+", "N", reason)
    r = re.sub(r"\s+", " ", r).strip()
    return r[:120]


def parse_log(path: str) -> list[tuple[int, str, str, str]]:
    """-> [(idx, input, reason, bucket)]"""
    out = []
    for line in open(path, encoding="utf-8", errors="replace"):
        m = SKIP_RE.match(line.rstrip("\n"))
        if not m:
            continue
        idx = int(m.group(1))
        rest = m.group(2)
        if " => " in rest:
            inp, reason = rest.split(" => ", 1)
        else:
            # forma `input (motivo)`
            k = rest.rfind(" (")
            inp, reason = (rest[:k], rest[k + 1:]) if k > 0 else (rest, "")
        out.append((idx, inp.strip(), reason.strip(), classify(reason)))
    return out


BUCKETS = ["runner", "engine", "parser", "oracle", "equality"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--area", help="una sola area (default: tutte)")
    ap.add_argument("--list", action="store_true", help="elenca ogni entry skippata")
    ap.add_argument("--bucket", choices=BUCKETS, help="filtra su un bucket")
    args = ap.parse_args()

    logs = sorted(glob.glob(LOG_GLOB))
    if args.area:
        logs = [p for p in logs if os.path.basename(p) == f"log_{args.area}.txt"]
    if not logs:
        print("nessun log in build-golden/log_*.txt — esegui prima "
              "bash scripts/run_golden_measurement.sh [--area <a>] [--skip-maxima]",
              file=sys.stderr)
        return 2

    per_area: dict[str, list] = {}
    for p in logs:
        area = os.path.basename(p)[len("log_"):-len(".txt")]
        per_area[area] = parse_log(p)

    width = max(len(a) for a in per_area) + 1
    print(f"{'area':<{width}} {'skip':>5} " + " ".join(f"{b:>9}" for b in BUCKETS))
    print("-" * (width + 6 + 10 * len(BUCKETS)))
    totals: Counter[str] = Counter()
    for area, rows in sorted(per_area.items()):
        c = Counter(b for _, _, _, b in rows)
        totals.update(c)
        totals["skip"] += len(rows)
        print(f"{area:<{width}} {len(rows):>5} " +
              " ".join(f"{c.get(b, 0):>9}" for b in BUCKETS))
    print("-" * (width + 6 + 10 * len(BUCKETS)))
    print(f"{'TOTALE':<{width}} {totals['skip']:>5} " +
          " ".join(f"{totals.get(b, 0):>9}" for b in BUCKETS))

    # Motivi ricorrenti per bucket: e' la lista di lavoro effettiva.
    print("\nMotivi ricorrenti (bucket · occorrenze · motivo normalizzato):")
    grouped: dict[tuple[str, str], int] = defaultdict(int)
    for rows in per_area.values():
        for _, _, reason, bucket in rows:
            if args.bucket and bucket != args.bucket:
                continue
            grouped[(bucket, normalize(reason))] += 1
    for (bucket, reason), n in sorted(grouped.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {bucket:<9} {n:>4}  {reason}")

    if args.list:
        print("\nEntry skippate:")
        for area, rows in sorted(per_area.items()):
            for idx, inp, reason, bucket in rows:
                if args.bucket and bucket != args.bucket:
                    continue
                print(f"  [{area}/{idx:>3}] {bucket:<9} {inp}\n"
                      f"        -> {reason}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
