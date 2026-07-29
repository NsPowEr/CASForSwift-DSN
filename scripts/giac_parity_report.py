#!/usr/bin/env python3
"""giac_parity_report.py — merge giac coverage scans with CAS golden results
into a per-area parity scoreboard (PARITY_GIAC.md).

Phase 1 (this script): area-level parity.
  - giac side: build-golden/giac_out/<area>/manifest.json + per-entry
    <i>.giac.out verdicts (from scripts/run_golden_giac.sh).
  - CAS side: build-golden/golden_<area>.json aggregates (pass/fail/skip vs
    Maxima) produced by cas_golden_runner.

Phase 2 (task in TASKLIST_MASTER.md): entry-level cross-diff once
cas_golden_runner grows a --per-entry-json dump. This script already reads
per-entry giac verdicts, so it will consume that dump without rework.

The scoreboard is a MEASUREMENT ARTIFACT, not a task tracker: any gap worth
fixing must be elaborated as a task inside TASKLIST_MASTER.md (single source
of truth). This script prints ready-to-adapt candidate task lines for areas
whose parity delta exceeds the reporting threshold, but never writes to
TASKLIST_MASTER.md itself.

Usage:
  python3 scripts/giac_parity_report.py \
      [--giac-dir build-golden/giac_out] \
      [--golden-dir build-golden] \
      [--corpus-dir test/golden/corpus] \
      [--out PARITY_GIAC.md] \
      [--delta-threshold-pct 10]
"""

import argparse
import datetime
import json
import os
import sys


def read_giac_area(giac_area_dir):
    """Return (manifest dict or None, list of (idx, verdict))."""
    manifest_path = os.path.join(giac_area_dir, "manifest.json")
    manifest = None
    if os.path.isfile(manifest_path):
        with open(manifest_path) as f:
            manifest = json.load(f)
    verdicts = []
    if os.path.isdir(giac_area_dir):
        for name in sorted(os.listdir(giac_area_dir)):
            if not name.endswith(".giac.out"):
                continue
            idx = int(name.split(".")[0])
            with open(os.path.join(giac_area_dir, name)) as f:
                verdict = f.readline().strip()
            verdicts.append((idx, verdict))
    verdicts.sort()
    return manifest, verdicts


def read_cas_area(golden_dir, area):
    """Return aggregate dict {pass, fail, skip} or None."""
    path = os.path.join(golden_dir, f"golden_{area}.json")
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        data = json.load(f)
    agg = {"pass": 0, "fail": 0, "skip": 0}
    for v in data.values():
        if isinstance(v, dict):
            for k in agg:
                agg[k] += v.get(k, 0)
    return agg


def read_corpus_inputs(corpus_dir, area):
    """Return list of input strings for an area (index-aligned with giac outs)."""
    if area == "bronstein":
        path = os.path.join(corpus_dir, area, "integrals.jsonl")
    else:
        path = os.path.join(corpus_dir, area, "basic.jsonl")
    inputs = []
    if os.path.isfile(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    inputs.append(json.loads(line).get("input", ""))
                except json.JSONDecodeError:
                    inputs.append("")
    return inputs


def pct(n, d):
    return 100.0 * n / d if d else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--giac-dir", default="build-golden/giac_out")
    ap.add_argument("--golden-dir", default="build-golden")
    ap.add_argument("--corpus-dir", default="test/golden/corpus")
    ap.add_argument("--out", default="PARITY_GIAC.md")
    ap.add_argument("--delta-threshold-pct", type=float,
                    default=float(os.environ.get("PARITY_DELTA_THRESHOLD_PCT", 10)))
    args = ap.parse_args()

    if not os.path.isdir(args.giac_dir):
        print(f"ERROR: giac scan dir not found: {args.giac_dir}\n"
              f"Run scripts/run_golden_giac.sh per area first "
              f"(see .claude/skills/giac-parity-scan).", file=sys.stderr)
        return 1

    areas = sorted(d for d in os.listdir(args.giac_dir)
                   if os.path.isdir(os.path.join(args.giac_dir, d)))
    if not areas:
        print(f"ERROR: no area subdirectories in {args.giac_dir}", file=sys.stderr)
        return 1

    rows = []
    candidate_tasks = []
    cas_mtimes = {}
    for area in areas:
        gp = os.path.join(args.golden_dir, f"golden_{area}.json")
        if os.path.isfile(gp):
            cas_mtimes[area] = os.path.getmtime(gp)
    for area in areas:
        manifest, verdicts = read_giac_area(os.path.join(args.giac_dir, area))
        if manifest is None:
            print(f"WARN: no manifest for area '{area}' — skipped", file=sys.stderr)
            continue
        cas = read_cas_area(args.golden_dir, area)
        total = manifest["total"]
        giac_ok = manifest["answered"]
        giac_rate = pct(giac_ok, total)
        if cas:
            cas_total = cas["pass"] + cas["fail"] + cas["skip"]
            cas_rate = pct(cas["pass"], cas_total)
            delta = giac_rate - cas_rate
        else:
            cas_rate = None
            delta = None
        rows.append((area, total, manifest, cas, giac_rate, cas_rate, delta))

        if delta is not None and delta > args.delta_threshold_pct:
            # Phase 1 has no per-entry CAS verdicts (task A35): entry-level
            # attribution ("giac solves it AND we don't") is impossible here.
            # Listing giac-ANSWERED examples would be misleading (they include
            # entries the CAS also passes), so candidates carry only the
            # area-level delta plus the CAS-side aggregate breakdown.
            candidate_tasks.append((area, delta, cas))

    now = datetime.datetime.now().isoformat(timespec="seconds")
    lines = []
    lines.append("# PARITY GIAC — Scoreboard CAS vs Giac 2.0.0")
    lines.append("")
    lines.append(f"> Generato: {now} · `python3 scripts/giac_parity_report.py`")
    lines.append("> **Artefatto di MISURA, non tracker**: ogni gap da chiudere va elaborato")
    lines.append("> come task in `TASKLIST_MASTER.md` (single source of truth).")
    lines.append("> giac% = risposte in forma chiusa di giac sul corpus (copertura oracle).")
    lines.append("> CAS% = pass del CAS vs Maxima sullo stesso corpus (correttezza+copertura).")
    lines.append("> Δ = giac% − CAS%: positivo grande ⇒ area dove giac ci batte.")
    lines.append("")
    lines.append("| Area | Entries | giac answered | giac uneval | giac timeout/err | giac% | CAS pass/fail/skip | CAS% | Δ pp |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for (area, total, m, cas, giac_rate, cas_rate, delta) in rows:
        te = m["timeout"] + m["error"]
        cas_cell = f"{cas['pass']}/{cas['fail']}/{cas['skip']}" if cas else "n/a"
        cas_rate_cell = f"{cas_rate:.0f}%" if cas_rate is not None else "n/a"
        delta_cell = f"{delta:+.0f}" if delta is not None else "n/a"
        lines.append(f"| {area} | {total} | {m['answered']} | {m['unevaluated']} "
                     f"| {te} | {giac_rate:.0f}% | {cas_cell} | {cas_rate_cell} | {delta_cell} |")
    lines.append("")

    # Staleness warning: CAS golden data much older than the freshest area is
    # measuring an old engine — deltas there are not comparable.
    if cas_mtimes:
        newest = max(cas_mtimes.values())
        stale_days = float(os.environ.get("PARITY_STALE_DAYS", 7))
        stale = [(a, (newest - t) / 86400.0) for a, t in sorted(cas_mtimes.items())
                 if (newest - t) > stale_days * 86400.0]
        if stale:
            lines.append(f"⚠ **Dati CAS stantii** (mtime `golden_<area>.json` più vecchio di "
                         f"{stale_days:.0f}g rispetto al più recente — il Δ misura un motore vecchio, "
                         f"rigenerare con `run_golden_measurement.sh --area <a> --skip-maxima` "
                         f"dopo rebuild del runner):")
            for a, days in stale:
                lines.append(f"- `{a}`: {days:.0f} giorni indietro")
            lines.append("")

    if candidate_tasks:
        lines.append(f"## Aree oltre soglia (Δ > {args.delta_threshold_pct:.0f} pp) — candidati task")
        lines.append("")
        lines.append("Da elaborare in `TASKLIST_MASTER.md` (formato `### A<N> · titolo — [E·C·S·R]`),")
        lines.append("previa verifica a codice della causa (REGOLA memoria: mai fidarsi del solo report).")
        lines.append("Caveat fase 1 (fino ad A35): niente attribuzione per-entry — il Δ è area-level;")
        lines.append("parte degli skip CAS sono limiti del RUNNER golden (matrix literal, gcd multivariato),")
        lines.append("non del motore: distinguere a codice prima di aprire task sul motore.")
        lines.append("")
        for area, delta, cas in candidate_tasks:
            cas_cell = (f"pass {cas['pass']} / fail {cas['fail']} / skip {cas['skip']}"
                        if cas else "n/a")
            lines.append(f"- **{area}** (Δ {delta:+.0f} pp) — lato CAS: {cas_cell} "
                         f"→ indagare prima i FAIL (possibili silent-wrong), poi gli SKIP.")
        lines.append("")

    with open(args.out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Scoreboard scritto: {args.out} ({len(rows)} aree)")
    if candidate_tasks:
        print(f"Candidati task: {', '.join(a for a, _, _ in candidate_tasks)} "
              f"→ elaborare in TASKLIST_MASTER.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
