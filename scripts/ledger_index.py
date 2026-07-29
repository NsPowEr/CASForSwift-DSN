#!/usr/bin/env python3
"""ledger_index.py — queryable index over the large markdown trackers.

Why: HARDCODE_LEDGER.md (~200 KB) and the task trackers (TASKLIST_MASTER.md +
archived CAS_TASKS.md / PLAN_*.md …) are linear markdown. Loading a whole
100-200 KB file to read ONE entry wastes the context budget; as the project
grows these files dominate it. This tool keeps the markdown as the editable
source of truth and builds a disposable SQLite+FTS index in .cache/ (git-ignored,
rebuilt on demand), so agents query targeted slices.

Permanent design choices (so it keeps working as the project grows):
  * Source set is DISCOVERED via globs (SOURCE_SPECS) → new trackers are indexed
    with zero code change. SUPERSEDED files self-mark with a top banner and are
    auto-tagged `archived`, so current state is distinguishable from history.
  * Ledger status is derived by an EXPLAINABLE 3-tier rule (body `**Stato**:`
    line → header verdict → section) and every row carries status_source +
    status_raw, so the classification is auditable. Ambiguity surfaces as
    `unknown` — never a silent wrong CLOSED. (Fixes the historical
    "ledger_index status flags unreliable" footgun.)
  * Search never crashes and never silently drops punctuation-bearing queries:
    ranked FTS5 + a substring (LIKE) pass are merged.

Commands: build · hc <id> · task <id> · search <text> · open · stats · doctor ·
sources.  Global flags: --json (machine output), --all (include archived task
files / unknown-status entries). Non-destructive: never writes to the markdown.
"""
from __future__ import annotations

import argparse
import contextlib
import glob
import io
import json
import os
import re
import sqlite3
import sys

SCHEMA_VERSION = 3
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_DIR = os.path.join(ROOT, ".cache")
DB_PATH = os.path.join(DB_DIR, "ledger_index.db")

# role 'ledger' = `### <ID> …` entries with status; 'tasks' = trackers (headings
# + id refs). Glob-driven so future PLAN_*/TASKLIST_*/_LEDGER files auto-index.
SOURCE_SPECS: list[tuple[str, str]] = [
    ("ledger", "HARDCODE_LEDGER.md"), ("ledger", "*_LEDGER.md"),
    ("tasks", "TASKLIST_MASTER.md"), ("tasks", "TASKLIST_*.md"),
    ("tasks", "CAS_TASKS.md"), ("tasks", "PLAN_*.md"), ("tasks", "TODO*.md"),
    ("tasks", "STATE.md"), ("tasks", "HANDOFF_*.md"), ("tasks", "*_SESSION_*.md"),
    # Tracker superseded archiviati (2026-07-19): restano indicizzati per la
    # tracciabilità Refs; il banner SUPERSEDED li classifica già come storici.
    ("tasks", "docs/archive/CAS_TASKS.md"), ("tasks", "docs/archive/PLAN_*.md"),
    ("tasks", "docs/archive/TODO*.md"), ("tasks", "docs/archive/STATE.md"),
    ("tasks", "docs/archive/HANDOFF_*.md"), ("tasks", "docs/archive/*_SESSION_*.md"),
]
# A tracker marks ITSELF historical with a top banner, e.g.
#   `> # ⚠️ SUPERSEDED — file storico (non aggiornare)`
# Match only the self-banner (first non-empty line, leading markdown/emoji
# stripped) so a CURRENT file that merely *names* superseded files (e.g.
# TASKLIST_MASTER.md, the single source of truth) is not mis-tagged.
ARCHIVED_BANNER = re.compile(r"^(superseded|deprecat|file storico)", re.I)

STRONG_ID_RES = [re.compile(p) for p in (
    r"\bCAS-L\d+-\d+[A-Za-z]?\b", r"\bT-\d{2,4}\b", r"\bHPP-[0-9A-Za-z.]+",
    r"\bHC-[0-9A-Za-z._-]+", r"\bBUG-[0-9A-Za-z._-]+",
)]
# Leading token of a heading, kept only if id-shaped — so TASKLIST_MASTER local
# ids (`### A2 · …`) are queryable without matching prose headings.
HEADING_ID_RE = re.compile(
    r"^(?:A\d{1,3}|T-\d{2,4}|CAS-L\d+-\d+[A-Za-z]?|HPP-[\w.]+|HC-[\w.-]+|"
    r"BUG-[\w.-]+|F\d[\w.-]*)$")

# Status vocabulary (Italian + English + emoji), evaluated in order.
STATUS_RULES: list[tuple[str, re.Pattern]] = [
    ("wontfix", re.compile(r"ACCETTAT[OA]|WON'?T[\s-]?FIX|BY[\s-]?DESIGN|"
                           r"permanente(?:mente)?|NON\s+debt", re.I)),
    ("partial", re.compile(r"PARZIAL|MITIGAT[OA]|\bresidu[oa]l?\b|in\s+corso|"
                           r"\bWIP\b|\bpartial\b", re.I)),
    ("closed",  re.compile(r"CHIUS[OA]|RISOLT[OA]|\bDONE\b|✅|\bclosed\b|\bfixed\b", re.I)),
    ("open",    re.compile(r"APERT[OA]|\bPENDING\b|\bdeferred\b|\bTODO\b|"
                           r"\bopen\b|da\s+fare", re.I)),
]
STATO_LINE_RE = re.compile(r"\*\*Stat[ou]s?\*\*\s*:?\s*(.+)", re.I)
DASH_SPLIT_RE = re.compile(r"\s[—\-]\s")  # spaced em-dash / hyphen only
MD_STRIP = "~*_`# "


def classify_token(text: str) -> str | None:
    for status, rx in STATUS_RULES:
        if rx.search(text):
            return status
    return None


def classify_status(header: str, body: str, section: str) -> tuple[str, str, str]:
    """Return (status, raw_evidence, tier). Explainable; never guesses CLOSED."""
    m = STATO_LINE_RE.search(body)  # tier 1: body verdict, most reliable
    if m and (st := classify_token(m.group(1).strip())):
        return st, m.group(1).strip()[:140], "body"
    tail = DASH_SPLIT_RE.split(header)[-1]  # tier 2: header (usually last segment)
    if st := (classify_token(tail) or classify_token(header)):
        return st, header[:140], "header"
    sl = section.lower()  # tier 3: section context
    if sl.startswith("storico"):
        return "closed", section, "section"
    if "voci aperte" in sl:
        return "open", section, "section"
    if "future enhancement" in sl:
        return "wontfix", section, "section"
    if "pending" in sl or "deferred" in sl:
        return "open", section, "section"
    if "fidelity audit" in sl:
        return "info", section, "section"
    return "unknown", "", "none"


# ── source discovery ─────────────────────────────────────────────────────────
class Source:
    __slots__ = ("path", "rel", "role", "current")

    def __init__(self, path: str, role: str):
        self.path, self.role = path, role
        self.rel = os.path.relpath(path, ROOT)
        self.current = not self._archived(path)

    @staticmethod
    def _archived(path: str) -> bool:
        try:
            with open(path, encoding="utf-8") as f:
                for raw in f:
                    if raw.strip():
                        core = re.sub(r"^[^A-Za-z]+", "", raw.strip())
                        return bool(ARCHIVED_BANNER.match(core))
        except OSError:
            pass
        return False


def discover_sources() -> list[Source]:
    seen: set[str] = set()
    out: list[Source] = []
    for role, pattern in SOURCE_SPECS:
        for path in sorted(glob.glob(os.path.join(ROOT, pattern))):
            if path not in seen and os.path.isfile(path):
                seen.add(path)
                out.append(Source(path, role))
    return out


# ── parsing ──────────────────────────────────────────────────────────────────
def parse_ledger(text: str) -> list[dict]:
    """Split a ledger into `### <ID> …` entries with explainable status."""
    lines = text.splitlines()
    entries: list[dict] = []
    section = ""
    cur: dict | None = None
    buf: list[str] = []

    def flush() -> None:
        if cur is not None:
            cur["body"] = "\n".join(buf).strip()
            s, raw, src = classify_status(cur["title"], cur["body"], cur["section"])
            cur.update(status=s, status_raw=raw, status_source=src)
            entries.append(cur)

    for i, ln in enumerate(lines):
        if ln.startswith("## ") and not ln.startswith("### "):
            section = ln[3:].strip()
        if ln.startswith("### "):
            flush()
            header = ln[4:].strip()
            tok = header.split()[0] if header else f"ENTRY-{i}"
            cur = {"id": tok.strip(MD_STRIP) or tok, "title": header,
                   "section": section, "line": i + 1, "body": ""}
            buf = [ln]
        elif cur is not None:
            buf.append(ln)
    flush()
    return entries


def parse_task_file(text: str) -> list[dict]:
    """Index task headings (id-shaped) + strong-id refs with ±2 line context."""
    lines = text.splitlines()
    rows: list[dict] = []
    seen: set[tuple[str, int]] = set()

    def add(tid: str, i: int, kind: str) -> None:
        if (tid, i) not in seen:
            seen.add((tid, i))
            rows.append({"id": tid, "kind": kind, "line": i + 1,
                         "context": "\n".join(lines[max(0, i - 2): i + 3])})

    for i, ln in enumerate(lines):
        h = re.match(r"^#{1,6}\s+(\S+)", ln)
        if h:
            first = h.group(1).strip(MD_STRIP)
            if first and HEADING_ID_RE.match(first):
                add(first, i, "heading")
        for rx in STRONG_ID_RES:
            for m in rx.finditer(ln):
                add(m.group(0), i, "ref")
    return rows


# ── index build / staleness ──────────────────────────────────────────────────
def _stamp(sources: list[Source]) -> str:
    sig: list = [["schema", SCHEMA_VERSION]]
    for s in sorted(sources, key=lambda x: x.rel):
        try:
            st = os.stat(s.path)
            sig.append([s.rel, st.st_size, int(st.st_mtime)])
        except OSError:
            sig.append([s.rel, -1, -1])
    return json.dumps(sig, sort_keys=True)


def build(verbose: bool = False) -> sqlite3.Connection:
    os.makedirs(DB_DIR, exist_ok=True)
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    con = sqlite3.connect(DB_PATH)
    con.executescript(
        "CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);"
        "CREATE TABLE sources(rel TEXT, role TEXT, current INTEGER);"
        "CREATE TABLE ledger(id TEXT, title TEXT, status TEXT, status_raw TEXT,"
        " status_source TEXT, section TEXT, source TEXT, line INTEGER, body TEXT);"
        "CREATE TABLE tasks(id TEXT, kind TEXT, source TEXT, current INTEGER,"
        " line INTEGER, context TEXT);"
        "CREATE VIRTUAL TABLE fts USING fts5(role, source, ref, line UNINDEXED, body);"
    )
    sources = discover_sources()
    n_led = n_task = 0
    for s in sources:
        con.execute("INSERT INTO sources VALUES(?,?,?)", (s.rel, s.role, int(s.current)))
        try:
            with open(s.path, encoding="utf-8") as f:
                text = f.read()
        except OSError:
            continue
        if s.role == "ledger":
            for e in parse_ledger(text):
                n_led += 1
                con.execute("INSERT INTO ledger VALUES(?,?,?,?,?,?,?,?,?)",
                            (e["id"], e["title"], e["status"], e["status_raw"],
                             e["status_source"], e["section"], s.rel, e["line"], e["body"]))
                con.execute("INSERT INTO fts VALUES('ledger',?,?,?,?)",
                            (s.rel, e["id"], e["line"], e["body"]))
        else:
            for t in parse_task_file(text):
                n_task += 1
                con.execute("INSERT INTO tasks VALUES(?,?,?,?,?,?)",
                            (t["id"], t["kind"], s.rel, int(s.current), t["line"], t["context"]))
                con.execute("INSERT INTO fts VALUES('tasks',?,?,?,?)",
                            (s.rel, t["id"], t["line"], t["context"]))
    con.execute("INSERT INTO meta VALUES('schema',?)", (str(SCHEMA_VERSION),))
    con.execute("INSERT INTO meta VALUES('stamp',?)", (_stamp(sources),))
    con.commit()
    if verbose:
        print(f"indexed {len(sources)} sources → {n_led} ledger entries, {n_task} task rows")
    return con


def connect() -> sqlite3.Connection:
    if os.path.exists(DB_PATH):
        try:
            con = sqlite3.connect(DB_PATH)
            meta = dict(con.execute("SELECT key,value FROM meta").fetchall())
            if meta.get("schema") == str(SCHEMA_VERSION) and \
                    meta.get("stamp") == _stamp(discover_sources()):
                con.row_factory = sqlite3.Row
                return con
            con.close()
        except sqlite3.DatabaseError:
            pass
    con = build()
    con.row_factory = sqlite3.Row
    return con


# ── search helpers ───────────────────────────────────────────────────────────
def _fts_query(raw: str) -> str | None:
    toks = re.findall(r"\w+", raw, re.UNICODE)
    return " ".join(f'"{t}"' for t in toks) if toks else None


def _excerpt(text: str, needle: str, width: int = 80) -> str:
    flat = " ".join(text.split())
    pos = flat.lower().find(needle.lower()) if needle else -1
    if pos < 0:
        return flat[:width] + ("…" if len(flat) > width else "")
    a, b = max(0, pos - width // 2), min(len(flat), pos + len(needle) + width // 2)
    seg = flat[a:b].replace(needle, f">>{needle}<<")
    return ("…" if a else "") + seg + ("…" if b < len(flat) else "")


def _src_label(con, rel: str) -> str:
    r = con.execute("SELECT current FROM sources WHERE rel=? LIMIT 1", (rel,)).fetchone()
    return rel + ("" if (r and r[0]) else "  [archived]")


def _emit(rows) -> int:
    print(json.dumps([dict(r) for r in rows], ensure_ascii=False, indent=2))
    return 0


def _none(msg: str, as_json: bool) -> int:
    print("[]" if as_json else msg)
    return 1


# ── commands ─────────────────────────────────────────────────────────────────
def cmd_hc(con, term, as_json=False, **_):
    rows = con.execute(
        "SELECT id,status,status_source,source,line,body FROM ledger "
        "WHERE id LIKE ? OR title LIKE ? ORDER BY source,line",
        (f"%{term}%", f"%{term}%")).fetchall()
    if not rows:
        return _none(f"no hardcode entry matching '{term}'", as_json)
    if as_json:
        return _emit(rows)
    for rid, status, ssrc, src, line, body in rows:
        print(f"── {rid}  [{status.upper()} via {ssrc}]  {_src_label(con, src)}:{line}")
        print(body, end="\n\n")
    return 0


def cmd_task(con, term, as_json=False, include_all=False, **_):
    # Exact id first, then current trackers, then headings (definitions) before
    # refs — so `task A2` leads with the `### A2` entry, not `HC-…A2…` refs.
    order = ("ORDER BY (id<>?), current DESC, "
             "CASE kind WHEN 'heading' THEN 0 ELSE 1 END, source, line")
    base = "SELECT id,kind,source,current,line,context FROM tasks WHERE id LIKE ? "
    rows = con.execute(base + ("" if include_all else "AND current=1 ") + order,
                       (f"%{term}%", term)).fetchall()
    if not rows and not include_all:  # fall back to archived so CAS-L* still resolve
        rows = con.execute(base + order, (f"%{term}%", term)).fetchall()
    if not rows:
        return _none(f"no task matching '{term}'", as_json)
    if as_json:
        return _emit(rows)
    for tid, kind, src, _cur, line, ctx in rows:
        print(f"── {tid}  [{kind}]  {_src_label(con, src)}:{line}")
        print(ctx, end="\n\n")
    return 0


def cmd_search(con, query, as_json=False, **_):
    results: list[dict] = []
    keys: set[tuple[str, int]] = set()
    first = (query.split() or [query])[0]
    fq = _fts_query(query)
    if fq:
        try:
            ranked = con.execute(
                "SELECT source,ref,line,body FROM fts WHERE fts MATCH ? "
                "ORDER BY rank LIMIT 40", (fq,)).fetchall()
        except sqlite3.OperationalError:
            ranked = []  # malformed FTS expr → rely on LIKE pass below
        for r in ranked:
            if (r["source"], r["line"]) not in keys:
                keys.add((r["source"], r["line"]))
                results.append({"source": r["source"], "ref": r["ref"], "line": r["line"],
                                "snippet": _excerpt(r["body"], first), "via": "fts"})
    esc = query.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
    for r in con.execute("SELECT source,ref,line,body FROM fts "
                         "WHERE body LIKE ? ESCAPE '\\' LIMIT 40", (f"%{esc}%",)).fetchall():
        if (r["source"], r["line"]) not in keys:
            keys.add((r["source"], r["line"]))
            results.append({"source": r["source"], "ref": r["ref"], "line": r["line"],
                            "snippet": _excerpt(r["body"], query), "via": "like"})
    if not results:
        return _none(f"no match for '{query}'", as_json)
    if as_json:
        print(json.dumps(results, ensure_ascii=False, indent=2))
        return 0
    for r in results[:40]:
        print(f"[{r['via']}] {r['ref']}  {_src_label(con, r['source'])}:{r['line']}")
        print(f"   {r['snippet']}")
    return 0


def cmd_open(con, as_json=False, include_all=False, **_):
    wanted = ("open", "partial") + (("unknown",) if include_all else ())
    rows = con.execute(
        f"SELECT id,status,status_source,source,line,title FROM ledger "
        f"WHERE status IN ({','.join('?' * len(wanted))}) ORDER BY status,source,line",
        wanted).fetchall()
    if as_json:
        return _emit(rows)
    print(f"{len(rows)} actionable hardcode entries (open/partial"
          f"{'/unknown' if include_all else ''}):")
    for rid, status, ssrc, src, line, title in rows:
        print(f"  [{status:<7}] {rid}  ({_src_label(con, src)}:{line})  {title[:80]}")
    return 0


def cmd_stats(con, as_json=False, **_):
    dist = {r[0]: r[1] for r in con.execute(
        "SELECT status,COUNT(*) FROM ledger GROUP BY status").fetchall()}
    tasks = con.execute("SELECT COUNT(DISTINCT id) FROM tasks").fetchone()[0]
    cur = con.execute("SELECT COUNT(DISTINCT id) FROM tasks WHERE current=1").fetchone()[0]
    if as_json:
        print(json.dumps({"ledger_status": dist, "task_ids": tasks,
                          "task_ids_current": cur}, ensure_ascii=False, indent=2))
        return 0
    order = ["open", "partial", "unknown", "wontfix", "info", "closed"]
    parts = "  ".join(f"{k}={dist[k]}" for k in order if dist.get(k))
    print(f"hardcode ledger : {sum(dist.values())} entries — {parts}")
    print(f"  actionable (open+partial) = {dist.get('open', 0) + dist.get('partial', 0)}"
          f"  · ambiguous (unknown) = {dist.get('unknown', 0)}")
    print(f"task ids        : {tasks} distinct ({cur} in current trackers)")
    return 0


def cmd_sources(con, as_json=False, **_):
    rows = con.execute("SELECT rel,role,current FROM sources ORDER BY role,rel").fetchall()
    if as_json:
        return _emit(rows)
    for rel, role, current in rows:
        print(f"  {role:<7} {'current ' if current else 'archived'}  {rel}")
    return 0


def cmd_doctor(con, as_json=False, **_):
    specs = discover_sources()
    issues: list[str] = []
    if not any(s.role == "ledger" for s in specs):
        issues.append("no ledger source discovered")
    if not any(s.role == "tasks" and s.current for s in specs):
        issues.append("no CURRENT task tracker discovered")
    unknown = con.execute("SELECT id,source,line FROM ledger WHERE status='unknown' "
                          "ORDER BY source,line").fetchall()
    search_ok = True
    for probe in ("x2-2", "*16U", "matrix±matrix", 'a"b'):
        try:  # probe must not raise on punctuation; silence its output
            with contextlib.redirect_stdout(io.StringIO()):
                cmd_search(con, probe, as_json=True)
        except sqlite3.Error as e:  # pragma: no cover
            search_ok = False
            issues.append(f"search crashed on {probe!r}: {e}")
    if as_json:
        print(json.dumps({"sources": len(specs), "search_robust": search_ok,
                          "ledger_unknown_status": [f"{r[0]} @ {r[1]}:{r[2]}" for r in unknown],
                          "issues": issues}, ensure_ascii=False, indent=2))
        return 0 if not issues else 1
    print(f"sources discovered : {len(specs)}")
    print(f"search robustness  : {'OK' if search_ok else 'FAIL'}")
    print(f"unknown-status entries (add a **Stato**: line): {len(unknown)}")
    for r in unknown[:20]:
        print(f"    {r[0]}  ({r[1]}:{r[2]})")
    if issues:
        print("ISSUES:\n  ✗ " + "\n  ✗ ".join(issues))
        return 1
    print("doctor: OK")
    return 0


# ── cli ──────────────────────────────────────────────────────────────────────
def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="ledger_index.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.add_argument("--all", dest="include_all", action="store_true",
                   help="include archived task files / unknown-status entries")
    sub = p.add_subparsers(dest="cmd")
    sub.add_parser("build")
    for name in ("hc", "task"):
        sub.add_parser(name).add_argument("term")
    sub.add_parser("search").add_argument("term", nargs="+")
    for name in ("open", "stats", "doctor", "sources"):
        sub.add_parser(name)

    args = p.parse_args(argv)
    if not args.cmd:
        p.print_help()
        return 2
    if args.cmd == "build":
        build(verbose=True)
        return 0
    con = connect()
    kw = {"as_json": args.json, "include_all": args.include_all}
    dispatch = {"hc": cmd_hc, "task": cmd_task, "search": cmd_search, "open": cmd_open,
                "stats": cmd_stats, "doctor": cmd_doctor, "sources": cmd_sources}
    fn = dispatch[args.cmd]
    if args.cmd in ("hc", "task"):
        return fn(con, args.term, **kw)
    if args.cmd == "search":
        return fn(con, " ".join(args.term), **kw)
    return fn(con, **kw)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
