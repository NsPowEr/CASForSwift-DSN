#!/usr/bin/env python3
"""ledger_index.py — queryable index over the large markdown ledgers.

Why this exists (scaling to 4x):
    HARDCODE_LEDGER.md (~183 KB, 142 entries) and CAS_TASKS.md (~108 KB, ~199
    task ids) are linear markdown. An agent that wants ONE entry must otherwise
    load the whole file into context. At 4x project size these files dominate the
    context budget. This tool keeps the markdown as the human-editable source of
    truth and builds a disposable SQLite index (with full-text search) in
    .cache/, so agents query targeted slices instead of reading 300 KB.

Non-destructive: never writes to the markdown. The DB is rebuilt on demand and
is git-ignored (.cache/). If the markdown changes, just `build` again.

Usage:
    ledger_index.py build                 # (re)build the index from markdown
    ledger_index.py hc <id|substr>        # show hardcode ledger entry/entries
    ledger_index.py task <id|substr>      # show task lines for an id (e.g. CAS-L1-07)
    ledger_index.py search <text...>      # full-text search across both files
    ledger_index.py open                  # list OPEN hardcode entries
    ledger_index.py stats                 # counts (open/closed/tasks)

Auto-build: query commands rebuild the index automatically if the DB is missing
or older than either source file, so callers never see stale data.
"""
from __future__ import annotations

import os
import re
import sqlite3
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEDGER_MD = os.path.join(ROOT, "HARDCODE_LEDGER.md")
TASKS_MD = os.path.join(ROOT, "CAS_TASKS.md")
DB_DIR = os.path.join(ROOT, ".cache")
DB_PATH = os.path.join(DB_DIR, "ledger_index.db")

CLOSED_RE = re.compile(r"\b(CHIUS[OA]|RISOLTA|RISOLTO|DONE)\b", re.IGNORECASE)
TASK_ID_RE = re.compile(r"CAS-L\d+-\d+[A-Za-z]?")


# ── parsing ──────────────────────────────────────────────────────────────────
def parse_ledger(text: str) -> list[dict]:
    """Split HARDCODE_LEDGER.md into entries keyed by `### <ID> — ...` headers.

    Open vs closed is decided by (a) the `## Storico (risolti)` section boundary
    and (b) status keywords in the header, whichever marks it closed.
    """
    lines = text.splitlines()
    # Find the "Storico (risolti)" boundary line index, if any.
    storico_idx = next(
        (i for i, ln in enumerate(lines) if ln.strip().lower().startswith("## storico")),
        len(lines),
    )
    entries: list[dict] = []
    cur: dict | None = None
    buf: list[str] = []
    for i, ln in enumerate(lines):
        if ln.startswith("### "):
            if cur is not None:
                cur["body"] = "\n".join(buf).strip()
                entries.append(cur)
            header = ln[4:].strip()
            entry_id = header.split()[0] if header else f"ENTRY-{i}"
            in_storico = i >= storico_idx
            closed = in_storico or bool(CLOSED_RE.search(header))
            cur = {
                "id": entry_id,
                "title": header,
                "status": "CLOSED" if closed else "OPEN",
                "line": i + 1,
                "body": "",
            }
            buf = [ln]
        elif cur is not None:
            buf.append(ln)
    if cur is not None:
        cur["body"] = "\n".join(buf).strip()
        entries.append(cur)
    return entries


def parse_tasks(text: str) -> list[dict]:
    """Index every task-id occurrence with a small context window.

    CAS_TASKS.md has no per-id heading structure, so we capture each line that
    mentions a CAS-L id together with the surrounding ±2 lines, deduped per id.
    """
    lines = text.splitlines()
    rows: list[dict] = []
    seen: set[tuple[str, int]] = set()
    for i, ln in enumerate(lines):
        for m in TASK_ID_RE.finditer(ln):
            tid = m.group(0)
            key = (tid, i)
            if key in seen:
                continue
            seen.add(key)
            ctx = "\n".join(lines[max(0, i - 2): i + 3])
            rows.append({"id": tid, "line": i + 1, "context": ctx})
    return rows


# ── index build ──────────────────────────────────────────────────────────────
def build() -> sqlite3.Connection:
    os.makedirs(DB_DIR, exist_ok=True)
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
    con = sqlite3.connect(DB_PATH)
    con.executescript(
        """
        CREATE TABLE ledger(id TEXT, title TEXT, status TEXT, line INTEGER, body TEXT);
        CREATE TABLE tasks(id TEXT, line INTEGER, context TEXT);
        CREATE VIRTUAL TABLE fts USING fts5(source, ref, line, body);
        """
    )
    if os.path.exists(LEDGER_MD):
        with open(LEDGER_MD, encoding="utf-8") as f:
            for e in parse_ledger(f.read()):
                con.execute(
                    "INSERT INTO ledger VALUES(?,?,?,?,?)",
                    (e["id"], e["title"], e["status"], e["line"], e["body"]),
                )
                con.execute(
                    "INSERT INTO fts VALUES('ledger',?,?,?)",
                    (e["id"], e["line"], e["body"]),
                )
    if os.path.exists(TASKS_MD):
        with open(TASKS_MD, encoding="utf-8") as f:
            for t in parse_tasks(f.read()):
                con.execute(
                    "INSERT INTO tasks VALUES(?,?,?)",
                    (t["id"], t["line"], t["context"]),
                )
                con.execute(
                    "INSERT INTO fts VALUES('tasks',?,?,?)",
                    (t["id"], t["line"], t["context"]),
                )
    con.commit()
    return con


def _fresh(path: str) -> bool:
    if not os.path.exists(DB_PATH):
        return False
    db_m = os.path.getmtime(DB_PATH)
    return all(
        (not os.path.exists(p)) or os.path.getmtime(p) <= db_m
        for p in (LEDGER_MD, TASKS_MD, __file__)
    )


def connect() -> sqlite3.Connection:
    if _fresh(DB_PATH):
        return sqlite3.connect(DB_PATH)
    return build()


# ── commands ─────────────────────────────────────────────────────────────────
def cmd_hc(con: sqlite3.Connection, term: str) -> int:
    rows = con.execute(
        "SELECT id,status,line,body FROM ledger WHERE id LIKE ? OR title LIKE ? ORDER BY line",
        (f"%{term}%", f"%{term}%"),
    ).fetchall()
    if not rows:
        print(f"no hardcode entry matching '{term}'")
        return 1
    for rid, status, line, body in rows:
        print(f"── {rid}  [{status}]  HARDCODE_LEDGER.md:{line}")
        print(body)
        print()
    return 0


def cmd_task(con: sqlite3.Connection, term: str) -> int:
    rows = con.execute(
        "SELECT id,line,context FROM tasks WHERE id LIKE ? ORDER BY line",
        (f"%{term}%",),
    ).fetchall()
    if not rows:
        print(f"no task matching '{term}'")
        return 1
    for tid, line, ctx in rows:
        print(f"── {tid}  CAS_TASKS.md:{line}")
        print(ctx)
        print()
    return 0


def cmd_search(con: sqlite3.Connection, query: str) -> int:
    try:
        rows = con.execute(
            "SELECT source,ref,line,snippet(fts,3,'>>','<<','…',12) "
            "FROM fts WHERE fts MATCH ? LIMIT 40",
            (query,),
        ).fetchall()
    except sqlite3.OperationalError:
        rows = con.execute(
            "SELECT source,ref,line,snippet(fts,3,'>>','<<','…',12) "
            "FROM fts WHERE body LIKE ? LIMIT 40",
            (f"%{query}%",),
        ).fetchall()
    if not rows:
        print(f"no match for '{query}'")
        return 1
    for source, ref, line, snip in rows:
        src_file = "HARDCODE_LEDGER.md" if source == "ledger" else "CAS_TASKS.md"
        print(f"[{source}] {ref}  {src_file}:{line}")
        print(f"   {snip}")
    return 0


def cmd_open(con: sqlite3.Connection) -> int:
    rows = con.execute(
        "SELECT id,line,title FROM ledger WHERE status='OPEN' ORDER BY line"
    ).fetchall()
    print(f"{len(rows)} OPEN hardcode entries:")
    for rid, line, title in rows:
        print(f"  {rid}  (:{line})  {title[:90]}")
    return 0


def cmd_stats(con: sqlite3.Connection) -> int:
    opened = con.execute("SELECT COUNT(*) FROM ledger WHERE status='OPEN'").fetchone()[0]
    closed = con.execute("SELECT COUNT(*) FROM ledger WHERE status='CLOSED'").fetchone()[0]
    tasks = con.execute("SELECT COUNT(DISTINCT id) FROM tasks").fetchone()[0]
    print(f"hardcode ledger : {opened} open, {closed} closed ({opened + closed} total)")
    print(f"task ids        : {tasks} distinct")
    return 0


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    cmd, rest = argv[0], argv[1:]
    if cmd == "build":
        build()
        print(f"index built → {os.path.relpath(DB_PATH, ROOT)}")
        return 0
    con = connect()
    if cmd == "hc" and rest:
        return cmd_hc(con, rest[0])
    if cmd == "task" and rest:
        return cmd_task(con, rest[0])
    if cmd == "search" and rest:
        return cmd_search(con, " ".join(rest))
    if cmd == "open":
        return cmd_open(con)
    if cmd == "stats":
        return cmd_stats(con)
    print(__doc__)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
