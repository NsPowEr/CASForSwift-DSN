#!/usr/bin/env bash
# tasks_audit.sh — Unified parser for all task list files in the repo.
#
# Reads the canonical task lists, normalizes their heterogeneous formats
# (Markdown checkboxes, pipe tables, ad-hoc ID lists) into a single row
# stream:  <file>|<format>|<id>|<status>|<section>
#
# Status is normalized to a controlled vocabulary:
#   DONE, PARTIAL, OPEN, PENDING, UNKNOWN
#
# Recognized inputs (all optional — missing files silently skipped):
#   TODO.md                       — checkbox only
#   TODO_PH8.md                   — checkbox + F8.x/CAS-L/HC-F ids
#   CAS_TASKS.md                  — pipe table with CAS-L-* ids
#   PLAN_HP_PRIME_PARITY.md       — pipe table with CAS-L-*/HC-F* ids
#   PLAN_TASKS_REMAINING.md       — pipe table with Task N + ad-hoc status
#   PLAN_NEXT_SESSIONS.md         — checkbox + F8.x/W*.x + pipe table
#   PLAN_F3_F8_GAP_CLOSURE.md     — pipe table with CAS-L-*/HC-F*/F* ids
#
# Output modes:
#   (default)    human-readable summary
#   --tsv        tab-separated rows for piping
#   --json       JSON object for tooling
#
# Exit code: 0 always (informational).

set -u
shopt -s nullglob

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

MODE="${1:-summary}"

# ─── Status normalization map ────────────────────────────────────────────────
# Maps every observed status string to one of: DONE, PARTIAL, OPEN, PENDING, UNKNOWN
norm_status() {
  local s="${1:-}"
  # lowercase + trim (POSIX-safe: tr instead of ${s,,})
  s="$(printf '%s' "$s" | tr '[:upper:]' '[:lower:]')"
  s="${s## }"; s="${s%% }"
  # Prefix/substring matching: real ledgers append parenthetical scope notes,
  # e.g. "Parziale (solo Z, no PID)" or "Risolta (ComplexLit integrato)".
  case "$s" in
    # DONE — note: "parziale" is checked FIRST below, so "risolta" here is safe
    risolta*|done*|"✓ done"*|closed*|complete*|completed*|verified*|verificata*|"[x]"|x|"✅"*) echo "DONE" ;;
    # PARTIAL
    parziale*|partial*|"⚠ partial"*|mvp*|"in progress"*|in_progress*|wip*) echo "PARTIAL" ;;
    # OPEN
    aperta*|open*|blocked*) echo "OPEN" ;;
    # PENDING
    pending*|todo*|planned*|queued*) echo "PENDING" ;;
    # Unknown
    *) echo "UNKNOWN" ;;
  esac
}

# ─── Per-file extractors ─────────────────────────────────────────────────────
# Each extractor writes rows to the global ROWS array in the form:
#   file<TAB>format<TAB>id<TAB>status<TAB>section

declare -a ROWS=()

# TODO.md — flat checkbox list
extract_todo_md() {
  local file="$1" section="(root)" line id status
  while IFS= read -r line; do
    # heading → section
    if [[ "$line" =~ ^##+\ +(.+)$ ]]; then
      section="${BASH_REMATCH[1]}"
      continue
    fi
    # checkbox: "- [x] text" or "- [ ] text"
    if [[ "$line" =~ ^[[:space:]]*-\ \[([x ])\]\ *(.*)$ ]]; then
      local mark="${BASH_REMATCH[1]}" text="${BASH_REMATCH[2]}"
      if [[ "$mark" == "x" ]]; then status="DONE"; else status="PENDING"; fi
      id="$text"
      ROWS+=("${file}"$'	'"checkbox"$'	'"$id"$'	'"$status"$'	'"$section")
    fi
  done < "$file"
}

# TODO_PH8.md — same checkbox parser (also has F8.x/CAS-L/HC-F ids inline)
extract_todo_ph8() {
  extract_todo_md "$1"
}

# PLAN_NEXT_SESSIONS.md — mixed checkbox + pipe table
extract_plan_next_sessions() {
  local file="$1" line section="(root)"
  while IFS= read -r line; do
    if [[ "$line" =~ ^##+\ +(.+)$ ]]; then
      section="${BASH_REMATCH[1]}"; continue
    fi
    # checkbox
    if [[ "$line" =~ ^[[:space:]]*-\ \[([x ])\]\ *(.*)$ ]]; then
      local mark="${BASH_REMATCH[1]}" text="${BASH_REMATCH[2]}"
      local status="PENDING"; [[ "$mark" == "x" ]] && status="DONE"
      ROWS+=("${file}"$'	'"checkbox"$'	'"$text"$'	'"$status"$'	'"$section")
      continue
    fi
    # pipe-table status row: | 4 | ✓ DONE | note |
    if [[ "$line" =~ ^\|([^\|]+)\|([^\|]+)\|([^\|]+)\| ]]; then
      local id="${BASH_REMATCH[1]// /}" st="${BASH_REMATCH[2]}" note="${BASH_REMATCH[3]}"
      st="${st## }"; st="${st%% }"
      # skip header / separator
      [[ "$id" == "Task" || "$id" == "----" || "$id" == ---* ]] && continue
      ROWS+=("${file}"$'	'"table"$'	'"$id"$'	'"$(norm_status "$st")"$'	'"$note")
    fi
  done < "$file"
}

# PLAN_TASKS_REMAINING.md — pipe table with Task N + ad-hoc status
extract_plan_tasks_remaining() {
  local file="$1"
  local line
  while IFS= read -r line; do
    if [[ "$line" =~ ^\|([^\|]+)\|([^\|]+)\|([^\|]+)\| ]]; then
      local id="${BASH_REMATCH[1]// /}" st="${BASH_REMATCH[2]}" note="${BASH_REMATCH[3]}"
      st="${st## }"; st="${st%% }"
      [[ "$id" == "Task" || "$id" == ---* ]] && continue
      ROWS+=("${file}"$'	'"table"$'	'"$id"$'	'"$(norm_status "$st")"$'	'"$note")
    fi
  done < "$file"
}

# CAS_TASKS.md / PLAN_HP_PRIME_PARITY.md / PLAN_F3_F8_GAP_CLOSURE.md
# — pipe table with id in col 2, status in col 5
# Schema (verified empirically): | ID | Domain | Name | Level | Status | Note | Priority | Verification |
extract_pipe_id_status() {
  local file="$1" col_id="${2:-2}" col_status="${3:-5}"
  local line
  while IFS= read -r line; do
    [[ "$line" != \|* ]] && continue
    IFS='|' read -ra parts <<< "$line"
    # bash pads with empty leading/trailing entry
    local n=${#parts[@]}
    [[ $n -lt $(( col_status + 1 )) ]] && continue
    local id="${parts[$col_id]}"  st="${parts[$col_status]}"
    id="${id## }"; id="${id%% }"
    st="${st## }"; st="${st%% }"
    [[ -z "$id" || "$id" == "ID" || "$id" == ---* ]] && continue
    ROWS+=("${file}"$'	'"table"$'	'"$id"$'	'"$(norm_status "$st")"$'	'"")
  done < "$file"
}

# Narrative scan — files where tasks live in prose with inline IDs + status
# words (PLAN_HP_PRIME_PARITY.md, PLAN_F3_F8_GAP_CLOSURE.md). Extracts every
# recognized task ID and infers status from status keywords on the same line.
# IDs recognized: HPP-xxx, HC-xxx, CAS-Lx-yy, Fn.m, Wn.m
extract_narrative_ids() {
  local file="$1" line section="(root)"
  while IFS= read -r line; do
    if [[ "$line" =~ ^##+\ +(.+)$ ]]; then
      section="${BASH_REMATCH[1]}"; continue
    fi
    # find IDs on this line
    local ids
    ids=$(printf '%s' "$line" | grep -oE '(HPP-[A-Za-z0-9.]+|HC-[A-Za-z0-9.-]+|CAS-L[0-9]+-[0-9]+|\bF[0-9]+\.[0-9]+[A-Za-z]?|\bW[0-9]+\.[0-9]+)' | sort -u)
    [[ -z "$ids" ]] && continue
    # infer status from line keywords (case-insensitive)
    local lc st="UNKNOWN"
    lc="$(printf '%s' "$line" | tr '[:upper:]' '[:lower:]')"
    if   [[ "$lc" == *aperta*permanente* || "$lc" == *aperte*permanenti* ]]; then st="OPEN"
    elif [[ "$lc" == *aperta* || "$lc" == *open* ]]; then st="OPEN"
    elif [[ "$lc" == *parziale* || "$lc" == *partial* ]]; then st="PARTIAL"
    elif [[ "$lc" == *risolta* || "$lc" == *done* || "$lc" == *chiusa* || "$lc" == *closed* ]]; then st="DONE"
    fi
    # Narrative scan only emits rows when a status keyword was found on the
    # line. Lines that merely mention an ID without status (cross-references,
    # coverage tables) would otherwise flood the report with UNKNOWN noise.
    [[ "$st" == "UNKNOWN" ]] && continue
    local id
    while IFS= read -r id; do
      [[ -z "$id" ]] && continue
      ROWS+=("${file}"$'	'"narrative"$'	'"$id"$'	'"$st"$'	'"$section")
    done <<< "$ids"
  done < "$file"
}

# ─── Dispatch ────────────────────────────────────────────────────────────────
[[ -f TODO.md ]]                  && extract_todo_md                TODO.md
[[ -f TODO_PH8.md ]]              && extract_todo_ph8              TODO_PH8.md
# CAS_TASKS schema: | <empty> | ID | Domain | Name | Level | Status | ... |
# After IFS='|' split: parts[1]=ID, parts[5]=Status
[[ -f CAS_TASKS.md ]]             && extract_pipe_id_status        CAS_TASKS.md 1 5
[[ -f PLAN_TASKS_REMAINING.md ]]  && extract_plan_tasks_remaining  PLAN_TASKS_REMAINING.md
[[ -f PLAN_NEXT_SESSIONS.md ]]    && extract_plan_next_sessions    PLAN_NEXT_SESSIONS.md
[[ -f PLAN_HP_PRIME_PARITY.md ]]  && extract_narrative_ids         PLAN_HP_PRIME_PARITY.md
[[ -f PLAN_F3_F8_GAP_CLOSURE.md ]]&& extract_narrative_ids         PLAN_F3_F8_GAP_CLOSURE.md

# ─── Output ──────────────────────────────────────────────────────────────────
case "$MODE" in
  --tsv)
    printf "file\tformat\tid\tstatus\tnote\n"
    for r in "${ROWS[@]}"; do
      printf "%s\n" "$r"
    done
    ;;
  --json)
    # Use python3 to emit a clean JSON object.
    python3 - "${ROWS[@]}" <<'PY'
import sys, json
RANK = {"OPEN":0, "PENDING":1, "UNKNOWN":2, "PARTIAL":3, "DONE":4}
raw = []
for line in sys.argv[1:]:
    parts = line.split('\t', 4)
    if len(parts) < 5: parts += [''] * (5 - len(parts))
    raw.append({
        "file": parts[0],
        "format": parts[1],
        "id": parts[2],
        "status": parts[3],
        "note": parts[4],
    })
best = {}
for r in raw:
    key = (r["file"], r["id"])
    if key not in best or RANK[r["status"]] < RANK[best[key]["status"]]:
        best[key] = r
rows = list(best.values())
# Aggregate
from collections import Counter
status_total = Counter(r["status"] for r in rows)
file_total   = Counter(r["file"] for r in rows)
print(json.dumps({
    "total_rows": len(rows),
    "by_status": dict(status_total),
    "by_file":   dict(file_total),
    "rows":      rows,
}, indent=2, ensure_ascii=False))
PY
    ;;
  *)
    # Human-readable summary
    python3 - "${ROWS[@]}" <<'PY'
import sys
from collections import Counter, defaultdict
# Status precedence: when the same (file,id) appears multiple times across
# narrative lines, keep the WORST (most actionable) status.
RANK = {"OPEN":0, "PENDING":1, "UNKNOWN":2, "PARTIAL":3, "DONE":4}
raw = []
for line in sys.argv[1:]:
    parts = line.split('\t', 4)
    if len(parts) < 5: parts += [''] * (5 - len(parts))
    raw.append({
        "file": parts[0],
        "format": parts[1],
        "id": parts[2],
        "status": parts[3],
        "note": parts[4].strip(),
    })
# Dedup by (file, id) keeping worst status
best = {}
for r in raw:
    key = (r["file"], r["id"])
    if key not in best or RANK[r["status"]] < RANK[best[key]["status"]]:
        best[key] = r
rows = list(best.values())
by_file_status = defaultdict(Counter)
for r in rows:
    by_file_status[r["file"]][r["status"]] += 1
overall = Counter(r["status"] for r in rows)
print("═══════════════════════════════════════════════════════════")
print("  CAS Engine — Tasks Audit (unified across all ledgers)")
print(f"  total rows parsed: {len(rows)}")
print("═══════════════════════════════════════════════════════════")
print(f"{'STATUS':<10} {'COUNT':>6}")
print("─" * 18)
for s in ("DONE","PARTIAL","OPEN","PENDING","UNKNOWN"):
    print(f"{s:<10} {overall.get(s,0):>6}")
print("─" * 18)
print()
print(f"{'FILE':<32} {'DONE':>5} {'PART':>5} {'OPEN':>5} {'PEND':>5} {'UNK':>5} {'TOT':>5}")
print("─" * 62)
for f in sorted(by_file_status):
    c = by_file_status[f]
    print(f"{f:<32} {c.get('DONE',0):>5} {c.get('PARTIAL',0):>5} "
          f"{c.get('OPEN',0):>5} {c.get('PENDING',0):>5} "
          f"{c.get('UNKNOWN',0):>5} {sum(c.values()):>5}")
print("─" * 62)
# Open + unknown items: actionable
actionable = [r for r in rows if r["status"] in ("OPEN","UNKNOWN","PENDING")]
if actionable:
    print()
    print(f"⚠ Actionable items (OPEN/UNKNOWN/PENDING): {len(actionable)}")
    # Cap to 30 most relevant
    for r in sorted(actionable, key=lambda x:(x["file"],x["id"]))[:30]:
        note = f" — {r['note']}" if r["note"] else ""
        print(f"  [{r['status']:8}] {r['file']:30} {r['id']}{note}")
    if len(actionable) > 30:
        print(f"  ... and {len(actionable)-30} more (use --tsv for full list)")
PY
    ;;
esac
