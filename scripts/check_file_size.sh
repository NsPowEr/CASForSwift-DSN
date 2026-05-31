#!/usr/bin/env bash
# F0.6 Anti-monolith CI gate.
#
# Scans src/ and include/ for .cpp/.hpp files exceeding MAX_LINES (default 500).
# Files listed in scripts/file_size_whitelist.txt are exempt (pre-existing debt
# with mandatory split tickets).
#
# Usage:
#   scripts/check_file_size.sh [--max-lines 500] [--whitelist path] [--verbose]
#
# Exit codes:
#   0 — all non-whitelisted files are within the line limit
#   1 — one or more files violate the limit

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

MAX_LINES=500
WHITELIST_FILE="$script_dir/file_size_whitelist.txt"
VERBOSE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-lines)  MAX_LINES="$2";       shift 2 ;;
        --whitelist)  WHITELIST_FILE="$2";  shift 2 ;;
        --verbose)    VERBOSE=1;            shift   ;;
        *) echo "Unknown argument: $1" >&2; exit 1  ;;
    esac
done

# Delegate to Python3 (bash 3.2 on macOS lacks associative arrays)
python3 - "$repo_root" "$WHITELIST_FILE" "$MAX_LINES" "$VERBOSE" << 'PYEOF'
import sys
import os
import subprocess

repo_root      = sys.argv[1]
whitelist_file = sys.argv[2]
max_lines      = int(sys.argv[3])
verbose        = sys.argv[4] == "1"

# Load whitelist (repo-relative paths)
whitelist = set()
if os.path.isfile(whitelist_file):
    with open(whitelist_file) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if line:
                whitelist.add(line)

src_dir     = os.path.join(repo_root, "src")
include_dir = os.path.join(repo_root, "include")

violations = 0
checked    = 0
skipped    = 0
all_files  = []

for base in [src_dir, include_dir]:
    if not os.path.isdir(base):
        continue
    for dirpath, _dirs, files in os.walk(base):
        for fname in files:
            if fname.endswith(".cpp") or fname.endswith(".hpp"):
                all_files.append(os.path.join(dirpath, fname))

all_files.sort()

for abs_path in all_files:
    rel_path = os.path.relpath(abs_path, repo_root)
    # Skip CMake generated files
    if "_deps" in rel_path or "CMakeFiles" in rel_path:
        continue

    with open(abs_path, "rb") as f:
        line_count = sum(1 for _ in f)

    checked += 1

    if rel_path in whitelist:
        if verbose:
            print(f"EXEMPT  {line_count:5d}  {rel_path}")
        skipped += 1
        continue

    if line_count > max_lines:
        print(f"FAIL    {line_count:5d}  {rel_path}  (limit: {max_lines})")
        violations += 1
    else:
        if verbose:
            print(f"OK      {line_count:5d}  {rel_path}")

print(f"---")
print(f"Anti-monolith scan: {checked} files checked, {skipped} whitelisted, {violations} violation(s) (limit: {max_lines} lines).")

if violations > 0:
    print(f"FAIL: {violations} file(s) exceed {max_lines} lines. Split required per CLAUDE.md §ANTI-MONOLITH.", file=sys.stderr)
    print(f"  Add a split task to CAS_TASKS.md and (if genuinely unavoidable right now)", file=sys.stderr)
    print(f"  add to scripts/file_size_whitelist.txt with a mandatory split ticket reference.", file=sys.stderr)
    sys.exit(1)

print("OK: no new anti-monolith violations.")
sys.exit(0)
PYEOF
