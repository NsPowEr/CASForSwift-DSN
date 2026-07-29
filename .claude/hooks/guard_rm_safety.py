#!/usr/bin/env python3
# Logic for the PreToolUse "REGOLA EVIDENCE-FIRST" deletion guard, invoked by
# guard_rm_safety.sh. Philosophy (CLAUDE.md):
#
#   - regenerable build output (build*/, .cache/, untracked root stragglers)
#     stays freely deletable — deleting it is lossless by construction;
#   - git-TRACKED paths and protected source/infra dirs are denied: deletion
#     there must go through `git rm`/`git mv` (history keeps a copy) or the
#     attic protocol, with evidence collected first;
#   - absolute paths outside the repo are denied except temp/attic locations
#     (the agent has no business deleting the user's files, and the oracle
#     installs live outside the repo).
#
# The human can always run the command directly in a terminal — this only
# gates the agent's Bash tool. Fires only when rm/find is the *executed*
# command of a pipeline segment, so mentions inside grep/echo do not trip it.
import json
import os
import re
import shlex
import subprocess
import sys

try:
    data = json.load(sys.stdin)
except Exception:
    sys.exit(0)

cmd = (data.get("tool_input", {}) or {}).get("command", "") or ""
if not cmd.strip():
    sys.exit(0)
# Fast path: nothing deletion-shaped anywhere in the command.
if not re.search(r"\b(rm|find)\b", cmd):
    sys.exit(0)

HOOK_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HOOK_DIR))

# First path component (repo-relative) → protected. Everything an agent could
# not regenerate: sources, tests, docs, rules, references, agent config, git.
PROTECTED_ROOTS = {
    "src", "include", "test", "docs", "scripts", "tools",
    ".APROJECT_REFERENCES", ".claude", ".git",
    "!_VECCHI DOCUMENTI",
}
# Absolute prefixes where deletion is fine (temp, scratchpads, attic).
_tmp = os.environ.get("TMPDIR", "/tmp")
SAFE_ABS_PREFIXES = tuple(
    os.path.realpath(p) for p in ("/tmp", "/private/tmp", _tmp)
) + (os.path.expanduser("~/cas-attic"),)

segments = re.split(r"&&|\|\||;|\||\n", cmd)


def tokenize(seg):
    try:
        return shlex.split(seg, comments=False, posix=True)
    except ValueError:
        return seg.split()


def first_exec(seg):
    toks = tokenize(seg)
    i = 0
    while i < len(toks) and re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", toks[i]):
        i += 1  # leading env assignment
    if i < len(toks) and os.path.basename(toks[i]) in ("timeout", "gtimeout"):
        i += 1
        while i < len(toks) and (toks[i].startswith("-")
                                 or re.match(r"^[0-9]+[smhd]?$", toks[i])):
            i += 1
    if i < len(toks):
        return os.path.basename(toks[i]), toks[i:]
    return None, []


def deny(target, why):
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason":
                "REGOLA EVIDENCE-FIRST (CLAUDE.md): cancellazione bloccata — "
                f"{target!r}: {why}\n"
                "Percorsi ammessi: (a) output rigenerabile (build*/, .cache/) "
                "→ rm libero; (b) file TRACKED → `git rm`/`git mv` con prove "
                "(refs grep, rigenerabilità, reversibilità); (c) untracked "
                "unico → attic `~/cas-attic-<data>/`, mai rm. Se l'azione è "
                "davvero necessaria la esegue l'utente a mano.",
        }
    }))
    sys.exit(0)


def is_tracked(rel):
    try:
        out = subprocess.run(
            ["git", "ls-files", "--", rel], cwd=REPO,
            capture_output=True, text=True, timeout=5)
        return bool(out.stdout.strip())
    except Exception:
        return False  # never wedge deletion of clearly-unprotected paths


def check_target(tok):
    t = os.path.expanduser(tok)
    if os.path.isabs(t):
        rt = os.path.realpath(t)
        if any(rt == p or rt.startswith(p + os.sep) for p in SAFE_ABS_PREFIXES):
            return
        repo_real = os.path.realpath(REPO)
        if rt == repo_real:
            deny(tok, "è la root del repository")
        if rt.startswith(repo_real + os.sep):
            rel = os.path.relpath(rt, repo_real)
        else:
            deny(tok, "path assoluto fuori dal repository (non temp/attic)")
    else:
        rel = os.path.normpath(t)
    if rel in (".", "..") or rel.startswith(".." + os.sep):
        deny(tok, "esce dallo scope del repository")
    first = rel.split(os.sep, 1)[0]
    if first in PROTECTED_ROOTS:
        deny(tok, f"directory protetta `{first}/` (sorgenti/infra non rigenerabili)")
    if is_tracked(rel):
        deny(tok, "path tracciato da git (usa `git rm`/`git mv`, la storia recupera)")


for seg in segments:
    base, toks = first_exec(seg)
    if base == "rm":
        for tok in toks[1:]:
            if tok.startswith("-") or tok == "--":
                continue
            check_target(tok)
    elif base == "find":
        rest = toks[1:]
        destructive = "-delete" in rest or (
            "-exec" in rest and any(os.path.basename(t) == "rm" for t in rest))
        if not destructive:
            continue
        for tok in rest:
            if tok.startswith("-"):
                break  # find roots come before the first predicate
            check_target(tok)

sys.exit(0)
