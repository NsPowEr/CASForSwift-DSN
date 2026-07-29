#!/usr/bin/env bash
# PreToolUse Bash guard: mechanize the CLAUDE.md "Git safety" rule.
#
# The constitution (PROTOCOLLO ANTI-LOOP → Git safety) says: only `git stash
# push` for backups; NO `git reset --hard`, NO `git restore --source`. That
# rule was prose-only — this hook gives it teeth, same philosophy as
# guard_anti_loop.sh. It also denies the close cousins that destroy work the
# same way (worktree-discarding restore/checkout, clean -f, stash drop/clear,
# force push) because an agent reaching for any of them is in the same failure
# mode: throwing away state instead of preserving it.
#
# The human can always run these commands directly in a terminal — the hook
# only gates the agent's Bash tool.
#
# Reads tool_input JSON from stdin (Claude Code hook protocol). Never throws.

set -euo pipefail

CMD=$(python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('tool_input',{}).get('command',''))" 2>/dev/null || echo "")
[[ -z "$CMD" ]] && exit 0

# Fast path: not a git command at all.
echo "$CMD" | grep -q "git" || exit 0

deny() {
    local reason="$1"
    REASON="$reason" python3 -c '
import json, os
print(json.dumps({
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "deny",
    "permissionDecisionReason":
      "GIT SAFETY (CLAUDE.md): " + os.environ["REASON"] + "\n"
      "Backup consentito: SOLO `git stash push`. Se il comando è davvero "
      "necessario, deve eseguirlo l’utente a mano nel proprio terminale."
  }
}))'
    exit 0
}

# Normalize whitespace for matching; keep original for messages.
NORM=$(echo "$CMD" | tr -s '[:space:]' ' ')

# Patterns are prefix-free (no leading "git ") so `git -C <dir> …`,
# `/usr/bin/git …` and `command git …` are covered too; the grep gate above
# ensures we only inspect commands that mention git at all. A rare false
# positive (e.g. echoing one of these strings in a git-mentioning command) is
# an acceptable price: the deny message explains itself and the agent rephrases.
case "$NORM" in
    *"reset --hard"*|*"reset --merge"*)
        deny "\`git reset --hard/--merge\` distrugge il working tree. Usa \`git stash push\` per salvare lo stato, \`git reset --soft\`/\`--mixed\` per muovere HEAD." ;;
    *"restore --source"*|*"restore --worktree"*)
        deny "\`git restore --source/--worktree\` sovrascrive file non committati (vietato esplicitamente). Leggi il file storico con \`git show <ref>:<path>\` e applica le differenze a mano." ;;
    *"checkout -- "*|*"checkout ."*|*"checkout -f"*|*"checkout --force"*)
        deny "\`git checkout -- <path>\`/\`-f\` scarta modifiche non committate in silenzio. Se serve davvero, prima \`git stash push\`." ;;
    *"clean -f"*|*"clean -xf"*|*"clean -fx"*|*"clean -fd"*|*"clean -df"*|*"clean --force"*)
        deny "\`git clean -f*\` cancella file non tracciati (possibile lavoro in corso non ancora aggiunto)." ;;
    *"stash drop"*|*"stash clear"*)
        deny "\`git stash drop/clear\` distrugge i backup: lo stash è la rete di sicurezza del protocollo anti-loop." ;;
    *"push --force"*|*"push -f "*|*"push -f")
        deny "force-push riscrive storia remota. Decisione dell'utente, mai dell'agente." ;;
    *"branch -D"*|*"branch --delete --force"*)
        deny "\`git branch -D\` cancella un branch anche se non merged (lavoro potenzialmente perso). Usa \`git branch -d\` (safe) o lascia decidere l'utente." ;;
    *"git rm -f"*|*"git rm -rf"*|*"git rm --force"*)
        deny "\`git rm -f\` scarta anche le modifiche non committate del file. Usa \`git rm\` semplice (rifiuta se dirty) dopo aver raccolto le prove (REGOLA EVIDENCE-FIRST)." ;;
    *"worktree remove --force"*|*"worktree remove -f"*)
        deny "\`git worktree remove --force\` cancella un worktree anche se dirty." ;;
esac

exit 0
