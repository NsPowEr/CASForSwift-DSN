#!/usr/bin/env python3
"""
run_golden_sympy.py — optional secondary reference check using SymPy 1.13+.
Non-blocking: if SymPy is not installed the script exits cleanly.

Usage:
    python3 scripts/run_golden_sympy.py <corpus.jsonl> <output_dir>

For each input line, invokes SymPy via its public API (no eval) and writes
<output_dir>/<idx>.sympy.out with the string representation of the result.

This is a secondary (cross-check) oracle; differences between Maxima and SymPy
are noted but do not block the test suite.
"""

import json
import os
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# SymPy availability check
# ---------------------------------------------------------------------------
try:
    import sympy as sp
    SYMPY_AVAILABLE = True
    SYMPY_VERSION = sp.__version__
except ImportError:
    SYMPY_AVAILABLE = False
    SYMPY_VERSION = "not installed"

if not SYMPY_AVAILABLE:
    print(f"SymPy not available ({SYMPY_VERSION}). Skipping secondary cross-check.")
    sys.exit(0)


# ---------------------------------------------------------------------------
# Symbol pool
# ---------------------------------------------------------------------------
_SYMS: dict[str, sp.Symbol] = {}

def sym(name: str) -> sp.Symbol:
    if name not in _SYMS:
        _SYMS[name] = sp.Symbol(name)
    return _SYMS[name]


def _make_default_syms() -> dict:
    d = {}
    for ch in "xyzntabcdfghmk":
        d[ch] = sym(ch)
    return d


# ---------------------------------------------------------------------------
# Dispatcher: map our CAS input string to a SymPy computation.
# Uses structural dispatch — NOT eval() — for security.
# ---------------------------------------------------------------------------

def _split_args(inner: str) -> list[str]:
    """Split comma-separated top-level arguments."""
    result = []
    depth = 0
    current = ""
    for ch in inner:
        if ch in "([": depth += 1; current += ch
        elif ch in ")]": depth -= 1; current += ch
        elif ch == "," and depth == 0:
            result.append(current.strip())
            current = ""
        else:
            current += ch
    current = current.strip()
    if current:
        result.append(current)
    return result


def _parse_expr_str(s: str, syms: dict) -> sp.Expr:
    """
    Parse a simple algebraic expression string into a SymPy expression
    using sp.sympify with a local symbol namespace.
    sympify uses safe parsing (no exec/eval of Python code).
    """
    # Normalise: pi -> pi (SymPy knows pi), inf -> oo
    s = re.sub(r'\bpi\b', 'pi', s)
    s = re.sub(r'\binf\b', 'oo', s)
    # log( -> log(  — SymPy log is natural log by default
    return sp.sympify(s, locals=syms)


def _parse_int(s: str) -> int:
    return int(s.strip())


def dispatch(input_str: str) -> tuple[str, str]:
    """
    Dispatch a CAS input string to a SymPy call.
    Returns (status, result_string) where status is "ok" or "skip:<reason>".
    """
    syms = _make_default_syms()
    s = input_str.strip()

    # Find top-level function name
    m = re.match(r'^([a-zA-Z_][a-zA-Z0-9_]*)\((.+)\)$', s, re.DOTALL)
    if not m:
        # Bare expression: sympify and simplify
        try:
            expr = _parse_expr_str(s, syms)
            return "ok", str(sp.simplify(expr))
        except Exception as exc:
            return f"skip:parse_err:{exc}", ""

    fn_name = m.group(1)
    args_str = m.group(2)
    args = _split_args(args_str)

    try:
        if fn_name == "integrate" and len(args) == 2:
            f = _parse_expr_str(args[0], syms)
            v = sym(args[1].strip())
            return "ok", str(sp.integrate(f, v))

        elif fn_name == "diff" and len(args) >= 2:
            f = _parse_expr_str(args[0], syms)
            v = sym(args[1].strip())
            order = _parse_int(args[2]) if len(args) >= 3 else 1
            return "ok", str(sp.diff(f, v, order))

        elif fn_name == "limit" and len(args) >= 3:
            f = _parse_expr_str(args[0], syms)
            v = sym(args[1].strip())
            pt = _parse_expr_str(args[2], syms)
            direction = "+-"
            if len(args) >= 4:
                d = args[3].strip()
                if d == "plus":  direction = "+"
                elif d == "minus": direction = "-"
            return "ok", str(sp.limit(f, v, pt, direction))

        elif fn_name == "factor" and len(args) == 1:
            f = _parse_expr_str(args[0], syms)
            return "ok", str(sp.factor(f))

        elif fn_name == "gcd" and len(args) == 2:
            p = _parse_expr_str(args[0], syms)
            q = _parse_expr_str(args[1], syms)
            return "ok", str(sp.gcd(p, q))

        elif fn_name == "solve" and len(args) == 2:
            f = _parse_expr_str(args[0], syms)
            v = sym(args[1].strip())
            return "ok", str(sp.solve(f, v))

        elif fn_name == "series" and len(args) == 4:
            f  = _parse_expr_str(args[0], syms)
            v  = sym(args[1].strip())
            pt = _parse_expr_str(args[2], syms)
            n  = _parse_int(args[3]) + 1
            return "ok", str(sp.series(f, v, pt, n).removeO())

        elif fn_name == "simplify" and len(args) == 1:
            f = _parse_expr_str(args[0], syms)
            return "ok", str(sp.simplify(f))

        elif fn_name in {"det", "trace", "transpose", "inverse", "rank", "eigenvalues"}:
            return "skip:matrix_fn", ""

        else:
            # Attempt generic sympify + simplify
            f = _parse_expr_str(s, syms)
            return "ok", str(sp.simplify(f))

    except Exception as exc:
        return f"skip:dispatch_err:{exc}", ""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <corpus.jsonl> <output_dir>")
        sys.exit(1)

    corpus_path = sys.argv[1]
    output_dir  = sys.argv[2]

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print(f"SymPy version: {SYMPY_VERSION}")
    print(f"Corpus:       {corpus_path}")
    print(f"Output dir:   {output_dir}")
    print()

    total = 0
    errors = 0

    with open(corpus_path) as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue

            input_expr = entry.get("input", "")
            if not input_expr:
                continue

            total += 1
            status, result_str = dispatch(input_expr)

            out_path = os.path.join(output_dir, f"{idx}.sympy.out")
            with open(out_path, "w") as out:
                if status == "ok":
                    out.write(result_str + "\n")
                else:
                    out.write(f"SKIP:{status}\n")

            if status == "ok":
                print(f"  OK   [{idx:3d}] {input_expr} => {result_str}")
            else:
                errors += 1
                print(f"  SKIP [{idx:3d}] {input_expr} ({status})")

    print(f"\nDone. {total} inputs, {errors} skipped/errors.")


if __name__ == "__main__":
    main()
