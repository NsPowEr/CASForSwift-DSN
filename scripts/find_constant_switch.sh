#!/usr/bin/env bash
# F7.5.F1 — Scan tutti i siti switch (Constant::value) o switch su MathConstant.
# Ogni switch deve gestire NegInfinity/ComplexInfinity/Indeterminate/NaN
# esplicitamente o ritornare Indeterminate diagnostico (REGOLA 0.2).
#
# Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Extended_Real_AST.md
# Ledger: HC-F70-A43-EXTENDED-REAL.

set -euo pipefail

ROOT="${1:-.}"

echo "=== switch su MathConstant ==="
grep -rnE 'switch\s*\(.*->value\)|switch\s*\(.*\.value\)|switch\s*\(.*MathConstant.*\)' \
  --include='*.cpp' --include='*.hpp' "${ROOT}/src" "${ROOT}/include" \
  2>/dev/null | sort -u || true

echo
echo "=== case MathConstant:: ==="
grep -rnE 'case MathConstant::' \
  --include='*.cpp' --include='*.hpp' "${ROOT}/src" "${ROOT}/include" \
  2>/dev/null | sort -u || true

echo
echo "=== case bare (Pi/E/I/Infinity senza prefisso enum) ==="
grep -rnE '^\s*case (Pi|E|I|Infinity|NegInfinity|ComplexInfinity|Indeterminate|NaN|EulerGamma):' \
  --include='*.cpp' --include='*.hpp' "${ROOT}/src" "${ROOT}/include" \
  2>/dev/null | sort -u || true
