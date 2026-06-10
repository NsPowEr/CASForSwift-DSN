#!/usr/bin/env bash
# F7.5.G1 — Giac binary integrity check
#
# Computes sha256 of the giac binary and verifies it matches the pinned
# manifest in scripts/giac_X.Y.Z_manifest.sha256.
#
# Workflow:
#   1. Install giac via scripts/install_giac.sh
#   2. Run scripts/giac_integrity.sh --pin  → writes manifest
#   3. Subsequent calls verify against manifest; build fails on mismatch.

set -euo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

GIAC_BIN="$(command -v giac || command -v icas || true)"
if [[ -z "${GIAC_BIN}" ]]; then
    echo "ERROR: giac/icas not on PATH. Run scripts/install_giac.sh first." >&2
    exit 1
fi

GIAC_VERSION="$("${GIAC_BIN}" --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
if [[ -z "${GIAC_VERSION}" ]]; then
    echo "ERROR: unable to determine giac version from --version output." >&2
    exit 1
fi

MANIFEST="${PROJECT_ROOT}/scripts/giac_${GIAC_VERSION}_manifest.sha256"

CURRENT_HASH="$(shasum -a 256 "${GIAC_BIN}" | awk '{print $1}')"

case "${1:-verify}" in
    --pin)
        printf '%s  %s\n' "${CURRENT_HASH}" "$(basename "${GIAC_BIN}")" > "${MANIFEST}"
        printf 'giac_version=%s\n'  "${GIAC_VERSION}" >> "${MANIFEST}"
        printf 'pinned_path=%s\n'   "${GIAC_BIN}"    >> "${MANIFEST}"
        echo "Pinned manifest written: ${MANIFEST}"
        ;;
    verify|*)
        if [[ ! -f "${MANIFEST}" ]]; then
            echo "ERROR: no manifest at ${MANIFEST}. Run with --pin first." >&2
            exit 1
        fi
        EXPECTED_HASH="$(head -1 "${MANIFEST}" | awk '{print $1}')"
        if [[ "${CURRENT_HASH}" != "${EXPECTED_HASH}" ]]; then
            echo "ERROR: giac binary sha256 mismatch!" >&2
            echo "  Expected: ${EXPECTED_HASH}" >&2
            echo "  Actual  : ${CURRENT_HASH}"   >&2
            echo "  Binary  : ${GIAC_BIN}"        >&2
            echo "Re-pin if intentional upgrade: scripts/giac_integrity.sh --pin" >&2
            exit 2
        fi
        echo "Giac integrity OK (sha256 matches pinned manifest, version ${GIAC_VERSION})."
        ;;
esac
