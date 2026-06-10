#!/usr/bin/env bash
# F7.5.G1 — Giac install scaffold
#
# Giac is GPL-3.0-or-later. We use it as a SECOND REFERENCE ORACLE,
# fork/exec only (same rule as Maxima — see CLAUDE.md §6).
#
# NO source modifications, NO embedding source/binaries in CAS Engine.
#
# Brew has no `giac` formula (only `ginac`, the C++ library — different).
# Two install paths:
#   (1) Pre-built .pkg from xcas.fr/install/ (manual download required)
#   (2) Source build from sourceforge giacxcas/giac-X.Y.Z-source.tar.gz
#
# This script: detects existing install, otherwise prints download instructions.
# Actual install is INTERACTIVE — user accepts GPL-3.0 + downloads payload.

set -euo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
GIAC_TARGET_VERSION="${GIAC_TARGET_VERSION:-1.9.0}"

color_red()   { printf '\033[31m%s\033[0m\n' "$*"; }
color_green() { printf '\033[32m%s\033[0m\n' "$*"; }
color_amber() { printf '\033[33m%s\033[0m\n' "$*"; }

# (1) Detect existing binary.
GIAC_BIN="$(command -v giac || true)"
ICAS_BIN="$(command -v icas || true)"

if [[ -n "${GIAC_BIN}" ]]; then
    INSTALLED_VER="$(giac --version 2>&1 | head -1 || echo "unknown")"
    color_green "Giac found: ${GIAC_BIN}"
    color_green "Version    : ${INSTALLED_VER}"
    echo
    echo "Next step: scripts/giac_integrity.sh to pin sha256."
    exit 0
fi

if [[ -n "${ICAS_BIN}" ]]; then
    color_green "icas (Giac CLI) found: ${ICAS_BIN}"
    color_amber "Note: use 'icas' as oracle binary; create 'giac' symlink if needed."
    exit 0
fi

# (2) Not installed → instructions.
color_red "Giac NOT installed on this system."
echo
echo "Install options (manual, user-driven — GPL-3.0-or-later):"
echo
echo "  (a) macOS .pkg installer:"
echo "      https://www-fourier.ujf-grenoble.fr/~parisse/giac/giac.dmg"
echo "      Mounts a DMG with Xcas.app; binaries land in:"
echo "        /Applications/Xcas.app/Contents/MacOS/"
echo "      Add to PATH:"
echo "        export PATH=\"/Applications/Xcas.app/Contents/MacOS:\$PATH\""
echo
echo "  (b) Source build (preferred for CI reproducibility):"
echo "      curl -O https://www-fourier.ujf-grenoble.fr/~parisse/debian/dists/stable/main/source/giac_${GIAC_TARGET_VERSION}.tar.gz"
echo "      tar xzf giac_${GIAC_TARGET_VERSION}.tar.gz"
echo "      cd giac-${GIAC_TARGET_VERSION}"
echo "      ./configure --prefix=/usr/local --disable-gui"
echo "      make -j\$(sysctl -n hw.ncpu)"
echo "      sudo make install"
echo
echo "  (c) MacPorts (if MacPorts installed):"
echo "      sudo port install giac"
echo
echo "After install, re-run: scripts/install_giac.sh"
echo "Then run: scripts/giac_integrity.sh to pin the binary sha256."
echo
color_amber "GPL-3.0 compliance: Giac is used ONLY as fork/exec oracle. NO source"
color_amber "embedding, NO derivative work, NO redistribution of Giac binaries"
color_amber "with CAS Engine. See CLAUDE.md §6 (Maxima oracle rule applies symmetrically)."

exit 2
