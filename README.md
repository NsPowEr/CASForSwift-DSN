# REAL CAS ENGINE C++

Industrial-grade computer algebra system in modern C++.

Current focus: calculus foundations (`F4`) in pure C++, on top of completed repository/parser refactors (`F0`, `F1`, `F1b`) and a partially completed symbolic core (`F2`).
Swift/C API boundary is explicitly suspended until the mathematical engine is complete.

## Status Snapshot

- `F0`: completed
- `F1`: completed with parser/lexer/round-trip suite green
- `F1b`: completed at implementation level; benchmark discipline is active, but deeper memory/perf characterization remains open
- `F2`: partially completed; canonical simplification is usable, property-based coverage is still missing
- `F4`: in progress with base differentiation, a dedicated `implicit_diff()` API for relations `F(x,y)=0`, integration, definite integrals via TFC, limits (including rational poles, logarithmic lateral singularities, basic growth comparisons, and first signed `-inf` cases), and Taylor series already available with an initial precomputed Maclaurin library plus generalized binomial support for `(1+x)^n`
- `F8`: suspended by project policy until the C++ engine is mature

## Build

Recommended generator:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
./build/tools/cas_ui/cas_ui
```

`Debug` and `RelWithDebInfo` enable AddressSanitizer and UndefinedBehaviorSanitizer on Clang/GCC.

## Optional Manual GUI Lab

The `GUI/` directory is an optional, detachable macOS Qt/QML workspace for manual CAS testing (sidebar, notebook-style cell history, inspector, command palette, and 2D plot preview). It is not part of the mathematical core and is disabled by default.

```bash
cmake -S . -B build-gui -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCAS_ENABLE_GUI=ON
cmake --build build-gui --target cas_gui
```

The Qt/QML frontend requires Qt 6.7+ with `Quick` and `QuickControls2`. Without Qt, the core build remains unaffected.

## Milestone Verification

```bash
./scripts/verify_milestone.sh M0
```

Milestones con gate prestazionale:

```bash
./scripts/verify_milestone.sh M1b
```

## Benchmark Control

Esecuzione benchmark `Release`:

```bash
bash scripts/benchmark.sh
```

Verifica regressioni contro la baseline versionata:

```bash
bash scripts/benchmark.sh --check
```

Aggiornamento intenzionale della baseline dopo una modifica approvata:

```bash
bash scripts/benchmark.sh --update-baseline
```

Output persistente e report JSON:

```bash
bash scripts/benchmark.sh \
  --build-dir build-bench-report \
  --output build-bench-report/current.txt \
  --report-json build-bench-report/report.json \
  --check
```

Policy e baseline correnti:

- `test/benchmarks/policy.json`
- `test/benchmarks/baseline_release.txt`

## Scope Rules

- Pure C++ engine only in this phase
- No Swift or C API implementation work
- No floating-point arithmetic in symbolic foundations
- Tests validate structure and typed behavior, never formatter output
