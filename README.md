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

## Coverage (F0.3)

Line coverage via gcov/lcov. Enable with:

```bash
cmake -S . -B build_cov -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_COVERAGE=ON \
  -DCAS_ENABLE_SANITIZERS=OFF
cmake --build build_cov
cmake --build build_cov --target coverage
# Report: build_cov/coverage/index.html
```

Requires `lcov` and `genhtml` (`brew install lcov` on macOS).
Coverage target: ≥90% line coverage per area marked `Risolta` (F0 exit gate).

## Property-Based Tests (F0.4)

Property tests use [rapidcheck](https://github.com/emil-e/rapidcheck) (fetched automatically via CMake FetchContent).

Run all property tests:

```bash
ctest --test-dir build -R property --output-on-failure
```

Six test areas: GCD/LCM identity, D(∫f)=f corpus, roots substitution, matrix inverse, factor reconstruction, partial fractions roundtrip.

## Sanitizers (F0.6)

ASan + UBSan are enabled in `Debug` and `RelWithDebInfo` builds by default (`CAS_ENABLE_SANITIZERS=ON`).

Dedicated sanitizer build:

```bash
cmake -S . -B build_san -G Ninja \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCAS_ENABLE_SANITIZERS=ON
cmake --build build_san
ASAN_OPTIONS=detect_leaks=0 ctest --test-dir build_san --output-on-failure
```

TSan (thread sanitizer) availability is detected at configure time and logged.

## Build Cache (F0.6)

ccache is the recommended compiler launcher to speed up incremental rebuilds:

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
```

On macOS: `brew install ccache`. Typically provides 5–20× speedup on rebuilds during active development.

## Anti-Monolith Gate (F0.6)

All `.cpp` and `.hpp` files under `src/` and `include/` must be ≤ 500 lines.
Pre-existing violations are whitelisted in `scripts/file_size_whitelist.txt` with mandatory split tickets.

Check locally:

```bash
bash scripts/check_file_size.sh --verbose
```

This check also runs as a CI job on every push and PR (`anti-monolith` job).

## Benchmark Regression Gate (F0.6)

After building the benchmark, compare against the pinned baseline:

```bash
bash scripts/check_bench_regression.sh \
  --current build-bench/current.txt \
  --threshold 10
```

Exit 1 if any metric regresses by more than 10%. The baseline is `test/benchmarks/baseline_release.txt`.

## Mutation Testing (F0.6, optional)

Mutation testing via [mull](https://github.com/mull-project/mull) (optional, runs weekly in CI):

```bash
brew install mull          # macOS
bash scripts/run_mutation.sh
# Report: mutation-report/index.html
```

Target mutation score: ≥70% per module. Non-blocking on PR; weekly scheduled CI job only.

## Scope Rules

- Pure C++ engine only in this phase
- No Swift or C API implementation work
- No floating-point arithmetic in symbolic foundations
- Tests validate structure and typed behavior, never formatter output
