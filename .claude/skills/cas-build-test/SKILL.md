---
name: cas-build-test
description: Build CAS con Ninja e lancia cas_foundation_tests con gtest_filter opzionale. Usa ASAN_OPTIONS=detect_leaks=0. Argomenti opzionali: build dir (default: build) e filtro GTest.
---

Build the CAS engine and run targeted GoogleTest suite.

## Usage

`/cas-build-test [build_dir] [gtest_filter]`

- `build_dir`: default `build`. Alternatives: `build-bench` (release/benchmark), `build-golden` (golden runner)
- `gtest_filter`: default `AcidTest.*:SupremeTest.*`. Use GTest filter syntax.

## Steps

1. Parse arguments from user message (build_dir and filter if provided)
2. Run build:
   ```bash
   ninja -C {build_dir}
   ```
3. If build fails → report error, stop.
4. Run tests:
   ```bash
   ASAN_OPTIONS=detect_leaks=0 gtimeout 120 ./{build_dir}/cas_foundation_tests --gtest_filter="{filter}"
   ```
5. Report: pass count, fail count, failed test names, total time.

## Examples

- `/cas-build-test` → builds `build/`, runs AcidTest + SupremeTest
- `/cas-build-test build-check AcidTest.Test6*` → specific build dir + filter
- `/cas-build-test build CalculusDiffTest*:CalculusIntegrateTest*` → calculus suite
