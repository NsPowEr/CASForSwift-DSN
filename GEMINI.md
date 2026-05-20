# REAL CAS ENGINE C++ — Project Instructions

This document provides essential context and instructions for developing the REAL CAS ENGINE C++. All agents and contributors must strictly adhere to the architecture and workflows defined here.

## Project Overview
An industrial-grade Computer Algebra System (CAS) implemented in modern C++20. The project prioritizes mathematical correctness, performance via structural sharing, and a zero-hardcode policy.

- **Current Focus:** Calculus foundations (Module `F4`), including differentiation, integration (Risch/Hermite), limits (Gruntz/MRV), and Taylor series.
- **Key Modules:**
  - `foundation`: BigInt (limb-based) and Rational arithmetic.
  - `ast`: Immutable expression tree with interning and Arena allocation.
  - `symbolic`: Simplification core, assumptions, and rewrite engine.
  - `algebra`: Polynomial arithmetic, GCD (modular/heuristic), and factorization.
  - `calculus`: Symbolic differentiation and integration.

## Architectural Mandates (The "Technical Constitution")
Refer to `CLAUDE.md` for the full "Technical Constitution". Key highlights:

1.  **NO HARDCODE:** Every constant must have a formal mathematical justification. Computational limits must be configurable via `CASContext`.
2.  **Structural Sharing:** Expressions are immutable. Functions must return the original `ExprPtr` if no changes are made. Pointer identity is used for $O(1)$ equality checks.
3.  **Arena Allocation:** All AST nodes must be allocated in an `AstArena`. `std::make_unique` or `new` for individual nodes is strictly forbidden.
4.  **Exact Arithmetic:** Use `BigInt` and `Rational`. `int64_t` or `double` are forbidden in the symbolic core.
5.  **Zero Warning Policy:** The project must compile without warnings (`-Wall -Wextra -Wpedantic`).

## Building and Running

### Build Commands
Recommended generator is Ninja with Clang/GCC for sanitizer support.

```bash
# Configure and build the core engine
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build

# Build with GUI enabled (optional lab)
cmake -S . -B build-gui -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCAS_ENABLE_GUI=ON
cmake --build build-gui --target cas_gui
```

### Testing
Testing is mandatory for every change.

```bash
# Run all tests
ctest --test-dir build --output-on-failure

# Run specific milestone verification
./scripts/verify_milestone.sh M1b
```

### Benchmarking
Performance gates are enforced via benchmarking.

```bash
# Run release benchmarks
bash scripts/benchmark.sh

# Check against baseline
bash scripts/benchmark.sh --check
```

## Development Conventions

- **Language:** C++20.
- **Error Handling:** No `throw/catch`. Use `Result<T>` (monadic error handling).
- **File Length:** Maximum 500 lines per file. Split into specialized modules if exceeded.
- **Naming:** Follow existing conventions (snake_case for functions/variables, PascalCase for classes).
- **Documentation:** Use Doxygen-style comments for public APIs in `include/cas/`.

## Key Files
- `CLAUDE.md`: The "Supreme Law" - architectural and anti-hardcode rules.
- `CAS_TASKS.md`: Current development tasks and roadmap.
- `include/cas/ast.hpp`: Core expression definitions and `AstArena`.
- `include/cas/symbolic.hpp`: `CASContext` and simplification interfaces.
- `CMakeLists.txt`: Build system configuration.

## Workflow: Before Submitting
1.  **Mathematical Integrity:** Ensure 100% test pass rate.
2.  **Benchmark Gate:** No performance regressions against `baseline_release.txt`.
3.  **Sanitizers:** ASan and UBSan must be clean.
4.  **Zero Warnings:** Compilation must be warning-free.
