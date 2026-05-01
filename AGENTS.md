# Repository Guidelines

## Project Structure & Module Organization
`.APROJECT_REFERENCES/` is the source of truth for architecture, module boundaries, and testing rules. Read it before changing code. Core engine code belongs in `src/` and public headers in `include/`; both directories are planned by the project contract even if they are still sparse. Keep manual diagnostics isolated in `tools/cas_ui/`. Do not add Swift, iOS, or bridge code here.

Key references:
- `CLAUDE.md`: project-wide operating rules for agents and contributors
- `CONTRIBUTING.md`: non-negotiable development constraints
- `.APROJECT_REFERENCES/03_ENGINE_MODULES/`: module specs
- `.APROJECT_REFERENCES/04_TESTING_STRATEGY/`: testing philosophy and regression rules

## Build, Test, and Development Commands
Use CMake only.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
./build/tools/cas_ui/cas_ui
ctest --test-dir build --output-on-failure
```

`Debug` and `RelWithDebInfo` enable AddressSanitizer and UBSanitizer on Clang/GCC. The root `CMakeLists.txt` enforces `-Wall -Wextra -Wpedantic -Werror` or `/W4 /WX`; fix warnings instead of weakening flags.

## Coding Style & Naming Conventions
Target `C++20` and keep code/comments in English. Use 4-space indentation and avoid compiler extensions. Follow the module contracts in `.APROJECT_REFERENCES` exactly: no circular dependencies, no implicit fallback from symbolic/calculus code into numeric code, and no floating point in the symbolic core. Prefer explicit result/error types in public APIs as required by `CLAUDE.md`.

Name files and types by domain responsibility, for example `lexer.cpp`, `parser.hpp`, `AlgebraEngine`, `NumericEvaluator`. Keep tool-only code under `tools/`, never mixed into engine modules.

## Testing Guidelines
GoogleTest is the intended framework; the `test/` subtree is reserved for it. Tests must validate AST structure or mathematical equivalence, never `toString()` output. Add regression coverage for every bug fix and run sanitizer-enabled builds locally before closing a task.

## Commit & Pull Request Guidelines
Git history is not available in this checkout, so use short imperative commit messages such as `Add parser error propagation`. Keep each commit focused on one module or behavior. PRs should include: purpose, affected module, references consulted in `.APROJECT_REFERENCES`, test evidence (`ctest`, sanitizer runs, manual `cas_ui` checks), and screenshots only if `tools/cas_ui` UI behavior changes.
## Mandato di Guardia Architetturale
Gli agenti sono i custodi dell'integrità del motore. È obbligatorio:
1. **Rilevamento Preventivo Monoliti**: Se un file supera le 400 righe, l'agente deve proporre uno split preventivo in moduli verticali (es. `simplify_arithmetic.cpp`).
2. **Hardening Error Path**: Ogni funzione deve tornare `Result<T>`. Trovati `throw` o `assert` in path simbolici vanno segnalati come debiti critici.
3. **Validazione AST**: Ogni nuovo nodo deve usare `AstArena::make` e rispettare l'Hash-consing (identità dei puntatori).
4. **Coerenza Matematica**: Bloccare l'uso di "cerotti" numerici o euristiche fragili nel core simbolico. Preferire sempre il sistema di rewrite universale.

## Agent-Specific Instructions
...
If a required spec file is missing or marked as future work, stop and ask for the specification instead of inventing architecture. This repository is documentation-driven; implementation follows approved references, not ad hoc design.
