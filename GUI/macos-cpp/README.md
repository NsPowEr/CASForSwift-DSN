# macOS GUI in C++ — Pacchetto completo

Stack scelto: **Qt 6.7 + QML** con bridge Objective-C++ per le feature Apple-specifiche.
Motivazioni in [`01-architecture.md`](./01-architecture.md).

## Indice
1. [`01-architecture.md`](./01-architecture.md) — confronto Qt vs ImGui vs Slint vs AppKit, scelta finale
2. [`02-project-structure.md`](./02-project-structure.md) — repo, CMake, vcpkg
3. [`03-native-integration.md`](./03-native-integration.md) — bridge AppKit (vibrancy, Touch Bar, Quick Look)
4. [`AGENTS.md`](./AGENTS.md) — task list per agente AI

## Codice scaffold incluso
| File | Cosa | Stato |
|---|---|---|
| `src/viewmodels/AppCore.h` | Singleton root, esposto a QML | scaffold |
| `src/viewmodels/ThemeVM.h/.cpp` | Token engine (porting da `palette()`) | implementato |
| `src/qml-items/MathRenderer.h/.cpp` | Math renderer custom QQuickPaintedItem | scaffold |
| `src/platform/macos/NativeWindow.mm` | Vibrancy + traffic lights + appearance | implementato |
| `app/main.cpp` | Bootstrap | implementato |
| `qml/Main.qml` | Window root | scaffold |
| `qml/workspace/CellView.qml` | Cella In[n]/Out[n] | scaffold |

## Sintesi della scelta

| Esigenza CAS | Soluzione Qt | Soluzione native |
|---|---|---|
| Math rendering 2D | `QQuickPaintedItem` + microtex | NSBezierPath custom |
| Plot 2D/3D | `QtCharts` + `QQuick3D` | Charts.framework limitata |
| Editor 2D formule | `MathInput.qml` + AST C++ | NSTextView custom |
| Code/REPL | `Qt Quick` + KSyntaxHighlighting | NSTextView |
| Command palette | QML overlay + fuzzy C++ | NSPanel |
| Traffic lights / vibrancy | bridge `NativeWindow.mm` | nativo |
| Sync con backend | `QNetworkAccessManager` + `QWebSocket` | URLSession |
| Kernel embedded | link diretto libgiac (C++) | stesso |
| Cross-platform futuro | gratis | da riscrivere |

## Punto chiave
**Giac è C++ → linkalo direttamente.** Niente FFI, niente overhead. Stesso codice del kernel server-side. Per sandbox: XPC service separato (equivalente macOS-nativo del subprocess + seccomp che ho usato sul backend).
