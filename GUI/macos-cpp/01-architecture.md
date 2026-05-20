# macOS GUI in C++ — Architettura

> Scelta tecnica: **Qt 6.7 + QML (Qt Quick)** con bridge nativo macOS via Objective-C++ per le funzioni Apple-specifiche (NSWindow vibrancy, Touch Bar, Apple Pencil su iPad sidecar, Services menu).

## Perché Qt6 + QML (ragionamento)

### Cosa serve a un CAS desktop pro
1. **Rendering matematico 2D dinamico** — frazioni, integrali, matrici, sommatorie
2. **Editor di codice/testo** con syntax highlighting, multi-cursor, undo
3. **Plotting 2D/3D** interattivo con pan/zoom/rotate
4. **Command palette** stile Raycast/VSCode
5. **Layout responsivo** (notebook + inspector + sidebar collassabili)
6. **Theming dark/light** istantaneo, accent color customizzabile
7. **Native feel** macOS (vibrancy, toolbar unificata, traffic lights)

### Confronto serio dei 4 candidati

| Criterio | **Qt 6 + QML** | Dear ImGui | Slint | AppKit C++ |
|---|---|---|---|---|
| Math rendering custom | ⭐⭐⭐⭐⭐ `QQuickPaintedItem` + `QPainter` | ⭐⭐⭐ canvas immediato | ⭐⭐⭐⭐ `Path` API | ⭐⭐⭐⭐ NSBezierPath |
| Code editor pronto | ⭐⭐⭐⭐ KSyntaxHighlighting, QSci | ⭐⭐ TextEditor base | ⭐⭐ da implementare | ⭐⭐⭐⭐ NSTextView |
| Plot 2D/3D | ⭐⭐⭐⭐⭐ QtCharts + Qt3D + QQuick3D | ⭐⭐⭐ ImPlot | ⭐⭐ da fare | ⭐⭐ Charts.framework limitata |
| Native macOS feel | ⭐⭐⭐⭐ con un po' di QSS + native bridge | ⭐ no | ⭐⭐⭐ default | ⭐⭐⭐⭐⭐ |
| Productivity (DSL UI) | ⭐⭐⭐⭐⭐ QML | ⭐⭐ tutto in C++ | ⭐⭐⭐⭐⭐ DSL dedicato | ⭐⭐ XIB/storyboard o codice |
| Hot reload UI | ⭐⭐⭐⭐⭐ `qml` runtime | ⭐⭐ rebuild | ⭐⭐⭐⭐ | ⭐⭐ |
| Cross-platform futuro | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐ solo Apple |
| Maturità ecosistema | ⭐⭐⭐⭐⭐ 30 anni | ⭐⭐⭐⭐ | ⭐⭐ giovane | ⭐⭐⭐⭐⭐ |
| Licenza | LGPLv3 / commercial | MIT | GPL/commercial/Royalty | Apple SDK |
| Accessibility | ⭐⭐⭐⭐ QAccessible + VoiceOver via Cocoa | ⭐ scarsa | ⭐⭐ | ⭐⭐⭐⭐⭐ nativa |

### Verdetto
**Qt 6 + QML** vince perché:
- QML è dichiarativo (vicino a SwiftUI come modello mentale, ma in C++/JS-like)
- Il C++ resta puro nei reducer/business logic; QML gestisce solo presentazione
- Il bridge `QObject` → QML è dichiarativo via `Q_PROPERTY`
- Posso fare il math renderer custom in C++ con `QPainter` (vector, super performante) e usarlo da QML
- Conservo l'opzione cross-platform "gratis"

Le poche cose non native (vibrancy NSWindow, Touch Bar, Services menu, Quick Look, dock badge) le risolvo con un thin layer **Objective-C++** in `platform/macos/`. Qt6 espone già `nativeInterface()` per accesso a `NSView*`/`NSWindow*`.

---

## Stack scelto

```
┌────────────────────────────────────────────────┐
│  QML (Qt Quick 2)         → presentation        │
│  ─ Window, Pages, Sidebar, Notebook, Inspector  │
│  ─ MathRenderer item (C++ backed)               │
│  ─ Theme singleton (token-driven)               │
└──────────────────┬─────────────────────────────┘
                   │ Q_PROPERTY / Signal-Slot
┌──────────────────▼─────────────────────────────┐
│  C++ ViewModels (QObject)                       │
│  ─ SessionVM, NotebookVM, KernelVM, ThemeVM     │
└──────────────────┬─────────────────────────────┘
                   │
┌──────────────────▼─────────────────────────────┐
│  C++ Core (no Qt)         → portable            │
│  ─ Models (Cell, Output, Variable, Session)     │
│  ─ KernelClient (libcurl + websocketpp)         │
│  ─ MathParser (own AST)                         │
│  ─ EmbeddedKernel (libgiac C++ direttamente!)   │
└──────────────────┬─────────────────────────────┘
                   │ Obj-C++ bridge dove serve
┌──────────────────▼─────────────────────────────┐
│  Platform layer (macOS)                         │
│  ─ NSWindow vibrancy + traffic lights           │
│  ─ Touch Bar items                              │
│  ─ Services / Quick Look                        │
│  ─ Notarization, code signing                   │
└─────────────────────────────────────────────────┘
```

## Vantaggio chiave per il CAS
**Giac è già C++**: lo linko direttamente nel processo (no FFI, no overhead). Il kernel embedded è la stessa libreria che gira sul backend. Per la sandbox uso `XPC service` separato (ogni evaluate gira in un processo XPC isolato — equivalente macOS-nativo del subprocess + seccomp che ho descritto per il backend).

## Build system
- **CMake 3.27+** con preset per macOS arm64 + universal binary
- Qt 6.7+ via `find_package(Qt6 COMPONENTS Quick QuickControls2 Charts ...)`
- `vcpkg` o `Conan` per dipendenze C++ third-party (libgiac, libcurl, websocketpp, fmt, spdlog)
- Output: `.app` bundle code-signed + notarized

## Stati / dependency injection
Userò **single-source-of-truth pattern**:
- `AppCore` singleton C++ in `main.cpp` espone i ViewModel come `qmlRegisterSingletonInstance`
- ViewModels sono `QObject` con `Q_PROPERTY NOTIFY` — QML reagisce automaticamente
- Niente Redux/MVU: per Qt è overkill. Property binding + signal/slot bastano.

## Threading
- UI thread: solo Qt Quick
- Kernel thread: pool di `QThread` con `QFuture`/`QPromise` per job
- Network thread: `QNetworkAccessManager` async (default)
- WebSocket: `QWebSocket` su thread dedicato

## Compatibilità target
- **macOS 12 Monterey** minimo (universal binary x86_64 + arm64)
- Test su Sonoma 14, Sequoia 15, Tahoe 16
- Qt 6.7+ richiesto per `Qt Quick Controls 2` macOS style
