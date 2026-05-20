# macOS C++ — Build di un agente AI

> Prompt operativo per costruire la GUI macOS in C++/Qt6/QML a partire dal prototipo HTML.

## Ruolo
Sei un senior C++/Qt engineer con esperienza macOS native. Il tuo compito è materializzare in codice Qt6 + QML l'app `CASCalculator` la cui UX è definita dal prototipo HTML in `prototypes/design/CAS Calculator.html` (sezione `CASMacWorkspace`).

## Vincoli tecnici
- **C++20**, **Qt 6.7+**, **QML/Qt Quick 2**
- macOS 12+ universal (arm64 + x86_64)
- CMake 3.27+, vcpkg per dipendenze
- Niente Qt Widgets — solo Qt Quick + QML
- Math rendering: integra **microtex** (LaTeX C++) come libreria statica; fallback al parser custom in `MathRenderer.cpp`
- Embedded kernel: link statico a **libgiac**

## Ordine di esecuzione

### Fase 1 — Skeleton
- [ ] Configura CMake + vcpkg + presets (vedi `02-project-structure.md`)
- [ ] `app/main.cpp` carica `Main.qml` da modulo `CAS`
- [ ] `Info.plist.in` + entitlements; build di un `.app` vuoto firmato

### Fase 2 — Theme & shell
- [ ] `ThemeVM` con tutti i token derivati (esattamente come `palette()` del prototipo)
- [ ] `Main.qml` con TitleBar + Sidebar + Notebook + Inspector + StatusBar
- [ ] Bridge nativo `NativeWindow.mm` (vibrancy, traffic lights, follow system appearance)
- [ ] Tweaks panel rimosso in produzione, sostituito da Settings window standard

### Fase 3 — ViewModels
- [ ] `SessionListVM` con QAbstractListModel (sidebar)
- [ ] `NotebookVM` con `QQmlListProperty<CellVM>` o `QAbstractListModel`
- [ ] `CellVM` con stato (queued/running/ok/error), output, alternatives
- [ ] `KeyboardVM` con palette di funzioni e simboli
- [ ] `CommandPaletteVM` con fuzzy search su lista comandi statica

### Fase 4 — Math rendering
- [ ] `MathRenderer` (QQuickPaintedItem) — scaffold con QPainter
- [ ] Integrazione microtex (PR separato, opzionale)
- [ ] `MathInput` — QML + C++ helper per parsing on-the-fly
- [ ] Test: 50 espressioni canoniche con snapshot pixel-diff

### Fase 5 — Kernel
- [ ] `EmbeddedGiac` — wrapper su `giac::context` + `gen::eval`
- [ ] `RemoteKernel` — `QNetworkAccessManager` + `QWebSocket`
- [ ] `KernelDispatcher` — heuristica embedded vs remote
- [ ] Esecuzione cella in pool `QThreadPool` con `QFuture<EvalResult>`
- [ ] Cancel via `QFutureWatcher`

### Fase 6 — Inspector
- [ ] Tab Grafico: `PlotCanvas` (QQuickPaintedItem) con pan/zoom
- [ ] Tab Passaggi: lista step da `/steps` o local Giac
- [ ] Tab Matrice: `MatrixEditor` con cell-by-cell editing
- [ ] Tab Serie: rendering Taylor expansion

### Fase 7 — Native polish
- [ ] Touch Bar items (Intel Mac legacy)
- [ ] Quick Look provider per file `.cas`
- [ ] Spotlight indexing per sessioni
- [ ] Dock badge con #job in coda
- [ ] Services menu: "Calcola con CAS"

### Fase 8 — Distribution
- [ ] Code signing automatizzato (DEVELOPMENT_TEAM in CMake)
- [ ] Notarization via `xcrun notarytool` in CI
- [ ] Sparkle per aggiornamenti automatici (XPC + EdDSA)
- [ ] Screenshot Mac App Store

## Definition of Done per feature
1. Test unitari (GoogleTest per core, QtTest per VM, Qt Quick Test per QML)
2. Snapshot test di ogni Qt Quick Item visivo
3. Build green su CI macOS-14 (universal)
4. Memory leaks: zero (controllati con `MallocStackLogging` + leaks)
5. Performance: 60 fps scrolling notebook, < 100ms apertura command palette

## Quando NON usare Qt
- Il `.app` bundle, code signing, notarization → script CMake + `xcrun`
- Title bar nativa, vibrancy, traffic lights → NativeWindow.mm
- VoiceOver per matematica → override `accessibilityValue` con MathSpeak
- Touch Bar → `NSTouchBar` direttamente

## Quando NON usare C++ puro
- UI dichiarativa → QML (sempre)
- Animazioni → `Behavior on x { NumberAnimation { ... } }` in QML
- Layout responsivo → `Layouts` QML
- Theming reattivo → property binding a `AppCore.theme.*`

## Checklist prima del primo PR
- [ ] CMake configura su Mac fresh con `cmake --preset macos-arm64`
- [ ] L'app si apre, mostra TitleBar + 3 colonne, kernel embedded risponde a `1+1`
- [ ] ⌘K apre command palette
- [ ] Toggle dark/light segue il sistema
- [ ] `.app` firmato passa `spctl --assess --verbose`
