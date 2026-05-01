# GUI Architecture Reduction

Questo documento riduce `.APROJECT_REFERENCES/08_handoff-GUI` al caso reale di questo repository.

## Decisione

La cartella `GUI/` e' un modulo opzionale di test manuale. Non implementa la CAS Calculator completa descritta dal handoff. Implementa solo un laboratorio locale collegato al core C++.

## Cosa viene preso dal handoff

- macOS in C++/Qt/QML, non SwiftUI.
- Notebook-like workspace: input, output, inspector, plot pane.
- Math output multiplo: plain text, LaTeX, rendering 2D testuale.
- Plot 2D/3D come feature visuale, ma con priorita al 2D gia supportabile dal sampler C++.
- Shortcuts e command palette come obiettivo futuro.

## Cosa viene escluso qui

- Backend FastAPI, REST, WebSocket, Redis, PostgreSQL, MinIO.
- Auth, Apple Sign-In, Cloud sync, multi-device state.
- Claude/AI client lato UI.
- Swift bridge, XCFramework, SwiftData, SwiftUI.
- Packaging notarizzato e auto-update.

## Moduli

```
GUI/
├── CMakeLists.txt              # entrypoint opzionale
├── src/
│   ├── CasGuiSession.hpp       # adapter C++ puro verso il CAS
│   ├── CasGuiSession.cpp
│   └── QtGuiMain.cpp           # frontend Qt/QML, solo se CAS_ENABLE_GUI=ON
└── qml/
    └── Main.qml                # UI test lab
```

## Dependency rule

```
QML -> QtGuiMain.cpp -> CasGuiSession -> cas_core
```

Il verso inverso e' vietato. `cas_core` non conosce Qt, QML, AppKit o dettagli GUI.

## Attacca e stacca

Il root `CMakeLists.txt` contiene solo:

```cmake
option(CAS_ENABLE_GUI "Build the optional detachable manual GUI lab." OFF)
if(CAS_ENABLE_GUI)
    add_subdirectory(GUI)
endif()
```

Quindi:

- default: nessun target GUI, nessuna dipendenza Qt;
- abilitato: si costruisce `cas_gui`;
- rimozione: eliminare `GUI/` e togliere il blocco opzionale non altera il core.

## Regole matematiche

- La GUI non semplifica da sola: chiama sempre il motore.
- La GUI non introduce euristiche numeriche nel symbolic core.
- Il plot usa solo percorsi esplicitamente numerici (`cas::numeric`).
- Ogni nuova operazione manuale deve essere prima un adapter sottile, poi eventualmente un test.
