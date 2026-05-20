# macOS C++ — Struttura del progetto

```
macos-cpp/
├── CMakeLists.txt
├── CMakePresets.json
├── vcpkg.json
├── README.md
│
├── app/                          # Eseguibile .app bundle
│   ├── CMakeLists.txt
│   ├── main.cpp
│   ├── Info.plist.in
│   ├── CASCalculator.entitlements
│   └── resources/
│       ├── Assets.xcassets/      # icone .icns
│       └── qml.qrc               # bundle QML
│
├── src/
│   ├── core/                     # C++ puro, niente Qt
│   │   ├── models/
│   │   │   ├── Cell.h/.cpp
│   │   │   ├── Output.h/.cpp
│   │   │   ├── Session.h/.cpp
│   │   │   ├── Variable.h/.cpp
│   │   │   └── Plot.h/.cpp
│   │   ├── parser/
│   │   │   ├── MathLexer.h/.cpp
│   │   │   ├── MathParser.h/.cpp
│   │   │   └── AstNode.h/.cpp
│   │   ├── kernel/
│   │   │   ├── IKernel.h          # interfaccia astratta
│   │   │   ├── EmbeddedGiac.h/.cpp # link diretto a libgiac
│   │   │   ├── RemoteKernel.h/.cpp # client REST/WS verso backend
│   │   │   └── KernelDispatcher.h/.cpp
│   │   └── util/
│   │       ├── Result.h           # tagged union ok/err
│   │       └── Logger.h
│   │
│   ├── viewmodels/                # QObject — bridge a QML
│   │   ├── AppCore.h/.cpp         # root singleton
│   │   ├── ThemeVM.h/.cpp
│   │   ├── SessionListVM.h/.cpp
│   │   ├── NotebookVM.h/.cpp
│   │   ├── CellVM.h/.cpp
│   │   ├── KeyboardVM.h/.cpp
│   │   ├── CommandPaletteVM.h/.cpp
│   │   └── PlotVM.h/.cpp
│   │
│   ├── qml-items/                 # custom Qt Quick items in C++
│   │   ├── MathRenderer.h/.cpp    # QQuickPaintedItem — disegna LaTeX
│   │   ├── PlotCanvas.h/.cpp      # 2D plotting interattivo
│   │   ├── MatrixEditor.h/.cpp
│   │   └── KaTeXItem.h/.cpp       # alternativa: WebEngine + KaTeX
│   │
│   └── platform/macos/            # Objective-C++ (.mm)
│       ├── NativeWindow.mm        # vibrancy, traffic light positioning
│       ├── TouchBar.mm
│       ├── ServicesProvider.mm
│       └── DockBadge.mm
│
├── qml/
│   ├── Main.qml                   # ApplicationWindow root
│   ├── Theme.qml                  # singleton tokens
│   ├── pages/
│   │   ├── NotebookPage.qml
│   │   ├── ConverterPage.qml
│   │   └── FormularyPage.qml
│   ├── workspace/
│   │   ├── Sidebar.qml
│   │   ├── Notebook.qml
│   │   ├── CellView.qml
│   │   ├── Inspector.qml
│   │   ├── StatusBar.qml
│   │   └── TitleBar.qml
│   ├── controls/
│   │   ├── TabPill.qml
│   │   ├── KeyButton.qml
│   │   ├── CommandPalette.qml
│   │   ├── TweakPanel.qml
│   │   └── IconButton.qml
│   └── math/
│       ├── MathView.qml           # wrapper di MathRenderer
│       ├── MathInput.qml          # editor 2D
│       └── MatrixView.qml
│
├── tests/
│   ├── core/                      # GoogleTest
│   ├── parser/
│   └── viewmodels/                # QtTest + Qt Quick Test
│
└── third_party/                   # vendored se necessario
    ├── giac/                      # CAS engine (link statico)
    └── microtex/                  # rendering LaTeX in C++ (opzione)
```

## CMake — struttura essenziale

`CMakeLists.txt` (root):

```cmake
cmake_minimum_required(VERSION 3.27)
project(CASCalculator VERSION 0.1.0 LANGUAGES CXX OBJCXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET 12.0)
set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64")

find_package(Qt6 6.7 REQUIRED COMPONENTS
    Core Gui Qml Quick QuickControls2
    Network WebSockets Charts Concurrent Test)

qt_standard_project_setup(REQUIRES 6.7)

add_subdirectory(src/core)
add_subdirectory(src/viewmodels)
add_subdirectory(src/qml-items)
add_subdirectory(src/platform/macos)
add_subdirectory(app)

enable_testing()
add_subdirectory(tests)
```

`app/CMakeLists.txt`:

```cmake
qt_add_executable(CASCalculator MACOSX_BUNDLE
    main.cpp
)

qt_add_qml_module(CASCalculator
    URI CAS
    VERSION 1.0
    QML_FILES
        ../qml/Main.qml
        ../qml/Theme.qml
        ../qml/workspace/Sidebar.qml
        # ... tutti gli altri
    RESOURCES
        resources/Assets.xcassets
)

target_link_libraries(CASCalculator PRIVATE
    cas_core cas_viewmodels cas_qml_items cas_platform_macos
    Qt6::Quick Qt6::QuickControls2 Qt6::Charts Qt6::WebSockets
)

set_target_properties(CASCalculator PROPERTIES
    MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist.in
    MACOSX_BUNDLE_GUI_IDENTIFIER app.cascalc.desktop
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    XCODE_ATTRIBUTE_CODE_SIGN_ENTITLEMENTS
        ${CMAKE_CURRENT_SOURCE_DIR}/CASCalculator.entitlements
)
```

`vcpkg.json`:

```json
{
  "name": "cas-calculator",
  "version": "0.1.0",
  "dependencies": [
    "fmt", "spdlog", "nlohmann-json",
    {"name": "qt", "version>=": "6.7.0",
     "features": ["quick", "quickcontrols2", "websockets", "charts"]}
  ]
}
```

## Build

```bash
cmake --preset macos-arm64
cmake --build --preset macos-arm64 --config Release
# Universal binary
cmake --preset macos-universal
```

Output: `build/macos-arm64/app/CASCalculator.app` — pronto per code-sign e notarization.
