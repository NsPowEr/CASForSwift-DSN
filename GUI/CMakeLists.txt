cmake_minimum_required(VERSION 3.25)

# ─── Platform-native layer (macOS only) ──────────────────────────────────────
# Compiled separately so that Objective-C++ does not contaminate the main
# C++20 translation units.  Wrapped in if(APPLE) so the library is a no-op
# on hypothetical future non-Apple platforms.
if(APPLE)
    add_library(cas_platform_macos STATIC
        macos-cpp/src/platform/macos/NativeWindow.mm
        macos-cpp/src/platform/macos/NativeWindow.h
    )
    set_source_files_properties(
        macos-cpp/src/platform/macos/NativeWindow.mm
        PROPERTIES COMPILE_FLAGS "-fobjc-arc -fobjc-weak"
    )
    target_include_directories(cas_platform_macos
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/macos-cpp/src/platform/macos
    )
    target_link_libraries(cas_platform_macos
        PUBLIC
            Qt6::Gui
            "-framework AppKit"
            "-framework Cocoa"
            "-framework CoreServices"
            "-framework QuartzCore"
    )
    target_compile_features(cas_platform_macos PUBLIC cxx_std_20)
endif()

# Libreria riusabile del layer GUI locale (AppCore/ViewModel/QML items)
qt_add_library(cas_gui_vm STATIC
    macos-cpp/src/core/kernel/KernelDispatcher.h
    macos-cpp/src/core/storage/SessionStore.cpp
    macos-cpp/src/core/storage/SessionStore.h
    macos-cpp/src/viewmodels/AppCore.cpp
    macos-cpp/src/viewmodels/AppCoreSessionActions.cpp
    macos-cpp/src/viewmodels/AppCoreCellActions.cpp
    macos-cpp/src/viewmodels/AppCorePresentation.cpp
    macos-cpp/src/viewmodels/AppCoreUndo.cpp
    macos-cpp/src/viewmodels/AppCoreWorkspace.cpp
    macos-cpp/src/viewmodels/AppCoreLayoutState.cpp
    macos-cpp/src/viewmodels/AppCore.h
    macos-cpp/src/viewmodels/CellVM.h
    macos-cpp/src/viewmodels/CommandPaletteVM.h
    macos-cpp/src/viewmodels/KeyboardVM.h
    macos-cpp/src/viewmodels/NotebookVM.h
    macos-cpp/src/viewmodels/PlotVM.h
    macos-cpp/src/viewmodels/SessionListVM.h
    macos-cpp/src/viewmodels/ThemeVM.cpp
    macos-cpp/src/viewmodels/MathHighlighter.cpp
    macos-cpp/src/viewmodels/MathHighlighter.h
    macos-cpp/src/viewmodels/AutocompleteVM.cpp
    macos-cpp/src/viewmodels/AutocompleteVM.h
    macos-cpp/src/qml-items/MathRenderer.cpp
    macos-cpp/src/qml-items/MathRenderer.h
    macos-cpp/src/qml-items/PlotCanvas.cpp
    macos-cpp/src/qml-items/PlotCanvas.h
)

target_include_directories(cas_gui_vm
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/macos-cpp/src
        ${CMAKE_CURRENT_SOURCE_DIR}/macos-cpp/src/viewmodels
        ${CMAKE_CURRENT_SOURCE_DIR}/macos-cpp/src/qml-items
        ${CMAKE_CURRENT_SOURCE_DIR}/macos-cpp/src/platform/macos
)

# Modulo QML
qt_add_qml_module(cas_gui_vm
    URI CAS
    VERSION 1.0
    OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CAS
    QML_FILES
        macos-cpp/qml/Main.qml
        macos-cpp/qml/workspace/Sidebar.qml
        macos-cpp/qml/workspace/Notebook.qml
        macos-cpp/qml/workspace/Inspector.qml
        macos-cpp/qml/workspace/InspectorResultTab.qml
        macos-cpp/qml/workspace/InspectorRepresentationsTab.qml
        macos-cpp/qml/workspace/InspectorStepsTab.qml
        macos-cpp/qml/workspace/InspectorVariablesTab.qml
        macos-cpp/qml/workspace/InspectorPlotTab.qml
        macos-cpp/qml/workspace/StatusBar.qml
        macos-cpp/qml/workspace/TitleBar.qml
        macos-cpp/qml/workspace/CellView.qml
        macos-cpp/qml/controls/CommandPalette.qml
        macos-cpp/qml/math/MathView.qml
        macos-cpp/qml/math/MathInput.qml
        macos-cpp/qml/math/AutocompletePopup.qml
)

target_link_libraries(cas_gui_vm
    PUBLIC
        cas_gui_bridge
        Qt6::Core
        Qt6::Gui
        Qt6::Qml
        Qt6::Quick
        Qt6::QuickControls2
        Qt6::Concurrent
)

# Link the platform library only when it was actually built (macOS).
if(APPLE)
    target_link_libraries(cas_gui_vm PUBLIC cas_platform_macos)
endif()

target_compile_features(cas_gui_vm PUBLIC cxx_std_20)

# Eseguibile principale
qt_add_executable(cas_gui
    macos-cpp/app/main.cpp
)

target_link_libraries(cas_gui
    PRIVATE
        cas_gui_vm
)

target_compile_features(cas_gui PRIVATE cxx_std_20)
