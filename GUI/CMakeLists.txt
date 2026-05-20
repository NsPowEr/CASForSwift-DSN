cmake_minimum_required(VERSION 3.25)

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
