cmake_minimum_required(VERSION 3.25)

add_library(cas_gui_session STATIC
    src/CasGuiSession.cpp
)

target_include_directories(cas_gui_session
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(cas_gui_session PUBLIC cas_core)
target_compile_features(cas_gui_session PUBLIC cxx_std_20)

find_package(Qt6 6.7 QUIET COMPONENTS Quick QuickControls2)

if(NOT Qt6_FOUND)
    message(STATUS "Qt6 Quick/QuickControls2 not found; optional target cas_gui is not built.")
    return()
endif()

set(CMAKE_AUTOMOC ON)

qt_standard_project_setup()

qt_add_executable(cas_gui
    src/QtGuiMain.cpp
)

qt_add_qml_module(cas_gui
    URI CASGui
    VERSION 1.0
    QML_FILES
        qml/Main.qml
)

target_link_libraries(cas_gui
    PRIVATE
        cas_gui_session
        Qt6::Quick
        Qt6::QuickControls2
)

target_compile_features(cas_gui PRIVATE cxx_std_20)
