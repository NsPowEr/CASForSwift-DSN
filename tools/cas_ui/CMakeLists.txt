add_executable(cas_ui src/main.cpp)

# Link contro il core engine C++ locale
target_link_libraries(cas_ui PRIVATE cas_core)

# Niente bridge Swift o riferimenti incrociati qui
set_target_properties(cas_ui PROPERTIES FOLDER "Tools")
