# test/golden/CMakeLists.txt — Golden test runner (FASE 0.5)
#
# Integration instructions for Block C (top-level CMakeLists.txt):
#   1. Add at the end of the top-level CMakeLists.txt:
#        add_subdirectory(test/golden)
#   2. Optionally register the runner as a CTest label:
#        enable_testing()
#        add_test(NAME golden_runner_simplify
#                 COMMAND golden_runner
#                   ${CMAKE_SOURCE_DIR}/test/golden/corpus/simplify/basic.jsonl
#                   ${CMAKE_BINARY_DIR}/golden_maxima_out/simplify
#                   --json ${CMAKE_BINARY_DIR}/golden_simplify.json)
#        set_tests_properties(golden_runner_simplify PROPERTIES LABELS "golden")
#
# The runner is standalone (no GTest); link against cas_core only.

cmake_minimum_required(VERSION 3.20)

add_executable(golden_runner
    main.cpp
)

target_include_directories(golden_runner PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(golden_runner PRIVATE cas_core)

target_compile_features(golden_runner PRIVATE cxx_std_20)

# Optional: register CTest targets if testing is enabled
if(BUILD_TESTING)
    set(GOLDEN_CORPUS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/corpus")
    set(GOLDEN_MAXIMA_OUT "${CMAKE_BINARY_DIR}/golden_maxima_out")

    foreach(AREA simplify factor gcd integrate diff limit solve series special_fn)
        set(CORPUS_FILE "${GOLDEN_CORPUS_DIR}/${AREA}/basic.jsonl")
        set(MAXIMA_OUT  "${GOLDEN_MAXIMA_OUT}/${AREA}")
        set(JSON_OUT    "${CMAKE_BINARY_DIR}/golden_${AREA}.json")

        if(EXISTS "${CORPUS_FILE}")
            add_test(NAME "golden_${AREA}"
                COMMAND golden_runner
                    "${CORPUS_FILE}"
                    "${MAXIMA_OUT}"
                    --json "${JSON_OUT}"
            )
            set_tests_properties("golden_${AREA}" PROPERTIES
                LABELS "golden"
                TIMEOUT 120
            )
        endif()
    endforeach()

    # Bronstein sub-corpus
    if(EXISTS "${GOLDEN_CORPUS_DIR}/bronstein/integrals.jsonl")
        add_test(NAME "golden_bronstein"
            COMMAND golden_runner
                "${GOLDEN_CORPUS_DIR}/bronstein/integrals.jsonl"
                "${GOLDEN_MAXIMA_OUT}/bronstein"
                --json "${CMAKE_BINARY_DIR}/golden_bronstein.json"
        )
        set_tests_properties("golden_bronstein" PROPERTIES
            LABELS "golden"
            TIMEOUT 300
        )
    endif()
endif()
