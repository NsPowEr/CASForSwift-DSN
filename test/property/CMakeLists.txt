# F0.4 — Property-based tests via rapidcheck + GTest
# Each test file is small and focused on one algebraic area.
# Run with: ctest -R property

set(CAS_PROPERTY_SOURCES
    test_gcd_lcm_property.cpp
    test_diff_integrate_inverse.cpp
    test_roots_substitute.cpp
    test_matrix_inverse.cpp
    test_factor_product.cpp
    test_partial_fractions_roundtrip.cpp
)

add_executable(cas_property_tests ${CAS_PROPERTY_SOURCES})

target_link_libraries(cas_property_tests
    PRIVATE
        cas_core
        rapidcheck
        rapidcheck_gtest
        GTest::gtest_main
)

target_include_directories(cas_property_tests
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/test/helpers
)

target_compile_features(cas_property_tests PRIVATE cxx_std_20)

include(GoogleTest)
gtest_discover_tests(cas_property_tests
    DISCOVERY_TIMEOUT 120
    PROPERTIES
        ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0"
        LABELS "property"
)
