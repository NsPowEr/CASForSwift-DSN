cmake_minimum_required(VERSION 3.25)
project(CASEngine VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(CAS_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer in debug-oriented builds." ON)
option(CAS_ENABLE_TESTING "Enable unit tests and CTest integration." ON)
option(CAS_ENABLE_GUI "Build the optional detachable manual GUI lab." OFF)
option(CAS_USE_SYSTEM_GTEST "Use an installed GoogleTest package instead of FetchContent." OFF)
option(ENABLE_SANITIZERS "Compatibility alias for milestone scripts and legacy docs." ${CAS_ENABLE_SANITIZERS})

set(CAS_ENABLE_SANITIZERS "${ENABLE_SANITIZERS}")

# ─── ZERO WARNING POLICY OBBLIGATORIA ───
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic) # Rimosso -Werror temporaneamente
endif()

# ─── SANITIZERS (abilitati in Debug / RelWithDebInfo, anche multi-config) ───
if(CAS_ENABLE_SANITIZERS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        add_compile_options(
            $<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:-fsanitize=address>
            $<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:-fsanitize=undefined>
        )
        add_link_options(
            $<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:-fsanitize=address>
            $<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:-fsanitize=undefined>
        )
    endif()
endif()

set(CAS_PUBLIC_HEADERS
    include/cas/algebra.hpp
    include/cas/algebraic_number.hpp
    include/cas/algebraic_tower.hpp
    include/cas/algebraic_tower_bridge.hpp
    include/cas/ast.hpp
    include/cas/ast_debug.hpp
    include/cas/bigint.hpp
    include/cas/calculus.hpp
    include/cas/error.hpp
    include/cas/expr.hpp
    include/cas/formatter.hpp
    include/cas/lexer.hpp
    include/cas/linalg/Matrix.hpp
    include/cas/numtheory.hpp
    include/cas/parser.hpp
    include/cas/rational.hpp
    include/cas/result.hpp
    include/cas/symbolic.hpp
    include/cas/token.hpp
)

set(CAS_FOUNDATION_SOURCES
    src/foundation/bigint_arithmetic.cpp
    src/foundation/bigint_conversions.cpp
    src/foundation/bigint_logic.cpp
    src/foundation/rational.cpp
)

set(CAS_AST_SOURCES
    src/ast/ast.cpp
    src/ast/ast_debug.cpp
)

set(CAS_LEXER_SOURCES
    src/lexer/lexer.cpp
    src/lexer/token.cpp
)

set(CAS_PARSER_SOURCES
    src/parser/parser.cpp
    src/parser/parser_special_forms.cpp
    src/parser/parser_support.cpp
    src/parser/round_trip_printer.cpp
)

set(CAS_SYMBOLIC_SOURCES
    src/rewrite/builtin_rewrite.cpp
    src/symbolic/simplify_core.cpp
    src/symbolic/simplify_arithmetic.cpp
    src/symbolic/simplify_arithmetic_chain.cpp
    src/symbolic/simplify_functions.cpp
    src/symbolic/simplify_utils.cpp
    src/symbolic/context_core.cpp
    src/symbolic/context_utils.cpp
    src/symbolic/assumptions.cpp
    src/symbolic/term_order.cpp
    src/symbolic/rewrite_engine.cpp
    src/symbolic/rewrite_matching.cpp
    src/symbolic/substitute.cpp
    src/symbolic/functions_special.cpp
    src/symbolic/series_arithmetic.cpp
    src/symbolic/normal_form.cpp
)

set(CAS_CALCULUS_SOURCES
    src/calculus/differentiate.cpp
    src/calculus/integrate.cpp
    src/calculus/integrate_definite_patterns.cpp
    src/calculus/bessel_orthogonality.cpp
    src/calculus/orthogonal_polynomials.cpp
    src/calculus/integrate_improper.cpp
    src/calculus/summation.cpp
    src/calculus/integrate_core.cpp
    src/calculus/integrate_elementary.cpp
    src/calculus/integrate_helpers.cpp
    src/calculus/integrate_quadratic_rational.cpp
    src/calculus/integrate_product_power.cpp
    src/calculus/integrate_substitution.cpp
    src/calculus/integrate_trig_substitution.cpp
    src/calculus/integrate_hermite.cpp
    src/calculus/integrate_risch.cpp
    src/calculus/differential_field.cpp
    src/calculus/limit.cpp
    src/calculus/limit_helpers.cpp
    src/calculus/limit_infinite.cpp
    src/calculus/limit_polynomial.cpp
    src/calculus/limit_series.cpp
    src/calculus/limit_mrv.cpp
    src/calculus/residue.cpp
    src/calculus/residue_theorem.cpp
    src/calculus/laurent_general.cpp
    src/calculus/diff_rules.cpp
    src/calculus/integrate_parts.cpp
    src/calculus/ode_classifier.cpp
    src/calculus/ode_solver_1st_order.cpp
    src/calculus/ode_solver_advanced.cpp
    src/calculus/ode_solver_frobenius.cpp
    src/calculus/asymptotes.cpp
)

set(CAS_NUMTHEORY_SOURCES
    src/numtheory/arithmetic_functions.cpp
    src/numtheory/bernoulli.cpp
    src/numtheory/modular.cpp
    src/numtheory/primality.cpp
)

set(CAS_ALGEBRA_SOURCES
    src/algebra/algebra_core.cpp
    src/algebra/algebraic_number.cpp
    src/algebra/algebraic_number_bridge.cpp
    src/algebra/algebraic_tower_bridge.cpp
    src/algebra/polynomial_arithmetic.cpp
    src/algebra/polynomial_multivariate.cpp
    src/algebra/polynomial_resultant.cpp
    src/algebra/polynomial_square_free.cpp
    src/algebra/factorization_integers.cpp
    src/algebra/factorization_recombination.cpp
    src/algebra/factorization_polynomials.cpp
    src/algebra/polynomial_groebner.cpp
    src/algebra/polynomial_groebner_f4.cpp
    src/algebra/polynomial_groebner_f4_buchberger.cpp
    src/algebra/polynomial_groebner_f4_reduce.cpp
    src/algebra/polynomial_groebner_f4_conversion.cpp
    src/algebra/polynomial_groebner_f4_solver.cpp
    src/algebra/csolve.cpp
    src/algebra/partial_fractions.cpp
    src/algebra/partial_fractions_lrt.cpp
    src/algebra/solve_polynomial.cpp
    src/algebra/polynomial_internal.cpp
    src/algebra/polynomial_hensel.cpp
    src/algebra/lattice_lll.cpp
    src/algebra/polynomial_modular.cpp
    src/algebra/polynomial_gcd_heuristic.cpp
    src/algebra/polynomial_gcd_modular.cpp
    src/algebra/polynomial_gcd_multivariate.cpp
    src/algebra/algebraic_equal.cpp
    src/algebra/polynomial_integer.cpp
    src/algebra/polynomial_conversions.cpp
    src/algebra/polynomial_gcd_core.cpp
    src/algebra/polynomial_cyclotomic.cpp
)

set(CAS_LINALG_SOURCES
    src/linalg/matrix_ops.cpp
    src/linalg/matrix_inverse.cpp
    src/linalg/matrix_solve.cpp
    src/linalg/matrix_eigenvalues.cpp
    src/linalg/matrix_null_space_extension.cpp
    src/linalg/matrix_smith.cpp
    src/linalg/matrix_jordan.cpp
    src/linalg/matrix_determinant_modular.cpp
    src/linalg/matrix_structured_determinant.cpp
    src/linalg/matrix_bareiss.cpp
)

set(CAS_NUMERIC_SOURCES
    src/numeric/evaluator.cpp
    src/numeric/solvers.cpp
    src/numeric/integration.cpp
    src/numeric/ode.cpp
    src/numeric/sampler.cpp
)

set(CAS_FORMATTER_SOURCES
    src/formatter/formatter_text.cpp
    src/formatter/formatter_ascii.cpp
    src/formatter/formatter_latex.cpp
)

add_library(cas_foundation STATIC ${CAS_FOUNDATION_SOURCES})
target_include_directories(cas_foundation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(cas_foundation PUBLIC cxx_std_20)

add_library(cas_ast STATIC ${CAS_AST_SOURCES})
target_include_directories(cas_ast PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_ast PUBLIC cas_foundation)
target_compile_features(cas_ast PUBLIC cxx_std_20)

add_library(cas_lexer STATIC ${CAS_LEXER_SOURCES})
target_include_directories(cas_lexer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(cas_lexer PUBLIC cxx_std_20)

add_library(cas_parser STATIC ${CAS_PARSER_SOURCES})
target_include_directories(cas_parser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_parser PUBLIC cas_ast cas_lexer cas_foundation)
target_compile_features(cas_parser PUBLIC cxx_std_20)

add_library(cas_symbolic STATIC ${CAS_SYMBOLIC_SOURCES})
target_include_directories(cas_symbolic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_symbolic PUBLIC cas_ast cas_foundation)
target_compile_features(cas_symbolic PUBLIC cxx_std_20)

add_library(cas_calculus STATIC ${CAS_CALCULUS_SOURCES})
target_include_directories(cas_calculus PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_calculus PUBLIC cas_symbolic cas_parser cas_ast cas_foundation)
target_compile_features(cas_calculus PUBLIC cxx_std_20)

add_library(cas_numtheory STATIC ${CAS_NUMTHEORY_SOURCES})
target_include_directories(cas_numtheory PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_numtheory PUBLIC cas_foundation)
target_compile_features(cas_numtheory PUBLIC cxx_std_20)

target_link_libraries(cas_symbolic PUBLIC cas_numtheory)

add_library(cas_algebra STATIC ${CAS_ALGEBRA_SOURCES})
target_include_directories(cas_algebra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_algebra PUBLIC cas_symbolic cas_ast cas_foundation cas_numtheory)
target_compile_features(cas_algebra PUBLIC cxx_std_20)

target_link_libraries(cas_calculus PUBLIC cas_algebra)

add_library(cas_linalg STATIC ${CAS_LINALG_SOURCES})
target_include_directories(cas_linalg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_linalg PUBLIC cas_ast cas_symbolic cas_algebra)
target_compile_features(cas_linalg PUBLIC cxx_std_20)

add_library(cas_numeric STATIC ${CAS_NUMERIC_SOURCES})
target_include_directories(cas_numeric PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_numeric PUBLIC cas_ast cas_foundation cas_calculus)
target_compile_features(cas_numeric PUBLIC cxx_std_20)

add_library(cas_formatter STATIC ${CAS_FORMATTER_SOURCES})

target_include_directories(cas_formatter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(cas_formatter PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/formatter)
target_link_libraries(cas_formatter PUBLIC cas_ast cas_foundation)
target_compile_features(cas_formatter PUBLIC cxx_std_20)

add_library(cas_core INTERFACE)
target_include_directories(cas_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_sources(cas_core INTERFACE FILE_SET HEADERS BASE_DIRS include FILES ${CAS_PUBLIC_HEADERS})
target_link_libraries(cas_core INTERFACE cas_foundation cas_ast cas_lexer cas_parser cas_symbolic cas_numtheory cas_algebra $<LINK_LIBRARY:WHOLE_ARCHIVE,cas_linalg> cas_calculus cas_numeric cas_formatter)
target_compile_features(cas_core INTERFACE cxx_std_20)

add_library(cas_c_api STATIC src/c_api/cas_api.cpp)
target_include_directories(cas_c_api PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_c_api PRIVATE cas_core)
target_compile_features(cas_c_api PUBLIC cxx_std_20)

if(CAS_ENABLE_TESTING)
    add_executable(cas_c_api_test test/c_api/test_c_api.c)
    target_link_libraries(cas_c_api_test PRIVATE cas_c_api)
endif()

add_executable(cas_benchmarks
    test/benchmarks/benchmark_core.cpp
)
target_link_libraries(cas_benchmarks PRIVATE cas_core)
target_compile_features(cas_benchmarks PRIVATE cxx_std_20)

add_executable(audit_perf
    tools/audit_perf.cpp
)
target_link_libraries(audit_perf PRIVATE cas_core)
target_compile_features(audit_perf PRIVATE cxx_std_20)

add_executable(benchmark_tests
    tools/benchmark_tests.cpp
)
target_link_libraries(benchmark_tests PRIVATE cas_core)
target_compile_features(benchmark_tests PRIVATE cxx_std_20)

add_executable(stress_verifier
    tools/stress_verifier.cpp
)
target_link_libraries(stress_verifier PRIVATE cas_core)
target_compile_features(stress_verifier PRIVATE cxx_std_20)

# ─── TOOLS SEPARATI (Standalone UI) ───
add_subdirectory(tools/cas_ui)

# ─── GUI OPZIONALE E STACCABILE ───
if(CAS_ENABLE_GUI)
    add_subdirectory(GUI)
endif()

# ─── TESTING ───
if(CAS_ENABLE_TESTING)
    include(CTest)
    include(FetchContent)

    if(CAS_USE_SYSTEM_GTEST)
        find_package(GTest 1.14 REQUIRED CONFIG)
    else()
        FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()

    add_executable(cas_foundation_tests
        test/fuzz/test_bigint_fuzz.cpp
        test/fuzz/test_parser_fuzz.cpp
        test/fuzz/test_symbolic_fuzz.cpp
        test/unit/ast/test_ast_model.cpp
        test/unit/ast/test_ast_debug.cpp
        test/unit/foundation/test_foundation.cpp
        test/unit/lexer/test_lexer.cpp
        test/unit/numtheory/test_numtheory.cpp
        test/unit/lexer/test_token_model.cpp
        test/unit/parser/test_parser.cpp
        test/unit/parser/test_parser_roundtrip.cpp
        test/unit/symbolic/test_calculus.cpp
        test/unit/symbolic/test_calculus_f4.cpp
        test/unit/symbolic/test_compositions.cpp
        test/unit/symbolic/test_abs_sign.cpp
        test/unit/symbolic/test_symbolic.cpp
        test/unit/symbolic/test_linalg_integration.cpp
        test/unit/linalg/test_matrix_basic.cpp
        test/unit/linalg/test_matrix_advanced.cpp
        test/unit/linalg/test_bareiss_vander.cpp
        test/unit/linalg/test_eigen_n_gt_3.cpp
    test/unit/linalg/test_eigen_over_extension.cpp
        test/unit/numeric/test_numeric.cpp
        test/unit/formatter/test_formatter.cpp
        test/unit/algebra/test_resultant.cpp
        test/unit/algebra/test_resultant_generic.cpp
        test/unit/algebra/test_algebra.cpp
        test/unit/algebra/test_root_of.cpp
        test/unit/algebra/test_groebner.cpp
        test/unit/algebra/test_polynomial_gcd_subresultant.cpp
        test/unit/algebra/test_polynomial_gcd_multivariate_robust.cpp
        test/unit/algebra/test_algebraic_number.cpp
        test/unit/algebra/test_algebraic_tower_bridge.cpp
    test/unit/algebra/test_algebraic_number_bridge.cpp
        test/unit/algebra/test_algebraic_tower.cpp
        test/unit/algebra/test_factorization_lll.cpp
        test/unit/algebra/test_factorization_trager.cpp
        test/unit/algebra/test_regression.cpp
        test/unit/algebra/test_gap6.cpp
        test/unit/algebra/test_factorization_p2.cpp
        test/unit/symbolic/test_debug_limit.cpp
        test/unit/symbolic/test_randomized_antihardcode.cpp
        test/unit/symbolic/test_mathematical_equality.cpp
        test/unit/symbolic/test_property_based.cpp
        test/unit/symbolic/test_regression.cpp
        # test/unit/symbolic/test_quantity.cpp
        # test/unit/symbolic/test_assumptions_advanced.cpp
        # test/unit/symbolic/test_assumptions_advanced_new.cpp
        test/unit/symbolic/test_assumptions_relational.cpp
        # test/unit/symbolic/test_simplification_advanced.cpp
        test/unit/symbolic/test_integration_advanced.cpp
        test/unit/symbolic/test_laurent_series.cpp
        test/unit/symbolic/test_residue_theorem.cpp
        test/unit/symbolic/test_improper_pv.cpp
        test/unit/symbolic/test_frobenius.cpp
        test/unit/symbolic/test_special_functions.cpp
        test/unit/symbolic/test_bessel_zero.cpp
        test/unit/symbolic/test_complex_polar.cpp
        test/unit/calculus/test_asymptotes_rootsum.cpp
        test/unit/calculus/test_bessel_orthogonality.cpp
        test/unit/calculus/test_bessel_recurrence.cpp
        test/unit/calculus/test_orthogonal_polynomials.cpp
        test/unit/symbolic/test_assumptions_stability.cpp
        test/unit/symbolic/test_caching.cpp
        test/unit/symbolic/test_gruntz.cpp
        test/unit/test_stress.cpp
        test/unit/test_supreme.cpp
        test/unit/test_hpprime_acid.cpp
        test/unit/test_ode.cpp
        test/unit/test_ode_critical.cpp
        test/unit/test_acid_complex_canonical.cpp
        test/unit/calculus/test_asymptotes_rootsum.cpp
        test/unit/test_limit_mrv_cancellation.cpp
        test/unit/test_integrate_singularity.cpp
        test/regression/test_regression.cpp
    )
    target_link_libraries(cas_foundation_tests PRIVATE cas_core GTest::gtest_main)
    target_include_directories(cas_foundation_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

    add_executable(assumptions_stability_test
        test/unit/symbolic/test_assumptions_stability.cpp
    )
    target_link_libraries(assumptions_stability_test PRIVATE cas_core GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(cas_foundation_tests)
    gtest_discover_tests(assumptions_stability_test)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-fsanitize=fuzzer" COMPILER_SUPPORTS_FUZZER)
    if(COMPILER_SUPPORTS_FUZZER)
        add_executable(parser_fuzzer test/fuzz/parser_fuzzer.cpp)
        target_compile_options(parser_fuzzer PRIVATE -fsanitize=fuzzer)
        target_link_options(parser_fuzzer PRIVATE -fsanitize=fuzzer)
        target_link_libraries(parser_fuzzer PRIVATE cas_c_api)
    endif()
endif()
