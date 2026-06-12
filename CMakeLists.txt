cmake_minimum_required(VERSION 3.25)
project(CASEngine VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ─── APPLE SILICON E IOS ARCHITECTURE SUPPORT ───
if(APPLE AND NOT CMAKE_OSX_ARCHITECTURES)
    # Se CMake gira sotto Rosetta 2 (uname restituisce x86_64), force ARM se c'è l'Homebrew di Apple Silicon
    # o se è una compilazione standard moderna per garantire il funzionamento su iPadOS/iOS/macOS M-series.
    if(EXISTS "/opt/homebrew/bin/brew")
        set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for Apple platforms" FORCE)
    else()
        execute_process(COMMAND uname -m OUTPUT_VARIABLE HOST_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CMAKE_OSX_ARCHITECTURES "${HOST_ARCH}" CACHE STRING "Build architectures for Apple platforms" FORCE)
    endif()
endif()
option(CAS_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer in debug-oriented builds." ON)
option(CAS_ENABLE_TESTING "Enable unit tests and CTest integration." ON)
option(CAS_ENABLE_GUI "Build the optional detachable manual GUI lab." OFF)
option(CAS_USE_SYSTEM_GTEST "Use an installed GoogleTest package instead of FetchContent." OFF)
option(ENABLE_SANITIZERS "Compatibility alias for milestone scripts and legacy docs." ${CAS_ENABLE_SANITIZERS})
option(ENABLE_COVERAGE "Enable gcov/lcov line coverage instrumentation (gcc/clang only). Disables optimizations." OFF)

set(CAS_ENABLE_SANITIZERS "${ENABLE_SANITIZERS}")

# ─── COVERAGE INSTRUMENTATION ───────────────────────────────────────────────
# Enable with: cmake -DENABLE_COVERAGE=ON -DCMAKE_BUILD_TYPE=Debug
# Then build and run: cmake --build <dir> --target coverage
if(ENABLE_COVERAGE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        add_compile_options(--coverage -fprofile-arcs -ftest-coverage -O0 -g)
        add_link_options(--coverage)
    else()
        message(WARNING "ENABLE_COVERAGE is ON but compiler is not Clang or GCC — ignored.")
    endif()
endif()

# ─── ZERO WARNING POLICY OBBLIGATORIA ───
# CLAUDE.md § REGOLA ZERO requires zero warning compilation.
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
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
    include/cas/gaussian_int.hpp
    include/cas/complex_rational.hpp
    include/cas/galois.hpp
    include/cas/interval.hpp
    include/cas/lexer.hpp
    include/cas/linalg/Matrix.hpp
    include/cas/numtheory.hpp
    include/cas/parser.hpp
    include/cas/rational.hpp
    include/cas/result.hpp
    include/cas/solve_inequality.hpp
    include/cas/statistics.hpp
    include/cas/symbolic.hpp
    include/cas/token.hpp
    include/cas/units.hpp
)

set(CAS_FOUNDATION_SOURCES
    src/foundation/bigint_arithmetic.cpp
    src/foundation/bigint_conversions.cpp
    src/foundation/bigint_logic.cpp
    # F1.1 production-grade BigInt algorithms
    src/foundation/bigint_mul_toom3.cpp
    src/foundation/bigint_div_knuth_d.cpp
    src/foundation/bigint_div_burnikel_ziegler.cpp
    src/foundation/bigint_gcd_lehmer.cpp
    src/foundation/bigint_numtheory.cpp
    src/foundation/bigint_factor_pollard.cpp
    src/foundation/rational.cpp
    src/foundation/complex_rational.cpp

    src/foundation/error.cpp
)

set(CAS_AST_SOURCES
    src/ast/ast.cpp          # AstArena lifecycle only (F8.0 split)
    src/ast/ast_compare.cpp  # structural_equal, expr_hash, expr_kind_name (F8.0 split)
    src/ast/ast_clone.cpp    # clone_into_arena (F8.0 split)
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
    src/rewrite/discrimination_net.cpp
    src/symbolic/simplify_core.cpp
    src/symbolic/simplify_arithmetic.cpp
    src/symbolic/simplify_arithmetic_power.cpp
    src/symbolic/simplify_arithmetic_chain_liketerm.cpp
    src/symbolic/summation_gosper.cpp
    src/symbolic/summation_zeilberger.cpp
    src/symbolic/summation_zeilberger_helpers.cpp
    src/symbolic/simplify_arithmetic_chain_sum.cpp
    src/symbolic/simplify_arithmetic_chain.cpp
    src/symbolic/simplify_arithmetic_chain_gamma.cpp
    src/symbolic/simplify_arithmetic_chain_sqrt.cpp
    src/symbolic/simplify_arithmetic_chain_exp.cpp
    src/symbolic/simplify_extended_real.cpp
    src/symbolic/simplify_functions.cpp
    src/symbolic/simplify_trig_chebyshev.cpp
    src/symbolic/simplify_trig_tables.cpp
    src/symbolic/simplify_trig.cpp
    src/symbolic/simplify_trig_inverse.cpp
    src/symbolic/simplify_hyperbolic.cpp
    src/symbolic/simplify_combinatorial.cpp
    src/symbolic/simplify_exp_log.cpp
    src/symbolic/simplify_special_fn.cpp
    src/symbolic/simplify_special_fn_hyper.cpp
    src/symbolic/simplify_orthogonal_polys.cpp
    src/symbolic/simplify_bessel_orthogonal.cpp
    src/symbolic/simplify_complex.cpp
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
    src/symbolic/units.cpp
    src/symbolic/complex_qi.cpp
)

set(CAS_CALCULUS_SOURCES
    src/calculus/differentiate.cpp
    src/calculus/integrate.cpp
    src/calculus/integrate_definite_patterns.cpp
    src/calculus/bessel_orthogonality.cpp
    src/calculus/orthogonal_polynomials.cpp
    src/calculus/integrate_improper.cpp
    src/calculus/summation.cpp
    src/calculus/summation_abramov.cpp
    src/calculus/summation_zeilberger_driver.cpp
    src/calculus/integrate_core.cpp
    src/calculus/integrate_elementary.cpp
    src/calculus/integrate_helpers.cpp
    src/calculus/integrate_quadratic_rational.cpp
    src/calculus/integrate_product_power.cpp
    src/calculus/integrate_substitution.cpp
    src/calculus/integrate_trig_substitution.cpp
    src/calculus/integrate_weierstrass.cpp
    src/calculus/integrate_hermite.cpp
    src/calculus/integrate_risch.cpp
    src/calculus/integrate_risch_rde.cpp
    src/calculus/integrate_risch_hermite.cpp
    src/calculus/risch_parametric.cpp
    src/calculus/risch_de_extensions.cpp
    src/calculus/risch_rde_bronstein.cpp
    src/calculus/risch_log_expand.cpp
    src/calculus/risch_exp_expand.cpp
    src/calculus/z_transform.cpp
    src/calculus/mellin_transform.cpp
    src/calculus/inverse_laplace_residue.cpp
    src/calculus/integrate_dirac.cpp
    src/calculus/fourier_transform.cpp
    src/calculus/differential_field.cpp
    src/calculus/limit.cpp
    src/calculus/limit_helpers.cpp
    src/calculus/limit_infinite.cpp
    src/calculus/limit_polynomial.cpp
    src/calculus/limit_series.cpp
    src/calculus/limit_mrv.cpp
    src/calculus/limit_mrv_compare.cpp
    src/calculus/limit_mrv_set.cpp
    src/calculus/limit_mrv_leading.cpp
    src/calculus/limit_mrv_exp.cpp
    src/calculus/limit_quotient_d2.cpp
    src/calculus/residue.cpp
    src/calculus/residue_theorem.cpp
    src/calculus/laurent_general.cpp
    src/calculus/pade.cpp
    src/calculus/puiseux.cpp
    src/calculus/multiple_integral.cpp
    src/calculus/diff_rules.cpp
    src/calculus/integrate_parts.cpp
    src/calculus/ode_classifier.cpp
    src/calculus/ode_classifier_first_order.cpp
    src/calculus/ode_classifier_higher_order.cpp
    src/calculus/ode_solver_1st_order.cpp
    src/calculus/ode_solver_advanced.cpp
    src/calculus/ode_kovacic_case1.cpp
    src/calculus/ode_kovacic_laurent.cpp
    src/calculus/ode_kovacic.cpp
    src/calculus/ode_solver_frobenius.cpp
    src/calculus/ode_solver_nonlinear.cpp
    src/calculus/asymptotes.cpp
    src/calculus/numeric_diff.cpp
    src/calculus/laplace_transform.cpp
    src/calculus/ode_laplace.cpp
)

set(CAS_NUMTHEORY_SOURCES
    src/numtheory/arithmetic_functions.cpp
    src/numtheory/bernoulli.cpp
    src/numtheory/modular.cpp
    src/numtheory/primality.cpp
)

set(CAS_ALGEBRA_SOURCES
    src/algebra/algebra_core.cpp
    src/algebra/algebra_polynomial_api.cpp
    src/algebra/polynomial_sparse_interpolation.cpp
    src/algebra/algebraic_number.cpp
    src/algebra/algebraic_number_qx.cpp
    src/algebra/algebraic_number_bridge.cpp
    src/algebra/algebraic_tower_bridge.cpp
    src/algebra/algebraic_tower_primitive.cpp
    src/algebra/algebraic_tower_primitive_nested.cpp
    src/algebra/polynomial_arithmetic.cpp
    src/algebra/polynomial_multivariate.cpp
    src/algebra/factor_multivariate_mpoly.cpp
    src/algebra/factor_multivariate_convert.cpp
    src/algebra/factor_multivariate_lc.cpp
    src/algebra/factor_multivariate_hensel.cpp
    src/algebra/factor_multivariate_wang.cpp
    src/algebra/polynomial_resultant.cpp
    src/algebra/polynomial_square_free.cpp
    src/algebra/factorization_integers.cpp
    src/algebra/factorization_recombination.cpp
    src/algebra/van_hoeij_factor.cpp
    src/algebra/factorization_num_den.cpp
    src/algebra/factorization_poly_gcd.cpp
    src/algebra/factorization_trager.cpp
    src/algebra/factorization_berlekamp.cpp
    src/algebra/polynomial_half_gcd.cpp
    src/algebra/factorization_tower.cpp
    src/algebra/factorization_tower_helpers.cpp
    src/algebra/factorization_tower_n.cpp
    src/algebra/algebraic_number_rootof.cpp
    src/algebra/polynomial_groebner.cpp
    src/algebra/polynomial_groebner_f4.cpp
    src/algebra/polynomial_groebner_f4_buchberger.cpp
    src/algebra/polynomial_groebner_f4_reduce.cpp
    src/algebra/polynomial_groebner_f4_conversion.cpp
    src/algebra/polynomial_groebner_f4_solver.cpp
    src/algebra/polynomial_groebner_fglm.cpp
    src/algebra/polynomial_groebner_f5.cpp
    src/algebra/csolve.cpp
    src/algebra/fsolve.cpp
    src/algebra/solve_inequality.cpp
    src/algebra/partial_fractions.cpp
    src/algebra/partial_fractions_lrt.cpp
    src/algebra/solve_polynomial.cpp
    src/algebra/polynomial_internal.cpp
    src/algebra/polynomial_hensel.cpp
    src/algebra/lattice_lll.cpp
    src/algebra/polynomial_modular.cpp
    src/algebra/polynomial_gcd_heuristic.cpp
    src/algebra/polynomial_gcd_modular.cpp
    src/algebra/polynomial_gcd_brown.cpp
    src/algebra/polynomial_gcd_brown_modular.cpp
    src/algebra/polynomial_gcd_brown_lc_scaling.cpp
    src/algebra/polynomial_gcd_fp_recursive.cpp
    src/algebra/polynomial_gcd_zippel_prony.cpp
    src/algebra/polynomial_gcd_multivariate_helpers.cpp
    src/algebra/polynomial_gcd_multivariate_interp.cpp
    src/algebra/polynomial_gcd_multivariate_sparse.cpp
    src/algebra/polynomial_gcd_multivariate_linear.cpp
    src/algebra/polynomial_gcd_multivariate.cpp
    src/algebra/algebraic_equal.cpp
    src/algebra/algebraic_equal_cyclotomic.cpp
    src/algebra/hyperbolic_normalize.cpp
    src/algebra/polynomial_integer.cpp
    src/algebra/polynomial_conversions.cpp
    src/algebra/polynomial_gcd_core.cpp
    src/algebra/polynomial_gcd_crt.cpp
    src/algebra/polynomial_resultant_crt.cpp
    src/algebra/polynomial_cyclotomic.cpp
    src/algebra/gaussian_int.cpp
    src/algebra/gaussian_factor.cpp
    src/algebra/galois.cpp
    src/algebra/galois_deg5.cpp
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
    src/linalg/matrix_lu.cpp
    src/linalg/matrix_qr.cpp
    src/linalg/matrix_expr_helpers.cpp
    src/linalg/matrix_cholesky.cpp
    src/linalg/matrix_determinant.cpp
    src/linalg/matrix_companion.cpp
    src/linalg/matrix_hermite.cpp
)

set(CAS_NUMERIC_SOURCES
    src/numeric/evaluator.cpp
    src/numeric/solvers.cpp
    src/numeric/sturm.cpp
    src/numeric/interpolation_spline.cpp
    src/numeric/ode_rkf45.cpp
    src/numeric/aberth_root_isolator.cpp
    src/numeric/lipschitz.cpp
    src/numeric/integration.cpp
    src/numeric/adaptive_integration.cpp
    src/numeric/interpolation.cpp
    src/statistics/normal.cpp
    src/statistics/regression.cpp
    src/statistics/discrete.cpp
    src/statistics/continuous.cpp
    src/statistics/hypothesis_tests.cpp
    src/statistics/multivariate_regression.cpp
    src/numeric/ode.cpp
    src/numeric/sampler.cpp
    src/numeric/bigfloat.cpp
    src/numeric/bigfloat_eval.cpp
    src/numeric/interval.cpp
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

# HPP-005 closure: cas_calculus uses BigFloat (MPFR) for boundary
# singularity checks in integrate.cpp.  Locate MPFR/GMP headers before
# the cas_calculus target so its compile flags can reference them.
find_path(MPFR_INCLUDE_DIR_EARLY mpfr.h
    HINTS /opt/homebrew/opt/mpfr/include
    REQUIRED)
find_path(GMP_INCLUDE_DIR_EARLY gmp.h
    HINTS /opt/homebrew/opt/gmp/include
    REQUIRED)

add_library(cas_calculus STATIC ${CAS_CALCULUS_SOURCES})
target_include_directories(cas_calculus PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_calculus PUBLIC cas_symbolic cas_parser cas_ast cas_foundation)
target_compile_features(cas_calculus PUBLIC cxx_std_20)
target_include_directories(cas_calculus PRIVATE
    ${MPFR_INCLUDE_DIR_EARLY}
    ${GMP_INCLUDE_DIR_EARLY})

add_library(cas_numtheory STATIC ${CAS_NUMTHEORY_SOURCES})
target_include_directories(cas_numtheory PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_numtheory PUBLIC cas_foundation)
target_compile_features(cas_numtheory PUBLIC cxx_std_20)

target_link_libraries(cas_symbolic PUBLIC cas_numtheory)

add_library(cas_algebra STATIC ${CAS_ALGEBRA_SOURCES})
target_include_directories(cas_algebra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(cas_algebra PRIVATE
    ${MPFR_INCLUDE_DIR_EARLY}
    ${GMP_INCLUDE_DIR_EARLY})
target_link_libraries(cas_algebra PUBLIC cas_symbolic cas_ast cas_foundation cas_numtheory)
target_compile_features(cas_algebra PUBLIC cxx_std_20)

target_link_libraries(cas_calculus PUBLIC cas_algebra)

add_library(cas_linalg STATIC ${CAS_LINALG_SOURCES})
target_include_directories(cas_linalg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cas_linalg PUBLIC cas_ast cas_symbolic cas_algebra)
target_compile_features(cas_linalg PUBLIC cxx_std_20)

# ── MPFR / GMP (required for BigFloat / L3-01) ────────────────────────────────
find_library(MPFR_LIB mpfr
    HINTS /opt/homebrew/opt/mpfr/lib
    REQUIRED)
find_library(GMP_LIB gmp
    HINTS /opt/homebrew/opt/gmp/lib
    REQUIRED)
find_path(MPFR_INCLUDE_DIR mpfr.h
    HINTS /opt/homebrew/opt/mpfr/include
    REQUIRED)
find_path(GMP_INCLUDE_DIR gmp.h
    HINTS /opt/homebrew/opt/gmp/include
    REQUIRED)

add_library(cas_numeric STATIC ${CAS_NUMERIC_SOURCES})
target_include_directories(cas_numeric PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MPFR_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR})
target_link_libraries(cas_numeric PUBLIC cas_ast cas_foundation cas_calculus
    ${MPFR_LIB} ${GMP_LIB})
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

add_library(cas_gui_bridge STATIC GUI/src/CasGuiSession.cpp)
target_include_directories(cas_gui_bridge
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/GUI/src
)
target_link_libraries(cas_gui_bridge PUBLIC cas_core)
target_compile_features(cas_gui_bridge PUBLIC cxx_std_20)

if(CAS_ENABLE_TESTING)
    add_executable(cas_c_api_test test/c_api/test_c_api.c)
    target_link_libraries(cas_c_api_test PRIVATE cas_c_api)
endif()

add_executable(cas_benchmarks
    test/benchmarks/benchmark_core.cpp
)
target_link_libraries(cas_benchmarks PRIVATE cas_core)
target_compile_features(cas_benchmarks PRIVATE cxx_std_20)


# F0.6 Benchmark Regression CI Gate
add_custom_target(check_bench_regression
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_bench_regression.sh --current benchmark_results.txt --baseline ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/baseline_release.json
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Checking for benchmark regressions..."
)

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

# F7.4-C1 — Golden corpus runner (CAS engine vs Maxima diff harness).
add_executable(cas_golden_runner
    test/golden/main.cpp
)
target_link_libraries(cas_golden_runner PRIVATE cas_core)
target_include_directories(cas_golden_runner PRIVATE test/golden)
target_compile_features(cas_golden_runner PRIVATE cxx_std_20)

# ─── TOOLS SEPARATI (Standalone UI) ───
add_subdirectory(tools/cas_ui)

# ─── GUI OPZIONALE E STACCABILE ───
if(CAS_ENABLE_GUI)
    find_package(Qt6 6.7 REQUIRED COMPONENTS Core Gui Qml Quick QuickControls2 Concurrent Widgets)
    set(CMAKE_AUTOMOC ON)
    qt_standard_project_setup()
    if(COMMAND qt_policy)
        qt_policy(SET QTP0001 NEW)
        qt_policy(SET QTP0004 NEW)
    endif()

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
        # Suppress warnings on gtest/gmock targets: -Werror is for our code only.
        foreach(_gtest_target gtest gtest_main gmock gmock_main)
            if(TARGET ${_gtest_target})
                target_compile_options(${_gtest_target} PRIVATE -w)
            endif()
        endforeach()
    endif()

    add_executable(cas_foundation_tests
        test/fuzz/test_bigint_fuzz.cpp
        test/fuzz/test_parser_fuzz.cpp
        test/fuzz/test_symbolic_fuzz.cpp
        test/unit/error/test_diagnostic.cpp
        test/unit/error/test_unimplemented_payload.cpp

        test/unit/ast/test_ast_model.cpp
        test/unit/ast/test_ast_debug.cpp
        test/unit/foundation/test_foundation.cpp
        test/unit/foundation/test_bigint_production.cpp
        test/unit/lexer/test_lexer.cpp
        test/unit/numtheory/test_numtheory.cpp
        test/unit/lexer/test_token_model.cpp
        test/unit/parser/test_parser.cpp
        test/unit/parser/test_parser_roundtrip.cpp
        test/unit/symbolic/test_calculus.cpp
test/unit/symbolic/test_unwinding_number.cpp
        test/unit/symbolic/test_complex_simplification.cpp
        test/unit/symbolic/test_extended_real.cpp
        test/unit/symbolic/test_simplify_hyperbolic_zero.cpp
        test/unit/symbolic/test_simplify_combinatorial.cpp
        test/unit/symbolic/test_special_fn_e1_identities.cpp
        test/unit/symbolic/test_bessel_e2_identities.cpp
        test/unit/symbolic/test_calculus_f4.cpp
        test/unit/symbolic/test_multiple_integral.cpp
        test/unit/symbolic/test_limit_tower_adaptive.cpp
        test/unit/symbolic/test_mrv_cancellation_probe.cpp
        test/unit/calculus/test_risch_logarithmic_probe.cpp
        test/unit/calculus/test_differential_field_probe.cpp
        test/unit/calculus/test_differentiate_smoke.cpp
        test/unit/calculus/test_limit_smoke.cpp
        test/unit/calculus/test_limit_gruntz_nested_log.cpp
        test/unit/calculus/test_integrate_interrupt.cpp
        test/unit/calculus/test_integrate_risch_smoke.cpp
        test/unit/calculus/test_integrate_substitution.cpp
        test/unit/algebra/test_polynomial_sparse_interpolation.cpp
        test/unit/algebra/test_algebraic_number_rootof.cpp
test/unit/algebra/test_cyclotomic_rootof_d2.cpp
        test/unit/algebra/test_csolve_linear.cpp
        test/unit/symbolic/test_summation_gosper.cpp
        test/unit/algebra/test_factorization_polynomials_smoke.cpp
        test/unit/algebra/test_polynomial_gcd_multivariate_smoke.cpp
        test/unit/parser/test_decimal_to_rational.cpp
        test/unit/symbolic/test_denesting_recursive.cpp
        test/unit/symbolic/test_rootof_auto_trigger.cpp
        test/unit/calculus/test_integrate_risch_exp_mix.cpp
        test/unit/calculus/test_risch_parametric.cpp
        test/unit/calculus/test_risch_de_extensions.cpp
        test/unit/calculus/test_risch_exp_expand.cpp
        test/unit/calculus/test_z_transform.cpp
        test/unit/calculus/test_mellin_transform.cpp
        test/unit/calculus/test_inverse_laplace_residue.cpp
        test/unit/calculus/test_dirac_heaviside_smoke.cpp
        test/unit/calculus/test_integrate_dirac.cpp
        test/unit/calculus/test_fourier_transform.cpp
        test/unit/calculus/test_hypergeometric_pfq.cpp
        test/unit/calculus/test_special_closed_forms.cpp
        test/unit/calculus/test_elliptic_integrals.cpp
        test/unit/numeric/test_mpfr_foundation.cpp
test/unit/numeric/test_sturm_bigfloat.cpp
        test/unit/algebra/test_gaussian_int.cpp
        test/unit/algebra/test_gaussian_factor.cpp
        test/unit/foundation/test_rational_cf.cpp
        test/unit/symbolic/test_complex_qi.cpp
        test/unit/ast/test_arena_shard_locks.cpp
        test/unit/ast/test_arena_reset.cpp
        test/unit/ast/test_arena_memory_budget.cpp
        test/unit/foundation/test_bigint_limb_budget.cpp
        test/unit/ast/test_arena_hash_dos.cpp
        test/unit/numeric/test_implicit_contour_vector.cpp
        test/unit/statistics/test_hypothesis_tests.cpp
        test/unit/numeric/test_spline_hermite_rkf45.cpp
        test/unit/symbolic/test_async_depth_scope.cpp
        test/unit/symbolic/test_cancellation.cpp
        test/unit/symbolic/test_cache_assumption_invalidation.cpp
        test/unit/symbolic/test_pow_div_simplify.cpp
        test/unit/algebra/test_swell_guard.cpp
        test/unit/calculus/test_weierstrass_substitution.cpp
        test/unit/symbolic/test_complex_log_branch.cpp
        test/unit/numeric/test_float_contextual.cpp
        test/unit/symbolic/test_units.cpp
        test/unit/symbolic/test_piecewise.cpp
        test/unit/algebra/test_monomial_orders.cpp
        test/unit/numeric/test_interval.cpp
        test/unit/numeric/test_aberth_root_isolator.cpp
        test/unit/calculus/test_numeric_diff.cpp
        test/unit/symbolic/test_branch_cuts_global.cpp
        test/unit/linalg/test_lu_decomposition.cpp
        test/unit/calculus/test_laplace.cpp
        test/unit/algebra/test_galois.cpp
        test/unit/algebra/test_galois_deg5_f3.cpp
        test/unit/algebra/test_solve_closure.cpp
        test/unit/calculus/test_ode_laplace.cpp
        test/unit/symbolic/test_jacobi.cpp
        test/unit/linalg/test_qr_decomposition.cpp
test/unit/linalg/test_cholesky_ldlt.cpp
test/unit/linalg/test_companion_matrix.cpp
test/unit/linalg/test_jordan_certificator.cpp
test/unit/linalg/test_special_determinants.cpp
test/unit/linalg/test_smith_qx.cpp
test/unit/linalg/test_hermite.cpp
test/unit/linalg/test_f4_stress.cpp
        test/unit/calculus/test_pfd_inverse_laplace.cpp
        test/unit/symbolic/test_lambert_w.cpp
        test/unit/algebra/test_cyclotomic_mobius.cpp
        test/unit/algebra/test_berlekamp_probe.cpp
        test/unit/algebra/test_hensel_quadratic.cpp
        test/unit/algebra/test_berlekamp_factorization.cpp
test/unit/algebra/test_half_gcd.cpp
        test/unit/symbolic/test_math_correctness.cpp
        test/unit/symbolic/test_compositions.cpp
        test/unit/symbolic/test_abs_sign.cpp
        test/unit/symbolic/test_trig_identities.cpp
        test/unit/symbolic/test_chebyshev_trig.cpp
        test/unit/rewrite/test_discrimination_net.cpp
        test/unit/symbolic/test_complex_log.cpp
        test/unit/symbolic/test_symbolic.cpp
        test/unit/symbolic/test_linalg_integration.cpp
        test/unit/golden/test_matrix_adapter_d2.cpp
        test/unit/linalg/test_matrix_basic.cpp
        test/unit/linalg/test_matrix_advanced.cpp
        test/unit/linalg/test_bareiss_vander.cpp
        test/unit/linalg/test_eigen_n_gt_3.cpp
    test/unit/linalg/test_eigen_over_extension.cpp
        test/unit/numeric/test_numeric.cpp
        test/unit/numeric/test_bigfloat.cpp
        test/unit/numeric/test_gauss_kronrod.cpp
        test/unit/numeric/test_interpolation.cpp
        test/unit/numeric/test_parametric_sampler.cpp
        test/unit/statistics/test_normal.cpp
        test/unit/statistics/test_regression.cpp
        test/unit/statistics/test_discrete.cpp
        test/unit/statistics/test_continuous.cpp
        test/unit/formatter/test_formatter.cpp
        test/unit/gui/test_cas_gui_session.cpp
        test/unit/algebra/test_resultant.cpp
        test/unit/algebra/test_resultant_generic.cpp
        test/unit/algebra/test_algebra.cpp
        test/unit/algebra/test_root_of.cpp
        test/unit/algebra/test_groebner.cpp
        test/unit/algebra/test_groebner_f5_pruning_f3.cpp
        test/unit/algebra/test_polynomial_gcd_subresultant.cpp
        test/unit/algebra/test_polynomial_gcd_multivariate_robust.cpp
        test/unit/algebra/test_gcd_brown_f3.cpp
        test/unit/algebra/test_factor_multivariate_f3.cpp
        test/unit/algebra/test_modular_crt.cpp
        test/unit/algebra/test_algebraic_number.cpp
        test/unit/algebra/test_algebraic_number_qx.cpp
        test/unit/algebra/test_algebraic_tower_bridge.cpp
    test/unit/algebra/test_algebraic_number_bridge.cpp
        test/unit/algebra/test_algebraic_tower.cpp
        test/unit/algebra/test_fsolve.cpp
        test/unit/algebra/test_solve_inequality.cpp
        test/unit/algebra/test_f2_gate_benchmark.cpp
        test/unit/algebra/test_factorization_lll.cpp
        test/unit/algebra/test_van_hoeij.cpp
        test/unit/algebra/test_factorization_trager.cpp
        test/unit/algebra/test_factorization_tower.cpp
    test/unit/algebra/test_factorization_tower_n_f3.cpp
        test/unit/algebra/test_primitive_element_f3.cpp
        test/unit/symbolic/test_equivalence_subset.cpp
        test/unit/algebra/test_regression.cpp
        test/unit/algebra/test_gap6.cpp
        test/unit/algebra/test_factorization_p2.cpp
        test/unit/symbolic/test_debug_limit.cpp
        test/unit/symbolic/test_randomized_antihardcode.cpp
        test/unit/symbolic/test_mathematical_equality.cpp
        test/unit/symbolic/test_property_based.cpp
        test/unit/symbolic/test_regression.cpp
        test/unit/symbolic/test_quantity.cpp
        test/unit/symbolic/test_assumptions_advanced.cpp
        test/unit/symbolic/test_assumptions_advanced_new.cpp
        test/unit/symbolic/test_assumptions_relational.cpp
        # test/unit/symbolic/test_assumptions_canonical_50.cpp  # file missing from disk (pre-existing stash issue)
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
        test/unit/calculus/test_pade_approximant.cpp
        test/unit/calculus/test_puiseux.cpp
        test/unit/calculus/test_residue_theorem_aberth.cpp
        test/unit/calculus/test_definite_summation.cpp
        test/unit/calculus/test_zeilberger_summation.cpp
        test/unit/symbolic/test_assumptions_stability.cpp
        test/unit/symbolic/test_caching.cpp
        test/unit/symbolic/test_gruntz.cpp
        test/unit/test_stress.cpp
        test/unit/test_supreme.cpp
        test/unit/test_hpprime_acid.cpp
        test/unit/test_ode.cpp
        test/unit/test_ode_critical.cpp
        test/unit/test_ode_class_f5.cpp
        test/unit/test_ode_solv_f5.cpp
        test/unit/test_acid_complex_canonical.cpp
        test/unit/calculus/test_asymptotes_rootsum.cpp
        test/unit/test_limit_mrv_cancellation.cpp
        test/unit/test_integrate_singularity.cpp
        test/regression/test_regression.cpp
        test/unit/symbolic/test_coverage_f1_debt5.cpp
        test/unit/symbolic/test_context_core_coverage.cpp
        test/unit/symbolic/test_orthogonal_coverage.cpp
        test/unit/symbolic/test_coverage_symbolic_v2.cpp
        test/unit/symbolic/test_coverage_symbolic_v3.cpp
    )
    target_link_libraries(cas_foundation_tests PRIVATE cas_core cas_gui_bridge GTest::gtest_main)
    target_include_directories(cas_foundation_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/GUI/src
        ${CMAKE_CURRENT_SOURCE_DIR}/test/golden
    )

    if(CAS_ENABLE_GUI)
        qt_add_executable(cas_gui_vm_tests
            GUI/test/test_appcore_undo_redo.cpp
            GUI/test/test_workspace_recovery.cpp
            GUI/test/test_workspace_layout_state.cpp
            GUI/test/test_autocomplete_vm.cpp
            GUI/test/test_variables_actions.cpp
            GUI/test/test_notebook_interaction_state.cpp
            GUI/test/test_command_availability.cpp
            GUI/test/test_command_palette_behavior.cpp
            GUI/test/test_session_sidebar_model.cpp
        )
        target_link_libraries(cas_gui_vm_tests
            PRIVATE
                cas_gui_vm
                GTest::gtest
        )
        target_include_directories(cas_gui_vm_tests PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/macos-cpp/src
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/macos-cpp/src/viewmodels
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/macos-cpp/src/core
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/src
        )
        add_test(NAME cas_gui_vm_tests COMMAND cas_gui_vm_tests)

        add_executable(cas_gui_qml_smoke_tests
            GUI/test/test_qml_smoke.cpp
        )
        target_link_libraries(cas_gui_qml_smoke_tests
            PRIVATE
                cas_gui_vmplugin
                cas_gui_vmplugin_init
                Qt6::Quick
                Qt6::Gui
                Qt6::Qml
                Qt6::Concurrent
                Qt6::Widgets
                GTest::gtest
        )
        target_include_directories(cas_gui_qml_smoke_tests PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/macos-cpp/src
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/macos-cpp/src/viewmodels
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/macos-cpp/src/qml-items
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/macos-cpp/src/core
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/GUI/src
        )
        add_test(NAME cas_gui_qml_smoke_tests COMMAND cas_gui_qml_smoke_tests)
    endif()

    add_executable(assumptions_stability_test
        test/unit/symbolic/test_assumptions_stability.cpp
    )
    target_link_libraries(assumptions_stability_test PRIVATE cas_core GTest::gtest_main)

    include(GoogleTest)
    # Per-test wall-clock cap (CLAUDE.md REGOLA TIMEOUT TEST): ogni test gtest
    # individuale viene killato da CTest oltre 60s. Protezione di rete contro
    # loop infiniti introdotti accidentalmente nel core algoritmico.
    gtest_discover_tests(cas_foundation_tests DISCOVERY_TIMEOUT 60
        PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0" TIMEOUT 60)
    gtest_discover_tests(assumptions_stability_test DISCOVERY_TIMEOUT 60
        PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0" TIMEOUT 60)
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

# ─── SANITIZER TARGET: cas_tests_san ────────────────────────────────────────
# Dedicated sanitizer build: compile cas_foundation_tests sources with
# ASan + UBSan always-on (independent of CAS_ENABLE_SANITIZERS build flag).
# TSan is added if the compiler supports it.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-fsanitize=thread" COMPILER_SUPPORTS_TSAN)

if(CAS_ENABLE_TESTING AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    # Interface library that carries ASan+UBSan flags.
    add_library(cas_san_flags INTERFACE)
    target_compile_options(cas_san_flags INTERFACE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
        -g -O1
    )
    target_link_options(cas_san_flags INTERFACE
        -fsanitize=address,undefined
    )
    if(COMPILER_SUPPORTS_TSAN)
        # TSan is incompatible with ASan; provide separate tsan variant if needed.
        # For now document availability without auto-enabling.
        message(STATUS "TSan supported by compiler (use -fsanitize=thread separately if needed).")
    endif()
endif()

# ─── RAPIDCHECK + PROPERTY TESTS (F0.4) ──────────────────────────────────────
if(CAS_ENABLE_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        rapidcheck
        GIT_REPOSITORY https://github.com/emil-e/rapidcheck.git
        GIT_TAG        b2d9ed2dddefc4b84318d664b4f221eb792d89c7  # 2025-01-latest HEAD
    )
    set(RC_ENABLE_GTEST ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(rapidcheck)
    # Suppress warnings on rapidcheck: -Werror is for our code only.
    if(TARGET rapidcheck)
        target_compile_options(rapidcheck PRIVATE -w)
    endif()

    add_subdirectory(test/property)
endif()

# ─── COVERAGE TARGET ─────────────────────────────────────────────────────────
# Requires ENABLE_COVERAGE=ON build. Run: cmake --build <dir> --target coverage
if(ENABLE_COVERAGE)
    find_program(LCOV_PATH lcov)
    find_program(GENHTML_PATH genhtml)
    if(LCOV_PATH AND GENHTML_PATH)
        add_custom_target(coverage
            COMMENT "Generating HTML coverage report in ${CMAKE_BINARY_DIR}/coverage/"
            COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR}
                    --zerocounters
            COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
                    --output-on-failure
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR}
                    --capture
                    --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    --ignore-errors mismatch
            COMMAND ${LCOV_PATH} --remove ${CMAKE_BINARY_DIR}/coverage.info
                    '/usr/*' '*/test/*' '*/_deps/*' '*/tools/*' '*/GUI/*'
                    --output-file ${CMAKE_BINARY_DIR}/coverage_filtered.info
                    --ignore-errors unused
            COMMAND ${GENHTML_PATH} ${CMAKE_BINARY_DIR}/coverage_filtered.info
                    --output-directory ${CMAKE_BINARY_DIR}/coverage
                    --title "CAS Engine Coverage"
                    --legend --demangle-cpp
            COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan
                    "Coverage report: ${CMAKE_BINARY_DIR}/coverage/index.html"
        )
    else()
        message(WARNING "lcov/genhtml not found — 'coverage' target unavailable. Install lcov.")
    endif()
endif()

# F0.6 / F7.0-A1.5 Anti-monolith CI Gate (HARD-FAIL on every build).
# Wired into ALL target so a violation breaks the build, not just an opt-in check.
# Whitelist managed in scripts/file_size_whitelist.txt — pre-existing debt only,
# no new entries without a corresponding split task in CAS_TASKS.md.
add_custom_target(check_anti_monolith ALL
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_file_size.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "F7.0-A1.5: Anti-monolith scan (files > 500 lines, whitelist-aware)..."
    VERBATIM
)
