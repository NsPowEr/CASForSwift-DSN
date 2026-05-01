# Test Report - 2026-04-24

Role: A-07 CAS Testing Agent.

Scope:
- Verify current CMake build and CTest suite.
- Do not modify production code.
- Do not weaken structural/math test oracles.

Commands:

```bash
cmake --build build
ctest --test-dir build --output-on-failure
ctest --test-dir build --rerun-failed --output-on-failure
```

Results:
- `cmake --build build`: PASS.
- Full `ctest`: FAIL, 357/365 passed, 8 failed.
- `ctest --rerun-failed`: FAIL, 8/8 failed again.

Confirmed failures:
- `CalculusIntegrateTest.IntegratesExtendedElementaryTable`
  - Error: `Only powers of the integration variable with constant exponent are implemented`.
- `CalculusIntegrateTest.MatchesPrimitiveFormsRequiringAbsoluteValuesOrReciprocalNormalization`
  - Error: `Only powers of the integration variable with constant exponent are implemented`.
- `CalculusLimitTest.ComputesDirectSubstitutionAndLHopitalCases`
  - Error: `Rational division by zero` for `limit(sin(x)/x, x, 0)`.
- `CalculusLimitTest.ComputesBasicInfiniteGrowthComparisons`
  - Returns unresolved expressions at infinity instead of expected constants/infinities.
- `CalculusLimitTest.RejectsLogarithmicGrowthAtNegativeInfinityOutsideDomain`
  - Expected domain error, actual result is success.
- `CalculusLimitTest.ComputesLateralPoleLimitsForRationalFunctions`
  - Error: `Rational division by zero` for lateral pole handling.
- `SymbolicSimplifyTest.ExpandsSqrtOfPositiveQuotient`
  - Structural mismatch for positive-assumption `sqrt(a / b) -> sqrt(a) / sqrt(b)`.
- `SymbolicRewriteTest.ContextUsesBuiltinProviderForPositiveSqrtQuotientRule`
  - Structural mismatch for provider-backed positive-assumption sqrt quotient rewrite.

Assessment:
- No local test-only fix was applied. The failing assertions are valid structural or mathematical contracts, not formatter oracles.
- Failures point to production behavior in calculus limit/integrate and symbolic rewrite/simplify, outside this task's safe ownership.
- The result is consistent with the roadmap/TODO hardening debt around calculus closure and rewrite assumptions.
