---
name: hardcode-audit
description: Scansiona src/ e include/ per hardcode vietati dalle 10 categorie di CLAUDE.md. Produce report con file:riga, categoria e suggerimento fix. Claude-only.
user-invocable: false
---

Audit the CAS Engine codebase for forbidden hardcodes as defined in CLAUDE.md categories 1-10.

## Scan targets

Directories: `src/`, `include/`
File types: `*.cpp`, `*.hpp`, `*.h`

## Categories to check

For each finding, report: **file:line | category | snippet | fix**

### Cat 1 — Budget computazionali fissi
Grep patterns:
- `depth >= [0-9]` in integrate/simplify paths
- `constexpr.*MAX_.*= [0-9]` (non configurabile via CASContext)
- `subset_count >= [0-9]` in factorization
- `max_steps` con valore letterale

### Cat 2 — Costanti magiche algebriche
- `\+ [0-9]+U?;` vicino a `max_samples` o `required_samples`
- `= [0-9]{3,}` in algoritmi probabilistici (non in test)
- `score = [0-9]+` hardcoded (non derivato da assumptions)

### Cat 3 — Set e range fissi
- Array statici `kCandidates` o simili con primes/valori fissi
- `for.*n <= 100` o loop con bound letterale in algoritmi matematici
- Serie Taylor con liste fisse di funzioni

### Cat 6 — Seed e randomness
- `seed.*42` o `seed.*13` o `seed = [0-9]+`  
- Fallback a primo fisso (es. `p = 13`, `p = 17`) in Cantor-Zassenhaus

### Cat 7 — Nomi variabili interni fissi
- String literals `"__` prefixed in limit/ODE code
- `"C1"`, `"C2"` come nomi costanti ODE

### Cat 9 — Intervalli polling fissi
- `kTimeoutCheckInterval = [0-9]+` non configurabile

## Output format

```
HARDCODE AUDIT REPORT — CAS Engine
Date: {today}

CRITICAL (hot path: integrate/groebner/factorize):
  src/algebra/factorization_integers.cpp:232 | Cat3 | kCandidates[] fisso | Hash-based cyclic selection
  ...

HIGH:
  ...

MEDIUM:
  ...

Total: N findings (C critical, H high, M medium)
Already in HARDCODE_LEDGER.md: (list if found)
```

Mark as CRITICAL if pattern is in: `integrate_*.cpp`, `factorization_*.cpp`, `polynomial_groebner*.cpp`, `polynomial_gcd*.cpp`.
Mark as HIGH if in: `simplify_*.cpp`, `limit*.cpp`, `ode_*.cpp`.
Mark as MEDIUM otherwise.
