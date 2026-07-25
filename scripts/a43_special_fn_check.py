#!/usr/bin/env python3
"""Banco di verifica per Nonelementary_Antiderivatives.md (task A43).

Ogni formula della spec e' verificata qui PRIMA di entrare nel motore, secondo
il metodo stabilito in A45/A46: simulare e certificare in mpmath, non dedurre.

Uso:  python3 scripts/a43_special_fn_check.py
Exit: 0 se tutte le verifiche passano, 1 altrimenti (usabile in CI).
"""

import sys

from mpmath import (chi, ci, cos, diff, ei, erf, erfi, exp, j, li, log, mp,
                    mpf, pi, polylog, quad, shi, si, sin, sqrt)

mp.dps = 30

FAILURES = []


def check(label, lhs, rhs, rel=mpf("1e-20")):
    """Confronto a errore RELATIVO: erfi(x) cresce come e^{x^2}, un confronto
    assoluto fallirebbe per pura scala anche con formula esatta."""
    denom = abs(rhs) if abs(rhs) > 0 else mpf(1)
    err = abs(lhs - rhs) / denom
    ok = err < rel
    if not ok:
        FAILURES.append(f"{label}: rel_err={err}")
    print(f"  {'OK ' if ok else 'FAIL'}  {label}")
    return ok


print("== derivate (spec §3) ==")
for x in (mpf("0.7"), mpf("2.3"), mpf("5.1")):
    check(f"d/dx Ei({x})", diff(ei, x), exp(x) / x)
    check(f"d/dx Si({x})", diff(si, x), sin(x) / x)
    check(f"d/dx Ci({x})", diff(ci, x), cos(x) / x)
    check(f"d/dx erfi({x})", diff(erfi, x), 2 / sqrt(pi) * exp(x**2))
    check(f"d/dx Shi({x})", diff(shi, x), mp.sinh(x) / x)
    check(f"d/dx Chi({x})", diff(chi, x), mp.cosh(x) / x)
    if x > 1:
        check(f"d/dx li({x})", diff(li, x), 1 / log(x))
for z in (mpf("0.3"), mpf("-0.8"), mpf("0.55")):
    check(f"d/dz Li2({z})", diff(lambda t: polylog(2, t), z), -log(1 - z) / z)

print("== identita' di riduzione (spec §4) ==")
for x in (mpf("1.7"), mpf("3.2")):
    check(f"li({x}) = Ei(ln {x})", li(x), ei(log(x)))
for x in (mpf("0.4"), mpf("1.1")):
    check(f"erfi({x}) = -i erf(i {x})", erfi(x), (-j * erf(j * x)).real)
    check(f"Shi({x}) = (Ei(x)-Ei(-x))/2", shi(x), (ei(x) - ei(-x)) / 2)
    check(f"Chi({x}) = (Ei(x)+Ei(-x))/2", chi(x), (ei(x) + ei(-x)) / 2)

print("== primitive (spec §5) ==")
check("int_1^3 e^t/t = Ei(3)-Ei(1)",
      quad(lambda t: exp(t) / t, [1, 3]), ei(3) - ei(1))
check("int_1^3 sin(t)/t = Si(3)-Si(1)",
      quad(lambda t: sin(t) / t, [1, 3]), si(3) - si(1))
check("int_1^3 cos(t)/t = Ci(3)-Ci(1)",
      quad(lambda t: cos(t) / t, [1, 3]), ci(3) - ci(1))
check("int_2^4 dt/ln t = li(4)-li(2)",
      quad(lambda t: 1 / log(t), [2, 4]), li(4) - li(2))
check("int_0^0.6 e^{t^2} = sqrt(pi)/2 erfi(0.6)",
      quad(lambda t: exp(t**2), [0, mpf("0.6")]),
      sqrt(pi) / 2 * erfi(mpf("0.6")))
check("int_0^0.5 ln(1-t)/t = -Li2(0.5)",
      quad(lambda t: log(1 - t) / t, [0, mpf("0.5")]),
      -polylog(2, mpf("0.5")))

print("== primitive per parti (spec §5, righe derivate) ==")
# int e^x/x^2 = Ei(x) - e^x/x ; verificata come differenza su [1,3]
F = lambda t: ei(t) - exp(t) / t
check("int_1^3 e^t/t^2 = [Ei(t)-e^t/t]",
      quad(lambda t: exp(t) / t**2, [1, 3]), F(3) - F(1))
# int e^x ln x = e^x ln x - Ei(x)
G = lambda t: exp(t) * log(t) - ei(t)
check("int_1^3 e^t ln t = [e^t ln t - Ei(t)]",
      quad(lambda t: exp(t) * log(t), [1, 3]), G(3) - G(1))
# int e^-x ln x = -e^-x ln x + Ei(-x)
# NB: il segno di Ei(-x) e' +, non -. La prima stesura della spec aveva -, ed e'
# stato questo banco a smentirla: d/dx[-e^-x ln x] = e^-x ln x - e^-x/x, e
# int e^-x/x dx = Ei(-x), quindi il termine di compensazione si somma.
H = lambda t: -exp(-t) * log(t) + ei(-t)
check("int_1^3 e^-t ln t = [-e^-t ln t + Ei(-t)]",
      quad(lambda t: exp(-t) * log(t), [1, 3]), H(3) - H(1))

if FAILURES:
    print(f"\n{len(FAILURES)} VERIFICHE FALLITE:")
    for f in FAILURES:
        print("  " + f)
    sys.exit(1)
print("\nTutte le verifiche passano.")
