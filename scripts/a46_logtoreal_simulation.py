#!/usr/bin/env python3
"""A46 — simulazione fedele di LogToReal / LogToAtan (Bronstein §2.8) in sympy.

Serve a fissare l'algoritmo PRIMA di scriverlo in C++ (lezione A45: simulare
riga-per-riga invece di debuggare da integrate()). Riproduce:

  1. i dati che `integrate_rational_lrt` passa a `rioboo_conversion`
     (R_ij = fattore irriducibile della parte squarefree del resultant di
     Rothstein-Trager, S_i = subresultante di grado i);
  2. LogToAtan(A, B)          -> spec righe 2188-2200;
  3. LogToReal(R, S)          -> spec righe 2381-2395;
  4. certificato d/dx F == f (simbolico + numerico multi-punto).

Uso:  python3 scripts/a46_logtoreal_simulation.py [--verbose]
"""

from __future__ import annotations

import argparse
import sys

import sympy as sp

x = sp.Symbol("x", real=True)
t = sp.Symbol("t")


# --------------------------------------------------------------------------
# 1. dati LRT (quello che il C++ ha già in mano quando chiama rioboo_conversion)
# --------------------------------------------------------------------------
def lrt_data(num, den):
    """Ritorna [(R_ij, S_i)] con R_ij irriducibile in Q[t], S_i in Q[t][x]."""
    P = sp.Poly(num, x)
    Q = sp.Poly(den, x)
    target = sp.Poly(P.as_expr() - t * Q.diff(x).as_expr(), x)
    chain = sp.subresultants(Q, target, x)
    R = sp.Poly(sp.resultant(Q, target, x), t)

    out = []
    _, sqf = sp.sqf_list(R.as_expr(), t)
    for factor, mult in sqf:
        fac_poly = sp.Poly(factor, t)
        if fac_poly.degree() == 0:
            continue
        i = mult
        if i == Q.degree():
            S_i = Q
        else:
            S_i = None
            for step in chain:
                sp_step = sp.Poly(step, x)
                if sp_step.degree() == i:
                    S_i = sp_step
                    break
        if S_i is None:
            continue
        # presentazione per fattori razionali (spec riga 2378): stessa somma
        # formale, ma i fattori irriducibili tengono le forme chiuse separate
        for irr, _m in sp.factor_list(factor, t)[1]:
            irr_poly = sp.Poly(irr, t)
            if irr_poly.degree() > 0:
                out.append((irr_poly, S_i))
    return out


# --------------------------------------------------------------------------
# 2. LogToAtan (spec righe 2188-2200)
# --------------------------------------------------------------------------
def as_alg_poly(expr) -> sp.Poly:
    """Poly in x sull'estensione algebrica generata dai coefficienti.

    `extension=True` è indispensabile: senza, sympy prende sqrt(2) come
    *generatore* e ogni gcdex fallisce con "univariate polynomial expected"
    (l'equivalente C++ è: i coefficienti vivono in K = Q(a,b), non in Q).
    """
    return sp.Poly(sp.expand(expr), x, extension=True)


def log_to_atan(A: sp.Poly, B: sp.Poly, depth: int = 0):
    """f tale che df/dx = d/dx [ i·log((A+iB)/(A-iB)) ]."""
    if depth > 64:
        raise RuntimeError("LogToAtan: ricorsione oltre il bound")
    if B.is_zero:
        raise ValueError("LogToAtan richiede B != 0")

    quot, rem = sp.div(A, B)
    if rem.is_zero:                                     # if B | A
        return 2 * sp.atan(sp.simplify(quot.as_expr()))
    if A.degree() < B.degree():                         # LogToAtan(-B, A)
        return log_to_atan(-B, A, depth + 1)

    # (D, C, G) <- ExtendedEuclidean(B, -A):  B·D - A·C = G
    D, C, G = sp.gcdex(B, -A)
    arg = sp.cancel(((A * D + B * C).as_expr()) / G.as_expr())
    return 2 * sp.atan(sp.radsimp(arg)) + log_to_atan(D, C, depth + 1)


# --------------------------------------------------------------------------
# 3. LogToReal (spec righe 2381-2395), variante "radici calcolate"
#    (sanzionata dalla spec: «whenever the real roots ... can be computed
#     efficiently ... it can be more efficient to first compute the roots,
#     and then call LogToAtan»)
# --------------------------------------------------------------------------
def log_to_real(R: sp.Poly, S: sp.Poly, verbose=False):
    terms = []
    roots = sp.roots(R, multiple=False)
    if not roots or sum(roots.values()) != R.degree():
        raise RuntimeError(f"radici non esprimibili per R = {R.as_expr()}")

    for root, mult in roots.items():
        root = sp.radsimp(sp.simplify(sp.expand(root)))
        a, b = sp.simplify(sp.re(root)), sp.simplify(sp.im(root))
        if b == 0:                                       # somma sulle radici reali
            S_a = as_alg_poly(S.as_expr().subs(t, root))
            terms.append(mult * a * sp.log(S_a.as_expr()))
            continue
        if sp.simplify(b).is_negative:                   # solo b > 0: la coppia
            continue                                     # coniugata è già coperta
        S_ab = sp.expand(S.as_expr().subs(t, a + sp.I * b))
        A = as_alg_poly(sp.re(S_ab))
        B = as_alg_poly(sp.im(S_ab))
        if verbose:
            print(f"    root a={a}  b={b}\n      A={A.as_expr()}\n      B={B.as_expr()}")
        norm = sp.expand(A.as_expr() ** 2 + B.as_expr() ** 2)
        terms.append(mult * (a * sp.log(norm) + b * log_to_atan(A, B)))
    return sp.Add(*terms)


# --------------------------------------------------------------------------
# 4. certificato
# --------------------------------------------------------------------------
def certify(F, f, label, verbose=False):
    d = sp.simplify(sp.diff(F, x) - f)
    ok_sym = d == 0
    pts = [sp.Rational(1, 3), sp.Rational(3, 2), sp.Rational(-7, 5), sp.Rational(11, 4), sp.Rational(-1, 2)]
    worst = 0
    for p in pts:
        val = sp.N(sp.diff(F, x).subs(x, p) - f.subs(x, p), 30)
        worst = max(worst, abs(complex(val)))
    ok_num = worst < 1e-25
    status = "OK " if (ok_sym or ok_num) else "FAIL"
    print(f"  [{status}] {label}: D(F)=f  simbolico={ok_sym}  num_max_err={sp.N(worst,3)}")
    if verbose or not (ok_sym or ok_num):
        print(f"        F = {F}")
    return ok_sym or ok_num


CASES = [
    ("1/(x**2+1)", 1, x**2 + 1),
    ("1/(x**4+1)", 1, x**4 + 1),
    ("x**2/(x**4+1)", x**2, x**4 + 1),
    ("x**3/(x**4+1)", x**3, x**4 + 1),
    ("1/(x**6+1)", 1, x**6 + 1),
    ("1/(x**5+x+1)", 1, x**5 + x + 1),
    ("x/(x**5+x+1)", x, x**5 + x + 1),
    ("1/(x**8+1)", 1, x**8 + 1),
    ("(x**2-1)/(x**4+1)", x**2 - 1, x**4 + 1),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--only", default=None)
    args = ap.parse_args()

    failures = []
    for label, num, den in CASES:
        if args.only and args.only not in label:
            continue
        f = num / den
        print(f"\n=== {label}")
        try:
            data = lrt_data(num, den)
            parts = []
            for R_ij, S_i in data:
                print(f"  R={sp.factor(R_ij.as_expr())}  deg={R_ij.degree()}  "
                      f"S deg_x={S_i.degree()}")
                parts.append(log_to_real(R_ij, S_i, verbose=args.verbose))
            F = sp.Add(*parts)
            if not certify(F, f, label, verbose=args.verbose):
                failures.append(label)
        except Exception as exc:  # noqa: BLE001 — la sim deve dire *perché* cede
            print(f"  [SKIP] {label}: {type(exc).__name__}: {exc}")
            failures.append(f"{label} ({type(exc).__name__})")

    print("\n--- riepilogo ---")
    print("falliti/non coperti:", failures if failures else "nessuno")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())


# --------------------------------------------------------------------------
# 5. variante spec "sistema P=Q=0" (righe 2381-2395 alla lettera): misura quali
#    gradi vanno effettivamente risolti, per sapere che solver serve in C++.
# --------------------------------------------------------------------------
def system_route_report(R: sp.Poly, S: sp.Poly):
    u, v = sp.symbols("u v", real=True)
    Ruv = sp.expand(R.as_expr().subs(t, u + sp.I * v))
    P = sp.Poly(sp.expand(sp.re(Ruv)), u, v)
    Q = sp.Poly(sp.expand(sp.im(Ruv)), u, v)
    H = sp.Poly(sp.resultant(P, Q, v), u)
    print(f"    P = {P.as_expr()}")
    print(f"    Q = {Q.as_expr()}")
    print(f"    H = res_v(P,Q) = {sp.factor(H.as_expr())}")
    for fac, mult in sp.factor_list(H.as_expr(), u)[1]:
        fp = sp.Poly(fac, u)
        if fp.degree() == 0:
            continue
        print(f"      fattore di H: {fac}   deg_u={fp.degree()}  mult={mult}")
        for a in sp.roots(fp, multiple=True):
            if not sp.im(a) == 0:
                continue
            Pa = sp.Poly(sp.expand(P.as_expr().subs(u, a)), v)
            Qa = sp.Poly(sp.expand(Q.as_expr().subs(u, a)), v)
            g = sp.gcd(Pa, Qa)
            print(f"        a={a}  gcd_v(P(a,v),Q(a,v)) = {g.as_expr()}  deg_v={g.degree()}")
            for fac2, _m2 in sp.factor_list(g.as_expr(), v, extension=True)[1]:
                print(f"          fattore in v: {fac2}  deg={sp.Poly(fac2, v).degree()}")


# --------------------------------------------------------------------------
# 6. VIA F — quartico irriducibile su Q -> due quadratici reali su Q(sqrt(D)),
#    poi il ramo grado-2 di rioboo_conversion (gia' nel motore) su ciascuno.
#
#    Tutte le decisioni (radice razionale del cubico risolvente, realta' dei
#    due discriminanti) avvengono sui RAZIONALI: nessun test di segno su
#    espressioni radicali, nessun solver di sistemi bivariati.
# --------------------------------------------------------------------------
def rioboo_degree_two(a, b, G: sp.Poly, z):
    """Trascrizione fedele del ramo grado-2 del C++ (partial_fractions_lrt.cpp
    righe 167-234) con a, b simbolici invece che razionali."""
    monic_r = sp.Poly(z**2 + a * z + b, z)
    Q1_coeffs, Q0_coeffs = [], []
    for i in range(G.degree() + 1):
        g_i = sp.Poly(G.as_expr().coeff(x, i), z)
        rem = sp.rem(g_i, monic_r)
        rc = sp.Poly(rem, z).all_coeffs()[::-1] + [0, 0]
        Q0_coeffs.append(sp.simplify(rc[0]))
        Q1_coeffs.append(sp.simplify(rc[1]))
    Q0 = sum(c * x**i for i, c in enumerate(Q0_coeffs))
    Q1 = sum(c * x**i for i, c in enumerate(Q1_coeffs))
    disc = sp.simplify(4 * b - a**2)
    sqrt_disc = sp.sqrt(disc)
    norm = sp.simplify(sp.expand(Q0**2 - a * Q1 * Q0 + b * Q1**2))
    ln_term = -a / 2 * sp.log(norm)
    atan_term = -sqrt_disc * sp.atan(sp.simplify(sqrt_disc * Q1 / (2 * Q0 - a * Q1)))
    return sp.simplify(ln_term + atan_term)


def quartic_to_real_quadratics(R: sp.Poly, z):
    """R quartico -> [(a1,b1),(a2,b2)] con z^2+a_i z+b_i su Q(sqrt(D)).

    Cubico risolvente: y^3 - q y^2 + (p r - 4 s) y - (p^2 s - 4 q s + r^2).
    Serve una radice RAZIONALE y0 (esiste esattamente quando il quartico
    fattorizza in due quadratici su un'estensione quadratica), poi
      e,g = (p +- sqrt(p^2-4q+4y0))/2      f,h = (y0 +- sqrt(y0^2-4s))/2
    con l'accoppiamento verificato per espansione.
    """
    c = R.all_coeffs()
    lc = c[0]
    p, q, r, s = [sp.nsimplify(v / lc) for v in c[1:]]
    y = sp.Symbol("y")
    cubic = sp.Poly(y**3 - q * y**2 + (p * r - 4 * s) * y - (p**2 * s - 4 * q * s + r**2), y)
    rational_roots = [rt for rt in sp.roots(cubic, multiple=True) if rt.is_rational]
    out = []
    for y0 in rational_roots:
        d1 = sp.nsimplify(p**2 - 4 * q + 4 * y0)
        d2 = sp.nsimplify(y0**2 - 4 * s)
        if d1 < 0 or d2 < 0:
            continue
        e = (p + sp.sqrt(d1)) / 2
        g = (p - sp.sqrt(d1)) / 2
        for f, h in [((y0 + sp.sqrt(d2)) / 2, (y0 - sp.sqrt(d2)) / 2),
                     ((y0 - sp.sqrt(d2)) / 2, (y0 + sp.sqrt(d2)) / 2)]:
            prod = sp.expand((z**2 + e * z + f) * (z**2 + g * z + h))
            if sp.simplify(prod - sp.expand(R.as_expr() / lc)) == 0:
                out = [(sp.radsimp(e), sp.radsimp(f)), (sp.radsimp(g), sp.radsimp(h))]
                return out, y0
    return None, None


def route_f(num, den, verbose=False):
    """rioboo_conversion via VIA F su tutti i fattori del resultant."""
    z = t
    parts = []
    for R_ij, S_i in lrt_data(num, den):
        G = sp.Poly(S_i.as_expr(), x)
        deg = R_ij.degree()
        if deg == 1:
            root = sp.roots(R_ij, multiple=True)[0]
            parts.append(root * sp.log(sp.expand(S_i.as_expr().subs(z, root))))
        elif deg == 2:
            cc = R_ij.all_coeffs()
            parts.append(rioboo_degree_two(cc[1] / cc[0], cc[2] / cc[0], G, z))
        elif deg == 4:
            if sp.real_roots(R_ij):
                raise RuntimeError("quartico con radici reali: fuori dalla classe supportata")
            quads, y0 = quartic_to_real_quadratics(R_ij, z)
            if quads is None:
                raise RuntimeError("nessuna radice razionale utilizzabile del cubico risolvente")
            if verbose:
                print(f"      y0={y0}  quadratici: {quads}")
            for a, b in quads:
                parts.append(rioboo_degree_two(a, b, G, z))
        else:
            raise RuntimeError(f"grado {deg} fuori dalla classe supportata (radicali non garantiti)")
    return sp.Add(*parts)
