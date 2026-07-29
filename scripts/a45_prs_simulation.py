#!/usr/bin/env python3
"""A45 — banco di prova della catena subresultante (LRT), fuori dal motore.

Perche' esiste: i due tentativi falliti su A45 hanno diagnosticato il bug
*attraverso* integrate(), con troppe variabili in gioco. Questo script simula
FEDELMENTE, riga per riga, le due implementazioni C++ reali con coefficienti
sympy esatti, e le confronta con l'oracolo (sympy subresultants/resultant, a
sua volta verificato contro Maxima). Risponde a "quale ricorrenza sbaglia, e
su quali sequenze" in secondi, senza compilare il motore.

  parte 1: resultant_generic_sim()  <- src/algebra/polynomial_resultant_generic.hpp:229-372
  parte 2: subresultant_prs_head()  <- src/algebra/partial_fractions_lrt.cpp:97-155

Misura 2026-07-24 (vedi A45 in TASKLIST_MASTER.md):
  - resultant_generic: corretto su tutto il corpus, sequenze difettose incluse
    (gap 2 e 3) — il tracking di `c` (gamma) applica tau_k implicitamente.
  - partial_fractions_lrt: sbaglia esattamente e solo sulle sequenze difettose
    — manca il fattore tau_k (Bronstein, Symbolic_Integration_I.md righe 1035
    e 1066-1070: resultant = R_k solo se deg(R_{k-1}) = 1, altrimenti tau_k*R_k).

Uso:  python3 scripts/a45_prs_simulation.py     (exit 1 se il riferimento regredisce)
"""

import sympy as sp
from sympy import Poly, resultant, subresultants, symbols

z, x = symbols("z x")


# --------------------------------------------------------------------------
# Helper comuni. Un polinomio e' una lista di coefficienti a grado crescente,
# [] = polinomio nullo — stessa convenzione di PolyExpr / std::vector<Coeff>.
# --------------------------------------------------------------------------
def deg(p):
    return 0 if not p else len(p) - 1


def strip(p):
    while p and sp.expand(p[-1]) == 0:
        p.pop()
    return p


def is_zero(p):
    return all(sp.expand(c) == 0 for c in p)


def to_list(expr, var):
    p = Poly(expr, var)
    return [sp.expand(p.coeff_monomial(var**i)) for i in range(p.degree() + 1)]


# --------------------------------------------------------------------------
# PARTE 1 — polynomial_resultant_generic.hpp (il riferimento CORRETTO)
# --------------------------------------------------------------------------
def prem_generic(a, b):
    """pseudo_remainder_generic, righe 169-227."""
    a = strip(list(a))
    divisor = strip(list(b))
    if not divisor:
        raise ValueError("zero divisor")
    m, n = deg(a), deg(divisor)
    if not a or m < n:
        return a
    lead_b = divisor[-1]
    r = list(a)
    for step in range(0, m - n + 1):
        r = strip(r)
        if not r:
            break
        deg_r = deg(r)
        if deg_r == m - step:
            lead_r = r[-1]
            r = strip([c * lead_b for c in r])
            shift = deg_r - n
            if len(r) < len(divisor) + shift:
                r += [sp.Integer(0)] * (len(divisor) + shift - len(r))
            for i in range(len(divisor)):
                r[shift + i] = sp.expand(r[shift + i] - divisor[i] * lead_r)
        else:
            r = strip([c * lead_b for c in r])
    return strip(r)


def resultant_generic_sim(a, b):
    """resultant_generic, righe 229-372 — traduzione 1:1."""
    a, b = strip(list(a)), strip(list(b))
    if not a or not b:
        return sp.Integer(0)
    n, m = deg(a), deg(b)
    sign_correction = 1
    if n < m:
        a, b = b, a
        n, m = m, n
        if (n * m) % 2 != 0:
            sign_correction = -1
    minus_one = sp.Integer(-1)
    b_scale = minus_one if ((n - m + 1) % 2 != 0) else sp.Integer(1)
    h = strip([sp.expand(c * b_scale) for c in prem_generic(a, b)])
    lc = b[-1]
    s_last = sp.expand(minus_one * lc ** (n - m))
    c = s_last
    guard = 0
    while not is_zero(h):
        guard += 1
        if guard > 60:
            raise RuntimeError("chain did not terminate")
        k = deg(h)
        a, b = b, h
        d_next = deg(a) - k
        b_scale = sp.expand(minus_one * lc * c**d_next)
        h = strip([sp.cancel(sp.expand(cc) / b_scale) for cc in prem_generic(a, b)])
        # NB: lc va riletto QUI, dopo lo swap. Lo stash A45 stall usava il
        # valore precedente (-lc(R_{i-1}) invece di -lc(R_i)): e' quel bug.
        lc = b[-1]
        if d_next > 1:
            c = sp.cancel(sp.expand((minus_one * lc) ** d_next) / sp.expand(c ** (d_next - 1)))
        else:
            c = sp.expand(minus_one * lc)
        s_last = sp.expand(minus_one * c)
    if deg(b) > 0 and not is_zero(b):
        return sp.Integer(0)
    return sp.expand(minus_one * s_last) if sign_correction == -1 else sp.expand(s_last)


# --------------------------------------------------------------------------
# PARTE 2 — partial_fractions_lrt.cpp (la versione con il BUG)
# --------------------------------------------------------------------------
def prem_lrt(A, B):
    """pseudo_remainder_lrt, righe 52-90."""
    B = strip(list(B))
    if not B:
        raise ValueError("zero divisor")
    R = strip(list(A))
    m, n = deg(R), deg(B)
    if m < n:
        return R
    b_n = B[-1]
    for step in range(0, m - n + 1):
        if is_zero(R):
            break
        deg_r = deg(R)
        lc_r_before = R[-1]
        R = [sp.expand(c * b_n) for c in R]
        if deg_r == m - step:
            shift = deg_r - n
            prod = [sp.Integer(0)] * (len(B) + shift)
            for i in range(len(B)):
                prod[i + shift] = sp.expand(B[i] * lc_r_before)
            sz = max(len(R), len(prod))
            nxt = [sp.Integer(0)] * sz
            for i in range(len(R)):
                nxt[i] = R[i]
            for i in range(len(prod)):
                nxt[i] = sp.expand(nxt[i] - prod[i])
            R = strip(nxt)
    return strip(R)


def subresultant_prs_head(P1, P2):
    """subresultant_prs_lrt, righe 97-155 — ricorrenza g/h, SENZA tau_k.

    Ritorna la catena come lista di (coefficienti, grado), nello stesso ordine
    in cui integrate_rational_lrt la consuma (righe 324-361).
    """
    P1, P2 = strip(list(P1)), strip(list(P2))
    steps = [(list(P1), deg(P1)), (list(P2), deg(P2))]
    g = sp.Integer(1)
    h = sp.Integer(1)
    guard = 0
    while True:
        guard += 1
        if guard > 40:
            raise RuntimeError("chain did not terminate")
        d1, d2 = deg(P1), deg(P2)
        if d1 < d2:
            P1, P2, d1, d2 = P2, P1, d2, d1
        d = d1 - d2
        P3 = prem_lrt(P1, P2)
        if is_zero(P3):
            break
        divisor = sp.expand(g * h**d)
        P3 = strip([sp.cancel(sp.expand(c) / divisor) for c in P3])
        P1, P2 = P2, P3
        g = P1[-1]
        h = sp.cancel(sp.expand(g**d) / sp.expand(h ** (d - 1))) if d > 1 else g
        steps.append((list(P2), deg(P2)))
    return steps


# --------------------------------------------------------------------------
# Corpus: i casi A45 reali (tutti nel golden, integrate/basic.jsonl:117-120)
# piu' altri costruiti per coprire sequenze regolari e difettose.
# --------------------------------------------------------------------------
def build_corpus():
    out = []
    for P in [x**2, x**2 + 1, x**2 - 1, x**3, x**3 + x, sp.Integer(1), x, x**2 + x + 1]:
        out.append((x**4 + 1, P))
    for P in [x**2, x**2 + 1, x, x**3 + 1]:
        out.append((x**4 + x**2 + 1, P))
    for P in [x, x**2, x**4, x**3 - 2]:
        out.append((x**5 + x + 1, P))
    for P in [sp.Integer(1), x**2]:
        out.append((x**3 - x, P))
    out.append((x**2 + 1, sp.Integer(1)))
    out.append((x**2 - 1, sp.Integer(1)))
    # Bronstein, esempio 2.5.1 (Symbolic_Integration_I.md riga 1841)
    out.append((x**6 - 5 * x**4 + 5 * x**2 + 4, x**4 - 3 * x**2 + 6))
    return out


def main():
    print("A45 — catene subresultanti del motore vs oracolo sympy (verificato su Maxima)\n")
    hdr = f"{'integrale':<34} {'gradi PRS':<22} {'difett':<7} {'generic':<8} {'lrt(HEAD)':<9}"
    print(hdr)
    print("-" * len(hdr))

    corpus = build_corpus()
    bad_generic, bad_head, missing_sub = [], [], []
    for Q, P in corpus:
        label = f"({P})/({Q})"
        T = sp.expand(P - z * sp.diff(Q, x))
        A, B = to_list(Q, x), to_list(T, x)
        true_res = sp.expand(resultant(Q, T, x))
        oracle_chain = subresultants(Poly(Q, x), Poly(T, x), x)
        degs = [Poly(c, x).degree() if c != 0 else -1 for c in oracle_chain]
        defective = any(degs[i] - degs[i + 1] > 1 for i in range(1, len(degs) - 1))

        try:
            ok_gen = sp.simplify(sp.expand(resultant_generic_sim(A, B) - true_res)) == 0
        except Exception as exc:  # noqa: BLE001 — diagnostico, non controllo di flusso
            ok_gen = f"EXC({exc})"

        try:
            steps = subresultant_prs_head(A, B)
            head_degs = [d for _, d in steps]
            r_head = next((p[0] if p else sp.Integer(0) for p, d in reversed(steps) if d == 0), None)
            ok_head = r_head is not None and sp.simplify(sp.expand(r_head - true_res)) == 0
            # Bronstein riga 1828: serve R_m con deg_x(R_m) = i, dove i e' la
            # molteplicita' nella fattorizzazione squarefree del resultant.
            mults = {m for f, m in sp.factor_list(Poly(true_res, z))[1] if Poly(f, z).degree() > 0}
            deg_Q = Poly(Q, x).degree()
            if not all(m == deg_Q or m in head_degs[1:] for m in mults):
                missing_sub.append(label)
        except Exception as exc:  # noqa: BLE001
            ok_head = f"EXC({exc})"

        if ok_gen is not True:
            bad_generic.append(label)
        if ok_head is not True:
            bad_head.append(label)
        print(f"{label:<34} {str(degs):<22} {str(defective):<7} {str(ok_gen):<8} {str(ok_head):<9}")

    total = len(corpus)
    print(f"\nresultant_generic (riferimento) : {total - len(bad_generic)}/{total} corretti")
    print(f"partial_fractions_lrt (HEAD)    : {total - len(bad_head)}/{total} corretti")
    print(f"subresultante di grado i assente: {len(missing_sub)}/{total}")
    if bad_head:
        print("\nSBAGLIATI in partial_fractions_lrt (attesi: solo le sequenze difettose):")
        for label in bad_head:
            print(f"  - {label}")
    if bad_generic:
        print("\nATTENZIONE: il riferimento resultant_generic ha REGREDITO:")
        for label in bad_generic:
            print(f"  - {label}")
    return 1 if bad_generic else 0


if __name__ == "__main__":
    raise SystemExit(main())
