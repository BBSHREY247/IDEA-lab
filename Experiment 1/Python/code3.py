"""
Laplace Transform Solver - with step-by-step working
=====================================================
Type a function of t and get L{f(t)} as a function of s, with every step shown.
Also does inverse transforms (start the line with "inv").

Install:  pip install sympy
Run:      python laplace_solver.py                (interactive)
          python laplace_solver.py "t^2 sin3t"    (single problem)

Input examples:
    t^2 sin3t            t2sin3t            L{t^2 sin3t}
    e^-2t cos4t          t e^-3t            t^3 + 4t - 5
    sin t cos t          sin^2 3t           (t+1)^2 e^t
    cosh 2t              5                  3e^t
    inv 1/(s^2+1)        inv (3s+2)/(s^2-4)
"""
import re
import sys

import sympy as sp
from sympy.simplify.fu import TR8
from sympy.parsing.sympy_parser import (
    parse_expr, standard_transformations,
    implicit_multiplication_application, convert_xor,
)

t, s = sp.symbols("t s", positive=True)
TRANSFORMS = standard_transformations + (implicit_multiplication_application, convert_xor)
FUNCS = "sinh|cosh|tanh|sin|cos|tan|exp|ln|log"
ARG = r"-?\d*\.?\d*\*?t"            # 3t, -2t, t, 0.5*t
TRIG = (sp.sin, sp.cos, sp.sinh, sp.cosh)


# --------------------------------------------------------------------------
# Reading the input
# --------------------------------------------------------------------------
def clean(text: str) -> str:
    """Turn textbook-style input (t2sin3t, e^-2t, sin3t ...) into valid sympy text."""
    text = text.strip().replace("−", "-").replace("^", "**")
    text = text.replace("×", "*").replace("·", "*")
    # t2 -> t**2, s2 -> s**2   (power typed without ^)
    text = re.sub(r"(?<![A-Za-z])([ts])(\d+)(?![\d.])", r"\1**\2", text)
    # separate glued function names:  t**2sin3t -> t**2 sin3t
    text = re.sub(rf"(?<=[0-9t)])({FUNCS})", r" \1", text)
    # sin^2 3t -> (sin(3t))**2
    text = re.sub(rf"(sinh|cosh|sin|cos|tan)\s*\*\*\s*(\d+)\s*({ARG})", r"(\1(\3))**\2", text)
    # e^-2t -> exp(-2t)
    text = re.sub(rf"(?<![A-Za-z])e\s*\*\*\s*({ARG})", r"exp(\1)", text)
    # sin3t -> sin(3t)
    text = re.sub(rf"({FUNCS})\s*({ARG})", r"\1(\2)", text)
    # t**2(...) / t**2sin(...) -> explicit multiplication
    text = re.sub(r"(\*\*\d+)(?=[A-Za-z(])", r"\1*", text)
    return text


def parse(text: str):
    return parse_expr(clean(text), local_dict={"t": t, "s": s, "e": sp.E},
                      transformations=TRANSFORMS)


# --------------------------------------------------------------------------
# Printing helpers (ASCII only, so it works in any console)
# --------------------------------------------------------------------------
def fmt(e) -> str:
    return str(e).replace("**", "^")


def Lf(e) -> str:
    return "L{" + fmt(e) + "}"


def tidy(e):
    return sp.factor(sp.simplify(e))


class Log:
    def __init__(self):
        self.n = 0

    def step(self, title):
        self.n += 1
        print(f"\nStep {self.n}: {title}")

    def line(self, text, indent=1):
        print("    " * indent + text)


def lin_coeff(arg):
    """Return k if arg == k*t (k free of t), else None."""
    k = sp.simplify(arg / t)
    return None if k.has(t) else k


# --------------------------------------------------------------------------
# Forward transform
# --------------------------------------------------------------------------
def fallback(term, log, why):
    log.step(f"{why} - use the definition  L{{f}} = integral from 0 to infinity of e^(-st) f(t) dt")
    F = sp.laplace_transform(term, t, s, noconds=True)
    if F.has(sp.LaplaceTransform):
        raise ValueError(f"no rule found for {fmt(term)}")
    F = tidy(F)
    log.line(f"{Lf(term)} = {fmt(F)}")
    return F


def solve_term(term, log):
    """Transform ONE term (no top-level sum), showing the steps."""
    c, rest = term.as_independent(t, as_Add=False)

    # ---- constant ---------------------------------------------------------
    if rest == 1:
        log.step("Constant:  L{k} = k/s")
        F = c / s
        log.line(f"{Lf(term)} = {fmt(F)}")
        return F

    # ---- split the term into  c * t^n * e^(a t) * trig(b t)  ---------------
    n, a, trig = 0, sp.Integer(0), []
    need_rewrite, ok = False, True
    for fac in sp.Mul.make_args(rest):
        if fac == t:
            n += 1
        elif fac.is_Pow and fac.base == t and fac.exp.is_Integer and fac.exp > 0:
            n += int(fac.exp)
        elif fac.func == sp.exp:
            k = lin_coeff(fac.args[0])
            if k is None:
                ok = False
            else:
                a += k
        elif fac.func in TRIG:
            if lin_coeff(fac.args[0]) is None:
                ok = False
            trig.append(fac)
        elif fac.is_Pow and fac.base.func in TRIG and fac.exp.is_Integer and fac.exp > 1:
            need_rewrite = True
        else:
            ok = False

    # ---- products / powers of sin, cos: product-to-sum first --------------
    if ok and (len(trig) > 1 or need_rewrite):
        new = sp.expand(TR8(rest))
        if new != sp.expand(rest) and not new.has(sp.sinh, sp.cosh):
            log.step("Product/power of trig functions - rewrite as a sum (product-to-sum identities)")
            log.line(f"{fmt(rest)} = {fmt(new)}")
            F = transform(new, log)
            if c != 1:
                log.step(f"Constant multiple: multiply by {fmt(c)}")
                F = tidy(c * F)
                log.line(f"{Lf(term)} = {fmt(F)}")
            return F
        ok = False

    if not ok:
        return fallback(term, log, "Pattern not covered by the standard rules")

    # ---- describe the pieces ---------------------------------------------
    log.step(f"Identify the pieces of  {fmt(term)}")
    parts = []
    if c != 1:
        parts.append(f"constant factor {fmt(c)}")
    if n:
        parts.append("t" if n == 1 else f"t^{n}")
    if a != 0:
        parts.append(f"exponential e^({fmt(a * t)})")
    if trig:
        parts.append(fmt(trig[0]))
    log.line("Built from: " + ",  ".join(parts))

    # ---- basic table formula ---------------------------------------------
    n_diff = 0
    log.step("Start from the basic table formula")
    if trig:
        fn, b = trig[0].func, lin_coeff(trig[0].args[0])
        rules = {
            sp.sin: ("L{sin(bt)} = b/(s^2 + b^2)", b / (s**2 + b**2)),
            sp.cos: ("L{cos(bt)} = s/(s^2 + b^2)", s / (s**2 + b**2)),
            sp.sinh: ("L{sinh(bt)} = b/(s^2 - b^2)", b / (s**2 - b**2)),
            sp.cosh: ("L{cosh(bt)} = s/(s^2 - b^2)", s / (s**2 - b**2)),
        }
        rule, G = rules[fn]
        log.line(rule)
        log.line(f"with b = {fmt(b)}:  {Lf(trig[0])} = {fmt(G)}")
        n_diff = n
    elif n:
        G = sp.factorial(n) / s ** (n + 1)
        log.line("L{t^n} = n!/s^(n+1)")
        log.line(f"with n = {n}:  {Lf(t**n)} = {fmt(G)}")
    else:
        G = 1 / s
        log.line("L{1} = 1/s")
    F = G

    # ---- first shifting theorem ------------------------------------------
    if a != 0:
        log.step("First shifting theorem:  L{e^(at) f(t)} = F(s - a)")
        log.line(f"here a = {fmt(a)}, so replace s by  {fmt(s - a)}")
        shifted = G.subs(s, s - a)
        log.line(f"F(s) = {fmt(shifted)}")
        F = tidy(shifted)
        if sp.srepr(F) != sp.srepr(shifted):
            log.line(f"     = {fmt(F)}")

    # ---- multiplication by t^n -------------------------------------------
    if n_diff:
        log.step("Multiplication by t^n:  L{t^n f(t)} = (-1)^n * d^n/ds^n F(s)")
        log.line(f"here n = {n_diff},  F(s) = {fmt(F)}")
        D = F
        for k in range(1, n_diff + 1):
            D = tidy(sp.diff(D, s))
            label = "dF/ds" if k == 1 else f"d^{k}F/ds^{k}"
            log.line(f"{label} = {fmt(D)}")
        sign = (-1) ** n_diff
        F = tidy(sign * D)
        log.line(f"multiply by (-1)^{n_diff} = {sign}:   {fmt(F)}")

    # ---- constant ---------------------------------------------------------
    if c != 1:
        log.step(f"Constant multiple: multiply by {fmt(c)}")
        F = tidy(c * F)
        log.line(f"{Lf(term)} = {fmt(F)}")
    return F


def transform(expr, log):
    """Transform a full expression: split sums (linearity), then each term."""
    expr = sp.expand(expr)
    terms = sp.Add.make_args(expr)
    if len(terms) == 1:
        return solve_term(terms[0], log)

    log.step("Linearity:  L{f + g} = L{f} + L{g}  and  L{k f} = k L{f}")
    for i, tm in enumerate(terms, 1):
        log.line(f"term {i}: {fmt(tm)}")
    results = []
    for i, tm in enumerate(terms, 1):
        print(f"\n  --- Term {i}: {Lf(tm)} ---")
        Fi = solve_term(tm, log)
        results.append(Fi)
        log.line(f">> {Lf(tm)} = {fmt(Fi)}")
    log.step("Add the results of all terms")
    raw = sp.Add(*results, evaluate=False)
    log.line(f"= {fmt(raw)}")
    total = tidy(sum(results))
    log.line(f"= {fmt(total)}   (combined over a common denominator)")
    return total


def verify(expr, result):
    """True / False / None (could not check) against sympy's own transform."""
    try:
        ref = sp.laplace_transform(expr, t, s, noconds=True)
        if ref.has(sp.LaplaceTransform):
            return None
        diff = sp.simplify(sp.together(ref - result))
        if diff == 0:
            return True
        if diff.free_symbols <= {s}:
            return abs(complex(diff.subs(s, 3.7).evalf())) < 1e-9
        return False
    except Exception:
        return None


# --------------------------------------------------------------------------
# Inverse transform (lighter steps: partial fractions, then term by term)
# --------------------------------------------------------------------------
def inverse_steps(F, log):
    log.step("Break F(s) into simple (partial) fractions")
    pf = sp.apart(F, s)
    log.line(f"F(s) = {fmt(pf)}")
    terms = sp.Add.make_args(pf)

    log.step("Invert each term with the table  (L^-1 is linear)")
    total = 0
    for tm in terms:
        ft = sp.inverse_laplace_transform(tm, s, t)
        ft = sp.simplify(ft.subs(sp.Heaviside(t), 1))
        log.line(f"L^-1{{{fmt(tm)}}} = {fmt(ft)}")
        total += ft

    log.step("Add the results")
    total = sp.expand(sp.simplify(total))
    log.line(f"f(t) = {fmt(total)}")
    return total


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def banner(text):
    bar = "=" * max(50, len(text) + 4)
    print(f"\n{bar}\n  {text}\n{bar}")


def final(label, result):
    print(f"\nFINAL ANSWER\n    {label} = {fmt(result)}")
    print(f"    LaTeX: {sp.latex(result)}")


def solve(line: str):
    line = line.strip()
    log = Log()
    try:
        inv = re.match(r"^(?:inv(?:erse)?|l\s*\^?\s*-1)\s*", line, re.I)
        if inv:
            F = parse(line[inv.end():].strip().strip("{}").strip())
            banner(f"Find  L^-1{{ {fmt(F)} }}")
            res = inverse_steps(F, log)
            final(f"L^-1{{{fmt(F)}}}", res)
        else:
            m = re.match(r"^l\s*\{(.*)\}\s*$", line, re.I)
            f = parse(m.group(1) if m else line)
            banner(f"Find  {Lf(f)}")
            res = transform(f, log)
            final(Lf(f), res)
            ok = verify(f, res)
            if ok is True:
                print("    Check: matches SymPy's built-in Laplace transform.")
            elif ok is False:
                print("    Check: differs from SymPy's built-in result - please double-check.")
    except Exception as e:
        print(f"\nCould not solve: {e}")


def main():
    if len(sys.argv) > 1:
        solve(" ".join(sys.argv[1:]))
        return
    print(__doc__)
    while True:
        try:
            line = input("\nf(t) > ").strip()
        except (EOFError, KeyboardInterrupt):
            break
        if not line or line.lower() in {"q", "quit", "exit"}:
            break
        solve(line)


if __name__ == "__main__":
    main()