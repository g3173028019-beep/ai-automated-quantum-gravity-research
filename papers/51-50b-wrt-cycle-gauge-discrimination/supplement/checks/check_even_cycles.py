"""Reproduce the even-cycle extension independently from modular S,T.

No network calls.  Run with OPENBLAS_NUM_THREADS=1 python check_even_cycles.py.
Finite-level identities are numerical checks, not a substitute for the proof.
"""
from math import comb, log, pi
import json
from pathlib import Path
import numpy as np


def reduced_eigenvalues(r):
    a = np.arange(1, r, 2, dtype=float)
    s = np.minimum(a, r-a)
    f = (pi / (2*r*np.sin(pi*a/(2*r))))**2
    matrix = np.sqrt(f[:, None]*f[None, :])*np.minimum(s[:, None], s[None, :])
    return np.linalg.eigvalsh(matrix)


def modular_eigenvalues(r):
    a = np.arange(1, r, dtype=float)
    smat = np.sqrt(2/r)*np.sin(pi*np.outer(a, a)/r)
    twist_squared = np.exp(1j*pi*(a*a-1)/r)
    cap = np.abs(smat @ (twist_squared*smat[:, 0]))**2
    weight = cap/smat[0, :]**2
    transfer = smat*np.sqrt(weight[:, None]*weight[None, :])
    return np.linalg.eigvalsh(transfer)


def run():
    checks = []
    for r in range(3, 27):
        raw, red = modular_eigenvalues(r), reduced_eigenvalues(r)
        for m in (2, 3, 4):
            actual = float(np.sum(raw**(2*m)))
            reduced = float((2*r*r/pi**2)**m * np.sum(red**m))
            relative = abs(actual-reduced)/max(1., abs(actual))
            assert relative < 3e-11, (r, m, actual, reduced, relative)
            checks.append({'r': r, 'cycle_length': 2*m, 'relative_error': relative})
    rows = []
    for r in (31, 32, 101, 102, 301, 302, 801, 802, 1601, 1602):
        red = reduced_eigenvalues(r)
        row = {'r': r}
        for m in (2, 3, 4, 5):
            trace = float(np.sum(red**m))
            leading = comb(2*m-2, m-1)
            row[f'C{2*m}_ratio_to_leading'] = 2**m*trace/(leading*log(r))
            row[f'C{2*m}_bounded_remainder'] = 2**m*trace-leading*log(r)
        rows.append(row)
    output = {'identity_max_relative_error': max(x['relative_error'] for x in checks),
              'identity_checks': len(checks), 'asymptotics': rows}
    destination = Path(__file__).with_name('even_cycle_checks.json')
    destination.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    run()
