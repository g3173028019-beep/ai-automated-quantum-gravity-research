"""Compare the exact H-tree radial contraction and cotangent formula."""
import mpmath as mp

mp.mp.dps = 45
for r in (17, 18, 37, 38):
    odd = range(1, r, 2)
    theta = {a: mp.pi * a / r for a in odd}
    w = {a: mp.cot(theta[a] / 2) for a in odd}
    total = mp.fsum(w[a] * w[b] * mp.sin(theta[a] * b) for a in odd for b in odd)
    tau = (-1)**r * (mp.mpf(2) / r) * total / mp.sin(mp.pi / r)
    closed = (-1)**r * mp.fsum(w.values()) / mp.sin(mp.pi / r)
    error = abs(tau - closed)
    assert error < mp.mpf("1e-38")
    print(f"r={r}: absolute error {mp.nstr(error, 5)}")
