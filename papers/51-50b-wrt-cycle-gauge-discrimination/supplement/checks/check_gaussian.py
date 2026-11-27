"""Independent quadrature check of the exact P4 Gaussian normal-model formula."""
import numpy as np
from numpy.polynomial.legendre import leggauss

x, w = leggauss(512)
x = (x + 1.0) / 2.0
w = w / 2.0
s, t = x[:, None], x[None, :]
weights = w[:, None] * w[None, :]
constant = 4 * np.pi**2 / 3 * (4 - np.pi + 2 * np.log(2))
print(f"Predicted constant term after multiplying by R^2: {constant:.12f}")
print("R, R^2 I/(8 pi^2 log R), R^2 I - 8 pi^2 log R")
for R in (10, 100, 1000, 10000, 100000, 1000000):
    integrand = (
        R * t**2 / np.sqrt(1 + R**2 * t**2)
        * np.arcsinh(
            R * np.sqrt(1 + R**2 * t**2)
            / np.sqrt(1 + R**2 * (s**2 + t**2))
        )
    )
    scaled = 16 * np.pi**2 * np.sum(weights * integrand)
    print(R, scaled / (8 * np.pi**2 * np.log(R)),
          scaled - 8 * np.pi**2 * np.log(R))
