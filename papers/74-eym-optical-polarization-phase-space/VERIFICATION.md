# Verification

## Analytic checks

- The compact differential preparations satisfy the linear Einstein and Gauss constraints exactly in an incoming flat subpatch; extension by zero preserves the constraints on the larger laboratory slice. Match both the configuration and its derivative at subleading order using the two eikonal branches so that the error is small enough for the wavelength-rescaled action current.
- The Einstein TT kinetic factor is 1/(4κ) in the momentum, and the color kinetic factor is one. The physical amplitude hTT/(2√κ) makes them equal. The imaginary Hermitian ray flux follows with the stated coordinate–momentum sign, without retaining a double-carrier term. A fixed finite packet space suffices for convergence of the bracket inverse.
- The independent background and detector computations retain the off-axis color transport and use rotation/reflection of the registered aperture. The equality is of the full leading 2×2 gravity covariance at one cut, not of spatially resolved amplitudes. The full two-cut Poisson bracket and the phase-coherent two-cut product are distinct; both reverse under the reflecting conjugation. Check the Weyl calibration and the exact compact TT drive, including the sign relating the drive response to the bracket.
- Earlier work on null symplectic geometry, projective rays, wave conversion and non-Abelian plane-wave backgrounds is credited in the manuscript. Papers 29B, 29C and 60 provide distinct inputs. No unrestricted Poisson inverse or quantum model follows from the bounded packet calculation.

## Build

Run ./build.sh with pdfLaTeX and latexmk installed. No external figures, dependencies or data files are needed. The release PDF is eight US-letter pages, with a blank author/affiliation and a populated title field. The source was compiled with no unresolved citations or references, and all pages were rendered for visual inspection. The results are analytical, not numerical.
