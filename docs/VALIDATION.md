# Validation notes

Numeric verification details pulled out of the main README to keep it
readable - the README says what's supported; this says how it was
checked and against what.

## Spherical (puream) vs Cartesian integrals

Engine originally shipped Cartesian-Gaussians-only (no spherical-harmonic
transform in `src/integrals/`); that mode is still available
(`spherical = .false.`) and was verified against Psi4 1.10 on HF/water
with `set puream false` (forcing Psi4 to Cartesian too), agreement 2e-9
to 3e-8 Hartree in every case tested: cc-pVDZ
-76.027074711/-76.027074711 (Engine/Psi4-cart), cc-pVQZ
-76.065004703/-76.065004699, def2-SVP -75.962181653/-75.962181658,
6-31+G* -76.017424679/-76.017424653.

Spherical (puream) support was added later (`spherical = .true.`, now
the default - commits `369396ae` ECP+puream energy path, `cc6364b4`
puream force/gradient path, `f5b9eb79` f/g (Lt=3,4) shells for both) and
is what conventionally-spherical basis families (Dunning cc-pVxZ/
aug-cc-pVxZ, Ahlrichs def2-*) actually mean by their published numbers -
matches Psi4's own default (`puream true`). Verified against Psi4 across
HF/PBE/LDA/TPSS/CAM-B3LYP, exact and DF/RI-K, with and without ECPs, and
with f/g-shell basis sets - core-Hamiltonian/overlap eigenvalues checked
for rotation invariance as the oracle for systems where total-energy
comparison alone is ambiguous (e.g. lanthanide multi-minima cases - see
git log around the commits above for the individual per-functional/
per-system numbers, not reproduced here to avoid drift from the actual
verification runs). For s/p-only bases (Pople sets: 6-31G*, 6-311G**
etc.) spherical vs Cartesian makes no difference at all - both
conventions agree exactly at l<=1.

## ECP import verification

`lanl2dz`/`sdd`'s extended element coverage (added via
`scripts/import_basis.py --ecp-only --nwchem-ecp`) was verified two
ways: (1) NaCl and TiCl4 under the newly-added lanl2dz Na/Cl/Ti agree
with Psi4 to 1e-9/4e-7 Hartree; (2) re-importing Pd's LANL2DZ ECP
through the new nwchem-ecp parser under a throwaway label reproduced
the original Gaussian94-format-parsed Pd energy bit-for-bit (-125.864922660
both ways) - Psi4 doesn't bundle SDD at all, so this cross-parser check
substitutes for a direct SDD/Psi4 comparison. Bare atomic Fe/Ti HF
energies do NOT reliably cross-validate this way (both codes can
converge to different, both-valid electronic states for open-shell
d-block atoms) - use a closed-shell molecule instead.

## `def2universaljkfit` as a general-purpose aux basis

Verified against exact (non-DF) 2e integrals on water/6-31G
(-76.298236 vs -76.298205 Hartree) and on Fe/LANL2DZ with an ECP
present (-122.539205 vs -122.539231 Hartree).

## Solvation (COSMO/CPCM + SMD)

Plain CPCM/gepol energy+force verified against Psi4 (CPCM/GePol/
Bondi, ~0.02-0.03% agreement). SMD's CDS term verified bit-exact (7+
significant figures, energy AND gradient) against pyscf's own SMD
module across several bonded-environment test cases (H-C, C-O, C-N
amide, C-Br). Full SMD (electrostatics+CDS) verified across water/
toluene/n-hexane (eps 78 down to 1.9): agreement ~0.02-0.04 kcal/mol in
energy, ~1e-4 to 1e-6 Hartree/Bohr in force - Engine deliberately
reuses the same conductor-scaled (C-PCM-style) electrostatics as plain
CPCM rather than true IEF-PCM, confirmed not to cost meaningful
accuracy across that full polarity range.

## `-D3`/`-D3BJ` swap to dynamically-linked simple-dftd3

Energy and analytic gradient verified bit-for-bit against the old
table-based implementation (both ZERO and BJ damping, water/PBE0/
6-31G) - see `docs/RELEASE_EXPORT.md` for why the swap was made and
`src/dispersion/mod_dispersion_d3.f90`'s header for the implementation.
