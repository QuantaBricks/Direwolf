# Engine development log

A high-level timeline, not a line-by-line changelog (see `git log` for
that). Status entries below are only for what's actually shipped and
verified against a reference code (Psi4/ORCA/pyscf), not aspirational.

## Time-line

- **Jul 2017** - alpha-1.0, initial creation.
- **Mar 2021** - alpha-1.1, recreated.
- **2025 - present** - active, AI-assisted continuous development.

## Milestones (2025 - present)

- **SCF core**: DirectSCF, density fitting (RI-J/RI-K), COSX (chain-of-
  spheres exchange, full-range K and RSH K_LR), independent J/K build-
  mode dispatch (exact/RI/COSX combinable per term).
- **DFT**: dynamic two-stage XC grid (coarse-then-production, DIRECT
  and STORE evaluation paths), VV10/B97M-V/WB97X-V/WB97X-3C and other
  meta-GGA/RSH functionals, B97-3C/R2SCAN-3C/WB97X-3C composite
  methods, D3/D4 dispersion and gCP/SRB corrections (dynamically
  linked against the vendored simple-dftd3/dftd4/gcp libraries).
- **Integrals**: spherical harmonics (puream) through f/g shells, ECP
  (effective core potentials, energy + analytic force), analytic
  Hessian (1e + CPHF + response-Hessian), analytic force for DF/COSX/
  VV10 paths.
- **Solvation**: COSMO (GePol cavity, analytic force, charge
  conservation) and SMD, cross-verified against ORCA/NWChem/Gaussian.
- **I/O**: native namelist/XYZ/column `&atoms` formats, ORCA/Gaussian
  input-compatibility readers, TOML export, checkpoint/restart,
  Molden output.
- **Correctness/perf infrastructure**: bit-reproducible parallel
  reductions (fixed-order OpenMP merges), DSYEVD diagonalization,
  disk-backed (mmap) density-fitting tensors for large systems.

## Not started yet

- Fast Multipole 2e integrals (MARI-J)
- DF (RI-J/RI-K) analytic Hessian
- meta-GGA analytic force, TDDFT, GHO QM/MM boundary
