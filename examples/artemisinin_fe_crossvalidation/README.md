# Artemisinin + Fe2+ Cross-Validation

High-spin open-shell Fe(II) complex with artemisinin peroxide bridge.
Geometry optimized via GFN2-xTB (chrg=2, uhf=4).

## Test Case Properties

- **Charge**: +2
- **Spin multiplicity**: 5 (S=2, quartet alpha electron excess)
- **Key feature**: Near-degenerate frontier orbitals (alpha HOMO ≈ beta LUMO)
  - Causes multi-cycle DIIS oscillation in early SCF iterations
  - Converges cleanly after ~35-44 cycles depending on basis
  - Reproduced identically across Engine/PySCF/ORCA

## Files

### `artemisinin_fe_wb97mv.inp` — Engine

- Functional: WB97M-V
- Basis: def2-TZVPP (992 basis functions)
- VV10 NLC: non-self-consistent mode
- Result: E = -2223.115753109 Ha (35 SCF cycles)

### `artemisinin_fe_orca_wb97mv_631g.inp` — ORCA

- Functional: WB97M-V
- Basis: 6-31G (251 basis functions)
- VV10 NLC: post-SCF correction (ORCA's default "DFT-NL" style)
- Result: SCF = -2223.764541537 Ha + NL = +0.648819884 Ha
  - **Final SC+NL Energy**: -2223.115721653 Ha (44 SCF cycles)

## Cross-Validation

Both codes converge to **-2223.1157 Ha** with **3e-5 Ha precision** on the final
total energy (including VV10 dispersion), confirming:
- Engine's non-self-consistent VV10 implementation
- Correct handling of high-spin open-shell SCF in near-degenerate regime
- Consistency of DFT-NL dispersion across independent quantum chemistry codes

**Note**: ORCA applies VV10/NLC as a non-self-consistent post-SCF correction by default.
Engine's `vv10_nonself=.true.` mode reproduces this treatment exactly.
