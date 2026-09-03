# Artemisinin + Fe2+ Cross-Validation

High-spin open-shell Fe(II) complex with artemisinin peroxide bridge.
Geometry optimized via GFN2-xTB (chrg=2, uhf=4).

## Test Case Properties

- **Charge**: +2
- **Spin multiplicity**: 5 (S=2)
- **Key feature**: near-degenerate frontier orbitals (alpha HOMO ≈ beta LUMO)
  cause multi-cycle DIIS oscillation before the SCF settles. Direwolf
  auto-widens the DIIS window to 15 for any d/f-block metal; this system
  then converges in ~37 cycles (RI-JK) / ~51 (RIJCOSX).

## Files

### `artemisinin_fe_wb97mv.inp` — Direwolf, WB97M-V / def2-TZVPP (992 bf)

- RI-JK (`J = 'RI'`, `K = 'RI'`, `ri_aux_basis = 'def2universaljkfit'`)
- VV10 NLC: non-self-consistent (`vv10_nonself = .true.`)
- Result: **E = -2223.9316 Ha** (SCF -2224.57897, VV10 NLC +0.64739),
  37 SCF cycles. Reproduced to ~1e-4 Ha across a 20-commit span
  (2026-08-11 build -2223.931579, current HEAD RIJCOSX variant
  -2223.931719).

### `artemisinin_fe_wb97mv_def2tzvp.inp` — Direwolf, WB97M-V / def2-TZVP

Same setup, smaller basis. Run for a basis-sensitivity check.

### `artemisinin_fe_orca_wb97mv_631g.inp` — ORCA, WB97M-V / 6-31G (251 bf)

Independent code, much smaller basis. NOT directly comparable to the
def2-TZVPP total energy (basis-set incompleteness alone is many Ha);
useful only as a same-geometry SCF-stability / DFT-NL sanity check.

## History

An earlier version of this README recorded E = -2223.1158 Ha for the
def2-TZVPP Direwolf run and claimed 3e-5 Ha agreement with the 6-31G ORCA
number. Both were wrong: total energies in two different basis sets
cannot agree to 3e-5, and no Direwolf build (checked back to 2026-08-11)
converges this input anywhere near -2223.1158 — it converges to
-2223.9316. The stale number was corrected 2026-08-31.
