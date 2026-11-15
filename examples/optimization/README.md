# Geometry optimization examples

`runtype = 'opt'` switches a run from a single point to a geometry
minimization; the optional `&opt` namelist tunes the optimizer (every
setting has a default; giving the block also means optimize). Run these from the repository root, e.g.
`./Direwolf examples/optimization/water_opt.inp`.

Every cycle is appended to `<output>.opt.xyz` (multi-frame XYZ); the
converged geometry is also printed in the `.out` as `[FINAL GEOMETRY]`.
`restart` (default on) resumes from the last frame of an existing
`<output>.opt.xyz`, so delete that file to start fresh. The optimizer's current
Hessian is saved every cycle in `<output>.opt.hess` (Cartesian, Hartree/bohr^2) and
a resumed run with no Hessian source given starts from it instead of the model
(paracetamol, interrupted after 6 cycles: 8 more cycles instead of 15). A fresh
start removes it; a TS search does not use it.

| file | shows |
|---|---|
| `water_opt.inp` | the minimal case - `runtype = 'opt'` with an empty `&opt`, redundant internal coordinates, Lindh model Hessian |
| `water_opt_cart.inp` | `coord = 'cart'` - optimize in Cartesians instead |
| `water_opt_tight.inp` | `convergence = 'tight'` plus `max_steps` / `max_step_size`; tight also pulls the SCF up to `'tight'` automatically |
| `water_opt_thresholds.inp` | `convergence = 'loose'` (loose thresholds; also `'normal'`, `'tight'`) and per-criterion overrides `max_force`, `rms_force` (Hartree/bohr), `max_disp`, `rms_disp` (bohr), `energy_change` (Hartree); converged = gradient below both and (step below both or \|dE\| < energy_change). The thresholds in use are printed at the start |
| `water_opt_xtbhess.inp` | `hessian = 'read'` + `hessian_file` - seed from an external Cartesian Hessian (`water_gfn2.hessian`, produced by `xtb w.xyz --hess`) instead of the model (Gaussian readfc) |
| `water_opt_calcfc.inp` | `hessian = 'calc'` - analytic Hessian of the method itself at the first step (Gaussian calcfc) |
| `water_opt_xtb.inp` | `hessian = 'xtb'` - xtb Hessian at the first step, run on the current geometry with the job's charge and spin |
| `water_opt_recalc.inp` | `recomp_hess = 3` - analytic Hessian again every 3 cycles (Gaussian recalc=3) |
| `water_opt_calcall.inp` | `recomp_hess = 1` - analytic Hessian at every cycle (Gaussian calcall) |
| `water_opt_recalc_xtb.inp` | `hessian = 'xtb'` + `recomp_hess = 2` - xtb Hessian refreshed on the current geometry; every xtb field shown |
| `water_freeze.inp` | `constraints = 'atoms 1 2'` - atoms 1 and 2 stay exactly where they are, only atom 3 relaxes |
| `water_bond_constraint.inp` | `constraints = 'bond 1 2 1.10'` - the O-H bond is driven to 1.10 A and held there while everything else relaxes (`angle i j k deg`, `dihedral i j k l deg` likewise; no value = keep the starting value) |
| `water_scan.inp` | `constraints = 'bond 1 2 0.90 1.30 5'` - relaxed scan: the O-H bond from 0.90 to 1.30 A in 5 points (both ends included), everything else relaxes at each point; writes the energy table and `<base>_scan.xyz` (one frame per point). `angle` and `dihedral` likewise (degrees); one scan coordinate per run |
| `water_freeze_flag.inp` | a trailing `F` on a coordinate line freezes that atom (same as an `atoms` item in `constraints`; works in the column-format `&atoms` block) |
| `water_freeze_xyzfile.inp` + `water_freeze.xyz` | a trailing `F` on an xyz line (`xyzfile` input) freezes that atom |
| `water_freeze_bond.inp` | both together |
| `nh3_ts.inp` | `ts` - transition state (NH3 inversion, P-RFO), analytic start Hessian, frequencies at the end |
| `nh3_ts_follow.inp` | `ts_follow` - choose the mode to climb by the atoms/bond/angle/dihedral that move in it |
| `nh3_ts_xtb.inp` | `ts` with an xtb start Hessian (its negative curvature is kept) |
| `sn2_cl_ts_recalc.inp` | `ts` + `recomp_hess = 3` - Cl- + CH3Cl SN2 transition state, analytic Hessian refreshed every 3 cycles |
| `paracetamol_opt.inp` | a realistic 20-atom organic: r2SCAN-3c, 8 threads |

`&opt` fields (all optional): `max_steps` (100), `convergence`
(`'normal'` | `'tight'`), `max_step_size` (0.3 bohr), `coord`
(`'ric'` | `'cart'`), `restart` (`.true.`), `hessian_file`
(`''`), `write_ric` (`.false.`). See the top-level README's
"Geometry optimization" section for the convergence-threshold values.

Hessian fields:

| field | values | meaning |
|---|---|---|
| `hessian` | `'model'` (default), `'calc'`, `'read'`, `'xtb'` | Hessian source: Lindh model / analytic Hessian of this method / `hessian_file` / xtb. Recomputed Hessians (`recomp_hess`) use xtb when this is `'xtb'`, else the analytic one (external codes: xtb only). A TS search with `'model'` uses `'calc'` |
| `recomp_hess` | not set (default), `first`, N, `all` | `first`: analytic Hessian at the start only; N: recompute every N cycles; `all`: every cycle. A TS search takes the same values |
| `xtb_method` | `'gfn2'` (default), `'gfn1'`, `'gfn0'`, `'gfnff'` | xtb Hamiltonian |
| `xtb_cmd` | `'xtb'` | xtb executable (full path if it is not on PATH) |
| `xtb_args` | `''` | extra xtb arguments, e.g. `'--alpb water'` for a solvated run |
| `constraints` | `'atoms 1-6 9 12; bond 1 2; angle 3 1 2 109.5; dihedral 4 3 1 2 60'` | `atoms`: atoms that do not move (1-based). bond / angle / dihedral held at the starting value, or driven to the given end value (Angstrom, degrees) at most 0.1 bohr or 0.1 rad per cycle; `;` separated. Works in both `ric` and `cart`. A bond/angle/dihedral with three numbers after the atoms, `bond 1 2 0.9 1.3 5`, is a scan: start, end, number of points (see `water_scan.inp`) |
| `ts` | `.false.` | transition-state search |

xtb runs in `<output stem>.xtb/` (`geom.xyz`, `xtb.out`, `hessian`) on
the current geometry with `--chrg` and `--uhf` taken from the job. In a
minimization a computed Hessian gets the gradient term of the
Cartesian-to-internal transformation, and negative eigenvalues are
flipped (an xtb Hessian is also floored at 5e-3 in Cartesians); in a TS
search the negative curvature is kept.

Measured (8 threads; paracetamol r2SCAN-3c, TS cases PBE0/def2-SVP RIJCOSX):

| case | model | calc | xtb | recalc 3 (calc) | recalc 3 (xtb) | calcall |
|---|---|---|---|---|---|---|
| paracetamol, cycles / s | 20 / 138 | 12 / 163 | 18 / 136 | 17 / 489 | 21 / 157 | 9 / 553 |
| butadiene + ethylene TS, cycles / s | - | 10 / 199 | 11 / 130 | | | |
| Cl- + CH3Cl TS, cycles / s | - | 5 / 18 | 6 / 14 | 5 / 23 | | |

For an ordinary minimization the model Hessian is usually the fastest in
wall time; `'calc'` saves cycles and pays for one Hessian. `'xtb'` is
cheap and a good start Hessian for a TS search.

To make your own seed Hessian for `hessian_file`: any code that
writes a 3N x 3N Cartesian Hessian in Hartree/bohr² works;
`xtb <geom>.xyz --hess` writes one (`$hessian`, 5 values per line) that
is read as-is.
