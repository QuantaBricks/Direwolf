# Geometry optimization examples

`runtype = 'opt'` switches a run from a single point to a geometry
minimization; the optional `&opt` namelist tunes the optimizer (every
setting has a default; giving the block also means optimize). Run these from the repository root, e.g.
`./Direwolf examples/optimization/water_opt.inp`.

Every cycle is appended to `<output>.opt.xyz` (multi-frame XYZ); the
converged geometry is also printed in the `.out` as `[FINAL GEOMETRY]`.
`opt_restart` (default on) resumes from the last frame of an existing
`<output>.opt.xyz`, so delete that file to start fresh. The optimizer's current
Hessian is saved every cycle in `<output>.opt.hess` (Cartesian, Hartree/bohr^2) and
a resumed run with no Hessian source given starts from it instead of the model
(paracetamol, interrupted after 6 cycles: 8 more cycles instead of 15). A fresh
start removes it; a TS search does not use it.

| file | shows |
|---|---|
| `water_opt.inp` | the minimal case - `runtype = 'opt'` with an empty `&opt`, redundant internal coordinates, Lindh model Hessian |
| `water_opt_cart.inp` | `opt_coord = 'cart'` - optimize in Cartesians instead |
| `water_opt_tight.inp` | `opt_conv = 'tight'` plus `opt_maxcyc` / `opt_trust`; tight also pulls the SCF up to `'tight'` automatically |
| `water_opt_thresholds.inp` | `opt_conv = 'loose'` (loose thresholds; also `'normal'`, `'tight'`) and per-criterion overrides `opt_gmax`, `opt_grms` (Hartree/bohr), `opt_smax`, `opt_srms` (bohr), `opt_de` (Hartree); converged = gradient below both and (step below both or \|dE\| < opt_de). The thresholds in use are printed at the start |
| `water_opt_xtbhess.inp` | `opt_hess = 'read'` + `opt_hessian_file` - seed from an external Cartesian Hessian (`water_gfn2.hessian`, produced by `xtb w.xyz --hess`) instead of the model (Gaussian readfc) |
| `water_opt_calcfc.inp` | `opt_hess = 'calc'` - analytic Hessian of the method itself at the first step (Gaussian calcfc) |
| `water_opt_xtb.inp` | `opt_hess = 'xtb'` - xtb Hessian at the first step, run on the current geometry with the job's charge and spin |
| `water_opt_recalc.inp` | `opt_hess_every = 3` - analytic Hessian again every 3 cycles (Gaussian recalc=3) |
| `water_opt_calcall.inp` | `opt_hess_every = 1` - analytic Hessian at every cycle (Gaussian calcall) |
| `water_opt_recalc_xtb.inp` | `opt_hess = 'xtb'` + `opt_hess_every = 2` - xtb Hessian refreshed on the current geometry; every xtb field shown |
| `water_freeze.inp` | `opt_freeze = '1 2'` - atoms 1 and 2 stay exactly where they are, only atom 3 relaxes |
| `water_bond_constraint.inp` | `opt_constraints = 'bond 1 2 1.10'` - the O-H bond is driven to 1.10 A and held there while everything else relaxes (`angle i j k deg`, `dihedral i j k l deg` likewise; no value = keep the starting value) |
| `water_freeze_flag.inp` | a trailing `F` on a coordinate line freezes that atom (same as listing it in `opt_freeze`; works in the column-format `&atoms` block) |
| `water_freeze_xyzfile.inp` + `water_freeze.xyz` | a trailing `F` on an xyz line (`xyzfile` input) freezes that atom |
| `water_freeze_gaussian.gjf` | Gaussian input: freeze code `-1` after the element freezes it; with `Opt=ModRedundant` the section after the coordinates is read: `X n F` freezes an atom, `B i j F` / `A i j k F` / `D i j k l F` freeze a bond / angle / dihedral (a value before the `F` sets it). `S` scans and other actions are ignored with a note |
| `water_freeze_bond.inp` | both together |
| `nh3_ts.inp` | `opt_ts` - transition state (NH3 inversion, P-RFO), analytic start Hessian, frequencies at the end |
| `nh3_ts_xtb.inp` | `opt_ts` with an xtb start Hessian (its negative curvature is kept) |
| `sn2_cl_ts_recalc.inp` | `opt_ts` + `opt_hess_every = 3` - Cl- + CH3Cl SN2 transition state, analytic Hessian refreshed every 3 cycles |
| `paracetamol_opt.inp` | a realistic 20-atom organic: r2SCAN-3c, 8 threads |

`&opt` fields (all optional): `opt_maxcyc` (100), `opt_conv`
(`'normal'` | `'tight'`), `opt_trust` (0.3 bohr), `opt_coord`
(`'ric'` | `'cart'`), `opt_restart` (`.true.`), `opt_hessian_file`
(`''`), `opt_write_ric` (`.false.`). See the top-level README's
"Geometry optimization" section for the convergence-threshold values.

Hessian fields:

| field | values | meaning |
|---|---|---|
| `opt_hess` | `'model'` (default), `'calc'`, `'read'`, `'xtb'` | start Hessian: Lindh model / analytic Hessian of this method / `opt_hessian_file` / xtb. A TS search with `'model'` uses `'calc'` |
| `opt_hess_every` | 0 (default), N | recompute the Hessian every N cycles; 1 = every cycle |
| `opt_hess_recalc` | `'calc'`, `'xtb'` | source for the recomputed Hessians; default `'xtb'` when `opt_hess = 'xtb'`, else `'calc'`. External codes: xtb only |
| `opt_xtb_method` | `'gfn2'` (default), `'gfn1'`, `'gfn0'`, `'gfnff'` | xtb Hamiltonian |
| `opt_xtb_cmd` | `'xtb'` | xtb executable (full path if it is not on PATH) |
| `opt_xtb_args` | `''` | extra xtb arguments, e.g. `'--alpb water'` for a solvated run |
| `opt_freeze` | `'1-6 9 12'` | atoms that do not move (1-based) |
| `opt_constraints` | `'bond 1 2; angle 3 1 2 109.5; dihedral 4 3 1 2 60'` | bond / angle / dihedral held at the starting value, or driven to the given end value (Angstrom, degrees) at most 0.1 bohr or 0.1 rad per cycle; `;` separated. Works in both `ric` and `cart` |
| `opt_ts`, `opt_ts_hess_every` | `.false.`, 0 | transition-state search; `opt_ts_hess_every` is the older name of `opt_hess_every` |

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

To make your own seed Hessian for `opt_hessian_file`: any code that
writes a 3N x 3N Cartesian Hessian in Hartree/bohr² works;
`xtb <geom>.xyz --hess` writes one (`$hessian`, 5 values per line) that
is read as-is.
