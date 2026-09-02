# Engine environment variables

Every knob below is read from the process environment at run time. They
fall into two groups: a handful of **tuning switches** that change how a
production calculation runs, and a larger set of **debug/diagnostic**
switches that dump intermediates or force a code path for testing.

Calculation *settings* — basis, functional, charge, multiplicity, RI on/off,
aux basis, force on/off, memory cap, thread count, COSMO/SMD, checkpoint
and molden files — are **not** environment variables. They belong in the
`&molecule` namelist of the input file. A few env vars that used to
control those (`ENGINE_USE_DF`, `ENGINE_DF_AUX_BASIS`, `ENGINE_HF_FRAC`,
`ENGINE_SKIP_FORCE`, `ENGINE_MAX_TWOEI_GB`, `ENGINE_MAX_GRID_CACHE_GB`)
are **no longer read** — they survive only in source comments.

Unset is always a valid, sensible default. Anything marked *(default)*
below is what you get with the variable unset.

Every tuning switch below can also be set from the input file itself, in
an optional `&professional` namelist, instead of the shell environment -
see the "Input-file equivalent" section at the end of this doc. The two
are equivalent; `&professional` just sets the same env var for you before
the run starts (and wins over anything already in the shell environment).

---

## XC quadrature: STORE vs DIRECT, and the two-stage grid

`DFT_calc` can either cache every basis function's value at every grid
point (**STORE** — fast per iteration, large memory) or recompute them
per batch each iteration (**DIRECT** — O(1) in grid size, historically
slower). The two-stage grid closed most of DIRECT's speed gap, so DIRECT
is now the default for larger systems.

| Variable | Values | Meaning |
|---|---|---|
| `ENGINE_GRID_CACHE` | `0` / `1` | Force DIRECT (`0`) or STORE (`1`). Unset *(default)*: DIRECT when `nconts > 500`, STORE below. The threshold is a hard basis-function count, chosen to reserve memory headroom for J/K predictably. |
| `ENGINE_XC_DYNGRID` | `0` / `1` | Two-stage grid: run early SCF iterations on a radially reduced grid, then promote to the production grid. Unset *(default)*: on whenever DIRECT is active. Ignored in STORE mode (there the switch must free and refill the value cache, which costs more than the coarse stage saves). |
| `ENGINE_XC_DYNGRID_PRMS` | real | Density-RMS threshold at which to promote. Default `3e-4`. Switching *later* is monotonically worse — the coarse grid's converged solution is a fixed point of the wrong operator. |
| `ENGINE_XC_DYNGRID_RSCALE` | real | Coarse-stage radial scale factor. Default `0.60` (production is `1.848`). |
| `ENGINE_XC_DYNGRID_SPH` | `434`/`302`/`230`/`170` | Coarse-stage base Lebedev order. Default `434`, i.e. **unchanged from production**. Reducing it is a trap: 434→302 costs a fixed 5-6 extra SCF iterations at every radial scale tried. Take the saving from the radial axis only. |

The promotion happens *before* that iteration's `DFT_calc`, so the
converged energy is defined by the production grid alone (verified to
within ~1e-10 Ha against a fixed-fine-grid run).

## Alternative XC build paths (experimental, off by default)

| Variable | Values | Meaning |
|---|---|---|
| `ENGINE_XC_SHELLPAIR` | `1` | Shell-pair-driven XC build instead of the default batch-driven one. |
| `ENGINE_XC_PSCREEN` | real | Density-change screening threshold for the shell-pair path. |
| `ENGINE_XC_INCR` | `1` | Incremental Vxc against a stored reference density (shell-pair path only). |

## SCF

| Variable | Values | Meaning |
|---|---|---|
| `ENGINE_HARRIS_GUESS` | `0` / `1` | Force the SAD/atomic-HF initial guess off or on. Unset uses the driver's `harris_guess` namelist setting. |
| `ENGINE_INCREMENTAL_FOCK` | `0` / `1` | Incremental (Almlöf) Fock builds. On by default; `0` disables it on the RI-J path. |
| `ENGINE_DIIS_WINDOW` | integer | DIIS history depth. Default `10`. |

## Density fitting / force mode overrides

These override the automatic memory-based STORE-vs-DIRECT decision for
the DF tensors. Useful for reproducing a mode on a machine where the
estimate would pick the other one.

| Variable | Values | Meaning |
|---|---|---|
| `ENGINE_FORCE_DF_DIRECT` | `1` | Force DF energy path to DIRECT. |
| `ENGINE_FORCE_DF_STORE` | `1` | Force DF energy path to STORE. |
| `ENGINE_FORCE_DF_FORCE_DIRECT` | `1` | Force the DF *gradient* path to DIRECT. |

## wB97M-V and the rest of the Head-Gordon VV10 family (WB97X-V, B97M-V)

`functional = 'WB97M-V'` in the `&molecule` namelist selects Mardirossian &
Head-Gordon's wB97M-V (J. Chem. Phys. 144, 214110 (2016)): libxc's
`XC_HYB_MGGA_XC_WB97M_V` supplies the semilocal exchange-correlation and
range-separation (omega=0.3, 15% SR-HF/100% LR-HF, taken from libxc's own
`hyb_cam_coef`), and Engine's own `mod_vv10` supplies the VV10 nonlocal
correlation term with the functional's own fit **b=6.0, C=0.01** (not the
original Vydrov/Van Voorhis b=5.9/C=0.0093, which is only correct for the
plain "VV10"/rPW86PBE+VV10 case a few lines above it in `mod_xc.f90`).
Verified against real Psi4 (`scf_type df`, `df_basis_scf cc-pvdz-jkfit` —
**must** match Engine's own `J='RI'/K='RI'`/aux-basis settings, PK/exact
integrals on the Psi4 side are not a valid comparison) on caffeine/6-311G*:
total energy agrees to 3.2e-5 Ha.

Two siblings share this exact same VV10 machinery (same b=6.0/C=0.01 fit,
same `ENGINE_VV10_*` tuning below, same `J='RI'/K='RI'`/aux-basis requirement for
a valid Psi4/ORCA comparison) and were added the same way:

- **`WB97X-V`** (`XC_HYB_GGA_XC_WB97X_V`) — the earlier (2014), non-meta-GGA
  member: range-separated hybrid GGA (no tau) + VV10. Verified vs real
  Psi4 (self-consistent VV10) on caffeine/6-311G*: total energy agrees to
  2.4e-5 Ha at each code's own default XC grid, tightening to 2.3e-6 Ha
  once both sides' grids are pushed well past their defaults — the
  residual at default settings is ordinary cross-code quadrature noise for
  this "highly-parameterized" functional family (same category as the
  M06-2X/WB97X-D residual documented in `mod_xc.f90`), not an
  implementation bug.
- **`B97M-V`** (`XC_MGGA_XC_B97M_V`) — the non-range-separated sibling: pure
  meta-GGA (0% exact exchange, no RSH) + VV10, same b=6.0/C=0.01 fit taken
  without further adjustment per the 2015 paper. Verified vs real Psi4 on
  caffeine/6-311G*: total energy agrees to 1.3e-5 Ha. Its total energy
  differs from wB97M-V's by ~0.16 Ha on the same molecule/geometry — this
  is NOT the isolated cost of dropping the long-range correction (the two
  functionals are independently fit, not the same parameterization with
  omega toggled), it's two different functionals. Use wB97M-V/wB97X-V (not
  B97M-V) whenever the calculation is sensitive to long-range
  self-interaction error — charge-transfer/Rydberg excitations, reaction
  barrier heights, stretched/breaking bonds, radicals/anions; for a plain
  closed-shell equilibrium-geometry single point, B97M-V is a reasonable,
  cheaper choice.

Example namelist:

    functional   = 'WB97M-V'
    J            = 'RI'
    K            = 'RI'
    ri_aux_basis = 'ccpvdzjkfit'

The VV10 term needs its own nonlocal-correlation (NLC) quadrature grid —
separate from, and much coarser than, the main XC grid — for the O(N_pts^2)
kernel double sum. It's built once per run (atoms are fixed during SCF) and
cached.

| Variable | Values | Meaning |
|---|---|---|
| `ENGINE_VV10_RADIAL` | integer | NLC grid radial points per atom. Default `40`. |
| `ENGINE_VV10_ANGULAR` | Lebedev order (6,14,26,...,974) | NLC grid angular points per atom. Default `170` (163200 pts on caffeine's 24 atoms). |
| `ENGINE_VV10_RCUT` | bohr | Atom-block pair cutoff for the kernel sum; blocks farther apart are skipped. Default `20.0`. |
| `ENGINE_VV10_RHOCUT` | real | Points with total density below this are dropped from the kernel sum. Default `1e-9`. |

Default 40/170 was picked by direct comparison against Psi4's own wB97M-V,
which reuses its full production XC grid (75 radial x 302 angular, ~500K
points on caffeine) for the NLC pass rather than a separate coarse one —
copying that density is ~50x more expensive than 40/170 for a further
~1e-5 Ha of agreement, and isn't actually necessary to clear the accuracy
bar. The kernel sum is O(N_pts^2), so grid density trades accuracy for
speed steeply: the old 20/50 default (24000 pts) gave only 1.7e-4 Ha
agreement; halving 40/170 to 40/86 (82560 pts) gave 4.67e-5 Ha — passes a
<5e-5 Ha bar but with very thin margin, so it was not adopted as the
default. A separate closed-source implementation's own NL grid
(Lebedev-110, ~37K points) was checked too and is *not* a valid target
to copy: its own result on that grid is 2.5e-4 Ha off from Psi4, i.e.
less accurate than Engine's old 20/50 default.

Two speedup attempts on the kernel sum were tried and rejected (measured,
not guessed): an explicit `!$omp simd` hint on the inner loop (zero
effect — the compiler already auto-vectorizes it under the module's
`-funsafe-math-optimizations -march=native` flags), and a per-pair
distance early-exit inside the `allscreen` branch, analogous to Schwarz
screening (measured a ~1.7x *regression* — for a compact molecule whose
diameter is close to `ENGINE_VV10_RCUT`, few pairs are actually prunable,
and the added branch broke auto-vectorization for all of them). The O(N^2)
kernel is at its practical floor for this grid representation; a real
further speedup would need a different algorithm (FMM/multipole), which
is out of scope as a quick optimization.

## Input-file equivalent: the `&professional` namelist

Every tuning switch above (not the debug/diagnostic ones below - those
stay env-var-only) can be set directly in the input file instead of the
shell environment, in an optional `&professional` block. `run_engine.f90`
parses it and pushes each field you actually set to the corresponding
`ENGINE_*` environment variable (via a `setenv()` C binding, `overwrite=1`)
before the run starts - it does not change how any of those variables are
read; every module still reads its own `ENGINE_*` var exactly as
documented above. A field left out of the block is left alone (same as
the env var being unset); there is no "0 means default" ambiguity because
Fortran namelist I/O has no way to tell "not present" from "set to zero"
otherwise, so this is why the mapping goes through explicit sentinels
internally.

| `&professional` field | Env var it sets |
|---|---|
| `prof_grid_cache` (`0`/`1`) | `ENGINE_GRID_CACHE` |
| `prof_xc_dyngrid` (`0`/`1`) | `ENGINE_XC_DYNGRID` |
| `prof_xc_dyngrid_prms` (real) | `ENGINE_XC_DYNGRID_PRMS` |
| `prof_xc_dyngrid_rscale` (real) | `ENGINE_XC_DYNGRID_RSCALE` |
| `prof_xc_dyngrid_sph` (integer) | `ENGINE_XC_DYNGRID_SPH` |
| `prof_xc_shellpair` (`1`) | `ENGINE_XC_SHELLPAIR` |
| `prof_xc_pscreen` (real) | `ENGINE_XC_PSCREEN` |
| `prof_xc_incr` (`1`) | `ENGINE_XC_INCR` |
| `prof_incremental_fock` (`0`/`1`) | `ENGINE_INCREMENTAL_FOCK` |
| `prof_diis_window` (integer) | `ENGINE_DIIS_WINDOW` |
| `prof_force_df_direct` (`1`) | `ENGINE_FORCE_DF_DIRECT` |
| `prof_force_df_store` (`1`) | `ENGINE_FORCE_DF_STORE` |
| `prof_force_df_force_direct` (`1`) | `ENGINE_FORCE_DF_FORCE_DIRECT` |
| `prof_vv10_radial` (integer) | `ENGINE_VV10_RADIAL` |
| `prof_vv10_angular` (Lebedev order) | `ENGINE_VV10_ANGULAR` |
| `prof_vv10_rcut` (bohr) | `ENGINE_VV10_RCUT` |
| `prof_vv10_rhocut` (real) | `ENGINE_VV10_RHOCUT` |
| `prof_xc_grid_level` (`2`-`7`) | `ENGINE_XC_GRID_LEVEL` (picks a matched coarse+fine PAIR) |
| `prof_xc_grid_level_coarse` (`2`-`7`) | `ENGINE_XC_GRID_LEVEL_COARSE` (coarse side only, independent of fine) |
| `prof_xc_grid_level_fine` (`2`-`7`) | `ENGINE_XC_GRID_LEVEL_FINE` (fine/production side only, independent of coarse) |
| `prof_cosx_hi_level` (`1`/`2`/`3`) | sets `ENGINE_COSX_NRAD`/`NSPH` from a preset (coarser/**default**/finer than COSX's "hi" grid, 18/194) |
| `prof_cosx_md_level` (`1`/`2`/`3`) | sets `ENGINE_COSX_NRAD_MD`/`NSPH_MD` from a preset (coarser/**default**/finer than "md", 14/110) |
| `prof_cosx_lo_level` (`1`/`2`/`3`) | sets `ENGINE_COSX_NRAD_LO`/`NSPH_LO` from a preset (coarser/**default**/finer than "lo", 14/50) |

`ENGINE_HARRIS_GUESS` has no `&professional` field: `&molecule`'s own
`harris_guess` field already covers it directly.

Example - reproduce the old always-STORE behaviour and a shorter DIIS
window from the input file, no shell env needed:

    &professional
     prof_grid_cache  = 1
     prof_diis_window = 6
    &end

## Debug and diagnostics

None of these change results; they print or dump intermediates.

| Variable | Meaning |
|---|---|
| `ENGINE_DIIS_DEBUG=1` | Per-iteration DIIS/ADIIS trace. |
| `ENGINE_GTO_PROGRESS=1` | Progress reporting while filling the grid value cache. |
| `ENGINE_FORCE_DEBUG_PARTS=1` | Per-term gradient breakdown. |
| `ENGINE_HESSIAN_DEBUG1E=1` | One-electron Hessian terms. |
| `ENGINE_HESSIAN_DEBUG2E_DIRECT=1` | Two-electron Hessian, direct path. |
| `ENGINE_HESSIAN_TOTAL_DEBUG=1` | Total analytic Hessian. |
| `ENGINE_CPHF_DEBUG=1` | CPHF solver trace. |
| `ENGINE_CPHF_DUMP_DP=1` | Dump CPHF dP/dR. |
| `ENGINE_CPHF_GRAD_CHECK=1` | CPHF gradient consistency check. |
| `ENGINE_DUMP_S=1` | Dump the overlap matrix. |
| `ENGINE_DUMP_PTOT=1` | Dump the converged total density. |
| `ENGINE_DUMP_E1FIXED=1` | Dump the fixed one-electron energy terms. |
| `ENGINE_DF_RAW=1` | Dump raw DF three-index quantities (converged runs only). |
| `ENGINE_CCBK_RAW=1` | Dump raw CCBK quantities (converged runs only). |
| `ENGINE_RESP_DUMP=<file>` | With `resp_charges_on`, write the RESP fit's Merz-Kollman grid to `<file>.grid` (Angstrom, one `x y z` per line) and the fitted ESP to `<file>.esp` (a.u., one value per line). For cross-code checks: comparing only the final charges cannot separate a grid-generator difference from a density difference from a fitting-algorithm difference, but feeding this dump to the other code's fitter can. |

---

## Recipes

Reproduce the old always-STORE behaviour:

    ENGINE_GRID_CACHE=1 ./Engine job.inp job.out

Run DIRECT on a small system, with the two-stage grid (it follows
DIRECT automatically):

    ENGINE_GRID_CACHE=0 ./Engine job.inp job.out

Run DIRECT but with a fixed production grid throughout, e.g. to measure
what the two-stage grid is worth:

    ENGINE_GRID_CACHE=0 ENGINE_XC_DYNGRID=0 ./Engine job.inp job.out

Scan the coarse stage:

    ENGINE_GRID_CACHE=0 ENGINE_XC_DYNGRID_RSCALE=0.40 \
    ENGINE_XC_DYNGRID_PRMS=1e-3 ./Engine job.inp job.out
