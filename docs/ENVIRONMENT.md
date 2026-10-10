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
is now the unconditional default (2026-09-19 — dropped the earlier
`nconts > 500` gate that left every smaller system on STORE).

| Variable | Values | Meaning |
|---|---|---|
| `ENGINE_GRID_CACHE` | `0` / `1` | Force DIRECT (`0`) or STORE (`1`). Unset *(default)*: DIRECT, regardless of system size. |
| `ENGINE_XC_DYNGRID` | `0` / `1` | Two-stage grid: run early SCF iterations on a radially reduced grid, then promote to the production grid. Unset *(default)*: on when `nconts > 500`, independent of the STORE/DIRECT choice above (decoupled from it 2026-08-26 — see the comment at its call site in engine.f90 for why: a small system's STORE cache is cheap enough to just build once at full density, so the coarse-then-promote dance isn't worth it there regardless of which mode is active). A small system now runs plain single-stage DIRECT by default; DYNGRID kicking in is what actually changes above the 500-function line. |
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

## XC grid theory-table (default ON for most functionals since 2026-09-19)

A second, per-functional XC grid construction that replaces the classic
COARSE/FINE/L4-L7 tiers with a LEVEL (functional family) x TIER
(LO/REGULAR/HIGH accuracy) lookup, calibrated against reference point
budgets on a 33-molecule period-2 set plus a handful of individually
validated heavier elements.

**This is now the default construction for most functionals** (flipped
2026-09-19, once `xcgrid_zone_calibrated_for()` closed the last
correctness gap - see "Safe fallback" below). `ENGINE_XC_USE_THEORY_TABLE`
unset entirely means: ON for anything with a real, measured LEVEL
(everything under SMOOTHDFT's case-default - B3LYP/PBE/PBE0/BLYP/BP86/
TPSS/CAM-B3LYP/etc - plus R2SCAN/R2SCAN-3C/WB97X-D/WB97M-V and, since
2026-09-19/30, M06/M06-2X/M06-L); `xcgrid_theory_default_on()`
(mod_data.f90) is the authoritative list and currently has no off-list.
R2SCAN-3C has its own LEVEL `R2SCAN3C` (2026-10-06): LO and HIGH are
R2SCAN's, REGULAR keeps R2SCAN's angular zones with R2SCAN/HIGH's radial
scales, because the mTZVPP basis needs a finer radial grid than the def2-SVP
the R2SCAN cells were calibrated on (paracetamol 7.3e-5 -> 2.8e-6 Ha against
an L7 grid, omeprazole 1.5e-4 -> 1.0e-5, +26% wall time).
An explicit `ENGINE_XC_USE_THEORY_TABLE=1`/`=0` always overrides this
either direction, for anyone who wants theory-table on M06 anyway (not
unsafe, just unmeasured) or wants it off entirely regardless of
functional.

| Variable | Values | Meaning |
|---|---|---|
| `ENGINE_XC_USE_THEORY_TABLE` | `1` / `0` | Force on or off. Unset: on for a measured LEVEL, off otherwise - see above. |
| `ENGINE_XC_THEORY_LEVEL` | `SMOOTHDFT`/`R2SCAN`/`R2SCAN3C`/`WB97XD`/`WB97MV`/... | Functional family. Auto-inferred from the input file's `functional` if unset. |
| `ENGINE_XC_THEORY_TIER` | `LO`/`REGULAR`/`HIGH` | Accuracy tier. Default `REGULAR`. |
| `ENGINE_XC_THEORY_DYNAMIC` | `1` | Promote LO→REGULAR mid-SCF (mirrors `ENGINE_XC_DYNGRID` above, theory-table-specific). Off by default. |

**Coverage as of 2026-09-19** (see `data/xcgrid_theory_levels.h` and
`src/core/mod_data.f90`'s `xcgrid_theory_apply()` for the authoritative,
up-to-date list — this table rots faster than the code):

| | row1-2 (main group) | Fe | Ru | other row3-7 |
|---|---|---|---|---|
| SMOOTHDFT (GGA/hybrid-GGA) | calibrated | calibrated | calibrated | safe fallback (= production) |
| R2SCAN | calibrated | calibrated | calibrated (via safe fallback, see below) | safe fallback (= production) |
| WB97X-D | calibrated | calibrated (1 geometry) | calibrated (1 geometry) | safe fallback (= production) |
| WB97M-V | calibrated | calibrated (2 geometries) | calibrated (1 geometry) | safe fallback (= production) |

"Calibrated" beats production default on the geometry(ies) tested but
is not guaranteed to generalize — Fe's WB97M-V/R2SCAN fixes, for
instance, meet their own absolute accuracy target on a second geometry
(`fe_hexaaqua`) but lose to production default there on some tiers even
though they win clearly on the first (`ferrocene`).

**"Safe fallback" (added 2026-09-19, `xcgrid_zone_calibrated_for()` in
`mod_data.f90`)**: any element with no per-element zone/radial override
under the active LEVEL is angular-pruned through the exact same
production-default construction as theory-table-off, instead of
inheriting row1-2's calibrated zone table (which was never validated to
extend past period 2). This is what closed R2SCAN+Ru specifically — four
different direct overrides all made it worse than production, because
the real defect was structural (borrowing row1-2's zone breakpoints for
a 4d metal), not a tuning problem. Fixed at the mechanism level, Ru's
REGULAR tier now meets target and HIGH beats production by 228x with
*zero* manual per-element numbers. Net effect: turning theory-table on
can no longer regress an element nobody has individually measured yet -
worst case it behaves exactly like production default.

## SCF

| Variable | Values | Meaning |
|---|---|---|
| `ENGINE_HARRIS_GUESS` | `0` / `1` | Force the SAD/atomic-HF initial guess off or on. Unset uses the driver's `harris_guess` namelist setting. |
| `ENGINE_INCREMENTAL_FOCK` | `0` / `1` | Incremental (Almlöf) Fock builds. On by default; `0` disables it on the RI-J path. |
| `ENGINE_DIIS_WINDOW` | integer | DIIS history depth. Default `10`. |
| `ENGINE_DIIS_EPS` | real | Width of the ADIIS-to-DIIS blend, in terms of the DIIS error. Default `1.0`; `0.05` gives the earlier, slower hand-over (more iterations on metal-organic systems). |
| `ENGINE_NO_ADIIS=1` | flag | Use plain DIIS only (several times more iterations on hard cases). |

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
a valid Psi4 comparison) and were added the same way:

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

## COSX analytic Hessian and the CPSCF response (2026-09-28)

Defaults are the validated fast path. Reference (omeprazole HF/6-311G*
RIJCOSX Freq, 497 bf, 16 threads): 8:41 total, peak 5.28 GB, no scratch on
disk.
Everything below is optional.

**Memory opt-ins** (off by default; each trades time for memory):

| Variable | Effect | Measured on omeprazole |
|---|---|---|
| `ENGINE_COSX_HESS_LRFIT=1` | S^y / P^y (n^2 x 3N each) never stored: an extra grid pass for S^y.R | COSX phase -0.5..0.76 GB, +3.3 s |
| `ENGINE_COSX_HESS_PSUB=1` | main-pass post-processing (dK/dR.C, T2-T4, fit Z) per point sub-batch | main-pass peak 5.22 -> 4.84 GB, +7.3 s |
| `ENGINE_CPSCF_COSX_COOP=1` | response K: all threads on one grid block (one A / R / Hc) | 129-density call 1.96 -> 0.93 GB process peak, 21.0 -> 25.4 s |
| `ENGINE_FORCE_DF_DIRECT=1` | RI-J 3-centre store never built (also in the CPSCF) | CPSCF phase -1.5 GB, about +21 s |
| `ENGINE_COSX_HESS_KVCAP`, `ENGINE_COSX_HESS_ACAP` | doubles for the main pass's k^y/m^y and A buffers per thread (default 8e6 each) | with `ENGINE_CPSCF_COSX_RPTS=16` and all of the above: 3.44 GB, 10:50 |

The memory/time trade is step-like: the first items are cheap (about 5-20
s per GB), shrinking the per-thread batch buffers is not (about 66 s/GB),
and the global peak only drops when both the COSX pass and the CPSCF drop.

**Tuning knobs** (defaults in parentheses):

| Variable | Meaning |
|---|---|
| `ENGINE_COSX_HESS_KSCREEN` (1e-10), `..._KSCREEN2` (10x KSCREEN), `..._KSCREENK` (1e-8) | pair screens: Hessian terms, second-derivative integrals, dK/dR |
| `ENGINE_COSX_HESS_COLTHR` (1e-7) | b^y / a^y columns dropped from the dense A.R below this sqrt|w| max |
| `ENGINE_COSX_HESS_SNMEM` (1e9) | bytes over all threads for the per-thread S^y arrays (atom ranges) |
| `ENGINE_COSX_HESS_D2SEG`, `..._BLKPTS` (200), `..._SYMEM` (5.12e8, LRFIT) | d2-integral segment and grid chunk sizes (smaller chunks measured 20-25 s slower); LRFIT's atom-range budget |
| `ENGINE_CPSCF_COSX_GRID` (`lo`) | response-K grid: `lo` / `md` / `hi` |
| `ENGINE_CPSCF_COSX_RPTS` (64), `..._NBS` | response K: points per D_m x group; sub-batch (coop only) |
| `ENGINE_CPHF_CONV` (5e-5), `ENGINE_CPHF_SOLVER` (block CG; `pcg`, `diis`) | CPSCF convergence and solver (block CG since 2026-09-29: all perturbations share each step's search block; omeprazole applies RKS PBE 968 -> 390, UHF 2026 -> 656) |
| `ENGINE_CPHF_VERBOSE`, `ENGINE_CPHF_DIAGGAP` (0.05), `ENGINE_CPHF_PFLOOR` | solver log; UHF PCG: exact diagonal preconditioner below that gap, optional gap floor |

**Archived paths.** The alternative COSX Hessian paths (dense tangents,
full / mirror / two-leg dK/dR, whole-block k^y and b^y, full-n A, per-thread
x-leg, fused integral kernel, the dense and ZT7 Z^y forms, PSTREAM, the
shared S^y accumulator) were removed from the build on 2026-09-28 and are
kept, with the switch each used and why it was replaced, in
`src/integrals/cosx/legacy/` (last buildable at commit `ffdd90d1`). Their
variables (`ENGINE_COSX_HESS_BFLY`, `..._XKSHARED`, `..._KCOMP`, `..._ACMP`,
`..._APRUNE`, `..._NOA1`, `..._T3LR`, `..._KXX`, `..._LRZ`, `..._D2GATE`,
`..._FUSED`, `..._ZT7`, `..._SNUM`, `..._BCOMP`, `..._PSTREAM`, `..._DIAG`,
`..._SCRDIAG`, `..._KXDIAG`, `..._NOKX`, `ENGINE_COSX_HESS34_CHECK`,
`ENGINE_COSX_KX_MIRROR`, `ENGINE_COSX_KX_NONVAR`, `ENGINE_COSX_KX_FULL`)
are no longer read. Still read on the CPSCF response side:
`ENGINE_CPSCF_COSX_RSUB`, `..._KSHARED`, `..._NONVAR`, `..._MIRROR`,
`..._LOOP`.

**Diagnostics:** `ENGINE_COSX_HESS_VERBOSE=1` (phase timings, and
per-phase RSS/peak; note it resets VmHWM, so /usr/bin/time's peak is then
per-phase), `ENGINE_COSX_HESS_XCHK=1` (COSX vs exact fixed-density
exchange Hessian and dK/dR.C), `ENGINE_COSX_HESS_DUMP=1`, `ENGINE_COSX_BATCH_TIME=1`,
`ENGINE_COSX_BATCH_CHECK=1`, `ENGINE_DFJ_BATCH_VERBOSE=1`,
`ENGINE_SMAPS_DUMP=1` (smaps + malloc_info at the SCF end),
`ENGINE_DFJ_BENCH=nd` / `ENGINE_COSXK_BENCH=nd` (time the batched RI-J /
response K for nd densities after the SCF, then stop),
`ENGINE_FORCE_DF_STORE=1`.

Memory units (PROFILE lines appear only at `verbose = 2`): PROFILE lines print GiB, `/usr/bin/time -v` KB, the
`COSX HESS mem` lines decimal GB.

---

## Meta-GGA analytic Hessian, closed shell (2026-09-30)

The fixed-density term, dVxc/dR.C_occ, the CPSCF response kernel and its occ-occ
projections carry the tau channel (tau = 1/2 sum P grad phi.grad phi; kernels
v_tau, f_rt, f_st, f_tt). R2SCAN and hybrids (R2SCAN0) work for RHF; UKS
meta-GGA stops with a message. Ethylene 6-31G*, exact J/K, vs pyscf on the
identical dumped grid (`ENGINE_XC_GRID_DUMP`), off-diagonal atom blocks:
R2SCAN 8.3e-7, R2SCAN0 4.6e-6 (diagonal blocks are not comparable: pyscf's fixed
grid is not translation invariant, ours builds them from the sum rule).

Grid sensitivity: the tau terms make the Hessian much more grid-sensitive than
the energy. Water R2SCAN 6-311G*, sum_A H_AB (should be 0), response part:
theory-table REGULAR 1.0e-3, HIGH 2.4e-4, level 7 (232838 pts) 2.2e-5 (PBE at
REGULAR: 1.7e-5). At level 7 the analytic Hessian agrees with a central
difference of Direwolf's own analytic gradient on the same grid to 3.9e-5
(max element, |H| up to 0.67), frequencies to 0.22 cm-1.
Recommendation: for meta-GGA Freq set `ENGINE_XC_THEORY_TIER=HIGH`
for the whole job (SCF included, so energy, gradient and Hessian share one
grid). The default tier is unchanged.
## Range-separated hybrids, analytic Hessian (2026-09-30)

wB97X-D, CAM-B3LYP etc. (K_eff = frac K + RS_beta K_LR) work in the analytic
Hessian with K='cosx', RHF and UHF: the COSX Hessian is called a second time
with omega = RS_omega and the LR pieces are folded into the arrays that the
callers already scale by frac (direct term, dK/dR . C_occ), and the batched
response K (integrals_build_exchange_cosx_batch / _batch_c) adds
(RS_beta/frac) K_LR. Before this the Hessian silently dropped K_LR (water
wB97X-D 6-311G* frequencies 893 / 3643 / 3767 vs 1731 / 3746 / 3913 cm-1 from
differencing the analytic gradient). Now, RI-J + COSX, water, central
difference of the analytic gradient on locked COSX grids: max element
difference 1.02e-3 (B3LYP on the same setup: 1.01e-3, so this is the COSX
numerical floor, not K_LR), frequencies within 3 cm-1; water cation (UHF)
1.2e-3. Exact / RI-K with a range-separated functional stops with a message.
The COSX Hessian needs D = C C^T to 1e-10: use scf_conv = 'tight' for small
test molecules (the default convergence stops with "hess_main_pass: needs
dK/dR.C").

## Hessian findings and known gaps (2026-09-30)

- **VV10 (wB97M-V, wB97X-V, B97M-V) Hessian is not implemented.** The analytic
  Hessian runs without error but omits the second derivative of E_VV10, and now
  prints a WARNING. Methylamine, wB97M-V/def2-SVP, RIJCOSX: analytic vs finite
  difference of the analytic gradient differs by up to 13.9 cm-1 (mid-range
  modes 1300-1650); the no-VV10 wB97X-D noise baseline is 7.8 cm-1, so the VV10
  term is roughly 5-10 cm-1.
- **COSMO/CPCM electrostatic Hessian (gepol cavity, RHF and UHF)** is analytic:
  fixed-density second derivative of the stationary functional (charge response
  through the projected A^-1, tessera-centre and area Jacobians from the cavity
  derivative, second derivatives of the cavity by a central difference of those
  Jacobians, step ENGINE_COSMO_HESS_STEP = 2e-4 bohr), the derivative Fock matrix
  and the solvent kernel in the CPSCF. Methylamine, PBE/def2-SVP, water: analytic
  vs finite difference of the analytic COSMO gradient within 0.6 cm-1 (725.87 /
  892.10 / 1053.34 / 1106.69 / ... / 3633.66). It stops if the tessera count
  changes under the displacement (smaller step) or the cavity is not gepol.
  Neutral solutes only (same envelope-theorem limit as the COSMO force). The COSMO force
  also used to stop with "tessera count changed" at geometries where tesserae
  smaller than 1e-4 Ang^2 exist: the energy side drops them, the force side did
  not; both now drop them.
- **K='RI' (RI-K) has no analytic Hessian and stops with an explanation**; use
  K='cosx' or K='exact'. Range-separated hybrids (RHF and UHF) run with K='cosx'
  or with exact J + exact K: the erf(omega r)/r exchange enters the direct term
  (`integrals_hessian_2e_direct(..., omega=)`), the Fock derivative (four-index
  and density-role dK_LR) and the CPSCF response (dense K_LR per trial density,
  no batching yet, so the response is the slow part). wB97X/6-311G*, distorted
  water and its cation, on the dumped grid vs pyscf: off-diagonal atom blocks
  within 1.4e-6 (RKS) and 2.0e-6 (UKS). Exact K with RI-J and RSH is refused.
- **Classical terms** (D2/CHG/D3/D4/gCP/SMD-CDS) enter the Hessian as the central
  difference of their own analytic gradient; `-D` corrections need no further
  testing (D3BJ methylamine dimer: max 2.10 cm-1, grid noise).
- **Default-grid accuracy, meta-GGA**: Direwolf's default grid is not converged
  for R2SCAN (about +3 cm-1); the M06 family is grid sensitive in all codes
  (Direwolf-pyscf 1.4-3.6 cm-1 at defaults).
- **UKS meta-GGA**: the tau terms cost about 1.6-1.75x per grid point; the
  R2SCAN/PBE Hessian ratio of 2.3x on the caffeine cation is mostly the denser
  SCF/response grid the theory table picks for meta-GGA.

## DFT (LDA/GGA) analytic Hessian: XC terms and RI-J (2026-09-28)

Reference: omeprazole PBE/6-311G* RI-J Freq, 497 bf, 16 threads. Before:
489 s, peak 6.7 GB. Now: 196 s, peak 1.97 GB (the SCF's RI-J store; about
6 s more and 1.57 GB with `ENGINE_FORCE_DF_DIRECT=1`).
Frequencies are unchanged by all of the below.

Defaults (each old route kept behind a switch, all verified to 1e-11 on
distorted naphthalene):

| Piece | Default | Old route |
|---|---|---|
| fixed-density XC Hessian | blocked, shell-screened chunks, dgemm (68.8 -> 3.1 s) | `ENGINE_XC_HESS_PT=1` (also used with `ENGINE_XC_HESS_TI=0`) |
| dVxc/dR in the CPSCF right-hand side | dVxc/dR . C_occ in one blocked, value-screened pass, no n^2 x 3N per thread (58 -> 11 s, -4 GB) | `ENGINE_XC_VXCDERIV_FULL=1` |
| CPSCF XC kernel (fxc) | non-occupied index in the AO basis on the grid, value-screened AO rows, one dgemm per perturbation group per step (213 -> 71 s) | `ENGINE_XC_PROJ_GRP=1` (MO projection, grid-innermost), `..._CHK=1`, `..._PT=1` |
| RI-J dJ/dR | one dgemm per shell pair for all 3N (23 -> 5.5 s) | - |
| RI-J in the CPSCF | stays DIRECT after the store release, factored densities, one b- and one j-pass per build (129 densities 1.7 s) | `ENGINE_CPSCF_DF_STORE=1` (refill the store), `ENGINE_DFJ_BATCHC_DENSE=1` (16-density dense chunks) |
| pure functionals' CPSCF J | factored (J C_occ), no n^2 x nbatch arrays | `ENGINE_CPSCF_DENSE_J=1` |

CPSCF XC kernel screens (`xc_fxc_c`, also used by dVxc/dR . C_occ):
`ENGINE_XC_PROJ_VTHR` (1e-7) drops AO rows whose value / gradient (/ second
derivative) stays below it on a 128-point chunk (omeprazole 274 -> 154 rows,
naphthalene effect 2.7e-11). `ENGINE_XC_PROJ_LOC=1` (off) switches the kernel
to Pipek-Mezey localized occupied orbitals and drops, per chunk, those below
`ENGINE_XC_PROJ_OTHR` (default = VTHR): omeprazole keeps 90.9 / 85.0 / 69.4
of 91 at 1e-7 / 1e-5 / 1e-4 (CPSCF fxc 78 / 70 / 61 s, Hessian change 2e-9 /
1.7e-7), so it is off at this size.

Diagnostics: `ENGINE_DFJ_VERBOSE=1` (dJ/dR phase timings),
`ENGINE_DFJ_BATCH_VERBOSE=1` (per RI-J batch build), `ENGINE_DFJC_BENCH=nd`
(time the factored RI-J batch for nd densities after the SCF, then stop).
Open shell (UKS, LDA/GGA, since 2026-09-29): spin-polarized fixed-density
XC Hessian, dV^s_xc/dR . C^s_occ and the alpha/beta-coupled kernel in the UHF
CPSCF (on the same coarser CPSCF grid). Before this the UHF Hessian had no XC
terms at all. `ENGINE_HESSIAN_FORCE_UHF=1` runs a closed shell through the
UHF/UKS code (check: naphthalene PBE vs RKS 8.8e-10 with
`ENGINE_CPHF_CONV=1e-9`); `ENGINE_XC_GRID_DUMP=file` writes the SCF grid
(x y z w, bohr) for running pyscf on the identical quadrature (vinyl radical
UKS vs pyscf, off-diagonal atom blocks 1e-6). The UHF CPSCF runs factored like
the closed shell (RI-J with COSX or no exact exchange; PCG; no n^2 x
perturbation arrays); `ENGINE_CPSCF_UHF_DENSE=1` keeps the old dense DIIS
route (whose COSX form differs from RHF by up to 5e-5 on the `lo` grid).
Omeprazole cation (doublet) PBE/6-311G* RI-J Freq: 344 s with the block
CG (574 s with PCG; dense route > 20 min, 5.4 GB), 2.5 GB.
Omeprazole with the block CG: PBE 147 s, B3LYP
519 s.
VV10 Hessian is not implemented (see "Hessian findings and known gaps").

ECPs in the Hessian (since 2026-09-29): ECP second derivatives and the ECP
operator-centre first derivative were missing (CH3I def2-TZVPP vs pyscf
2.35e3, now 2e-5). Also fixed the same day: def2-SVP / TZVP / QZVP / TZVPD
silently had no ECP for Z >= 37 (only def2tzvpp.ecp exists in data/bases).

---

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

## Newton-CG SCF

Newton-CG orbital steps (finite-difference Fock Hessian-vector products, preconditioned CG) replace DIIS from a given SCF cycle onward, once the orbital gradient RMS is below `ENGINE_NEWTON_START` (default 5e-3).

- Gaussian-style route token `Newton-CG=N` (or TOML `newton_cg = N`, env `ENGINE_NEWTON_ITER`): first cycle that uses Newton-CG. Default 64; `0` disables completely.
- `ENGINE_NEWTON=1`: force Newton-CG from cycle 3.
- `ENGINE_NEWTON_CGTOL`: CG relative residual tolerance (default 0.1).
- `ENGINE_NEWTON_START`: orbital-gradient threshold below which a Newton step is taken (default 5e-3). Larger values start earlier and are much less stable.
- `ENGINE_NEWTON_H`: relative size of the finite-difference step (default 0.01 of the largest element of the trial density change).
- With COSX the trial Fock builds keep the SCF's incremental exchange state intact; `ENGINE_NEWTON_FDINC=0` uses a separate full rebuild instead.
- `ENGINE_NEWTON_ANALYTIC=1`: analytic Hessian-vector product (J, exact or COSX K and the XC kernel built on the density change alone) instead of finite differences; falls back to finite differences for range-separated hybrids, VV10, solvation and RI-K. `ENGINE_NEWTON_CHECK=1` prints its difference to the finite-difference product and the J/K/XC times for the first products of each step; `ENGINE_NEWTON_STOP=n` ends the run after n such products.

Newton-CG is an opt-in fallback: it lowers the number of cycles but each trial product costs about one Fock build, so the SCF takes longer than DIIS on the systems tested.

## Wavefunction stability

`ENGINE_STABILITY=1` (Gaussian-style route token `Stable`, TOML `stability = 1`): after SCF converges, find the lowest eigenvalue of the orbital Hessian (Lanczos on the finite-difference Fock response, so J/COSX-K/XC are all included) and report whether it is negative (< -1e-3). `ENGINE_STABILITY=2` (`Stable=Opt`, `stability = 2`): on an instability, rotate the orbitals along the lowest eigenvector (step norm 0.5) and re-converge, up to 3 rounds. Closed-shell (RKS/RHF) runs test restricted-to-restricted stability only; RHF->UHF instabilities are not tested. Costs up to 40 Fock builds per check.
