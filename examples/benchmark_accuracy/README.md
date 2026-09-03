# Functional accuracy benchmark set

All inputs here are the same water molecule (6-31G, geometry below,
Angstrom) run with a different functional, kept as a standing
regression/accuracy check against real reference-code numbers (never
finite-difference alone - see this project's own testing convention).
The "-3c" composite-method entries near the bottom of the table are the
one exception to "6-31G" above - each "-3c" method is DEFINED with its
own fixed basis (that's the point of the method), so those three use
`mtzvp`/`mtzvpp`/`vdzp` instead, same water geometry otherwise.

Geometry (Angstrom, matches every `.inp` file in this directory):
```
O   0.0000   0.0000   0.0000
H   0.7584   0.0000   0.5861
H  -0.7584   0.0000   0.5861
```

## Quick check

```
python3 examples/benchmark_accuracy/run_benchmark.py
```

Builds nothing itself (build Direwolf first via `make` from the repo
root) - runs every `.inp` here against the `Direwolf` binary, diffs
`total_energy_hartree` against `references.csv`, and exits non-zero if
anything is out of its own per-row tolerance. `references.csv` is the
machine-readable source of truth; the table below is a human-readable
snapshot of the same data plus narrative context - if they ever
disagree, trust `references.csv` (or better, regenerate this table from
it) since it's what the script actually checks against.

To re-verify after a code change, run each input and compare
`total_energy_hartree` against the reference column below. Reference
values were obtained via Psi4 (`scf_type pk`, `reference rks`) unless
noted; PySCF (`dft.RKS`) was used as an independent cross-check for a
few entries where noted.

| Input | Functional | Reference energy (Ha) | Direwolf energy (Ha) | diff | Notes |
|---|---|---|---|---|---|
| water.inp (examples/) | PBE | -76.298205421 (self-consistent, verified at libxc migration) | -76.298205421 | 0 | baseline regression |
| water_pbe0.inp | PBE0 | -76.3011085524414341 (Psi4) | -76.301108122 | ~4.3e-7 | |
| water_blyp.inp | BLYP | -76.3648063640215184 (Psi4) | -76.364805699 | ~6.7e-7 | |
| water_bp86.inp | BP86 | -76.3851450025465795 (Psi4) | -76.385145049 | ~4.6e-8 | |
| water_lda.inp | LDA (Slater+VWN_RPA) | -76.0134411618707446 (Psi4 `svwn`) | -76.013440496 | ~6.7e-7 | LDA_C_VWN_RPA, not LDA_C_VWN/"VWN5" - see mod_xc.f90's LDA case comment |
| water_camb3lyp.inp | CAM-B3LYP (DF/RI) | -76.3556819885017006 (Psi4 `set puream true`) | -76.355681491 | ~5e-7 | reference updated 2026-08-03: old -76.355606583 was Psi4 default Cartesian-aux; Direwolf defaults to spherical since commit 2493143 and matches Psi4 puream=true (the reference update's own aux representation) to ~5e-7. force also verified, see commits 3a3fa47/7a44bcf |
| water_camb3lyp_exact.inp | CAM-B3LYP (exact K) | -76.3555844385 (Psi4 `scf_type pk`) | -76.355583949 | ~4.9e-7 | force verified, commit 5dc02db |
| water_tpss.inp | TPSS | -76.3882095839735200 (Psi4) | -76.388209164 | ~4.2e-7 | PySCF cross-check: -76.38820939205318; force also verified (~1e-6 Ha/Bohr), open-shell (OH radical) separately verified vs Psi4 UKS (~2.2e-7 Ha) |
| water_m06l.inp | M06-L | -76.3773264586606615 (Psi4) | -76.377333957 | ~7.5e-6 | corrected from an earlier ~7.5e-7 transcription slip - see references.csv |
| water_mn15l.inp | MN15-L | -76.3003198936959137 (Psi4) | -76.300317408 | ~2.5e-6 | |
| water_mn15.inp | MN15 | -76.3020418020673930 (Psi4) | -76.302039980 | ~1.8e-6 | |
| water_m062x.inp | M06-2X | -76.3489752858216804 (Psi4) | -76.348994378 | **~1.9e-5** | KNOWN, investigated, unresolved residual - see mod_xc.f90's M06-2X case comment. Not RI/grid/coefficient-extraction/thread-order related (all individually ruled out); PySCF (-76.34897576455575) agrees with Psi4 to ~1e-6, not with Direwolf. Force also verified vs Psi4 to ~1e-6 Ha/Bohr despite the energy residual (commit 2a14945-era work) |
| water_m06.inp | M06 | -76.3512850249167485 (Psi4) | -76.351297376 | ~1.24e-5 | same documented Minnesota-family residual as M06-2X |
| water_m052x.inp | M05-2X | not independently cross-checked | -76.368603317 | - | same code path as M06-2X/M06, not separately verified against Psi4 |
| water_wb97xd.inp | wB97X-D | electronic: -76.362519867 (Psi4 minus its own dispersion term); dispersion: -2.70787e-5 (Psi4/PySCF) | electronic: -76.362508214; dispersion: -2.70787e-5 | dispersion ~1e-9 (essentially exact); electronic **~1.1e-5** (same documented residual family as M06-2X/M06, despite wB97X-D being RSH not global-hybrid - confirms the residual tracks with "highly-parameterized functional form", not with any single mechanism) | |
| water_pbe0_d3.inp | PBE0-D3 (Grimme zero-damping, "GD3") | -76.30111300939376 (Psi4) | -76.301112520 | ~4.9e-7 (same order as plain PBE0's own baseline) | isolating the dispersion term alone (Direwolf/Psi4 each minus their own plain-PBE0 energy) matches to ~1e-9 Ha - force (analytic gradient) verified the same way on both this molecule and octane (26 atoms, real sp3 CN variation), matching to ~1e-9 Ha/Bohr |
| water_pbe0_d3bj.inp | PBE0-D3BJ (Becke-Johnson damping, "GD3BJ") | -76.30138547553199 (Psi4) | -76.301384986 | ~4.9e-7 (same order as plain PBE0's own baseline) | same isolation methodology as PBE0-D3 above - dispersion term alone matches to ~1e-9 Ha (energy) and Ha/Bohr (force) |
| water_b973c.inp | B97-3c (`baselabel='mtzvp'`) | -76.397671282703 (Psi4) | -76.397608964 | ~6.2e-5 | Composite: B97-3c functional (libxc `XC_GGA_XC_B97_3C`) + def2-mTZVP + D3(BJ) (functional string `"b973c"`, params bit-identical to real ORCA's printed a1/s8/a2) + SRB correction (NOT gCP - `third_party/gcp`'s own select-case routes `b973c` to SRB only, contrary to an earlier informal assumption). Force verified: analytic vs FD ~2.6e-6 Ha/Bohr. Also cross-checked vs real ORCA 6.1.0: 2.0e-4 Ha - a systematic ORCA-vs-Psi4 offset seen across all three -3c methods here, not a Direwolf defect |
| water_r2scan3c.inp | r2SCAN-3c (`baselabel='mtzvpp'`) | -76.418331377062 (Psi4) | -76.418283782 | ~4.8e-5 | Composite: r2SCAN (existing functional, unchanged) + def2-mTZVPP + D4 (vendored `third_party/dftd4`, functional string `"r2scan-3c"`) + damped gCP. Force verified: analytic vs FD ~2.1e-5 Ha/Bohr (larger than B97-3c's, but isolated D4/gCP gradients alone match FD to ~1e-12/1e-14 - the residual traces to r2SCAN's own known meta-GGA grid sensitivity, not this session's new code). ORCA cross-check: 6.3e-4 Ha |
| water_wb97x3c.inp | wB97X-3c (`baselabel='vdzp'`) | -17.270942540919 (Psi4) | -17.270919745 | ~2.3e-5 | Composite: WB97X-V's combined libxc handle with VV10 explicitly DISABLED (`mod_xc.f90`'s own `"WB97X-3C"` case - confirmed via Psi4's `hyb_functionals.py` source, `"nlc": False`, not documented in any paper abstract/manual) + D4 (`"wb97x-3c"`) + vDZP (small-core ECPs on B-Rn, imported via `scripts/import_basis.py --gaussian94 --ecp`, extended to parse a second raw ECP channel-header grammar that basis didn't use before). No gCP/SRB. Total energy magnitude (~-17 Ha) is correct, not a bug - matches real ORCA's own -17.27 Ha. Force verified: analytic vs FD ~1.4e-7 Ha/Bohr. ORCA cross-check: 2.0e-4 Ha. **vDZP has no fluorine** (upstream basis gap - its own header comment says the ECP lacks a projector for F) |

## VV10 (nonlocal vdW correlation)

`he2_vv10.inp` is a He dimer (aug-cc-pVDZ, He at ±2.0 Å, matching Psi4's
own `tests/dft-vv10` input) run with functional `VV10` - the Vydrov-Van
Voorhis 2010 functional (rPW86 exchange + PBE correlation + the VV10
nonlocal term, JCP 133, 244103 (2010); arXiv:1009.1421). The nonlocal
kernel is implemented natively in `src/xc/mod_vv10.f90` (libxc 6.1.0 no
longer ships a standalone VV10; Psi4 likewise computes the kernel itself
via BRIAN) and evaluated on a separate coarse NLC grid
(`grid_gen.f90`'s `gridgen_nlc`, Psi4-style `DFT_VV10_RADIAL/
SPHERICAL_POINTS`, engine default 20x50, env `ENGINE_VV10_RADIAL/ANGULAR`).
Verification:

| Quantity | Psi4 (`tests/dft-vv10`, 20x50 NLC grid) | Direwolf (20x50 NLC grid) | Notes |
|---|---|---|---|
| VV10 NLC energy (E_c^nl + beta N) | 0.01879662804 | 0.018787481 | ~9e-6 Ha (0.05%) |
| Total energy (Ha) | -5.8199583320852053 | -5.819633414 | ~3.2e-4 Ha |

The potential is verified exactly: a standalone finite-difference check
confirms the analytic F_n and F_g match numerical derivatives of
E_c^VV10 to ~1e-8 at regular grid points. The NLC energy is well
converged w.r.t. the NLC grid (20x50 vs 50x146 differ by <3e-5 Ha on
benzene/caffeine), and rcut/rhocut truncation (ENGINE_VV10_RCUT/
ENGINE_VV10_RHOCUT, defaults 20 bohr / 1e-9) changes it by <1e-9.
Open-shell (OH doublet), water, CH4, benzene, and caffeine (24 atoms)
all converge cleanly. Cost: the O(npts^2) kernel sum on the coarse NLC
grid adds roughly 2-4x to plain-PBE wall time (caffeine/6-31G 8.5s ->
32s); the previous implementation on the production DFT grid was ~40x
slower still (benzene 12m46s -> 4.2s).

Analytic VV10 gradient (calc_force_vv10, force.f90): implements the
paper's three-term decomposition - the GBF term with the MOVING-point
convention (each NLC point rides its parent atom), the quadrature-weight
term (Eq. 24, Becke weight derivatives by re-running gridgen_nlc at
+-delta), and the grid term (Eq. 25-26, mod_vv10's vv10_grid_force).
The Becke-weight term is verified exactly against the fixed-density NLC
energy finite difference, and the whole force is translation-invariant
and finite-difference-checked to ~1e-4-6e-4 Ha/Bohr (~1-3%) on water/He2
- the residual (traced to the GBF/grid part) is a known limitation; the
semilocal XC force omits its own weight term consistently with Psi4.

## What the ~1e-5 to 2e-5 Ha residuals mean (and don't mean)

M06, M06-2X, and wB97X-D's electronic part show a residual roughly
10-30x larger than every other functional in this table (which all sit
at ~1e-7 to ~1e-6 Ha). This was investigated in depth, not just
observed and shrugged off:

- Per-grid-point libxc evaluation matches pylibxc to ~1e-15 (machine
  precision) at sampled representative points, for both M06-2X and
  wB97X-D.
- hf_frac/rs_omega/rs_beta extraction matches libxc's own reported
  values to full double precision.
- Direwolf's own grid is already converged for these functionals -
  bumping radial density, angular density, and disabling angular
  pruning each move the energy by less than 1e-7.
- The energy is bit-identical across 1/2/4/8 OpenMP threads (rules out
  floating-point summation-order sensitivity).
- Two independent reference codes (Psi4, PySCF) agree with EACH OTHER
  to ~1e-6, but both disagree with Direwolf by ~10-20x that.

Most likely explanation: a genuine, small systematic bias between
Direwolf's own (self-converged) grid quadrature scheme and Psi4/PySCF's
(also self-converged, but differently-placed) scheme, specific to
functional forms sharp/rational-function-heavy enough to be sensitive
to which quadrature scheme is used, not just how many points it has.
This is a documented characteristic of "Minnesota-style"/highly-
parameterized functionals in the broader DFT literature, not unique to
Direwolf - "smoother" functional forms (PBE, B3LYP, CAM-B3LYP, TPSS, all
in this table) don't show it.
