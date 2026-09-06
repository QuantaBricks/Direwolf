# EngineUp — the Direwolf library entry point

`EngineUp` (defined in `src/main/engine.f90`) runs one complete calculation:
initial guess → SCF → properties → optional analytic force. Every run-mode
choice is an explicit dummy argument — no environment variables, no global
config object. `src/harness/run_engine.f90` is just a thin file-driven front
end that parses a `&molecule` namelist into these arguments and calls
`EngineUp` once; a geometry optimizer / MD driver / QM-MM loop should call
`EngineUp` directly instead, mutating coordinates between steps.

---

## Calling it from Fortran

`EngineUp` is a plain external subroutine (not a module procedure), but it
has assumed-shape (`dens_in_*`) and allocatable-`intent(out)` (`dens_out_*`)
dummy arguments, so the Fortran standard **requires an explicit interface at
every call site**. Use the one that ships with the tree:

```fortran
use mod_engineup_interface, only: EngineUp
```

`mod_engineup_interface.f90` is a hand-maintained interface block that must
stay in sync with `engine.f90`'s actual declarations — if you change
`EngineUp`'s signature, change both.

Without that interface the `allocate` done inside `EngineUp` for
`dens_out_a/dens_out_b` silently fails to propagate back to the caller
(`gfortran -fcheck=all` catches it downstream as "Allocatable actual
argument … is not allocated").

---

## Argument reference

Order matches the actual signature. `n = ncenters`, `p = npc`.

### System

| arg | type | in/out | meaning |
|---|---|---|---|
| `ncenters` | `integer` | in | number of atoms |
| `imult` | `integer` | in | spin multiplicity 2S+1 (1 = closed-shell singlet) |
| `icharge` | `integer` | in | total molecular charge |
| `functional_in` | `character*30` | in | functional name, e.g. `'HF'`, `'B3LYP'`, `'WB97M-V'`, `'PBE_PBE'` — resolved by `mod_xc.f90:xc_select_functional` (`'HF'` is the real way to ask for plain Hartree-Fock) |
| `coord` | `real(8)(n,3)` | in | Cartesian coordinates, **Angstrom** |
| `atomchg` | `integer(n)` | in | nuclear charge (Z) per atom |
| `baselable` | `character*30` | in | global basis label — must match an `<Z>-<label>` block in `data/bases` (e.g. `'def2svp'`, `'6311gs'`) |
| `ecplabel` | `character*30` | in | global ECP label; `''` = same as `baselable` (no ECP unless the basis file carries one) |
| `atombasis` | `character(n)*30` | in | per-atom basis override; `''` in a slot = use `baselable` |
| `atomecp` | `character(n)*30` | in | per-atom ECP override; `''` = use `ecplabel` |

### Method / build mode

| arg | type | in/out | meaning |
|---|---|---|---|
| `j_mode` | `character(*)` | in | Coulomb build: `'exact'` or `'RI'` |
| `k_mode` | `character(*)` | in | exchange build: `'exact'`, `'RI'`, or `'cosx'` — independent of `j_mode`. `J='RI' K='cosx'` is RIJCOSX; `J='RI' K='RI'` is RIJK; `J='exact' K='exact'` is conventional. Case-insensitive. |
| `ri_aux_basis` | `character(*)` | in | auxiliary basis label for the RI path; `''` = auto-generated even-tempered aux. Overriding with e.g. `'ccpvqzjkfit'` reduces RI-J incompleteness error at higher cost. |
| `puream` | `logical` | in | spherical (`.true.`) vs Cartesian (`.false.`) Gaussians |
| `harris_guess` | `logical` | in | Harris/SAD-style initial guess (`.true.` recommended) |
| `do_force` | `logical` | in | also compute the analytic nuclear gradient |
| `vv10_nonself` | `logical` | in | evaluate VV10 non-local correlation non-self-consistently (post-SCF) instead of in every cycle |
| `mem_cap_gb` | `real(8)` | in | soft memory budget in GB; `0.0` = auto-detect. Drives the STORE-vs-DIRECT decision for the 2e / DF tensors. |
| `estimate_only` | `logical` | in | skip SCF and force entirely; return only `mem_grid_gb` / `mem_2e_gb` |
| `n_threads` | `integer` | in | OpenMP thread count; `0` = leave `OMP_NUM_THREADS` / the runtime default alone |
| `scf_conv_level` | `integer` | in | SCF convergence preset: `0` regular (dE 1e-6 / dRMS 1e-5), `1` fine (1e-8 / 1e-7), `2` tight (1e-9 / 1e-8). The `&molecule` namelist's `scf_conv = 'regular'/'fine'/'tight'` string maps onto this. |
| `basedir` | `character(*)` | in | directory to resolve `data/bases` etc. against; `''` = built-in default |

### Point charges (QM/MM)

| arg | type | in/out | meaning |
|---|---|---|---|
| `npc` | `integer` | in | number of external point charges (`0` = none) |
| `pc_charge` | `real(8)(p)` | in | point-charge values (e) |
| `pc_coord` | `real(8)(p,3)` | in | point-charge positions, **Angstrom** |

### Molden I/O

| arg | type | in/out | meaning |
|---|---|---|---|
| `molden_write` | `logical` | in | write a Molden file of the converged orbitals |
| `molden_file` | `character(*)` | in | path for the above |
| `molden_read` | `logical` | in | seed the guess from a Molden file |
| `molden_read_file` | `character(*)` | in | path for the above |

### COSMO / SMD solvation

All ignored unless `cosmo_on = .true.`.

| arg | type | in/out | meaning |
|---|---|---|---|
| `cosmo_on` | `logical` | in | enable implicit solvation |
| `cosmo_epsilon` | `real(8)` | in | solvent dielectric (ignored if `cosmo_solvent` is set) |
| `cosmo_radii_scale` | `real(8)` | in | Bondi-radius multiplier for the cavity |
| `cosmo_avg_area` | `real(8)` | in | target tessera area, Å² |
| `cosmo_sigma_rav` | `real(8)` | in | σ-profile averaging radius, Å (only with `cosmo_sigma_profile_file`) |
| `cosmo_cavity_type` | `character(*)` | in | `'gepol'` (analytic, has an analytic force) is the default and only type |
| `cosmo_rsolv` | `real(8)` | in | solvent probe radius |
| `cosmo_ks_nseg`, `cosmo_ks_nface` | `integer` | in | Klamt-surface segmentation controls |
| `cosmo_sigma_profile_file` | `character(*)` | in | write a σ-profile to this path; `''` = off |
| `cosmo_smd` | `logical` | in | add the SMD CDS (non-electrostatic) term; needs `cosmo_cavity_type='gepol'` |
| `cosmo_solvent` | `character(*)` | in | named solvent, e.g. `'water'`, `'methanol'` |

### Outputs

| arg | type | in/out | meaning |
|---|---|---|---|
| `force_out` | `real(8)(n,3)` | out | nuclear gradient, **Hartree/Bohr** (valid only if `do_force`) |
| `energy_out` | `real(8)` | out | total energy, **Hartree** |
| `MLcharge_out` | `real(8)(n)` | out | atomic charges (RESP if enabled, else the fallback population charge) |
| `iconv` | `integer` | out | `1` = SCF converged, `0` = not |
| `econv` | `real(8)` | out | final SCF convergence measure |
| `mem_grid_gb` | `real(8)` | out | predicted XC-grid-cache footprint, GB |
| `mem_2e_gb` | `real(8)` | out | predicted 2e / DF-tensor footprint, GB |

### Density-matrix warm start

| arg | type | in/out | meaning |
|---|---|---|---|
| `dens_in_a`, `dens_in_b` | `real(8)(:,:)` | in | α/β density guess; pass shape `(0,0)` for "no warm start" |
| `dens_out_a`, `dens_out_b` | `real(8)(:,:)` allocatable | out | converged α/β density — allocated by `EngineUp`, caller owns it afterward |

Feed a previous call's `dens_out_*` into the next call's `dens_in_*` (e.g. via
`move_alloc`) to warm-start a geometry step. `run_engine.f90` bridges this to
on-disk checkpoint files, but that is harness policy, not part of `EngineUp`.

---

## Units and conventions

- Input coordinates and point-charge positions: **Angstrom**. (`run_engine`'s
  `&molecule` parser converts a `unit='bohr'` input to Angstrom before the
  call.)
- Energy: **Hartree**. Gradient: **Hartree/Bohr**.
- `imult` is 2S+1 (so 1 for a closed-shell singlet, 2 for a doublet).
- `iconv == 1` is the only "results are trustworthy" signal — always check it.

---

## `estimate_only` mode

Set `estimate_only = .true.` to get `mem_grid_gb` / `mem_2e_gb` back without
running SCF or force. Everything else (`energy_out`, `force_out`, …) is
undefined in that mode. Useful for a driver that wants to pick `j_mode` /
`k_mode` / `mem_cap_gb` before committing to a full run.

---

## Minimal call

```fortran
use mod_engineup_interface, only: EngineUp
implicit none
integer :: iconv
real(8) :: e, econv, mem_g, mem_2e
real(8) :: coord(3,3), force(3,3), q(3)
integer :: Z(3)
character(30) :: ab(3), ae(3)
real(8), allocatable :: din_a(:,:), din_b(:,:), dout_a(:,:), dout_b(:,:)

Z = [8, 1, 1]
coord = reshape([0.d0,0.d0,0.d0,  0.7584d0,0.d0,0.5861d0,  -0.7584d0,0.d0,0.5861d0], [3,3], order=[2,1])
ab = ''  ;  ae = ''
allocate(din_a(0,0), din_b(0,0))

call EngineUp( 3, 1, 0, 'B3LYP', &
     coord, Z, 'def2svp', '', ab, ae, &
     'RI', 'cosx', '', .true., .true., .false., .false., 0.d0, &
     .false., 8, &
     0, [real(8)::], reshape([real(8)::],[0,3]), &
     .false., '', .false., '', &
     .false., 78.4d0, 1.2d0, 0.3d0, 0.5d0, &
     'gepol', 1.3d0, 0, 0, '', &
     .false., '', &
     force, e, q, iconv, econv, &
     mem_g, mem_2e, 0, '', &
     din_a, din_b, dout_a, dout_b )

if (iconv == 1) print *, 'E =', e
```

Link against `libengine.a` plus the third-party libs the top-level `Makefile`
already lists on the `Direwolf` link line (libcint, OpenBLAS, libxc,
simple-dftd3, gcp, dftd4, multicharge, toml-f).

---

## In-process integration (Python etc.)

`EngineUp`'s all-explicit-scalar/explicit-shape argument list (apart from the
four density arrays) is deliberately `iso_c_binding`-friendly. For QM/MM or
fragment methods, call it in a loop from the host process rather than
spawning `run_engine` per step — that avoids process startup + input
re-parsing on every energy/gradient evaluation. Provide your own thin C-
interoperable shim if calling from C/Python; keep the density arrays on the
Fortran side or marshal them explicitly.

---

## Caveats

- **One call at a time.** `EngineUp` uses module-level state and OpenMP
  regions internally; do not call it concurrently from multiple threads of
  the same process.
- A few knobs are still environment variables, read inside the run rather
  than passed here — notably `ENGINE_GRID_CACHE` (XC-grid STORE vs DIRECT)
  and the COSX grid-density overrides. See `docs/ENVIRONMENT.md`.
- `mod_engineup_interface.f90` is maintained by hand. A signature change that
  is not mirrored there produces wrong-looking runtime failures (bad
  descriptors, unallocated `dens_out_*`), not a compile error.
