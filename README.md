# Direwolf

An independent quantum chemistry (QM) engine providing QM calculation
results for other software, built for speed.

Developer: Xin Chen

**Status**: under active development since 2025 (AI-assisted); behavior
and interfaces may still change between commits. See
[docs/DEVLOG.md](docs/DEVLOG.md) for the project history and what's
shipped vs. not started yet.

## License

Copyright © 2025-2026 QuantaBricks. Developed by Xin Chen.

Direwolf is licensed under the **GNU Affero General Public License v3.0
or later** (`AGPL-3.0-or-later`) - see [LICENSE](LICENSE). Note AGPL
§13: running a modified version to provide a service over a network
obliges you to offer that version's complete source to its users.

For use that AGPL-3.0 does not permit - e.g. embedding Direwolf in a
closed-source product or service - a separate commercial license is
available from QuantaBricks.

Bundled third-party dependencies keep their own licenses (table at the
end of this file).

## Build

```
make
```
(requires gfortran, gcc, and cmake. Produces `./Direwolf` and
`lib/libengine.a`; `build/` holds `.o`/`.mod` files.)

Every dependency - libcint, libxc, OpenBLAS/LAPACK, and the `-D3`/`-D3BJ`
dispersion chain (simple-dftd3 + mctc-lib + toml-f) - is vendored
SOURCE under `third_party/` and built automatically by this same `make`
call the first time (each has its own Makefile rule; nothing is
prebuilt or checked in as a binary). Expect the first `make` to take a
few minutes (OpenBLAS alone runs its own test suite as part of its
build); subsequent builds skip anything already built. See "Third-party
dependencies" below for what's in `third_party/`, each one's license,
and how it's linked.

## Run

```
./Direwolf <input.inp> [output.out]
```
(output defaults to `<input.inp>.out` if omitted. Single-threaded unless
`OMP_NUM_THREADS` is set in the environment or `n_threads` is set in the
input file - see below.)

## Docker

```
docker build -t engine .
docker run --rm -v "$PWD":/data engine /data/input.inp /data/output.out
```

Multi-stage build: the first stage compiles Direwolf and every vendored
`third_party/` dependency from source (same as a plain `make`); the
final runtime image keeps only the `Direwolf` binary, `data/` (basis
sets - resolved relative to Direwolf's own path, so it must ship
alongside the binary), the 3 dynamically-linked LGPL `.so`s
(simple-dftd3/gcp/dftd4), and the glibc/libgfortran/libgomp runtime
libraries - about 130MB. `Dockerfile`/`.dockerignore` are tree-generic
(no path is specific to this repo vs. an exported `release/Direwolf-<ver>/`
tree - see [docs/RELEASE_EXPORT.md](docs/RELEASE_EXPORT.md)), so the
same `Dockerfile` builds either one unmodified.

Mount input and output under the **same** directory: Direwolf resolves
the output path's basename in its own current working directory
(`/app` inside the container), not the directory component of the path
you passed - see `run_engine.f90`. Put both files under one `-v` mount
(as in the example above) rather than under two different host
directories, or the output will land in `/app` inside the container
instead of where you expected.

## Supported functionals and methods

Set via `functional` in `&molecule` (case-insensitive). Any unrecognized
name falls back to plain PBE. See `docs/VALIDATION.md` for cross-code
verification numbers.

- **Hartree-Fock**: `HF`
- **LDA**: `LDA` (Slater + VWN-RPA, = Psi4 "SVWN")
- **GGA**: `BLYP`, `BP86`, `PBE`
- **Hybrid GGA**: `PBE0`, `B3LYP`
- **Range-separated hybrid**: `CAM-B3LYP`, `WB97X`, `WB97X-D`
- **meta-GGA**: `TPSS`, `R2SCAN`, `M06-L`, `MN15-L`
- **Hybrid meta-GGA**: `R2SCAN0`, `M06`, `M06-2X`, `M05-2X`, `MN15`
- **VV10 nonlocal correlation**: `VV10`, `B97M-V`, `WB97X-V`, `WB97M-V`
- **"-3c" composites** (basis + dispersion fixed by the method name, see
  `docs/3C_METHODS.md`): `B97-3c`, `R2SCAN-3c`, `WB97X-3c`

Append `-D2`, `-D3`, or `-D3BJ` to any name above for Grimme empirical
dispersion, e.g. `functional = 'B3LYP-D3BJ'`. The "-3c" methods are the
exception - their dispersion and basis are fixed, not user-selectable.

Other methods, documented in their own sections below:
- **Density fitting (RI-J/RI-K) and COSX**: `J`/`K`/`ri_aux_basis`, any functional above.
- **ECP** (effective core potentials): see "ECP" section below.
- **Implicit solvation** (COSMO/CPCM + SMD): see "Solvation" section below.
- **QM/MM point charges**: `&pointcharges` namelist (`npc`, `pc_q`, `pc_x`, `pc_y`, `pc_z`).
- **Analytic gradients**: `calc_force = .true.` (default), for every functional/method combination above.

## Input file format


Fortran namelist; see `examples/friendly_format/`, `examples/drug_molecules/`,
`examples/benchmark_accuracy/`, `examples/solvation/`, and `examples/misc/`
for working examples spanning water up to a ~100-basis-function alkane -
`examples/README.md` explains how the directories are grouped. The
namelist terminator can be either `/` or `&end`.

```
&molecule
 ncenters   = <number of atoms>
 imult      = <spin multiplicity, 2S+1 - 1 for closed-shell singlet>
 icharge    = <total molecular charge>
 functional = 'HF'                 plain Hartree-Fock, no XC functional
              'B3LYP' / 'B3LYP_HYB' B3LYP hybrid GGA
              'WB97M-V'             range-separated hybrid meta-GGA + VV10
                                    (see docs/ENVIRONMENT.md)
              anything else (e.g. 'PBE_PBE') falls through to a plain
              PBE GGA (PBE exchange + PBE correlation, no HF exchange)
 baselabel  = '6-31g'              basis set label - must match a
              "<atomic-number>-<baselabel>" block in data/bases
 J          = 'exact' / 'RI'       [optional, default 'exact'] Coulomb
              build mode, for BOTH the SCF energy and the analytic
              force/gradient
 K          = 'exact' / 'RI' / 'cosx'  [optional, default 'exact']
              exchange build mode, independent of J - the real RIJCOSX
              recipe is J='RI', K='cosx', not a package deal
 ri_aux_basis = 'ccpvdzjkfit'      [optional, default '' = auto-generated
              even-tempered auxiliary basis] set to load a real published
              fitting basis from data/bases instead - only cc-pVDZ-JKFIT
              (H/C/N/O/S) is bundled today
 calc_force = .true. / .false.     [optional, default .true.]
 spherical  = .true. / .false.     [optional, default .true.] spherical
              (puream) vs Cartesian Gaussians - see "Basis sets" below
 mem_cap_gb = <real>               [optional, default 0.0 = auto memory
              budget, see mod_meminfo.f90]
 estimate_only = .true. / .false.  [optional, default .false.] skip SCF/
              force entirely, just report the mem_grid_gb/mem_2e_gb
              memory estimates below
 n_threads  = <integer>            [optional, default 0 = leave
              OMP_NUM_THREADS/the OpenMP default thread count alone]
 resp_charges_on = .true./.false.  [optional, default .false.] compute
              two-stage RESP (restrained ESP-fit) atomic charges after
              the SCF converges - see "RESP charges" below
&end
&atoms
 O    0.000000    0.000000    0.000000
 H    0.758400    0.000000    0.586100
 H   -0.758400    0.000000    0.586100
&end
```

(`&atoms` column format, recommended: one atom per line, `element x y z`,
coordinates in Angstrom, no commas/quotes. `ncenters` is inferred from
the line count if left unset. `examples/run_engine.f90`'s own header
comment documents two more input styles: the original namelist array
form (`atomchg`/`x`/`y`/`z` as comma-separated arrays, still
backward-compatible) and an external `xyzfile`. It also documents the
optional `&professional` namelist for advanced performance tuning - see
docs/ENVIRONMENT.md.)


## ECP (effective core potential) and per-atom basis/ECP override

```
 ecplabel   = <ECP name>  [optional, default '' = same as baselabel]
              data/bases lookup name for the ECP block, independent of
              baselabel - pairs one valence basis with a differently-
              named core potential. An atom with no
              matching "<Z>-<ecp-name>.ecp" block is just all-electron.
```

Set `ecplabel` explicitly even when it equals `baselabel` (e.g.
`baselabel='lanl2dz'` + `ecplabel='lanl2dz'`) - makes the pairing
visible instead of relying on the default fallback, and is required if
you want a *different* valence basis on top of the same core potential.
See `examples/friendly_format/pd_hydrate_lanl2dz.inp`.

Per-atom override (`&atoms` column format only): add up to two optional
trailing columns to any atom line to mix basis families within one
molecule - `element x y z [basis] [ecp]`. Either column may be `-` to
skip it (fall back to the global default for just that slot); a line
with no extra columns behaves exactly as before.

```
&molecule
 baselabel = '6-31g'
&end
&atoms
 Fe   0.0000    0.0000    0.0000   lanl2dz   lanl2dz   ! both overridden
 O    2.0000    0.0000    0.0000   -         lanl2dz   ! ecp only
 H    2.5877    0.0000    0.7591                       ! global 6-31g, no ECP
&end
```

Resolution order per atom: basis = per-atom column, else global
baselabel. ECP = per-atom column, else global ecplabel, else (if
neither given) this atom's own resolved basis name.

## Basis sets

Bundled in `data/bases/` (imported via `scripts/import_basis.py` from
basissetexchange.org - see that script's header to add more; nothing
is auto-generated, so an unlisted element/basis combo errors clearly
at startup instead of silently running wrong). Split into one file per
basis-set family (`631g.bas`, `def2.bas`, `mtzvp.bas`, `sdd.bas`, ...) -
`mod_basis_files.f90` resolves which file holds a given label at lookup
time, so mixing basis families across atoms in the same molecule works
transparently. Elements 1-18 (H-Ar) unless noted otherwise.

Conventional basis-set names work directly: `baselabel`/`ecplabel`
accept either the compact label (left column below) or the common
keyword (right column), case-insensitively - e.g. `'6-31G*'` and
`'631gs'` are the same basis (translation table: `resolve_basis_name`,
`engine.f90`).

Direwolf supports both **spherical (puream) and Cartesian** Gaussians -
`spherical = .true. / .false.` in the `&molecule` namelist (default
`.true.`). Spherical is what conventionally-spherical basis families
(Dunning cc-pVxZ/aug-cc-pVxZ, Ahlrichs def2-*) actually mean by their
published numbers, and matches Psi4's own default (`puream true`); Pople
sets (6-31G*, 6-311G** etc.) are conventionally Cartesian either way, so
`spherical` doesn't move their energy. Covers energy, analytic gradients,
and f/g shells - see `docs/VALIDATION.md` for the verification numbers.
Cartesian mode (`spherical = .false.`) is still available and is what
Direwolf originally shipped with; it differs from the spherical/published
number by ~0.01-0.5 mHartree for d-and-higher conventionally-spherical
families, same magnitude either direction depending which one you're
diffing against a Psi4 default run.

### Orbital (valence) bases

Use as `baselabel` (compact label | common keyword aliases that
resolve to it):

| Compact label | Common keyword aliases | Notes |
|---|---|---|
| `6-31g` | 6-31G | H-Kr (Z=1-36) |
| `631gs` | 6-31G*, 6-31G(d) | |
| `631gdp` | 6-31G**, 6-31G(d,p) | |
| `631pgs` | 6-31+G*, 6-31+G(d) | |
| `631pgss` | 6-31+G**, 6-31+G(d,p) | |
| `6311g` | 6-311G | |
| `6311gs` | 6-311G*, 6-311G(d) | H-Ca (Z=1-20) + Ga-Kr (Z=31-36) - skips Sc-Zn (older entry) |
| `6311gdp` | 6-311G**, 6-311G(d,p) | |
| `6311pgss` | 6-311+G**, 6-311+G(d,p) | |
| `6311ppgss` | 6-311++G**, 6-311++G(d,p) | no He - not defined upstream |
| `321g` | 3-21G | |
| `ccpvdz` | cc-pVDZ | |
| `ccpvtz` | cc-pVTZ | |
| `ccpvqz` | cc-pVQZ | |
| `augccpvdz` | aug-cc-pVDZ | |
| `augccpvtz` | aug-cc-pVTZ | |
| `def2svp` | def2-SVP, SVP | |
| `def2tzvp` | def2-TZVP, TZVP | |
| `def2tzvpp` | def2-TZVPP, TZVPP | H-Rn (Z=1-86); all-electron through Kr(36), def2-ECP-paired from Rb(37) on - set `ecplabel='def2tzvpp'` explicitly for the ECP range |
| `def2tzvpd` | def2-TZVPD | |
| `def2qzvp` | def2-QZVP, QZVP | |

Bare `SVP`/`TZVP`/`TZVPP`/`QZVP` (no `def2-` prefix) alias to the def2
entries above. This is NOT guaranteed bit-identical to the original
1992/1994 Ahlrichs SVP/TZVP parameterization (a distinct, older dataset
BSE keeps separately as "Ahlrichs VDZ/pVDZ/TZV" - not imported here); ask
if exact reproduction of that specific legacy data is ever needed.

### ECP (core potential) families

Use as `ecplabel` (or `baselabel`, in the default coupled case):

| Compact label | Common keyword | Coverage |
|---|---|---|
| `lanl2dz` | LANL2DZ | H, Li-Bi (skips lanthanides Ce-Lu) plus U, Np, Pu - 71 elements, 62 ECP-bearing |
| `sdd` | SDD | K-Lr (Stuttgart RSC 1997) - complements lanl2dz's lanthanide gap - 63 elements, all ECP-bearing |

Import verification numbers: see `docs/VALIDATION.md`.

### DF auxiliary (fitting) basis

Use as `ri_aux_basis`:

| Compact label | Common keyword | Coverage |
|---|---|---|
| `ccpvdzjkfit` | cc-pVDZ-JKFIT | H, B-Ne, Al-Ar, Ga-Kr (Z=1,5-10,13-18,31-36) |
| `ccpvtzjkfit` | cc-pVTZ-JKFIT | same element set as ccpvdzjkfit (BSE publishes no other elements for this family) |
| `ccpvqzjkfit` | cc-pVQZ-JKFIT | same element set as ccpvdzjkfit |
| `def2universaljkfit` | def2-universal-JKFIT | H-Rn (Z=1-86); general-purpose RI-J/RI-K aux basis for ANY orbital basis, not just def2 - BSE's own metadata names this the recommended `jkfit` pairing for the WHOLE def2 orbital-basis family (SVP through QZVPP), so it's the right `ri_aux_basis` choice for `def2tzvpp`/`def2tzvp`/`def2svp` above, not just a fallback |

**Not bundled yet**: sdd for main-group elements past Ar, aug-cc-pVQZ,
def2-QZVPP, other ECP families (CRENBL/CRENBS, LANL08) - `def2tzvpp`'s
own def2-ECP (Rb-Rn) IS bundled now, see its row above. Elements past
Rn (87+) exist in
data/bases but aren't reachable by symbol through `run_engine.f90`
(only direct `atomchg` via in-process `EngineUp` can address them).

To add more: pull from Basis Set Exchange (NWChem format, or the
`--gaussian94`/`--ecp`/`--ecp-only --nwchem-ecp` flags of
`scripts/import_basis.py` for an ECP), run the import script, then add
a `resolve_basis_name` (`engine.f90`) case for a common-keyword alias.

## Solvation (implicit continuum: COSMO/CPCM electrostatics + SMD)

See `docs/ARCHITECTURE.md`'s `src/solvent/` entry for the module
layout. `&cosmo` namelist (all optional; `cosmo_on=.false.` is the
default vacuum behavior):

```
 cosmo_on          = .true./.false.  [default .false.]
 cosmo_solvent     = '<name>'        [default ''] e.g. 'water', 'methanol',
                     'dmso', 'toluene', 'hexane' - overrides cosmo_epsilon
 cosmo_epsilon     = <real>          [default 78.4] ignored if cosmo_solvent set
 cosmo_cavity_type = 'gepol'         [default] analytic, has an analytic force -
                     use this for energy+force
                   | 'ks1993'        COSMO-RS sigma-profiles (energy-only)
                   | 'yk1999'        smooth polynomial switching (energy-only)
                   | 'iswig'         Lange-Herbert erf switching (energy-only)
 cosmo_radii_scale = <real>          [default 1.2] Bondi-radius multiplier
 cosmo_avg_area    = <real>          [default 0.3, Ang^2] gepol tessera area
 cosmo_sigma_rav   = <real>          [default 0.5, Ang] only if cosmo_sigma_profile_file set
 cosmo_sigma_profile_file = '<path>' [default '' = off]
 cosmo_smd         = .true./.false.  [default .false.] needs cosmo_cavity_type='gepol'
                     + non-empty cosmo_solvent
```

Verified against Psi4 (CPCM) and pyscf (SMD's CDS term, bit-exact) -
see `docs/VALIDATION.md`.

### Example: water in implicit water (COSMO)

```fortran
&molecule
 imult      = 1
 icharge    = 0
 functional = 'B3LYP'
 baselabel  = '6-31g'
&end
&atoms
O           0.0000         0.0000         0.0000
H           0.7584         0.0000         0.5861
H          -0.7584         0.0000         0.5861
&end
&cosmo
 cosmo_on          = .true.
 cosmo_cavity_type = 'gepol'
 cosmo_solvent     = 'water'
&end
```

`cosmo_solvent = 'water'` sets `cosmo_epsilon` to water's dielectric
constant automatically (no need to also set `cosmo_epsilon` yourself);
`cosmo_cavity_type = 'gepol'` is the default and the only cavity type
with an analytic force, so it's the right choice unless you're doing an
energy-only COSMO-RS-style sigma-profile run (`'ks1993'`/`'yk1999'`/
`'iswig'` above). See `examples/solvation/water_sigma_profile.inp` for
that variant.

### Example: water in implicit toluene, with the SMD CDS correction

```fortran
&molecule
 imult      = 1
 icharge    = 0
 functional = 'B3LYP'
 baselabel  = '6-31g'
&end
&atoms
O           0.0000         0.0000         0.0000
H           0.7584         0.0000         0.5861
H          -0.7584         0.0000         0.5861
&end
&cosmo
 cosmo_on          = .true.
 cosmo_cavity_type = 'gepol'
 cosmo_solvent     = 'toluene'
 cosmo_smd         = .true.
&end
```

`cosmo_smd = .true.` adds the SMD cavity-dispersion-solvent-structure
(CDS) term on top of the GePol/CPCM electrostatics above - it requires
`cosmo_cavity_type = 'gepol'` and a non-empty `cosmo_solvent` (SMD's
CDS parameterization is per-solvent, unlike the electrostatic part
which only needs a dielectric constant). See `examples/solvation/water_smd_toluene.inp`.

## RESP charges

On by default (`resp_charges_on = .true.`) for any converged run - the
per-atom table in the output gains a third `RESP_q` column alongside
Mulliken/Lowdin, no separate table. Set `resp_charges_on = .false.` in
`&molecule` to skip the fit if you don't need it (see "Cost" below - it's
now cheap enough that this is mainly about a tidier two-column table, not
saved time). The implementation follows Psi4's RESP plugin
([cdsgroup/resp](https://github.com/cdsgroup/resp), BSD-3, which is what
Psi4 users actually run - Psi4 core itself ships no RESP) and through it
[Bayly:93:10269], matching its defaults exactly:

- Merz-Kollman grid, 4 shells at 1.4/1.6/1.8/2.0x VDW radius, 1.0 point
  per Ų, **GAMESS** VDW radii (not Bondi), GAMESS latitude-band
  unit-sphere sampling, and per-shell exclusion (a point on the 2.0x
  shell is rejected against other atoms' 2.0x radii).
- Stage 1: hyperbolic restraint `a=0.0005`, `b=0.1`, `IHFREE` (hydrogens
  are not restrained), iterated to `toler=1e-5` from the unrestrained ESP
  solution.
- Stage 2: `a=0.001`, sp3 carbons carrying at least one hydrogen are
  re-fit together with their own hydrogens (those hydrogens constrained
  equal to each other); every other atom is pinned to its stage-1 charge.

Verified on ibuprofen (33 atoms, B3LYP/def2-SVP) against Psi4 + the
plugin at the same level of theory, separating the three things a
final-charge comparison alone cannot tell apart (`ENGINE_RESP_DUMP`
writes the grid and ESP out for exactly this):

Both codes at `J='RI'`/`K='RI'` with `def2universaljkfit`, matching Psi4's
`scf_type df` (which already defaults to def2-universal-JKFIT here -
naming it explicitly changes Psi4's charges by exactly zero):

| | difference |
|---|---|
| Grid: point count, and each point's position | 1328 = 1328, max 5e-9 Å |
| Total SCF energy | 3.5e-5 Ha |
| ESP at those grid points | RMS 3.3e-6 a.u., 2.1e-4 relative |
| Fitting algorithm, both fitters run on Direwolf's own ESP | max 4.8e-7 e |
| End to end, each code on its own density | RMS 4.5e-5 e, max 1.2e-4 e |

Two independent checks say the residual is the SCF density and not the
fit. Feeding Direwolf's ESP through Psi4's own fitter reproduces the
end-to-end difference to every printed digit. And switching Direwolf from
RIJCOSX to matched RI-JK shrinks the energy gap 16.9x (5.9e-4 → 3.5e-5
Ha) and the charge gap 7.6x (9.1e-4 → 1.2e-4 e) together - the charges
track the density, as they must. The 4.8e-7 e residual in the algorithm
row is at the level the dump file's precision and the ESP fit's own
conditioning can resolve. See `examples/drug_molecules/ibuprofen_resp.inp`.

Repeated on three more drug molecules at the same matched RI-JK level, to
check the ibuprofen number wasn't a lucky cancellation:

| molecule | atoms | grid pts | max diff | RMS diff |
|---|---|---|---|---|
| paracetamol | 20 | 957 | 8.1e-5 e | 3.0e-5 e |
| caffeine | 24 | 1108 | 1.5e-4 e | 6.0e-5 e |
| fluoxetine | 40 | 1651 | 1.1e-4 e | 3.9e-5 e |
| ibuprofen | 33 | 1328 | 1.2e-4 e | 4.5e-5 e |

All four land in the same 1e-4 e band with no size trend, which is what
"the residual is SCF-density noise" predicts and a real algorithmic bug
would not.

**Cost.** The ESP-at-a-grid-point kernel (`resp_grid_and_esp`) originally
built a full nConts×nConts matrix per point purely to reduce it to one
scalar - the same integral libcint's `cint1e_grids` (COSX already used
it) evaluates for a whole block of points per call, with the shell-pair
loop outside the point loop and nothing of size nConts² ever
materialized. Same integrals, reordered: verified bit-identical against
the original point-at-a-time path (`ENGINE_RESP_ESP_MODE=legacy` keeps it
around as that reference) on both the ESP values (1e-12 a.u., the dump
file's own precision) and the final charges (max diff 0.0e+00 e).

| | `resp_grid_and_esp` before | after | speedup |
|---|---|---|---|
| ibuprofen, 33 atoms, 300 bf, 1328 pts | 2.76 s | 0.21 s | 13.1x |
| imatinib, 68 atoms, 673 bf, 2539 pts | 33.61 s | 0.60 s | 56.2x |

At 48 threads this puts the whole property at roughly 0.2-0.6 s even on a
70-atom drug molecule - under 1% of a typical job's wall time, down from
the ~12% it cost before this kernel existed. `resp_two_stage_fit` (the
linear solves) stays negligible (≤0.03 s) regardless of system size, so
`resp_charges_on = .false.` is now mostly useful for keeping the output
table to two columns, not for saving time.

(A separate idea - reusing the RI-J fitting coefficients the SCF already
solved for, instead of the exact density, to skip the integral pass
entirely - was tried and rejected: it needs a full RI-J re-solve anyway
whenever incremental Fock leaves the cache pointing at a delta density
(costing more than the exact path's entire ESP pass on a 68-atom test),
and the small RI-J fitting error on the raw electronic ESP gets amplified
~800x by the near-total cancellation between the nuclear and electronic
terms in the total ESP - final charges came out 18x further from Psi4.
See `esp_at_grid_batch_df` in `src/solvent/mod_integrals_cosmo.f90`,
unreachable from `resp_charges_on` but kept rather than deleted.)

## Output file (`<input>.out`)

```
estimate_only        = <T/F>
mem_grid_gb          = <DFT grid memory estimate, GB>
mem_2e_gb            = <2-electron integral memory estimate, GB>
converged            = <SCF convergence flag, 0/1>            (skipped if estimate_only)
total_energy_hartree = <total electronic+nuclear energy>      (skipped if estimate_only)
forces_hartree_per_bohr:
<Fx> <Fy> <Fz>   (one line per atom, same order as ncenters/atomchg -
                   this is dE/dR, NOT the physical force -dE/dR; see
                   the force_out sign-convention note above)
```

## Layout

```
data/               basis set library + physical-constant INCLUDE file
lib/                prebuilt third-party static libraries
build/              object/.mod output (gitignored)
src/main/           program entry point + EngineUp (library entry point)
src/harness/        run_engine.f90 (file-driven front end, see "Example"
                    below); src/harness/io/ parses &molecule/&atoms/&cosmo
                    (namelist/xyz/column &atoms formats, TOML export)
src/core/           MOL_info / GRID_info module state, basis-file lookup,
                    memory-budget prediction
src/integrals/      sole boundary to the electron-integral backend (libcint) -
                    mod_exchange.f90/mod_density_fitting.f90 (thin
                    public wrappers for exact-K/RI-J) sit here directly;
                    core/ shared setup, exact/ conventional J/K, df/
                    RI-J/RI-K, cosx/ chain-of-spheres exchange, force/
                    and hessian/ analytic derivatives
src/localization/   orbital localization (Pipek-Mezey)
src/dispersion/     empirical dispersion (-D2/-D3/-D3BJ/-D4, gCP/SRB)
src/xc/             sole boundary to the XC functional backend (libxc);
                    also DFT grid generation + GTO evaluation
                    (grid_gen.f90/gto_eval.f90/grid_dynamic.f90,
                    Lebedev.F) - shared with src/integrals/cosx/ and
                    src/force/, not XC-exclusive despite living here
src/scf/            SCF loop + Roothaan-equation solver, CPHF
src/force/          post-SCF force/gradient assembly
src/guess/          SCF initial guess, checkpoint/restart
src/util/           small matrix/IO debug helpers
src/solvent/        implicit continuum solvation (COSMO/CPCM + SMD) -
                    see "Solvation" above; cavity/ (GePol geometry),
                    cosmo_impl/ (SCF coupling, force), smd/ (CDS term)
```

EngineUp is safe to call more than once in the same process (resets all
module state via `reset_engine_state()` first).

Example: `examples/run_engine.f90` is a thin, file-driven front end
(namelist input, see `examples/friendly_format/water_gau.inp`) for
one-off/scripted runs and regression checks. It is NOT the intended
integration path for QM/MM or fragment-based-method use - link
`EngineUp` directly (e.g. via `iso_c_binding`) and call it in-process, in
a loop, instead.

**Note**: `EngineUp`'s `force_out` is the energy GRADIENT (dE/dR), not
the physical force (-dE/dR) - same sign convention as most QM codes'
"gradient" output, just named "force" here. Mind this if wiring it into
an optimizer or MD driver that expects -dE/dR.

**Known gap**: Mulliken charge analysis (`MLcharge_out`) is not
implemented; `EngineUp` reports it as zero rather than leaving it
undefined.

## Further reading

- [docs/DEVLOG.md](docs/DEVLOG.md) - project time-line and shipped-milestone
  history.
- [docs/ENVIRONMENT.md](docs/ENVIRONMENT.md) - every `ENGINE_*` tuning/
  debug environment variable, the `&professional` namelist that sets
  them from the input file, and wB97M-V/VV10 usage.
- [ARCHITECTURE.md](docs/ARCHITECTURE.md) - module boundaries, one
  `EngineUp` call's control flow, and the `src/solvent/` file layout.
- [VALIDATION.md](docs/VALIDATION.md) - cross-code verification numbers
  for basis import, ECP, solvation, dispersion.

## Third-party dependencies

| Dependency | License | Linking |
|---|---|---|
| `third_party/libcint` | Apache-2.0 | Static (built from source by `make`) |
| `third_party/libcint_ecp` | Apache-2.0 (vendored from PySCF) | Static (`lib/libcint_ecp.a`) |
| `third_party/simple-dftd3` | **LGPL-3.0-or-later** | **Dynamic** (`libs-dftd3.so`, `-Wl,-rpath` at link time) |
| `third_party/gcp` | **LGPL-3.0-or-later** | **Dynamic** (`libgcp.so`, `-Wl,-rpath` at link time) |
| `third_party/dftd4` | **LGPL-3.0-or-later** | **Dynamic** (`libdftd4.so`, `-Wl,-rpath` at link time) |
| `third_party/multicharge` | Apache-2.0 | Static (built from source by `make`) |
| `third_party/mctc-lib` | Apache-2.0 | Static (built from source by `make`) |
| `third_party/toml-f` | MIT / Apache-2.0 | Static (built from source by `make`) |
| `third_party/OpenBLAS` | BSD-3-Clause | Static (built from source by `make`, pthread-threaded build) |
| `third_party/libxc` | MPL-2.0 | Static (built from source by `make`) |

`simple-dftd3`/`gcp`/`dftd4` are the only copyleft ones, and each must
stay dynamically linked for that reason (LGPL permits linking from a
closed-source caller only when dynamic). `gcp` and `dftd4` back the
gCP/SRB basis-superposition-error correction and D4 dispersion used by
the "-3c" composite methods (B97-3c, r2SCAN-3c, wB97X-3c - see
`examples/benchmark_accuracy/README.md`).
