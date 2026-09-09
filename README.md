# Direwolf

<img src="imgs/quantabricks_logo.png" alt="QuantaBricks" width="220">

*A [QuantaBricks](https://www.quantabricks.xyz) product.*

An independent quantum-chemistry engine - Hartree-Fock and DFT energies,
analytic gradients, and properties - built for speed and for feeding
results to other software.

**Status**: in active development (AI-assisted); interfaces may still
change between commits. See [docs/DEVLOG.md](docs/DEVLOG.md) for history
and what is shipped.

## Install

```
make
```

Needs `gfortran`, `gcc`, `cmake`. Every dependency (libcint, libxc,
OpenBLAS, the D3/D4 dispersion chain) is vendored as source under
`third_party/` and built automatically on the first `make` - expect a
few minutes the first time, fast afterwards. Produces `./Direwolf`.

Or with Docker (`Dockerfile` in the repo root):

```
docker build -t direwolf .
docker run --rm -v "$PWD":/data direwolf /data/input.inp /data/output.out
```

(mount input and output under the **same** directory - Direwolf writes
the output next to its own working directory, not next to the input.)

## Run

```
./Direwolf <input.inp> [output.out]
```

Output defaults to `<input.inp>.out`. Single-threaded unless
`OMP_NUM_THREADS` is set or `n_threads` is given in the input file.
A minimal input:

```
&molecule
 ncenters   = 3
 imult      = 1
 icharge    = 0
 functional = 'B3LYP'
 baselabel  = '6-31g'
&end
&atoms
 O   0.0000   0.0000   0.0000
 H   0.7584   0.0000   0.5861
 H  -0.7584   0.0000   0.5861
&end
```

More worked examples in [`examples/`](examples/) (water up to ~100-basis
-function drug molecules); the full input grammar is under "Input file
format" below.

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


Fortran namelist. Working examples: `examples/friendly_format/`,
`examples/drug_molecules/`, `examples/benchmark_accuracy/`,
`examples/solvation/`, `examples/misc/` (`examples/README.md` groups them).
Namelist terminator is `/` or `&end`.

```
&molecule
 functional = 'B3LYP'
 baselabel  = 'def2svp'
&end
&atoms
 O    0.000000   0.000000   0.000000
 H    0.758400   0.000000   0.586100
 H   -0.758400   0.000000   0.586100
&end
```

Everything except a geometry has a default (shown below), so a minimal
input is just `functional` + `baselabel` + `&atoms`.

**`&molecule`**

| tag | default | meaning |
|---|---|---|
| `functional` | `'PBE_PBE'` | `'HF'`, `'B3LYP'`, `'PBE0'`, `'WB97M-V'`, … (full list under "Supported functionals and methods" above); unknown name falls back to plain PBE |
| `baselabel` | `'def2svp'` | global basis; must match a `<Z>-<label>` block in `data/bases` |
| `ecplabel` | `''` | global ECP; `''` = whatever the basis file carries |
| `imult` | `1` | spin multiplicity 2S+1 (1 = closed-shell singlet) |
| `icharge` | `0` | total charge |
| `J` | `'RI'` | Coulomb build: `'exact'` or `'RI'` |
| `K` | `'cosx'` | exchange build: `'exact'`, `'RI'`, or `'cosx'`; independent of `J` (default `RI`+`cosx` = RIJCOSX) |
| `ri_aux_basis` | `''` | aux basis for the RI path; `''` = auto even-tempered. Bundled: cc-pVDZ/TZ/QZ-JKFIT, def2-universal-JFIT/JKFIT |
| `calc_force` | `.true.` | also compute the analytic gradient |
| `spherical` | `.true.` | spherical (puream) vs Cartesian Gaussians |
| `scf_conv` | `'regular'` | `'regular'` (dE 1e-6 / dRMS 1e-5), `'fine'` (1e-8 / 1e-7), `'tight'` (1e-9 / 1e-8) |
| `vv10_nonself` | `.false.` | evaluate VV10 non-local correlation post-SCF instead of every cycle |
| `mem_cap_gb` | `0.0` | soft memory budget, GB; `0.0` = auto-detect |
| `estimate_only` | `.false.` | skip SCF/force, just print the memory estimate |
| `n_threads` | `0` | OpenMP threads; `0` = leave `OMP_NUM_THREADS` alone |
| `verbose` | `1` | `0` quiet, `1` normal, `2` debug |
| `resp_charges_on` | `.true.` | two-stage RESP atomic charges after SCF (see "RESP charges") |
| `ncenters` | `0` | atom count; `0` = infer from `&atoms` line count |
| `unit` | `'angstrom'` | `'angstrom'` or `'bohr'` (applies to `&atoms`) |
| `basedir` | `''` | directory to resolve `data/` against; `''` = built-in |
| `xyzfile` | `''` | read geometry from this `.xyz` instead of `&atoms` |

**`&atoms`** — one atom per line `element x y z` (no commas/quotes). Two
other accepted forms: namelist arrays `atomchg` / `x` / `y` / `z`, or an
external `xyzfile`. Per-atom overrides: `atom_basis(i)` / `atom_ecp(i)`
(`''` in a slot = use the global label).

**`&pointcharges`** (QM/MM) — `npc` (count), then arrays `pc_q` / `pc_x` /
`pc_y` / `pc_z`.

**`&checkpoint`** — `chk_read` / `chk_write` (`.false.`), `chk_file` (`''`):
density-matrix warm start to/from a file.

**`&molden`** — `molden_write` / `molden_read` (`.false.`), `molden_file` /
`molden_read_file` (`''`).

**`&cosmo`** (implicit solvation; all inert unless `cosmo_on`)

| tag | default | meaning |
|---|---|---|
| `cosmo_on` | `.false.` | enable COSMO |
| `cosmo_solvent` | `''` | named solvent (`'water'`, …); sets the dielectric |
| `cosmo_epsilon` | `78.4` | dielectric (ignored if `cosmo_solvent` set) |
| `cosmo_radii_scale` | `1.2` | Bondi-radius multiplier for the cavity |
| `cosmo_avg_area` | `0.3` | target tessera area, Å² |
| `cosmo_cavity_type` | `'gepol'` | only `'gepol'` (has an analytic force) |
| `cosmo_smd` | `.false.` | add the SMD non-electrostatic term (needs `gepol`) |
| `cosmo_rsolv` | `0.0` | solvent probe radius |
| `cosmo_sigma_rav` | `0.5` | σ-profile averaging radius, Å |
| `cosmo_sigma_profile_file` | `''` | write a σ-profile here |
| `cosmo_ks_nseg` / `cosmo_ks_nface` | `92` / `1082` | Klamt-surface segmentation |


**`&professional`** — advanced performance tuning (XC grid cache, COSX grid
density, …); see `docs/ENVIRONMENT.md`.


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
Mulliken/Lowdin. Set `resp_charges_on = .false.` in `&molecule` to skip
the fit. Full details, defaults, and cross-code verification are in
[docs/RESP_CHARGES.md](docs/RESP_CHARGES.md).

## Geometry optimization

Add an `&opt` namelist (its mere presence is enough) to minimize the
geometry instead of doing a single point. An ORCA/Gaussian input with an
`Opt` route token maps to the same path.

```
&molecule
 functional = 'PBE'
 baselabel  = 'def2svp'
 xyzfile    = 'start.xyz'
&end
&opt
 opt_coord       = 'ric'
 opt_hessian_file = 'hessian'   ! optional: seed from e.g. `xtb --hess`
&end
```

Quasi-Newton: RFO step with a trust radius and BFGS Hessian updates, run
in **redundant internal coordinates** by default (`opt_coord = 'cart'`
for Cartesians). The Hessian is seeded from the Lindh model; supply
`opt_hessian_file` to start from an external Cartesian Hessian
(Hartree/bohr², free-form - `xtb --hess` output is read as-is) instead.
Every cycle is appended to `<output>.opt.xyz` (multi-frame XYZ); the
converged geometry is also printed in the output as `[FINAL GEOMETRY]`.
`opt_restart` resumes from the last frame of `<output>.opt.xyz`.

**`&opt`** (its presence alone enables optimization; every field optional)

| tag | default | meaning |
|---|---|---|
| `opt_maxcyc` | `100` | max optimization cycles |
| `opt_conv` | `'normal'` | `'normal'` (gmax 4.5e-4 / grms 3.0e-4 / dmax 1.8e-3 / drms 1.2e-3 / dE 1e-6, the Gaussian set) or `'tight'` (gmax 1.5e-5 / grms 1.0e-5 / dE 1e-8) |
| `opt_trust` | `0.3` | initial trust radius, bohr |
| `opt_coord` | `'ric'` | `'ric'` (redundant internals) or `'cart'` |
| `opt_restart` | `.true.` | resume from `<output>.opt.xyz` if present |
| `opt_hessian_file` | `''` | initial Cartesian Hessian to seed from instead of the Lindh model |
| `opt_write_ric` | `.false.` | dump the generated internal-coordinate set to `<output>.ric` |

`opt_conv = 'tight'` also raises the SCF to `'tight'` automatically (an
optimizer needs gradients cleaner than its own convergence threshold),
and any `&opt` run raises the SCF to at least `'fine'`.

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

## How to cite

If you use Direwolf in published work, please cite it as:

> Y. Ma and X. Chen, Direwolf, Revision V0.99, QuantaBricks, Newark NJ, 2026.

```bibtex
@misc{direwolf,
  author       = {Ma, Yizhou and Chen, Xin},
  title        = {Direwolf},
  howpublished = {Revision V0.99, QuantaBricks, Newark NJ},
  year         = {2026}
}
```

### Please also cite the methods you use

Direwolf implements published methods. Citing Direwolf does not replace
citing the work that a given calculation actually relies on. The list
below covers the methods most commonly used through Direwolf; it is
**not exhaustive** - functionals, basis sets, dispersion models and
algorithms not named here still carry their own original references, and
those should be cited too.

**Semi-numerical exchange (COSX / chain-of-spheres), used by `K = 'cosx'`
and therefore by the default RIJCOSX path**

> F. Neese, F. Wennmohs, A. Hansen, U. Becke, *Chem. Phys.* **356**, 98 (2009).
> R. Izsák, F. Neese, *J. Chem. Phys.* **135**, 144105 (2011).
> B. Helmich-Paris, B. de Souza, F. Neese, R. Izsák, *J. Chem. Phys.* **155**, 104109 (2021).

**Numerical integration grid** - every DFT calculation in Direwolf (the
XC quadrature, the VV10 non-local grid and the COSX grid alike)
partitions space with Becke's multicenter scheme, combined with Lebedev
angular quadrature and Treutler-Ahlrichs radial mapping

> A. D. Becke, *J. Chem. Phys.* **88**, 2547 (1988).  (multicenter numerical integration / Becke partitioning)
> V. I. Lebedev, D. N. Laikov, *Dokl. Math.* **59**, 477 (1999).  (Lebedev angular grids)
> O. Treutler, R. Ahlrichs, *J. Chem. Phys.* **102**, 346 (1995).  (radial mapping / M3-M4 grids)

**Becke exchange and the B3LYP hybrid**

> A. D. Becke, *Phys. Rev. A* **38**, 3098 (1988).  (B88 exchange, used by BLYP/B3LYP)
> A. D. Becke, *J. Chem. Phys.* **98**, 5648 (1993).  (three-parameter hybrid, B3LYP)
> C. Lee, W. Yang, R. G. Parr, *Phys. Rev. B* **37**, 785 (1988).  (LYP correlation)

**Composite "-3c" methods (`B97-3c`, `r2SCAN-3c`, `wB97X-3c`)**

> S. Grimme, J. G. Brandenburg, C. Bannwarth, A. Hansen, *J. Chem. Phys.* **143**, 054107 (2015).  (PBEh-3c, the family's original)
> J. G. Brandenburg, C. Bannwarth, A. Hansen, S. Grimme, *J. Chem. Phys.* **148**, 064104 (2018).  (B97-3c)
> S. Grimme, A. Hansen, S. Ehlert, J.-M. Mewes, *J. Chem. Phys.* **154**, 064103 (2021).  (r2SCAN-3c)
> M. Müller, A. Hansen, S. Grimme, *J. Chem. Phys.* **158**, 014103 (2023).  (wB97X-3c)

**Dispersion and basis-set-superposition corrections**

> S. Grimme, J. Antony, S. Ehrlich, H. Krieg, *J. Chem. Phys.* **132**, 154104 (2010).  (D3)
> S. Grimme, S. Ehrlich, L. Goerigk, *J. Comput. Chem.* **32**, 1456 (2011).  (D3 Becke-Johnson damping)
> E. Caldeweyher, C. Bannwarth, S. Grimme, *J. Chem. Phys.* **147**, 034112 (2017).  (D4)
> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth, S. Grimme, *J. Chem. Phys.* **150**, 154122 (2019).  (D4)
> H. Kruse, S. Grimme, *J. Chem. Phys.* **136**, 154101 (2012).  (gCP geometrical counterpoise)

**VV10 non-local correlation and the VV10-containing functionals
(`wB97M-V`, `wB97X-V`, `B97M-V`)**

> O. A. Vydrov, T. Van Voorhis, *J. Chem. Phys.* **133**, 244103 (2010).  (VV10)
> N. Mardirossian, M. Head-Gordon, *J. Chem. Phys.* **144**, 214110 (2016).  (wB97M-V)
> N. Mardirossian, M. Head-Gordon, *Phys. Chem. Chem. Phys.* **16**, 9904 (2014).  (wB97X-V)
> N. Mardirossian, M. Head-Gordon, *J. Chem. Phys.* **142**, 074111 (2015).  (B97M-V)

**Other functionals** (r2SCAN, TPSS, M05/M06/MN15 family, B3LYP, PBE0,
CAM-B3LYP, wB97X-D, ...) are evaluated through **libxc**; cite the
functional's own original paper and

> S. Lehtola, C. Steigemann, M. J. T. Oliveira, M. A. L. Marques, *SoftwareX* **7**, 1 (2018).  (libxc)

**Underlying libraries** - integrals via **libcint** (Q. Sun,
*J. Comput. Chem.* **36**, 1664 (2015)); D3/D4/gCP through the
Grimme-group `simple-dftd3`, `dftd4` and `gcp` libraries; linear algebra
via OpenBLAS. Basis sets carry their own references (see the Basis Set
Exchange entry for the set you use).

## License

Copyright © 2025-2026 QuantaBricks.

This open version of Direwolf is licensed under the **GNU Affero General
Public License v3.0 or later** (`AGPL-3.0-or-later`) - see
[LICENSE](LICENSE). Under AGPL §13, offering a modified Direwolf as a
network service obliges you to publish that version's complete source.
Bundled third-party components keep their own licenses (table above).

**Commercial and production use → Direwolf Advanced.** The Advanced
version is proprietary-licensed, carries the extra input/interoperability
and performance modules, and comes without the AGPL copyleft obligation.
Contact [contact@quanta-bricks.com](mailto:contact@quanta-bricks.com) /
[www.quantabricks.xyz](https://www.quantabricks.xyz).

**Data and model training.** Numerical output produced by Direwolf
(energies, forces, properties, and datasets built from them) is yours.
Using Direwolf to generate data and training machine-learning models on
that data is expressly permitted; the resulting datasets and models are
not covered by the AGPL and are not derivative works of Direwolf.
