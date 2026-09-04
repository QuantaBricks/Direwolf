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
Mulliken/Lowdin. Set `resp_charges_on = .false.` in `&molecule` to skip
the fit. Full details, defaults, and cross-code verification are in
[docs/RESP_CHARGES.md](docs/RESP_CHARGES.md).

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
