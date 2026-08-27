All examples use the friendly `&atoms` form - a plain "element x y z"
block, one atom per line, no commas or quotes. `ncenters` is inferred
from the atom count; `imult`/`icharge` default to 1/0 and are set
explicitly in `&molecule` only when non-default (see any lanl2dz
example below). The old comma-separated `atomchg`/`x`/`y`/`z` array
form is no longer used anywhere in this tree, though the parser
(`examples/run_engine.f90`) still accepts it for backward
compatibility.

Two atom-block variants exist:
- element/x/y/z lines directly in `&atoms` (most examples).
- an external `xyzfile` pointing at a standard 2-line-header .xyz file
  (see `friendly_format/water_xyz.inp` / `water.xyz`).

Both are parsed by the same `&molecule`/`&atoms` namelist front end
and accept the same `&molecule` keys (functional, baselabel, J, K,
calc_force, ...). J picks the Coulomb build mode ('exact', the
default, or 'RI'); K picks exchange ('exact' default, 'RI', or 'cosx')
independently of J - the real RIJCOSX recipe is J='RI', K='cosx'.

Directories are grouped by purpose:

- `friendly_format/` - plain-basis water inputs in every atom-block
  variant, plus every ECP (LANL2DZ) example: bare metal atoms
  (`fe_bare_lanl2dz.inp`, `fe3_bare_lanl2dz.inp`,
  `pd_bare_lanl2dz.inp`, `pt_bare_lanl2dz.inp` - single heavy atom, no
  light atoms, so a global `baselabel='lanl2dz'` is correct as-is),
  and mixed-element molecules (`fe_hydrate_lanl2dz.inp`,
  `pd_hydrate_lanl2dz.inp`, `pt_hydrate_lanl2dz.inp`,
  `fe_plus_farH.inp`, `pd_plus_farH.inp`). The mixed ones use the
  documented per-atom override pattern - global `baselabel='6-31g'`
  for the light atoms (O/H), with `lanl2dz lanl2dz` override columns
  on only the metal atom's line - not a global `baselabel='lanl2dz'`,
  which would incorrectly push the light atoms onto LANL2DZ's own
  (small/dated) main-group basis too. See the top-level README's "ECP"
  section for the override-column syntax.
- `drug_molecules/` - drug-like test molecules (aspirin, caffeine,
  imatinib, ...), mostly for performance/scale testing.
- `benchmark_accuracy/` - one functional per file, all on the same
  water geometry, cross-checked against real Psi4 output
  (`references.csv`/`run_benchmark.py`) - see that directory's own
  README. Includes `water_wb97mv.inp` (range-separated hybrid meta-GGA
  + VV10 nonlocal correlation - see docs/ENVIRONMENT.md for the VV10
  grid-density knobs) and `he2_vv10.inp` (plain VV10).
- `solvation/` - implicit continuum solvation (`&cosmo` namelist):
  - `water_cosmo.inp` - plain CPCM/GePol electrostatics in water.
  - `water_smd_toluene.inp` - full SMD (electrostatics + CDS
    non-electrostatic term) in toluene.
  - `water_sigma_profile.inp` - COSMO-RS-style sigma-profile export
    (`cosmo_sigma_profile_file`, `cosmo_cavity_type='ks1993'`) - an
    advanced/non-OSS feature, see docs/ADVANCED_FEATURES.md.
    `ks1993`/`iswig` cavities are energy-only (no analytic gradient),
    so this example also sets `calc_force = .false.`.
- `misc/` - feature demos that aren't drug molecules or ECP:
  `c8h18.inp`/`c8h18_rijk.inp` and `sh2.inp`/`sh2_rijk.inp` (plain vs
  DF/RI-JK), `water_bohr.inp` (`unit='bohr'`), `water_element.inp`
  (explicit `element` array), `water_chkread.inp`/`water_chkwrite.inp`
  (`&checkpoint`), `water_pc.inp` (`&pointcharges`, QM/MM), and
  `water_rijk.inp` (plain water with DF/RI-JK).
