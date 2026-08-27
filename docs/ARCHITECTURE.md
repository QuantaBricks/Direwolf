# Engine — Program Architecture

This document describes how the codebase is organized: the module
boundaries, the call flow of one `EngineUp` invocation, and the design
principles that were established (and preserved through several refactors)
so future changes stay consistent with them.

## Top-level flow

`src/main/engine.f90` exposes the library's only two entry points,
`EngineUp` and `reset_engine_state`. One `EngineUp` call runs, in order:

```
reset_engine_state()          ! clears all module state (safe to call EngineUp repeatedly)
xc_select_functional(...)     ! picks the DFT functional, HF_exchange_frac ('HF' = plain HF, frac=1)
engine_use_df_j/engine_use_df_k/cosx_enabled/engine_df_aux_basis = j_mode/k_mode/ri_aux_basis
                               ! J: exact or RI. K: exact, RI, or cosx - independent of J
                               ! (real RIJCOSX is J='RI', K='cosx')
engine_mem_cap_bytes = mem_cap_gb    ! shared memory budget for every STORE-vs-DIRECT/cache decision below
integrals_init(info)          ! read basis, build S/Hcore, set up the 2e-integral path (STORE or DIRECT)
mem_2e_gb = nRec or nConts^2*nContsAux (whichever J/K pick) * 8 bytes   ! informational, always computed
gridgen(info)                 ! build the DFT quadrature grid (Becke partitioning + Lebedev)
mem_grid_gb = size(Grids)*nconts*4*8 bytes   ! informational, always computed when the grid runs at all
GTOeval(info)                 ! evaluate basis functions/gradients on every grid point
                               ! (gridgen/GTOeval both SKIPPED for plain HF - HF_exchange_frac>=1 -
                               ! since scf_build_fock never touches the grid in that case either;
                               ! GTOeval ALSO skipped, even when gridgen runs, in estimate_only mode)
[estimate_only mode returns HERE, before guess/SCF/force - see below]
guess(info, Gtype)            ! initial density guess
SCFcycle(info, ...)           ! the SCF loop itself (see below)
calc_force(info)              ! only if the driver's do_force argument is true
prof_report()                 ! print accumulated wall-clock profile
```

`EngineUp`'s calculation-mode parameters (J/K build mode + aux basis
name, whether to compute the gradient, a memory budget cap, and a
memory-estimation mode) are ordinary arguments, not env vars -
`j_mode`/`k_mode`/`ri_aux_basis`/`do_force`/`mem_cap_gb`/`estimate_only`,
set by the driver from the `&molecule` namelist's `J`/`K`/`ri_aux_basis`/
`calc_force`/`mem_cap_gb`/`estimate_only` fields (`examples/run_engine.f90`),
all optional with backward-compatible defaults (J/K exact, force on, no
memory cap, not estimate-only) so existing `.inp` files are unaffected.
Plain HF is requested via `functional = 'HF'` (mod_xc.f90's
`xc_select_functional`), same as any other functional string - no
separate flag needed.

`mem_cap_gb` (`mod_meminfo.f90`'s `engine_mem_cap_bytes`/
`get_memory_budget_bytes`) is the ONE budget every memory-aware
STORE-vs-DIRECT/cache decision in the codebase checks against - exact
TWOEI (`integrals_init`), DF's `dfB` (`build_df_integrals`), and the
grid `val0`/`val1` cache (`GTOeval`) - replacing three formerly
independent env-var overrides (`ENGINE_MAX_TWOEI_GB`,
`ENGINE_MAX_GRID_CACHE_GB`, and DF's lack of any cap at all). `<=0`
means "no cap, use whatever memory is actually available" (unchanged
default behavior).

`estimate_only=.true.` is a distinct calculation mode: `EngineUp`
computes `mem_grid_gb`/`mem_2e_gb` (how much RAM the grid cache and the
2e-integral store WOULD need for this molecule/basis/aux-basis, without
actually allocating them at full size - `integrals_init` already skips
materializing `TWOEI` via the same `skip_exact_2e` flag `engine_use_df`
uses, `integrals_estimate_aux_size` builds just enough of the aux basis
to know `nContsAux` without the expensive 3-center-integral/metric
build, and `GTOeval`'s per-point tabulation is skipped entirely) and
returns immediately - no `guess`/`SCFcycle`/`calc_force` ever runs, so
`energy_out`/`force_out`/`iconv`/`econv` stay at their zero defaults.
`mem_grid_gb`/`mem_2e_gb` are filled on ordinary (non-estimate) runs
too, since by the time they're computed the real numbers are already
known for free.

`examples/run_engine.f90` is the standalone driver (reads a namelist input
file, calls `EngineUp`, writes results) — the library itself has no
`program` statement, so it archives into `build/libengine.a` and any driver
can link against it.

## Directory / module map

```
src/core/       mod_data.f90        MOL_info: the shared module-state hub (geometry, basis,
                                     S/Hcore/Fa/Fb/Pa/Pb, nconts, ...). Nearly every other
                                     module does `use MOL_info` rather than passing this state
                                     as arguments.
                mod_meminfo.f90     get_available_memory_bytes() — reads /proc/meminfo's
                                     MemAvailable. Leaf module, no MOL_info dependency. Backs
                                     every STORE-vs-DIRECT memory decision in the codebase.

src/integrals/  mod_integrals.f90            Module state + public interface only (F2008
                                              "interface ... end interface" block, dummy-args
                                              only, no bodies). THE SOLE LIBCINT BOUNDARY —
                                              no file outside this module family calls libcint.
                mod_integrals_core.f90       submodule: integrals_init, 1e integrals (S/Hcore),
                                              basis/aux-basis file parsing, eri_get/INDEX_2E.
                mod_integrals_coulomb.f90    submodule: exact (non-DF) Coulomb build, incl.
                                              incremental-Fock bookkeeping.
                mod_integrals_exchange.f90   submodule: exact HF exchange build.
                mod_integrals_force.f90      submodule: force/gradient two-electron derivatives.
                mod_integrals_df.f90         submodule: density-fitting (RI-J) Coulomb build,
                                              both STORE and DIRECT modes; screening lists.
                                              (mod_density_fitting.f90 below is DF's public
                                              wrapper; the DF math itself lives here so it can
                                              share libcint access and module state with the
                                              rest of the integrals family.)

src/exchange/   mod_exchange.f90    Public wrapper around mod_integrals_exchange (HF_exchange_frac,
                                     exchange_build) — SCF/force code goes through this, not
                                     mod_integrals directly.

src/density_fitting/
                mod_density_fitting.f90  Public wrapper: df_build_coulomb(nconts, Ptot, J),
                                     dispatches into mod_integrals_df's STORE/DIRECT paths.

src/xc/         mod_xc.f90          THE SOLE LIBXC BOUNDARY — functional selection
                                     (xc_select_functional), libxc calls.
                DFT.f90             DFT_calc: XC quadrature over the grid (uses mod_xc + the
                                     grid/GTO data from src/grid).

src/grid/       grid_gen.f90        Becke-weight grid construction (gridgen), with
                                     distance-pruned (PAIR_CUTOFF) partitioning.
                gto_eval.f90        Basis-function values/gradients on every grid point (GTOeval).
                Lebedev.f            Lebedev angular quadrature tables (F77 source).

src/guess/      guess.f90           Initial density guess (core Hamiltonian / SAD-style).

src/scf/        SCF.f90             SCFcycle: the iteration loop skeleton — density build,
                                     damping, convergence check, timing/print. Delegates the
                                     three heaviest/most distinct pieces below rather than
                                     inlining them.
                scf_diis.f90        diis_extrapolate: DIIS/ADIIS Fock extrapolation (history
                                     buffers owned by SCFcycle, passed in/out).
                scf_adiis.f90       ADIIS building blocks (build_adiis_terms, adiis_powell,
                                     adiis_fx) used by scf_diis.f90.
                scf_fock.f90        scf_build_fock: builds Fa/Fb and the energy contribution
                                     E_n for the current density (DFT_calc + Coulomb build,
                                     dispatching STORE/DIRECT/DF via use_df + exact exchange).
                solveKS.f90         solHFR_KS: Roothaan-Hall generalized eigenproblem solve
                                     (F C = S C E).

src/force/      force.f90           Nuclear-gradient assembly (1e + 2e derivative contributions,
                                     calls into mod_integrals_force via mod_integrals).

src/util/       mod_profile.f90     Wall-clock (system_clock) named-section timers:
                                     prof_start/prof_stop/prof_report/prof_reset.
                tool_io.f90         Small I/O helpers.
                tool_mat.f90        Small matrix helpers.

src/solvent/    mod_integrals_cosmo.f90  The one integral solvation needs: a tessera's
                                     electrostatic-potential matrix over the AO basis
                                     (ghost-point-charge trick reusing the existing
                                     nuclear-attraction kernel), energy and force variants.
                cavity/             Pure geometry - atom-centered-sphere cavity construction
                                     (4 tessellation schemes), Gauss-Bonnet tessera area/
                                     centroid + analytic Jacobians for gepol's force path.
                                     No libcint dependency. Dependency order top to bottom:
                                     mod_cosmo_constants -> mod_cosmo_dual (dual-number AD
                                     scaffolding) -> mod_cosmo_radii (per-element radius
                                     tables) -> mod_cosmo_geodesic (sphere discretization)
                                     -> mod_cosmo_clip (GePol multi-sphere clipping) ->
                                     mod_cosmo_gaubon (exact tessera area/centroid) ->
                                     mod_cosmo_cavity (orchestrator; also holds the 3
                                     alternative energy-only cavity types and the
                                     sigma-profile histogram).
                smd/                SMD's CDS (Cavitation/Dispersion/Solvent-structure)
                                     term - direct Fortran port of NWChem's MNSOL/MN-GSM
                                     reference (mnsol.F, Truhlar group; same code pyscf's
                                     own SMD module calls via its compiled libsolvent.so).
                                     Pure function of geometry, added directly to E/atmForce
                                     (engine.f90) alongside dispersion, not in the SCF Fock
                                     build. mod_smd_tables (constants/tables) ->
                                     mod_smd_sts (per-atom surface tension + gradient) ->
                                     mod_smd_daareal (analytic solid-angle SASA + gradient,
                                     the largest piece) -> mod_smd_cds (orchestrator).
                cosmo_impl/         Top-level driver, split by lifecycle stage, sharing
                                     state via mod_cosmo_state.f90's module variables:
                                     mod_cosmo_solvents (named-solvent lookup, stateless) ->
                                     mod_cosmo_init (cosmo_init, once before SCF) ->
                                     mod_cosmo_scf (cosmo_scf_step per iteration,
                                     cosmo_report_sigma_profile - see docs/
                                     ADVANCED_FEATURES.md) -> mod_cosmo_force
                                     (post-SCF analytic gradient) -> cosmo.f90
                                     (mod_cosmo orchestrator, re-exports everything).
```

## Design principles established this session (preserve these)

- **`mod_integrals.f90` is the sole libcint boundary.** When it grew past
  2800 lines, it was split by *file*, not by *responsibility*: the parent
  module still owns all state and the public interface; five Fortran
  submodules (`mod_integrals_core/coulomb/exchange/force/df.f90`) hold the
  actual bodies, split by feature (1e-setup, exact Coulomb, exact exchange,
  force, density fitting). Submodules see the parent's full state
  (public *and* private) automatically — no `use` statement needed — and
  compile against `mod_integrals.o`'s `.mod` file (Makefile enforces this
  ordering). This was a deliberate choice over giving each feature its own
  independent module with direct libcint access, specifically to keep
  "only one file ever calls libcint" true. Any new integral-level feature
  belongs in this family, as a new submodule if it's substantial, not as
  a new libcint call site elsewhere.
- **`mod_xc.f90` is the sole libxc boundary**, by the same logic — DFT.f90
  goes through it rather than calling libxc directly.
- **`MOL_info` (`mod_data.f90`) is the shared state hub.** Convention
  (matching `DFT_calc`/`solHFR_KS`/`scf_build_fock`): subroutines that need
  `nconts`/`Fa`/`Fb`/`Hcore`/`Pa`/`Pb` do `use MOL_info` and reference the
  module-level variables directly, rather than threading them through the
  argument list. (Passing one of these as an explicit dummy argument *and*
  also doing `use MOL_info` in the same scope is a compile error —
  "ambiguous reference" — since Fortran can't disambiguate two entities
  with the same name in one scope.)
- **STORE-vs-DIRECT memory-aware dispatch**, used identically in two
  independent places: the exact 4-center integral path (`integrals_init`,
  governed by `ENGINE_MAX_TWOEI_GB`) and the DF 3-center path
  (`build_df_integrals`'s `df_direct_mode`, governed by
  `ENGINE_FORCE_DF_DIRECT`). Both estimate the tensor's memory footprint,
  compare against `0.5 * available - 4GB reserve` (from
  `mod_meminfo.f90`), and pick between storing the full tensor once vs.
  recomputing integrals on the fly every SCF iteration. New large-tensor
  features should follow this same pattern rather than assuming either
  mode unconditionally.
- **Incremental Fock** (`J(Ptot) = J(Ptot_prev) + J(ΔP)`, exact by
  linearity) is implemented once per path (exact: `mod_integrals_coulomb`;
  DF: `mod_integrals_df`) via shared module-level `incr_*` state
  (`incr_Ptot_prev`, `incr_J_prev`, `incr_has_prev`, `incr_since_full`,
  `incr_full_period=20`). The two paths share the same state variables
  because they're mutually exclusive per run (a run is either exact or DF,
  never both) — do not assume both could be "active" simultaneously if
  extending this.
- **`SCF.f90` is a thin loop skeleton, not a monolith.** It was split
  incrementally (three separate, individually-verified extractions) into:
  `scf_adiis.f90` (self-contained ADIIS math), `scf_diis.f90`
  (DIIS/ADIIS extrapolation, history buffers still owned/allocated by
  `SCFcycle` and passed in/out), `scf_fock.f90` (Fock-matrix build +
  energy for the current density). `SCFcycle` itself keeps only: DIIS
  history allocation, the iteration loop calling into the three pieces
  above in order, density-matrix construction + damping, and the
  convergence check. This was an explicit, incremental request ("一点点拆开")
  — further splits should follow the same pattern (extract one
  self-contained, verifiable piece at a time, confirm bit-exact energies
  against an established baseline, commit, then continue) rather than a
  single large rewrite.
- **Wall-clock profiling via `mod_profile.f90`** (`system_clock`, not
  `cpu_time` — `cpu_time` sums across threads, which would misreport if
  OpenMP is ever added) wraps every major phase
  (`basis_and_integrals`, `grid_gen`, `gto_eval`, `scf_total`,
  `coulomb_build`, `exchange_build`, `dft_xc`, `df_direct_bpass`,
  `df_direct_solve`, `df_direct_jpass`, `force`, ...), printed once via
  `prof_report()` at the end of `EngineUp`. Any new hot path should be
  wrapped in a named `prof_start`/`prof_stop` pair rather than left
  unprofiled.

## Build

`Makefile`'s `VPATH` lists all `src/*` subdirectories plus `examples`, so
source files are referenced by filename only in `objects`/dependency rules
regardless of which subdirectory they live in. Module (`.mod`) dependency
ordering is expressed via explicit prerequisite rules (e.g.
`mod_integrals_core.o : mod_integrals.o`) since `make` can't infer Fortran
module dependencies on its own. The five `mod_integrals_*.f90` submodules
have no ordering constraint *among themselves* — only against their parent
`mod_integrals.o`.
