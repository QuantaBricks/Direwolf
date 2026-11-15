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

These are `EngineUp` arguments; all are ignored unless `cosmo_on = .true.`. In an input file the same settings are the keys of the
`&cosmo` block (`[cosmo]` in TOML) without the `cosmo_` prefix (`solvent`, `smd`, `epsilon`, `radii_scale`, ...); the block alone
switches the solvent on and `on = .false.` inside it switches it off. The prefixed spelling is still read. `convention = 'cpcm'` (`opts%cosmo_cpcm`) selects the C-PCM convention: f(eps) = (eps-1)/eps, hydrogen
cavity radius 1.10 A instead of Bondi's 1.20, and the same f for the optical response of a TDDFT run (`tddft_eps_opt`); the
default `'cosmo'` keeps Klamt's f = (eps-1)/(eps+0.5) and Bondi radii.

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

### Post-SCF tasks (`opts`, last argument)

`EngineUp`'s basic task is: initial guess, SCF, energy (and the force, when `do_force`). Everything else it can do
after a converged SCF is a **post-SCF task**, and which of them to run is asked through one argument,
`type(engine_options_t), intent(inout) :: opts` (`use mod_engine_options`). It is the only channel for these
requests and for what comes back from them (apart from `ENGINE_*` environment variables): `EngineUp` reads no state
that the caller set in other modules, and a caller (`run_engine`, an optimizer, an MD or QM/MM loop) only fills `opts`
and calls. Every field defaults to "off", so `type(engine_options_t) :: opts` left alone gives the basic task.

What belongs in `opts`: which post-SCF tasks to run, their options, how their report is written, and what they hand
back. What does not: the molecule, the method, the basis and the numerical thresholds (the arguments above), and
`do_force`, which stays a separate argument. A new post-SCF task adds fields to the type and its code to
`engine.f90`; the signature and the other callers stay as they are.

| task | fields (in unless noted) | meaning |
|---|---|---|
| Excited states | `tddft_nroots` | number of states; 0 = none |
| | (forces) | a single-point run with `tddft_nroots` and no `state` does not compute the ground-state forces (`calc_force` is turned off); `state`, `opt`, `irc` and `freq` runs keep them |
| | `tddft_robust`, `tddft_triplets`, `tddft_tda` | spare roots; triplets of a closed-shell reference; Tamm-Dancoff instead of the full problem |
| | `tddft_eps_opt` | optical dielectric constant (n^2, default 1.78) of the solvent response when COSMO/SMD is on |
| Polarizability | `polar`, `polar_omega` | static (and dynamic) polarizability; frequencies of the dynamic one, as in the input (`'0.05, 589nm'`) |
| Hessian, frequencies | `freq` | analytic Hessian and harmonic frequencies |
| | `hessian_return` | keep the Hessian (in `<freq_base>.hess`) for the caller and continue (a transition-state optimizer asks for it mid-run) |
| | `freq_thermo`, `freq_unit`, `freq_mult`, `freq_base` | thermochemistry block on/off; unit of the report; multiplicity for the thermochemistry; stem of `<stem>.hess`, `.freq.molden`, `.modes.csv` |
| Bond orders | `mayer` | Mayer bond orders of the converged density |
| | `mayer_bo` (out, `real(8)(:,:)` allocatable) | (natoms,natoms), filled when `mayer` is set and the SCF converged |
| | `progress` (in, `procedure(report_hook_i), pointer`) | called with the record so far (`opts%rep` content) whenever a part of it is complete: the result energy of a polarizability/TDDFT/frequency run (`energy_only`), the static tensor, each TDDFT state, the frequencies, the thermochemistry. `run_engine` sets it to `report_progress`, which rewrites `<base>.json`, so an interrupted long call leaves what was done |
| | `done` (out) | a Hessian/frequency run (`freq` without `hessian_return`): its output is complete, the caller ends the run (no results block); `run_engine` and the optimizer do, after writing the report |
| | `rep` (out, `engine_report_t`) | what the call computed besides energy and forces; each part has a `has_*`/`used` flag: `has_props` (dipole `dip_au`, Mulliken/Lowdin/RESP charges, orbital energies), `has_polar` (`stat`, `dyn(1:n_dyn)`), `tddft(1:3)` (singlets, triplets, open-shell; per state `w`, `f`, `mu`, main transitions), `has_freq` (`freq`, `inten`, `redmass`, `fc`, `disp`), `has_thermo` (`th_*`), `has_dipder` (`dmu`), `energy_only` (the run left EngineUp before the results block). Cleared at the start of every call |
| Run policy | `verbose` | 1 standard report, 2 debug extras |
| | `quiet` | a multi-step driver (optimizer, IRC) silences the per-call report |
| | `force_dense_allow` | `.false.`: the force stays on the SCF grid (an optimizer must differentiate the function it minimizes); default `.true.` |
| | `cosx_force_loose` | COSX gradient at the energy K screen (early optimizer steps) |
| Charges, COSMO-RS | `resp_charges` | RESP charges after the SCF |
| | `cosmors_on`, `cosmors_file`, `cosmors_compound` | write the COSMO-RS `.cosmo` file (advanced tier) |
| SCF report (out) | `scf_niter`, `scf_dP`, `scf_dE` | iterations, last density rms change and energy change |
| | `e_hcore`, `e_coul`, `e_exx`, `e_xc`, `e_nuc`, `vv10_used`, `e_vv10` | the energy parts of the converged SCF (Ha) |
| Reference | `uks` | unrestricted orbitals for a closed-shell molecule |
| | `uks_restart` (out) | `stability = 2` found the restricted solution unstable toward an unrestricted one: call again with `opts%uks = .true.`, passing `dens_out_*` as `dens_in_*` (`run_engine` does) |

EngineUp sets its internal state (verbosity, quiet, grid policy, the COSMO-RS request) from these fields at every call and
reads nothing a caller set in another module. It writes no result file of its own: everything it produces besides the energy
and force arguments comes back in `opts%rep` (`engine_report_t`, `mod_engine_report`), and the driver turns that into
`<base>.json` (`json_write_report`, called through `finish_engine_call` after every call). The one process-wide sink left is
the timing report (`mod_profile`).

`engine_output_stem(outfile)` (same module) turns an output file name into the `freq_base` stem. The densities a
caller wants back are the `dens_out_*` arguments, not part of `opts`.

---

## NMR shieldings (`runtype = 'nmr'`)

Magnetic shielding tensors sigma_ab = d2E/dB_a dm_b (ppm) of every nucleus with London orbitals, Hartree-Fock and LDA/GGA density functionals, pure or global hybrid, restricted closed-shell or unrestricted (`nmr_xc_giao` adds the London-phase derivative of the xc matrix; the xc kernel does not respond to the imaginary first-order density). Range-separated hybrids add the erf-attenuated exchange (`integrals_giao_jk` / `integrals_k_general` take the attenuation `omega`; K = cx K + beta K_LR). Meta-GGAs use the gauge-invariant kinetic-energy density tau~ = tau_p - |j_p|^2/(2 rho) (tau_p = 1/2 sum |grad psi|^2, j_p the paramagnetic current), which equals the covariant tau - |j|^2/(2 rho), is exact for one-electron systems in a field and does not contain the vector potential explicitly: the field enters only through the orbitals. With London orbitals the first-order Fock matrix gets the phase and a_nu terms of tau_p plus the explicit-current source of the j^2 term (`nmr_tau_block`), and the coupled-perturbed operator gets the current kernel -(v_tau/rho) J[d] . dj/dD (`nmr_xc_resp`), applied to the orbital response and to the orthonormality density; there is no explicit xc term for the nuclear moment. The result is independent of the gauge origin (a rigid translation changes nothing). Against the reference (def2-SVP): TPSS water 341.331 / 31.993 (reference 341.324 / 31.994), r2SCAN water 343.049 / 31.719 (343.045 / 31.718), TPSS formaldehyde C 5.640, O -401.680 (5.654, -401.663), M06-2X water H 31.567 (31.567), O 346.041 (346.125; M06-2X is grid sensitive). Leaving the current terms out (a common simplification) shifts the oxygen of water by about 1.4 ppm. VV10 (wB97M-V) depends only on rho and grad rho, so like a GGA it enters through the London phase of its potential on the NLC grid (water wB97M-V 348.424 / 31.683, reference 348.424 / 31.683).
`opts%nmr`; the GIAO integrals of the magnetic field and the nuclear moment come from libcint (`integrals_giao_1e`, `integrals_nmr_nucleus`,
`integrals_giao_jk`, the exact exchange of the antisymmetric first-order densities `integrals_k_general`). RIJCOSX (the default `J = 'RI'`, `K = 'cosx'`)
is supported: the Coulomb part of the field derivative is fitted (`integrals_giao_ri_j`: (ig mu nu|P) from libcint's `int3c2e_ig1` with the fitted
density of the SCF's auxiliary basis), the exchange part is COSX (`integrals_giao_cosx_k`): the London phase of the grid-side pair is a function of the
point, so K^B = K1 - K1^T is assembled from nine raw COSX builds whose grid-side function is scaled by a Cartesian component of the point
(`cosx_set_bra_kind`; no overlap fit for those), and the response exchange of the antisymmetric densities is the raw COSX build
(K[d] = (M^T - M)/2). Both use the hi COSX grid (`cosx_set_resp_grid`) with a screen proportional to the size of the densities
(`cosx_set_resp_relscreen`; an absolute screen stalls the conjugate gradients near 1e-7). The SCF of `runtype nmr` re-converges on the hi COSX grid
after the lo-grid convergence (`SCFcycle(..., hi_reconverge)`): the lo-grid density moved the oxygen of HF water by 0.12 ppm. Any other J/K choice falls
back to exact integrals. Against exact integrals (water, def2-SVP): HF O 347.717 (347.696), PBE0 345.660 (345.644), M06-2X, CAM-B3LYP, TPSS, wB97X,
wB97M-V within 0.06 ppm, hydrogens within 0.005 ppm; the remaining 0.02 ppm is the RI-J fit.
The coupled-perturbed equations have a symmetric positive definite operator (preconditioned conjugate gradients, the three field directions
sharing each exchange build). Output `[NMR]` and `nmr` in the JSON (`atoms[].isotropic`, `anisotropy`, `diamagnetic_isotropic`,
`paramagnetic_isotropic`, `tensor`, `diamagnetic`, `paramagnetic`). The anisotropy is sigma_3 - (sigma_1 + sigma_2)/2 of the principal values of
the symmetric part, ordered by absolute value. Unrestricted references: the first-order density is imaginary, so neither the Coulomb nor the semilocal kernel couples
the spins; each spin solves the closed-shell equations with D' = 2 P_sigma (its own exchange, d_oo and, for meta-GGAs, its own current kernel
v_tau,sigma / rho_sigma), J^B and the diamagnetic term take the total density and `nmr_xc_giao` evaluates the functional spin-polarized
(G_a = 2 v_sigma_aa grad rho_a + v_sigma_ab grad rho_b). Methyl radical, def2-SVP, against the reference: UHF C 127.920 / 97.087 (127.916 /
97.092), H 28.094 / 28.074 (same); UKS PBE0 C 108.835 (108.837), TPSS C 114.505 (114.504), CAM-B3LYP C 101.262 (101.250), RIJCOSX PBE0 C
108.836 / 119.413 (108.833 / 119.429), hydrogens within 0.001. Against a published-program reference (HF, def2-SVP and def2-TZVP, water, NH3, CH4, HF, C2H4,
HCN, H2CO, acetaldehyde, methanol, fluoroethane, benzene): isotropic values within 0.014 ppm, anisotropies within 0.042 ppm (the carbonyl oxygen
of formaldehyde, -458 ppm, included). Against the reference (def2-SVP; H2CO, NH3, HF, CH4, C2H4): PBE and PBE0 isotropic values within 0.03 ppm, CAM-B3LYP (H2CO, NH3, HF, water) within 0.003 ppm (B3LYP water: the VWN variant differs, 0.09 ppm at O). Environment: external point charges enter the London phase of the nuclear attraction as ghost nuclei (`atm_nuc`), the converged
solvent surface charges (COSMO/C-PCM, SMD electrostatics) the same way (`integrals_giao_charges`); neither responds (the first-order
density vanishes). Effective core potentials: a semilocal projector at centre C is given the London phase of C,
-1/2 [(R_mu - C) x <mu|r U|nu> - (R_nu - C) x <mu|U r|nu>] (origin independent; `<r mu|U|nu>` from `ECPscalar_rbra`, added to the
vendored ECP kernel - the kernel's own ignuc uses <mu|r U|nu> on both sides). Against the reference (def2-SVP): water PBE with two
point charges O 348.039 (348.037), in C-PCM water O 352.186 (352.187); CH3I HF C 195.199, I 460.660, H 30.155 (identical).
Not implemented: open-shell references.

---

## Spin-spin couplings (`runtype = 'nmr'`, `jcoupling = .true.`)

Indirect nuclear spin-spin coupling constants J (Hz) of a closed-shell molecule, Hartree-Fock or density functional, exact integrals or RI-J with COSX (`src/nmr/mod_jcoupling.f90`). The four Ramsey terms are computed separately and added: Fermi contact (FC), spin-dipole (SD), paramagnetic spin-orbit (PSO) and diamagnetic spin-orbit (DSO). `jcoupling_terms` selects them: `'all'` (default) or any of `fc`, `sd`, `pso`, `dso` separated by commas (`'fc'` is the cheap Fermi-contact-only run: one triplet response per nucleus instead of ten right-hand sides). The Jensen J-coupling basis sets are included as `baselabel = 'pcj0'` ... `'pcj4'` and `'augpcj0'` ... `'augpcj4'` (H to Ar; their s functions are dense at the nucleus, which the Fermi-contact term needs; the def2 sets are not dense enough there for 1H-1H). `jcoupling_atoms` restricts the nuclei that take part: `'all'` (default) or atom numbers and ranges such as `'1,3-5'`; the response equations are solved only for those (cost is linear in their number), and nuclei without a gyromagnetic ratio in the table (not 1H, 7Li, 13C, 15N, 17O, 19F, 29Si, 31P, 33S, 35Cl) are always left out.
FC and SD are triplet responses with the (A+B) operator of the triplet TDDFT (SD: the traceless tensor (3 r r^T - r^2 1)/r^5, six components per nucleus; the isotropic coupling is the sum over the nine tensor components, the FC/SD cross term has no isotropic part). PSO is the singlet response of the purely imaginary operator r x grad / r^3 with (A-B): a functional without exact exchange has (A-B) = orbital-energy difference, so no iteration is needed; the xc kernel of an LDA/GGA does not respond to an imaginary perturbation, a meta-GGA responds through the paramagnetic current (the same kernel as the shielding). DSO is an expectation value, integrated on the xc grid (a grid is built when the functional has none). The SD tensor integral is the second derivative of the nuclear-attraction matrix with respect to the nucleus plus (4 pi/3) delta_ab phi_mu(R) phi_nu(R) (the principal value of the dipole operator). The response uses the Galerkin subspace solver `jc_solve` (the COSX exchange is not symmetric); a negative eigenvalue of (A+B) (triplet instability, Hartree-Fock on a multiple bond) stops the run.
Output `[JCOUPLING]`: per atom pair J_total, the reduced coupling K and the four terms; `nmr.jcoupling.pairs[]` in the JSON (`J`, `K_reduced_1e19_T2_per_J`, `FC`, `SD`, `PSO`, `DSO`). Isotopes: 1H 13C 15N 17O 19F 29Si 31P 33S 35Cl 7Li. Tests: `test/jcoupling/run_tests.py` (each term against an independent code: H2, HF, H2O with Hartree-Fock, ethylene PBE0; FC 0.01 Hz, SD/PSO/DSO 0.05 Hz, a hybrid functional's FC 0.5 Hz because the xc-kernel grids differ).

## Hole and electron of the TDDFT states

`tddft_cube = 'ct'` (or any `hole:`, `elec:`, `diff:`, `trans:` entry) ends the `[TDDFT]` block of a closed-shell singlet or unrestricted TDDFT run with a table (and `tddft.states[].hole_electron` in the JSON). It is off by default: it costs about as much as the response itself on a drug-sized molecule.

```
state   D_CT (A)   q_CT (e)  sigma_h (A)  sigma_e (A)   Lambda
```

from the unrelaxed hole and electron densities of the response vectors, rho_h = sum_ij P^h_ij psi_i psi_j with
P^h = 1/2 [(X+Y)^T (X+Y) + (X-Y)^T (X-Y)] and rho_e = sum_ab P^e_ab chi_a chi_b with P^e = 1/2 [(X+Y)(X+Y)^T + (X-Y)(X-Y)^T] (spin summed,
both integrate to sum (X^2 + Y^2), 1 for a pure single excitation), evaluated on a regular grid (spacing 0.35 bohr or coarser,
points further than 6 bohr from every nucleus left out): `D_CT` the distance between the centroids of the normalized hole and
electron (Angstrom), `q_CT` the integral of the positive part of rho_e - rho_h (electrons), `sigma_h`, `sigma_e` the root mean
square extent of each around its centroid (Angstrom), `Lambda` the overlap of the normalized densities, the integral of
sqrt(rho_h rho_e) (about 0.02 uncertain on this grid, the others 1e-3). Small `Lambda`, large `D_CT` and `q_CT` near 1 mark a
charge-transfer state. Against pyscf's response vectors on the same molecule (formaldehyde, PBE0/def2-SVP): D_CT, sigma and q_CT to 1e-3,
the densities to 2e-4 (hole, difference) and 3e-3 (transition density) relative rms at grid points; the first moment of the
transition density equals the transition dipole. The cost is a few grid sweeps (40 atoms, 10 states: 16 s on 16 threads, a quarter
of the TDDFT run); without a request it is not computed. An unrestricted reference gives the same table (the densities are the sums over the two
spins, the response vectors being normalized together; against pyscf on the CH3 radical D_CT and sigma agree to 1e-3, q_CT to 0.006,
Lambda to 0.01); no triplets.

`cube = 'diff:1, hole:2, elec:2, trans:1'` (or `tddft_cube`) also writes the densities as cube files `<base>.s<state>.<kind>.cube`
on a 0.2 bohr grid (electron minus hole, hole, electron, transition density). With `opt` or `irc` they belong to the last frame.

---

## Cube files (`cube = 'den, esp'`)

`opts%cube_den` / `opts%cube_esp` / `opts%cube_esp_exact` (input key `cube` in `&molecule` / `[molecule]`: `'den'`, `'esp'`, `'esp_exact'` or a list; with RI-J on `'esp'` is the potential of the RI-J fitted density, `'esp_exact'` the exact one) write
`<stem>.den.cube` (electron density, e/bohr^3) and `<stem>.esp.cube` (electrostatic potential, Eh/e; `<stem>` is `opts%freq_base`)
after a converged single-point calculation. The grid extends 6 bohr beyond the atoms with 0.2 bohr spacing (coarser when the
box would hold more than 3e6 points). The density is the occupied-orbital form of the total density matrix at the grid
points, the potential `esp_at_grid_batch` plus the point nuclei (effective charge Z - ECP core electrons, cut off 0.05 bohr from
a nucleus). Against pyscf on the same molden file and random grid points: density 6e-6 relative, potential 3e-5 absolute
(H2O, CH3 UKS, CH3I with an ECP).

`cube = 'den, esp, all'` adds the cube files of every step: `<base>.c<cycle>.den.cube` for each optimization cycle, `<base>.ts`, `<base>.f<k>` and
`<base>.b<k>` (forward / backward point k) for an IRC (one more single point per point), the last frame keeping `<base>.den.cube`; a note on the
cost is printed at the start (tens of MB and seconds per file for a drug-sized molecule). With `opt` or `irc`, `cube` and `tddft_nroots` apply to the last frame only: after the optimization converges (before the
frequencies if `freq` is also set), or at the last frame of `<base>.irc.xyz` (the forward end point; the backward end
for a backward-only run). The optimization steps and IRC points run without TDDFT and without cube files.

---

## Citations (`[CITATIONS]`, `<base>.bib`)

At the end of the output file: `[CITATIONS]` ... `[CITATIONSEND]`, listing the methods this run used (grouped:
integrals, grid, COSX, density fitting, basis sets, dispersion, non-local correlation, composite methods, solvation, SCF,
geometry optimization / TS search / IRC, population analysis) and, for the exchange-correlation part, the Libxc version, its
reference and each selected functional's references, all read from Libxc's API. `<base>.bib` has the same entries in BibLaTeX.
The table is `src/harness/mod_citations.f90`: a new method with a literature reference gets one `add(...)` line and
one `use_ref(...)` condition there.

---

## Intrinsic reaction coordinate (`runtype = 'irc'`)

The input geometry is a transition state (one imaginary mode); the path is followed down both sides in
mass-weighted coordinates and written to `<base>.irc.xyz` (backward end ... TS ... forward end) and
`<base>.irc.dat` (s, E, E - E(TS), gmax). The output file has `[IRC]` (settings), `[IRCSTART]` (TS energy, lowest
mode, imaginary-mode count), one `[IRCWALK forward|backward]` table per direction ending with `[IRCWALKEND]`,
`[IRCRESULT]` (points, end energies, barriers, why each direction stopped; ends with `[IRCRESULTEND]`) and
`[IRCFILES]`. With `corrector` (default on) each predictor point is followed by a corrector on a fitted surface. Stops: gradient (`max_force`), Mayer bond orders no longer changing (`bond_order_change`), the quadratic
model's minimum, `max_points`, or a failure. Settings and examples: `examples/irc/README.md`.

---

## Static polarizability (`runtype = 'polar'`)

After the SCF the output file gets `[POLARIZABILITY]`: the symmetric 3 x 3 tensor alpha_ij = -d2E/dF_i dF_j
in atomic units (4 decimals), the isotropic value (xx+yy+zz)/3, the anisotropy (both also in A^3) and the
number of coupled-perturbed iterations; the block ends with `[POLARIZABILITYEND]`. Closed shell only. The response
uses the SCF's own Coulomb and exchange operators (RI-J, COSX) and the xc kernel on the SCF grid, converged to 1e-8.
For VV10 functionals the SCF includes VV10 but the response omits its nonlocal kernel; the block then ends with a note
line saying so (water wB97M-V/def2-SVP: -0.036 % of alpha; other systems not quantified).

With a solvent model (COSMO or SMD) the static tensor is the second derivative of the energy with the solvent
equilibrated: the response products carry the ground-state surface-charge operator with f(eps), as the frequencies do; the
dynamic tensors use only the electronic part of the solvent, f(eps_opt) with `tddft_eps_opt` (default 1.78), as the vertical
TDDFT excitation. The SMD CDS term does not depend on the density and is not part of the response. Water PBE0/def2-SVP
(`convention = 'cpcm'`, eps 78.4) against the finite-field second derivative of the pyscf C-PCM energy: xx/yy/zz
3.3353 / 7.6189 / 5.5979 against 3.3358 / 7.6195 / 5.5985 a.u. (gas phase 3.1256 / 7.0262 / 5.2587 in both); regression
test `test/polar_solvent` (all checks below are in it). SMD against pyscf SMD, same finite-field test: 3.4450 / 7.7274 / 5.7405
against 3.4480 / 7.7353 / 5.7447 (the two SMD cavities differ; total energies 9e-5 Ha apart). With eps_opt = eps and a small
omega the dynamic tensor equals the static one. Dynamic tensors (omega 0.05 and 0.1, PBE and PBE0, gas and solvent) against the
sum over all 95 TDDFT states of the same run: within 7e-4 a.u. (4-decimal output). RI-J (+4e-3 a.u. in the isotropic value, the
same in gas phase and solvent) and COSX (-2e-4 more) do not interact with the solvent term. The field acts on the solute only
(no cavity-field factor).

External point charges are allowed too (they change only the one-electron Hamiltonian): water PBE0/def2-SVP with two charges,
diagonal 3.1087 / 6.9738 / 5.2159 and xz -0.0485 against the pyscf finite-field values (identical to 4 decimals); with C-PCM
3.3160 / 7.5530 / 5.5481 against 3.3164 / 7.5535 / 5.5486; the dynamic tensor with the charges against the sum over states: 2.6e-4.

`polar_omega = '0.05, 589nm, 2.1ev'` (a.u., or with nm / ev units, up to 16 values) adds the dynamic polarizability
alpha(omega) at each frequency, written below the static tensor and to `polarizability.dynamic` in `<base>.json`.
Pure functionals solve [(A+B) - w^2 (A-B)^-1](X+Y) = -2v by PCG; hybrids and range-separated hybrids (K='exact'
or 'cosx') solve the coupled pair (A+B)U - wW = -v, (A-B)W - wU = 0 in two subspaces, with the kernel on the TDDFT grid
(`ENGINE_POLAR_GRID=scf`: the SCF grid; caffeine/def2-SVP alpha changes by <= 1.2e-3 a.u., the solve is 8-10 % faster)
(the static value of a range-separated hybrid with K='exact' goes through the same solver at w = 0). Frequencies at or
above the HOMO-LUMO orbital gap are rejected, and so are those beyond the lowest dipole-active excitation: the hybrid
solver compares omega with the lowest root of its converged spaces, the pure-functional CG stops on negative
curvature (water PBE0: 0.28 runs, 0.30 stops; its first excitation is 0.2916). Water/def2-SVP isotropic alpha at w = 0, 0.05, 0.15 against pyscf:
PBE0 5.1684 / 5.2048 / 5.5444 (pyscf 5.1684 / 5.2048 / 5.5444), CAM-B3LYP 5.1828 / 5.2209 / 5.5807 (5.1828 / 5.2209 /
5.5807); with K='cosx' within 3e-4.
Each tensor is followed by `Wall time (s):` for its solve. Caffeine/def2-SVP, 24 threads, RIJCOSX: PBE0 static 13.5 s,
w = 0.05 13.1 s, w = 0.15 16.0 s; CAM-B3LYP 18.0 / 20.1 / 23.7 s; PBE 5.5 / 5.4 / 5.9 s (w = 0.05 / 0.10).

---

## Excited states (`tddft_nroots = n`)

`tddft_nroots = 5` (in `&molecule`, or `[molecule]` in TOML) adds the lowest five singlet excited states after the
SCF: full linear-response TDDFT (not TDA), adiabatic kernel, closed shell, pure functionals only (LDA, GGA,
meta-GGA, global and meta-GGA hybrids such as B3LYP, PBE0, M06-2X and R2SCAN0, and range-separated hybrids such as CAM-B3LYP and wB97X,
with K='cosx' or K='exact'). With `opt` or `irc` it is run on the last frame only (see "Cube files"); external point
charges (`&pointcharges`) are allowed (they enter only the one-electron Hamiltonian, so the orbital energies carry them; checked against pyscf's `qmmm.mm_charge` TDDFT, water PBE0/def2-SVP 5 states within 2e-4 eV, OH radical UKS within 3e-4 eV; with RI-J + COSX within 5e-3 eV of exact J/K; together with C-PCM water 3-4e-3 eV from pyscf, the cavity difference); a COSMO/SMD solvent is treated as a non-equilibrium vertical excitation (below). The output file gets `[TDDFT]`: per state the energy (eV, Ha), wavelength (nm), oscillator
strength f and transition dipole (a.u.), then the main occ -> vir pairs (weight X^2 - Y^2 >= 0.05, orbitals
numbered from 1); `<base>.json` gets `tddft.states`. The xc kernel is integrated on one TDDFT grid for every functional (radial 0.4 / 0.8 of the M3 counts for H /
other rows, zone orders 14/26/86/110/86, H peak 50; imatinib 129k points; `ENGINE_CPSCF_XC_GRID=full`
uses the SCF grid, `=lo` the theory-table LO tier). Kernel-grid error against the SCF grid, max over 5 states:
caffeine PBE 3e-7, PBE0 7e-7, CAM-B3LYP 2e-6, TPSS 3e-6, wB97X 2.6e-5, M06-2X 4.7e-5 Ha; water (diffuse states)
PBE0 9e-6, PBE 2.2e-5, CAM-B3LYP 2.4e-5, R2SCAN0 8.9e-5, wB97X 9.8e-5, wB97M-V 1.4e-4, M06-2X 3.7e-4, M06 8.0e-4 Ha
(Minnesota kernels on water need ~29k points per atom for 1e-5). Minnesota functionals (M06 family) have their
own kernel grid (radial 0.55/1.1, zones 26/74/194/266/194, H peak 194; 156k points on caffeine, 15k on water): M06-2X
/ M06 are 9e-6 / 7.6e-6 Ha off the SCF-grid kernel on caffeine (was 4.7e-5 / 1.2e-4) and 9e-6 / 3.2e-5 on water (was
3.7e-4 / 8.0e-4). Imatinib 10 roots TDDFT: PBE 22 s, PBE0 105 s, M06-2X 186 s (was 145 s;
kernel 429k points against 129k before). The response uses the SCF's RI-J when the SCF did, which
shifts excitation energies slightly (water PBE/def2-SVP: up to 1.9e-4 Ha, 5e-3 eV, below exact J). Checked against pyscf
TDDFT with exact J/K (water PBE, BLYP, TPSS: energies to 1e-5 Ha, f to 4 decimals; benzene PBE: 12 states incl.
the degenerate pairs). Hybrids: water PBE0 and B3LYP with exact J/K to 2e-6 Ha, f to 4 decimals; with COSX
(J exact) 1.6e-5 .. 4.3e-5 Ha. Range-separated: water CAM-B3LYP exact J/K 9e-6 Ha; wB97X 5e-5 Ha on the LO kernel grid,
1.7e-5 Ha with ENGINE_CPSCF_XC_GRID=full (B97-type kernels are more grid-sensitive). Meta-GGA hybrids (no code of
their own: the kernel handles tau): water M06-2X / R2SCAN0 / M06 / wB97M-V against pyscf 1.2e-5 / 1.4e-5 / 6e-6 /
4e-6 Ha on the SCF grid, up to 6.8e-5 Ha (M06) on the LO kernel grid. Imatinib M06-2X 10 roots: TDDFT 207 s
(our LO grid for the M06-2X kernel: 412k points). VV10: the ground state includes it, the response kernel does not.
Davidson (residual 1e-5 on Omega = D^1/2 (A+B) D^1/2) can miss a low state that is pushed up strongly by the
Coulomb coupling in a small trial space: benzene PBE with `tddft_nroots = 4` gives the three lowest states and then
7.25 eV, not the bright E1u pair at 7.19 eV (`tddft_nroots = 10` finds it).
`tddft_robust = .true.` tracks max(n, 6) spare roots, which finds it, at about 1.5x the response vectors.

Open-shell references (UKS, any multiplicity) are supported: one list of states, each main transition
labelled with its spin (`from a -> to a`, JSON `"spin"`), f = 2/3 w |mu|^2 with mu summed over both spins; the Davidson
vectors are the stacked (alpha, beta) vir-occ amplitudes and the products come from `cpscf_apb_amb_uhf` (RI-J or
exact J, COSX or exact K, range-separated hybrids, kernel from `xc_fxc_c_u`). Against pyscf UKS/UHF TDDFT on NH2 (doublet) and CH2
(triplet), def2-SVP, exact J/K: HF 1e-6 Ha, LDA 2e-6; on the classical grids (retired 2026-10-21) NH2 PBE
1e-6 and TPSS 6e-6, CH2 PBE 8e-5, PBE0 6e-5, CAM-B3LYP 3e-5, TPSS 2.7e-4 (largest over 6 states; pyscf's own level 5 and
9 grids differ by 1e-5 in the first CH2 state). On the default (theory-table) grid the open-shell frontier orbital energies are within 1e-5 Ha for PBE/B3LYP and
4e-5 for r2SCAN (see docs/ENVIRONMENT.md, "XC grids of open-shell references"); the first CH2 excitation is
4e-5 Ha (PBE) to 1.3e-4 (TPSS) off pyscf's level-9 result. `tddft_triplets` is for closed shells only. A nearly
degenerate pair (slightly bent OH, 7e-4 Ha) is found, with the COSX error (4e-4) of the same size as the energy.

**Stability of a closed-shell reference.** `stability = 1` (namelist, TOML, or ENGINE_STABILITY) checks after the SCF,
besides the old restricted -> restricted test (finite-difference orbital Hessian, Lanczos, 40 Fock builds), the
restricted -> unrestricted one: the lowest eigenvalue of the triplet-type (A+B), from the TDDFT products themselves
(Davidson, 5-8 steps, SCF grid kernel). A negative value (< -1e-3) means the restricted determinant is a saddle point
toward a lower unrestricted one; the triplet excitation energies of such a reference are imaginary.
`stability = 2` then redoes the run unrestricted (`uks = .true.`: alpha and beta orbitals of a closed-shell molecule are
separate) and its own stability step follows the negative mode to the symmetry-broken solution; TDDFT runs on that
reference. H2 at 2.0 A, def2-SVP, exact J/K: triplet-type eigenvalue -5.77e-2 (PBE) and -9.46e-2 (PBE0) against
-5.76e-2 / -9.48e-2 from the unrestricted orbital Hessian; symmetry-broken energies -1.009740 / -1.011179 Ha against
pyscf's broken-symmetry UKS -1.009738 / -1.011178 (RKS: -0.995619 / -0.984744). Water and caffeine (PBE, PBE0): stable.
Every driver (single point, optimization, IRC, Hessian) redoes the call unrestricted, starting from the density the
restricted pass ended with, and keeps `uks` on for the remaining steps of an optimization or IRC. The output keeps only
the unrestricted pass (the restricted pass's SCF listing is cut off; a line states the eigenvalue that caused the
restart). H2 from 2.0 A optimized with `stability = 2` ends at 0.761 A. `uks = .true.` alone gives the unrestricted
solution that SCF reaches from the (symmetric) guess; add `stability = 2` to break the symmetry.

`tddft_tda = .true.` solves the Tamm-Dancoff problem A X = w X instead of the full one (X+Y = X-Y = X, sum X^2 = 1; every
kind of state: singlets, triplets, open shells; pure, hybrid, range-separated, meta-GGA). Water/def2-SVP against pyscf
`TDA` (exact J/K, SCF-grid kernel): PBE/PBE0/CAM-B3LYP 2e-6, TPSS 1e-5 Ha, singlets and triplets; NH2/CH2 UKS 1e-5
(CH2 second state 7e-5, grid). Caffeine/def2-SVP RIJCOSX, 10 roots: TDDFT part PBE 2.9 s (full 2.9), PBE0 7.2 s (full 8.6).
Energies differ from the full ones by 0.006-0.15 eV (PBE) and 0.02-0.19 eV (PBE0) on the first states, f of the first
PBE0 state 0.163 against 0.142: TDA is an approximation of the linear response, not a free speed-up. The output
says which method was used (`Method:` line, `"method"` in the JSON).

A full problem with an imaginary excitation energy (w^2 < 0: the restricted determinant is a saddle point toward a
lower unrestricted one; typical for the triplets of stretched bonds, diradicals) used to stop the run. It is now
detected (pure functionals: lowest Ritz w^2 < 0; hybrids: a negative 1/w^2 larger in size than the nroots-th
positive one), reported with the lowest w^2 found and that block is redone in the Tamm-Dancoff approximation, whose
roots are real. H2 at 2.0 A, def2-SVP: PBE lowest triplet 0.0182 Ha (pyscf TDA 0.01820; pyscf's full TDDFT silently drops
it and starts at 0.6517), PBE0 0.00247 (pyscf TDA 0.00247). The ground state itself is not re-optimised.

`tddft_triplets = .true.` adds the lowest n triplet states of the same closed-shell reference after the singlets: a
second `[TDDFT]` block headed "Triplet excited states" (f and transition dipole 0) and `tddft_triplet` in
`<base>.json`. A+B loses the Coulomb term and takes the kernel f_aa - f_ab (the open-shell kernel at rho_a = rho_b),
A-B is the singlet one. Water/def2-SVP, 5 roots, exact J/K, SCF-grid kernel, against pyscf (`singlet=False`):
PBE 3e-6, PBE0 2e-6, CAM-B3LYP 2e-6, TPSS 8e-6 Ha. Caffeine PBE0 RIJCOSX, 10 roots: singlets 9.1 s, triplets 8.1 s.
Davidson convergence: residual 1e-4 (energies agree with a 1e-5 run to 1e-7 Ha, f within 1e-4). With exact
exchange X+Y and X-Y have their own trial spaces, and each X+Y / X-Y pair shares one exchange density (the
symmetric and antisymmetric parts of its K); `ENGINE_TDDFT_SOLVER=single` uses one space for both instead.
Imatinib (673 functions, 10 roots, 24 threads, RIJCOSX): CAM-B3LYP 164 s for the TDDFT part.



---

## Frequencies, IR intensities and thermochemistry (`freq=.true.`)

Input key `freq = .true.` (`&molecule` namelist or TOML) runs the analytic Hessian after the
SCF and writes, to the output file at any `verbose` level:

- `[FREQUENCIES]` (cm-1; translation/rotation
  modes first, imaginary modes appear as negative values among the vibrations);
- `[NORMALMODES]`: groups of three modes, each with frequency, IR intensity (km/mol), reduced mass,
  force constant and the unit-norm Cartesian displacement of every atom (dx dy dz, 4 decimals); `[IR]` summary;
  `[DIPOLEDERIV]`: dipole derivatives dmu_c/dR (e), one row per atom and dipole component;
  the Cartesian Hessian (Hartree/bohr^2) is only in `<base>.hess`, not in the output file;
- `[THERMO]` (conditions `freq_temp` K, `freq_pres` atm, `freq_conc` mol/L - the last replaces the pressure by a concentration
  standard state; `opts%freq_temp`, `opts%freq_pres` in Pa): ZPE, thermal E/H/G corrections and total/translational/rotational/vibrational/electronic
  entropy (quasi-RRHO) at 298.15 K, 1 atm (ideal gas, rigid rotor, harmonic oscillator; rotational symmetry number
  found automatically, electronic term R ln(multiplicity)); sums with the SCF energy;
- `<base>.modes.csv` (per-mode, per-atom displacements for plotting) and `<base>.freq.molden`, where `<base>` is the output file name without `.out` (in-process EngineUp callers that don't set it get `direwolf`).

With `runtype='opt, freq'` the Hessian is evaluated once at the converged geometry. IR intensities use the
CPSCF response already computed for the Hessian (no extra solve) plus the dipole Pulay term.
External point charges are supported; they break the translational and rotational invariance, so no mode is projected
out then (all 3N modes are vibrations, the hindered translations and rotations included). Not supported / not applicable: VV10
functionals (nonlocal second derivative omitted), K='RI'. Frequencies are harmonic; the vibrational entropy uses Grimme's quasi-RRHO
(reference 100 cm-1).
The symmetry-number tolerance is loose (0.2 bohr), so slightly distorted geometries may get a higher
symmetry number than the exact point group.

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

### TDDFT in a solvent (COSMO / SMD)

With a solvent model on, the excitation is treated as vertical (non-equilibrium): the ground state is converged with the
static dielectric constant, and only the electronic part of the solvent follows the transition density. The (A+B) products
add the ground-state CPCM operator with f(eps_opt) = (eps_opt-1)/(eps_opt+0.5) in place of f(eps); `tddft_eps_opt` (default
1.78, water's n^2) sets eps_opt. Triplets of a closed-shell reference have no electrostatic response (no solvent term); UKS
references use the response of the total transition density. Formaldehyde PBE0/def2-SVP in water against pyscf C-PCM
`equilibrium_solvation=False`: excitation energies 4.0092/8.8477/9.1482 against 4.0030/8.8298/9.1462 eV (gas phase
3.8928/8.5213/8.9569 against 3.8929/8.5212/8.9569); the residual is the cavity difference, not the response. The SMD CDS term
belongs to the ground state and is not part of the excitation energy.
