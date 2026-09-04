# RESP charges

On by default (`resp_charges_on = .true.`) for any converged run - the
per-atom table in the output gains a third `RESP_q` column alongside
Mulliken/Lowdin, no separate table. Set `resp_charges_on = .false.` in
`&molecule` to skip the fit if you don't need it (see "Cost" below - it's
now cheap enough that this is mainly about a tidier two-column table, not
saved time). The implementation follows Psi4's RESP plugin
([cdsgroup/resp](https://github.com/cdsgroup/resp), BSD-3, which is what
Psi4 users actually run - Psi4 core itself ships no RESP) and through it
the original RESP paper ([Bayly, Cieplak, Cornell, Kollman, *J. Phys.
Chem.* **1993**, 97, 10269](https://doi.org/10.1021/j100142a004)),
matching its defaults exactly:

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
