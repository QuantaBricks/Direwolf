# Excited states (TDDFT)

Vertical singlet excitation energies, wavelengths and oscillator strengths after the SCF (`tddft_nroots = n`, closed
shell or unrestricted reference; full linear response, `tddft_tda = .true.` for the Tamm-Dancoff approximation). The
reference is in `docs/INTERFACE.md`, "Excited states". Run with `./Direwolf <file>`.

- `formaldehyde_uv.inp` - PBE0/def2-SVP, RIJCOSX, 5 states.
- `formaldehyde_uv_water.inp` - the same in water: a `&cosmo` block switches the solvent on; the excitation is treated as
  vertical (non-equilibrium), the optical response using the optical dielectric constant `tddft_eps_opt` (default 1.78).
  `convention = 'cpcm'` selects the C-PCM convention (f = (eps-1)/eps, hydrogen radius 1.10 A) that the other quantum
  chemistry programs use; without it the default `'cosmo'` convention is used.
- `formaldehyde_opt_uv.inp` - `runtype = 'opt'` with `tddft_nroots` and `cube`: the optimization runs without TDDFT, and
  the excited states and the cube files (`<base>.den.cube`, `<base>.esp.cube`) are for the optimized geometry only.

- `formaldehyde_hole_electron.inp` - the same states with `tddft_cube`: cube files of the electron-minus-hole density of
  state 1 and of the hole and electron of state 2 (`formaldehyde_hole_electron.s1.diff.cube`, `.s2.hole.cube`,
  `.s2.elec.cube`; show them as isosurfaces in any cube viewer).

## Output

`[TDDFT]` in the output file lists, per state, the energy (eV and Ha), the wavelength (nm), the oscillator strength f and
the transition dipole (a.u.), followed by the main occupied -> virtual pairs. `<base>.json` carries the same in
`tddft.states` (`energy_ev`, `wavelength_nm`, `oscillator_strength`, `transition_dipole_au`, `transitions`).

Formaldehyde, PBE0/def2-SVP, gas phase (`formaldehyde_uv.inp`):

```
state   E (eV)     nm        f
  1      4.0040   309.65   0.0000     n -> pi*, symmetry forbidden
  2      8.5191   145.54   0.1516
  3      9.1641   135.29   0.0013
  4      9.7027   127.78   0.0057
  5     10.3540   119.75   0.0000
```

In water (`formaldehyde_uv_water.inp`) the n -> pi* state moves to 4.1049 eV (+0.10 eV, the usual blue shift) and the
bright state to 8.8044 eV.

Each singlet state also gets a hole / electron line (distance of the centroids `D_CT`, moved charge `q_CT`, extents
`sigma_h`, `sigma_e`, overlap `Lambda`; see `docs/INTERFACE.md`): for formaldehyde

```
state   D_CT (A)   q_CT (e)  sigma_h (A)  sigma_e (A)   Lambda
    1      0.365      0.741        1.209        1.186    0.508     HOMO -> LUMO (n -> pi*)
    2      1.449      0.811        1.229        1.924    0.433
```

A spectrum is the sum of one broadening function per state, with a height proportional to f; the output is the stick
spectrum, and drawing the broadened curve is left to the plotting tool.

- `formaldehyde_nto.inp` - `nto = '1, 3'` in `&tddft`: natural transition orbitals of states 1 and 3. Each state's file `<base>.s<state>.nto.molden` holds the hole orbitals H1, H2, ... (occupation = weight of the pair) and the particle orbitals P1, P2, ...; the pair weights are also listed in the `[NTO]` block of the output. The pairs with weight >= 0.01 also get cube files (`<base>.s<state>.nto.h1.cube`, `.p1.cube`, ...: orbital values) for any cube viewer. Open the molden file in Multiwfn, Avogadro or any orbital viewer.