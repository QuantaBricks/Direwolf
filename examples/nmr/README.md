# NMR shielding tensors (`runtype = 'nmr'`)

Magnetic shielding tensors of every nucleus with gauge-including atomic orbitals (GIAO), so the result does not depend on the gauge
origin. Run with `./Direwolf water_nmr.inp`.

- `water_nmr.inp` - water, Hartree-Fock / def2-SVP.
- `benzene_nics.inp` - benzene, PBE0 / def2-SVP, NICS of the ring (`nics = 'rings'`).
- `ethylene_jcoupling.inp` - ethylene, PBE0 / def2-SVP, shieldings and spin-spin couplings (`jcoupling = .true.`).

RIJCOSX (the default) and exact integrals (`J = 'exact'`, `K = 'exact'`) are supported; other choices fall back to exact integrals (a notice is printed). Hartree-Fock and all LDA/GGA/meta-GGA functionals: pure, global hybrid, range-separated hybrid, and VV10 (e.g. wB97M-V). Effective core potentials, solvent models (COSMO/C-PCM, SMD) and external point charges are supported. Closed-shell and
unrestricted (open-shell, `imult` > 1 or `uks = .true.`) references.

## Output

`[NMR]` lists, per atom, the isotropic shielding, the anisotropy (sigma_3 - (sigma_1 + sigma_2)/2 of the principal values ordered by
absolute value) and the diamagnetic and paramagnetic parts of the isotropic value, all in ppm. `verbose = 2` prints the full tensors,
`<base>.json` holds them for every atom (`nmr.atoms[].tensor`, `diamagnetic`, `paramagnetic`; `tensor[a][b]` = d2E/dB_a dm_b).

Water, HF/def2-SVP:

```
atom  elem   isotropic   anisotropy   diamagnetic   paramagnetic
   1    O       347.696       37.733       379.876        -32.180
   2    H        31.181       20.150        31.157          0.024
```

A chemical shift is the shielding of a reference compound minus the shielding of the nucleus, computed with the same method and basis.


## Spin-spin couplings (`jcoupling = .true.`)

With `runtype = 'nmr'`, `jcoupling = .true.` adds the indirect nuclear spin-spin coupling constants J (Hz) and the reduced couplings
K (1e19 T^2/J) of every pair of nuclei, as `[JCOUPLING]` in the output and `nmr.jcoupling.pairs` in `<base>.json`.

- **Fermi-contact term only.** The paramagnetic and diamagnetic spin-orbit and the spin-dipole terms are not included. The contact
  term dominates the one-bond couplings to hydrogen and the vicinal H-H couplings; for couplings between heavier nuclei and for
  geminal H-H the other terms matter, and the total J can differ in size and sign.
- Isotopes (gyromagnetic ratio of the most abundant spin-1/2 nucleus): 1H, 13C, 15N, 19F, 29Si, 31P, and 17O, 33S, 35Cl, 7Li. Pairs
  with other elements are not listed. The electron g factor is 2 (the value of the Ramsey derivation).
- Closed-shell molecules; Hartree-Fock, LDA/GGA/meta-GGA, hybrids, RIJCOSX (default) and exact integrals. Solvent and point charges
  enter through the orbitals.
- The triplet response needs a stable reference. Hartree-Fock on a multiple bond (ethylene, C=O) has a triplet instability: the
  run stops with a message instead of printing the values an unstable reference gives. A hybrid functional (PBE0, B3LYP) is stable there.
- The basis needs tight s functions at the nuclei for converged contact couplings (def2-TZVP, cc-pVTZ or larger; no pcJ-n sets are
  installed yet).


## NICS (`nics`)

With `runtype = 'nmr'`, `nics = 'rings'` puts probes in every ring of the molecule; `nics = 'x y z; x y z'` puts them at points in Angstrom.
A probe is a point with no nucleus and no basis function, so it changes nothing in the wave function and costs no extra response
(the magnetic response is solved once; the probes only enter the last contraction). NICS = -sigma (ppm): negative above an aromatic
ring, positive above an antiaromatic one. Output: `[NICS]` and `nmr.nics.probes` in `<base>.json`.

- `rings`: rings of 3 to 8 members found from the covalent radii (the smallest ring through each bond, so a fused system gives its
  faces). Per ring NICS(0) at the centroid, NICS(+1) and NICS(-1) 1 Angstrom above and below the mean plane, their average, and the zz
  components along the ring normal (NICS(0)zz, NICS(1)zz). Non-planar rings use the least-squares plane.
- Points: the isotropic value and the component along the lab z axis.
- Benzene, PBE0/def2-SVP: NICS(0) -8.87, NICS(1) -11.24, NICS(1)zz -29.50; cyclobutadiene NICS(0) +28.19. A probe on a
  nucleus gives exactly -sigma of that nucleus.
