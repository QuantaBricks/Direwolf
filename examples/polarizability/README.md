Static dipole polarizability (`runtype = 'polar'`), water PBE0/def2-SVP, RIJCOSX:

- `water_polar.inp` - native format; the output file gets a `[POLARIZABILITY]` block with the
  3 x 3 tensor alpha_ij = -d2E/dF_i dF_j (a.u., 4 decimals), the isotropic value (xx+yy+zz)/3 and the
  anisotropy, both also in A^3.

Run with `./Direwolf water_polar.inp`. Closed-shell HF, global hybrids and (meta-)GGA with J = exact or
RI and K = exact or cosx are supported. Not implemented: open shell, range-separated hybrids with exact
exchange, solvent models, point charges, and `polar` together with `opt` (optimize first, then run
`polar` on the result). The coupled-perturbed equations are solved to 1e-8 (the Hessian uses 5e-5), on the
SCF grid.
