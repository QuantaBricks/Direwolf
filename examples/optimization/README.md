# Geometry optimization examples

An `&opt` namelist (its presence alone is enough) switches a run from a
single point to a geometry minimization. Run these from the repository
root, e.g. `./Direwolf examples/optimization/water_opt.inp`.

Every cycle is appended to `<output>.opt.xyz` (multi-frame XYZ); the
converged geometry is also printed in the `.out` as `[FINAL GEOMETRY]`.
`opt_restart` (default on) resumes from the last frame of an existing
`<output>.opt.xyz`, so delete that file to start fresh.

| file | shows |
|---|---|
| `water_opt.inp` | the minimal case - `&opt` with no fields, redundant internal coordinates, Lindh model Hessian |
| `water_opt_cart.inp` | `opt_coord = 'cart'` - optimize in Cartesians instead |
| `water_opt_tight.inp` | `opt_conv = 'tight'` plus `opt_maxcyc` / `opt_trust`; tight also pulls the SCF up to `'tight'` automatically |
| `water_opt_xtbhess.inp` | `opt_hessian_file` - seed from an external Cartesian Hessian (`water_gfn2.hessian`, produced by `xtb w.xyz --hess`) instead of the model |
| `paracetamol_opt.inp` | a realistic 20-atom organic: PBE/def2-SVP, RI-J + COSX, 8 threads |

`&opt` fields (all optional): `opt_maxcyc` (100), `opt_conv`
(`'normal'` | `'tight'`), `opt_trust` (0.3 bohr), `opt_coord`
(`'ric'` | `'cart'`), `opt_restart` (`.true.`), `opt_hessian_file`
(`''`), `opt_write_ric` (`.false.`). See the top-level README's
"Geometry optimization" section for the convergence-threshold values.

To make your own seed Hessian for `opt_hessian_file`: any code that
writes a 3N x 3N Cartesian Hessian in Hartree/bohr² works;
`xtb <geom>.xyz --hess` writes one (`$hessian`, 5 values per line) that
is read as-is.
