Harmonic frequencies, IR intensities and thermochemistry (`runtype = 'freq'`),
all on water with B3LYP/def2-SVP, RIJCOSX:

- `water_freq.inp` - native format, frequencies at the given geometry.
- `water_opt_freq.inp` - native format, `runtype = 'opt, freq'`: optimize, then
  frequencies at the converged geometry.
- `water_toml_freq.toml` / `water_toml_opt_freq.toml` - TOML format, `runtype = "freq"` and
  `runtype = "opt, freq"`.
- `water_freq.gjf` - Gaussian format, `# B3LYP/def2SVP Opt Freq`.

Run with `./Direwolf <file>`. The output file gets `[FREQUENCIES]`,
`[NORMALMODES]`, `[IR]`, `[DIPOLEDERIV]` and `[THERMO]`; the Hessian matrix is not printed there, it is in `<base>.hess`; also written: `<base>.hess`, `<base>.modes.csv` and
`<base>.freq.molden` are written next to the output (`water_freq.out` -> `water_freq.hess`, ...). See the top-level
README's "Frequencies" section for limits (no `K = 'RI'`, no point charges).
