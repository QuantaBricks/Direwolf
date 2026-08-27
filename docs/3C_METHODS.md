# 3C Methods (Composite Methods)

Three composite methods with fixed basis, dispersion, and corrections:

| Method | Functional | Basis | Dispersion | Correction |
|--------|-----------|-------|-----------|-----------|
| B97-3C | B97-3C | def2-mTZVP | D3(BJ) | SRB |
| R2SCAN-3C | R2SCAN | def2-mTZVPP | D4 | gCP |
| WB97X-3C | WB97X-V (no VV10) | vDZP | D4 | none |

**Important**: Do NOT set `baselabel` or add `-D3`/`-D4` suffixes. The basis/dispersion/corrections are all baked into the method name.

## Example 1: Water with B97-3C

```fortran
&molecule
 functional = 'B97-3C'
&end
&atoms
O    0.0000   0.0000   0.0000
H    0.7584   0.0000   0.5861
H   -0.7584   0.0000   0.5861
&end
```

## Example 2: Benzene with R2SCAN-3C

```fortran
&molecule
 functional = 'R2SCAN-3C'
 calc_force = .true.
&end
&atoms
C   -1.4008    0.8088    0.0000
C   -1.4008   -0.8088    0.0000
C    0.0000   -1.6175    0.0000
C    1.4008   -0.8088    0.0000
C    1.4008    0.8088    0.0000
C    0.0000    1.6175    0.0000
H   -2.4871    1.4340    0.0000
H   -2.4871   -1.4340    0.0000
H    0.0000   -2.8680    0.0000
H    2.4871   -1.4340    0.0000
H    2.4871    1.4340    0.0000
H    0.0000    2.8680    0.0000
&end
```

## Example 3: Alkane with WB97X-3C

```fortran
&molecule
 functional = 'WB97X-3C'
&end
&atoms
C    0.0000    0.0000    0.0000
C    1.5400    0.0000    0.0000
H   -0.5100   -0.8900    0.0000
H   -0.5100    0.4450    0.8700
H   -0.5100    0.4450   -0.8700
H    2.0500   -0.8900    0.0000
H    2.0500    0.4450    0.8700
H    2.0500    0.4450   -0.8700
&end
```

## Limitations

- **WB97X-3C has no Fluorine** (upstream basis gap) — use `def2tzvp` + ECP for F-containing molecules
- All three methods have fixed bases — cannot override with `baselabel` or `ecplabel`
- RI/DF (`j_mode`/`k_mode` = `'RI'`) works normally with all three. R2SCAN-3C's
  def2-mTZVPP basis gets a purpose-built RI-J aux basis (`def2mtzvpprij`,
  Z=1-18+26); B97-3C (def2-mTZVP) and WB97X-3C (vDZP) fall back to the
  generic `def2universaljfit`/`def2universaljkfit` aux basis, same as any
  other basis with no specialized aux-basis case.
