# COSX analytic Hessian: archived alternative paths (not compiled)

These files are the COSX Hessian sources as of commit `ffdd90d1`
(2026-09-28), the last version in which every alternative path below was
still selectable by an environment variable and built with the code.
They are kept for reference only: `legacy/` is not on the Makefile's
VPATH. To run one of these paths, check out `ffdd90d1` and set the
variable listed.

The active sources keep the default path plus the memory opt-ins
`ENGINE_COSX_HESS_LRFIT=1` and `ENGINE_COSX_HESS_PSUB=1` and the
diagnostics `ENGINE_COSX_HESS_VERBOSE`, `..._DUMP`, `..._XCHK`.

Reference case for the numbers: omeprazole HF/6-311G* RIJCOSX Freq,
497 bf, 16 threads, unless noted; "naphthalene" = the distorted
naphthalene element-wise check (run-to-run noise ~1e-10).

| Path (switch at `ffdd90d1`) | What it did | Replaced by / why dropped |
|---|---|---|
| dense tangents (`lowrank` false, `ENGINE_COSX_HESS_T3LR=0`, `..._NOA1=0`) | b^y, a^y, k^y, m^y held as nb x 3N x n per block | occupied-basis factors (B1C, A1C, k^y C); also valid without the overlap fit since `ffdd90d1` |
| full dK/dR (`ENGINE_COSX_KX_FULL`, no occupied coefficients) | n x n x 3N Kx | Kx C (n x nocc x 3N), the only form the CPSCF reads |
| mirror dK/dR (`ENGINE_COSX_KX_MIRROR=1`) | symmetrized M^y + M^yT | single free leg (`70cd360b`): ibuprofen 163 -> 154 s, frequencies <= 0.10 cm-1 |
| two-leg variational dK/dR (`ENGINE_COSX_KX_NONVAR=0`) | second leg x.Gt, kav^y = A a^y + m^y | single leg (`e70e7ea9`): ibuprofen 154 -> 139 s, <= 0.19 cm-1 |
| whole-block k^y / m^y (`ENGINE_COSX_HESS_KCOMP=0`) | nb x 3N x n per thread | per point sub-batch, mapped to k^y C (`4126e612`): peak 12.4 -> 8.9 GB, main 179 -> 172 s |
| stored b^y (`ENGINE_COSX_HESS_BFLY=0`), b^y per sub-batch (`..._BCOMP=1`) | nb x 3N x n per thread / per sub-batch | expanded per point from B1C in the A.R gather (`9ab449e3`): main 170 -> 165 s |
| a^y bound instead of exact maxima (`ENGINE_COSX_HESS_ABOUND=1`, `..._APRUNE=0`) | pair-screen maxima of a^y from a product bound / full a^y | exact pruned maxima (`687e399f`) |
| full-n A (`ENGINE_COSX_HESS_ACMP=0`) | A(n,n,nbs) | A over the pair-active AOs within a 64 MB budget (`71ff3d75`) |
| per-thread x-leg accumulator (`ENGINE_COSX_HESS_XKSHARED=0`, `..._KXX=0`) | XKloc n x nocc x 3N per thread | shared XKsum under column-group locks (`9ab449e3`) |
| fused 34-component integral kernel (`ENGINE_COSX_HESS_FUSED=1`) | one Rys pass for value, ip, ipip, ipvip | slower on caffeine (ESP 756 -> 926 thread-s) |
| Z^y forms (`ENGINE_COSX_HESS_ZT7=0/2`, `..._LRZ=0`) | per-thread n x n x 3N Z^y / T7 per block / shared n x n x 3N Z^y | Z^y as n x nocc x 3N (`42c9e860`) |
| streamed S^y / P^y (`ENGINE_COSX_HESS_PSTREAM=1`, `..._SYMEM`) | S^y, P^y by atom range | fit 9.9 -> 37.7 s; `LRFIT` does the same memory job faster |
| shared S^y accumulator (`ENGINE_COSX_HESS_SNUM=rng`) | one S^y under per-atom locks | +6 s; per-thread atom ranges (`ENGINE_COSX_HESS_SNMEM`) |
| diagnostics (`ENGINE_COSX_HESS_DIAG`, `..._SCRDIAG`, `..._KXDIAG`, `ENGINE_COSX_HESS34_CHECK`, `..._NOKX`) | block/AO fill statistics, pointwise screen audit, Kx leg statistics, fused-kernel check, Kx skipped (timing only) | used while tuning the screens |
