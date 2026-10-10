# Intrinsic reaction coordinate (IRC) examples

`runtype = 'irc'` follows the reaction path from a **transition state** down both
sides, to the reactant and the product, so a TS can be checked to connect the
intended minima. The input geometry must be a TS (exactly one imaginary mode); get it
with `opt_ts = .true.` first (see `examples/optimization/`). `irc` stands alone: it is not
combined with `opt`, `freq` or `polar`. Run from the repository root, e.g.
`./Direwolf examples/irc/nh3_irc.inp`.

Method: mass-weighted Cartesian coordinates (arc length in amu^1/2 bohr), one local
quadratic (LQA) step per point, one energy + gradient per point. The analytic Hessian is
computed at the TS, then updated with the Bofill formula; it is recomputed every
`irc_hess_every` points and whenever a step went uphill (Hratchian and Schlegel,
J. Chem. Theory Comput. 1, 61 (2005)). With `irc_corrector` (default on) the predictor point is followed by the
paper's corrector: a distance-weighted interpolant is built from the two second-order Taylor expansions at both ends of the
step (energy, gradient and Hessian already at hand), steepest descent on it is integrated over the same arc length with Euler steps and
Bulirsch-Stoer extrapolation, and the next step starts from that corrected point with the fitted gradient: no extra
energy or gradient calculation. It is skipped for a point where the extrapolation did not converge (error above 1e-3) or the
shift is more than a quarter of the step (near a minimum, where the fitted surface is flat). The first step from the TS has no corrector,
and the energy printed for a corrected point is the one on the fitted surface. HCN -> HNC at step 0.3, against a step-0.05 path:
mean distance of the path points to the reference path 0.0110 -> 0.0052 A, same number of gradients.

The optional `&irc` namelist (or `[irc]` TOML table):

| setting | default | meaning |
|---|---|---|
| `irc_direction` | `both` | `both`, `forward` (along +v, v the imaginary mode) or `backward` |
| `irc_step` | 0.10 | arc length per point, amu^1/2 bohr |
| `irc_maxpoints` | 50 | points per direction |
| `irc_hess_every` | 5 | analytic Hessian every N points; 0 = only at the TS |
| `irc_gmax` | 3e-4 | stop when the largest Cartesian gradient component is below this (Hartree/bohr) |
| `irc_corrector` | `.true.` | the corrector step described above |
| `irc_dbo` | 0.01 | stop when no Mayer bond order changed by more than this over the last 3 points, once the path is 1 kcal/mol below the TS and the gradient is below 2e-3; 0 = off |

A direction also stops when the quadratic model reaches its minimum, at `irc_maxpoints`, or on a failure
(the `[IRCRESULT]` block says why). `irc_dbo` ends a path early, at a point where the molecule
is recognizable but not yet at the minimum: optimize the end point (`runtype = 'opt'`) to
get the minimum itself.

Output: `<base>.irc.xyz` (backward end ... TS ... forward end) and `<base>.irc.dat`
(s, E, E - E(TS) in kcal/mol, gmax). `plot_irc.py <base>` draws the energy profile and the
distances that change most:

```
python3 examples/irc/plot_irc.py nh3_irc
```

| file | shows |
|---|---|
| `nh3_irc.inp` | NH3 inversion, 13 points a side, both ends the same minimum |
| `hcn_hnc_irc.inp` | HCN <-> HNC, an asymmetric path; barriers 47.3 / 33.9 kcal/mol |
| `h3_irc.inp` | H + H2 exchange, open shell (UKS, doublet), 17 points a side, barrier 4.25 kcal/mol |
| `pd_h2_irc.inp` | H2 + Pd(PH3)2, Pd with the def2 ECP: TS (one imaginary mode, -452 cm-1) to the eta2-H2 complex on one side (6.4 kcal/mol down) and to the cis-PdH2(PH3)2 minimum on the other (only 0.1 kcal/mol down - a very shallow well) |

The Pd TS did not come from a guess: P-RFO on the first guesses ended in the eta2-H2 minimum or did not converge. A relaxed scan found it:
`opt_constraints = 'bond 4 5 1.10'` (the H-H distance) at 0.90, 1.00 ... 1.60 A, each point started from the
previous one; the energy flattens at 1.45-1.60 A, and `opt_ts` from the 1.60 A point converged to the TS in a few cycles.
| `sn2_cosmo_irc.inp` | Cl- + CH3Cl in water (`&cosmo`, COSMO); 31 points a side, barrier 16.8 kcal/mol, ends equal to 1.7e-4 Ha |
| `sn2_iodide_irc.inp` | I- + CH3I, iodine with the def2 ECP; 27 points a side, barrier 5.26 kcal/mol, ends equal to 1e-7 Ha |

Finding the TS is the hard part, not the IRC: P-RFO from a rough guess often ends in a minimum or a flat
region. `opt_ts` therefore checks the structure it converged to with one analytic Hessian and prints either
`[OPT] Hessian check: one imaginary mode - a transition state` or a WARNING with the number of imaginary modes,
and `irc` refuses a structure that does not have exactly one. Tried and not usable as IRC tests, because
they have no separate TS on the electronic energy surface or a floppy ion-molecule surface: I + H2 -> HI + H (the TS is
the products), I + HI (a bound I-H-I intermediate), I + CH4 (the converged structure had 3 imaginary modes),
Cl- + CH3I from a guess (the search ran off to the products).
| `pd_pme3_h2_irc.inp` | H2 + Pd(PMe3)2, 29 atoms with the Pd ECP: TS (-719 cm-1) to the eta2-H2 complex (13 points, 3.1 kcal/mol down) and the cis-PdH2(PMe3)2 side (8 points, 2.0 kcal/mol down) |

| `snh3_c2h4_irc.inp` | SnH3 radical + C2H4, open shell (doublet) with the Sn ECP: TS (-176 cm-1), 27 + 30 points, product side 8.09 kcal/mol down (gradient 4.9e-4, flat), reactant side stopped at 1.31 kcal/mol down |

Open-shell reactions with an ECP atom took several tries to find a usable TS: with iodine radicals the 5p hole
is degenerate in a linear geometry and the unrestricted CPHF response blows up (Hessian eigenvalues of -6e9: I + CH3I, I + HI),
and I + H2, I + CH4, H + CH3I, H + SnH4, SnH3 + HSnH3 are (nearly) barrierless on this surface. HI and HI+ frequencies are fine
(2341 and 2243 cm-1), so the open-shell ECP Hessian itself works. A constrained scan in `opt_coord = 'ric'` of a floppy
ion-molecule complex (Cl- ... CH3I) can stall on COSX noise; `opt_coord = 'cart'` converged in 15 cycles.

