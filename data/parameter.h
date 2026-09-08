! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later


! NOTE: every literal in this file carries an explicit d0 suffix. Without
! it a bare decimal like 1.88972612456506 is a DEFAULT-KIND (single
! precision) literal that gets truncated to ~7 significant digits before
! being promoted into the real(8) parameter - measured 1.98e-8 relative
! error on ans2bohr and 2.78e-8 on PI, which is what was left of the
! nuclear-repulsion mismatch against pyscf after the constant itself was
! corrected. PI feeds every grid weight (radw ~ 2*PI/(nrad+1)*...*4*PI),
! so the error propagated into the XC quadrature too.
!
! CODATA 2018 Bohr radius, a0 = 0.52917721092 Angstrom - the value pyscf
! uses (pyscf.data.nist.BOHR), so cross-code geometries and hence nuclear
! repulsion energies agree. The previous 1.889725989 was the CODATA 1986
! value (a0 = 0.5291772489), off by 7.2e-8 relative: enough to shift the
! nuclear repulsion of a 42-atom molecule by ~9e-5 Ha and to account for
! nearly all of a 1e-4 Ha total-energy gap against ORCA/pyscf.
real(8),parameter   :: ans2bohr = 1.88972612456506d0
real(8),parameter   :: PI = 3.14159265358979324d0

! covrad: per-element radius parameter (Angstrom) feeding grid_gen.f90's
! Becke-partition atomic-size adjustment (chi=covrad(i)/covrad(j)) and
! the radial-grid scale "parm" (=covrad(Z)/2*ans2bohr). Originally sized
! 26 (H-Fe) - any element past Fe (Z>26) indexed PAST THE END of this
! compile-time constant array: undefined-behavior read, empirically
! landing on a near-zero value that silently zeroed radw and therefore
! the ENTIRE DFT XC grid contribution for that atom (Exc exactly
! 0.000000000 every SCF cycle, no crash, no warning) - found via a bare
! Pd2+/lanl2dz/B3LYP atom test whose energy was ~10 Ha off from Psi4,
! traced by bisecting Fe (Z=26, in-bounds, correct) vs Pd (Z=46,
! out-of-bounds, broken) with an otherwise-identical ECP/basis
! structure. This is why every prior ECP validation (Fe/Pd/Pt) used
! HF, never DFT - HF's Fock build never touches the XC grid at all, so
! it never exercised this array. Extended Z=27-103 with standard
! (Cordero et al. 2008-style) single-bond covalent radii in Angstrom,
! same convention/scale as the original 26 entries; H-Fe values
! (indices 1-26) are UNCHANGED from before this fix.
real(8),parameter   :: covrad(103) = [0.38d0,0.32d0,1.34d0,0.9d0,0.82d0,0.77d0,0.75d0,&
                                     0.73d0,0.71d0,0.69d0,1.54d0,1.30d0,1.18d0,1.11d0,&
                                     1.06d0,1.02d0,0.99d0,0.97d0,1.96d0,1.74d0,1.44d0,&
                                     1.36d0,1.25d0,1.27d0,1.39d0,1.25d0,&
                                     1.11d0,1.10d0,1.12d0,1.18d0,1.24d0,1.21d0,1.21d0,&
                                     1.16d0,1.14d0,1.17d0,2.10d0,1.85d0,1.63d0,1.54d0,&
                                     1.47d0,1.38d0,1.28d0,1.25d0,1.25d0,1.20d0,1.28d0,&
                                     1.36d0,1.42d0,1.40d0,1.40d0,1.36d0,1.33d0,1.31d0,&
                                     2.32d0,1.96d0,1.80d0,1.63d0,1.76d0,1.74d0,1.73d0,&
                                     1.72d0,1.68d0,1.69d0,1.68d0,1.67d0,1.66d0,1.65d0,&
                                     1.64d0,1.70d0,1.62d0,1.52d0,1.46d0,1.37d0,1.31d0,&
                                     1.29d0,1.22d0,1.23d0,1.24d0,1.33d0,1.44d0,1.44d0,&
                                     1.51d0,1.45d0,1.47d0,1.42d0,2.23d0,2.01d0,1.86d0,&
                                     1.75d0,1.69d0,1.70d0,1.71d0,1.72d0,1.66d0,1.66d0,&
                                     1.68d0,1.68d0,1.65d0,1.67d0,1.73d0,1.76d0,1.61d0 ]

