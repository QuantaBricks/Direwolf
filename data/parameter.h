! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0


real(8),parameter   :: ans2bohr = 1.889725989
real(8),parameter   :: PI = 3.14159265359

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
real(8),parameter   :: covrad(103) = [0.38,0.32,1.34,0.9,0.82,0.77,0.75,&
                                     0.73,0.71,0.69,1.54,1.30,1.18,1.11,&
                                     1.06,1.02,0.99,0.97,1.96,1.74,1.44,&
                                     1.36,1.25,1.27,1.39,1.25,&
                                     1.11,1.10,1.12,1.18,1.24,1.21,1.21,&
                                     1.16,1.14,1.17,2.10,1.85,1.63,1.54,&
                                     1.47,1.38,1.28,1.25,1.25,1.20,1.28,&
                                     1.36,1.42,1.40,1.40,1.36,1.33,1.31,&
                                     2.32,1.96,1.80,1.63,1.76,1.74,1.73,&
                                     1.72,1.68,1.69,1.68,1.67,1.66,1.65,&
                                     1.64,1.70,1.62,1.52,1.46,1.37,1.31,&
                                     1.29,1.22,1.23,1.24,1.33,1.44,1.44,&
                                     1.51,1.45,1.47,1.42,2.23,2.01,1.86,&
                                     1.75,1.69,1.70,1.71,1.72,1.66,1.66,&
                                     1.68,1.68,1.65,1.67,1.73,1.76,1.61 ]

