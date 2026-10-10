! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! XC quadrature grid parameters, organized by THEORY LEVEL (functional
! family) x ACCURACY TIER, replacing the ad-hoc ENGINE_XC_USE_PYSCF_TABLE/
! ENGINE_XC_USE_TREUTLER_AIJ/ENGINE_XC_RADPOT_ROWSCALE*/ENGINE_XC_PYSCF_*
! environment-variable knobs used to find the numbers below. Those env
! vars still work (nothing here removes them) - this file is where a
! validated combination gets written down permanently instead of
! re-typed as a shell export.
!
! ONE NAMED ARRAY PER (LEVEL, TIER) CELL, deliberately not a packed
! multi-dimensional array: a reshape's flattening order is easy to get
! wrong (row-major vs column-major, which dimension is fastest) and
! that mistake is silent - it type-checks and compiles, it just reads
! back the wrong cell. A flat, named "XCGRID_<QUANTITY>_<LEVEL>_<TIER>"
! per cell cannot be misindexed: get the name right and the value is
! right by construction.
!
!
! METAGGA/M06's LO and REGULAR cells are NOT independently validated -
! they are seeded 2026-09-12 from the OLD seven-tier ladder
! (mod_xcgrid_defaults.f90/GRID_info's XCGRID_M3_RSCALE/XCGRID_SPH_IN/
! XCGRID_SPH_EDGE) as a documented starting point, not a fresh search:
! REGULAR = the OLD L5 tier (the production grid meta-GGAs actually run
! on today), LO = the OLD COARSE tier (the early-SCF grid xcgrid_dynamic
! already promotes FROM mid-run). See each array's own comment for the
! exact provenance and the caveat that "inherited from the old default"
! is not the same claim as "measured against the new e_exc targets".
! HIGH is still a real sentinel (-1) for METAGGA/M06 (their own HIGH
! cells) - there was no OLD tier obviously corresponding to "even more
! than L7 was ever used for", so it is left for an actual search rather
! than another inherited guess. Every sentinel cell is an explicit -1
! (real) / -1 (integer), not a guess sitting in a data slot pretending to
! be a result. The reading code (mod_data.f90's xcgrid_theory_cell/
! xcgrid_theory_apply) falls back to SMOOTHDFT x REGULAR whenever it hits
! one.
!
!
! HOW TO ADD A NEW LEVEL (beyond R2SCAN/M06): add an XCGRID_LVL_* constant
! in mod_data.f90, add that level's three tiers as new named-array blocks
! down here (start every array at the -1 sentinel), add its branch in
! xcgrid_theory_cell(), and add its ENGINE_XC_THEORY_LEVEL string to the
! parser in xcgrid_tune_from_env() - four places, always in step.
!

! ============================================================
! Per-ROW (1..7 periodic row, H/He=1 .. Fr-and-heavier=7) radial
! shell-count multiplier (xcgrid_radpot_rowscale's per-row generalization).
! ============================================================

real(8), parameter :: XCGRID_ROWSCALE_SMOOTHDFT_REGULAR(7) = &
   (/0.63d0, 0.81d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0/)

real(8), parameter :: XCGRID_ROWSCALE_SMOOTHDFT_LO(7) = &
   (/0.45d0, 0.68d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0/)

real(8), parameter :: XCGRID_ROWSCALE_SMOOTHDFT_HIGH(7) = &
   (/0.756d0, 0.912d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0/)

real(8), parameter :: XCGRID_ROWSCALE_METAGGA_LO(7) = &
   (/0.65d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0/)
real(8), parameter :: XCGRID_ROWSCALE_METAGGA_REGULAR(7) = &
   (/2.10d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0/)
real(8), parameter :: XCGRID_ROWSCALE_METAGGA_HIGH(7) = &
   (/-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0/)

real(8), parameter :: XCGRID_ROWSCALE_R2SCAN_LO(7) = &
   (/0.4d0, 1.15d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0/)
real(8), parameter :: XCGRID_ROWSCALE_R2SCAN_REGULAR(7) = &
   (/0.75d0, 1.45d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0/)
real(8), parameter :: XCGRID_ROWSCALE_R2SCAN_HIGH(7) = &
   (/0.8d0, 2.45d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0/)

real(8), parameter :: XCGRID_ROWSCALE_R2SCAN3C_REGULAR(7) = &
   (/0.8d0, 2.45d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0/)
integer, parameter :: XCGRID_ZONE434_R2SCAN3C_REGULAR(5) = (/14, 38, 110, 350, 146/)
integer, parameter :: XCGRID_ZONE302_R2SCAN3C_REGULAR(5) = (/14, 38, 110, 194, 146/)
integer, parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN3C_REGULAR_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN3C_REGULAR_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN3C_REGULAR_SCALE(4) = &
   (/2.45d0, 1.65d0, 1.75d0, 1.6d0/)

real(8), parameter :: XCGRID_ROWSCALE_M06_LO(7) = &
   (/0.5d0, 0.9d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0/)
real(8), parameter :: XCGRID_ROWSCALE_M06_REGULAR(7) = &
   (/0.9d0, 1.3d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0/)
real(8), parameter :: XCGRID_ROWSCALE_M06_HIGH(7) = &
   (/1.1d0, 1.3d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0, 2.10d0/)

real(8), parameter :: XCGRID_ROWSCALE_WB97XD_LO(7) = &
   (/0.45d0, 0.90d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0/)
real(8), parameter :: XCGRID_ROWSCALE_WB97XD_REGULAR(7) = &
   (/0.55d0, 1.10d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0/)
real(8), parameter :: XCGRID_ROWSCALE_WB97XD_HIGH(7) = &
   (/0.80d0, 1.55d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0/)

real(8), parameter :: XCGRID_ROWSCALE_WB97MV_LO(7) = &
   (/0.40d0, 0.60d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0, 0.65d0/)
real(8), parameter :: XCGRID_ROWSCALE_WB97MV_REGULAR(7) = &
   (/0.60d0, 0.90d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0/)
real(8), parameter :: XCGRID_ROWSCALE_WB97MV_HIGH(7) = &
   (/0.70d0, 1.80d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0, 0.95d0/)

! ============================================================
! Per-ELEMENT rowscale overrides (sparse - takes priority over the
! per-row table above when an element has one). Indexed by atomic
! number Z, 1..103 like data/parameter.h's own covrad table.
! ============================================================

integer, parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_REGULAR_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_REGULAR_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_REGULAR_SCALE(4) = &
   (/0.76d0, 0.85d0, 0.78d0, 0.74d0/)

! SMOOTHDFT / LO - VALIDATED. C=0.68, N=0.55, O=0.70, F=0.68. Same
! pattern as REGULAR: N is the accuracy bottleneck for N-rich cases
! (hydrazine/urea), F/O were raised relative to REGULAR's cut to buy
! back the angular truncation on the F/O-rich molecules
! (formyl_fluoride/difluoromethane/water) that bind the 1e-5 target.
integer, parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_LO_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_LO_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_LO_SCALE(4) = &
   (/0.68d0, 0.55d0, 0.70d0, 0.68d0/)

! SMOOTHDFT / HIGH - VALIDATED. C=0.912 (the row2+ base), N=0.98,
! O=0.936, F=0.888. Unlike REGULAR/LO, the extra accuracy comes from the
! angular peak (below), not from raising the radial rowscales - the
! radial shell count is a rounded integer (radpot = nint(NRAD*rowscale)),
! so the error is a step function of rowscale and 0.912 is the shell
! count that clears the 3e-7 bar for the C-rich worst cases. F/O are the
! worst-case point users (formyl_fluoride/difluoromethane at ~1.72x
! REGULAR), so they are not raised further.
integer, parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_HIGH_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_HIGH_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_SMOOTHDFT_HIGH_SCALE(4) = &
   (/0.912d0, 0.98d0, 0.936d0, 0.888d0/)

integer, parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_LO_N = 0
integer, parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_LO_Z(1) = (/0/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_LO_SCALE(1) = (/-1.0d0/)

integer, parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_REGULAR_N = 0
integer, parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_REGULAR_Z(1) = (/0/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_REGULAR_SCALE(1) = (/-1.0d0/)

integer, parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_HIGH_N = 0
integer, parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_HIGH_Z(1) = (/0/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_METAGGA_HIGH_SCALE(1) = (/-1.0d0/)

integer, parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_LO_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_LO_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_LO_SCALE(4) = &
   (/1.15d0, 0.85d0, 0.95d0, 1.05d0/)

! R2SCAN / REGULAR - VALIDATED 2026-09-13. C=1.45, N=1.4, O=1.35, F=1.45.
! N is highest (acetonitrile binds the 3e-6 target); see the ROWSCALE
! comment for the accuracy/point summary.
integer, parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_REGULAR_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_REGULAR_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_REGULAR_SCALE(4) = &
   (/1.45d0, 1.4d0, 1.35d0, 1.45d0/)

integer, parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_HIGH_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_HIGH_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_R2SCAN_HIGH_SCALE(4) = &
   (/2.45d0, 1.65d0, 1.75d0, 1.6d0/)

! M06 / LO - VALIDATED 2026-09-30 on M06-2X: C=0.9, N=0.7, O=1.0, F=0.8
! (H 0.5). Worst e_exc 1.63e-5 (octane), target 2e-5.
integer, parameter :: XCGRID_ELEM_OVERRIDE_M06_LO_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_M06_LO_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_M06_LO_SCALE(4) = &
   (/0.9d0, 0.7d0, 1.0d0, 0.8d0/)

! M06 / REGULAR - VALIDATED 2026-09-30 on M06-2X: C=1.3, N=1.0, O=1.0,
! F=1.1 (H 0.9). Worst e_exc 1.94e-6 (urea), target 2e-6; about half the
! points of the old L5 production grid (which had 1.65e-6 worst).
integer, parameter :: XCGRID_ELEM_OVERRIDE_M06_REGULAR_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_M06_REGULAR_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_M06_REGULAR_SCALE(4) = &
   (/1.3d0, 1.0d0, 1.0d0, 1.1d0/)

! M06 / HIGH - VALIDATED 2026-09-30 on M06-2X: C=1.3, N=1.2, O=1.3,
! F=1.3 (H 1.1). Worst e_exc 5.47e-7 (formyl_fluoride), target 6e-7.
integer, parameter :: XCGRID_ELEM_OVERRIDE_M06_HIGH_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_M06_HIGH_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_M06_HIGH_SCALE(4) = &
   (/1.3d0, 1.2d0, 1.3d0, 1.3d0/)

! WB97X-D / LO - VALIDATED 2026-09-16 (transcribed 2026-09-18). C=0.90,
! N=0.80, O=0.80, F=0.80 (H is row1=0.45).
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_LO_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_LO_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_LO_SCALE(4) = &
   (/0.90d0, 0.80d0, 0.80d0, 0.80d0/)

! WB97X-D / REGULAR - VALIDATED 2026-09-16 (transcribed 2026-09-18).
! C=1.10, N=0.80, O=0.90, F=1.05.
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_REGULAR_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_REGULAR_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_REGULAR_SCALE(4) = &
   (/1.10d0, 0.80d0, 0.90d0, 1.05d0/)

! WB97X-D / HIGH - VALIDATED 2026-09-16 (transcribed 2026-09-18). C=1.55,
! N=1.25, O=1.95, F=1.20. See this file's header block for the 5/33
! noise-floor caveat on this tier (formyl_fluoride/glyoxal/benzene/
! formaldehyde/formic_acid pass without a real sub-5e-7 measurement).
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_HIGH_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_HIGH_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_WB97XD_HIGH_SCALE(4) = &
   (/1.55d0, 1.25d0, 1.95d0, 1.20d0/)

integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_LO_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_LO_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_LO_SCALE(4) = &
   (/0.60d0, 0.50d0, 0.60d0, 0.70d0/)

integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_REGULAR_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_REGULAR_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_REGULAR_SCALE(4) = &
   (/0.90d0, 0.70d0, 0.90d0, 0.90d0/)

integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_HIGH_N = 4
integer, parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_HIGH_Z(4) = (/6, 7, 8, 9/)
real(8), parameter :: XCGRID_ELEM_OVERRIDE_WB97MV_HIGH_SCALE(4) = &
   (/1.80d0, 0.90d0, 1.00d0, 1.10d0/)

! ============================================================
! Angular zone orders (xcgrid_prune_sphpot_pyscf's pat434/pat302,
! zones 1..5: near-core, 2, pre-peak flank, peak, post-peak flank).
! ============================================================

integer, parameter :: XCGRID_ZONE434_SMOOTHDFT_REGULAR(5) = (/14, 38, 110, 350, 146/)
integer, parameter :: XCGRID_ZONE302_SMOOTHDFT_REGULAR(5) = (/14, 38, 110, 170, 146/)

! SMOOTHDFT / LO - VALIDATED. zone434 = (14,26,86,230,86), H/He peak
! (zone302(4)) = 110. The heavy-atom peak drops 350->230 and both flanks
! drop to 86/86; near-core stays 14/26. This is where the LO point saving
! comes from (radial can't drop much - see the rowscale comment), and it
! is the frontier at e_exc <= 1e-5: lowering any single entry one more
! Lebedev step pushes the worst case over 1e-5. The other zone302
! entries are unused (shared with pat434) but kept equal for coherence.
integer, parameter :: XCGRID_ZONE434_SMOOTHDFT_LO(5) = (/14, 26, 86, 230, 86/)
integer, parameter :: XCGRID_ZONE302_SMOOTHDFT_LO(5) = (/14, 26, 86, 110, 86/)

! SMOOTHDFT / HIGH - VALIDATED. zone434 = (14,38,110,590,170), H/He peak
! (zone302(4)) = 194. zone1/2/3 stay at REGULAR's 14/38/110 unchanged;
! the heavy-atom peak (zone4) moves 350->590 - benzene (the worst case
! at 2.74e-7) needs the peak order, raising the flanks or zone1 only
! adds points without lowering the worst case. zone5 (post-peak flank)
! also moves, 146->170, one Lebedev step up from REGULAR - NOT "stays
! put" despite sharing REGULAR's zone1/2/3. zone302's other entries are
! unused (shared with pat434) but kept equal for coherence.
integer, parameter :: XCGRID_ZONE434_SMOOTHDFT_HIGH(5) = (/14, 38, 110, 590, 170/)
integer, parameter :: XCGRID_ZONE302_SMOOTHDFT_HIGH(5) = (/14, 38, 110, 194, 170/)

! METAGGA/R2SCAN/M06 LO/REGULAR - seeded 2026-09-12 from the OLD seven-
! tier ladder, same inheritance and caveat as the ROWSCALE arrays above.
! LO = OLD COARSE tier: XCGRID_SPH_IN(2)=XCGRID_SPH_EDGE(2)=302, i.e. no
! inner/outer split at all, so all 5 zones get 302. REGULAR = OLD L5
! tier: XCGRID_SPH_IN(5)=974 (inner), XCGRID_SPH_EDGE(5)=590 (outer) -
! the OLD ladder only splits radial shells into "inner"/"outer" by index
! fraction, not this table's 5 physical r/r_atom zones, so the mapping
! used here is the simplest faithful one: zones 1-4 (near-core through
! the density peak) get the inner value, zone 5 (the outermost,
! post-peak flank) gets the outer value - not a claim that this split
! point is correct for the new zone scheme, just the least invented
! choice available without a real search.
integer, parameter :: XCGRID_ZONE434_METAGGA_LO(5) = (/302, 302, 302, 302, 302/)
integer, parameter :: XCGRID_ZONE302_METAGGA_LO(5) = (/302, 302, 302, 302, 302/)
integer, parameter :: XCGRID_ZONE434_METAGGA_REGULAR(5) = (/974, 974, 974, 974, 590/)
integer, parameter :: XCGRID_ZONE302_METAGGA_REGULAR(5) = (/974, 974, 974, 974, 590/)
integer, parameter :: XCGRID_ZONE434_METAGGA_HIGH(5) = (/-1,-1,-1,-1,-1/)
integer, parameter :: XCGRID_ZONE302_METAGGA_HIGH(5) = (/-1,-1,-1,-1,-1/)

! R2SCAN / LO - VALIDATED 2026-09-13. zone434 = (14,38,110,266,110), H/He
! peak (zone302(4)) = 194. Near-core stays cheap (14/38), the heavy-atom
! peak is 266 and both flanks 110. The angular is what buys the LO point
! cut here; radial rowscales are already near their own floor.
integer, parameter :: XCGRID_ZONE434_R2SCAN_LO(5) = (/14, 38, 110, 266, 110/)
integer, parameter :: XCGRID_ZONE302_R2SCAN_LO(5) = (/14, 38, 110, 194, 110/)

! R2SCAN / REGULAR - VALIDATED 2026-09-13. zone434 = (14,38,110,350,146),
! H/He peak = 194. Same near-core/flanks as SMOOTHDFT/REGULAR; the heavy-
! atom peak is 350. R2SCAN's radial dominance means this angular scheme
! is cheap relative to the rowscales above.
integer, parameter :: XCGRID_ZONE434_R2SCAN_REGULAR(5) = (/14, 38, 110, 350, 146/)
integer, parameter :: XCGRID_ZONE302_R2SCAN_REGULAR(5) = (/14, 38, 110, 194, 146/)

! R2SCAN / HIGH - VALIDATED 2026-09-13. zone434 = (50,146,146,434,146),
! H/He peak (zone302(4)) = 434. Unlike SMOOTHDFT/HIGH (where only the
! heavy-atom peak moved), R2SCAN's 5e-7 target forces the WHOLE zone
! scheme up - every zone is at 146/434, because the angular error (not
! radial) is the binding term for tau-dependent functionals. This is
! where HIGH's ~2x-REGULAR point cost comes from.
integer, parameter :: XCGRID_ZONE434_R2SCAN_HIGH(5) = (/50, 146, 146, 434, 146/)
integer, parameter :: XCGRID_ZONE302_R2SCAN_HIGH(5) = (/50, 146, 146, 434, 146/)

! M06 angular zones - VALIDATED 2026-09-30 on M06-2X (heavy atom peak /
! H-He peak): LO (26,74,194,350,194)/266, REGULAR (110,302,434,590,434)/
! 590, HIGH (110,302,434,770,434)/770.
integer, parameter :: XCGRID_ZONE434_M06_LO(5) = (/26, 74, 194, 350, 194/)
integer, parameter :: XCGRID_ZONE302_M06_LO(5) = (/26, 74, 194, 266, 194/)
integer, parameter :: XCGRID_ZONE434_M06_REGULAR(5) = (/110, 302, 434, 590, 434/)
integer, parameter :: XCGRID_ZONE302_M06_REGULAR(5) = (/110, 302, 434, 590, 434/)
integer, parameter :: XCGRID_ZONE434_M06_HIGH(5) = (/110, 302, 434, 770, 434/)
integer, parameter :: XCGRID_ZONE302_M06_HIGH(5) = (/110, 302, 434, 770, 434/)

! WB97X-D / LO - VALIDATED 2026-09-16 (transcribed 2026-09-18). zone434 =
! (14,38,194,266,86), H/He peak (zone302(4)) = 110.
integer, parameter :: XCGRID_ZONE434_WB97XD_LO(5) = (/14, 38, 194, 266, 86/)
integer, parameter :: XCGRID_ZONE302_WB97XD_LO(5) = (/14, 38, 194, 110, 86/)

! WB97X-D / REGULAR - VALIDATED 2026-09-16 (transcribed 2026-09-18).
! zone434 = (14,110,302,434,110), H/He peak = 194.
integer, parameter :: XCGRID_ZONE434_WB97XD_REGULAR(5) = (/14, 110, 302, 434, 110/)
integer, parameter :: XCGRID_ZONE302_WB97XD_REGULAR(5) = (/14, 110, 302, 194, 110/)

! WB97X-D / HIGH - VALIDATED 2026-09-16 (transcribed 2026-09-18). zone434
! = (50,194,434,770,194), H/He peak = 302.
integer, parameter :: XCGRID_ZONE434_WB97XD_HIGH(5) = (/50, 194, 434, 770, 194/)
integer, parameter :: XCGRID_ZONE302_WB97XD_HIGH(5) = (/50, 194, 434, 302, 194/)

! WB97M-V / LO and REGULAR - VALIDATED 2026-09-18. Both use zone434 =
! (14,38,110,350,146), H/He peak (zone302(4)) = 194. Radial-dominated
! (like R2SCAN): raising the NEARCORE orders measurably raised cost for
! little accuracy, so LO and REGULAR share the R2SCAN/REGULAR angular.
integer, parameter :: XCGRID_ZONE434_WB97MV_LO(5) = (/14, 38, 110, 350, 146/)
integer, parameter :: XCGRID_ZONE302_WB97MV_LO(5) = (/14, 38, 110, 194, 146/)
integer, parameter :: XCGRID_ZONE434_WB97MV_REGULAR(5) = (/14, 38, 110, 350, 146/)
integer, parameter :: XCGRID_ZONE302_WB97MV_REGULAR(5) = (/14, 38, 110, 194, 146/)

! WB97M-V / HIGH - VALIDATED 2026-09-18. zone434 = (50,146,146,434,146),
! H/He peak = 434. The 5e-7 target forces the whole zone scheme up, same
! as R2SCAN/HIGH.
integer, parameter :: XCGRID_ZONE434_WB97MV_HIGH(5) = (/50, 146, 146, 434, 146/)
integer, parameter :: XCGRID_ZONE302_WB97MV_HIGH(5) = (/50, 146, 146, 434, 146/)

! ============================================================
! Becke-partition atomic-size adjustment scheme: 1=Becke (JCP 88, 2547
! (1988), chi=r_i/r_j), 2=Treutler (JCP 102, 346 (1995),
! chi=sqrt(r_i)/sqrt(r_j)) - see xcgrid_use_treutler_aij's own comment
! in mod_data.f90 for the derivation. -1 = not populated.
! ============================================================

integer, parameter :: XCGRID_AIJ_SCHEME_SMOOTHDFT_REGULAR = 2   ! Treutler - VALIDATED
integer, parameter :: XCGRID_AIJ_SCHEME_SMOOTHDFT_LO = 2   ! Treutler - VALIDATED
integer, parameter :: XCGRID_AIJ_SCHEME_SMOOTHDFT_HIGH = 2   ! Treutler - VALIDATED
! METAGGA/R2SCAN/M06 LO/REGULAR: the OLD seven-tier ladder never had an
! aij scheme concept (it always used classic 1988 Becke), so there is
! nothing to "inherit" here - set to Treutler (2) anyway, carrying
! forward SMOOTHDFT's validated choice rather than reverting to classic
! Becke, which has never been shown to be better in this codebase.
integer, parameter :: XCGRID_AIJ_SCHEME_METAGGA_LO = 2
integer, parameter :: XCGRID_AIJ_SCHEME_METAGGA_REGULAR = 2
integer, parameter :: XCGRID_AIJ_SCHEME_METAGGA_HIGH = -1
integer, parameter :: XCGRID_AIJ_SCHEME_R2SCAN_LO = 2
integer, parameter :: XCGRID_AIJ_SCHEME_R2SCAN_REGULAR = 2
integer, parameter :: XCGRID_AIJ_SCHEME_R2SCAN_HIGH = 2
integer, parameter :: XCGRID_AIJ_SCHEME_M06_LO = 2
integer, parameter :: XCGRID_AIJ_SCHEME_M06_REGULAR = 2
integer, parameter :: XCGRID_AIJ_SCHEME_M06_HIGH = 2
! WB97X-D: not separately re-derived, carrying forward SMOOTHDFT/R2SCAN/
! M06's validated Treutler choice (same reasoning as those cells' own
! comments - never shown to be worse in this codebase).
integer, parameter :: XCGRID_AIJ_SCHEME_WB97XD_LO = 2
integer, parameter :: XCGRID_AIJ_SCHEME_WB97XD_REGULAR = 2
integer, parameter :: XCGRID_AIJ_SCHEME_WB97XD_HIGH = 2
! WB97M-V: same as WB97X-D, carrying forward Treutler (2).
integer, parameter :: XCGRID_AIJ_SCHEME_WB97MV_LO = 2
integer, parameter :: XCGRID_AIJ_SCHEME_WB97MV_REGULAR = 2
integer, parameter :: XCGRID_AIJ_SCHEME_WB97MV_HIGH = 2
