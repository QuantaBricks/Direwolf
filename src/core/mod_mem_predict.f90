! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_mem_predict: predicts this run's TOTAL peak resident memory (not

module mod_mem_predict
implicit none
private
public :: predict_energy_peak_bytes, predict_force_extra_bytes, predict_core_overhead_bytes, &
          TSTACK_BLOCK_BYTES, TSTACK_MIN_BLOCK_WIDTH, GRID_RAGGED_FRACTION_EST, MEM_SAFETY_MARGIN

real(8),parameter :: TSTACK_BLOCK_BYTES = 5.0d7
integer,parameter :: TSTACK_MIN_BLOCK_WIDTH = 100

real(8),parameter :: DF_COMPACT_FRACTION_EST = 0.8d0

real(8),parameter :: DF_RSH_EXTRA_FRACTION = 0.9d0

real(8),parameter :: GRID_RAGGED_FRACTION_EST = 0.6d0

real(8),parameter :: MEM_SAFETY_MARGIN = 1.2d0

contains

integer(8) function predict_core_overhead_bytes(nConts, natoms) result(bytes)
implicit none
integer,intent(in) :: nConts, natoms
integer(8) :: core_bytes, diis_bytes, transient_bytes, n8
integer :: diis_max
character(len=8) :: diis_env

n8 = int(nConts,8)
core_bytes = 9_8*n8*n8*8_8 + int(natoms,8)*int(natoms,8)*4_8

diis_max = 10
diis_env = ""
call get_environment_variable("ENGINE_DIIS_WINDOW", diis_env)
if (len_trim(diis_env) .gt. 0) read(diis_env,*) diis_max
if (diis_max .lt. 1) diis_max = 1
diis_bytes = 6_8*n8*n8*int(diis_max,8)*4_8

transient_bytes = 4_8*n8*n8*8_8
bytes = core_bytes + diis_bytes + transient_bytes
end function predict_core_overhead_bytes

integer(8) function predict_energy_peak_bytes(nConts, natoms, nContsAux, need_k, &
                                               want_store, grid_cache_bytes, rs_omega_active, &
                                               df_total_pairs_actual) result(bytes)
implicit none
integer,intent(in) :: nConts, natoms, nContsAux
integer(8),intent(in) :: grid_cache_bytes
logical,intent(in) :: need_k, want_store
logical,intent(in),optional :: rs_omega_active
integer(8),intent(in),optional :: df_total_pairs_actual
integer(8) :: core_bytes, diis_bytes, transient_bytes, store_bytes
integer(8) :: n8, nAux8, dfBK_slot_bytes
integer :: diis_max
character(len=8) :: diis_env

n8 = int(nConts,8)
nAux8 = int(nContsAux,8)
core_bytes = 9_8*n8*n8*8_8 + int(natoms,8)*int(natoms,8)*4_8

diis_max = 10
diis_env = ""
call get_environment_variable("ENGINE_DIIS_WINDOW", diis_env)
if (len_trim(diis_env) .gt. 0) read(diis_env,*) diis_max
if (diis_max .lt. 1) diis_max = 1
diis_bytes = 6_8*n8*n8*int(diis_max,8)*4_8

transient_bytes = 4_8*n8*n8*8_8

store_bytes = nAux8*nAux8*8_8 + nAux8*8_8 &
            + nAux8*nAux8*8_8 &
            + 2_8*nAux8*nAux8*8_8
if (want_store) then
   if (present(df_total_pairs_actual)) then
      dfBK_slot_bytes = df_total_pairs_actual*nAux8*8_8
      store_bytes = store_bytes + dfBK_slot_bytes
      if (present(rs_omega_active)) then
         if (rs_omega_active) &
            store_bytes = store_bytes + 2_8*dfBK_slot_bytes
      endif
   else
      dfBK_slot_bytes = n8*n8*nAux8*8_8
      store_bytes = store_bytes &
                  + int(real(dfBK_slot_bytes,8)*DF_COMPACT_FRACTION_EST,8)
      if (present(rs_omega_active)) then
         if (rs_omega_active) &
            store_bytes = store_bytes &
                        + int(real(dfBK_slot_bytes,8)*DF_COMPACT_FRACTION_EST*DF_RSH_EXTRA_FRACTION,8)
      endif
   endif
endif

bytes = int(real(core_bytes + diis_bytes + transient_bytes + grid_cache_bytes + store_bytes,8) &
            * MEM_SAFETY_MARGIN, 8)
end function predict_energy_peak_bytes

integer(8) function predict_force_extra_bytes(nConts, nContsAux, need_k, want_store) result(bytes)
implicit none
integer,intent(in) :: nConts, nContsAux
logical,intent(in) :: need_k, want_store
integer(8) :: force_local_bytes, tstack_bytes, n8, nAux8
integer :: nT
n8 = int(nConts,8)
nAux8 = int(nContsAux,8)
force_local_bytes = 6_8*n8*n8*8_8
tstack_bytes = 0_8
if (want_store) then
   nT = 1
   if (need_k) nT = 3
   block
      integer(8) :: row_cost_bytes, block_width
      row_cost_bytes = n8*n8*int(nT,8)*8_8
      block_width = max(int(TSTACK_BLOCK_BYTES/real(row_cost_bytes,8),8), int(TSTACK_MIN_BLOCK_WIDTH,8))
      block_width = min(block_width, nAux8)
      tstack_bytes = n8*n8*block_width*int(nT,8)*8_8 + nAux8*nAux8*int(nT,8)*8_8
   end block
endif
bytes = int(real(force_local_bytes + tstack_bytes,8) * MEM_SAFETY_MARGIN, 8)
end function predict_force_extra_bytes

end module mod_mem_predict
