! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Density-fitting Coulomb+exchange-force gradient, disk/memory-stored 3-center integral mode.

submodule (mod_integrals) df_force_store_impl
implicit none
contains

module subroutine integrals_force_df(nConts_in, natoms_in, Pa, Pb, Ca, nOccA, Cb, nOccB, &
                                      force_J, force_Ka, force_Kb, need_k)
use MOL_info, only: nAtoms
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in, natoms_in, nOccA, nOccB
real(8),intent(in) :: Pa(nConts_in,nConts_in), Pb(nConts_in,nConts_in)
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: force_J(natoms_in,3), force_Ka(natoms_in,3), force_Kb(natoms_in,3)
logical,intent(in) :: need_k

real(8),allocatable :: Ptot(:,:), bvec(:), tvec(:), cvec(:), bprime(:), Ptot_compact(:)
real(8),allocatable :: Wstack(:,:,:), Fstack(:,:,:)
real(8),allocatable :: h_ij_a(:,:,:), h_ij_b(:,:,:)
real(8),allocatable :: dmax_shell(:,:)
integer :: A, nT, it, idxKa, idxKb, slot
logical :: restricted_kb, pure_j_only

if (df_force_direct_mode) then
   if (.not. need_k) then
      call integrals_force_df_direct_J(nConts_in, natoms_in, Pa+Pb, force_J)
      force_Ka = 0.0d0
      force_Kb = 0.0d0
   else
      call integrals_force_df_direct_JK(nConts_in, natoms_in, Pa, Pb, Ca, nOccA, Cb, nOccB, &
                                         force_J, force_Ka, force_Kb)
   endif
   return
endif

call ensure_df_built()
if (.not. df_Minvhalf_built) call build_df_metric_invhalf()

restricted_kb = .false.
if (need_k) then
   restricted_kb = (nOccB .gt. 0) .and. (nOccB .eq. nOccA) .and. &
                   all(Ca(:,1:nOccA) .eq. Cb(:,1:nOccB))
endif
nT = 1
idxKa = 0
idxKb = 0
if (need_k .and. nOccA .gt. 0) then
   nT = nT + 1
   idxKa = nT
endif
if (need_k .and. (nOccB .gt. 0) .and. (.not. restricted_kb)) then
   nT = nT + 1
   idxKb = nT
endif
allocate(Wstack(nContsAux,nContsAux,nT))

allocate(Ptot(nConts_in,nConts_in))
Ptot = Pa + Pb

pure_j_only = (idxKa .eq. 0) .and. (idxKb .eq. 0)
if (pure_j_only) then
   allocate(dmax_shell(0:nBases-1,0:nBases-1))
   call compute_shell_density_bound(nConts_in, Ptot, dmax_shell)
else
   allocate(dmax_shell(0,0))
endif

allocate(Ptot_compact(df_total_pairs))
!$omp parallel do private(slot)
do slot = 1,df_total_pairs
   Ptot_compact(slot) = Ptot(df_pair_row(slot), df_pair_col(slot))
enddo
!$omp end parallel do
call openblas_set_num_threads(omp_get_max_threads())
if (dfK_built) then
   allocate(cvec(nContsAux))
   allocate(bprime(nContsAux))
   call dgemv('T', df_total_pairs, nContsAux, 1.0d0, dfK_compact, df_total_pairs, &
              Ptot_compact, 1, 0.0d0, bprime, 1)
   call dgemv('N', nContsAux, nContsAux, 1.0d0, df_Minvhalf, nContsAux, bprime, 1, 0.0d0, cvec, 1)
   deallocate(bprime)
else
   allocate(cvec(nContsAux))
   if (df_direct_cvec_cache_valid) then
      if (all(df_direct_cvec_cache_Ptot .eq. Ptot)) then
         cvec = df_direct_cvec_cache
      else
         call df_force_cvec_dfB_fallback(Ptot_compact, cvec)
      endif
   else
      call df_force_cvec_dfB_fallback(Ptot_compact, cvec)
   endif
endif
deallocate(Ptot_compact)
call openblas_set_num_threads(1)

!$omp parallel do
do A = 1,nContsAux
   Wstack(:,A,1) = cvec(:)*cvec(A)
enddo
!$omp end parallel do

if (need_k) then
   if (nOccA .gt. 0) then
      allocate(h_ij_a(nOccA,nOccA,nContsAux))
      call df_build_H_K(nConts_in, Ca, nOccA, h_ij_a, Wstack(:,:,idxKa))
   endif
   if ((nOccB .gt. 0) .and. (.not. restricted_kb)) then
      allocate(h_ij_b(nOccB,nOccB,nContsAux))
      call df_build_H_K(nConts_in, Cb, nOccB, h_ij_b, Wstack(:,:,idxKb))
   endif
endif

allocate(Fstack(natoms_in,3,nT))
Fstack = 0.0d0
call prof_start("df_force_contract")
block
   use mod_mem_predict, only: TSTACK_BLOCK_BYTES, TSTACK_MIN_BLOCK_WIDTH
   integer(8) :: ip1_opt_shared, ip2_opt_shared
   integer :: ksh, dk, maxP, max_rows, blk_lo, blk_hi, blk_width, blk_offset
   integer(8) :: row_cost_bytes
   real(8),allocatable :: Tstack_blk(:,:,:,:)

   maxP = 0
   do ksh = nBases,nBases+nBasesAux-1
      maxP = max(maxP, cgto_engine(ksh,basDF))
   enddo
   row_cost_bytes = int(nConts_in,8)*int(nConts_in,8)*int(nT,8)*8_8
   max_rows = int(TSTACK_BLOCK_BYTES/real(row_cost_bytes,8))
   max_rows = max(max_rows, maxP, TSTACK_MIN_BLOCK_WIDTH)
   max_rows = min(max_rows, nContsAux)

   call threec2e_ip1_optimizer_engine(ip1_opt_shared, atm, nAtoms, basDF, nBases+nBasesAux, envDF)
   call threec2e_ip2_optimizer_engine(ip2_opt_shared, atm, nAtoms, basDF, nBases+nBasesAux, envDF)

   blk_lo = nBases
   blk_offset = 0
   do while (blk_lo .le. nBases+nBasesAux-1)
      blk_hi = blk_lo
      blk_width = cgto_engine(blk_lo,basDF)
      do while (blk_hi+1 .le. nBases+nBasesAux-1)
         dk = cgto_engine(blk_hi+1,basDF)
         if (blk_width+dk .gt. max_rows) exit
         blk_hi = blk_hi+1
         blk_width = blk_width+dk
      enddo

      allocate(Tstack_blk(nConts_in,nConts_in,blk_width,nT))
      !$omp parallel do
      do A = 1,blk_width
         Tstack_blk(:,:,A,1) = Ptot(:,:)*cvec(blk_offset+A)
      enddo
      !$omp end parallel do
      if (need_k) then
         if (nOccA .gt. 0) call df_H_K_backtransform_block(nConts_in, nOccA, Ca, h_ij_a, &
                                 blk_offset+1, blk_offset+blk_width, Tstack_blk(:,:,:,idxKa))
         if ((nOccB .gt. 0) .and. (.not. restricted_kb)) then
            call df_H_K_backtransform_block(nConts_in, nOccB, Cb, h_ij_b, &
                     blk_offset+1, blk_offset+blk_width, Tstack_blk(:,:,:,idxKb))
         endif
      endif

      call df_force_3c2e_contract(nConts_in, nT, blk_width, blk_lo, blk_hi, blk_offset, &
                                   ip1_opt_shared, ip2_opt_shared, &
                                   Tstack_blk, natoms_in, Fstack, pure_j_only, dmax_shell)
      deallocate(Tstack_blk)

      blk_offset = blk_offset + blk_width
      blk_lo = blk_hi + 1
   enddo
   call cintdel_optimizer(ip1_opt_shared)
   call cintdel_optimizer(ip2_opt_shared)
end block
call df_force_2c2e_contract(nT, Wstack, natoms_in, Fstack)
call prof_stop("df_force_contract")
deallocate(Ptot, cvec, Wstack, dmax_shell)
if (allocated(h_ij_a)) deallocate(h_ij_a)
if (allocated(h_ij_b)) deallocate(h_ij_b)

force_J = Fstack(:,:,1)
force_Ka = 0.0d0
if (idxKa .gt. 0) force_Ka = Fstack(:,:,idxKa)
force_Kb = 0.0d0
if (need_k .and. nOccB .gt. 0) then
   if (restricted_kb) then
      force_Kb = force_Ka
   else
      force_Kb = Fstack(:,:,idxKb)
   endif
endif
deallocate(Fstack)

end subroutine integrals_force_df

subroutine df_force_cvec_dfB_fallback(Ptot_compact, cvec)
implicit none
real(8),intent(in) :: Ptot_compact(df_total_pairs)
real(8),intent(out) :: cvec(nContsAux)
real(8),allocatable :: bvec(:), tvec(:)
allocate(bvec(nContsAux), tvec(nContsAux))
call dgemv('T', df_total_pairs, nContsAux, 1.0d0, dfB_compact, df_total_pairs, &
           Ptot_compact, 1, 0.0d0, bvec, 1)
call dgemv('T', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, bvec, 1, 0.0d0, tvec, 1)
tvec = tvec * df_evalinv
call dgemv('N', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, tvec, 1, 0.0d0, cvec, 1)
deallocate(bvec,tvec)
end subroutine df_force_cvec_dfB_fallback

subroutine df_build_H_K(nConts_in, C, nOcc, h_ij, W_K)
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in, nOcc
real(8),intent(in) :: C(nConts_in,nConts_in)
real(8),intent(out) :: h_ij(nOcc,nOcc,nContsAux), W_K(nContsAux,nContsAux)

real(8),allocatable :: Bmi(:,:,:), g_ij(:,:,:)
integer :: A

call prof_start("dfK_halftransform")
allocate(Bmi(nOcc,nConts_in,nContsAux))
block
   real(8),allocatable :: C_gathered(:,:)
   integer :: m, rs, rc, nn, slot2
   !$omp parallel private(C_gathered,m,rs,rc,nn,slot2)
   allocate(C_gathered(df_max_row_count,nOcc))
   !$omp do schedule(guided)
   do m = 1,nConts_in
      rs = df_row_start(m)
      rc = df_row_count(m)
      if (rc .gt. 0) then
         do nn = 1,rc
            slot2 = rs+nn-1
            C_gathered(nn,1:nOcc) = C(df_pair_col(slot2), 1:nOcc)
         enddo
         call dgemm('T','N', nOcc, nContsAux, rc, 1.0d0, &
                    C_gathered, df_max_row_count, dfK_compact(rs,1), df_total_pairs, 0.0d0, &
                    Bmi(1,m,1), nOcc*nConts_in)
      else
         Bmi(:,m,:) = 0.0d0
      endif
   enddo
   !$omp end do
   deallocate(C_gathered)
   !$omp end parallel
end block

allocate(g_ij(nOcc,nOcc,nContsAux))
!$omp parallel do
do A = 1,nContsAux
   call dgemm('N','N', nOcc, nOcc, nConts_in, 1.0d0, &
              Bmi(:,:,A), nOcc, C(:,1:nOcc), nConts_in, 0.0d0, g_ij(:,:,A), nOcc)
enddo
!$omp end parallel do
deallocate(Bmi)

call openblas_set_num_threads(omp_get_max_threads())
call dgemm('N','N', nOcc*nOcc, nContsAux, nContsAux, 1.0d0, &
           g_ij, nOcc*nOcc, df_Minvhalf, nContsAux, 0.0d0, h_ij, nOcc*nOcc)
call openblas_set_num_threads(1)
deallocate(g_ij)

call openblas_set_num_threads(omp_get_max_threads())
call dgemm('T','N', nContsAux, nContsAux, nOcc*nOcc, 1.0d0, &
           h_ij, nOcc*nOcc, h_ij, nOcc*nOcc, 0.0d0, W_K, nContsAux)
call openblas_set_num_threads(1)
call prof_stop("dfK_halftransform")

end subroutine df_build_H_K

subroutine df_H_K_backtransform_block(nConts_in, nOcc, C, h_ij, A_lo, A_hi, H_ao_blk)
implicit none
integer,intent(in) :: nConts_in, nOcc, A_lo, A_hi
real(8),intent(in) :: C(nConts_in,nConts_in), h_ij(nOcc,nOcc,nContsAux)
real(8),intent(out) :: H_ao_blk(nConts_in,nConts_in,A_hi-A_lo+1)
real(8),allocatable :: temp(:,:)
integer :: A

!$omp parallel private(A,temp)
allocate(temp(nConts_in,nOcc))
!$omp do
do A = A_lo,A_hi
   call dgemm('N','N', nConts_in, nOcc, nOcc, 1.0d0, &
              C(:,1:nOcc), nConts_in, h_ij(:,:,A), nOcc, 0.0d0, temp, nConts_in)
   call dgemm('N','T', nConts_in, nConts_in, nOcc, 1.0d0, &
              temp, nConts_in, C(:,1:nOcc), nConts_in, 0.0d0, H_ao_blk(:,:,A-A_lo+1), nConts_in)
enddo
!$omp end do
deallocate(temp)
!$omp end parallel

end subroutine df_H_K_backtransform_block

subroutine df_force_3c2e_contract(nConts_in, nT, nAux_blk, ksh_lo, ksh_hi, aux_blk_offset, &
                                   ip1_opt, ip2_opt, T, natoms_in, force_out, pure_j_only, dmax_shell)
use MOL_info, only: nAtoms
implicit none
integer,intent(in) :: nConts_in, natoms_in, nT, nAux_blk, ksh_lo, ksh_hi, aux_blk_offset
integer(8),intent(in) :: ip1_opt, ip2_opt
real(8),intent(in) :: T(nConts_in,nConts_in,nAux_blk,nT)
real(8),intent(inout) :: force_out(natoms_in,3,nT)
logical,intent(in) :: pure_j_only
real(8),intent(in) :: dmax_shell(0:,0:)

integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP,offP_local,it
integer :: shls(4), atomI, atomJ, atomK
real(8),allocatable,target :: buf_pool(:)
real(8),pointer :: buf(:,:,:,:)
real(8) :: val(3), pair_weight
integer,allocatable :: ao_offset0(:), aux_offset0(:)
integer :: maxDao, maxDk_here

allocate(ao_offset0(0:nBases-1))
offi = 0
do ish = 0,nBases-1
   ao_offset0(ish) = offi
   offi = offi + cgto_engine(ish, basDF)
enddo
allocate(aux_offset0(nBases:nBases+nBasesAux-1))
offP = 0
do ksh = nBases,nBases+nBasesAux-1
   aux_offset0(ksh) = offP
   offP = offP + cgto_engine(ksh, basDF)
enddo

maxDao = 0
do ish = 0,nBases-1
   maxDao = max(maxDao, cgto_engine(ish, basDF))
enddo
maxDk_here = 0
do ksh = ksh_lo,ksh_hi
   maxDk_here = max(maxDk_here, cgto_engine(ksh, basDF))
enddo

!$omp parallel do private(ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP,offP_local,it) &
!$omp&   private(shls,atomI,atomJ,atomK,buf,buf_pool,val,pair_weight) &
!$omp&   reduction(+:force_out) schedule(dynamic)
do ish = 0,nBases-1
   if (.not. allocated(buf_pool)) allocate(buf_pool(maxDao*maxDao*maxDk_here*3))
   di = cgto_engine(ish, basDF)
   atomI = basDF(1,ish+1) + 1
   offi = ao_offset0(ish)
   do jsh = 0,ish
      dj = cgto_engine(jsh, basDF)
      atomJ = basDF(1,jsh+1) + 1
      offj = ao_offset0(jsh)
      pair_weight = merge(1.0d0, 2.0d0, ish .eq. jsh)
      do ksh = ksh_lo,ksh_hi
         dk = cgto_engine(ksh, basDF)
         atomK = basDF(1,ksh+1) + 1
         offP = aux_offset0(ksh)
         offP_local = offP - aux_blk_offset

         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) then
            cycle
         endif

         if (pure_j_only) then
            if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1)*dmax_shell(ish,jsh) &
                .lt. density_screen_cutoff) cycle
         endif

         shls(1) = ish
         shls(2) = jsh
         shls(3) = ksh
         buf(1:di,1:dj,1:dk,1:3) => buf_pool(1:di*dj*dk*3)
         call threec2e_ip1_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip1_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf(p,q,r,:)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  do it = 1,nT
                     force_out(atomI,:,it) = force_out(atomI,:,it) + 2.0d0*T(offi+p,offj+q,offP_local+r,it)*val
                  enddo
               enddo
            enddo
         enddo

         if (jsh .ne. ish) then
            shls(1) = jsh
            shls(2) = ish
            shls(3) = ksh
            buf(1:dj,1:di,1:dk,1:3) => buf_pool(1:dj*di*dk*3)
            call threec2e_ip1_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip1_opt)
            do r = 1,dk
               do p = 1,di
                  do q = 1,dj
                     val = buf(q,p,r,:)*NorVEC(offj+q)*NorVEC(offi+p)*NorVECAux(offP+r)
                     do it = 1,nT
                        force_out(atomJ,:,it) = force_out(atomJ,:,it) + pair_weight*T(offj+q,offi+p,offP_local+r,it)*val
                     enddo
                  enddo
               enddo
            enddo
         endif

         shls(1) = ish
         shls(2) = jsh
         shls(3) = ksh
         buf(1:di,1:dj,1:dk,1:3) => buf_pool(1:di*dj*dk*3)
         call threec2e_ip2_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip2_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf(p,q,r,:)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  do it = 1,nT
                     force_out(atomK,:,it) = force_out(atomK,:,it) + pair_weight*T(offi+p,offj+q,offP_local+r,it)*val
                  enddo
               enddo
            enddo
         enddo

      enddo
   enddo
enddo
!$omp end parallel do
deallocate(ao_offset0, aux_offset0)

end subroutine df_force_3c2e_contract

subroutine df_force_3c2e_contract_lr(nConts_in, nT, nAux_blk, ksh_lo, ksh_hi, aux_blk_offset, T, natoms_in, force_out)
use MOL_info, only: nAtoms
use mod_exchange, only: RS_omega
implicit none
integer,intent(in) :: nConts_in, natoms_in, nT, nAux_blk, ksh_lo, ksh_hi, aux_blk_offset
real(8),intent(in) :: T(nConts_in,nConts_in,nAux_blk,nT)
real(8),intent(inout) :: force_out(natoms_in,3,nT)

integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP,offP_local,it
integer :: shls(4), atomI, atomJ, atomK
real(8),allocatable :: buf(:,:,:,:)
real(8) :: val(3), pair_weight
integer(8) :: ip1_opt, ip2_opt
integer,allocatable :: ao_offset0(:), aux_offset0(:)

envDF(9) = RS_omega
call threec2e_ip1_optimizer_engine(ip1_opt, atm, nAtoms, basDF, nBases+nBasesAux, envDF)
call threec2e_ip2_optimizer_engine(ip2_opt, atm, nAtoms, basDF, nBases+nBasesAux, envDF)

allocate(ao_offset0(0:nBases-1))
offi = 0
do ish = 0,nBases-1
   ao_offset0(ish) = offi
   offi = offi + cgto_engine(ish, basDF)
enddo
allocate(aux_offset0(nBases:nBases+nBasesAux-1))
offP = 0
do ksh = nBases,nBases+nBasesAux-1
   aux_offset0(ksh) = offP
   offP = offP + cgto_engine(ksh, basDF)
enddo

!$omp parallel do private(ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP,offP_local,it) &
!$omp&   private(shls,atomI,atomJ,atomK,buf,val,pair_weight) &
!$omp&   reduction(+:force_out) schedule(dynamic)
do ish = 0,nBases-1
   di = cgto_engine(ish, basDF)
   atomI = basDF(1,ish+1) + 1
   offi = ao_offset0(ish)
   do jsh = 0,ish
      dj = cgto_engine(jsh, basDF)
      atomJ = basDF(1,jsh+1) + 1
      offj = ao_offset0(jsh)
      pair_weight = merge(1.0d0, 2.0d0, ish .eq. jsh)
      do ksh = ksh_lo,ksh_hi
         dk = cgto_engine(ksh, basDF)
         atomK = basDF(1,ksh+1) + 1
         offP = aux_offset0(ksh)
         offP_local = offP - aux_blk_offset

         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) then
            cycle
         endif

         shls(1) = ish
         shls(2) = jsh
         shls(3) = ksh
         allocate(buf(di,dj,dk,3))
         call threec2e_ip1_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip1_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf(p,q,r,:)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  do it = 1,nT
                     force_out(atomI,:,it) = force_out(atomI,:,it) + 2.0d0*T(offi+p,offj+q,offP_local+r,it)*val
                  enddo
               enddo
            enddo
         enddo
         deallocate(buf)

         if (jsh .ne. ish) then
            shls(1) = jsh
            shls(2) = ish
            shls(3) = ksh
            allocate(buf(dj,di,dk,3))
            call threec2e_ip1_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip1_opt)
            do r = 1,dk
               do p = 1,di
                  do q = 1,dj
                     val = buf(q,p,r,:)*NorVEC(offj+q)*NorVEC(offi+p)*NorVECAux(offP+r)
                     do it = 1,nT
                        force_out(atomJ,:,it) = force_out(atomJ,:,it) + pair_weight*T(offj+q,offi+p,offP_local+r,it)*val
                     enddo
                  enddo
               enddo
            enddo
            deallocate(buf)
         endif

         shls(1) = ish
         shls(2) = jsh
         shls(3) = ksh
         allocate(buf(di,dj,dk,3))
         call threec2e_ip2_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip2_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf(p,q,r,:)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  do it = 1,nT
                     force_out(atomK,:,it) = force_out(atomK,:,it) + pair_weight*T(offi+p,offj+q,offP_local+r,it)*val
                  enddo
               enddo
            enddo
         enddo
         deallocate(buf)

      enddo
   enddo
enddo
!$omp end parallel do
deallocate(ao_offset0, aux_offset0)
call cintdel_optimizer(ip1_opt)
call cintdel_optimizer(ip2_opt)
envDF(9) = 0.0d0
end subroutine df_force_3c2e_contract_lr

subroutine df_build_g_ij_raw(nConts_in, C, nOcc, Braw_compact, g_ij)
implicit none
integer,intent(in) :: nConts_in, nOcc
real(8),intent(in) :: C(nConts_in,nConts_in)
real(8),intent(in) :: Braw_compact(df_total_pairs,nContsAux)
real(8),intent(out) :: g_ij(nOcc,nOcc,nContsAux)

real(8),allocatable :: Bmi(:,:,:)
integer :: A

allocate(Bmi(nOcc,nConts_in,nContsAux))
block
   real(8),allocatable :: C_gathered(:,:)
   integer :: m, rs, rc, nn, slot2
   !$omp parallel private(C_gathered,m,rs,rc,nn,slot2)
   allocate(C_gathered(df_max_row_count,nOcc))
   !$omp do schedule(guided)
   do m = 1,nConts_in
      rs = df_row_start(m)
      rc = df_row_count(m)
      if (rc .gt. 0) then
         do nn = 1,rc
            slot2 = rs+nn-1
            C_gathered(nn,1:nOcc) = C(df_pair_col(slot2), 1:nOcc)
         enddo
         call dgemm('T','N', nOcc, nContsAux, rc, 1.0d0, &
                    C_gathered, df_max_row_count, Braw_compact(rs,1), df_total_pairs, 0.0d0, &
                    Bmi(1,m,1), nOcc*nConts_in)
      else
         Bmi(:,m,:) = 0.0d0
      endif
   enddo
   !$omp end do
   deallocate(C_gathered)
   !$omp end parallel
end block

!$omp parallel do
do A = 1,nContsAux
   call dgemm('N','N', nOcc, nOcc, nConts_in, 1.0d0, &
              Bmi(:,:,A), nOcc, C(:,1:nOcc), nConts_in, 0.0d0, g_ij(:,:,A), nOcc)
enddo
!$omp end parallel do
deallocate(Bmi)
end subroutine df_build_g_ij_raw

module subroutine integrals_force_df_lr(nConts_in, natoms_in, Ca, nOccA, Cb, nOccB, &
                                         force_Ka_lr, force_Kb_lr)
use MOL_info, only: nAtoms
use mod_exchange, only: RS_omega
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in, natoms_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: force_Ka_lr(natoms_in,3), force_Kb_lr(natoms_in,3)

logical :: restricted_kb
real(8),allocatable :: h_ij_a(:,:,:), h_ij_b(:,:,:)
real(8),allocatable :: c_ij_lr_a(:,:,:), c_ij_lr_b(:,:,:)
real(8),allocatable :: Fstack(:,:,:)
integer :: nT, idxA, idxB

force_Ka_lr = 0.0d0
force_Kb_lr = 0.0d0
if (.not. (RS_omega .gt. 0.0d0)) return

if (df_direct_mode) then
   print *, "ERROR: CAM-B3LYP RI-K DIRECT-mode analytic gradient not implemented ", &
            "(STORE-only scope limit - energy-only DIRECT-mode CAM-B3LYP/K=RI runs are fine, ", &
            "but this molecule needs DIRECT mode AND a force calc was requested)"
   stop 1
endif

call ensure_df_built()
if (.not. df_Minvhalf_built) call build_df_metric_invhalf()
if (.not. dfK_built) then
   call build_df_exchange_metric()
   dfK_built = .true.
endif

restricted_kb = (nOccB .gt. 0) .and. (nOccB .eq. nOccA) .and. &
                all(Ca(:,1:nOccA) .eq. Cb(:,1:nOccB))

nT = 0
idxA = 0
idxB = 0
if (nOccA .gt. 0) then
   nT = nT + 1
   idxA = nT
endif
if ((nOccB .gt. 0) .and. (.not. restricted_kb)) then
   nT = nT + 1
   idxB = nT
endif
if (nT .eq. 0) return

allocate(Fstack(natoms_in,3,nT))
Fstack = 0.0d0

call prof_start("df_force_contract_lr")

if (idxA .gt. 0) then
   allocate(h_ij_a(nOccA,nOccA,nContsAux))
   block
      real(8),allocatable :: Wscratch(:,:)
      allocate(Wscratch(nContsAux,nContsAux))
      call df_build_H_K(nConts_in, Ca, nOccA, h_ij_a, Wscratch)
      deallocate(Wscratch)
   end block

   allocate(c_ij_lr_a(nOccA,nOccA,nContsAux))
   call df_build_g_ij_raw(nConts_in, Ca, nOccA, dfB_compact_LR, c_ij_lr_a)
   block
      real(8),allocatable :: tmp(:,:,:)
      allocate(tmp(nOccA,nOccA,nContsAux))
      call openblas_set_num_threads(omp_get_max_threads())
      call dgemm('N','N', nOccA*nOccA, nContsAux, nContsAux, 1.0d0, &
                 c_ij_lr_a, nOccA*nOccA, df_Minvhalf, nContsAux, 0.0d0, tmp, nOccA*nOccA)
      call dgemm('N','N', nOccA*nOccA, nContsAux, nContsAux, 1.0d0, &
                 tmp, nOccA*nOccA, df_Minvhalf, nContsAux, 0.0d0, c_ij_lr_a, nOccA*nOccA)
      call openblas_set_num_threads(1)
      deallocate(tmp)
   end block
endif

if (idxB .gt. 0) then
   allocate(h_ij_b(nOccB,nOccB,nContsAux))
   block
      real(8),allocatable :: Wscratch(:,:)
      allocate(Wscratch(nContsAux,nContsAux))
      call df_build_H_K(nConts_in, Cb, nOccB, h_ij_b, Wscratch)
      deallocate(Wscratch)
   end block

   allocate(c_ij_lr_b(nOccB,nOccB,nContsAux))
   call df_build_g_ij_raw(nConts_in, Cb, nOccB, dfB_compact_LR, c_ij_lr_b)
   block
      real(8),allocatable :: tmp(:,:,:)
      allocate(tmp(nOccB,nOccB,nContsAux))
      call openblas_set_num_threads(omp_get_max_threads())
      call dgemm('N','N', nOccB*nOccB, nContsAux, nContsAux, 1.0d0, &
                 c_ij_lr_b, nOccB*nOccB, df_Minvhalf, nContsAux, 0.0d0, tmp, nOccB*nOccB)
      call dgemm('N','N', nOccB*nOccB, nContsAux, nContsAux, 1.0d0, &
                 tmp, nOccB*nOccB, df_Minvhalf, nContsAux, 0.0d0, c_ij_lr_b, nOccB*nOccB)
      call openblas_set_num_threads(1)
      deallocate(tmp)
   end block
endif

block
   real(8),allocatable :: Wstack2c(:,:,:), Praw(:,:)
   allocate(Wstack2c(nContsAux,nContsAux,nT))
   if (idxA .gt. 0) then
      allocate(Praw(nContsAux,nContsAux))
      call openblas_set_num_threads(omp_get_max_threads())
      call dgemm('T','N', nContsAux, nContsAux, nOccA*nOccA, 1.0d0, &
                 h_ij_a, nOccA*nOccA, c_ij_lr_a, nOccA*nOccA, 0.0d0, Praw, nContsAux)
      call openblas_set_num_threads(1)
      Wstack2c(:,:,idxA) = Praw + transpose(Praw)
      deallocate(Praw)
   endif
   if (idxB .gt. 0) then
      allocate(Praw(nContsAux,nContsAux))
      call openblas_set_num_threads(omp_get_max_threads())
      call dgemm('T','N', nContsAux, nContsAux, nOccB*nOccB, 1.0d0, &
                 h_ij_b, nOccB*nOccB, c_ij_lr_b, nOccB*nOccB, 0.0d0, Praw, nContsAux)
      call openblas_set_num_threads(1)
      Wstack2c(:,:,idxB) = Praw + transpose(Praw)
      deallocate(Praw)
   endif
   call df_force_2c2e_contract(nT, Wstack2c, natoms_in, Fstack)
   deallocate(Wstack2c)
end block

block
   use mod_mem_predict, only: TSTACK_BLOCK_BYTES
   integer :: ksh, dk, maxP, max_rows, blk_lo, blk_hi, blk_width, blk_offset, pass
   integer(8) :: row_cost_bytes
   integer(8) :: ip1_lr_opt, ip2_lr_opt
   real(8),allocatable :: Tblk(:,:,:,:)
   real(8) :: dmax_shell_dummy(0,0)
   logical :: lr_opt_built

   maxP = 0
   do ksh = nBases,nBases+nBasesAux-1
      maxP = max(maxP, cgto_engine(ksh,basDF))
   enddo
   call threec2e_ip1_optimizer_engine(ip1_lr_opt, atm, nAtoms, basDF, nBases+nBasesAux, envDF)
   call threec2e_ip2_optimizer_engine(ip2_lr_opt, atm, nAtoms, basDF, nBases+nBasesAux, envDF)
   row_cost_bytes = int(nConts_in,8)*int(nConts_in,8)*int(nT,8)*8_8
   max_rows = int(TSTACK_BLOCK_BYTES/real(row_cost_bytes,8))
   max_rows = max(max_rows, maxP)
   max_rows = min(max_rows, nContsAux)

   do pass = 1,2
      blk_lo = nBases
      blk_offset = 0
      do while (blk_lo .le. nBases+nBasesAux-1)
         blk_hi = blk_lo
         blk_width = cgto_engine(blk_lo,basDF)
         do while (blk_hi+1 .le. nBases+nBasesAux-1)
            dk = cgto_engine(blk_hi+1,basDF)
            if (blk_width+dk .gt. max_rows) exit
            blk_hi = blk_hi+1
            blk_width = blk_width+dk
         enddo

         allocate(Tblk(nConts_in,nConts_in,blk_width,nT))
         if (pass .eq. 1) then
            if (idxA .gt. 0) call df_H_K_backtransform_block(nConts_in, nOccA, Ca, c_ij_lr_a, &
                                   blk_offset+1, blk_offset+blk_width, Tblk(:,:,:,idxA))
            if (idxB .gt. 0) call df_H_K_backtransform_block(nConts_in, nOccB, Cb, c_ij_lr_b, &
                                   blk_offset+1, blk_offset+blk_width, Tblk(:,:,:,idxB))
            call df_force_3c2e_contract(nConts_in, nT, blk_width, blk_lo, blk_hi, blk_offset, &
                                         ip1_lr_opt, ip2_lr_opt, Tblk, natoms_in, Fstack, &
                                         .false., dmax_shell_dummy)
         else
            if (idxA .gt. 0) call df_H_K_backtransform_block(nConts_in, nOccA, Ca, h_ij_a, &
                                   blk_offset+1, blk_offset+blk_width, Tblk(:,:,:,idxA))
            if (idxB .gt. 0) call df_H_K_backtransform_block(nConts_in, nOccB, Cb, h_ij_b, &
                                   blk_offset+1, blk_offset+blk_width, Tblk(:,:,:,idxB))
            call df_force_3c2e_contract_lr(nConts_in, nT, blk_width, blk_lo, blk_hi, blk_offset, &
                                            Tblk, natoms_in, Fstack)
         endif
         deallocate(Tblk)

         blk_offset = blk_offset + blk_width
         blk_lo = blk_hi + 1
      enddo
   enddo
end block
call prof_stop("df_force_contract_lr")

if (idxA .gt. 0) force_Ka_lr = Fstack(:,:,idxA)
if (idxB .gt. 0) then
   force_Kb_lr = Fstack(:,:,idxB)
else if (restricted_kb) then
   force_Kb_lr = force_Ka_lr
endif

deallocate(Fstack)
if (allocated(h_ij_a)) deallocate(h_ij_a)
if (allocated(h_ij_b)) deallocate(h_ij_b)
if (allocated(c_ij_lr_a)) deallocate(c_ij_lr_a)
if (allocated(c_ij_lr_b)) deallocate(c_ij_lr_b)
end subroutine integrals_force_df_lr

module subroutine df_force_2c2e_contract(nT, W, natoms_in, force_out)
use MOL_info, only: nAtoms
implicit none
integer,intent(in) :: natoms_in, nT
real(8),intent(in) :: W(nContsAux,nContsAux,nT)
real(8),intent(inout) :: force_out(natoms_in,3,nT)

integer :: ish,jsh,p,q,di,dj,offP,offQ,atomI,atomJ,it
integer :: shls(4)
real(8),allocatable :: buf(:,:,:)
real(8) :: val(3), half_perm
integer(8) :: ip1_opt, ip2_opt
integer,allocatable :: aux_offset0(:)

call twoc2e_ip1_optimizer_engine(ip1_opt, atm, nAtoms, basDF, nBases+nBasesAux, envDF)
call twoc2e_ip2_optimizer_engine(ip2_opt, atm, nAtoms, basDF, nBases+nBasesAux, envDF)

allocate(aux_offset0(nBases:nBases+nBasesAux-1))
offP = 0
do ish = nBases,nBases+nBasesAux-1
   aux_offset0(ish) = offP
   offP = offP + cgto_engine(ish, basDF)
enddo

!$omp parallel do private(ish,jsh,p,q,di,dj,offP,offQ,atomI,atomJ,it) &
!$omp&   private(shls,buf,val,half_perm) reduction(+:force_out) schedule(dynamic)
do ish = nBases,nBases+nBasesAux-1
   di = cgto_engine(ish, basDF)
   atomI = basDF(1,ish+1) + 1
   offP = aux_offset0(ish)
   do jsh = nBases,ish
      dj = cgto_engine(jsh, basDF)
      atomJ = basDF(1,jsh+1) + 1
      offQ = aux_offset0(jsh)
      half_perm = merge(0.5d0, 1.0d0, ish .eq. jsh)
      shls(1) = ish
      shls(2) = jsh

      if (aux_shell_bound(ish-nBases+1)*aux_shell_bound(jsh-nBases+1) &
          .lt. SCHWARZ_CUTOFF) then
         cycle
      endif

      allocate(buf(di,dj,3))
      call twoc2e_ip1_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip1_opt)
      do q = 1,dj
         do p = 1,di
            val = buf(p,q,:)*NorVECAux(offP+p)*NorVECAux(offQ+q)
            do it = 1,nT
               force_out(atomI,:,it) = force_out(atomI,:,it) - half_perm*W(offP+p,offQ+q,it)*val
            enddo
         enddo
      enddo
      deallocate(buf)

      allocate(buf(di,dj,3))
      call twoc2e_ip2_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip2_opt)
      do q = 1,dj
         do p = 1,di
            val = buf(p,q,:)*NorVECAux(offP+p)*NorVECAux(offQ+q)
            do it = 1,nT
               force_out(atomJ,:,it) = force_out(atomJ,:,it) - half_perm*W(offP+p,offQ+q,it)*val
            enddo
         enddo
      enddo
      deallocate(buf)

   enddo
enddo
!$omp end parallel do
deallocate(aux_offset0)
call cintdel_optimizer(ip1_opt)
call cintdel_optimizer(ip2_opt)
end subroutine df_force_2c2e_contract

end submodule df_force_store_impl
