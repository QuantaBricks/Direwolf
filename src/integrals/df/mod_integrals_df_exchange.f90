! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Density-fitting (RI-K) exchange matrix build, including the long-range (range-separated) variant.

submodule (mod_integrals) df_exchange_impl
implicit none
contains

module subroutine build_df_exchange_metric()
use MOL_info
use omp_lib, only: omp_get_max_threads
implicit none

if (.not. df_Minvhalf_built) call build_df_metric_invhalf()

block
   integer,parameter :: DFK_XFORM_CHUNK = 4096
   real(8),allocatable :: dfB_blk(:,:), dfK_blk(:,:)
   integer,allocatable :: primary_list(:)
   integer :: slot, nprimary, p, cstart, cend, csize, blksz
   nprimary = count(df_pair_row .le. df_pair_col)
   allocate(primary_list(nprimary))
   p = 0
   do slot = 1,df_total_pairs
      if (df_pair_row(slot) .le. df_pair_col(slot)) then
         p = p+1
         primary_list(p) = slot
      endif
   enddo
   blksz = min(DFK_XFORM_CHUNK, max(nprimary,1))
   allocate(dfB_blk(blksz,nContsAux), dfK_blk(blksz,nContsAux))
   do cstart = 1,nprimary,DFK_XFORM_CHUNK
      cend = min(cstart+DFK_XFORM_CHUNK-1, nprimary)
      csize = cend-cstart+1
      !$omp parallel do private(p)
      do p = cstart,cend
         dfB_blk(p-cstart+1,:) = dfB_compact(primary_list(p),:)
      enddo
      !$omp end parallel do
      call openblas_set_num_threads(omp_get_max_threads())
      call dgemm('N','N', csize, nContsAux, nContsAux, 1.0d0, &
                 dfB_blk, blksz, df_Minvhalf, nContsAux, 0.0d0, dfK_blk, blksz)
      call openblas_set_num_threads(1)
      !$omp parallel do private(p,slot)
      do p = cstart,cend
         slot = primary_list(p)
         dfB_compact(slot,:) = dfK_blk(p-cstart+1,:)
         dfB_compact(df_pair_mirror(slot),:) = dfK_blk(p-cstart+1,:)
      enddo
      !$omp end parallel do
   enddo
   deallocate(dfB_blk,dfK_blk,primary_list)
   call move_alloc(dfB_compact, dfK_compact)
end block
end subroutine build_df_exchange_metric

module subroutine build_df_exchange_metric_full_lr()
use MOL_info
use omp_lib, only: omp_get_max_threads
implicit none
integer,parameter :: DFC_XFORM_CHUNK = 4096
integer :: cstart, cend, csize

if (.not. dfK_built) then
   call build_df_exchange_metric()
   dfK_built = .true.
endif
if (.not. df_Minvhalf_built) call build_df_metric_invhalf()

allocate(dfC_full_compact(df_total_pairs,nContsAux))
do cstart = 1,df_total_pairs,DFC_XFORM_CHUNK
   cend = min(cstart+DFC_XFORM_CHUNK-1, df_total_pairs)
   csize = cend-cstart+1
   call openblas_set_num_threads(omp_get_max_threads())
   call dgemm('N','N', csize, nContsAux, nContsAux, 1.0d0, &
              dfK_compact(cstart,1), df_total_pairs, df_Minvhalf, nContsAux, 0.0d0, &
              dfC_full_compact(cstart,1), df_total_pairs)
   call openblas_set_num_threads(1)
enddo
df_full_lr_metric_built = .true.
end subroutine build_df_exchange_metric_full_lr

module subroutine integrals_build_exchange_df(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Ka(nConts_in,nConts_in), Kb(nConts_in,nConts_in)
logical :: restricted

call ensure_df_built()

restricted = (nOccB .eq. nOccA) .and. all(Ca(:,1:nOccA) .eq. Cb(:,1:nOccB))

if (df_direct_mode) then
   if (restricted) then
      call df_build_exchange_direct(nConts_in, Ca, nOccA, Cb, 0, Ka, Kb)
      Kb = Ka
   else
      call df_build_exchange_direct(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
   endif
   return
endif

if (.not. dfK_built) then
   call prof_start("df_k_metric_build")
   call build_df_exchange_metric()
   call prof_stop("df_k_metric_build")
   dfK_built = .true.
endif

Ka = 0.0d0
if (nOccA .gt. 0) then
   block
      integer,parameter :: DF_K_BMI_CHUNK = 2000
      real(8),allocatable :: Bmi_chunk(:,:,:)
      real(8),allocatable :: Ca_gathered(:,:)
      integer :: m, rs, rc, nn, slot, qstart, qend, csize
      do qstart = 1,nContsAux,DF_K_BMI_CHUNK
         qend = min(qstart+DF_K_BMI_CHUNK-1, nContsAux)
         csize = qend-qstart+1
         allocate(Bmi_chunk(nOccA,csize,nConts_in))
         !$omp parallel private(Ca_gathered,m,rs,rc,nn,slot)
         allocate(Ca_gathered(df_max_row_count,nOccA))
         !$omp do schedule(guided)
         do m = 1,nConts_in
            rs = df_row_start(m)
            rc = df_row_count(m)
            if (rc .gt. 0) then
               do nn = 1,rc
                  slot = rs+nn-1
                  Ca_gathered(nn,1:nOccA) = Ca(df_pair_col(slot), 1:nOccA)
               enddo
               call dgemm('T','N', nOccA, csize, rc, 1.0d0, &
                          Ca_gathered, df_max_row_count, dfK_compact(rs,qstart), df_total_pairs, 0.0d0, &
                          Bmi_chunk(1,1,m), nOccA)
            else
               Bmi_chunk(1:nOccA,1:csize,m) = 0.0d0
            endif
         enddo
         !$omp end do
         deallocate(Ca_gathered)
         !$omp end parallel
         call openblas_set_num_threads(omp_get_max_threads())
         call dgemm('T','N', nConts_in, nConts_in, nOccA*csize, 1.0d0, &
                    Bmi_chunk, nOccA*csize, Bmi_chunk, nOccA*csize, &
                    merge(0.0d0,1.0d0,qstart.eq.1), Ka, nConts_in)
         call openblas_set_num_threads(1)
         deallocate(Bmi_chunk)
      enddo
   end block
endif

if (restricted) then
   Kb = Ka
else
   Kb = 0.0d0
   if (nOccB .gt. 0) then
      block
         integer,parameter :: DF_K_BMI_CHUNK = 2000
         real(8),allocatable :: Bmi_chunk(:,:,:)
         real(8),allocatable :: Cb_gathered(:,:)
         integer :: m, rs, rc, nn, slot, qstart, qend, csize
         do qstart = 1,nContsAux,DF_K_BMI_CHUNK
            qend = min(qstart+DF_K_BMI_CHUNK-1, nContsAux)
            csize = qend-qstart+1
            allocate(Bmi_chunk(nOccB,csize,nConts_in))
            !$omp parallel private(Cb_gathered,m,rs,rc,nn,slot)
            allocate(Cb_gathered(df_max_row_count,nOccB))
            !$omp do schedule(guided)
            do m = 1,nConts_in
               rs = df_row_start(m)
               rc = df_row_count(m)
               if (rc .gt. 0) then
                  do nn = 1,rc
                     slot = rs+nn-1
                     Cb_gathered(nn,1:nOccB) = Cb(df_pair_col(slot), 1:nOccB)
                  enddo
                  call dgemm('T','N', nOccB, csize, rc, 1.0d0, &
                             Cb_gathered, df_max_row_count, dfK_compact(rs,qstart), df_total_pairs, 0.0d0, &
                             Bmi_chunk(1,1,m), nOccB)
               else
                  Bmi_chunk(1:nOccB,1:csize,m) = 0.0d0
               endif
            enddo
            !$omp end do
            deallocate(Cb_gathered)
            !$omp end parallel
            call openblas_set_num_threads(omp_get_max_threads())
            call dgemm('T','N', nConts_in, nConts_in, nOccB*csize, 1.0d0, &
                       Bmi_chunk, nOccB*csize, Bmi_chunk, nOccB*csize, &
                       merge(0.0d0,1.0d0,qstart.eq.1), Kb, nConts_in)
            call openblas_set_num_threads(1)
            deallocate(Bmi_chunk)
         enddo
      end block
   endif
endif
end subroutine integrals_build_exchange_df

module subroutine integrals_build_exchange_df_lr(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Ka(nConts_in,nConts_in), Kb(nConts_in,nConts_in)
logical :: restricted

call ensure_df_built()

restricted = (nOccB .eq. nOccA) .and. all(Ca(:,1:nOccA) .eq. Cb(:,1:nOccB))

if (df_direct_mode) then
   if (restricted) then
      call df_build_exchange_direct_lr(nConts_in, Ca, nOccA, Cb, 0, Ka, Kb)
      Kb = Ka
   else
      call df_build_exchange_direct_lr(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
   endif
   return
endif

if (.not. df_full_lr_metric_built) call build_df_exchange_metric_full_lr()

Ka = 0.0d0
if (nOccA .gt. 0) then
   block
      integer,parameter :: DF_K_BMI_CHUNK = 2000
      real(8),allocatable :: Bmi_full_chunk(:,:,:), Bmi_lr_chunk(:,:,:)
      real(8),allocatable :: Ka_raw(:,:)
      real(8),allocatable :: Ca_gathered(:,:)
      integer :: m, rs, rc, nn, slot, qstart, qend, csize
      allocate(Ka_raw(nConts_in,nConts_in))
      do qstart = 1,nContsAux,DF_K_BMI_CHUNK
         qend = min(qstart+DF_K_BMI_CHUNK-1, nContsAux)
         csize = qend-qstart+1
         allocate(Bmi_full_chunk(nOccA,csize,nConts_in))
         allocate(Bmi_lr_chunk(nOccA,csize,nConts_in))
         !$omp parallel private(Ca_gathered,m,rs,rc,nn,slot)
         allocate(Ca_gathered(df_max_row_count,nOccA))
         !$omp do schedule(guided)
         do m = 1,nConts_in
            rs = df_row_start(m)
            rc = df_row_count(m)
            if (rc .gt. 0) then
               do nn = 1,rc
                  slot = rs+nn-1
                  Ca_gathered(nn,1:nOccA) = Ca(df_pair_col(slot), 1:nOccA)
               enddo
               call dgemm('T','N', nOccA, csize, rc, 1.0d0, &
                          Ca_gathered, df_max_row_count, dfC_full_compact(rs,qstart), df_total_pairs, 0.0d0, &
                          Bmi_full_chunk(1,1,m), nOccA)
               call dgemm('T','N', nOccA, csize, rc, 1.0d0, &
                          Ca_gathered, df_max_row_count, dfB_compact_LR(rs,qstart), df_total_pairs, 0.0d0, &
                          Bmi_lr_chunk(1,1,m), nOccA)
            else
               Bmi_full_chunk(1:nOccA,1:csize,m) = 0.0d0
               Bmi_lr_chunk(1:nOccA,1:csize,m) = 0.0d0
            endif
         enddo
         !$omp end do
         deallocate(Ca_gathered)
         !$omp end parallel
         call openblas_set_num_threads(omp_get_max_threads())
         call dgemm('T','N', nConts_in, nConts_in, nOccA*csize, 1.0d0, &
                    Bmi_full_chunk, nOccA*csize, Bmi_lr_chunk, nOccA*csize, &
                    merge(0.0d0,1.0d0,qstart.eq.1), Ka_raw, nConts_in)
         call openblas_set_num_threads(1)
         deallocate(Bmi_full_chunk,Bmi_lr_chunk)
      enddo
      Ka = 0.5d0*(Ka_raw + transpose(Ka_raw))
      deallocate(Ka_raw)
   end block
endif

if (restricted) then
   Kb = Ka
else
   Kb = 0.0d0
   if (nOccB .gt. 0) then
      block
         integer,parameter :: DF_K_BMI_CHUNK = 2000
         real(8),allocatable :: Bmi_full_chunk(:,:,:), Bmi_lr_chunk(:,:,:)
         real(8),allocatable :: Kb_raw(:,:)
         real(8),allocatable :: Cb_gathered(:,:)
         integer :: m, rs, rc, nn, slot, qstart, qend, csize
         allocate(Kb_raw(nConts_in,nConts_in))
         do qstart = 1,nContsAux,DF_K_BMI_CHUNK
            qend = min(qstart+DF_K_BMI_CHUNK-1, nContsAux)
            csize = qend-qstart+1
            allocate(Bmi_full_chunk(nOccB,csize,nConts_in))
            allocate(Bmi_lr_chunk(nOccB,csize,nConts_in))
            !$omp parallel private(Cb_gathered,m,rs,rc,nn,slot)
            allocate(Cb_gathered(df_max_row_count,nOccB))
            !$omp do schedule(guided)
            do m = 1,nConts_in
               rs = df_row_start(m)
               rc = df_row_count(m)
               if (rc .gt. 0) then
                  do nn = 1,rc
                     slot = rs+nn-1
                     Cb_gathered(nn,1:nOccB) = Cb(df_pair_col(slot), 1:nOccB)
                  enddo
                  call dgemm('T','N', nOccB, csize, rc, 1.0d0, &
                             Cb_gathered, df_max_row_count, dfC_full_compact(rs,qstart), df_total_pairs, 0.0d0, &
                             Bmi_full_chunk(1,1,m), nOccB)
                  call dgemm('T','N', nOccB, csize, rc, 1.0d0, &
                             Cb_gathered, df_max_row_count, dfB_compact_LR(rs,qstart), df_total_pairs, 0.0d0, &
                             Bmi_lr_chunk(1,1,m), nOccB)
               else
                  Bmi_full_chunk(1:nOccB,1:csize,m) = 0.0d0
                  Bmi_lr_chunk(1:nOccB,1:csize,m) = 0.0d0
               endif
            enddo
            !$omp end do
            deallocate(Cb_gathered)
            !$omp end parallel
            call openblas_set_num_threads(omp_get_max_threads())
            call dgemm('T','N', nConts_in, nConts_in, nOccB*csize, 1.0d0, &
                       Bmi_full_chunk, nOccB*csize, Bmi_lr_chunk, nOccB*csize, &
                       merge(0.0d0,1.0d0,qstart.eq.1), Kb_raw, nConts_in)
            call openblas_set_num_threads(1)
            deallocate(Bmi_full_chunk,Bmi_lr_chunk)
         enddo
         Kb = 0.5d0*(Kb_raw + transpose(Kb_raw))
         deallocate(Kb_raw)
      end block
   endif
endif
end subroutine integrals_build_exchange_df_lr

module subroutine df_build_exchange_direct(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
use MOL_info
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Ka(nConts_in,nConts_in), Kb(nConts_in,nConts_in)
real(8),allocatable :: Braw_a(:,:,:), Braw_b(:,:,:)
real(8),allocatable :: Bmi_a(:,:,:), Bmi_b(:,:,:)
real(8),allocatable :: buf3c(:,:,:)
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP,ii
integer :: shls(4)
real(8) :: val, maxdiag
integer,allocatable :: ao_offset0(:), aux_offset0(:)

if (.not. df_Minvhalf_built) call build_df_metric_invhalf()

allocate(Braw_a(nConts_in,nOccA,nContsAux))
Braw_a = 0.0d0
allocate(Braw_b(nConts_in,nOccB,nContsAux))
Braw_b = 0.0d0

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
maxdiag = maxval(aux_shell_bound)

call prof_start("dfK_direct_halftransform")
!$omp parallel do private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf3c,p,q,r,val) &
!$omp&   schedule(dynamic)
do ksh = nBases,nBases+nBasesAux-1
   dk = cgto_engine(ksh, basDF)
   offP = aux_offset0(ksh)
   do ish = 0,nBases-1
      di = cgto_engine(ish, basDF)
      offi = ao_offset0(ish)
      do jsh = 0,ish
         dj = cgto_engine(jsh, basDF)
         offj = ao_offset0(jsh)
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
         shls(1) = ish
         shls(2) = jsh
         shls(3) = ksh
         allocate(buf3c(di,dj,dk))
         call threec2e_engine(buf3c, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf3c(p,q,r)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  if (nOccA .gt. 0) then
                     Braw_a(offi+p,:,offP+r) = Braw_a(offi+p,:,offP+r) + val*Ca(offj+q,1:nOccA)
                     if (ish .ne. jsh) then
                        Braw_a(offj+q,:,offP+r) = Braw_a(offj+q,:,offP+r) + val*Ca(offi+p,1:nOccA)
                     endif
                  endif
                  if (nOccB .gt. 0) then
                     Braw_b(offi+p,:,offP+r) = Braw_b(offi+p,:,offP+r) + val*Cb(offj+q,1:nOccB)
                     if (ish .ne. jsh) then
                        Braw_b(offj+q,:,offP+r) = Braw_b(offj+q,:,offP+r) + val*Cb(offi+p,1:nOccB)
                     endif
                  endif
               enddo
            enddo
         enddo
         deallocate(buf3c)
      enddo
   enddo
enddo
!$omp end parallel do
deallocate(ao_offset0, aux_offset0)
call prof_stop("dfK_direct_halftransform")

call prof_start("dfK_direct_metric")
call openblas_set_num_threads(omp_get_max_threads())
if (nOccA .gt. 0) then
   allocate(Bmi_a(nConts_in,nOccA,nContsAux))
   call dgemm('N','N', nConts_in*nOccA, nContsAux, nContsAux, 1.0d0, &
              Braw_a, nConts_in*nOccA, df_Minvhalf, nContsAux, 0.0d0, Bmi_a, nConts_in*nOccA)
   deallocate(Braw_a)
endif
if (nOccB .gt. 0) then
   allocate(Bmi_b(nConts_in,nOccB,nContsAux))
   call dgemm('N','N', nConts_in*nOccB, nContsAux, nContsAux, 1.0d0, &
              Braw_b, nConts_in*nOccB, df_Minvhalf, nContsAux, 0.0d0, Bmi_b, nConts_in*nOccB)
   deallocate(Braw_b)
endif
call openblas_set_num_threads(1)
call prof_stop("dfK_direct_metric")

call prof_start("dfK_direct_contract")
Ka = 0.0d0
if (nOccA .gt. 0) then
   block
      real(8),allocatable :: Ka_local(:,:)
      !$omp parallel private(ii,Ka_local)
      allocate(Ka_local(nConts_in,nConts_in))
      Ka_local = 0.0d0
      !$omp do
      do ii = 1,nContsAux
         call dgemm('N','T', nConts_in, nConts_in, nOccA, 1.0d0, &
                    Bmi_a(:,:,ii), nConts_in, Bmi_a(:,:,ii), nConts_in, 1.0d0, Ka_local, nConts_in)
      enddo
      !$omp end do
      !$omp critical
      Ka = Ka + Ka_local
      !$omp end critical
      deallocate(Ka_local)
      !$omp end parallel
   end block
   deallocate(Bmi_a)
endif
Kb = 0.0d0
if (nOccB .gt. 0) then
   block
      real(8),allocatable :: Kb_local(:,:)
      !$omp parallel private(ii,Kb_local)
      allocate(Kb_local(nConts_in,nConts_in))
      Kb_local = 0.0d0
      !$omp do
      do ii = 1,nContsAux
         call dgemm('N','T', nConts_in, nConts_in, nOccB, 1.0d0, &
                    Bmi_b(:,:,ii), nConts_in, Bmi_b(:,:,ii), nConts_in, 1.0d0, Kb_local, nConts_in)
      enddo
      !$omp end do
      !$omp critical
      Kb = Kb + Kb_local
      !$omp end critical
      deallocate(Kb_local)
      !$omp end parallel
   end block
   deallocate(Bmi_b)
endif
call prof_stop("dfK_direct_contract")
end subroutine df_build_exchange_direct

subroutine df_direct_halftransform_3c2e(omega, opt, nConts_in, nOccA, Ca, nOccB, Cb, Braw_a, Braw_b)
use MOL_info
implicit none
real(8),intent(in) :: omega
integer(8),intent(in) :: opt
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Braw_a(nConts_in,nOccA,nContsAux), Braw_b(nConts_in,nOccB,nContsAux)
real(8),allocatable :: buf3c(:,:,:)
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP
integer :: shls(4)
real(8) :: val, maxdiag
integer,allocatable :: ao_offset0(:), aux_offset0(:)

Braw_a = 0.0d0
Braw_b = 0.0d0

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
maxdiag = maxval(aux_shell_bound)

if (omega .gt. 0.0d0) envDF(9) = omega

!$omp parallel do private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf3c,p,q,r,val) &
!$omp&   schedule(dynamic)
do ksh = nBases,nBases+nBasesAux-1
   dk = cgto_engine(ksh, basDF)
   offP = aux_offset0(ksh)
   do ish = 0,nBases-1
      di = cgto_engine(ish, basDF)
      offi = ao_offset0(ish)
      do jsh = 0,ish
         dj = cgto_engine(jsh, basDF)
         offj = ao_offset0(jsh)
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
         shls(1) = ish
         shls(2) = jsh
         shls(3) = ksh
         allocate(buf3c(di,dj,dk))
         call threec2e_engine(buf3c, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf3c(p,q,r)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  if (nOccA .gt. 0) then
                     Braw_a(offi+p,:,offP+r) = Braw_a(offi+p,:,offP+r) + val*Ca(offj+q,1:nOccA)
                     if (ish .ne. jsh) then
                        Braw_a(offj+q,:,offP+r) = Braw_a(offj+q,:,offP+r) + val*Ca(offi+p,1:nOccA)
                     endif
                  endif
                  if (nOccB .gt. 0) then
                     Braw_b(offi+p,:,offP+r) = Braw_b(offi+p,:,offP+r) + val*Cb(offj+q,1:nOccB)
                     if (ish .ne. jsh) then
                        Braw_b(offj+q,:,offP+r) = Braw_b(offj+q,:,offP+r) + val*Cb(offi+p,1:nOccB)
                     endif
                  endif
               enddo
            enddo
         enddo
         deallocate(buf3c)
      enddo
   enddo
enddo
!$omp end parallel do
deallocate(ao_offset0, aux_offset0)

if (omega .gt. 0.0d0) envDF(9) = 0.0d0
end subroutine df_direct_halftransform_3c2e

module subroutine df_build_exchange_direct_lr(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
use MOL_info
use mod_exchange, only: RS_omega
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Ka(nConts_in,nConts_in), Kb(nConts_in,nConts_in)
real(8),allocatable :: Braw_a(:,:,:), Braw_b(:,:,:)
real(8),allocatable :: Cmi_full_a(:,:,:), Cmi_full_b(:,:,:)
integer :: ii

if (.not. df_Minvhalf_built) call build_df_metric_invhalf()
if (.not. df_Minv_full_built) call build_df_minv_full()
call ensure_df_lr_optimizer_built()

call prof_start("dfK_direct_lr_full_halftransform")
allocate(Braw_a(nConts_in,nOccA,nContsAux), Braw_b(nConts_in,nOccB,nContsAux))
call df_direct_halftransform_3c2e(0.0d0, int3c2e_opt, nConts_in, nOccA, Ca, nOccB, Cb, Braw_a, Braw_b)
call prof_stop("dfK_direct_lr_full_halftransform")

call prof_start("dfK_direct_lr_full_metric")
allocate(Cmi_full_a(nConts_in,nOccA,nContsAux), Cmi_full_b(nConts_in,nOccB,nContsAux))
call openblas_set_num_threads(omp_get_max_threads())
if (nOccA .gt. 0) then
   call dgemm('N','N', nConts_in*nOccA, nContsAux, nContsAux, 1.0d0, &
              Braw_a, nConts_in*nOccA, df_Minv_full, nContsAux, 0.0d0, Cmi_full_a, nConts_in*nOccA)
endif
if (nOccB .gt. 0) then
   call dgemm('N','N', nConts_in*nOccB, nContsAux, nContsAux, 1.0d0, &
              Braw_b, nConts_in*nOccB, df_Minv_full, nContsAux, 0.0d0, Cmi_full_b, nConts_in*nOccB)
endif
call openblas_set_num_threads(1)
deallocate(Braw_a, Braw_b)
call prof_stop("dfK_direct_lr_full_metric")

call prof_start("dfK_direct_lr_attenuated_halftransform")
allocate(Braw_a(nConts_in,nOccA,nContsAux), Braw_b(nConts_in,nOccB,nContsAux))
call df_direct_halftransform_3c2e(RS_omega, int3c2e_opt_lr, nConts_in, nOccA, Ca, nOccB, Cb, Braw_a, Braw_b)
call prof_stop("dfK_direct_lr_attenuated_halftransform")

call prof_start("dfK_direct_lr_contract")
Ka = 0.0d0
if (nOccA .gt. 0) then
   block
      real(8),allocatable :: Ka_local(:,:)
      !$omp parallel private(ii,Ka_local)
      allocate(Ka_local(nConts_in,nConts_in))
      Ka_local = 0.0d0
      !$omp do
      do ii = 1,nContsAux
         call dgemm('N','T', nConts_in, nConts_in, nOccA, 1.0d0, &
                    Cmi_full_a(:,:,ii), nConts_in, Braw_a(:,:,ii), nConts_in, 1.0d0, Ka_local, nConts_in)
      enddo
      !$omp end do
      !$omp critical
      Ka = Ka + Ka_local
      !$omp end critical
      deallocate(Ka_local)
      !$omp end parallel
   end block
   Ka = 0.5d0*(Ka + transpose(Ka))
endif
Kb = 0.0d0
if (nOccB .gt. 0) then
   block
      real(8),allocatable :: Kb_local(:,:)
      !$omp parallel private(ii,Kb_local)
      allocate(Kb_local(nConts_in,nConts_in))
      Kb_local = 0.0d0
      !$omp do
      do ii = 1,nContsAux
         call dgemm('N','T', nConts_in, nConts_in, nOccB, 1.0d0, &
                    Cmi_full_b(:,:,ii), nConts_in, Braw_b(:,:,ii), nConts_in, 1.0d0, Kb_local, nConts_in)
      enddo
      !$omp end do
      !$omp critical
      Kb = Kb + Kb_local
      !$omp end critical
      deallocate(Kb_local)
      !$omp end parallel
   end block
   Kb = 0.5d0*(Kb + transpose(Kb))
endif
deallocate(Braw_a, Braw_b, Cmi_full_a, Cmi_full_b)
call prof_stop("dfK_direct_lr_contract")
end subroutine df_build_exchange_direct_lr

module subroutine integrals_estimate_aux_size()
implicit none
if (len_trim(engine_df_aux_basis) .gt. 0) then
   call build_df_aux_basis_fromfile(trim(engine_df_aux_basis))
else
   call build_df_aux_basis()
endif
if (allocated(basDF))     deallocate(basDF)
if (allocated(envDF))     deallocate(envDF)
if (allocated(NorVECAux)) deallocate(NorVECAux)
end subroutine integrals_estimate_aux_size

end submodule df_exchange_impl
