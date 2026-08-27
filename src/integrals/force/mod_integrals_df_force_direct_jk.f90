! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Density-fitting Coulomb+exchange-force gradient, direct (non-stored) 3-center integral mode.

submodule (mod_integrals) df_force_direct_jk_impl
implicit none
contains

subroutine build_h_ij_direct(nConts_in, C, nOcc, h_ij, W_K)
use MOL_info
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in, nOcc
real(8),intent(in) :: C(nConts_in,nConts_in)
real(8),intent(out),allocatable :: h_ij(:,:,:), W_K(:,:)
real(8),allocatable :: Braw(:,:,:), g_ij(:,:,:), t1(:,:,:)
real(8),allocatable :: buf3c(:,:,:)
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP,A
integer :: shls(4)
real(8) :: val
integer,allocatable :: ao_offset0(:), aux_offset0(:)

allocate(Braw(nConts_in,nOcc,nContsAux))
Braw = 0.0d0
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

call prof_start("dfK_direct_force_halftransform")
!$omp parallel do private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,p,q,r,val,buf3c) &
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
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         allocate(buf3c(di,dj,dk))
         call threec2e_engine(buf3c, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf3c(p,q,r)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  Braw(offi+p,:,offP+r) = Braw(offi+p,:,offP+r) + val*C(offj+q,1:nOcc)
                  if (ish .ne. jsh) then
                     Braw(offj+q,:,offP+r) = Braw(offj+q,:,offP+r) + val*C(offi+p,1:nOcc)
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
call prof_stop("dfK_direct_force_halftransform")

call prof_start("dfK_direct_force_metric")
allocate(g_ij(nOcc,nOcc,nContsAux))
!$omp parallel do
do A = 1,nContsAux
   call dgemm('T','N', nOcc, nOcc, nConts_in, 1.0d0, &
              Braw(:,:,A), nConts_in, C(:,1:nOcc), nConts_in, 0.0d0, g_ij(:,:,A), nOcc)
enddo
!$omp end parallel do
deallocate(Braw)

allocate(t1(nOcc,nOcc,nContsAux))
call openblas_set_num_threads(omp_get_max_threads())
call dgemm('N','N', nOcc*nOcc, nContsAux, nContsAux, 1.0d0, &
           g_ij, nOcc*nOcc, df_evec, nContsAux, 0.0d0, t1, nOcc*nOcc)
call openblas_set_num_threads(1)
!$omp parallel do
do A = 1,nContsAux
   t1(:,:,A) = t1(:,:,A)*df_evalinv(A)
enddo
!$omp end parallel do
allocate(h_ij(nOcc,nOcc,nContsAux))
call openblas_set_num_threads(omp_get_max_threads())
call dgemm('N','T', nOcc*nOcc, nContsAux, nContsAux, 1.0d0, &
           t1, nOcc*nOcc, df_evec, nContsAux, 0.0d0, h_ij, nOcc*nOcc)
call openblas_set_num_threads(1)
deallocate(t1, g_ij)

allocate(W_K(nContsAux,nContsAux))
call openblas_set_num_threads(omp_get_max_threads())
call dgemm('T','N', nContsAux, nContsAux, nOcc*nOcc, 1.0d0, &
           h_ij, nOcc*nOcc, h_ij, nOcc*nOcc, 0.0d0, W_K, nContsAux)
call openblas_set_num_threads(1)
call prof_stop("dfK_direct_force_metric")
end subroutine build_h_ij_direct

module subroutine integrals_force_df_direct_JK(nConts_in, natoms_in, Pa, Pb, Ca, nOccA, Cb, nOccB, &
                                                force_J, force_Ka, force_Kb)
use MOL_info, only: nAtoms
use mod_profile, only: prof_start, prof_stop
implicit none
integer,intent(in) :: nConts_in, natoms_in, nOccA, nOccB
real(8),intent(in) :: Pa(nConts_in,nConts_in), Pb(nConts_in,nConts_in)
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: force_J(natoms_in,3), force_Ka(natoms_in,3), force_Kb(natoms_in,3)
real(8),allocatable :: Ptot(:,:), bvec(:), tvec(:), cvec(:), W_J(:,:)
real(8),allocatable :: bvec_local(:), dmax_shell(:,:)
real(8),allocatable :: h_ij_a(:,:,:), W_Ka(:,:), h_ij_b(:,:,:), W_Kb(:,:)
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offP,offi,offj,t
integer :: shls(4)
real(8) :: val, maxdiag
logical :: restricted_kb

call ensure_df_built()

allocate(Ptot(nConts_in,nConts_in))
Ptot = Pa + Pb
allocate(bvec(nContsAux), tvec(nContsAux), cvec(nContsAux))
bvec = 0.0d0
maxdiag = maxval(aux_shell_bound)
allocate(dmax_shell(0:nBases-1,0:nBases-1))
call compute_shell_density_bound(nConts_in, Ptot, dmax_shell)
call prof_start("df_force_direct_bpass")
!$omp parallel private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,p,q,r,val,t) &
!$omp&   private(bvec_local)
block
   real(8),allocatable :: buf3c(:,:,:)
   allocate(bvec_local(nContsAux))
   bvec_local = 0.0d0
   !$omp do schedule(dynamic)
   do t = 1,df_n_triples
      ish = df_triple_ish(t); jsh = df_triple_jsh(t); ksh = df_triple_ksh(t)
      di = df_shell_ncgto(ish); dj = df_shell_ncgto(jsh); dk = df_shell_ncgto(ksh)
      offi = df_shell_offset(ish); offj = df_shell_offset(jsh); offP = df_shell_offset(ksh)
      if (schwarz_bound(ish,jsh)*maxdiag*dmax_shell(ish,jsh) &
          .lt. density_screen_cutoff) cycle
      shls(1) = ish; shls(2) = jsh; shls(3) = ksh
      allocate(buf3c(di,dj,dk))
      call threec2e_engine(buf3c, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
      do r = 1,dk
         do q = 1,dj
            do p = 1,di
               val = buf3c(p,q,r)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
               if (ish .ne. jsh) then
                  bvec_local(offP+r) = bvec_local(offP+r) + val*(Ptot(offi+p,offj+q)+Ptot(offj+q,offi+p))
               else
                  bvec_local(offP+r) = bvec_local(offP+r) + val*Ptot(offi+p,offj+q)
               endif
            enddo
         enddo
      enddo
      deallocate(buf3c)
   enddo
   !$omp end do
   !$omp critical
   bvec = bvec + bvec_local
   !$omp end critical
   deallocate(bvec_local)
end block
!$omp end parallel
call prof_stop("df_force_direct_bpass")
call dgemv('T', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, bvec, 1, 0.0d0, tvec, 1)
tvec = tvec * df_evalinv
call dgemv('N', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, tvec, 1, 0.0d0, cvec, 1)
deallocate(bvec, tvec, dmax_shell)
allocate(W_J(nContsAux,nContsAux))
!$omp parallel do private(t)
do t = 1,nContsAux
   W_J(:,t) = cvec(:)*cvec(t)
enddo
!$omp end parallel do

call build_h_ij_direct(nConts_in, Ca, nOccA, h_ij_a, W_Ka)
restricted_kb = (nOccB .gt. 0) .and. (nOccB .eq. nOccA) .and. all(Ca(:,1:nOccA) .eq. Cb(:,1:nOccB))
if ((nOccB .gt. 0) .and. (.not. restricted_kb)) then
   call build_h_ij_direct(nConts_in, Cb, nOccB, h_ij_b, W_Kb)
endif

force_J = 0.0d0
force_Ka = 0.0d0
force_Kb = 0.0d0
call prof_start("df_force_direct_contract")
call df_force_3c2e_contract_JKdirect(nConts_in, Ptot, cvec, &
     Ca, nOccA, h_ij_a, Cb, nOccB, h_ij_b, restricted_kb, natoms_in, &
     force_J, force_Ka, force_Kb)
block
   real(8) :: W_J3(nContsAux,nContsAux,1), force_J3(natoms_in,3,1)
   real(8) :: W_Ka3(nContsAux,nContsAux,1), force_Ka3(natoms_in,3,1)
   W_J3(:,:,1) = W_J
   force_J3(:,:,1) = force_J
   call df_force_2c2e_contract(1, W_J3, natoms_in, force_J3)
   force_J = force_J3(:,:,1)
   W_Ka3(:,:,1) = W_Ka
   force_Ka3(:,:,1) = force_Ka
   call df_force_2c2e_contract(1, W_Ka3, natoms_in, force_Ka3)
   force_Ka = force_Ka3(:,:,1)
end block
deallocate(Ptot)
if (nOccB .gt. 0) then
   if (restricted_kb) then
      force_Kb = force_Ka
   else
      block
         real(8) :: W_Kb3(nContsAux,nContsAux,1), force_Kb3(natoms_in,3,1)
         W_Kb3(:,:,1) = W_Kb
         force_Kb3(:,:,1) = force_Kb
         call df_force_2c2e_contract(1, W_Kb3, natoms_in, force_Kb3)
         force_Kb = force_Kb3(:,:,1)
      end block
   endif
endif
call prof_stop("df_force_direct_contract")

deallocate(cvec, W_J, h_ij_a, W_Ka)
if (allocated(h_ij_b)) deallocate(h_ij_b)
if (allocated(W_Kb)) deallocate(W_Kb)
end subroutine integrals_force_df_direct_JK

subroutine df_force_3c2e_contract_JKdirect(nConts_in, Ptot, cvec, Ca, nOccA, h_ij_a, &
                                            Cb, nOccB, h_ij_b, restricted_kb, natoms_in, &
                                            force_J, force_Ka, force_Kb)
use MOL_info, only: nAtoms
implicit none
integer,intent(in) :: nConts_in, natoms_in, nOccA, nOccB
real(8),intent(in) :: Ptot(nConts_in,nConts_in), cvec(nContsAux)
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(in) :: h_ij_a(nOccA,nOccA,nContsAux)
real(8),intent(in) :: h_ij_b(:,:,:)
logical,intent(in) :: restricted_kb
real(8),intent(inout) :: force_J(natoms_in,3), force_Ka(natoms_in,3), force_Kb(natoms_in,3)

integer :: ish,jsh,ksh,p,q,r,i,j,di,dj,dk,offi,offj,offP
integer :: shls(4), atomI, atomJ, atomK
real(8),allocatable :: buf(:,:,:,:)
real(8),allocatable :: Hao_a_local(:,:,:), Hao_b_local(:,:,:), vtmp(:,:)
real(8) :: val(3), pair_weight, TvalJ, TvalKa, TvalKb
integer(8) :: ip1_opt, ip2_opt
integer,allocatable :: ao_offset0(:), aux_offset0(:)
logical :: need_kb_local
integer :: max_shell_dim

need_kb_local = (nOccB .gt. 0) .and. (.not. restricted_kb)

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

block
   integer :: sh
   max_shell_dim = 1
   do sh = 0,nBases-1
      max_shell_dim = max(max_shell_dim, cgto_engine(sh, basDF))
   enddo
end block

!$omp parallel private(ish,jsh,ksh,p,q,r,i,j,di,dj,dk,offi,offj,offP) &
!$omp&   private(shls,atomI,atomJ,atomK,buf,val,pair_weight,TvalJ,TvalKa,TvalKb) &
!$omp&   private(Hao_a_local,Hao_b_local,vtmp) &
!$omp&   reduction(+:force_J,force_Ka,force_Kb)
allocate(Hao_a_local(max_shell_dim,max_shell_dim,nContsAux))
allocate(vtmp(max(nOccA,nOccB),nContsAux))
if (need_kb_local) allocate(Hao_b_local(max_shell_dim,max_shell_dim,nContsAux))
!$omp do schedule(dynamic)
do ish = 0,nBases-1
   di = cgto_engine(ish, basDF)
   atomI = basDF(1,ish+1) + 1
   offi = ao_offset0(ish)
   do jsh = 0,ish
      dj = cgto_engine(jsh, basDF)
      atomJ = basDF(1,jsh+1) + 1
      offj = ao_offset0(jsh)
      pair_weight = merge(1.0d0, 2.0d0, ish .eq. jsh)

      do q = 1,dj
         do r = 1,nContsAux
            do i = 1,nOccA
               vtmp(i,r) = sum(h_ij_a(i,1:nOccA,r)*Ca(offj+q,1:nOccA))
            enddo
         enddo
         do r = 1,nContsAux
            do p = 1,di
               Hao_a_local(p,q,r) = sum(Ca(offi+p,1:nOccA)*vtmp(1:nOccA,r))
            enddo
         enddo
      enddo
      if (need_kb_local) then
         do q = 1,dj
            do r = 1,nContsAux
               do i = 1,nOccB
                  vtmp(i,r) = sum(h_ij_b(i,1:nOccB,r)*Cb(offj+q,1:nOccB))
               enddo
            enddo
            do r = 1,nContsAux
               do p = 1,di
                  Hao_b_local(p,q,r) = sum(Cb(offi+p,1:nOccB)*vtmp(1:nOccB,r))
               enddo
            enddo
         enddo
      endif

      do ksh = nBases,nBases+nBasesAux-1
         dk = cgto_engine(ksh, basDF)
         atomK = basDF(1,ksh+1) + 1
         offP = aux_offset0(ksh)

         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle

         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         allocate(buf(di,dj,dk,3))
         call threec2e_ip1_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip1_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf(p,q,r,:)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  TvalJ = Ptot(offi+p,offj+q)*cvec(offP+r)
                  force_J(atomI,:) = force_J(atomI,:) + 2.0d0*TvalJ*val
                  TvalKa = Hao_a_local(p,q,offP+r)
                  force_Ka(atomI,:) = force_Ka(atomI,:) + 2.0d0*TvalKa*val
                  if (need_kb_local) then
                     TvalKb = Hao_b_local(p,q,offP+r)
                     force_Kb(atomI,:) = force_Kb(atomI,:) + 2.0d0*TvalKb*val
                  endif
               enddo
            enddo
         enddo
         deallocate(buf)

         if (jsh .ne. ish) then
            shls(1) = jsh; shls(2) = ish; shls(3) = ksh
            allocate(buf(dj,di,dk,3))
            call threec2e_ip1_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip1_opt)
            do r = 1,dk
               do p = 1,di
                  do q = 1,dj
                     val = buf(q,p,r,:)*NorVEC(offj+q)*NorVEC(offi+p)*NorVECAux(offP+r)
                     TvalJ = Ptot(offj+q,offi+p)*cvec(offP+r)
                     force_J(atomJ,:) = force_J(atomJ,:) + pair_weight*TvalJ*val
                     TvalKa = Hao_a_local(p,q,offP+r)
                     force_Ka(atomJ,:) = force_Ka(atomJ,:) + pair_weight*TvalKa*val
                     if (need_kb_local) then
                        TvalKb = Hao_b_local(p,q,offP+r)
                        force_Kb(atomJ,:) = force_Kb(atomJ,:) + pair_weight*TvalKb*val
                     endif
                  enddo
               enddo
            enddo
            deallocate(buf)
         endif

         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         allocate(buf(di,dj,dk,3))
         call threec2e_ip2_engine(buf, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, ip2_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val = buf(p,q,r,:)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  TvalJ = Ptot(offi+p,offj+q)*cvec(offP+r)
                  force_J(atomK,:) = force_J(atomK,:) + pair_weight*TvalJ*val
                  TvalKa = Hao_a_local(p,q,offP+r)
                  force_Ka(atomK,:) = force_Ka(atomK,:) + pair_weight*TvalKa*val
                  if (need_kb_local) then
                     TvalKb = Hao_b_local(p,q,offP+r)
                     force_Kb(atomK,:) = force_Kb(atomK,:) + pair_weight*TvalKb*val
                  endif
               enddo
            enddo
         enddo
         deallocate(buf)
      enddo
   enddo
enddo
!$omp end do
deallocate(Hao_a_local, vtmp)
if (allocated(Hao_b_local)) deallocate(Hao_b_local)
!$omp end parallel
deallocate(ao_offset0, aux_offset0)
call cintdel_optimizer(ip1_opt)
call cintdel_optimizer(ip2_opt)
end subroutine df_force_3c2e_contract_JKdirect

end submodule df_force_direct_jk_impl
