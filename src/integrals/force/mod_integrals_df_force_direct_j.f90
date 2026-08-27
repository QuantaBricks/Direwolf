! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Density-fitting Coulomb-force gradient, direct (non-stored) 3-center integral mode.

submodule (mod_integrals) df_force_direct_j_impl
implicit none
contains

module subroutine integrals_force_df_direct_J(nConts_in, natoms_in, Ptot, force_J)
use MOL_info, only: nAtoms
use mod_profile, only: prof_start, prof_stop
implicit none
integer,intent(in) :: nConts_in, natoms_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: force_J(natoms_in,3)
real(8),allocatable :: bvec(:), tvec(:), cvec(:), W_J(:,:)
real(8),allocatable :: bvec_local(:)
real(8),allocatable :: dmax_shell(:,:)
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offP,offi,offj,t,A
integer :: shls(4)
real(8) :: val, maxdiag

call ensure_df_built()
allocate(cvec(nContsAux))

if (df_direct_cvec_cache_valid) then
   if (all(df_direct_cvec_cache_Ptot .eq. Ptot)) then
      cvec = df_direct_cvec_cache
      allocate(W_J(nContsAux,nContsAux))
      !$omp parallel do private(A)
      do A = 1,nContsAux
         W_J(:,A) = cvec(:)*cvec(A)
      enddo
      !$omp end parallel do
      force_J = 0.0d0
      call prof_start("df_force_direct_contract")
      call df_force_3c2e_contract_Jdirect(nConts_in, Ptot, cvec, natoms_in, force_J)
      block
         real(8) :: W_J3(nContsAux,nContsAux,1), force_J3(natoms_in,3,1)
         W_J3(:,:,1) = W_J
         force_J3(:,:,1) = force_J
         call df_force_2c2e_contract(1, W_J3, natoms_in, force_J3)
         force_J = force_J3(:,:,1)
      end block
      call prof_stop("df_force_direct_contract")
      deallocate(cvec, W_J)
      return
   endif
endif

allocate(bvec(nContsAux), tvec(nContsAux))
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
!$omp parallel do private(A)
do A = 1,nContsAux
   W_J(:,A) = cvec(:)*cvec(A)
enddo
!$omp end parallel do

force_J = 0.0d0
call prof_start("df_force_direct_contract")
call df_force_3c2e_contract_Jdirect(nConts_in, Ptot, cvec, natoms_in, force_J)
block
   real(8) :: W_J3(nContsAux,nContsAux,1), force_J3(natoms_in,3,1)
   W_J3(:,:,1) = W_J
   force_J3(:,:,1) = force_J
   call df_force_2c2e_contract(1, W_J3, natoms_in, force_J3)
   force_J = force_J3(:,:,1)
end block
call prof_stop("df_force_direct_contract")

deallocate(cvec, W_J)
end subroutine integrals_force_df_direct_J

subroutine df_force_3c2e_contract_Jdirect(nConts_in, Ptot, cvec, natoms_in, force_out)
use MOL_info, only: nAtoms
implicit none
integer,intent(in) :: nConts_in, natoms_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in), cvec(nContsAux)
real(8),intent(inout) :: force_out(natoms_in,3)

integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP
integer :: shls(4), atomI, atomJ, atomK
real(8),allocatable :: buf(:,:,:,:)
real(8) :: val(3), pair_weight, Tval
integer(8) :: ip1_opt, ip2_opt
integer,allocatable :: ao_offset0(:), aux_offset0(:)

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

!$omp parallel do private(ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP) &
!$omp&   private(shls,atomI,atomJ,atomK,buf,val,pair_weight,Tval) &
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
      do ksh = nBases,nBases+nBasesAux-1
         dk = cgto_engine(ksh, basDF)
         atomK = basDF(1,ksh+1) + 1
         offP = aux_offset0(ksh)

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
                  Tval = Ptot(offi+p,offj+q)*cvec(offP+r)
                  force_out(atomI,:) = force_out(atomI,:) + 2.0d0*Tval*val
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
                     Tval = Ptot(offj+q,offi+p)*cvec(offP+r)
                     force_out(atomJ,:) = force_out(atomJ,:) + pair_weight*Tval*val
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
                  Tval = Ptot(offi+p,offj+q)*cvec(offP+r)
                  force_out(atomK,:) = force_out(atomK,:) + pair_weight*Tval*val
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
end subroutine df_force_3c2e_contract_Jdirect

end submodule df_force_direct_j_impl
