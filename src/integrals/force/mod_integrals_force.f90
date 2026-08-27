! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Assembles the total analytic energy gradient (force) from 1-electron and 2-electron integral derivatives.

submodule (mod_integrals) force_impl
implicit none
contains

module subroutine integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
    real(8),intent(out) :: dS(nConts,nConts,3), dHcore(nConts,nConts,3)
    real(8),intent(out) :: dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3)

    integer :: i,j
    integer :: shls1e(2),shls2e(4)
    integer :: di,dj,dk,dl
    real(8),allocatable :: buf1e(:,:,:),buf1eV(:,:,:),buf1eT(:,:,:)
    integer(8) :: int2e_ip1_opt
    integer,allocatable :: ao_offset(:), shell_dim(:)
    integer :: si_off, num_off
    integer :: ish,jsh,ksh,lsh,pair_ij,pair_kl,a,b,c,d
    integer :: n_pairs_all,flat_ij
    integer(8) :: fstat_visit,fstat_cut_sw,fstat_cut_dn,fstat_call
    real(8) :: force_density_cutoff
    real(8),allocatable :: dmax_shell(:,:)
    logical :: all_distinct
    real(8),allocatable :: buf_i(:,:,:,:,:),buf_j(:,:,:,:,:),buf_k(:,:,:,:,:),buf_l(:,:,:,:,:)
    real(8),allocatable :: dS_local(:,:,:), dHcore_local(:,:,:)
    real(8),allocatable :: dJi_local(:,:,:), dKa_local(:,:,:), dKb_local(:,:,:)
    real(8),allocatable :: buf1eECPg(:,:,:)
    integer :: dims_arr(2)

call twoe_ip1_optimizer_engine(int2e_ip1_opt, atm, size(atm,2), bas, nBases, env)

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo

dS = 0
dHcore = 0
!$omp parallel private(i,j,shls1e,di,dj,buf1e,buf1eT,buf1eV,buf1eECPg,dims_arr) &
!$omp&   private(dS_local,dHcore_local)
allocate(dS_local(nConts,nConts,3), dHcore_local(nConts,nConts,3))
dS_local = 0
dHcore_local = 0
!$omp do schedule(static,1)   ! see mod_integrals_coulomb.f90's build_coulomb_direct: ordered merge alone is not enough
do i = 0,nBases-1
   do j = 0,nBases-1
      shls1e(1)=i
      shls1e(2)=j
      di = shell_dim(i)
      dj = shell_dim(j)
      allocate (buf1e(di,dj,3))
      allocate (buf1eT(di,dj,3))
      allocate (buf1eV(di,dj,3))
      call ovlp1e_ip_engine(buf1e, shls1e , atm, size(atm,2),bas,nBases,env,0_8)
      call store1edrv(shls1e,di,dj,ao_offset,nBases,buf1e,NorVEC,nConts,dS_local)
      call nuc1e_ip_engine(buf1eV, shls1e , atm_nuc, nAtoms_nuc,bas,nBases,env_nuc,0_8)
      call kin1e_ip_engine(buf1eT, shls1e , atm, size(atm,2),bas,nBases,env,0_8)
      if (necpbas > 0) then
         allocate(buf1eECPg(di,dj,3))
         dims_arr(1) = di
         dims_arr(2) = dj
         call ecp1e_ipnuc_engine(buf1eECPg, dims_arr, shls1e, atm, &
                                  size(atm,2), basECP, nBases+necpbas, envECP, 0_8, 0_8)
         call store1edrv(shls1e,di,dj,ao_offset,nBases,buf1eT+buf1eV+buf1eECPg,NorVEC,nConts,dHcore_local)
         deallocate(buf1eECPg)
      else
         call store1edrv(shls1e,di,dj,ao_offset,nBases,buf1eT+buf1eV,NorVEC,nConts,dHcore_local)
      endif
      deallocate (buf1e)
      deallocate (buf1eT)
      deallocate (buf1eV)
   enddo
enddo
!$omp end do
block
   integer :: merge_tid, merge_nthreads
   merge_nthreads = omp_get_num_threads()
   do merge_tid = 0, merge_nthreads-1
      !$omp barrier
      if (omp_get_thread_num() .eq. merge_tid) then
         dS = dS + dS_local
         dHcore = dHcore + dHcore_local
      endif
   enddo
end block
deallocate(dS_local, dHcore_local)
!$omp end parallel

dJi = 0
dKa = 0
dKb = 0
allocate(dmax_shell(0:nBases-1,0:nBases-1))
call compute_shell_density_bound(nConts, abs(Pa)+abs(Pb), dmax_shell)
fstat_visit = 0; fstat_cut_sw = 0; fstat_cut_dn = 0; fstat_call = 0

force_density_cutoff = 1.0d-9
block
   character(len=32) :: fdcenv
   fdcenv = ""
   call get_environment_variable("ENGINE_FORCE_SCREEN_CUTOFF", fdcenv)
   if (len_trim(fdcenv) .gt. 0) read(fdcenv,*) force_density_cutoff
end block

n_pairs_all = nBases*(nBases+1)/2
!$omp parallel private(flat_ij,ish,jsh,ksh,lsh,pair_ij,pair_kl,shls2e,di,dj,dk,dl) &
!$omp&   private(buf_i,buf_j,buf_k,buf_l,all_distinct,a,b,c,d) &
!$omp&   private(dJi_local,dKa_local,dKb_local) reduction(+:fstat_visit,fstat_cut_sw,fstat_cut_dn,fstat_call)
allocate(dJi_local(nConts,nConts,3), dKa_local(nConts,nConts,3), dKb_local(nConts,nConts,3))
dJi_local = 0
dKa_local = 0
dKb_local = 0
!$omp do schedule(static,1)   ! see mod_integrals_coulomb.f90's build_coulomb_direct: ordered merge alone is not enough
do flat_ij = 0,n_pairs_all-1
   pair_ij = flat_ij
   ish = int((sqrt(8.0d0*dble(pair_ij)+1.0d0)-1.0d0)/2.0d0)
   do while (ish*(ish+1)/2 .gt. pair_ij)
      ish = ish - 1
   enddo
   do while ((ish+1)*(ish+2)/2 .le. pair_ij)
      ish = ish + 1
   enddo
   jsh = pair_ij - ish*(ish+1)/2
   do ksh = 0,nBases-1
   do lsh = 0,ksh
      pair_kl = ksh*(ksh+1)/2 + lsh
      if (pair_kl .gt. pair_ij) cycle
      fstat_visit = fstat_visit + 1
      if (schwarz_bound(ish,jsh)*schwarz_bound(ksh,lsh) .lt. SCHWARZ_CUTOFF) then
         fstat_cut_sw = fstat_cut_sw + 1
         cycle
      endif
      if (schwarz_bound(ish,jsh)*schwarz_bound(ksh,lsh) &
          *max(dmax_shell(ish,jsh),dmax_shell(ksh,lsh)) .lt. force_density_cutoff) then
         fstat_cut_dn = fstat_cut_dn + 1
         cycle
      endif
      fstat_call = fstat_call + 1
      di = shell_dim(ish); dj = shell_dim(jsh); dk = shell_dim(ksh); dl = shell_dim(lsh)

      shls2e = (/ish,jsh,ksh,lsh/)
      allocate(buf_i(di,dj,dk,dl,3))
      call twoe_ip1_engine(buf_i, shls2e, atm, size(atm,2),bas,nBases,env,int2e_ip1_opt)
      call store2edrv(shls2e,ao_offset,buf_i,dJi_local,dKa_local,dKb_local,Pa,Pb,di,dj,dk,dl,nConts,nBases,NorVEC)

      if (jsh .ne. ish) then
         shls2e = (/jsh,ish,ksh,lsh/)
         allocate(buf_j(dj,di,dk,dl,3))
         call twoe_ip1_engine(buf_j, shls2e, atm, size(atm,2),bas,nBases,env,int2e_ip1_opt)
         call store2edrv(shls2e,ao_offset,buf_j,dJi_local,dKa_local,dKb_local,Pa,Pb,dj,di,dk,dl,nConts,nBases,NorVEC)
      endif

      if (pair_kl .eq. pair_ij) then
      else
         shls2e = (/ksh,lsh,ish,jsh/)
         allocate(buf_k(dk,dl,di,dj,3))
         call twoe_ip1_engine(buf_k, shls2e, atm, size(atm,2),bas,nBases,env,int2e_ip1_opt)
         call store2edrv(shls2e,ao_offset,buf_k,dJi_local,dKa_local,dKb_local,Pa,Pb,dk,dl,di,dj,nConts,nBases,NorVEC)

         if (lsh .eq. ksh) then
         else
            all_distinct = (ish.ne.jsh) .and. (ish.ne.ksh) .and. (ish.ne.lsh) .and. &
                            (jsh.ne.ksh) .and. (jsh.ne.lsh)
            if (all_distinct) then
               allocate(buf_l(dl,dk,di,dj,3))
               do d = 1,dj
                  do c = 1,di
                     do b = 1,dk
                        do a = 1,dl
                           buf_l(a,b,c,d,:) = -( buf_i(c,d,b,a,:) &
                                               + buf_j(d,c,b,a,:) &
                                               + buf_k(b,a,c,d,:) )
                        enddo
                     enddo
                  enddo
               enddo
               shls2e = (/lsh,ksh,ish,jsh/)
               call store2edrv(shls2e,ao_offset,buf_l,dJi_local,dKa_local,dKb_local,Pa,Pb,dl,dk,di,dj,nConts,nBases,NorVEC)
               deallocate(buf_l)
            else
               shls2e = (/lsh,ksh,ish,jsh/)
               allocate(buf_l(dl,dk,di,dj,3))
               call twoe_ip1_engine(buf_l, shls2e, atm, size(atm,2),bas,nBases,env,int2e_ip1_opt)
               call store2edrv(shls2e,ao_offset,buf_l,dJi_local,dKa_local,dKb_local,Pa,Pb,dl,dk,di,dj,nConts,nBases,NorVEC)
               deallocate(buf_l)
            endif
         endif
         deallocate(buf_k)
      endif
      if (jsh .ne. ish) deallocate(buf_j)
      deallocate(buf_i)
   enddo
enddo
enddo
!$omp end do
block
   integer :: merge_tid, merge_nthreads
   merge_nthreads = omp_get_num_threads()
   do merge_tid = 0, merge_nthreads-1
      !$omp barrier
      if (omp_get_thread_num() .eq. merge_tid) then
         dJi = dJi + dJi_local
         dKa = dKa + dKa_local
         dKb = dKb + dKb_local
      endif
   enddo
end block
deallocate(dJi_local, dKa_local, dKb_local)
!$omp end parallel
block
   character(len=8) :: fstat_env
   integer :: fstat_stat
   call get_environment_variable('ENGINE_JK_STATS', fstat_env, status=fstat_stat)
   if (fstat_stat .eq. 0 .and. trim(fstat_env) .eq. '1') then
      write(6,'(a,i12,a,i12,a,i12,a,i12,a,i12,a,f6.2,a)') &
         ' FSTATS pairs=',n_pairs_all,' visited=',fstat_visit, &
         ' cut_schwarz=',fstat_cut_sw,' cut_density=',fstat_cut_dn, &
         ' calls=',fstat_call,' (',100.0d0*real(fstat_call,8)/max(1.0d0,real(fstat_visit,8)),'% survive)'
      flush(6)
   endif
end block

call cintdel_optimizer(int2e_ip1_opt)
deallocate(ao_offset, shell_dim, dmax_shell)

end subroutine integrals_compute_force

module subroutine integrals_compute_force_exchange_lr(nConts, Pa, Pb, dKa_lr, dKb_lr)
    use mod_exchange, only: RS_omega
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
    real(8),intent(out) :: dKa_lr(nConts,nConts,3), dKb_lr(nConts,nConts,3)

    integer :: shls2e(4)
    integer :: di,dj,dk,dl
    integer(8) :: int2e_ip1_opt_lr
    integer,allocatable :: ao_offset(:), shell_dim(:)
    integer :: si_off, num_off
    integer :: ish,jsh,ksh,lsh,pair_ij,pair_kl,a,b,c,d
    integer :: n_pairs_all,flat_ij
    real(8),allocatable :: dmax_shell(:,:)
    logical :: all_distinct
    real(8),allocatable :: buf_i(:,:,:,:,:),buf_j(:,:,:,:,:),buf_k(:,:,:,:,:),buf_l(:,:,:,:,:)
    real(8),allocatable :: dJi_local(:,:,:), dKa_local(:,:,:), dKb_local(:,:,:)

    allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
    num_off = 0
    do si_off = 0,nBases-1
       ao_offset(si_off) = num_off
       shell_dim(si_off) = cgto_engine(si_off, bas)
       num_off = num_off + shell_dim(si_off)
    enddo

    env(9) = RS_omega
    call twoe_ip1_optimizer_engine(int2e_ip1_opt_lr, atm, size(atm,2), bas, nBases, env)

    dKa_lr = 0
    dKb_lr = 0
    allocate(dmax_shell(0:nBases-1,0:nBases-1))
    call compute_shell_density_bound(nConts, abs(Pa)+abs(Pb), dmax_shell)

    n_pairs_all = nBases*(nBases+1)/2
    !$omp parallel private(flat_ij,ish,jsh,ksh,lsh,pair_ij,pair_kl,shls2e,di,dj,dk,dl) &
    !$omp&   private(buf_i,buf_j,buf_k,buf_l,all_distinct,a,b,c,d) &
    !$omp&   private(dJi_local,dKa_local,dKb_local)
    allocate(dJi_local(nConts,nConts,3), dKa_local(nConts,nConts,3), dKb_local(nConts,nConts,3))
    dJi_local = 0
    dKa_local = 0
    dKb_local = 0
    !$omp do schedule(static,1)   ! see mod_integrals_coulomb.f90's build_coulomb_direct: ordered merge alone is not enough
    do flat_ij = 0,n_pairs_all-1
       pair_ij = flat_ij
       ish = int((sqrt(8.0d0*dble(pair_ij)+1.0d0)-1.0d0)/2.0d0)
       do while (ish*(ish+1)/2 .gt. pair_ij)
          ish = ish - 1
       enddo
       do while ((ish+1)*(ish+2)/2 .le. pair_ij)
          ish = ish + 1
       enddo
       jsh = pair_ij - ish*(ish+1)/2
       do ksh = 0,nBases-1
       do lsh = 0,ksh
          pair_kl = ksh*(ksh+1)/2 + lsh
          if (pair_kl .gt. pair_ij) cycle
          if (schwarz_bound(ish,jsh)*schwarz_bound(ksh,lsh) .lt. SCHWARZ_CUTOFF) cycle
          if (schwarz_bound(ish,jsh)*schwarz_bound(ksh,lsh) &
              *max(dmax_shell(ish,jsh),dmax_shell(ksh,lsh)) .lt. density_screen_cutoff) cycle
          di = shell_dim(ish); dj = shell_dim(jsh); dk = shell_dim(ksh); dl = shell_dim(lsh)

          shls2e = (/ish,jsh,ksh,lsh/)
          allocate(buf_i(di,dj,dk,dl,3))
          call twoe_ip1_engine(buf_i, shls2e, atm, size(atm,2),bas,nBases,env,int2e_ip1_opt_lr)
          call store2edrv(shls2e,ao_offset,buf_i,dJi_local,dKa_local,dKb_local,Pa,Pb,di,dj,dk,dl,nConts,nBases,NorVEC)

          if (jsh .ne. ish) then
             shls2e = (/jsh,ish,ksh,lsh/)
             allocate(buf_j(dj,di,dk,dl,3))
             call twoe_ip1_engine(buf_j, shls2e, atm, size(atm,2),bas,nBases,env,int2e_ip1_opt_lr)
             call store2edrv(shls2e,ao_offset,buf_j,dJi_local,dKa_local,dKb_local,Pa,Pb,dj,di,dk,dl,nConts,nBases,NorVEC)
          endif

          if (pair_kl .eq. pair_ij) then
          else
             shls2e = (/ksh,lsh,ish,jsh/)
             allocate(buf_k(dk,dl,di,dj,3))
             call twoe_ip1_engine(buf_k, shls2e, atm, size(atm,2),bas,nBases,env,int2e_ip1_opt_lr)
             call store2edrv(shls2e,ao_offset,buf_k,dJi_local,dKa_local,dKb_local,Pa,Pb,dk,dl,di,dj,nConts,nBases,NorVEC)

             if (lsh .eq. ksh) then
             else
                all_distinct = (ish.ne.jsh) .and. (ish.ne.ksh) .and. (ish.ne.lsh) .and. &
                                (jsh.ne.ksh) .and. (jsh.ne.lsh)
                if (all_distinct) then
                   allocate(buf_l(dl,dk,di,dj,3))
                   do d = 1,dj
                      do c = 1,di
                         do b = 1,dk
                            do a = 1,dl
                               buf_l(a,b,c,d,:) = -( buf_i(c,d,b,a,:) &
                                                   + buf_j(d,c,b,a,:) &
                                                   + buf_k(b,a,c,d,:) )
                            enddo
                         enddo
                      enddo
                   enddo
                   shls2e = (/lsh,ksh,ish,jsh/)
                   call store2edrv(shls2e,ao_offset,buf_l,dJi_local,dKa_local,dKb_local,Pa,Pb,dl,dk,di,dj,nConts,nBases,NorVEC)
                   deallocate(buf_l)
                else
                   shls2e = (/lsh,ksh,ish,jsh/)
                   allocate(buf_l(dl,dk,di,dj,3))
                   call twoe_ip1_engine(buf_l, shls2e, atm, size(atm,2),bas,nBases,env,int2e_ip1_opt_lr)
                   call store2edrv(shls2e,ao_offset,buf_l,dJi_local,dKa_local,dKb_local,Pa,Pb,dl,dk,di,dj,nConts,nBases,NorVEC)
                   deallocate(buf_l)
                endif
             endif
             deallocate(buf_k)
          endif
          if (jsh .ne. ish) deallocate(buf_j)
          deallocate(buf_i)
       enddo
    enddo
    enddo
    !$omp end do
    block
       integer :: merge_tid, merge_nthreads
       merge_nthreads = omp_get_num_threads()
       do merge_tid = 0, merge_nthreads-1
          !$omp barrier
          if (omp_get_thread_num() .eq. merge_tid) then
             dKa_lr = dKa_lr + dKa_local
             dKb_lr = dKb_lr + dKb_local
          endif
       enddo
    end block
    deallocate(dJi_local, dKa_local, dKb_local)
    !$omp end parallel

    call cintdel_optimizer(int2e_ip1_opt_lr)
    env(9) = 0.0d0
    deallocate(ao_offset, shell_dim, dmax_shell)

end subroutine integrals_compute_force_exchange_lr

module subroutine integrals_compute_force_1e(nConts, dS, dHcore)
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(out) :: dS(nConts,nConts,3), dHcore(nConts,nConts,3)

    integer :: i,j
    integer :: shls1e(2)
    integer :: di,dj
    real(8),allocatable :: buf1e(:,:,:),buf1eV(:,:,:),buf1eT(:,:,:)
    integer,allocatable :: ao_offset(:), shell_dim(:)
    integer :: si_off, num_off
    real(8),allocatable :: dS_local(:,:,:), dHcore_local(:,:,:)
    real(8),allocatable :: buf1eECPg(:,:,:)
    integer :: dims_arr(2)

    allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
    num_off = 0
    do si_off = 0,nBases-1
       ao_offset(si_off) = num_off
       shell_dim(si_off) = cgto_engine(si_off, bas)
       num_off = num_off + shell_dim(si_off)
    enddo

    dS = 0
    dHcore = 0
    !$omp parallel private(i,j,shls1e,di,dj,buf1e,buf1eT,buf1eV,buf1eECPg,dims_arr) &
    !$omp&   private(dS_local,dHcore_local)
    allocate(dS_local(nConts,nConts,3), dHcore_local(nConts,nConts,3))
    dS_local = 0
    dHcore_local = 0
    !$omp do schedule(static,1)   ! see mod_integrals_coulomb.f90's build_coulomb_direct: ordered merge alone is not enough
    do i = 0,nBases-1
       do j = 0,nBases-1
          shls1e(1)=i
          shls1e(2)=j
          di = shell_dim(i)
          dj = shell_dim(j)
          allocate (buf1e(di,dj,3))
          allocate (buf1eT(di,dj,3))
          allocate (buf1eV(di,dj,3))
          call ovlp1e_ip_engine(buf1e, shls1e , atm, size(atm,2),bas,nBases,env,0_8)
          call store1edrv(shls1e,di,dj,ao_offset,nBases,buf1e,NorVEC,nConts,dS_local)
          call nuc1e_ip_engine(buf1eV, shls1e , atm_nuc, nAtoms_nuc,bas,nBases,env_nuc,0_8)
          call kin1e_ip_engine(buf1eT, shls1e , atm, size(atm,2),bas,nBases,env,0_8)
          if (necpbas > 0) then
             allocate(buf1eECPg(di,dj,3))
             dims_arr(1) = di
             dims_arr(2) = dj
             call ecp1e_ipnuc_engine(buf1eECPg, dims_arr, shls1e, atm, &
                                      size(atm,2), basECP, nBases+necpbas, envECP, 0_8, 0_8)
             call store1edrv(shls1e,di,dj,ao_offset,nBases,buf1eT+buf1eV+buf1eECPg,NorVEC,nConts,dHcore_local)
             deallocate(buf1eECPg)
          else
             call store1edrv(shls1e,di,dj,ao_offset,nBases,buf1eT+buf1eV,NorVEC,nConts,dHcore_local)
          endif
          deallocate (buf1e)
          deallocate (buf1eT)
          deallocate (buf1eV)
       enddo
    enddo
    !$omp end do
    block
       integer :: merge_tid, merge_nthreads
       merge_nthreads = omp_get_num_threads()
       do merge_tid = 0, merge_nthreads-1
          !$omp barrier
          if (omp_get_thread_num() .eq. merge_tid) then
             dS = dS + dS_local
             dHcore = dHcore + dHcore_local
          endif
       enddo
    end block
    deallocate(dS_local, dHcore_local)
    !$omp end parallel

    deallocate(ao_offset, shell_dim)

end subroutine integrals_compute_force_1e

module subroutine integrals_nuc_attraction_deriv(coor_bohr, nConts, Dr)
    implicit none
    real(8),intent(in) :: coor_bohr(3)
    integer,intent(in) :: nConts
    real(8),intent(out) :: Dr(nConts,nConts,3)

    integer :: i2,j2,di,dj
    integer :: shls1e(2)
    real(8),allocatable :: buf1e(:,:,:)
    integer,allocatable :: ao_offset(:), shell_dim(:)
    integer :: si_off, num_off

    env(5) = coor_bohr(1)
    env(6) = coor_bohr(2)
    env(7) = coor_bohr(3)

    allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
    num_off = 0
    do si_off = 0,nBases-1
       ao_offset(si_off) = num_off
       shell_dim(si_off) = cgto_engine(si_off, bas)
       num_off = num_off + shell_dim(si_off)
    enddo

    Dr = 0.0d0

    !$omp parallel do private(i2,j2,shls1e,di,dj,buf1e) schedule(dynamic)
    do i2 = 0,nBases-1
      do j2 = 0,nBases-1
       if (schwarz_bound(i2,j2) .lt. SCHWARZ_CUTOFF) cycle
       shls1e(1)=i2
       shls1e(2)=j2
       di = shell_dim(i2)
       dj = shell_dim(j2)
       allocate (buf1e(di,dj,3))
       call rinv1e_ip_engine(buf1e, shls1e , atm, size(atm,2),bas,nBases,env,0_8)
       call store1edrv(shls1e,di,dj,ao_offset,nBases,buf1e,NorVEC,nConts,Dr)
       deallocate(buf1e)
      enddo
    enddo
    !$omp end parallel do
    deallocate(ao_offset, shell_dim)

end subroutine integrals_nuc_attraction_deriv

module subroutine integrals_nuc_deriv_allatoms(nConts, natoms_in, Ptot, force_out)
use MOL_info, only: atoms, nAtoms
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: nConts, natoms_in
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: force_out(natoms_in,3)

integer :: i2,j2,di,dj,i,j,k,a
integer :: shls1e(4)
integer :: si_off, num_off, off_grids, maxd
real(8) :: w
real(8),allocatable :: buf(:), env_g(:)
integer,allocatable :: ao_offset(:), shell_dim(:)

allocate(env_g(size(env) + 3*natoms_in))
env_g = 0.0d0
env_g(1:size(env)) = env(1:size(env))
off_grids = size(env)
do a = 1,natoms_in
   env_g(off_grids + 3*(a-1) + 1) = atoms(a)%coor(1)*ans2bohr
   env_g(off_grids + 3*(a-1) + 2) = atoms(a)%coor(2)*ans2bohr
   env_g(off_grids + 3*(a-1) + 3) = atoms(a)%coor(3)*ans2bohr
enddo
env_g(12) = dble(natoms_in)
env_g(13) = dble(off_grids)

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
maxd = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
   maxd = max(maxd, shell_dim(si_off))
enddo

force_out = 0.0d0

!$omp parallel private(i2,j2,di,dj,i,j,k,a,shls1e,w,buf) &
!$omp&   reduction(+:force_out)
allocate(buf(natoms_in*maxd*maxd*3))
!$omp do schedule(dynamic)
do i2 = 0,nBases-1
  do j2 = 0,nBases-1
   if (schwarz_bound(i2,j2) .lt. SCHWARZ_CUTOFF) cycle
   di = shell_dim(i2)
   dj = shell_dim(j2)
   shls1e(1) = i2
   shls1e(2) = j2
   shls1e(3) = 0
   shls1e(4) = natoms_in
   call grids1e_ip_engine(buf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8)
   do k = 1,3
      do j = 1,dj
         do i = 1,di
            w = Ptot(ao_offset(i2)+i, ao_offset(j2)+j) &
                * NorVEC(ao_offset(i2)+i) * NorVEC(ao_offset(j2)+j)
            if (w .eq. 0.0d0) cycle
            do a = 1,natoms_in
               force_out(a,k) = force_out(a,k) + w * &
                  buf(a + natoms_in*((i-1) + di*((j-1) + dj*(k-1))))
            enddo
         enddo
      enddo
   enddo
  enddo
enddo
!$omp end do
deallocate(buf)
!$omp end parallel

deallocate(ao_offset, shell_dim, env_g)

end subroutine integrals_nuc_deriv_allatoms

module subroutine integrals_ecp_force(nConts_in, natoms_in, Ptot, ecpForce)
use MOL_info
    implicit none
INCLUDE 'parameter.h'
    integer,intent(in) :: nConts_in, natoms_in
    real(8),intent(in) :: Ptot(nConts_in,nConts_in)
    real(8),intent(out) :: ecpForce(natoms_in,3)

    integer :: i,j,iatm,ia
    integer :: shls1e(2)
    integer :: di,dj
    integer :: dims_arr(2)
    real(8),allocatable :: buf1e(:,:,:)
    integer,allocatable :: ao_offset(:), shell_dim(:)
    integer :: si_off, num_off
    logical,allocatable :: atom_has_ecp(:)
    real(8),allocatable :: Dr_local(:,:,:), Dr(:,:,:)
    real(8),parameter :: ECP_SKIP_R2 = 144.0d0
    real(8) :: r_iatm(3), r_i(3), r_j(3)
    integer :: iatm_of_i, iatm_of_j

    ecpForce = 0.0d0
    if (necpbas .le. 0) return

    allocate(atom_has_ecp(natoms_in))
    atom_has_ecp = .false.
    do i = nBases+1, nBases+necpbas
       ia = basECP(1,i) + 1
       atom_has_ecp(ia) = .true.
    enddo

    allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
    num_off = 0
    do si_off = 0,nBases-1
       ao_offset(si_off) = num_off
       shell_dim(si_off) = cgto_engine(si_off, bas)
       num_off = num_off + shell_dim(si_off)
    enddo

    allocate(Dr(nConts_in,nConts_in,3))
    do iatm = 1,natoms_in
       if (.not. atom_has_ecp(iatm)) cycle
       envECP(18) = dble(iatm-1)
       r_iatm = atoms(iatm)%coor*ans2bohr
       Dr = 0.0d0
       !$omp parallel private(i,j,shls1e,di,dj,buf1e,dims_arr) private(Dr_local) &
       !$omp&   private(r_i,r_j,iatm_of_i,iatm_of_j)
       allocate(Dr_local(nConts_in,nConts_in,3))
       Dr_local = 0.0d0
       !$omp do schedule(dynamic)
       do i = 0,nBases-1
          do j = 0,nBases-1
             shls1e(1) = i
             shls1e(2) = j
             di = shell_dim(i)
             dj = shell_dim(j)
             allocate(buf1e(di,dj,3))
             iatm_of_i = bas(1,i+1) + 1
             iatm_of_j = bas(1,j+1) + 1
             r_i = atoms(iatm_of_i)%coor*ans2bohr
             r_j = atoms(iatm_of_j)%coor*ans2bohr
             if (sum((r_i-r_iatm)**2) .gt. ECP_SKIP_R2 .and. &
                 sum((r_j-r_iatm)**2) .gt. ECP_SKIP_R2) then
                buf1e = 0.0d0
             else
                dims_arr(1) = di
                dims_arr(2) = dj
                call ecp1e_iprinv_engine(buf1e, dims_arr, shls1e, atm, &
                                          size(atm,2), basECP, nBases+necpbas, envECP, 0_8, 0_8)
             endif
             call store1edrv(shls1e,di,dj,ao_offset,nBases,buf1e,NorVEC,nConts_in,Dr_local)
             deallocate(buf1e)
          enddo
       enddo
       !$omp end do
       !$omp critical
       Dr = Dr + Dr_local
       !$omp end critical
       deallocate(Dr_local)
       !$omp end parallel

       do i = 1,nConts_in
          do j = 1,nConts_in
             ecpForce(iatm,:) = ecpForce(iatm,:) + 2.0d0*Ptot(i,j)*Dr(i,j,:)
          enddo
       enddo
    enddo
    deallocate(Dr, ao_offset, shell_dim, atom_has_ecp)

end subroutine integrals_ecp_force

end submodule force_impl
