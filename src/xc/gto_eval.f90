! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Evaluates Gaussian basis functions (values and derivatives) on DFT quadrature grid points.

subroutine GTOeval(info)
use MOL_info
use GRID_INFO
use mod_meminfo, only: get_memory_budget_bytes
use mod_integrals, only: engine_puream
use mod_profile, only: fmt_gb
    implicit none
INCLUDE 'parameter.h'
    integer    :: i,j,k,info
    integer(8) :: need_bytes,avail_bytes,bytes_per_point,max_cache_bytes
    integer(8) :: cache_pts
    real(8),parameter :: SCREEN_EPS = 1.0d-7, SCREEN_MARGIN = 2.0d0
    real(8),parameter :: SIG_CUTOFF = 1.0d-7
    logical :: gto_progress_debug
    integer(8) :: gto_wc_rate,gto_wc_start,gto_wc_now
    character(len=8) :: dbgenv,cacheenv

    maxShell_shared = 0
    maxGauss_shared = 0
    do i = 1,natoms
       maxShell_shared = max(maxShell_shared,size(atoms(i)%shell))
       do j = 1,size(atoms(i)%shell)
          maxGauss_shared = max(maxGauss_shared,atoms(i)%shell(j)%nGauss)
       enddo
    enddo
    allocate(rcut2_shared(maxShell_shared,natoms))
    rcut2_shared = 0
    do i = 1,natoms
       do j = 1,size(atoms(i)%shell)
          rcut2_shared(j,i) = (sqrt(-log(SCREEN_EPS)/minval(atoms(i)%shell(j)%exponents)) + SCREEN_MARGIN)**2
       enddo
    enddo

    n_cache_batches = (size(Grids) + CACHE_BATCH - 1) / CACHE_BATCH
    allocate(batch_nsig(n_cache_batches))
    block
       integer,allocatable :: batch_sig_idx_tmp(:,:)
       integer :: ib2,bstart,bnb,ip2,l2,nsig2,ndim2,kk2,ilvl
       real(8),allocatable :: rcut2_lvl(:,:)
       type(BatchShellList),pointer :: dft_shl_p(:)
       real(8) :: coorr2(3),dist2_2
       real(8),allocatable :: rcut2_sig(:,:)
       logical :: shell_sig
       allocate(rcut2_sig(maxShell_shared,natoms))
       rcut2_sig = 0
       do i = 1,natoms
          do j = 1,size(atoms(i)%shell)
             rcut2_sig(j,i) = (sqrt(-log(SIG_CUTOFF)/minval(atoms(i)%shell(j)%exponents)) + SCREEN_MARGIN)**2
          enddo
       enddo
       allocate(batch_sig_idx_tmp(nconts,n_cache_batches))
       !$omp parallel do private(ib2,bstart,bnb,ip2,i,j,l2,ndim2,coorr2,dist2_2,shell_sig,nsig2,kk2) &
       !$omp&            schedule(dynamic)
       do ib2 = 1,n_cache_batches
          bstart = (ib2-1)*CACHE_BATCH+1
          bnb = min(CACHE_BATCH, size(Grids)-bstart+1)
          nsig2 = 0
          l2 = 1
          do i = 1,natoms
             do j = 1,size(atoms(i)%shell)
                if (engine_puream) then
                   ndim2 = 2*atoms(i)%shell(j)%angMoment+1
                else
                   ndim2 = (atoms(i)%shell(j)%angMoment+1)*(atoms(i)%shell(j)%angMoment+2)/2
                endif
                shell_sig = .false.
                do ip2 = 1,bnb
                   coorr2 = Grids(bstart+ip2-1)%coor - atoms(i)%coor*ans2bohr
                   dist2_2 = sum(coorr2**2)
                   if (dist2_2 .le. rcut2_sig(j,i)) then
                      shell_sig = .true.
                      exit
                   endif
                enddo
                if (shell_sig) then
                   do kk2 = 1,ndim2
                      batch_sig_idx_tmp(nsig2+kk2,ib2) = l2+kk2-1
                   enddo
                   nsig2 = nsig2+ndim2
                endif
                l2 = l2+ndim2
             enddo
          enddo
          batch_nsig(ib2) = nsig2
       enddo
       !$omp end parallel do
       max_batch_nsig = maxval(batch_nsig)
       allocate(batch_sig_idx(max_batch_nsig,n_cache_batches))
       !$omp parallel do private(ib2)
       do ib2 = 1,n_cache_batches
          batch_sig_idx(1:batch_nsig(ib2),ib2) = batch_sig_idx_tmp(1:batch_nsig(ib2),ib2)
       enddo
       !$omp end parallel do
       n_dft_batches = (size(Grids) + DFT_BATCH - 1) / DFT_BATCH
       allocate(dft_batch_shells(n_dft_batches))
       allocate(dft_batch_shells_lvl(n_dft_batches,XC_NLEVEL-1))
       allocate(rcut2_lvl(maxShell_shared,natoms))
       do ilvl = 1,XC_NLEVEL
       if (ilvl .eq. XC_NLEVEL) then
          dft_shl_p => dft_batch_shells
          rcut2_lvl = rcut2_sig
       else
          dft_shl_p => dft_batch_shells_lvl(:,ilvl)
          do i = 1,natoms
             do j = 1,size(atoms(i)%shell)
                rcut2_lvl(j,i) = (sqrt(-log(XC_LEVEL_CUT(ilvl))/minval(atoms(i)%shell(j)%exponents)) &
                                  + SCREEN_MARGIN)**2
             enddo
          enddo
       endif

       !$omp parallel do private(ib2,bstart,bnb,ip2,i,j,l2,ndim2,coorr2,dist2_2,shell_sig,nsig2,kk2) &
       !$omp&            schedule(dynamic)
       do ib2 = 1,n_dft_batches
          bstart = (ib2-1)*DFT_BATCH+1
          bnb = min(DFT_BATCH, size(Grids)-bstart+1)
          nsig2 = 0
          kk2 = 0
          do i = 1,natoms
             do j = 1,size(atoms(i)%shell)
                if (engine_puream) then
                   ndim2 = 2*atoms(i)%shell(j)%angMoment+1
                else
                   ndim2 = (atoms(i)%shell(j)%angMoment+1)*(atoms(i)%shell(j)%angMoment+2)/2
                endif
                shell_sig = .false.
                do ip2 = 1,bnb
                   coorr2 = Grids(bstart+ip2-1)%coor - atoms(i)%coor*ans2bohr
                   dist2_2 = sum(coorr2**2)
                   if (dist2_2 .le. rcut2_lvl(j,i)) then
                      shell_sig = .true.
                      exit
                   endif
                enddo
                if (shell_sig) then
                   kk2 = kk2 + 1
                   nsig2 = nsig2 + ndim2
                endif
             enddo
          enddo
          dft_shl_p(ib2)%n_shells = kk2
          dft_shl_p(ib2)%n_sig = nsig2
          allocate(dft_shl_p(ib2)%sh_atom(max(kk2,1)))
          allocate(dft_shl_p(ib2)%sh_idx(max(kk2,1)))
          allocate(dft_shl_p(ib2)%sh_off(max(kk2,1)))
          allocate(dft_shl_p(ib2)%sh_ndim(max(kk2,1)))
          allocate(dft_shl_p(ib2)%sig_idx(max(nsig2,1)))
          nsig2 = 0
          kk2 = 0
          l2 = 1
          do i = 1,natoms
             do j = 1,size(atoms(i)%shell)
                if (engine_puream) then
                   ndim2 = 2*atoms(i)%shell(j)%angMoment+1
                else
                   ndim2 = (atoms(i)%shell(j)%angMoment+1)*(atoms(i)%shell(j)%angMoment+2)/2
                endif
                shell_sig = .false.
                do ip2 = 1,bnb
                   coorr2 = Grids(bstart+ip2-1)%coor - atoms(i)%coor*ans2bohr
                   dist2_2 = sum(coorr2**2)
                   if (dist2_2 .le. rcut2_lvl(j,i)) then
                      shell_sig = .true.
                      exit
                   endif
                enddo
                if (shell_sig) then
                   kk2 = kk2 + 1
                   dft_shl_p(ib2)%sh_atom(kk2) = i
                   dft_shl_p(ib2)%sh_idx(kk2)  = j
                   dft_shl_p(ib2)%sh_off(kk2)  = nsig2
                   dft_shl_p(ib2)%sh_ndim(kk2) = ndim2
                   do ip2 = 1,ndim2
                      dft_shl_p(ib2)%sig_idx(nsig2+ip2) = l2+ip2-1
                   enddo
                   nsig2 = nsig2 + ndim2
                endif
                l2 = l2+ndim2
             enddo
          enddo
       enddo
       !$omp end parallel do
       enddo
       deallocate(rcut2_lvl)
       max_dft_batch_nsig = 0
       do ib2 = 1,n_dft_batches
          max_dft_batch_nsig = max(max_dft_batch_nsig, dft_batch_shells(ib2)%n_sig)
       enddo
       max_dft_batch_nsig_loose = 0
       do ilvl = 1,XC_NLEVEL-1
          do ib2 = 1,n_dft_batches
             max_dft_batch_nsig_loose = max(max_dft_batch_nsig_loose, dft_batch_shells_lvl(ib2,ilvl)%n_sig)
          enddo
       enddo
       dft_loose_built = .true.
       block
          integer(8) :: tt(XC_NLEVEL)
          integer :: lv
          do lv = 1,XC_NLEVEL
             tt(lv) = 0
             do ib2 = 1,n_dft_batches
                if (lv .eq. XC_NLEVEL) then
                   tt(lv) = tt(lv) + int(dft_batch_shells(ib2)%n_sig,8)
                else
                   tt(lv) = tt(lv) + int(dft_batch_shells_lvl(ib2,lv)%n_sig,8)
                endif
             enddo
          enddo
          if (engine_verbose .ge. 2) &
             print '(" XC geo levels (AOs/batch): ",3F8.1)', &
                   (real(tt(lv),8)/real(n_dft_batches,8), lv=1,XC_NLEVEL)
       end block
       dft_batch_shells_built = .true.
       deallocate(batch_sig_idx_tmp,rcut2_sig)
    end block
    if (engine_verbose .ge. 2) then
       print *,"Compact grid cache: max_batch_nsig=",max_batch_nsig," of nconts=",nconts, &
               " (",n_cache_batches," cache-batches of ",CACHE_BATCH," points)"
       print *,"DIRECT shell lists: max n_sig=",max_dft_batch_nsig," over ",n_dft_batches, &
               " batches of ",DFT_BATCH," points"
    endif
    call flush(6)

    grid_cache_mode = .not. xc_direct_mode

    if (.not. grid_cache_mode) then
       print *,"grid_cache_mode = OFF (DIRECT XC): DFT_calc will "// &
               "recompute basis-function values per batch every SCF "// &
               "iteration instead of reading a precomputed dense array - "// &
               "trades speed for O(1)-in-grid-size memory (Psi4-style)."
       call flush(6)
       return
    endif

    block
       integer :: ib4,bnb4
       integer(8) :: total_pf
       total_pf = 0_8
       do ib4 = 1,n_cache_batches
          bnb4 = min(CACHE_BATCH, size(Grids)-(ib4-1)*CACHE_BATCH)
          total_pf = total_pf + int(batch_nsig(ib4),8)*int(bnb4,8)
       enddo
       need_bytes = total_pf * 4_8 * 4_8
       bytes_per_point = need_bytes / max(int(size(Grids),8),1_8)
    end block
    avail_bytes = get_memory_budget_bytes()
    max_cache_bytes = int(avail_bytes,8)*8_8/10_8
    if (engine_verbose .ge. 2) then
       print '(A)', "  Grid val0/val1 storage would need="//trim(fmt_gb(real(need_bytes,8)/1024**3))// &
               "GB  available="//trim(fmt_gb(real(avail_bytes,8)/1024**3))//"GB"
       call flush(6)
    endif

    cache_pts = min(int(size(Grids),8), max_cache_bytes/bytes_per_point)
    cache_pts = (cache_pts/CACHE_BATCH)*CACHE_BATCH
    if (need_bytes .le. max_cache_bytes) then
       cache_pts = size(Grids)
    else if (cache_pts .le. 0) then
       print *,"Grid cache budget too small to cache even one DFT batch - "// &
               "falling back to fully on-the-fly (grid_cache_mode=OFF)."
       call flush(6)
       grid_cache_mode = .false.
       grid_cache_capacity = 0
       return
    else
       print *,"Grid cache budget only fits",cache_pts," of",size(Grids), &
               " grid points - caching that prefix, recomputing the rest "// &
               "on-the-fly per DFT batch."
       call flush(6)
    endif
    grid_cache_capacity = int(cache_pts)

  block
     integer :: ib4,npts_ib4,n_blocks_cached
     n_blocks_cached = (grid_cache_capacity + CACHE_BATCH - 1) / CACHE_BATCH
     allocate(val_blocks(n_blocks_cached))
     do ib4 = 1,n_blocks_cached
        npts_ib4 = min(CACHE_BATCH, grid_cache_capacity-(ib4-1)*CACHE_BATCH)
        allocate(val_blocks(ib4)%val0(batch_nsig(ib4),npts_ib4))
        allocate(val_blocks(ib4)%val1(batch_nsig(ib4),3,npts_ib4))
     enddo
  end block

  dbgenv = ""
  call get_environment_variable("ENGINE_GTO_PROGRESS", dbgenv)
  gto_progress_debug = (trim(dbgenv) .eq. "1")
  call system_clock(count_rate=gto_wc_rate)
  call system_clock(gto_wc_start)
  !$omp parallel private(k) firstprivate(gto_wc_start,gto_wc_rate)
  block
     real(8),allocatable :: val0_tmp(:),val1_tmp(:,:)
     integer :: ib3,local_col3
     allocate(val0_tmp(nconts),val1_tmp(nconts,3))
     !$omp do schedule(dynamic,1024)
     do k = 1,grid_cache_capacity
       if (gto_progress_debug .and. mod(k,100000).eq.0) then
          call system_clock(gto_wc_now)
          print '("GTOPROGRESS k=",I9," /",I9,"  elapsed=",F8.2," s")', &
                k, grid_cache_capacity, real(gto_wc_now-gto_wc_start,8)/gto_wc_rate
          call flush(6)
       endif
       ib3 = (k-1)/CACHE_BATCH + 1
       local_col3 = k - (ib3-1)*CACHE_BATCH
       call GTOeval_point(Grids(k)%coor,val0_tmp,val1_tmp)
       val_blocks(ib3)%val0(1:batch_nsig(ib3),local_col3) = val0_tmp(batch_sig_idx(1:batch_nsig(ib3),ib3))
       val_blocks(ib3)%val1(1:batch_nsig(ib3),:,local_col3) = val1_tmp(batch_sig_idx(1:batch_nsig(ib3),ib3),:)
     enddo
     !$omp end do
     deallocate(val0_tmp,val1_tmp)
  end block
  !$omp end parallel

end subroutine

subroutine GTOeval_point(coor_grid,val0,val1)
use MOL_info
use GRID_info, only: rcut2_shared,maxShell_shared,maxGauss_shared
use mod_integrals, only: engine_puream
implicit none
INCLUDE 'parameter.h'
real(8),intent(in) :: coor_grid(3)
real(8),intent(out) :: val0(nconts),val1(nconts,3)
integer :: i,j,l,a,b,c,ndim,ndim_out,idim,gi,nG,Lt
real(8) :: COORR(3),dist2
real(8) :: ex_shell(maxGauss_shared)
real(8) :: val,val1p(3),ex,cn,ge
real(8) :: cxa,cxa1,cxam1,cyb,cyb1,cybm1,czc,czc1,czcm1
real(8) :: px(0:6),py(0:6),pz(0:6),S0,S1
real(8) :: val0_cart(15),val1_cart(15,3)

l = 1
do i = 1,natoms
   COORR = coor_grid-atoms(i)%coor*ans2bohr
   dist2 = sum(COORR**2)
   do j = 1,size(atoms(i)%shell)
      Lt = atoms(i)%shell(j)%angMoment
      ndim = (Lt+1)*(Lt+2)/2
      ndim_out = merge(2*Lt+1, ndim, engine_puream)
      if (dist2 .gt. rcut2_shared(j,i)) then
         val0(l:l+ndim_out-1) = 0
         val1(l:l+ndim_out-1,:) = 0
         l = l + ndim_out
         cycle
      endif
      nG = atoms(i)%shell(j)%nGauss
      ex_shell(1:nG) = DEXP(-atoms(i)%shell(j)%exponents(1:nG)*dist2)
      px(0) = 1.0d0
      py(0) = 1.0d0
      pz(0) = 1.0d0
      do idim = 1,Lt+1
         px(idim) = px(idim-1)*COORR(1)
         py(idim) = py(idim-1)*COORR(2)
         pz(idim) = pz(idim-1)*COORR(3)
      enddo
      idim = 0
      do a =atoms(i)%shell(j)%angMoment,0,-1
         do b = atoms(i)%shell(j)%angMoment-a,0,-1
            c = atoms(i)%shell(j)%angMoment-a-b
            idim = idim+1

            cxa = px(a)
            cxa1 = px(a+1)
            cxam1 = px(max(a-1,0))
            cyb = py(b)
            cyb1 = py(b+1)
            cybm1 = py(max(b-1,0))
            czc = pz(c)
            czc1 = pz(c+1)
            czcm1 = pz(max(c-1,0))

            S0 = 0.0d0
            S1 = 0.0d0
            do gi = 1,nG
               ex = ex_shell(gi)*atoms(i)%shell(j)%cnVal(gi,idim)
               S0 = S0 + ex
               S1 = S1 + atoms(i)%shell(j)%exponents(gi)*ex
            enddo
            val = cxa*cyb*czc*S0
            val1p(1) = cyb*czc*(dble(a)*cxam1*S0 - 2.0d0*cxa1*S1)
            val1p(2) = cxa*czc*(dble(b)*cybm1*S0 - 2.0d0*cyb1*S1)
            val1p(3) = cxa*cyb*(dble(c)*czcm1*S0 - 2.0d0*czc1*S1)

              val0_cart(idim) = val
              val1_cart(idim,1) = val1p(1)
              val1_cart(idim,2) = val1p(2)
              val1_cart(idim,3) = val1p(3)
         enddo
      enddo
      if (Lt .ge. 2 .and. engine_puream) then
         if (Lt .eq. 2) then
            val0(l)   = val0_cart(2)
            val0(l+1) = val0_cart(5)
            val0(l+2) = val0_cart(6) - 0.5d0*val0_cart(1) - 0.5d0*val0_cart(4)
            val0(l+3) = val0_cart(3)
            val0(l+4) = (sqrt(3.0d0)/2.0d0)*(val0_cart(1) - val0_cart(4))
            val1(l,:)   = val1_cart(2,:)
            val1(l+1,:) = val1_cart(5,:)
            val1(l+2,:) = val1_cart(6,:) - 0.5d0*val1_cart(1,:) - 0.5d0*val1_cart(4,:)
            val1(l+3,:) = val1_cart(3,:)
            val1(l+4,:) = (sqrt(3.0d0)/2.0d0)*(val1_cart(1,:) - val1_cart(4,:))
         else if (Lt .eq. 3) then
            val0(l)   = (3.0d0*sqrt(2.0d0)/4.0d0)*val0_cart(2) - (sqrt(10.0d0)/4.0d0)*val0_cart(7)
            val0(l+1) = val0_cart(5)
            val0(l+2) = (sqrt(30.0d0)/5.0d0)*val0_cart(9) - (sqrt(30.0d0)/20.0d0)*val0_cart(2) &
                        - (sqrt(6.0d0)/4.0d0)*val0_cart(7)
            val0(l+3) = val0_cart(10) - (3.0d0*sqrt(5.0d0)/10.0d0)*val0_cart(3) &
                        - (3.0d0*sqrt(5.0d0)/10.0d0)*val0_cart(8)
            val0(l+4) = (sqrt(30.0d0)/5.0d0)*val0_cart(6) - (sqrt(6.0d0)/4.0d0)*val0_cart(1) &
                        - (sqrt(30.0d0)/20.0d0)*val0_cart(4)
            val0(l+5) = (sqrt(3.0d0)/2.0d0)*val0_cart(3) - (sqrt(3.0d0)/2.0d0)*val0_cart(8)
            val0(l+6) = (sqrt(10.0d0)/4.0d0)*val0_cart(1) - (3.0d0*sqrt(2.0d0)/4.0d0)*val0_cart(4)

            val1(l,:)   = (3.0d0*sqrt(2.0d0)/4.0d0)*val1_cart(2,:) - (sqrt(10.0d0)/4.0d0)*val1_cart(7,:)
            val1(l+1,:) = val1_cart(5,:)
            val1(l+2,:) = (sqrt(30.0d0)/5.0d0)*val1_cart(9,:) - (sqrt(30.0d0)/20.0d0)*val1_cart(2,:) &
                          - (sqrt(6.0d0)/4.0d0)*val1_cart(7,:)
            val1(l+3,:) = val1_cart(10,:) - (3.0d0*sqrt(5.0d0)/10.0d0)*val1_cart(3,:) &
                          - (3.0d0*sqrt(5.0d0)/10.0d0)*val1_cart(8,:)
            val1(l+4,:) = (sqrt(30.0d0)/5.0d0)*val1_cart(6,:) - (sqrt(6.0d0)/4.0d0)*val1_cart(1,:) &
                          - (sqrt(30.0d0)/20.0d0)*val1_cart(4,:)
            val1(l+5,:) = (sqrt(3.0d0)/2.0d0)*val1_cart(3,:) - (sqrt(3.0d0)/2.0d0)*val1_cart(8,:)
            val1(l+6,:) = (sqrt(10.0d0)/4.0d0)*val1_cart(1,:) - (3.0d0*sqrt(2.0d0)/4.0d0)*val1_cart(4,:)
         else if (Lt .eq. 4) then
            val0(l)   = (sqrt(5.0d0)/2.0d0)*val0_cart(2) - (sqrt(5.0d0)/2.0d0)*val0_cart(7)
            val0(l+1) = (3.0d0*sqrt(2.0d0)/4.0d0)*val0_cart(5) - (sqrt(10.0d0)/4.0d0)*val0_cart(12)
            val0(l+2) = (3.0d0*sqrt(7.0d0)/7.0d0)*val0_cart(9) - (sqrt(35.0d0)/14.0d0)*val0_cart(2) &
                        - (sqrt(35.0d0)/14.0d0)*val0_cart(7)
            val0(l+3) = (sqrt(70.0d0)/7.0d0)*val0_cart(14) - (3.0d0*sqrt(14.0d0)/28.0d0)*val0_cart(5) &
                        - (3.0d0*sqrt(70.0d0)/28.0d0)*val0_cart(12)
            val0(l+4) = val0_cart(15) - (3.0d0*sqrt(105.0d0)/35.0d0)*val0_cart(6) &
                        - (3.0d0*sqrt(105.0d0)/35.0d0)*val0_cart(13) + (3.0d0/8.0d0)*val0_cart(1) &
                        + (3.0d0*sqrt(105.0d0)/140.0d0)*val0_cart(4) + (3.0d0/8.0d0)*val0_cart(11)
            val0(l+5) = (sqrt(70.0d0)/7.0d0)*val0_cart(10) - (3.0d0*sqrt(70.0d0)/28.0d0)*val0_cart(3) &
                        - (3.0d0*sqrt(14.0d0)/28.0d0)*val0_cart(8)
            val0(l+6) = (3.0d0*sqrt(21.0d0)/14.0d0)*val0_cart(6) - (3.0d0*sqrt(21.0d0)/14.0d0)*val0_cart(13) &
                        - (sqrt(5.0d0)/4.0d0)*val0_cart(1) + (sqrt(5.0d0)/4.0d0)*val0_cart(11)
            val0(l+7) = (sqrt(10.0d0)/4.0d0)*val0_cart(3) - (3.0d0*sqrt(2.0d0)/4.0d0)*val0_cart(8)
            val0(l+8) = (sqrt(35.0d0)/8.0d0)*val0_cart(1) - (3.0d0*sqrt(3.0d0)/4.0d0)*val0_cart(4) &
                        + (sqrt(35.0d0)/8.0d0)*val0_cart(11)

            val1(l,:)   = (sqrt(5.0d0)/2.0d0)*val1_cart(2,:) - (sqrt(5.0d0)/2.0d0)*val1_cart(7,:)
            val1(l+1,:) = (3.0d0*sqrt(2.0d0)/4.0d0)*val1_cart(5,:) - (sqrt(10.0d0)/4.0d0)*val1_cart(12,:)
            val1(l+2,:) = (3.0d0*sqrt(7.0d0)/7.0d0)*val1_cart(9,:) - (sqrt(35.0d0)/14.0d0)*val1_cart(2,:) &
                          - (sqrt(35.0d0)/14.0d0)*val1_cart(7,:)
            val1(l+3,:) = (sqrt(70.0d0)/7.0d0)*val1_cart(14,:) - (3.0d0*sqrt(14.0d0)/28.0d0)*val1_cart(5,:) &
                          - (3.0d0*sqrt(70.0d0)/28.0d0)*val1_cart(12,:)
            val1(l+4,:) = val1_cart(15,:) - (3.0d0*sqrt(105.0d0)/35.0d0)*val1_cart(6,:) &
                          - (3.0d0*sqrt(105.0d0)/35.0d0)*val1_cart(13,:) + (3.0d0/8.0d0)*val1_cart(1,:) &
                          + (3.0d0*sqrt(105.0d0)/140.0d0)*val1_cart(4,:) + (3.0d0/8.0d0)*val1_cart(11,:)
            val1(l+5,:) = (sqrt(70.0d0)/7.0d0)*val1_cart(10,:) - (3.0d0*sqrt(70.0d0)/28.0d0)*val1_cart(3,:) &
                          - (3.0d0*sqrt(14.0d0)/28.0d0)*val1_cart(8,:)
            val1(l+6,:) = (3.0d0*sqrt(21.0d0)/14.0d0)*val1_cart(6,:) - (3.0d0*sqrt(21.0d0)/14.0d0)*val1_cart(13,:) &
                          - (sqrt(5.0d0)/4.0d0)*val1_cart(1,:) + (sqrt(5.0d0)/4.0d0)*val1_cart(11,:)
            val1(l+7,:) = (sqrt(10.0d0)/4.0d0)*val1_cart(3,:) - (3.0d0*sqrt(2.0d0)/4.0d0)*val1_cart(8,:)
            val1(l+8,:) = (sqrt(35.0d0)/8.0d0)*val1_cart(1,:) - (3.0d0*sqrt(3.0d0)/4.0d0)*val1_cart(4,:) &
                          + (sqrt(35.0d0)/8.0d0)*val1_cart(11,:)
         else
            print *, 'ERROR: puream grid evaluation for Lt=',Lt,' not implemented'
            stop 1
         endif
      else
         val0(l:l+ndim_out-1) = val0_cart(1:ndim)
         val1(l:l+ndim_out-1,:) = val1_cart(1:ndim,:)
      endif
      l = l + ndim_out
   enddo
enddo

end subroutine

subroutine GTOeval01(coor,a,b,c,nGauss,cn_vec,expo,ex_vec,val,val1)
implicit none
INCLUDE 'parameter.h'
real(8)      :: coor(3)
integer      :: a,b,c,nGauss
real(8)      :: cn_vec(nGauss),expo(nGauss),ex_vec(nGauss)
real(8)      :: val
real(8)      :: val1(3)
real(8)      :: ex,cn
integer      :: i

val = 0
val1 = 0
do i = 1,nGauss
   ex = ex_vec(i)
   cn = cn_vec(i)

   val = coor(1)**a * coor(2)**b * coor(3)**c * ex * cn + val

   if (a .eq. 0) then
      val1(1) = -2*expo(i)*coor(1)**(a+1) * coor(2)**b * coor(3)**c * ex * cn + val1(1)
   else
      val1(1) = ( a*coor(1)**(a-1) - 2*expo(i)*coor(1)**(a+1) ) * &
                coor(2)**b * coor(3)**c * ex * cn + val1(1)
   endif

   if (b .eq. 0) then
      val1(2) = -2*expo(i)*coor(2)**(b+1) * coor(1)**a * coor(3)**c * ex * cn + val1(2)
   else
      val1(2) = ( b*coor(2)**(b-1) - 2*expo(i)*coor(2)**(b+1) ) * &
                coor(1)**a * coor(3)**c * ex * cn + val1(2)
   endif

   if (c .eq. 0) then
      val1(3) = -2*expo(i)*coor(3)**(c+1) * coor(2)**b * coor(1)**a * ex * cn + val1(3)
   else
      val1(3) = ( c*coor(3)**(c-1) - 2*expo(i)*coor(3)**(c+1) ) * &
                coor(2)**b * coor(1)**a * ex * cn + val1(3)
   endif
enddo
end subroutine GTOeval01

subroutine GTOeval0(coor,a,b,c,nGauss,coeff,expo,val)
implicit none
INCLUDE 'parameter.h'
real(8)      :: coor(3)
integer      :: a,b,c,nGauss
real(8)      :: coeff(nGauss),expo(nGauss)
real(8)      :: Nor1,val,Nor2
real(8)      :: rs
integer      :: Lx1,Ly1,Lz1,Lx2,Ly2,Lz2,Lt
integer,external ::factorial,factorial2
integer      :: i,j,k,l

Lt = a+b+c
Lx1 = factorial2(a)
Ly1 = factorial2(b)
Lz1 = factorial2(c)
val = 0
do i = 1,nGauss
   Nor2 = (2**(2*Lt+1.5)*expo(i)**(Lt+1.5))**0.5 * (Lx1 * Ly1 * Lz1 * PI **(1.5))**(-0.5)
   rs = coor(1)**2 + coor(2)**2 + coor(3)**2
   val =coor(1)**a * coor(2)**b * coor(3)**c * &
         DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val
enddo
end subroutine

subroutine GTOeval1(coor,a,b,c,nGauss,coeff,expo,val)
implicit none
INCLUDE 'parameter.h'
real(8)      :: coor(3)
integer      :: a,b,c,nGauss
real(8)      :: coeff(nGauss),expo(nGauss)
real(8)      :: Nor1,Nor2
real(8)      :: val(3)
real(8)      :: rs
integer      :: Lx1,Ly1,Lz1,Lx2,Ly2,Lz2,Lt
integer,external ::factorial,factorial2
integer      :: i,j,k,l

Lt = a+b+c
Lx1 = factorial2(a)
Ly1 = factorial2(b)
Lz1 = factorial2(c)

val = 0
do i = 1,nGauss
   Nor2 = (2**(2*Lt+1.5)*expo(i)**(Lt+1.5))**0.5 * (Lx1 * Ly1 * Lz1 * PI **(1.5))**(-0.5)
   rs = coor(1)**2 + coor(2)**2 + coor(3)**2

   if (a .eq. 0) then
   val(1) = -2*expo(i)*coor(1)**(a+1) *&
            coor(2)**b * coor(3)**c *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(1)
   else
   val(1) = ( a*coor(1)**(a-1) - 2*expo(i)*coor(1)**(a+1) )*&
            coor(2)**b * coor(3)**c *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(1)
   endif

   if (b .eq. 0) then
   val(2) = -2*expo(i)*coor(2)**(b+1)*&
            coor(1)**a * coor(3)**c *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(2)

   else
   val(2) = ( b*coor(2)**(b-1) - 2*expo(i)*coor(2)**(b+1) )*&
            coor(1)**a * coor(3)**c *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(2)
   endif

   if (c .eq. 0) then
   val(3) = -2*expo(i)*coor(3)**(c+1)*&
            coor(2)**b * coor(1)**a *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(3)

   else
   val(3) = ( c*coor(3)**(c-1) - 2*expo(i)*coor(3)**(c+1) )*&
            coor(2)**b * coor(1)**a *DEXP(-expo(i)*rs) * coeff(i) * Nor2+ val(3)
   endif
enddo
end subroutine

subroutine GTOHess(coor,a,b,c,nGauss,cnVal,expo,val)
implicit none
INCLUDE 'parameter.h'
real(8)      :: coor(3)
integer      :: a,b,c,nGauss
real(8)      :: cnVal(nGauss),expo(nGauss)
real(8)      :: val(6)
real(8)      :: rs
real(8)      :: Tempa,Tempb,Tempc
real(8)      :: ex,cn
integer      :: Lx1,Ly1,Lz1,Lx2,Ly2,Lz2,Lt
integer,external ::factorial,factorial2
integer      :: i,j,k,l

Lt = a+b+c
val = 0

rs = coor(1)**2 + coor(2)**2 + coor(3)**2
do i = 1,nGauss
   ex = DEXP(-expo(i)*rs)
   cn = cnVal(i)

   if (a .lt. 2) then
   val(1) = (-2*expo(i)*(2*a+1)*coor(1)**a + 4*expo(i)**2 * coor(1)**(a+2)) *&
            coor(2)**b * coor(3)**c *ex * cn+ val(1)
   else
   val(1) = ( a*(a-1)*coor(1)**(a-2) &
             -2*expo(i)*(2*a+1)*coor(1)**a + 4*expo(i)**2 * coor(1)**(a+2) )*&
             coor(2)**b * coor(3)**c *ex * cn+ val(1)
   endif
   if (a .eq. 0) then
   Tempa = -2*expo(i)*coor(1)**(a+1)
   else
   Tempa =   a*coor(1)**(a-1) - 2*expo(i)*coor(1)**(a+1)
   endif

   if (b .eq. 0) then
   Tempb = -2*expo(i)*coor(2)**(b+1)
   else
   Tempb =   b*coor(2)**(b-1) - 2*expo(i)*coor(2)**(b+1)
   endif

   Tempc = coor(3)**c *ex * cn
   val(2) = Tempa*Tempb*Tempc+val(2)
   if (a .eq. 0) then
   Tempa = -2*expo(i)*coor(1)**(a+1)
   else
   Tempa =   a*coor(1)**(a-1) - 2*expo(i)*coor(1)**(a+1)
   endif

   if (c .eq. 0) then
   Tempb = -2*expo(i)*coor(3)**(c+1)
   else
   Tempb =   c*coor(3)**(c-1) - 2*expo(i)*coor(3)**(c+1)
   endif

   Tempc = coor(2)**b *ex * cn
   val(4) = Tempa*Tempb*Tempc+val(4)
   if (b .lt. 2) then
   val(3) = (-2*expo(i)*(2*b+1)*coor(2)**b + 4*expo(i)**2 * coor(2)**(b+2)) *&
            coor(1)**a * coor(3)**c *ex * cn+ val(3)
   else
   val(3) = ( b*(b-1)*coor(2)**(b-2) &
             -2*expo(i)*(2*b+1)*coor(2)**b + 4*expo(i)**2 * coor(2)**(b+2) )*&
             coor(1)**a * coor(3)**c *ex * cn+ val(3)
   endif
   if (b .eq. 0) then
   Tempa = -2*expo(i)*coor(2)**(b+1)
   else
   Tempa =   b*coor(2)**(b-1) - 2*expo(i)*coor(2)**(b+1)
   endif

   if (c .eq. 0) then
   Tempb = -2*expo(i)*coor(3)**(c+1)
   else
   Tempb =   c*coor(3)**(c-1) - 2*expo(i)*coor(3)**(c+1)
   endif

   Tempc = coor(1)**a *ex * cn
   val(5) = Tempa*Tempb*Tempc+val(5)

   if (c .lt. 2) then
   val(6) = (-2*expo(i)*(2*c+1)*coor(3)**c + 4*expo(i)**2 * coor(3)**(c+2)) *&
            coor(1)**a * coor(2)**b *ex * cn+ val(6)
   else
   val(6) = ( c*(c-1)*coor(3)**(c-2) &
             -2*expo(i)*(2*c+1)*coor(3)**c + 4*expo(i)**2 * coor(3)**(c+2) )*&
             coor(1)**a * coor(2)**b *ex * cn+ val(6)
   endif
enddo

end subroutine

subroutine GTOHess_batch(coor_batch,nb,a,b,c,nGauss,cnVal,expo,val_batch)
implicit none
integer,intent(in) :: nb,a,b,c,nGauss
real(8),intent(in) :: coor_batch(3,nb),cnVal(nGauss),expo(nGauss)
real(8),intent(out) :: val_batch(6,nb)
integer :: L
real(8) :: px(0:max(a,b,c)+2,nb),py(0:max(a,b,c)+2,nb),pz(0:max(a,b,c)+2,nb)
real(8) :: ex_all(nb,nGauss)

L = max(a,b,c)
call GTOHess_build_tables(coor_batch,nb,L,nGauss,expo,L+2,nb,px,py,pz,ex_all)
call GTOHess_batch_tab(nb,a,b,c,nGauss,cnVal,expo,L+2,nb,px,py,pz,ex_all,val_batch)

end subroutine

subroutine GTOHess_build_tables(coor_batch,nb,L,nGauss,expo,ldp,ldex,px,py,pz,ex_all)
implicit none
integer,intent(in) :: nb,L,nGauss,ldp,ldex
real(8),intent(in) :: coor_batch(3,nb),expo(nGauss)
real(8),intent(out) :: px(0:ldp,nb),py(0:ldp,nb),pz(0:ldp,nb)
real(8),intent(out) :: ex_all(ldex,nGauss)
real(8) :: rs(nb)
integer :: ip,n,i

do ip = 1,nb
   px(0,ip) = 1.0d0
   py(0,ip) = 1.0d0
   pz(0,ip) = 1.0d0
   do n = 1,L+2
      px(n,ip) = px(n-1,ip)*coor_batch(1,ip)
      py(n,ip) = py(n-1,ip)*coor_batch(2,ip)
      pz(n,ip) = pz(n-1,ip)*coor_batch(3,ip)
   enddo
   rs(ip) = coor_batch(1,ip)**2+coor_batch(2,ip)**2+coor_batch(3,ip)**2
enddo
do i = 1,nGauss
   do ip = 1,nb
      ex_all(ip,i) = DEXP(-expo(i)*rs(ip))
   enddo
enddo

end subroutine

subroutine GTOHess_batch_tab(nb,a,b,c,nGauss,cnVal,expo,ldp,ldex,px,py,pz,ex_all,val_batch)
implicit none
integer,intent(in) :: nb,a,b,c,nGauss,ldp,ldex
real(8),intent(in) :: cnVal(nGauss),expo(nGauss)
real(8),intent(in) :: px(0:ldp,nb),py(0:ldp,nb),pz(0:ldp,nb),ex_all(ldex,nGauss)
real(8),intent(out) :: val_batch(6,nb)
real(8) :: D1x,D1y,D1z,D2x,D2y,D2z,exc
real(8) :: ca2,cb2,cc2,ca1,cb1,cc1,e,e2,cn
integer :: i,ip,am2,bm2,cm2,am1,bm1,cm1

am2 = max(a-2,0) ; bm2 = max(b-2,0) ; cm2 = max(c-2,0)
am1 = max(a-1,0) ; bm1 = max(b-1,0) ; cm1 = max(c-1,0)
ca2 = dble(a*(a-1)) ; cb2 = dble(b*(b-1)) ; cc2 = dble(c*(c-1))
ca1 = dble(a)       ; cb1 = dble(b)       ; cc1 = dble(c)

val_batch = 0
do i = 1,nGauss
   e  = expo(i)
   e2 = e*e
   cn = cnVal(i)
   do ip = 1,nb
      exc = ex_all(ip,i)*cn
      D1x = ca1*px(am1,ip) - 2.0d0*e*px(a+1,ip)
      D1y = cb1*py(bm1,ip) - 2.0d0*e*py(b+1,ip)
      D1z = cc1*pz(cm1,ip) - 2.0d0*e*pz(c+1,ip)
      D2x = ca2*px(am2,ip) - 2.0d0*e*dble(2*a+1)*px(a,ip) + 4.0d0*e2*px(a+2,ip)
      D2y = cb2*py(bm2,ip) - 2.0d0*e*dble(2*b+1)*py(b,ip) + 4.0d0*e2*py(b+2,ip)
      D2z = cc2*pz(cm2,ip) - 2.0d0*e*dble(2*c+1)*pz(c,ip) + 4.0d0*e2*pz(c+2,ip)
      val_batch(1,ip) = val_batch(1,ip) + D2x*py(b,ip)*pz(c,ip)*exc
      val_batch(2,ip) = val_batch(2,ip) + D1x*D1y*pz(c,ip)*exc
      val_batch(3,ip) = val_batch(3,ip) + D2y*px(a,ip)*pz(c,ip)*exc
      val_batch(4,ip) = val_batch(4,ip) + D1x*D1z*py(b,ip)*exc
      val_batch(5,ip) = val_batch(5,ip) + D1y*D1z*px(a,ip)*exc
      val_batch(6,ip) = val_batch(6,ip) + D2z*px(a,ip)*py(b,ip)*exc
   enddo
enddo

end subroutine

subroutine GTOeval_shell_batch(i,j,coor_batch,nb,off,ld,val0_out,val1_out)
use MOL_info
use GRID_info, only: rcut2_shared,maxGauss_shared
use mod_integrals, only: engine_puream
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: i,j,nb,off,ld
real(8),intent(in) :: coor_batch(3,nb)
real(8),intent(inout) :: val0_out(ld,nb),val1_out(ld,3,nb)
integer :: ip,a,b,c,ndim,ndim_out,idim,gi,nG,Lt
real(8) :: COORR(3),dist2
real(8) :: ex_shell(maxGauss_shared)
real(8) :: val,val1p(3),ex,cn,ge
real(8) :: cxa,cxa1,cxam1,cyb,cyb1,cybm1,czc,czc1,czcm1
real(8) :: S0,S1,acen(3)
real(8) :: px(0:6),py(0:6),pz(0:6)
real(8) :: val0_cart(15),val1_cart(15,3)

Lt = atoms(i)%shell(j)%angMoment
ndim = (Lt+1)*(Lt+2)/2
ndim_out = merge(2*Lt+1, ndim, engine_puream)
nG = atoms(i)%shell(j)%nGauss

acen = atoms(i)%coor*ans2bohr

do ip = 1,nb
   COORR = coor_batch(:,ip)-acen
   dist2 = sum(COORR**2)
   if (dist2 .gt. rcut2_shared(j,i)) then
      val0_out(off+1:off+ndim_out,ip) = 0
      val1_out(off+1:off+ndim_out,:,ip) = 0
      cycle
   endif
      ex_shell(1:nG) = DEXP(-atoms(i)%shell(j)%exponents(1:nG)*dist2)
      px(0) = 1.0d0
      py(0) = 1.0d0
      pz(0) = 1.0d0
      do idim = 1,Lt+1
         px(idim) = px(idim-1)*COORR(1)
         py(idim) = py(idim-1)*COORR(2)
         pz(idim) = pz(idim-1)*COORR(3)
      enddo
      idim = 0
      do a = Lt,0,-1
         do b = Lt-a,0,-1
            c = Lt-a-b
            idim = idim+1

            cxa = px(a)
            cxa1 = px(a+1)
            cxam1 = px(max(a-1,0))
            cyb = py(b)
            cyb1 = py(b+1)
            cybm1 = py(max(b-1,0))
            czc = pz(c)
            czc1 = pz(c+1)
            czcm1 = pz(max(c-1,0))

            S0 = 0.0d0
            S1 = 0.0d0
            do gi = 1,nG
               ex = ex_shell(gi)*atoms(i)%shell(j)%cnVal(gi,idim)
               S0 = S0 + ex
               S1 = S1 + atoms(i)%shell(j)%exponents(gi)*ex
            enddo

            val0_cart(idim)   = cxa*cyb*czc*S0
            val1_cart(idim,1) = cyb*czc*(dble(a)*cxam1*S0 - 2.0d0*cxa1*S1)
            val1_cart(idim,2) = cxa*czc*(dble(b)*cybm1*S0 - 2.0d0*cyb1*S1)
            val1_cart(idim,3) = cxa*cyb*(dble(c)*czcm1*S0 - 2.0d0*czc1*S1)
         enddo
      enddo
      if (Lt .ge. 2 .and. engine_puream) then
         if (Lt .eq. 2) then
            val0_out(off+1,ip)   = val0_cart(2)
            val0_out(off+1+1,ip) = val0_cart(5)
            val0_out(off+1+2,ip) = val0_cart(6) - 0.5d0*val0_cart(1) - 0.5d0*val0_cart(4)
            val0_out(off+1+3,ip) = val0_cart(3)
            val0_out(off+1+4,ip) = (sqrt(3.0d0)/2.0d0)*(val0_cart(1) - val0_cart(4))
            val1_out(off+1,:,ip)   = val1_cart(2,:)
            val1_out(off+1+1,:,ip) = val1_cart(5,:)
            val1_out(off+1+2,:,ip) = val1_cart(6,:) - 0.5d0*val1_cart(1,:) - 0.5d0*val1_cart(4,:)
            val1_out(off+1+3,:,ip) = val1_cart(3,:)
            val1_out(off+1+4,:,ip) = (sqrt(3.0d0)/2.0d0)*(val1_cart(1,:) - val1_cart(4,:))
         else if (Lt .eq. 3) then
            val0_out(off+1,ip)   = (3.0d0*sqrt(2.0d0)/4.0d0)*val0_cart(2) - (sqrt(10.0d0)/4.0d0)*val0_cart(7)
            val0_out(off+1+1,ip) = val0_cart(5)
            val0_out(off+1+2,ip) = (sqrt(30.0d0)/5.0d0)*val0_cart(9) - (sqrt(30.0d0)/20.0d0)*val0_cart(2) &
                        - (sqrt(6.0d0)/4.0d0)*val0_cart(7)
            val0_out(off+1+3,ip) = val0_cart(10) - (3.0d0*sqrt(5.0d0)/10.0d0)*val0_cart(3) &
                        - (3.0d0*sqrt(5.0d0)/10.0d0)*val0_cart(8)
            val0_out(off+1+4,ip) = (sqrt(30.0d0)/5.0d0)*val0_cart(6) - (sqrt(6.0d0)/4.0d0)*val0_cart(1) &
                        - (sqrt(30.0d0)/20.0d0)*val0_cart(4)
            val0_out(off+1+5,ip) = (sqrt(3.0d0)/2.0d0)*val0_cart(3) - (sqrt(3.0d0)/2.0d0)*val0_cart(8)
            val0_out(off+1+6,ip) = (sqrt(10.0d0)/4.0d0)*val0_cart(1) - (3.0d0*sqrt(2.0d0)/4.0d0)*val0_cart(4)

            val1_out(off+1,:,ip)   = (3.0d0*sqrt(2.0d0)/4.0d0)*val1_cart(2,:) - (sqrt(10.0d0)/4.0d0)*val1_cart(7,:)
            val1_out(off+1+1,:,ip) = val1_cart(5,:)
            val1_out(off+1+2,:,ip) = (sqrt(30.0d0)/5.0d0)*val1_cart(9,:) - (sqrt(30.0d0)/20.0d0)*val1_cart(2,:) &
                          - (sqrt(6.0d0)/4.0d0)*val1_cart(7,:)
            val1_out(off+1+3,:,ip) = val1_cart(10,:) - (3.0d0*sqrt(5.0d0)/10.0d0)*val1_cart(3,:) &
                          - (3.0d0*sqrt(5.0d0)/10.0d0)*val1_cart(8,:)
            val1_out(off+1+4,:,ip) = (sqrt(30.0d0)/5.0d0)*val1_cart(6,:) - (sqrt(6.0d0)/4.0d0)*val1_cart(1,:) &
                          - (sqrt(30.0d0)/20.0d0)*val1_cart(4,:)
            val1_out(off+1+5,:,ip) = (sqrt(3.0d0)/2.0d0)*val1_cart(3,:) - (sqrt(3.0d0)/2.0d0)*val1_cart(8,:)
            val1_out(off+1+6,:,ip) = (sqrt(10.0d0)/4.0d0)*val1_cart(1,:) - (3.0d0*sqrt(2.0d0)/4.0d0)*val1_cart(4,:)
         else if (Lt .eq. 4) then
            val0_out(off+1,ip)   = (sqrt(5.0d0)/2.0d0)*val0_cart(2) - (sqrt(5.0d0)/2.0d0)*val0_cart(7)
            val0_out(off+1+1,ip) = (3.0d0*sqrt(2.0d0)/4.0d0)*val0_cart(5) - (sqrt(10.0d0)/4.0d0)*val0_cart(12)
            val0_out(off+1+2,ip) = (3.0d0*sqrt(7.0d0)/7.0d0)*val0_cart(9) - (sqrt(35.0d0)/14.0d0)*val0_cart(2) &
                        - (sqrt(35.0d0)/14.0d0)*val0_cart(7)
            val0_out(off+1+3,ip) = (sqrt(70.0d0)/7.0d0)*val0_cart(14) - (3.0d0*sqrt(14.0d0)/28.0d0)*val0_cart(5) &
                        - (3.0d0*sqrt(70.0d0)/28.0d0)*val0_cart(12)
            val0_out(off+1+4,ip) = val0_cart(15) - (3.0d0*sqrt(105.0d0)/35.0d0)*val0_cart(6) &
                        - (3.0d0*sqrt(105.0d0)/35.0d0)*val0_cart(13) + (3.0d0/8.0d0)*val0_cart(1) &
                        + (3.0d0*sqrt(105.0d0)/140.0d0)*val0_cart(4) + (3.0d0/8.0d0)*val0_cart(11)
            val0_out(off+1+5,ip) = (sqrt(70.0d0)/7.0d0)*val0_cart(10) - (3.0d0*sqrt(70.0d0)/28.0d0)*val0_cart(3) &
                        - (3.0d0*sqrt(14.0d0)/28.0d0)*val0_cart(8)
            val0_out(off+1+6,ip) = (3.0d0*sqrt(21.0d0)/14.0d0)*val0_cart(6) - (3.0d0*sqrt(21.0d0)/14.0d0)*val0_cart(13) &
                        - (sqrt(5.0d0)/4.0d0)*val0_cart(1) + (sqrt(5.0d0)/4.0d0)*val0_cart(11)
            val0_out(off+1+7,ip) = (sqrt(10.0d0)/4.0d0)*val0_cart(3) - (3.0d0*sqrt(2.0d0)/4.0d0)*val0_cart(8)
            val0_out(off+1+8,ip) = (sqrt(35.0d0)/8.0d0)*val0_cart(1) - (3.0d0*sqrt(3.0d0)/4.0d0)*val0_cart(4) &
                        + (sqrt(35.0d0)/8.0d0)*val0_cart(11)

            val1_out(off+1,:,ip)   = (sqrt(5.0d0)/2.0d0)*val1_cart(2,:) - (sqrt(5.0d0)/2.0d0)*val1_cart(7,:)
            val1_out(off+1+1,:,ip) = (3.0d0*sqrt(2.0d0)/4.0d0)*val1_cart(5,:) - (sqrt(10.0d0)/4.0d0)*val1_cart(12,:)
            val1_out(off+1+2,:,ip) = (3.0d0*sqrt(7.0d0)/7.0d0)*val1_cart(9,:) - (sqrt(35.0d0)/14.0d0)*val1_cart(2,:) &
                          - (sqrt(35.0d0)/14.0d0)*val1_cart(7,:)
            val1_out(off+1+3,:,ip) = (sqrt(70.0d0)/7.0d0)*val1_cart(14,:) - (3.0d0*sqrt(14.0d0)/28.0d0)*val1_cart(5,:) &
                          - (3.0d0*sqrt(70.0d0)/28.0d0)*val1_cart(12,:)
            val1_out(off+1+4,:,ip) = val1_cart(15,:) - (3.0d0*sqrt(105.0d0)/35.0d0)*val1_cart(6,:) &
                          - (3.0d0*sqrt(105.0d0)/35.0d0)*val1_cart(13,:) + (3.0d0/8.0d0)*val1_cart(1,:) &
                          + (3.0d0*sqrt(105.0d0)/140.0d0)*val1_cart(4,:) + (3.0d0/8.0d0)*val1_cart(11,:)
            val1_out(off+1+5,:,ip) = (sqrt(70.0d0)/7.0d0)*val1_cart(10,:) - (3.0d0*sqrt(70.0d0)/28.0d0)*val1_cart(3,:) &
                          - (3.0d0*sqrt(14.0d0)/28.0d0)*val1_cart(8,:)
            val1_out(off+1+6,:,ip) = (3.0d0*sqrt(21.0d0)/14.0d0)*val1_cart(6,:) - (3.0d0*sqrt(21.0d0)/14.0d0)*val1_cart(13,:) &
                          - (sqrt(5.0d0)/4.0d0)*val1_cart(1,:) + (sqrt(5.0d0)/4.0d0)*val1_cart(11,:)
            val1_out(off+1+7,:,ip) = (sqrt(10.0d0)/4.0d0)*val1_cart(3,:) - (3.0d0*sqrt(2.0d0)/4.0d0)*val1_cart(8,:)
            val1_out(off+1+8,:,ip) = (sqrt(35.0d0)/8.0d0)*val1_cart(1,:) - (3.0d0*sqrt(3.0d0)/4.0d0)*val1_cart(4,:) &
                          + (sqrt(35.0d0)/8.0d0)*val1_cart(11,:)
         else
            print *, 'ERROR: puream grid evaluation for Lt=',Lt,' not implemented'
            stop 1
         endif
      else
         val0_out(off+1:off+ndim_out,ip) = val0_cart(1:ndim)
         val1_out(off+1:off+ndim_out,:,ip) = val1_cart(1:ndim,:)
      endif
enddo

end subroutine

subroutine GTOeval_batch_compact(ibatch,batch_start,nb,ld,val0_sig,val1_sig)
use GRID_info, only: Grids,dft_batch_shells,dft_batch_shells_lvl,BatchShellList, &
                     xc_geo_level,XC_NLEVEL,dft_loose_built
implicit none
integer,intent(in) :: ibatch,batch_start,nb,ld
real(8),intent(inout) :: val0_sig(ld,nb),val1_sig(ld,3,nb)
integer :: k
real(8) :: coor_batch(3,nb)
type(BatchShellList),pointer :: sl

if (dft_loose_built .and. xc_geo_level .lt. XC_NLEVEL) then
   sl => dft_batch_shells_lvl(ibatch,xc_geo_level)
else
   sl => dft_batch_shells(ibatch)
endif

do k = 1,nb
   coor_batch(:,k) = Grids(batch_start+k-1)%coor
enddo
do k = 1,sl%n_shells
   call GTOeval_shell_batch(sl%sh_atom(k),sl%sh_idx(k), &
                            coor_batch,nb,sl%sh_off(k), &
                            ld,val0_sig,val1_sig)
enddo

end subroutine

subroutine GTOeval_build_shell_tables()
use MOL_info
use GRID_info
use mod_integrals, only: engine_puream
implicit none
integer :: i,j,k,ib,ish,l,ndim
integer,allocatable :: cnt(:),fill(:)

if (shell_tables_built) return

ish = 0
do i = 1,natoms
   ish = ish + size(atoms(i)%shell)
enddo
vxc_nshell_tot = ish
allocate(vxc_sh_atom(ish),vxc_sh_local(ish),vxc_sh_ao0(ish),vxc_sh_ndim(ish))
ish = 0
l = 1
do i = 1,natoms
   do j = 1,size(atoms(i)%shell)
      ish = ish + 1
      ndim = merge(2*atoms(i)%shell(j)%angMoment+1, &
                   (atoms(i)%shell(j)%angMoment+1)*(atoms(i)%shell(j)%angMoment+2)/2, &
                   engine_puream)
      vxc_sh_atom(ish)  = i
      vxc_sh_local(ish) = j
      vxc_sh_ao0(ish)   = l
      vxc_sh_ndim(ish)  = ndim
      l = l + ndim
   enddo
enddo

allocate(batch_natom(n_dft_batches))
allocate(batch_atoms(natoms,n_dft_batches))
do ib = 1,n_dft_batches
   batch_natom(ib) = 0
   do k = 1,dft_batch_shells(ib)%n_shells
      if (batch_natom(ib) .eq. 0) then
         batch_natom(ib) = 1
         batch_atoms(1,ib) = dft_batch_shells(ib)%sh_atom(k)
      else if (dft_batch_shells(ib)%sh_atom(k) .ne. batch_atoms(batch_natom(ib),ib)) then
         batch_natom(ib) = batch_natom(ib) + 1
         batch_atoms(batch_natom(ib),ib) = dft_batch_shells(ib)%sh_atom(k)
      endif
   enddo
enddo
allocate(batch_phimax(n_dft_batches))
allocate(ao_atom(nconts))
do ish = 1,vxc_nshell_tot
   ao_atom(vxc_sh_ao0(ish):vxc_sh_ao0(ish)+vxc_sh_ndim(ish)-1) = vxc_sh_atom(ish)
enddo
allocate(cnt(vxc_nshell_tot),fill(vxc_nshell_tot))
allocate(sh_in_batch(vxc_nshell_tot,n_dft_batches))
sh_in_batch = 0
cnt = 0
do ib = 1,n_dft_batches
   do k = 1,dft_batch_shells(ib)%n_shells
      ish = 0
      do i = 1,dft_batch_shells(ib)%sh_atom(k)-1
         ish = ish + size(atoms(i)%shell)
      enddo
      ish = ish + dft_batch_shells(ib)%sh_idx(k)
      cnt(ish) = cnt(ish) + 1
      sh_in_batch(ish,ib) = dft_batch_shells(ib)%sh_off(k) + 1
   enddo
enddo
allocate(shell_batches(vxc_nshell_tot))
do ish = 1,vxc_nshell_tot
   shell_batches(ish)%n = cnt(ish)
   allocate(shell_batches(ish)%batch(max(cnt(ish),1)))
   allocate(shell_batches(ish)%off(max(cnt(ish),1)))
enddo
fill = 0
do ib = 1,n_dft_batches
   do k = 1,dft_batch_shells(ib)%n_shells
      ish = 0
      do i = 1,dft_batch_shells(ib)%sh_atom(k)-1
         ish = ish + size(atoms(i)%shell)
      enddo
      ish = ish + dft_batch_shells(ib)%sh_idx(k)
      fill(ish) = fill(ish) + 1
      shell_batches(ish)%batch(fill(ish)) = ib
      shell_batches(ish)%off(fill(ish))   = dft_batch_shells(ib)%sh_off(k)
   enddo
enddo
deallocate(cnt,fill)
shell_tables_built = .true.

end subroutine
