! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX analytic COSX exchange gradient builder - split out of mod_integrals_exchange_cosx.f90

submodule (mod_integrals:exchange_cosx_impl) exchange_cosx_force_impl
contains
module subroutine integrals_force_exchange_cosx(nConts, natoms, Da, Db, RSomega, forceKa, forceKb)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
real(8),intent(in) :: RSomega
real(8),intent(out) :: forceKa(natoms,3), forceKb(natoms,3)

logical,save :: cosx_grad_wt = .false.
logical,save :: cosx_grad_wt_checked = .false.
character(len=8) :: cosx_grad_wt_env
real(8),allocatable :: xed_a(:), xed_b(:), wt_a(:,:), wt_b(:,:), eps_zero(:)
integer :: gp
logical,save :: cosx_grad_resp = .false.
logical,save :: cosx_grad_resp_checked = .false.
real(8),save :: fgr_ao_own = 1.0d0, fgr_ao_grid = 2.0d0
real(8),save :: fgr_esp_own = 0.0d0, fgr_esp_grid = 2.0d0
real(8),save :: fgr_wt = 1.0d0
character(len=16) :: gr_env
real(8),allocatable :: dGa(:,:,:), dGb(:,:,:), dbuf(:)
real(8) :: t_ao, t_esp, wij_d(3)
integer :: gpar, kdir

integer :: ib, nblocks, bstart, nb, g
integer,pointer :: blkstart_p(:)
integer :: si, sj, di, dj, i, j, ao_i, ao_j, maxd, sj_idx
integer :: shls1e(4), off_grids
real(8) :: nij, wij, x_block_max, k_bound
real(8),allocatable :: coor_batch(:,:), w_batch(:)
real(8),allocatable :: Xb(:,:), Xb1(:,:,:)
real(8),allocatable :: Fa(:,:), Fb(:,:), Ga(:,:), Gb(:,:)
real(8),allocatable :: Ha(:,:), Hb(:,:)
real(8),allocatable :: env_g(:), buf(:)
real(8),allocatable :: int_cache(:)
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),allocatable :: F_gmax(:)
integer :: npts, kappa, dloc, atom_idx
real(8) :: sKa, sKb
real(8),allocatable :: forceKa_local(:,:), forceKb_local(:,:)
logical :: use_cache
real(8),pointer :: G_hi_a_p(:,:), G_hi_b_p(:,:)
real(8),pointer :: esp_bound_p(:,:)
integer,pointer :: ext_ptr_p(:), ext_idx_p(:)
real(8) :: kscreen_use
integer :: ni_b, nj_b, offi_b, offj_b, pi_b, pj_b
real(8) :: r2_b, ei_b, ci_b, ej_b, cj_b, p_b

if (.not. cosx_ready .or. cosx_nBases_cached .ne. nBases) then
   print *, "COSX force: called before the COSX grid/screening cache ", &
            "was built - SCF must converge with cosx_enabled first"
   call flush(6)
   stop 1
endif

if (RSomega .gt. 0.0d0) then
   if ((.not. cosx_esp_bound_lr_ready) .or. (cosx_esp_bound_lr_omega .ne. RSomega)) then
      if (allocated(cosx_esp_bound_lr)) deallocate(cosx_esp_bound_lr)
      allocate(cosx_esp_bound_lr(0:nBases-1,0:nBases-1))
      cosx_esp_bound_lr = 0.0d0
      !$omp parallel do default(shared) &
      !$omp&   private(si,sj,ni_b,nj_b,offi_b,offj_b,r2_b,pi_b,pj_b,ei_b,ci_b,ej_b,cj_b,p_b) schedule(dynamic)
      do si = 0,nBases-1
         ni_b = bas(3,si+1)
         offi_b = bas(6,si+1)
         do sj = si,nBases-1
            nj_b = bas(3,sj+1)
            offj_b = bas(6,sj+1)
            r2_b = sum((cosx_shell_cen(:,si)-cosx_shell_cen(:,sj))**2)
            do pi_b = 1,ni_b
               ei_b = env(offi_b+pi_b)
               ci_b = env(bas(7,si+1)+pi_b)
               do pj_b = 1,nj_b
                  ej_b = env(offj_b+pj_b)
                  cj_b = env(bas(7,sj+1)+pj_b)
                  p_b = ei_b + ej_b
                  cosx_esp_bound_lr(si,sj) = cosx_esp_bound_lr(si,sj) + &
                     abs(ci_b*cj_b)*exp(-r2_b*ei_b*ej_b/p_b) * &
                     2.0d0*acos(-1.0d0)*RSomega/(p_b*sqrt(p_b+RSomega*RSomega))
               enddo
            enddo
            cosx_esp_bound_lr(sj,si) = cosx_esp_bound_lr(si,sj)
         enddo
      enddo
      !$omp end parallel do
      cosx_esp_bound_lr_max = maxval(cosx_esp_bound_lr)
      call cosx_build_ext_lr(nBases)
      cosx_kscreen_lr = cosx_kscreen
      cosx_esp_bound_lr_ready = .true.
      cosx_esp_bound_lr_omega = RSomega
   endif
   esp_bound_p => cosx_esp_bound_lr
   ext_ptr_p => cosx_ext_ptr_lr;  ext_idx_p => cosx_ext_idx_lr
   kscreen_use = cosx_kscreen_lr
else
   esp_bound_p => cosx_esp_bound
   ext_ptr_p => cosx_ext_ptr;  ext_idx_p => cosx_ext_idx
   kscreen_use = cosx_kscreen
endif

if (RSomega .eq. 0.0d0) then
   use_cache = cosx_G_hi_valid
   if (use_cache) then
      G_hi_a_p => cosx_G_hi_a ; G_hi_b_p => cosx_G_hi_b
   endif
else
   use_cache = cosx_G_hi_valid_lr
   if (use_cache) then
      G_hi_a_p => cosx_G_hi_a_lr ; G_hi_b_p => cosx_G_hi_b_lr
   endif
endif

forceKa = 0.0d0
forceKb = 0.0d0

if (.not. cosx_grad_wt_checked) then
   cosx_grad_wt_checked = .true.
   call get_environment_variable("ENGINE_COSX_GRAD_WT", cosx_grad_wt_env)
   if (len_trim(cosx_grad_wt_env) .gt. 0) cosx_grad_wt = (trim(cosx_grad_wt_env) .eq. "1")
endif
if (.not. cosx_grad_resp_checked) then
   cosx_grad_resp_checked = .true.
   call get_environment_variable("ENGINE_COSX_GRAD_RESP", gr_env)
   if (len_trim(gr_env) .gt. 0) cosx_grad_resp = (trim(gr_env) .eq. "1")
   call get_environment_variable("ENGINE_COSX_GR_AO_OWN",  gr_env); if (len_trim(gr_env)>0) read(gr_env,*) fgr_ao_own
   call get_environment_variable("ENGINE_COSX_GR_AO_GRID", gr_env); if (len_trim(gr_env)>0) read(gr_env,*) fgr_ao_grid
   call get_environment_variable("ENGINE_COSX_GR_ESP_OWN", gr_env); if (len_trim(gr_env)>0) read(gr_env,*) fgr_esp_own
   call get_environment_variable("ENGINE_COSX_GR_ESP_GRID",gr_env); if (len_trim(gr_env)>0) read(gr_env,*) fgr_esp_grid
   call get_environment_variable("ENGINE_COSX_GR_WT",      gr_env); if (len_trim(gr_env)>0) read(gr_env,*) fgr_wt
endif
cosx_grad_resp = cosx_grad_resp .and. (RSomega .eq. 0.0d0) .and. allocated(cosx_atom_hi)
cosx_grad_wt = (cosx_grad_wt .or. cosx_grad_resp) .and. (RSomega .eq. 0.0d0) .and. allocated(cosx_atom_hi)
if (cosx_grad_resp) use_cache = .false.
if (cosx_grad_wt) then
   allocate(xed_a(cosx_npts_hi), xed_b(cosx_npts_hi))
   xed_a = 0.0d0
   xed_b = 0.0d0
endif

maxd = 0
do si = 0,nBases-1
   maxd = max(maxd, cosx_shell_dim(si))
enddo

npts = cosx_npts_hi

blkstart_p => cosx_blkstart_hi
nblocks = cosx_nblk_hi

!$omp parallel default(shared) &
!$omp&   private(ib,bstart,nb,g,si,sj,di,dj,i,j,ao_i,ao_j,shls1e,off_grids,nij,wij) &
!$omp&   private(coor_batch,w_batch,Xb,Xb1,Fa,Fb,Ga,Gb,Ha,Hb,env_g,buf,int_cache) &
!$omp&   private(F_gmax,x_block_max,k_bound,sj_idx,kappa,dloc,atom_idx,sKa,sKb) &
!$omp&   private(forceKa_local,forceKb_local) &
!$omp&   private(dGa,dGb,dbuf,t_ao,t_esp,wij_d,gpar,kdir)
allocate(buf(BATCH*maxd*maxd))
allocate(int_cache(INT_CACHE_SIZE))
allocate(forceKa_local(natoms,3), forceKb_local(natoms,3))
forceKa_local = 0.0d0
forceKb_local = 0.0d0
allocate(env_g(size(env) + 3*BATCH))
env_g(1:size(env)) = env(1:size(env))
off_grids = size(env)
env_g(9) = RSomega

!$omp do schedule(dynamic)
do ib = 0,nblocks-1
   bstart = blkstart_p(ib) + 1
   nb = blkstart_p(ib+1) - blkstart_p(ib)

   allocate(coor_batch(3,nb), w_batch(nb))
   coor_batch = cosx_coor_hi(:,bstart:bstart+nb-1)
   w_batch = cosx_w_hi(bstart:bstart+nb-1)

   allocate(Xb(nConts,nb), Xb1(nConts,3,nb))
   Xb = 0.0d0
   Xb1 = 0.0d0
   do sj_idx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      si = cosx_blkidx_hi(sj_idx)
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                                coor_batch, nb, cosx_ao_offset(si), nConts, Xb, Xb1)
   enddo
   do g = 1,nb
      Xb(:,g) = Xb(:,g) * sqrt(abs(w_batch(g)))
      Xb1(:,:,g) = Xb1(:,:,g) * sqrt(abs(w_batch(g)))
   enddo

   allocate(Ha(nb,nConts), Hb(nb,nConts))
   if (use_cache) then
      call dgemm('N','N',nb,nConts,nConts,1.0d0,G_hi_a_p(bstart:bstart+nb-1,:),nb,Da,nConts,0.0d0,Ha,nb)
      call dgemm('N','N',nb,nConts,nConts,1.0d0,G_hi_b_p(bstart:bstart+nb-1,:),nb,Db,nConts,0.0d0,Hb,nb)
   else
      allocate(Fa(nb,nConts), Fb(nb,nConts))
      call dgemm('T','N',nb,nConts,nConts,1.0d0,Xb,nConts,Da,nConts,0.0d0,Fa,nb)
      call dgemm('T','N',nb,nConts,nConts,1.0d0,Xb,nConts,Db,nConts,0.0d0,Fb,nb)

      x_block_max = maxval(abs(Xb))
      allocate(F_gmax(0:nBases-1))
      do si = 0,nBases-1
         F_gmax(si) = max(maxval(abs(Fa(:,cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si)))), &
                           maxval(abs(Fb(:,cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si)))))
      enddo

      allocate(Ga(nb,nConts), Gb(nb,nConts))
      Ga = 0.0d0
      Gb = 0.0d0
      if (cosx_grad_resp) then
         allocate(dGa(nb,nConts,3), dGb(nb,nConts,3), dbuf(nb*maxd*maxd*3))
         dGa = 0.0d0
         dGb = 0.0d0
      endif

      do g = 1,nb
         env_g(off_grids+3*(g-1)+1) = coor_batch(1,g)
         env_g(off_grids+3*(g-1)+2) = coor_batch(2,g)
         env_g(off_grids+3*(g-1)+3) = coor_batch(3,g)
      enddo
      env_g(12) = dble(nb)
      env_g(13) = dble(off_grids)

      do si = 0,nBases-1
         di = cosx_shell_dim(si)
         do sj_idx = ext_ptr_p(si)+1, ext_ptr_p(si+1)
            sj = ext_idx_p(sj_idx)
            if (schwarz_bound(si,sj) .lt. SCHWARZ_CUTOFF) cycle
            k_bound = x_block_max*esp_bound_p(si,sj)*max(F_gmax(si),F_gmax(sj))
            if (k_bound .lt. kscreen_use) cycle
            dj = cosx_shell_dim(sj)
            shls1e(1) = si
            shls1e(2) = sj
            shls1e(3) = 0
            shls1e(4) = nb
            call grids1e_engine_cached(buf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            do j = 1,dj
               ao_j = cosx_ao_offset(sj) + j
               do i = 1,di
                  ao_i = cosx_ao_offset(si) + i
                  nij = NorVEC(ao_i)*NorVEC(ao_j)
                  if (nij .eq. 0.0d0) cycle
                  do g = 1,nb
                     wij = sign(1.0d0,w_batch(g)) * nij * buf(g + nb*((i-1)+di*(j-1)))
                     Ga(g,ao_i) = Ga(g,ao_i) + wij*Fa(g,ao_j)
                     Gb(g,ao_i) = Gb(g,ao_i) + wij*Fb(g,ao_j)
                     if (sj .ne. si) then
                        Ga(g,ao_j) = Ga(g,ao_j) + wij*Fa(g,ao_i)
                        Gb(g,ao_j) = Gb(g,ao_j) + wij*Fb(g,ao_i)
                     endif
                  enddo
               enddo
            enddo
            if (cosx_grad_resp) then
               call grids1e_ip_engine(dbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8)
               do kdir = 1,3
                  do j = 1,dj
                     ao_j = cosx_ao_offset(sj) + j
                     do i = 1,di
                        ao_i = cosx_ao_offset(si) + i
                        nij = NorVEC(ao_i)*NorVEC(ao_j)
                        if (nij .eq. 0.0d0) cycle
                        do g = 1,nb
                           wij = sign(1.0d0,w_batch(g)) * nij * &
                                 dbuf(g + nb*((i-1) + di*((j-1) + dj*(kdir-1))))
                           dGa(g,ao_i,kdir) = dGa(g,ao_i,kdir) + wij*Fa(g,ao_j)
                           dGb(g,ao_i,kdir) = dGb(g,ao_i,kdir) + wij*Fb(g,ao_j)
                        enddo
                     enddo
                  enddo
               enddo
               if (sj .ne. si) then
                  shls1e(1) = sj ; shls1e(2) = si
                  call grids1e_ip_engine(dbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8)
                  shls1e(1) = si ; shls1e(2) = sj
                  do kdir = 1,3
                     do i = 1,di
                        ao_i = cosx_ao_offset(si) + i
                        do j = 1,dj
                           ao_j = cosx_ao_offset(sj) + j
                           nij = NorVEC(ao_i)*NorVEC(ao_j)
                           if (nij .eq. 0.0d0) cycle
                           do g = 1,nb
                              wij = sign(1.0d0,w_batch(g)) * nij * &
                                    dbuf(g + nb*((j-1) + dj*((i-1) + di*(kdir-1))))
                              dGa(g,ao_j,kdir) = dGa(g,ao_j,kdir) + wij*Fa(g,ao_i)
                              dGb(g,ao_j,kdir) = dGb(g,ao_j,kdir) + wij*Fb(g,ao_i)
                           enddo
                        enddo
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo

      call dgemm('N','N',nb,nConts,nConts,1.0d0,Ga,nb,Da,nConts,0.0d0,Ha,nb)
      call dgemm('N','N',nb,nConts,nConts,1.0d0,Gb,nb,Db,nConts,0.0d0,Hb,nb)
      if (cosx_grad_resp) deallocate(dbuf)
      deallocate(Ga,Gb,F_gmax)
      if (.not. cosx_grad_resp) deallocate(Fa,Fb)
   endif

   if (cosx_grad_wt) then
      do g = 1,nb
         xed_a(bstart+g-1) = dot_product(Xb(1:nConts,g), Ha(g,1:nConts)) / &
                             sign(max(abs(w_batch(g)),1.0d-300), w_batch(g))
         xed_b(bstart+g-1) = dot_product(Xb(1:nConts,g), Hb(g,1:nConts)) / &
                             sign(max(abs(w_batch(g)),1.0d-300), w_batch(g))
      enddo
   endif

   do sj_idx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      si = cosx_blkidx_hi(sj_idx)
      atom_idx = cosx_shell_atom(si)
      di = cosx_shell_dim(si)
      do dloc = 1,3
         if (.not. cosx_grad_resp) then
            sKa = 0.0d0
            sKb = 0.0d0
            do i = 1,di
               kappa = cosx_ao_offset(si) + i
               do g = 1,nb
                  sKa = sKa + Xb1(kappa,dloc,g)*Ha(g,kappa)
                  sKb = sKb + Xb1(kappa,dloc,g)*Hb(g,kappa)
               enddo
            enddo
            forceKa_local(atom_idx,dloc) = forceKa_local(atom_idx,dloc) + sKa
            forceKb_local(atom_idx,dloc) = forceKb_local(atom_idx,dloc) + sKb
         else
            do i = 1,di
               kappa = cosx_ao_offset(si) + i
               do g = 1,nb
                  gpar = cosx_atom_hi(bstart+g-1)
                  t_ao  = Xb1(kappa,dloc,g)*Ha(g,kappa)
                  t_esp = Fa(g,kappa)*dGa(g,kappa,dloc)
                  forceKa_local(atom_idx,dloc) = forceKa_local(atom_idx,dloc) &
                       + fgr_ao_own*t_ao + fgr_esp_own*t_esp
                  forceKa_local(gpar,dloc) = forceKa_local(gpar,dloc) &
                       + fgr_ao_grid*t_ao + fgr_esp_grid*t_esp
                  t_ao  = Xb1(kappa,dloc,g)*Hb(g,kappa)
                  t_esp = Fb(g,kappa)*dGb(g,kappa,dloc)
                  forceKb_local(atom_idx,dloc) = forceKb_local(atom_idx,dloc) &
                       + fgr_ao_own*t_ao + fgr_esp_own*t_esp
                  forceKb_local(gpar,dloc) = forceKb_local(gpar,dloc) &
                       + fgr_ao_grid*t_ao + fgr_esp_grid*t_esp
               enddo
            enddo
         endif
      enddo
   enddo

   deallocate(coor_batch,w_batch,Xb,Xb1,Ha,Hb)
   if (cosx_grad_resp) deallocate(dGa,dGb,Fa,Fb)
enddo
!$omp end do

block
   integer :: merge_tid, merge_nthreads
   merge_nthreads = omp_get_num_threads()
   do merge_tid = 0, merge_nthreads-1
      !$omp barrier
      if (omp_get_thread_num() .eq. merge_tid) then
      forceKa = forceKa + forceKa_local
      forceKb = forceKb + forceKb_local
      endif
   enddo
end block
deallocate(buf,int_cache,env_g,forceKa_local,forceKb_local)
!$omp end parallel

if (cosx_grad_wt) then
   allocate(wt_a(natoms,3), wt_b(natoms,3), eps_zero(cosx_npts_hi))
   eps_zero = 0.0d0
   call vv10_nlc_weight_gradient(cosx_npts_hi, cosx_coor_hi, cosx_atom_hi, &
        max(cosx_npts_hi/max(natoms,1),1), xed_a, eps_zero, cosx_w_hi, -1.0d0, wt_a)
   call vv10_nlc_weight_gradient(cosx_npts_hi, cosx_coor_hi, cosx_atom_hi, &
        max(cosx_npts_hi/max(natoms,1),1), xed_b, eps_zero, cosx_w_hi, -1.0d0, wt_b)
   forceKa = forceKa + fgr_wt*wt_a
   forceKb = forceKb + fgr_wt*wt_b
   deallocate(wt_a, wt_b, eps_zero, xed_a, xed_b)
endif

call cosx_malloc_trim()

end subroutine integrals_force_exchange_cosx
end submodule exchange_cosx_force_impl
