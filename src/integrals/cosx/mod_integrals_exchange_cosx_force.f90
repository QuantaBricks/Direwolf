! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX analytic COSX exchange gradient builder - split out of mod_integrals_exchange_cosx.f90

submodule (mod_integrals:exchange_cosx_impl) exchange_cosx_force_impl
interface
   subroutine vv10_nlc_weight_gradient(np, coor, atom_of, nper, rho, eps, wquad, beta_v, gv_wt, dwdr_out)
   use MOL_info, only: natoms
   implicit none
   integer,intent(in) :: np, nper
   real(8),intent(in) :: coor(3,np)
   integer,intent(in) :: atom_of(np)
   real(8),intent(in) :: rho(np), eps(np), wquad(np), beta_v
   real(8),intent(out) :: gv_wt(natoms,3)
   real(8),optional,intent(out) :: dwdr_out(3,natoms,np)
   end subroutine vv10_nlc_weight_gradient
   subroutine vv10_nlc_weight_hessian(np, coor, atom_of, wquad, fct, hess_out)
   use MOL_info, only: natoms
   implicit none
   integer,intent(in) :: np
   real(8),intent(in) :: coor(3,np)
   integer,intent(in) :: atom_of(np)
   real(8),intent(in) :: wquad(np), fct(np)
   real(8),intent(out) :: hess_out(natoms,3,natoms,3)
   end subroutine vv10_nlc_weight_hessian
end interface
real(8),allocatable,save :: cosx_eps_a_hi_save(:), cosx_eps_b_hi_save(:)
logical,save :: cosx_eps_saved = .false.
real(8),allocatable,save :: cosx_hess_a_save_a(:,:,:,:), cosx_hess_a_save_b(:,:,:,:)
logical,save :: cosx_hess_a_saved = .false.
contains
subroutine cosx_ensure_esp_bound_lr(nBases_in, RSomega_in)
implicit none
integer,intent(in) :: nBases_in
real(8),intent(in) :: RSomega_in
integer :: si, sj, ni_b, nj_b, offi_b, offj_b, pi_b, pj_b
real(8) :: r2_b, ei_b, ci_b, ej_b, cj_b, p_b

if (cosx_esp_bound_lr_ready .and. cosx_esp_bound_lr_omega .eq. RSomega_in) return
if (allocated(cosx_esp_bound_lr)) deallocate(cosx_esp_bound_lr)
allocate(cosx_esp_bound_lr(0:nBases_in-1,0:nBases_in-1))
cosx_esp_bound_lr = 0.0d0
!$omp parallel do default(shared) &
!$omp&   private(si,sj,ni_b,nj_b,offi_b,offj_b,r2_b,pi_b,pj_b,ei_b,ci_b,ej_b,cj_b,p_b) schedule(dynamic)
do si = 0,nBases_in-1
   ni_b = bas(3,si+1)
   offi_b = bas(6,si+1)
   do sj = si,nBases_in-1
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
               2.0d0*acos(-1.0d0)*RSomega_in/(p_b*sqrt(p_b+RSomega_in*RSomega_in))
         enddo
      enddo
      cosx_esp_bound_lr(sj,si) = cosx_esp_bound_lr(si,sj)
   enddo
enddo
!$omp end parallel do
cosx_esp_bound_lr_max = maxval(cosx_esp_bound_lr)
call cosx_build_ext_lr(nBases_in)
cosx_kscreen_lr = cosx_kscreen
cosx_esp_bound_lr_ready = .true.
cosx_esp_bound_lr_omega = RSomega_in
end subroutine cosx_ensure_esp_bound_lr

subroutine cosx_force_exchange_v2(nConts, natoms, Da, Db, RSomega, forceKa, forceKb)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_get_wtime
use iso_c_binding, only: c_int, c_double
implicit none
interface
   function cosx_grad_contract(shls, atm, natm, bas, nbas, env, sph, nspin, &
                               vi, vj, sgw, oi, oj, ot, cache) bind(c, name="cosx_grad_contract") result(r)
   import :: c_int, c_double
   integer(c_int) :: r
   integer(c_int) :: shls(*), atm(*), bas(*)
   integer(c_int),value :: natm, nbas, sph, nspin
   real(c_double) :: env(*), vi(*), vj(*), sgw(*), oi(*), oj(*), ot(*), cache(*)
   end function cosx_grad_contract
end interface
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
real(8),intent(in) :: RSomega
real(8),intent(out) :: forceKa(natoms,3), forceKb(natoms,3)

real(8),pointer,contiguous :: Pfit(:,:)
logical :: fit
integer :: npts, nblocks, maxd, ib, bstart, nb, g, gg, d, ee
integer :: si, sj, di, dj, i, j, ao_i, ao_j, sj_idx, shls1e(4), off_grids
integer :: gpar, atom_i, atom_j, kappa, aq, dirq, kq, lq, info_q
real(8) :: nij, wij, sgw, theta, thetb, tval, dbra_v, dket_v, x_block_max, k_bound
real(8) :: seg_ab(3), seg_cl(3), seg_len2, seg_t, r_far, esp_eff
real(8),pointer :: esp_bound_p(:,:)
integer,pointer :: ext_ptr_p(:), ext_idx_p(:)
integer(8) :: n_offer, n_keep
real(8) :: tm_ao, tm_gemm, tm_espv, tm_espd, tm_gacc, tm_r3, tm_r2, tm_eps, tm_r4, tm0
real(8) :: tw0, tw_loop, tw_r1
logical :: tmg
real(8) :: kscreen_grad
real(8),parameter :: COSX_GRAD_KSCREEN_SCALE = 10.0d0
real(8) :: kscreen_pair
real(8),parameter :: COSX_GRAD_KPAIR_RATIO = 1.0d-2
logical :: do_esp_d
real(8),allocatable :: coor_batch(:,:), w_batch(:)
real(8),allocatable :: Xb(:,:), Xb1(:,:,:), Xb2(:,:,:,:), Qb(:,:)
real(8),allocatable :: Fa(:,:), Fb(:,:), Pha(:,:), Phb(:,:)
real(8),allocatable :: Ga(:,:), Gb(:,:), Gta(:,:), Gtb(:,:)
real(8),allocatable :: Ha(:,:), Hb(:,:), Hta(:,:), Htb(:,:)
real(8),allocatable :: Psia(:,:), Psib(:,:)
real(8),allocatable :: env_g(:), buf(:), dbuf(:), dbuf2(:), int_cache(:), vbuf(:)
logical :: grad_fused
integer :: nn, ii, jj
logical :: grad_contract
integer :: nsp, isp, rc_c
real(8),allocatable :: cvi(:), cvj(:), coi(:), coj(:), cot(:), csg(:)
real(8),allocatable :: F_gmax(:), P_gmax(:)
real(8),allocatable :: forceKa_local(:,:), forceKb_local(:,:)
real(8),allocatable :: pt_a(:,:), pt_b(:,:)
real(8),allocatable :: tha(:), thb(:)
logical :: same_spin
real(8) :: sKa, sKb, sJa, sJb
integer :: ob, okk
real(8),allocatable :: Va(:,:), Vb(:,:), Va_local(:,:), Vb_local(:,:)
real(8),allocatable :: Ma(:,:), Mb(:,:), DPa(:,:), DPb(:,:)
real(8),allocatable :: Xc(:,:), Dc(:,:), Mc(:,:), Tc(:,:)
integer,allocatable :: aoidx(:)
integer :: nnz
real(8),allocatable :: eps_a(:), eps_b(:), zero_eps(:), wt_a(:,:), wt_b(:,:)
real(8),allocatable :: Xr(:,:), Xr1(:,:,:), dsnum(:,:,:,:), dSan(:,:)
real(8),allocatable :: Yq(:,:), Wq(:,:), WX(:,:), WTX(:,:)
real(8),allocatable :: termA(:,:), termB(:,:), uq(:), wtq(:,:), fk_save(:,:)
integer :: ispin, gq
real(8),allocatable :: rhs_q(:,:), Zq(:,:), dS1e(:,:,:), dH1e(:,:,:)
real(8),allocatable :: wq_dw(:,:,:), wq_zero(:), wt_dummy(:,:)
integer,allocatable :: ao_atom(:)
integer,parameter :: INT_CACHE_SIZE = 4000000
logical :: hess_v2
real(8),allocatable :: dbpp(:), dbvp(:), dbpp2(:)
real(8),allocatable :: ha_local(:,:,:,:), hb_local(:,:,:,:)
real(8),allocatable :: ha_glob(:,:,:,:), hb_glob(:,:,:,:)
real(8) :: t_h, ipp, ivp, ivps, ikk
integer :: obpp, obpp2, obvp, obvps
real(8),allocatable :: SegAi(:), SegAj(:), SegBi(:), SegBj(:), Wmat_i(:,:), Wmat_j(:,:)
real(8) :: Gfdi, Gfdj, dTh, sgij, ordf, dAbr(3), dAkt(3)
real(8) :: Gx_i(3), Gx_j(3), dGai, dGaj, dGtai, dGtaj, dFaC, dFaiC, dFajC, dAijC, wij_v, dAbr_dd, dAkt_dd, sgnij_v
integer :: l2, a2, dd2, cc2, cc3, kappa2, sjj2, sk2, kk2
real(8) :: xlq
integer :: r3p_atom, r3p_dir
real(8) :: psi_dev, psi_scale
logical :: psi_chk
logical :: r3p_on
real(8) :: r3p_num, r3p_an, r3p_h
real(8),allocatable :: envp1(:), envp2(:), bufp1(:), bufp2(:)
real(8) :: c_x, c_a, c_w, c_q
character(len=16) :: v2env
logical :: frozen

forceKa = 0.0d0
forceKb = 0.0d0
c_x = 1.0d0 ; c_a = 1.0d0 ; c_w = 1.0d0 ; c_q = 1.0d0
call get_environment_variable("ENGINE_COSX_GRAD_FROZEN", v2env)
frozen = (trim(v2env) .eq. "1")
if (frozen) c_w = 0.0d0
call get_environment_variable("ENGINE_COSX_V2_CX", v2env); if (len_trim(v2env)>0) read(v2env,*) c_x
call get_environment_variable("ENGINE_COSX_V2_CA", v2env); if (len_trim(v2env)>0) read(v2env,*) c_a
call get_environment_variable("ENGINE_COSX_V2_CW", v2env); if (len_trim(v2env)>0) read(v2env,*) c_w
call get_environment_variable("ENGINE_COSX_V2_CQ", v2env); if (len_trim(v2env)>0) read(v2env,*) c_q
psi_dev = 0.0d0 ; psi_scale = 0.0d0
n_offer = 0_8 ; n_keep = 0_8
tm_ao=0d0; tm_gemm=0d0; tm_espv=0d0; tm_espd=0d0; tm_gacc=0d0
tm_r3=0d0; tm_r2=0d0; tm_eps=0d0; tm_r4=0d0
call get_environment_variable("ENGINE_COSX_V2_TIMING", v2env)
tmg = (trim(v2env) .eq. "1")
kscreen_grad = COSX_GRAD_KSCREEN_SCALE*cosx_kscreen
call get_environment_variable("ENGINE_COSX_GRAD_KSCREEN", v2env)
if (len_trim(v2env) .gt. 0) read(v2env,*) kscreen_grad
kscreen_pair = COSX_GRAD_KPAIR_RATIO*kscreen_grad
call get_environment_variable("ENGINE_COSX_GRAD_KPAIR", v2env)
if (len_trim(v2env) .gt. 0) read(v2env,*) kscreen_pair
if (RSomega .gt. 0.0d0) then
   call get_environment_variable("ENGINE_COSX_GRAD_KSCREEN_LR", v2env)
   if (len_trim(v2env) .gt. 0) read(v2env,*) kscreen_grad
   call get_environment_variable("ENGINE_COSX_GRAD_KPAIR_LR", v2env)
   if (len_trim(v2env) .gt. 0) then
      read(v2env,*) kscreen_pair
   else
      kscreen_pair = COSX_GRAD_KPAIR_RATIO*kscreen_grad
   endif
endif
call get_environment_variable("ENGINE_COSX_V2_ECHK", v2env)
psi_chk = (trim(v2env) .eq. "1")
r3p_on = .false.
r3p_atom = 0 ; r3p_dir = 0
r3p_num = 0.0d0 ; r3p_an = 0.0d0 ; r3p_h = 1.0d-4
call get_environment_variable("ENGINE_COSX_V2_R3PROBE", v2env)
if (len_trim(v2env) .gt. 0) then
   read(v2env,*) r3p_atom, r3p_dir
   r3p_on = (r3p_atom .ge. 1 .and. r3p_dir .ge. 1 .and. r3p_dir .le. 3)
endif
npts = cosx_npts_hi
nblocks = cosx_nblk_hi
grad_fused = .true.
grad_contract = .not. r3p_on

same_spin = all(Da .eq. Db)

fit = cosx_overlap_fit .and. allocated(cosx_QfullT_hi)
if (fit) fit = (size(cosx_QfullT_hi,1) .eq. nConts)
if (fit) Pfit => cosx_QfullT_hi

maxd = 0
do si = 0,nBases-1
   maxd = max(maxd, cosx_shell_dim(si))
enddo

hess_v2 = .false.
call get_environment_variable("ENGINE_COSX_HESS_V2", v2env)
if (trim(v2env) .eq. "1") hess_v2 = .true.
if (hess_v2) then
   allocate(ha_glob(3,natoms,3,natoms), hb_glob(3,natoms,3,natoms))
   ha_glob = 0.0d0 ; hb_glob = 0.0d0
endif

allocate(ao_atom(nConts))
do si = 0,nBases-1
   do i = 1,cosx_shell_dim(si)
      ao_atom(cosx_ao_offset(si)+i) = cosx_shell_atom(si)
   enddo
enddo

allocate(eps_a(npts), eps_b(npts))
eps_a = 0.0d0
eps_b = 0.0d0
if (fit) then
   allocate(Va(nConts,nConts), Vb(nConts,nConts))
   Va = 0.0d0
   Vb = 0.0d0
   allocate(Ma(nConts,nConts), DPa(nConts,nConts))
   call dgemm('T','N',nConts,nConts,nConts,1.0d0,Pfit,nConts,Da,nConts,0.0d0,Ma,nConts)
   call dgemm('N','N',nConts,nConts,nConts,1.0d0,Da,nConts,Pfit,nConts,0.0d0,DPa,nConts)
   if (.not. same_spin) then
      allocate(Mb(nConts,nConts), DPb(nConts,nConts))
      call dgemm('T','N',nConts,nConts,nConts,1.0d0,Pfit,nConts,Db,nConts,0.0d0,Mb,nConts)
      call dgemm('N','N',nConts,nConts,nConts,1.0d0,Db,nConts,Pfit,nConts,0.0d0,DPb,nConts)
   endif
endif

tw0 = omp_get_wtime()
!$omp parallel default(shared) &
!$omp&   private(ib,bstart,nb,g,gg,d,si,sj,di,dj,i,j,ao_i,ao_j,sj_idx,shls1e,off_grids) &
!$omp&   private(gpar,atom_i,atom_j,kappa,nij,wij,sgw,theta,tval,dbra_v,dket_v) &
!$omp&   private(x_block_max,k_bound,coor_batch,w_batch,Xb,Xb1,Xb2,Qb,Fa,Fb,Pha,Phb) &
!$omp&   private(Ga,Gb,Gta,Gtb,Ha,Hb,Hta,Htb,Psia,Psib,env_g,buf,dbuf,dbuf2,vbuf,nn,ii,jj) &
!$omp&   private(cvi,cvj,coi,coj,cot,csg,isp,rc_c,Xc,Dc,Mc,Tc,aoidx,nnz) &
!$omp&   private(int_cache,F_gmax,P_gmax,forceKa_local,forceKb_local,Va_local,Vb_local,thetb) &
!$omp&   private(dbpp,dbvp,dbpp2,ha_local,hb_local,SegAi,SegAj,SegBi,SegBj,Wmat_i,Wmat_j) &
!$omp&   private(dAbr_dd,dAkt_dd,dFaiC,dFajC,sgnij_v,ordf,Gfdi,Gfdj,dTh,sgij,dAbr,dAkt) &
!$omp&   private(l2,a2,dd2,cc2,cc3,kappa2,sjj2,sk2,kk2,Gx_i,Gx_j,dGai,dGaj,dGtai,dGtaj,dFaC,dAijC,wij_v) &
!$omp&   private(envp1,envp2,bufp1,bufp2,pt_a,pt_b,tha,thb) &
!$omp&   private(sKa,sKb,sJa,sJb,ob,okk,ee,seg_ab,seg_cl,seg_len2,seg_t,r_far,esp_eff) &
!$omp&   reduction(+:r3p_num,r3p_an) &
!$omp&   reduction(max:psi_dev,psi_scale) reduction(+:n_offer,n_keep) shared(psi_chk,same_spin) &
!$omp&   reduction(+:tm_ao,tm_gemm,tm_espv,tm_espd,tm_gacc,tm_r3,tm_r2,tm_eps) &
!$omp&   private(tm0,do_esp_d) shared(tmg,kscreen_grad,kscreen_pair)
allocate(buf(BATCH*maxd*maxd))
allocate(dbuf(BATCH*maxd*maxd*3), dbuf2(BATCH*maxd*maxd*3))
if (grad_fused) allocate(vbuf(BATCH*maxd*maxd*7))
if (grad_contract) allocate(cvi(BATCH*maxd*4), cvj(BATCH*maxd*4), coi(BATCH*maxd*4), &
                            coj(BATCH*maxd*4), cot(BATCH*12), csg(BATCH))
allocate(int_cache(INT_CACHE_SIZE))
allocate(forceKa_local(natoms,3), forceKb_local(natoms,3))
forceKa_local = 0.0d0
forceKb_local = 0.0d0
if (hess_v2) then
   allocate(dbpp(BATCH*maxd*maxd*9), dbvp(BATCH*maxd*maxd*9), dbpp2(BATCH*maxd*maxd*9))
   allocate(ha_local(3,natoms,3,natoms), hb_local(3,natoms,3,natoms))
   allocate(SegAi(natoms), SegAj(natoms), SegBi(natoms), SegBj(natoms), Wmat_i(natoms,3), Wmat_j(natoms,3))
   ha_local = 0.0d0 ; hb_local = 0.0d0
endif
if (fit) then
   allocate(Va_local(nConts,nConts), Vb_local(nConts,nConts))
   Va_local = 0.0d0
   Vb_local = 0.0d0
endif
if (RSomega .gt. 0.0d0) then
   call cosx_ensure_esp_bound_lr(nBases, RSomega)
   esp_bound_p => cosx_esp_bound_lr
   ext_ptr_p => cosx_ext_ptr_lr
   ext_idx_p => cosx_ext_idx_lr
else
   esp_bound_p => cosx_esp_bound
   ext_ptr_p => cosx_ext_ptr
   ext_idx_p => cosx_ext_idx
endif
allocate(env_g(size(env) + 3*BATCH))
env_g(1:size(env)) = env(1:size(env))
off_grids = size(env)
env_g(9) = RSomega
allocate(pt_a(BATCH,3), pt_b(BATCH,3), tha(BATCH), thb(BATCH))
if (r3p_on) then
   allocate(envp1(size(env_g)), envp2(size(env_g)))
   allocate(bufp1(BATCH*maxd*maxd), bufp2(BATCH*maxd*maxd))
endif

!$omp do schedule(dynamic)
do ib = 0,nblocks-1
   bstart = cosx_blkstart_hi(ib) + 1
   nb = cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib)

   allocate(coor_batch(3,nb), w_batch(nb))
   coor_batch = cosx_coor_hi(:,bstart:bstart+nb-1)
   w_batch = cosx_w_hi(bstart:bstart+nb-1)
   if (grad_contract) then
      do g = 1,nb
         csg(g) = sign(1.0d0,w_batch(g))
      enddo
   endif

   if (tmg) tm0 = omp_get_wtime()
   allocate(Xb(nConts,nb), Xb1(nConts,3,nb))
   Xb = 0.0d0
   Xb1 = 0.0d0
   if (hess_v2) then
      allocate(Xb2(nConts,3,3,nb))
      Xb2 = 0.0d0
   endif
   do sj_idx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      si = cosx_blkidx_hi(sj_idx)
      if (hess_v2) then
         call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                                   coor_batch, nb, cosx_ao_offset(si), nConts, Xb, Xb1, Xb2)
      else
         call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                                   coor_batch, nb, cosx_ao_offset(si), nConts, Xb, Xb1)
      endif
   enddo
   do g = 1,nb
      Xb(:,g) = Xb(:,g) * sqrt(abs(w_batch(g)))
      Xb1(:,:,g) = Xb1(:,:,g) * sqrt(abs(w_batch(g)))
      if (hess_v2) Xb2(:,:,:,g) = Xb2(:,:,:,g) * sqrt(abs(w_batch(g)))
   enddo
   if (tmg) tm_ao = tm_ao + omp_get_wtime() - tm0

   if (tmg) tm0 = omp_get_wtime()
   allocate(Fa(nb,nConts), Fb(nb,nConts), Pha(nb,nConts), Phb(nb,nConts))
   nnz = 0
   do sj_idx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      nnz = nnz + cosx_shell_dim(cosx_blkidx_hi(sj_idx))
   enddo
   allocate(aoidx(nnz), Xc(nnz,nb), Dc(nnz,nConts))
   nnz = 0
   do sj_idx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      si = cosx_blkidx_hi(sj_idx)
      do i = 1,cosx_shell_dim(si)
         aoidx(nnz+i) = cosx_ao_offset(si) + i
      enddo
      nnz = nnz + cosx_shell_dim(si)
   enddo
   do g = 1,nb
      Xc(:,g) = Xb(aoidx,g)
   enddo
   do j = 1,nConts
      Dc(:,j) = Da(aoidx,j)
   enddo
   call dgemm('T','N',nb,nConts,nnz,1.0d0,Xc,nnz,Dc,nnz,0.0d0,Fa,nb)
   if (.not. same_spin) then
      do j = 1,nConts
         Dc(:,j) = Db(aoidx,j)
      enddo
      call dgemm('T','N',nb,nConts,nnz,1.0d0,Xc,nnz,Dc,nnz,0.0d0,Fb,nb)
   endif
   if (fit) then
      do j = 1,nConts
         Dc(:,j) = Ma(aoidx,j)
      enddo
      call dgemm('T','N',nb,nConts,nnz,1.0d0,Xc,nnz,Dc,nnz,0.0d0,Pha,nb)
      if (.not. same_spin) then
         do j = 1,nConts
            Dc(:,j) = Mb(aoidx,j)
         enddo
         call dgemm('T','N',nb,nConts,nnz,1.0d0,Xc,nnz,Dc,nnz,0.0d0,Phb,nb)
      endif
   else
      Pha = Fa
      if (.not. same_spin) Phb = Fb
   endif

   if (psi_chk .and. fit) then
      psi_dev = max(psi_dev, maxval(abs(Pha))/max(maxval(abs(Fa)),1.0d-300))
      psi_scale = max(psi_scale, maxval(abs(Pfit)))
   endif

   if (tmg) tm_gemm = tm_gemm + omp_get_wtime() - tm0
   x_block_max = maxval(abs(Xb))
   allocate(F_gmax(0:nBases-1), P_gmax(0:nBases-1))
   do si = 0,nBases-1
      F_gmax(si) = maxval(abs(Fa(:,cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si))))
      P_gmax(si) = maxval(abs(Pha(:,cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si))))
      if (.not. same_spin) then
         F_gmax(si) = max(F_gmax(si), &
            maxval(abs(Fb(:,cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si)))))
         P_gmax(si) = max(P_gmax(si), &
            maxval(abs(Phb(:,cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si)))))
      endif
      F_gmax(si) = max(F_gmax(si), P_gmax(si))
   enddo

   allocate(Ga(nb,nConts), Gb(nb,nConts), Gta(nb,nConts), Gtb(nb,nConts))
   Ga = 0.0d0
   Gb = 0.0d0
   Gta = 0.0d0
   Gtb = 0.0d0
   pt_a(1:nb,1:3) = 0.0d0
   pt_b(1:nb,1:3) = 0.0d0

   do g = 1,nb
      env_g(off_grids+3*(g-1)+1) = coor_batch(1,g)
      env_g(off_grids+3*(g-1)+2) = coor_batch(2,g)
      env_g(off_grids+3*(g-1)+3) = coor_batch(3,g)
   enddo
   env_g(12) = dble(nb)
   env_g(13) = dble(off_grids)

   do si = 0,nBases-1
      di = cosx_shell_dim(si)
      atom_i = cosx_shell_atom(si)
      do sj_idx = ext_ptr_p(si)+1, ext_ptr_p(si+1)
         sj = ext_idx_p(sj_idx)
         n_offer = n_offer + 1_8
         if (schwarz_bound(si,sj) .lt. SCHWARZ_CUTOFF) cycle
         seg_ab = cosx_shell_cen(:,sj) - cosx_shell_cen(:,si)
         seg_len2 = sum(seg_ab*seg_ab)
         if (seg_len2 .gt. 0.0d0) then
            seg_t = sum((cosx_blkcen_hi(:,ib) - cosx_shell_cen(:,si))*seg_ab)/seg_len2
            seg_t = max(0.0d0, min(1.0d0, seg_t))
            seg_cl = cosx_shell_cen(:,si) + seg_t*seg_ab
         else
            seg_cl = cosx_shell_cen(:,si)
         endif
         r_far = sqrt(sum((cosx_blkcen_hi(:,ib) - seg_cl)**2)) - cosx_blkrad_hi(ib)
         if (r_far .gt. 0.0d0) then
            esp_eff = min(esp_bound_p(si,sj), cosx_esp_mono(si,sj)/r_far)
         else
            esp_eff = esp_bound_p(si,sj)
         endif
         k_bound = x_block_max*esp_eff*max(F_gmax(si),F_gmax(sj))
         if (k_bound .lt. kscreen_grad) cycle
         n_keep = n_keep + 1_8
         dj = cosx_shell_dim(sj)
         atom_j = cosx_shell_atom(sj)
         shls1e(1) = si ; shls1e(2) = sj ; shls1e(3) = 0 ; shls1e(4) = nb
         do_esp_d = (esp_eff*(P_gmax(si)*F_gmax(sj) + P_gmax(sj)*F_gmax(si)) &
                     .ge. kscreen_pair)
         if (do_esp_d .and. grad_contract .and. .not. hess_v2) then
            if (tmg) tm0 = omp_get_wtime()
            nsp = merge(1, 2, same_spin)
            do isp = 1,nsp
               do i = 1,di
                  ao_i = cosx_ao_offset(si) + i
                  ii = nb*((i-1) + di*2*(isp-1))
                  if (isp .eq. 1) then
                     cvi(ii+1:ii+nb)         = NorVEC(ao_i)*Fa(1:nb,ao_i)
                     cvi(ii+nb*di+1:ii+nb*di+nb) = NorVEC(ao_i)*Pha(1:nb,ao_i)
                  else
                     cvi(ii+1:ii+nb)         = NorVEC(ao_i)*Fb(1:nb,ao_i)
                     cvi(ii+nb*di+1:ii+nb*di+nb) = NorVEC(ao_i)*Phb(1:nb,ao_i)
                  endif
               enddo
               do j = 1,dj
                  ao_j = cosx_ao_offset(sj) + j
                  jj = nb*((j-1) + dj*2*(isp-1))
                  if (isp .eq. 1) then
                     cvj(jj+1:jj+nb)         = NorVEC(ao_j)*Fa(1:nb,ao_j)
                     cvj(jj+nb*dj+1:jj+nb*dj+nb) = NorVEC(ao_j)*Pha(1:nb,ao_j)
                  else
                     cvj(jj+1:jj+nb)         = NorVEC(ao_j)*Fb(1:nb,ao_j)
                     cvj(jj+nb*dj+1:jj+nb*dj+nb) = NorVEC(ao_j)*Phb(1:nb,ao_j)
                  endif
               enddo
            enddo
            rc_c = cosx_grad_contract(shls1e, atm, size(atm,2), bas, nBases, env_g, &
                                      merge(1,0,engine_puream), nsp, cvi, cvj, csg, coi, coj, cot, int_cache)
            if (rc_c .ne. 0) then
               print *, "COSX force: cosx_grad_contract needs single-contraction shells (rc=", rc_c, ")"
               stop 1
            endif
            do isp = 1,nsp
               do i = 1,di
                  ao_i = cosx_ao_offset(si) + i
                  ii = nb*((i-1) + di*2*(isp-1))
                  if (isp .eq. 1) then
                     Ga(1:nb,ao_i)  = Ga(1:nb,ao_i)  + NorVEC(ao_i)*coi(ii+1:ii+nb)
                     Gta(1:nb,ao_i) = Gta(1:nb,ao_i) + NorVEC(ao_i)*coi(ii+nb*di+1:ii+nb*di+nb)
                  else
                     Gb(1:nb,ao_i)  = Gb(1:nb,ao_i)  + NorVEC(ao_i)*coi(ii+1:ii+nb)
                     Gtb(1:nb,ao_i) = Gtb(1:nb,ao_i) + NorVEC(ao_i)*coi(ii+nb*di+1:ii+nb*di+nb)
                  endif
               enddo
               if (sj .ne. si) then
                  do j = 1,dj
                     ao_j = cosx_ao_offset(sj) + j
                     jj = nb*((j-1) + dj*2*(isp-1))
                     if (isp .eq. 1) then
                        Ga(1:nb,ao_j)  = Ga(1:nb,ao_j)  + NorVEC(ao_j)*coj(jj+1:jj+nb)
                        Gta(1:nb,ao_j) = Gta(1:nb,ao_j) + NorVEC(ao_j)*coj(jj+nb*dj+1:jj+nb*dj+nb)
                     else
                        Gb(1:nb,ao_j)  = Gb(1:nb,ao_j)  + NorVEC(ao_j)*coj(jj+1:jj+nb)
                        Gtb(1:nb,ao_j) = Gtb(1:nb,ao_j) + NorVEC(ao_j)*coj(jj+nb*dj+1:jj+nb*dj+nb)
                     endif
                  enddo
               else
                  do j = 1,dj
                     ao_j = cosx_ao_offset(sj) + j
                     jj = nb*((j-1) + dj*2*(isp-1))
                     if (isp .eq. 1) then
                        Gta(1:nb,ao_j) = Gta(1:nb,ao_j) + NorVEC(ao_j)*coj(jj+nb*dj+1:jj+nb*dj+nb)
                     else
                        Gtb(1:nb,ao_j) = Gtb(1:nb,ao_j) + NorVEC(ao_j)*coj(jj+nb*dj+1:jj+nb*dj+nb)
                     endif
                  enddo
               endif
            enddo
            do d = 1,3
               sKa = sum(cot(nb*(d-1)+1:nb*d))
               sJa = sum(cot(nb*(d+2)+1:nb*(d+3)))
               pt_a(1:nb,d) = pt_a(1:nb,d) + c_a*(cot(nb*(d-1)+1:nb*d) + cot(nb*(d+2)+1:nb*(d+3)))
               if (same_spin) then
                  sKb = sKa ; sJb = sJa
                  pt_b(1:nb,d) = pt_b(1:nb,d) + c_a*(cot(nb*(d-1)+1:nb*d) + cot(nb*(d+2)+1:nb*(d+3)))
               else
                  sKb = sum(cot(nb*(d+5)+1:nb*(d+6)))
                  sJb = sum(cot(nb*(d+8)+1:nb*(d+9)))
                  pt_b(1:nb,d) = pt_b(1:nb,d) + c_a*(cot(nb*(d+5)+1:nb*(d+6)) + cot(nb*(d+8)+1:nb*(d+9)))
               endif
               forceKa_local(atom_i,d) = forceKa_local(atom_i,d) + 0.25d0*c_a*sKa
               forceKa_local(atom_j,d) = forceKa_local(atom_j,d) + 0.25d0*c_a*sJa
               forceKb_local(atom_i,d) = forceKb_local(atom_i,d) + 0.25d0*c_a*sKb
               forceKb_local(atom_j,d) = forceKb_local(atom_j,d) + 0.25d0*c_a*sJb
            enddo
            if (tmg) tm_espd = tm_espd + omp_get_wtime() - tm0
            cycle
         endif
         if (do_esp_d .and. grad_fused .and. .not. hess_v2) then
            if (tmg) tm0 = omp_get_wtime()
            call grids1e_vipj_engine_cached(vbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            nn = nb*di*dj
            buf(1:nn) = vbuf(1:nn)
            dbuf(1:3*nn) = vbuf(nn+1:4*nn)
            if (sj .ne. si) then
               do d = 1,3
                  do jj = 1,dj
                     do ii = 1,di
                        dbuf2(nb*((jj-1)+dj*((ii-1)+di*(d-1)))+1:nb*((jj-1)+dj*((ii-1)+di*(d-1)))+nb) = &
                           vbuf(nn*(3+d)+nb*((ii-1)+di*(jj-1))+1:nn*(3+d)+nb*((ii-1)+di*(jj-1))+nb)
                     enddo
                  enddo
               enddo
            endif
            if (tmg) tm_espd = tm_espd + omp_get_wtime() - tm0
         else
         if (tmg) tm0 = omp_get_wtime()
         call grids1e_engine_cached(buf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         if (tmg) tm_espv = tm_espv + omp_get_wtime() - tm0
         if (do_esp_d) then
         if (tmg) tm0 = omp_get_wtime()
         call grids1e_ip_engine_cached(dbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         if (sj .ne. si) then
            shls1e(1) = sj ; shls1e(2) = si
            call grids1e_ip_engine_cached(dbuf2, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            shls1e(1) = si ; shls1e(2) = sj
         endif
         if (hess_v2) then
            call grids1e_ipip_engine_cached(dbpp, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            call grids1e_ipvip_engine_cached(dbvp, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            if (sj .ne. si) then
               shls1e(1) = sj ; shls1e(2) = si
               call grids1e_ipip_engine_cached(dbpp2, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               shls1e(1) = si ; shls1e(2) = sj
            else
               dbpp2(1:nb*di*dj*9) = dbpp(1:nb*di*dj*9)
            endif
         endif
         if (tmg) tm_espd = tm_espd + omp_get_wtime() - tm0
         endif
         endif

         if (tmg) tm0 = omp_get_wtime()
         do j = 1,dj
            ao_j = cosx_ao_offset(sj) + j
            do i = 1,di
               ao_i = cosx_ao_offset(si) + i
               nij = NorVEC(ao_i)*NorVEC(ao_j)
               if (nij .eq. 0.0d0) cycle
               do g = 1,nb
                  sgw = sign(1.0d0,w_batch(g))
                  wij = sgw * nij * buf(g + nb*((i-1)+di*(j-1)))
                  Ga(g,ao_i)  = Ga(g,ao_i)  + wij*Fa(g,ao_j)
                  Gta(g,ao_j) = Gta(g,ao_j) + wij*Pha(g,ao_i)
                  if (sj .ne. si) then
                     Ga(g,ao_j)  = Ga(g,ao_j)  + wij*Fa(g,ao_i)
                     Gta(g,ao_i) = Gta(g,ao_i) + wij*Pha(g,ao_j)
                  endif
                  if (.not. same_spin) then
                     Gb(g,ao_i)  = Gb(g,ao_i)  + wij*Fb(g,ao_j)
                     Gtb(g,ao_j) = Gtb(g,ao_j) + wij*Phb(g,ao_i)
                     if (sj .ne. si) then
                        Gb(g,ao_j)  = Gb(g,ao_j)  + wij*Fb(g,ao_i)
                        Gtb(g,ao_i) = Gtb(g,ao_i) + wij*Phb(g,ao_j)
                     endif
                  endif
               enddo
            enddo
         enddo

         if (tmg) tm_gacc = tm_gacc + omp_get_wtime() - tm0
         if (.not. do_esp_d) cycle
         if (tmg) tm0 = omp_get_wtime()
         do j = 1,dj
            ao_j = cosx_ao_offset(sj) + j
            do i = 1,di
               ao_i = cosx_ao_offset(si) + i
               nij = NorVEC(ao_i)*NorVEC(ao_j)
               if (nij .eq. 0.0d0) cycle
               if (sj .ne. si) then
                  do g = 1,nb
                     sgw = sign(1.0d0,w_batch(g)) * nij
                     tha(g) = sgw*(Pha(g,ao_i)*Fa(g,ao_j) + Pha(g,ao_j)*Fa(g,ao_i))
                  enddo
                  if (same_spin) then
                     thb(1:nb) = tha(1:nb)
                  else
                     do g = 1,nb
                        sgw = sign(1.0d0,w_batch(g)) * nij
                        thb(g) = sgw*(Phb(g,ao_i)*Fb(g,ao_j) + Phb(g,ao_j)*Fb(g,ao_i))
                     enddo
                  endif
               else
                  do g = 1,nb
                     sgw = sign(1.0d0,w_batch(g)) * nij
                     tha(g) = sgw*Pha(g,ao_i)*Fa(g,ao_j)
                  enddo
                  if (same_spin) then
                     thb(1:nb) = tha(1:nb)
                  else
                     do g = 1,nb
                        sgw = sign(1.0d0,w_batch(g)) * nij
                        thb(g) = sgw*Phb(g,ao_i)*Fb(g,ao_j)
                     enddo
                  endif
               endif
               do d = 1,3
                  ob = nb*((i-1) + di*((j-1) + dj*(d-1)))
                  sKa = 0.0d0 ; sKb = 0.0d0 ; sJa = 0.0d0 ; sJb = 0.0d0
                  if (sj .ne. si) then
                     okk = nb*((j-1) + dj*((i-1) + di*(d-1)))
                     do g = 1,nb
                        dbra_v = dbuf(ob+g)
                        dket_v = dbuf2(okk+g)
                        sKa = sKa + tha(g)*dbra_v
                        sJa = sJa + tha(g)*dket_v
                        sKb = sKb + thb(g)*dbra_v
                        sJb = sJb + thb(g)*dket_v
                        pt_a(g,d) = pt_a(g,d) + c_a*tha(g)*(dbra_v+dket_v)
                        pt_b(g,d) = pt_b(g,d) + c_a*thb(g)*(dbra_v+dket_v)
                     enddo
                  else
                     okk = nb*((j-1) + di*((i-1) + dj*(d-1)))
                     do g = 1,nb
                        dbra_v = dbuf(ob+g)
                        dket_v = dbuf(okk+g)
                        sKa = sKa + tha(g)*dbra_v
                        sJa = sJa + tha(g)*dket_v
                        sKb = sKb + thb(g)*dbra_v
                        sJb = sJb + thb(g)*dket_v
                        pt_a(g,d) = pt_a(g,d) + c_a*tha(g)*(dbra_v+dket_v)
                        pt_b(g,d) = pt_b(g,d) + c_a*thb(g)*(dbra_v+dket_v)
                     enddo
                  endif
                  forceKa_local(atom_i,d) = forceKa_local(atom_i,d) + 0.25d0*c_a*sKa
                  forceKa_local(atom_j,d) = forceKa_local(atom_j,d) + 0.25d0*c_a*sJa
                  forceKb_local(atom_i,d) = forceKb_local(atom_i,d) + 0.25d0*c_a*sKb
                  forceKb_local(atom_j,d) = forceKb_local(atom_j,d) + 0.25d0*c_a*sJb
                   if (r3p_on .and. d .eq. r3p_dir) then
                      if (atom_i .eq. r3p_atom) r3p_an = r3p_an - (sKa+sKb)
                      if (atom_j .eq. r3p_atom) r3p_an = r3p_an - (sJa+sJb)
                   endif
                enddo
                if (hess_v2) then
                   do d = 1,3
                      do ee = 1,3
                         obpp  = nb*((i-1)+di*((j-1)+dj*((d-1)*3+(ee-1))))
                         obpp2 = nb*((j-1)+dj*((i-1)+di*((d-1)*3+(ee-1))))
                         obvp  = obpp
                         obvps = nb*((i-1)+di*((j-1)+dj*((ee-1)*3+(d-1))))
                         do g = 1,nb
                            ipp  = dbpp(obpp+g)
                            ivp  = dbvp(obvp+g)
                            ivps = dbvp(obvps+g)
                            ikk  = dbpp2(obpp2+g)
                            gpar = cosx_atom_hi(bstart+g-1)
                            t_h  = 0.25d0*c_a*tha(g)
                            ha_local(atom_i,d,atom_i,ee) = ha_local(atom_i,d,atom_i,ee) + t_h*ipp
                            ha_local(atom_i,d,atom_j,ee) = ha_local(atom_i,d,atom_j,ee) + t_h*ivp
                            ha_local(atom_i,d,gpar,ee)   = ha_local(atom_i,d,gpar,ee)   - t_h*(ipp+ivp)
                            ha_local(atom_j,d,atom_i,ee) = ha_local(atom_j,d,atom_i,ee) + t_h*ivps
                            ha_local(atom_j,d,atom_j,ee) = ha_local(atom_j,d,atom_j,ee) + t_h*ikk
                            ha_local(atom_j,d,gpar,ee)   = ha_local(atom_j,d,gpar,ee)   - t_h*(ivps+ikk)
                            ha_local(gpar,d,atom_i,ee)   = ha_local(gpar,d,atom_i,ee)   - t_h*(ipp+ivps)
                            ha_local(gpar,d,atom_j,ee)   = ha_local(gpar,d,atom_j,ee)   - t_h*(ivp+ikk)
                            ha_local(gpar,d,gpar,ee)     = ha_local(gpar,d,gpar,ee)     + t_h*(ipp+ivp+ivps+ikk)
                            t_h  = 0.25d0*c_a*thb(g)
                            hb_local(atom_i,d,atom_i,ee) = hb_local(atom_i,d,atom_i,ee) + t_h*ipp
                            hb_local(atom_i,d,atom_j,ee) = hb_local(atom_i,d,atom_j,ee) + t_h*ivp
                            hb_local(atom_i,d,gpar,ee)   = hb_local(atom_i,d,gpar,ee)   - t_h*(ipp+ivp)
                            hb_local(atom_j,d,atom_i,ee) = hb_local(atom_j,d,atom_i,ee) + t_h*ivps
                            hb_local(atom_j,d,atom_j,ee) = hb_local(atom_j,d,atom_j,ee) + t_h*ikk
                            hb_local(atom_j,d,gpar,ee)   = hb_local(atom_j,d,gpar,ee)   - t_h*(ivps+ikk)
                            hb_local(gpar,d,atom_i,ee)   = hb_local(gpar,d,atom_i,ee)   - t_h*(ipp+ivps)
                            hb_local(gpar,d,atom_j,ee)   = hb_local(gpar,d,atom_j,ee)   - t_h*(ivp+ikk)
                             hb_local(gpar,d,gpar,ee)     = hb_local(gpar,d,gpar,ee)     + t_h*(ipp+ivp+ivps+ikk)
                          enddo
                       enddo
                    enddo
                 endif
                 if (hess_v2 .and. .not. fit) then
                    if (same_spin) Fb = Fa
                    ordf = merge(2.0d0,1.0d0, sj .ne. si)
                    do dd2 = 1,3
                       do g = 1,nb
                          gpar = cosx_atom_hi(bstart+g-1)
                          sgij = sign(1.0d0,w_batch(g))*NorVEC(ao_i)*NorVEC(ao_j)
                          do d = 1,3
                             ob = nb*((i-1) + di*((j-1) + dj*(d-1)))
                             if (sj .ne. si) then
                                okk = nb*((j-1) + dj*((i-1) + di*(d-1)))
                             else
                                okk = nb*((j-1) + di*((i-1) + dj*(d-1)))
                             endif
                             dAbr(d) = dbuf(ob+g) ; dAkt(d) = dbuf2(okk+g)
                          enddo
                          do a2 = 1,2
                             if (a2 .eq. 1) then
                                Gfdi = 0.0d0 ; Gfdj = 0.0d0
                                SegAi(1:natoms) = 0.0d0 ; SegAj(1:natoms) = 0.0d0
                                do l2 = 1,nConts
                                   t_h = Xb1(l2,dd2,g)
                                   Gfdi = Gfdi + t_h*Da(l2,ao_i) ; Gfdj = Gfdj + t_h*Da(l2,ao_j)
                                   cc2 = ao_atom(l2)
                                   SegAi(cc2) = SegAi(cc2) + t_h*Da(l2,ao_i)
                                   SegAj(cc2) = SegAj(cc2) + t_h*Da(l2,ao_j)
                                enddo
                                dTh = ordf*sgij*(Gfdi*Fa(g,ao_j) + Fa(g,ao_i)*Gfdj)
                                do d = 1,3
                                   t_h = 0.25d0*c_a*dTh
                                   ha_local(atom_i,d,gpar,dd2) = ha_local(atom_i,d,gpar,dd2) + t_h*dAbr(d)
                                   ha_local(atom_j,d,gpar,dd2) = ha_local(atom_j,d,gpar,dd2) + t_h*dAkt(d)
                                   ha_local(gpar,d,gpar,dd2)   = ha_local(gpar,d,gpar,dd2)   - t_h*(dAbr(d)+dAkt(d))
                                enddo
                                do cc2 = 1,natoms
                                   dTh = -ordf*sgij*(SegAi(cc2)*Fa(g,ao_j) + Fa(g,ao_i)*SegAj(cc2))
                                   do d = 1,3
                                      t_h = 0.25d0*c_a*dTh
                                      ha_local(atom_i,d,cc2,dd2) = ha_local(atom_i,d,cc2,dd2) + t_h*dAbr(d)
                                      ha_local(atom_j,d,cc2,dd2) = ha_local(atom_j,d,cc2,dd2) + t_h*dAkt(d)
                                      ha_local(gpar,d,cc2,dd2)   = ha_local(gpar,d,cc2,dd2)   - t_h*(dAbr(d)+dAkt(d))
                                   enddo
                                enddo
                             else
                                Gfdi = 0.0d0 ; Gfdj = 0.0d0
                                SegAi(1:natoms) = 0.0d0 ; SegAj(1:natoms) = 0.0d0
                                do l2 = 1,nConts
                                   t_h = Xb1(l2,dd2,g)
                                   Gfdi = Gfdi + t_h*Db(l2,ao_i) ; Gfdj = Gfdj + t_h*Db(l2,ao_j)
                                   cc2 = ao_atom(l2)
                                   SegAi(cc2) = SegAi(cc2) + t_h*Db(l2,ao_i)
                                   SegAj(cc2) = SegAj(cc2) + t_h*Db(l2,ao_j)
                                enddo
                                dTh = ordf*sgij*(Gfdi*Fb(g,ao_j) + Fb(g,ao_i)*Gfdj)
                                do d = 1,3
                                   t_h = 0.25d0*c_a*dTh
                                   hb_local(atom_i,d,gpar,dd2) = hb_local(atom_i,d,gpar,dd2) + t_h*dAbr(d)
                                   hb_local(atom_j,d,gpar,dd2) = hb_local(atom_j,d,gpar,dd2) + t_h*dAkt(d)
                                   hb_local(gpar,d,gpar,dd2)   = hb_local(gpar,d,gpar,dd2)   - t_h*(dAbr(d)+dAkt(d))
                                enddo
                                 do cc2 = 1,natoms
                                    dTh = -ordf*sgij*(SegAi(cc2)*Fb(g,ao_j) + Fb(g,ao_i)*SegAj(cc2))
                                    do d = 1,3
                                       t_h = 0.25d0*c_a*dTh
                                       hb_local(atom_i,d,cc2,dd2) = hb_local(atom_i,d,cc2,dd2) + t_h*dAbr(d)
                                       hb_local(atom_j,d,cc2,dd2) = hb_local(atom_j,d,cc2,dd2) + t_h*dAkt(d)
                                       hb_local(gpar,d,cc2,dd2)   = hb_local(gpar,d,cc2,dd2)   - t_h*(dAbr(d)+dAkt(d))
                                    enddo
                                 enddo
                              endif
                           enddo
                        enddo
                     enddo
                  endif
                  if (hess_v2 .and. .not. fit) then
                     if (same_spin) Fb = Fa
                     do dd2 = 1,3
                        do g = 1,nb
                           gpar = cosx_atom_hi(bstart+g-1)
                           ob = nb*((i-1)+di*((j-1)+dj*(dd2-1)))
                           if (sj .ne. si) then
                              okk = nb*((j-1)+dj*((i-1)+di*(dd2-1)))
                           else
                              okk = nb*((j-1)+di*((i-1)+dj*(dd2-1)))
                           endif
                           dAbr_dd = dbuf(ob+g) ; dAkt_dd = dbuf2(okk+g)
                           sgnij_v = sign(1.0d0,w_batch(g))*NorVEC(ao_i)*NorVEC(ao_j)
                           wij_v = sgnij_v*buf(g+nb*((i-1)+di*(j-1)))
                           do a2 = 1,2
                              if (a2 .eq. 1) then
                                 Gfdi=0.0d0;Gfdj=0.0d0;SegAi=0.0d0;SegAj=0.0d0
                                 do l2=1,nConts
                                    t_h=Xb1(l2,dd2,g)
                                    Gfdi=Gfdi+t_h*Da(l2,ao_i); Gfdj=Gfdj+t_h*Da(l2,ao_j)
                                    cc2=ao_atom(l2)
                                    SegAi(cc2)=SegAi(cc2)+t_h*Da(l2,ao_i)
                                    SegAj(cc2)=SegAj(cc2)+t_h*Da(l2,ao_j)
                                 enddo
                              else
                                 Gfdi=0.0d0;Gfdj=0.0d0;SegAi=0.0d0;SegAj=0.0d0
                                 do l2=1,nConts
                                    t_h=Xb1(l2,dd2,g)
                                    Gfdi=Gfdi+t_h*Db(l2,ao_i); Gfdj=Gfdj+t_h*Db(l2,ao_j)
                                    cc2=ao_atom(l2)
                                    SegAi(cc2)=SegAi(cc2)+t_h*Db(l2,ao_i)
                                    SegAj(cc2)=SegAj(cc2)+t_h*Db(l2,ao_j)
                                 enddo
                              endif
                              Wmat_i=0.0d0;Wmat_j=0.0d0;Gx_i=0.0d0;Gx_j=0.0d0
                              do sjj2=cosx_blkptr_hi(ib)+1,cosx_blkptr_hi(ib+1)
                                 sk2=cosx_blkidx_hi(sjj2)
                                 do kk2=1,cosx_shell_dim(sk2)
                                    kappa2=cosx_ao_offset(sk2)+kk2; cc2=cosx_shell_atom(sk2)
                                    do d=1,3
                                       t_h=Xb1(kappa2,d,g)
                                       if (a2 .eq. 1) then
                                          Wmat_i(cc2,d)=Wmat_i(cc2,d)+t_h*Da(ao_i,kappa2)
                                          Wmat_j(cc2,d)=Wmat_j(cc2,d)+t_h*Da(ao_j,kappa2)
                                          Gx_i(d)=Gx_i(d)+t_h*Da(ao_i,kappa2)
                                          Gx_j(d)=Gx_j(d)+t_h*Da(ao_j,kappa2)
                                       else
                                          Wmat_i(cc2,d)=Wmat_i(cc2,d)+t_h*Db(ao_i,kappa2)
                                          Wmat_j(cc2,d)=Wmat_j(cc2,d)+t_h*Db(ao_j,kappa2)
                                          Gx_i(d)=Gx_i(d)+t_h*Db(ao_i,kappa2)
                                          Gx_j(d)=Gx_j(d)+t_h*Db(ao_j,kappa2)
                                       endif
                                    enddo
                                 enddo
                              enddo
                              do cc2=1,natoms
                                 dAijC=0.0d0
                                 if (cc2.eq.atom_i) dAijC=dAijC+dAbr_dd
                                 if (cc2.eq.atom_j) dAijC=dAijC+dAkt_dd
                                 if (cc2.eq.gpar)   dAijC=dAijC-(dAbr_dd+dAkt_dd)
                                 dFaiC=-SegAi(cc2) ; dFajC=-SegAj(cc2)
                                 if (cc2.eq.gpar) then
                                    dFaiC=dFaiC+Gfdi ; dFajC=dFajC+Gfdj
                                 endif
                                 dGai=0.0d0;dGaj=0.0d0;dGtai=0.0d0;dGtaj=0.0d0
                                 if (a2 .eq. 1) then
                                    dGai = sgnij_v*(dAijC*Fa(g,ao_j) + wij_v*dFajC)
                                    dGtaj= sgnij_v*(dAijC*Fa(g,ao_i) + wij_v*dFaiC)
                                    if (sj .ne. si) then
                                       dGaj = sgnij_v*(dAijC*Fa(g,ao_i) + wij_v*dFaiC)
                                       dGtai= sgnij_v*(dAijC*Fa(g,ao_j) + wij_v*dFajC)
                                    endif
                                 else
                                    dGai = sgnij_v*(dAijC*Fb(g,ao_j) + wij_v*dFajC)
                                    dGtaj= sgnij_v*(dAijC*Fb(g,ao_i) + wij_v*dFaiC)
                                    if (sj .ne. si) then
                                       dGaj = sgnij_v*(dAijC*Fb(g,ao_i) + wij_v*dFaiC)
                                       dGtai= sgnij_v*(dAijC*Fb(g,ao_j) + wij_v*dFajC)
                                    endif
                                 endif
                                 do cc3=1,natoms
                                    do d=1,3
                                       t_h = 0.25d0*c_x*((dGai+dGtai)*Wmat_i(cc3,d) + (dGaj+dGtaj)*Wmat_j(cc3,d))
                                       if (cc3.eq.gpar) t_h = t_h &
                                            - 0.25d0*c_x*((dGai+dGtai)*Gx_i(d) + (dGaj+dGtaj)*Gx_j(d))
                                       if (a2 .eq. 1) then
                                          ha_local(cc3,d,cc2,dd2) = ha_local(cc3,d,cc2,dd2) + t_h
                                       else
                                          hb_local(cc3,d,cc2,dd2) = hb_local(cc3,d,cc2,dd2) + t_h
                                       endif
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo

         if (r3p_on) then
            envp1 = env_g ; envp2 = env_g
            envp1(atm(2,r3p_atom)+r3p_dir) = envp1(atm(2,r3p_atom)+r3p_dir) + r3p_h
            envp2(atm(2,r3p_atom)+r3p_dir) = envp2(atm(2,r3p_atom)+r3p_dir) - r3p_h
            call grids1e_engine_cached(bufp1, shls1e, atm, size(atm,2), bas, nBases, envp1, 0_8, int_cache)
            call grids1e_engine_cached(bufp2, shls1e, atm, size(atm,2), bas, nBases, envp2, 0_8, int_cache)
            do j = 1,dj
               ao_j = cosx_ao_offset(sj) + j
               do i = 1,di
                  ao_i = cosx_ao_offset(si) + i
                  nij = NorVEC(ao_i)*NorVEC(ao_j)
                  if (nij .eq. 0.0d0) cycle
                  do g = 1,nb
                     sgw = sign(1.0d0,w_batch(g))
                     theta = sgw * nij * (Pha(g,ao_i)*Fa(g,ao_j) + Phb(g,ao_i)*Fb(g,ao_j))
                     if (sj .ne. si) theta = theta + sgw * nij * &
                          (Pha(g,ao_j)*Fa(g,ao_i) + Phb(g,ao_j)*Fb(g,ao_i))
                     r3p_num = r3p_num + theta * &
                        (bufp1(g + nb*((i-1)+di*(j-1))) - bufp2(g + nb*((i-1)+di*(j-1)))) / (2.0d0*r3p_h)
                  enddo
               enddo
            enddo
         endif
         if (tmg) tm_r3 = tm_r3 + omp_get_wtime() - tm0
      enddo
   enddo

   if (tmg) tm0 = omp_get_wtime()
   allocate(Ha(nb,nConts), Hb(nb,nConts), Hta(nb,nConts), Htb(nb,nConts))
   allocate(Psia(nb,nConts), Psib(nb,nConts))
   if (fit) then
      call dgemm('N','N',nb,nConts,nConts,1.0d0,Ga,nb,DPa,nConts,0.0d0,Psia,nb)
      call dgemm('N','N',nb,nConts,nConts,1.0d0,Gta,nb,Da,nConts,1.0d0,Psia,nb)
   else
      Hta = Ga + Gta
      call dgemm('N','N',nb,nConts,nConts,1.0d0,Hta,nb,Da,nConts,0.0d0,Psia,nb)
   endif
   if (psi_chk .and. .not. fit) &
      call dgemm('N','N',nb,nConts,nConts,1.0d0,Ga,nb,Da,nConts,0.0d0,Ha,nb)
   if (same_spin) then
      Psib = Psia
   else
      if (fit) then
         call dgemm('N','N',nb,nConts,nConts,1.0d0,Gb,nb,DPb,nConts,0.0d0,Psib,nb)
         call dgemm('N','N',nb,nConts,nConts,1.0d0,Gtb,nb,Db,nConts,1.0d0,Psib,nb)
      else
         Htb = Gb + Gtb
         call dgemm('N','N',nb,nConts,nConts,1.0d0,Htb,nb,Db,nConts,0.0d0,Psib,nb)
      endif
   endif

   if (psi_chk .and. .not. fit) then
      psi_dev = max(psi_dev, maxval(abs(Psia - 2.0d0*Ha)))
      psi_scale = max(psi_scale, maxval(abs(Psia)))
   endif

   if (tmg) tm_gemm = tm_gemm + omp_get_wtime() - tm0
   if (tmg) tm0 = omp_get_wtime()
   do sj_idx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      si = cosx_blkidx_hi(sj_idx)
      atom_i = cosx_shell_atom(si)
      di = cosx_shell_dim(si)
      do d = 1,3
         do i = 1,di
            kappa = cosx_ao_offset(si) + i
            sKa = 0.0d0
            sKb = 0.0d0
            do g = 1,nb
               tha(g) = Xb1(kappa,d,g)*Psia(g,kappa)
               sKa = sKa + tha(g)
               pt_a(g,d) = pt_a(g,d) + c_x*tha(g)
            enddo
            if (same_spin) then
               sKb = sKa
               do g = 1,nb
                  pt_b(g,d) = pt_b(g,d) + c_x*tha(g)
               enddo
            else
               do g = 1,nb
                  thb(g) = Xb1(kappa,d,g)*Psib(g,kappa)
                  sKb = sKb + thb(g)
                  pt_b(g,d) = pt_b(g,d) + c_x*thb(g)
               enddo
            endif
            forceKa_local(atom_i,d) = forceKa_local(atom_i,d) + 0.25d0*c_x*sKa
            forceKb_local(atom_i,d) = forceKb_local(atom_i,d) + 0.25d0*c_x*sKb
            if (hess_v2) then
               do ee = 1,3
                  do g = 1,nb
                     gpar = cosx_atom_hi(bstart+g-1)
                     t_h = 0.25d0*c_x*Xb2(kappa,d,ee,g)
                     ha_local(atom_i,d,atom_i,ee) = ha_local(atom_i,d,atom_i,ee) - t_h*Psia(g,kappa)
                     ha_local(atom_i,d,gpar,ee)   = ha_local(atom_i,d,gpar,ee)   + t_h*Psia(g,kappa)
                     ha_local(gpar,d,atom_i,ee)   = ha_local(gpar,d,atom_i,ee)   + t_h*Psia(g,kappa)
                     ha_local(gpar,d,gpar,ee)     = ha_local(gpar,d,gpar,ee)     - t_h*Psia(g,kappa)
                     hb_local(atom_i,d,atom_i,ee) = hb_local(atom_i,d,atom_i,ee) - t_h*Psib(g,kappa)
                     hb_local(atom_i,d,gpar,ee)   = hb_local(atom_i,d,gpar,ee)   + t_h*Psib(g,kappa)
                     hb_local(gpar,d,atom_i,ee)   = hb_local(gpar,d,atom_i,ee)   + t_h*Psib(g,kappa)
                     hb_local(gpar,d,gpar,ee)     = hb_local(gpar,d,gpar,ee)     - t_h*Psib(g,kappa)
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo

   if (.not. frozen) then
   do d = 1,3
      do g = 1,nb
         gpar = cosx_atom_hi(bstart+g-1)
         forceKa_local(gpar,d) = forceKa_local(gpar,d) - 0.25d0*pt_a(g,d)
         forceKb_local(gpar,d) = forceKb_local(gpar,d) - 0.25d0*pt_b(g,d)
      enddo
   enddo
   endif

   if (tmg) tm_r2 = tm_r2 + omp_get_wtime() - tm0
   if (tmg) tm0 = omp_get_wtime()
   do g = 1,nb
      gg = bstart + g - 1
      eps_a(gg) = 0.5d0*dot_product(Xc(1:nnz,g), Psia(g,aoidx)) / &
                  sign(max(abs(w_batch(g)),1.0d-300), w_batch(g))
      if (same_spin) then
         eps_b(gg) = eps_a(gg)
      else
         eps_b(gg) = 0.5d0*dot_product(Xc(1:nnz,g), Psib(g,aoidx)) / &
                     sign(max(abs(w_batch(g)),1.0d-300), w_batch(g))
      endif
   enddo
   if (fit) then
      allocate(Tc(nnz,nConts))
      call dgemm('N','N',nnz,nConts,nb,1.0d0,Xc,nnz,Ga,nb,0.0d0,Tc,nnz)
      do j = 1,nConts
         Va_local(aoidx,j) = Va_local(aoidx,j) + Tc(:,j)
      enddo
      if (.not. same_spin) then
         call dgemm('N','N',nnz,nConts,nb,1.0d0,Xc,nnz,Gb,nb,0.0d0,Tc,nnz)
         do j = 1,nConts
            Vb_local(aoidx,j) = Vb_local(aoidx,j) + Tc(:,j)
         enddo
      endif
      deallocate(Tc)
   endif

   if (tmg) tm_eps = tm_eps + omp_get_wtime() - tm0
   if (hess_v2) deallocate(Xb2)
   deallocate(coor_batch,w_batch,Xb,Xb1,Fa,Fb,Pha,Phb,F_gmax,P_gmax)
   deallocate(aoidx,Xc,Dc)
   deallocate(Ga,Gb,Gta,Gtb,Ha,Hb,Hta,Htb,Psia,Psib)
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
         if (fit) then
            Va = Va + Va_local
            Vb = Vb + Vb_local
         endif
         if (hess_v2) then
            ha_glob = ha_glob + ha_local
            hb_glob = hb_glob + hb_local
         endif
      endif
   enddo
end block
if (fit) deallocate(Va_local,Vb_local)
deallocate(pt_a,pt_b,tha,thb)
if (r3p_on) deallocate(envp1,envp2,bufp1,bufp2)
deallocate(buf,dbuf,dbuf2,int_cache,env_g,forceKa_local,forceKb_local)
if (grad_fused) deallocate(vbuf)
if (grad_contract) deallocate(cvi,cvj,coi,coj,cot,csg)
if (hess_v2) deallocate(dbpp,dbvp,dbpp2,ha_local,hb_local,SegAi,SegAj,SegBi,SegBj,Wmat_i,Wmat_j)
!$omp end parallel
tw_loop = omp_get_wtime() - tw0
if (fit) then
   allocate(Tc(nConts,nConts))
   call dgemm('N','N',nConts,nConts,nConts,1.0d0,Va,nConts,Da,nConts,0.0d0,Tc,nConts)
   Va = Tc
   if (same_spin) then
      Vb = Va
   else
      call dgemm('N','N',nConts,nConts,nConts,1.0d0,Vb,nConts,Db,nConts,0.0d0,Tc,nConts)
      Vb = Tc
   endif
   deallocate(Tc, Ma, DPa)
   if (allocated(Mb)) deallocate(Mb, DPb)
endif
if (hess_v2) then
   if (.not. allocated(cosx_hess_a_save_a)) then
      allocate(cosx_hess_a_save_a(3,natoms,3,natoms), cosx_hess_a_save_b(3,natoms,3,natoms))
   endif
   cosx_hess_a_save_a = ha_glob ; cosx_hess_a_save_b = hb_glob
   cosx_hess_a_saved = .true.
   print '(A,ES12.4)', " COSX HESS A-TERM: max|d2E_a+d2E_b| = ", maxval(abs(ha_glob+hb_glob))
   print '(A,3ES12.4)', " COSX HESS A-TERM: transl sum (atom1) = ", &
         sum(ha_glob(1,1,:,:)+hb_glob(1,1,:,:)), sum(ha_glob(2,1,:,:)+hb_glob(2,1,:,:)), &
         sum(ha_glob(3,1,:,:)+hb_glob(3,1,:,:))
   deallocate(ha_glob,hb_glob)
endif
if (r3p_on) then
   print '(A,I3,A,I2,A)', " COSX V2 R3PROBE (atom",r3p_atom,", dir",r3p_dir,"):"
   print '(A,ES22.14)', "   analytic nuclear-role dE/dR = ", r3p_an
   print '(A,ES22.14)', "   numerical (same Theta/grid) = ", r3p_num
   print '(A,ES12.4)',  "   ratio an/num                = ", r3p_an/max(abs(r3p_num),1.0d-300)*sign(1.0d0,r3p_num)
   call flush(6)
endif

if (psi_chk) then
   print '(A,ES12.4,A,ES12.4)', " COSX V2 SCREEN-DIAG: max(max|Pha|/max|Fa|) =", psi_dev, &
        "   max|Qfull| =", psi_scale
   print '(A,I12,A,I12,A,F8.4)', " COSX V2 SCREEN-DIAG: pairs offered =", n_offer, &
        "  kept =", n_keep, "  frac =", dble(n_keep)/max(dble(n_offer),1.0d0)
   call flush(6)
endif
if (psi_chk .and. .not. fit) then
   print '(A,2ES14.6,A,ES12.4)', " COSX V2 PSI-CHECK: max|Psi-2H|=", psi_dev, &
        psi_scale, "  rel=", psi_dev/max(psi_scale,1.0d-300)
   call flush(6)
endif
block
   character(len=8) :: echk_env
   real(8) :: e_from_eps_a, e_from_eps_b
   call get_environment_variable("ENGINE_COSX_V2_ECHK", echk_env)
   if (trim(echk_env) .eq. "1") then
      e_from_eps_a = 0.0d0
      e_from_eps_b = 0.0d0
      do gg = 1,npts
         e_from_eps_a = e_from_eps_a + cosx_w_hi(gg)*eps_a(gg)
         e_from_eps_b = e_from_eps_b + cosx_w_hi(gg)*eps_b(gg)
      enddo
      print '(A,2ES24.16)', " COSX V2 ECHK: sum_g w_g eps_g (a,b) =", &
            e_from_eps_a, e_from_eps_b
      print '(A,ES24.16)',  " COSX V2 ECHK: -0.5*(a+b)            =", &
            -0.5d0*(e_from_eps_a+e_from_eps_b)
      call flush(6)
   endif
end block

block
   character(len=8) :: dmp_env
   integer :: udmp
   call get_environment_variable("ENGINE_COSX_V2_DUMP", dmp_env)
   if (trim(dmp_env) .eq. "1") then
      open(newunit=udmp,file='direwolf_cosx_v2dump.txt',status='replace')
      write(udmp,*) npts, natoms, nConts
      do gg = 1,npts
         write(udmp,'(3ES24.16,I8,3ES24.16)') cosx_coor_hi(1,gg), cosx_coor_hi(2,gg), &
              cosx_coor_hi(3,gg), cosx_atom_hi(gg), cosx_w_hi(gg), eps_a(gg), eps_b(gg)
      enddo
      do kq = 1,nConts
         write(udmp,'(*(ES24.16))') Da(kq,1:nConts)
      enddo
      do kq = 1,nConts
         write(udmp,'(*(ES24.16))') Db(kq,1:nConts)
      enddo
      do kq = 1,nConts
         write(udmp,'(I8,ES24.16)') cosx_shell_atom(0) , NorVEC(kq)
      enddo
      close(udmp)
   endif
end block

tw0 = omp_get_wtime()
allocate(wt_a(natoms,3), wt_b(natoms,3), zero_eps(npts))
zero_eps = 0.0d0
call vv10_nlc_weight_gradient(npts, cosx_coor_hi, cosx_atom_hi, &
     max(npts/max(natoms,1),1), eps_a, zero_eps, cosx_w_hi, -1.0d0, wt_a)
if (same_spin) then
   wt_b = wt_a
else
   call vv10_nlc_weight_gradient(npts, cosx_coor_hi, cosx_atom_hi, &
        max(npts/max(natoms,1),1), eps_b, zero_eps, cosx_w_hi, -1.0d0, wt_b)
endif
forceKa = forceKa - 0.25d0*c_w*wt_a
forceKb = forceKb - 0.25d0*c_w*wt_b
if (.not. allocated(cosx_eps_a_hi_save)) &
   allocate(cosx_eps_a_hi_save(npts), cosx_eps_b_hi_save(npts))
cosx_eps_a_hi_save = eps_a
cosx_eps_b_hi_save = eps_b
cosx_eps_saved = .true.
deallocate(wt_a, wt_b, zero_eps, eps_a, eps_b, ao_atom)

tw_r1 = omp_get_wtime() - tw0
if (tmg) tm0 = omp_get_wtime()
if (fit .and. cosx_Snum_lu_ready) then
   allocate(dS1e(nConts,nConts,3), dH1e(nConts,nConts,3))
   call integrals_compute_force_1e(nConts, dS1e, dH1e)
   deallocate(dH1e)
   allocate(ao_atom(nConts))
   do si = 0,nBases-1
      do i = 1,cosx_shell_dim(si)
         ao_atom(cosx_ao_offset(si)+i) = cosx_shell_atom(si)
      enddo
   enddo

   allocate(Yq(nConts,nConts), Wq(nConts,nConts))
   allocate(termA(natoms,3), termB(natoms,3))
   allocate(uq(npts), wtq(natoms,3), wq_zero(npts))
   wq_zero = 0.0d0

   do ispin = 1,2
      if (ispin .eq. 1) then
         Yq = Va
      else
         if (same_spin) then
            forceKb = forceKb + (forceKa - fk_save)
            exit
         endif
         Yq = Vb
      endif
      if (ispin .eq. 1) fk_save = forceKa
      call dgetrs('T', nConts, nConts, cosx_Snum_lu_hi, nConts, &
                  cosx_ipiv_hi, Yq, nConts, info_q)
      if (info_q .ne. 0) then
         print *, "COSX force v2: dgetrs failed for the fit response (info=",info_q,")"
         call flush(6)
         exit
      endif
      call dgemm('N','T',nConts,nConts,nConts,1.0d0,Yq,nConts, &
                 cosx_Qfull_hi,nConts,0.0d0,Wq,nConts)

      termA = 0.0d0
      do dirq = 1,3
         do lq = 1,nConts
            do kq = 1,nConts
               termA(ao_atom(kq),dirq) = termA(ao_atom(kq),dirq) - dS1e(kq,lq,dirq)*Yq(kq,lq)
               termA(ao_atom(lq),dirq) = termA(ao_atom(lq),dirq) - dS1e(lq,kq,dirq)*Yq(kq,lq)
            enddo
         enddo
      enddo

      termB = 0.0d0
      uq = 0.0d0
      !$omp parallel default(shared) &
      !$omp&   private(ib,bstart,nb,gq,g,dirq,kq,sj_idx,si,Xr,Xr1,WX,WTX,tval,sKa) &
      !$omp&   reduction(+:termB)
      allocate(WX(nConts,BATCH), WTX(nConts,BATCH))
      !$omp do schedule(dynamic)
      do ib = 0,nblocks-1
         bstart = cosx_blkstart_hi(ib) + 1
         nb = cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib)
         allocate(Xr(nConts,nb), Xr1(nConts,3,nb))
         Xr = 0.0d0
         Xr1 = 0.0d0
         do sj_idx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
            si = cosx_blkidx_hi(sj_idx)
            call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                                      cosx_coor_hi(:,bstart:bstart+nb-1), nb, &
                                      cosx_ao_offset(si), nConts, Xr, Xr1)
         enddo
         call dgemm('N','N',nConts,nb,nConts,1.0d0,Wq,nConts,Xr,nConts,0.0d0,WX,nConts)
         call dgemm('T','N',nConts,nb,nConts,1.0d0,Wq,nConts,Xr,nConts,0.0d0,WTX,nConts)
         do gq = 1,nb
            g = bstart + gq - 1
            uq(g) = dot_product(Xr(1:nConts,gq), WX(1:nConts,gq))
            do dirq = 1,3
               sKa = 0.0d0
               do kq = 1,nConts
                  tval = cosx_w_hi(g)*Xr1(kq,dirq,gq)*(WX(kq,gq) + WTX(kq,gq))
                  sKa = sKa + tval
                  termB(ao_atom(kq),dirq) = termB(ao_atom(kq),dirq) - tval
               enddo
               if (.not. frozen) termB(cosx_atom_hi(g),dirq) = termB(cosx_atom_hi(g),dirq) + sKa
            enddo
         enddo
         deallocate(Xr,Xr1)
      enddo
      !$omp end do
      deallocate(WX,WTX)
      !$omp end parallel

      call vv10_nlc_weight_gradient(npts, cosx_coor_hi, cosx_atom_hi, &
           max(npts/max(natoms,1),1), uq, wq_zero, cosx_w_hi, -1.0d0, wtq)
      if (.not. frozen) termB = termB + wtq

      if (ispin .eq. 1) then
         forceKa = forceKa - 0.25d0*c_q*(termA - termB)
      else
         forceKb = forceKb - 0.25d0*c_q*(termA - termB)
      endif
   enddo

   deallocate(Yq, Wq, termA, termB, uq, wtq, wq_zero, dS1e, ao_atom)
endif

if (fit) deallocate(Va, Vb)
if (tmg) then
   tm_r4 = omp_get_wtime() - tm0
   print '(A)', " COSX V2 TIMING (thread-summed seconds):"
   print '(A,F10.3,A,F10.3,A,F10.3)', "   AO eval  ", tm_ao,  "   dgemm   ", tm_gemm, "   ESP val ", tm_espv
   print '(A,F10.3,A,F10.3,A,F10.3)', "   ESP deriv", tm_espd,"   G/Gt acc", tm_gacc, "   R3      ", tm_r3
   print '(A,F10.3,A,F10.3,A,F10.3)', "   R2       ", tm_r2,  "   eps/V   ", tm_eps,  "   R4(wall)", tm_r4
   print '(A,F10.3,A,F10.3)', "   loop(wall)", tw_loop, "   R1(wall)", tw_r1
   call flush(6)
endif

call cosx_malloc_trim()

end subroutine cosx_force_exchange_v2

module subroutine cosx_hess_weight_term(natoms_in, d2E)
use MOL_info, only: natoms
implicit none
integer,intent(in) :: natoms_in
real(8),intent(out) :: d2E(3,natoms,3,natoms)
real(8),allocatable :: h1(:,:,:,:)
integer :: nn, i, j, dd, ee
d2E = 0.0d0
if (.not. cosx_eps_saved) return
if (natoms_in .ne. natoms) return
nn = size(cosx_eps_a_hi_save)
allocate(h1(natoms,3,natoms,3))
call vv10_nlc_weight_hessian(nn, cosx_coor_hi, cosx_atom_hi, cosx_w_hi, cosx_eps_a_hi_save, h1)
do i = 1,natoms ; do dd = 1,3 ; do j = 1,natoms ; do ee = 1,3
   d2E(dd,i,ee,j) = d2E(dd,i,ee,j) - 0.25d0*h1(i,dd,j,ee)
enddo ; enddo ; enddo ; enddo
call vv10_nlc_weight_hessian(nn, cosx_coor_hi, cosx_atom_hi, cosx_w_hi, cosx_eps_b_hi_save, h1)
do i = 1,natoms ; do dd = 1,3 ; do j = 1,natoms ; do ee = 1,3
   d2E(dd,i,ee,j) = d2E(dd,i,ee,j) - 0.25d0*h1(i,dd,j,ee)
enddo ; enddo ; enddo ; enddo
deallocate(h1)
end subroutine cosx_hess_weight_term

module subroutine cosx_hess_a_term(natoms_in, d2E_a, d2E_b)
use MOL_info, only: natoms
implicit none
integer,intent(in) :: natoms_in
real(8),intent(out) :: d2E_a(3,natoms,3,natoms), d2E_b(3,natoms,3,natoms)
d2E_a = 0.0d0 ; d2E_b = 0.0d0
if (.not. cosx_hess_a_saved) return
if (natoms_in .ne. natoms) return
d2E_a = cosx_hess_a_save_a
d2E_b = cosx_hess_a_save_b
end subroutine cosx_hess_a_term

subroutine cosx_grad_convention_probe(nConts, RSomega)
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: RSomega
real(8),parameter :: HPROBE = 1.0d-5
real(8),allocatable :: v0p(:,:), v0m(:,:), v1(:,:,:), v0(:,:)
real(8),allocatable :: envp(:), abuf(:), abufp(:), abufm(:), dbra(:), dket(:)
real(8),allocatable :: icache(:)
real(8) :: rg(3), cpt(3,1), fd, an, num, den
integer :: si, sj, di, dj, shp(4), offg, d, i, j, idx

if (.not. allocated(cosx_coor_hi)) return
si = 0
sj = min(1, nBases-1)
di = cosx_shell_dim(si)
dj = cosx_shell_dim(sj)
block
   real(8),allocatable :: q0(:,:), q0p(:,:), q0m(:,:), q1(:,:,:)
   real(8) :: qpt(3,1), qfd, qan, wr, worst, rat
   integer :: ss, dd, ii, jj, widx
   allocate(q0(nConts,1), q0p(nConts,1), q0m(nConts,1), q1(nConts,3,1))
   worst = 0.0d0
   widx = -1
   do ss = 0,nBases-1
      qpt(:,1) = cosx_shell_cen(:,ss) + 0.3d0
      do dd = 1,3
         q0 = 0.0d0 ; q1 = 0.0d0
         call GTOeval_shell_batch(cosx_shell_atom(ss), cosx_shell_local(ss), &
                                   qpt, 1, cosx_ao_offset(ss), nConts, q0, q1)
         q0p = 0.0d0
         qpt(dd,1) = qpt(dd,1) + HPROBE
         call GTOeval_shell_batch(cosx_shell_atom(ss), cosx_shell_local(ss), &
                                   qpt, 1, cosx_ao_offset(ss), nConts, q0p, q1)
         qpt(dd,1) = qpt(dd,1) - 2.0d0*HPROBE
         q0m = 0.0d0
         call GTOeval_shell_batch(cosx_shell_atom(ss), cosx_shell_local(ss), &
                                   qpt, 1, cosx_ao_offset(ss), nConts, q0m, q1)
         qpt(dd,1) = qpt(dd,1) + HPROBE
         q0 = 0.0d0 ; q1 = 0.0d0
         call GTOeval_shell_batch(cosx_shell_atom(ss), cosx_shell_local(ss), &
                                   qpt, 1, cosx_ao_offset(ss), nConts, q0, q1)
         do ii = 1,cosx_shell_dim(ss)
            jj = cosx_ao_offset(ss) + ii
            qfd = (q0p(jj,1) - q0m(jj,1)) / (2.0d0*HPROBE)
            qan = q1(jj,dd,1)
            if (abs(qfd) .lt. 1.0d-8) cycle
            rat = qan/qfd
            if (abs(rat - 1.0d0) .gt. worst) then
               worst = abs(rat - 1.0d0)
               widx = ss*100 + ii*10 + dd
            endif
         enddo
      enddo
   enddo
   print '(A,ES12.4,A,I6)', " COSX CONV PROBE: worst |val1/grad_r(phi) - 1| over ALL shells =", &
         worst, "   at shell*100+ao*10+dir =", widx
   deallocate(q0,q0p,q0m,q1)
end block
rg = cosx_shell_cen(:,si) + 0.3d0

allocate(v0(nConts,1), v0p(nConts,1), v0m(nConts,1), v1(nConts,3,1))
num = 0.0d0
den = 0.0d0
do d = 1,3
   v0 = 0.0d0 ; v1 = 0.0d0
   cpt(:,1) = rg
   call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                             cpt, 1, cosx_ao_offset(si), nConts, v0, v1)
   v0p = 0.0d0
   cpt(:,1) = rg ; cpt(d,1) = rg(d) + HPROBE
   call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                             cpt, 1, cosx_ao_offset(si), nConts, v0p, v1)
   v0m = 0.0d0
   cpt(:,1) = rg ; cpt(d,1) = rg(d) - HPROBE
   call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                             cpt, 1, cosx_ao_offset(si), nConts, v0m, v1)
   v1 = 0.0d0 ; v0 = 0.0d0
   cpt(:,1) = rg
   call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                             cpt, 1, cosx_ao_offset(si), nConts, v0, v1)
   do i = 1,di
      idx = cosx_ao_offset(si) + i
      fd = (v0p(idx,1) - v0m(idx,1)) / (2.0d0*HPROBE)
      an = v1(idx,d,1)
      num = num + an*fd
      den = den + fd*fd
   enddo
enddo
print *, "COSX CONV PROBE: val1_out / grad_r(phi) =", num/max(den,1.0d-300), &
         "  (+1 => val1 IS grad_r phi => dX/dR_center = -val1)"

allocate(envp(size(env) + 3))
envp(1:size(env)) = env(1:size(env))
offg = size(env)
envp(9) = RSomega
envp(12) = 1.0d0
envp(13) = dble(offg)
allocate(abuf(di*dj), abufp(di*dj), abufm(di*dj))
allocate(dbra(di*dj*3), dket(di*dj*3))
allocate(icache(200000))
shp(1) = si ; shp(2) = sj ; shp(3) = 0 ; shp(4) = 1
envp(offg+1:offg+3) = rg
call grids1e_ip_engine(dbra, shp, atm, size(atm,2), bas, nBases, envp, 0_8)
shp(1) = sj ; shp(2) = si
call grids1e_ip_engine(dket, shp, atm, size(atm,2), bas, nBases, envp, 0_8)
shp(1) = si ; shp(2) = sj
num = 0.0d0
den = 0.0d0
do d = 1,3
   envp(offg+1:offg+3) = rg ; envp(offg+d) = rg(d) + HPROBE
   call grids1e_engine_cached(abufp, shp, atm, size(atm,2), bas, nBases, envp, 0_8, icache)
   envp(offg+1:offg+3) = rg ; envp(offg+d) = rg(d) - HPROBE
   call grids1e_engine_cached(abufm, shp, atm, size(atm,2), bas, nBases, envp, 0_8, icache)
   do j = 1,dj
      do i = 1,di
         fd = (abufp((i-1)+di*(j-1)+1) - abufm((i-1)+di*(j-1)+1)) / (2.0d0*HPROBE)
         an = dbra((i-1) + di*((j-1) + dj*(d-1)) + 1) &
            + dket((j-1) + dj*((i-1) + di*(d-1)) + 1)
         num = num + an*fd
         den = den + fd*fd
      enddo
   enddo
enddo
print *, "COSX CONV PROBE: (ip_bra+ip_ket) / grad_rg(A) =", num/max(den,1.0d-300), &
         "  (-1 => ip IS dA/dR_center; +1 => ip is MINUS it)"

block
   use mod_xc_hess, only: eval_bf_derivs
   real(8),allocatable :: ev0(:), ev1(:,:), ev2(:,:,:), ev3(:,:,:,:), gv0(:), gv1(:,:,:)
   real(8) :: rp(3,1), dev, scl
   integer :: mm, dd2, ss2
   allocate(ev0(nConts), ev1(nConts,3), ev2(nConts,3,3), ev3(nConts,3,3,3))
   allocate(gv0(nConts), gv1(nConts,3,1))
   rp(:,1) = cosx_shell_cen(:,si) + 0.3d0
   call eval_bf_derivs(rp(:,1), ev0, ev1, ev2, ev3, 2)
   gv0 = 0.0d0 ; gv1 = 0.0d0
   do ss2 = 0,nBases-1
      call GTOeval_shell_batch(cosx_shell_atom(ss2), cosx_shell_local(ss2), &
                               rp, 1, cosx_ao_offset(ss2), nConts, gv0, gv1)
   enddo
   dev = 0.0d0 ; scl = 0.0d0
   do mm = 1,nConts
      scl = max(scl, maxval(dabs(ev1(mm,:))))
      do dd2 = 1,3
         dev = max(dev, dabs(ev1(mm,dd2)-gv1(mm,dd2,1)))
      enddo
   enddo
   print '(A,2ES12.4)', " COSX CONV PROBE: eval_bf_derivs v1 vs GTOeval val1 dev/scale = ", dev, scl
   block
      real(8) :: d0, s0, rat
      integer :: kk2
      d0 = 0.0d0 ; s0 = 0.0d0 ; rat = 0.0d0
      do mm = 1,nConts
         s0 = max(s0, dabs(ev0(mm)))
         d0 = max(d0, dabs(ev0(mm)-gv0(mm)))
         if (dabs(gv0(mm)) .gt. 1.0d-6) rat = max(rat, dabs(ev0(mm)/gv0(mm)))
      enddo
      print '(A,3ES12.4)', " COSX CONV PROBE: eval_bf v0-vs-GTOeval dev/scale, max ratio = ", d0, s0, rat
      do kk2 = 1,min(6,nConts)
         print '(A,I3,2(ES13.5,1X),2(ES13.5,1X))', " COSX CONV PROBE: AO ", kk2, &
               ev0(kk2), gv0(kk2), ev1(kk2,1), gv1(kk2,1,1)
      enddo
   end block
   deallocate(ev0,ev1,ev2,ev3,gv0,gv1)
end block

block
   real(8),allocatable :: w0a(:,:), w1a(:,:,:), w2a(:,:,:,:)
   real(8),allocatable :: wp1(:,:,:), wm1(:,:,:), w0d(:,:,:)
   real(8) :: qp(3,1), hp2, dv, sc
   integer :: mm2, mu2, nu2, kk2
   allocate(w0a(nConts,1), w1a(nConts,3,1), w2a(nConts,3,3,1))
   allocate(wp1(nConts,3,1), wm1(nConts,3,1), w0d(nConts,3,1))
   hp2 = 1.0d-5
   w0a = 0.0d0 ; w1a = 0.0d0 ; w2a = 0.0d0
   qp(:,1) = cosx_shell_cen(:,si) + 0.3d0
   call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), qp, 1, &
                            cosx_ao_offset(si), nConts, w0a, w1a, w2a)
   dv = 0.0d0 ; sc = 0.0d0
   do mu2 = 1,3
      qp(:,1) = cosx_shell_cen(:,si) + 0.3d0 ; qp(mu2,1) = qp(mu2,1) + hp2
      wp1 = 0.0d0 ; w0d = 0.0d0
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), qp, 1, &
                               cosx_ao_offset(si), nConts, w0d, wp1)
      qp(mu2,1) = qp(mu2,1) - 2.0d0*hp2
      wm1 = 0.0d0 ; w0d = 0.0d0
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), qp, 1, &
                               cosx_ao_offset(si), nConts, w0d, wm1)
      do kk2 = 1,cosx_shell_dim(si)
         mm2 = cosx_ao_offset(si) + kk2
         do nu2 = 1,3
            sc = max(sc, dabs(w2a(mm2,mu2,nu2,1)))
            dv = max(dv, dabs(w2a(mm2,mu2,nu2,1) &
                 - (wp1(mm2,nu2,1)-wm1(mm2,nu2,1))/(2.0d0*hp2)))
         enddo
      enddo
   enddo
   print '(A,2ES12.4)', " COSX CONV PROBE: val2 vs FD(val1) dev/scale = ", dv, sc
   deallocate(w0a,w1a,w2a,wp1,wm1,w0d)
end block

block
   real(8),allocatable :: ipp(:), ipv(:), pp(:), pm(:)
   real(8) :: res, sc, tt
   integer :: a, b, da, db, mu, nu, ii, jj, cmp, shpr(4), role
   allocate(ipp(di*dj*9), ipv(di*dj*9), pp(di*dj*3), pm(di*dj*3))
   do role = 1,2
      if (role == 1) then
         a = si ; b = sj
      else
         a = sj ; b = si
      endif
      da = cosx_shell_dim(a) ; db = cosx_shell_dim(b)
      shpr(1) = a ; shpr(2) = b ; shpr(3) = 0 ; shpr(4) = 1
      envp(offg+1:offg+3) = rg
      call grids1e_ipip_engine_cached(ipp, shpr, atm, size(atm,2), bas, nBases, envp, 0_8, icache)
      call grids1e_ipvip_engine_cached(ipv, shpr, atm, size(atm,2), bas, nBases, envp, 0_8, icache)
      res = 0.0d0 ; sc = 0.0d0
      do nu = 1,3
         envp(offg+1:offg+3) = rg ; envp(offg+nu) = rg(nu) + HPROBE
         call grids1e_ip_engine_cached(pp, shpr, atm, size(atm,2), bas, nBases, envp, 0_8, icache)
         envp(offg+1:offg+3) = rg ; envp(offg+nu) = rg(nu) - HPROBE
         call grids1e_ip_engine_cached(pm, shpr, atm, size(atm,2), bas, nBases, envp, 0_8, icache)
         do mu = 1,3
            cmp = (mu-1)*3 + (nu-1) + 1
            do jj = 1,db
               do ii = 1,da
                  tt = ( pp((ii-1)+da*((jj-1)+db*(mu-1))+1) &
                       - pm((ii-1)+da*((jj-1)+db*(mu-1))+1) ) / (2.0d0*HPROBE)
                  sc = sc + tt*tt
                  res = res + ( ipp((ii-1)+da*((jj-1)+db*(cmp-1))+1) &
                              + ipv((ii-1)+da*((jj-1)+db*(cmp-1))+1) - tt )**2
               enddo
            enddo
         enddo
      enddo
      print '(A,I1,A,ES12.4)', " COSX CONV PROBE: 2nd-deriv identity (role ", role, &
            ") rel.residual =", sqrt(res/max(sc,1.0d-300))
   enddo
   deallocate(ipp,ipv,pp,pm)
end block
call flush(6)
deallocate(v0,v0p,v0m,v1,envp,abuf,abufp,abufm,dbra,dket,icache)
end subroutine cosx_grad_convention_probe

module subroutine integrals_force_exchange_cosx(nConts, natoms, Da, Db, RSomega, forceKa, forceKb)
use GRID_info, only: cosx_becke_treutler
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
real(8),intent(in) :: RSomega
real(8),intent(out) :: forceKa(natoms,3), forceKb(natoms,3)
character(len=8) :: probe_env

if (.not. cosx_ready .or. cosx_nBases_cached .ne. nBases) then
   print *, "COSX force: called before the COSX grid/screening cache ", &
            "was built - SCF must converge with cosx_enabled first"
   call flush(6)
   stop 1
endif

call get_environment_variable("ENGINE_COSX_CONV_PROBE", probe_env)
if (trim(probe_env) .eq. "1") call cosx_grad_convention_probe(nConts, RSomega)

cosx_becke_treutler = cosx_hi_treutler
call cosx_force_exchange_v2(nConts, natoms, Da, Db, RSomega, forceKa, forceKb)
cosx_becke_treutler = .false.

end subroutine integrals_force_exchange_cosx
end submodule exchange_cosx_force_impl
