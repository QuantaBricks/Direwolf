! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! read basis set and integral

subroutine calc_force(info)
use MOL_info
use GRID_info
use mod_integrals, only: integrals_compute_force, integrals_compute_force_1e, &
                          integrals_force_df, engine_use_df_j, engine_use_df_k, &
                          integrals_ecp_force, integrals_force_df_lr, integrals_compute_force_exchange_lr, &
                          integrals_nuc_deriv_allatoms, &
                          engine_puream
use mod_exchange, only: HF_exchange_frac, RS_omega, RS_beta, cosx_enabled, exchange_force_cosx
use mod_xc, only: xc_uses_tau, xc_uses_vv10
use mod_vv10, only: vv10_nlc_grid_build, vv10_evaluate, vv10_grid_force, &
                          vv10_beta_val, vv10_nlc_radial, vv10_nlc_angular, vv10_rcut_val
use mod_profile, only: prof_start, prof_stop
use mod_cosmo, only: cosmo_enabled, cosmo_force_step
    implicit none
INCLUDE 'parameter.h'
    integer    :: i,j,k,l,info,i2,j2
    integer    :: igrd
    integer    :: ii,n_sig,n_sig_wide,sig_idx(nconts)
    integer,allocatable :: narrow_idx(:)
    real(8) :: SIG_CUTOFF

    real(8),allocatable :: dS(:,:,:),dHcore(:,:,:),dJi(:,:,:),dKa(:,:,:),dKb(:,:,:)
    real(8),allocatable :: dKa_lr(:,:,:),dKb_lr(:,:,:)
    real(8),allocatable :: force_J_df(:,:),force_Ka_df(:,:),force_Kb_df(:,:)
    real(8),allocatable :: force_Ka_lr_df(:,:),force_Kb_lr_df(:,:)
    real(8),allocatable :: force_Ka_cosx(:,:),force_Kb_cosx(:,:)
    real(8),allocatable :: force_Ka_cosx_lr(:,:),force_Kb_cosx_lr(:,:)
    logical :: need_exact_family, need_df_family
    real(8),allocatable :: fxc_a(:,:),fxc_b(:,:)
    real(8),allocatable :: nucderiv(:,:)
    real(8),allocatable :: ecpForce(:,:)
    real(8),allocatable :: coor_ang_cosmo(:,:), force_cosmo(:,:)
    integer,allocatable :: ao_atom_cosmo(:)
    real(8)   :: Drforce(3)
    real(8),allocatable :: W(:,:)
    real(8)             :: NucVec(3),NucForce(3)
    integer,parameter :: BATCH = DFT_BATCH
    integer    :: batch_start,nb,ip,igrid
    integer    :: idim0,a,b,c,kk,ll,dd
    integer    :: nbatches,ib,ib_cache,local_start
    integer,parameter :: MAX_HESS_TERM = 6
    integer    :: bf_atom(nconts),bf_shell(nconts),bf_nterm(nconts)
    integer    :: bf_term_a(MAX_HESS_TERM,nconts),bf_term_b(MAX_HESS_TERM,nconts)
    integer    :: bf_term_c(MAX_HESS_TERM,nconts),bf_term_idim(MAX_HESS_TERM,nconts)
    real(8)    :: bf_term_coef(MAX_HESS_TERM,nconts)
    integer    :: kk2
    real(8)    :: coor_batch(3,BATCH),COORR_batch(3,BATCH),eval_batch(6,BATCH)
    integer,parameter :: MAX_LPOW = 6
    real(8),allocatable :: hpx(:,:),hpy(:,:),hpz(:,:),hex(:,:)
    integer :: cur_atom,cur_shell,maxGauss
    real(8),allocatable :: HessPacked(:,:,:)
    real(8),allocatable :: Xop(:,:),Xop_b(:,:),ABC(:,:),ABC_b(:,:)
    real(8),allocatable :: val1_sig(:,:,:)
    integer :: nrhs
    real(8),allocatable :: val0_wide(:,:),val1_wide(:,:,:)
    real(8),allocatable :: maxval0_bf(:),maxval1_bf(:)
    integer,allocatable :: wide_idx(:)
    integer :: wide_ld,saved_geo_level
    real(8),allocatable :: Pbatch(:,:,:)
    real(8),allocatable :: Pbatch_b(:,:,:)
    logical :: mgga_active
    logical :: vv10_active
    integer,parameter :: HK(3,3) = reshape([1,2,4, 2,3,5, 4,5,6],[3,3])
    real(8),allocatable :: fxc_a_local(:,:),fxc_b_local(:,:)
    real(8),allocatable :: Pgat(:,:),Pgat_b(:,:),accf(:,:),accf_b(:,:)
    real(8),allocatable :: tauw(:)
    integer :: jj

allocate( dS (nConts,nConts,3) )
allocate( dHcore (nConts,nConts,3))
need_exact_family = (.not. engine_use_df_j) .or. &
                     ((.not. engine_use_df_k) .and. (.not. cosx_enabled))
need_df_family = engine_use_df_j .or. (engine_use_df_k .and. .not. cosx_enabled)

if (need_exact_family) then
   allocate( dJi (nConts,nConts,3) )
   allocate( dKa (nConts,nConts,3))
   allocate( dKb (nConts,nConts,3))
   call prof_start("force_2e_deriv")
   call integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
   call prof_stop("force_2e_deriv")
   if (RS_omega .gt. 0.0d0 .and. .not. cosx_enabled .and. .not. engine_use_df_k) then
      allocate(dKa_lr(nConts,nConts,3), dKb_lr(nConts,nConts,3))
      call prof_start("force_2e_deriv_exact_lr")
      call integrals_compute_force_exchange_lr(nConts, Pa, Pb, dKa_lr, dKb_lr)
      call prof_stop("force_2e_deriv_exact_lr")
   endif
endif

if (need_df_family) then
   if (.not. need_exact_family) then
      call prof_start("force_2e_deriv")
      call integrals_compute_force_1e(nConts, dS, dHcore)
      call prof_stop("force_2e_deriv")
   endif
   allocate(force_J_df(natoms,3), force_Ka_df(natoms,3), force_Kb_df(natoms,3))
   call prof_start("force_2e_deriv_df")
   call integrals_force_df(nConts, natoms, Pa, Pb, C_a, n_alpha, C_b, n_beta, &
                            force_J_df, force_Ka_df, force_Kb_df, &
                            engine_use_df_k .and. HF_exchange_frac .gt. 0.0d0)
   call prof_stop("force_2e_deriv_df")
   if (engine_use_df_k) then
      allocate(force_Ka_lr_df(natoms,3), force_Kb_lr_df(natoms,3))
      call prof_start("force_2e_deriv_df_lr")
      call integrals_force_df_lr(nConts, natoms, C_a, n_alpha, C_b, n_beta, &
                                  force_Ka_lr_df, force_Kb_lr_df)
      call prof_stop("force_2e_deriv_df_lr")
   endif
endif

if (cosx_enabled .and. HF_exchange_frac .gt. 0.0d0) then
   allocate(force_Ka_cosx(natoms,3), force_Kb_cosx(natoms,3))
   call prof_start("force_2e_deriv_cosx")
   call exchange_force_cosx(nConts, natoms, Pa, Pb, 0.0d0, force_Ka_cosx, force_Kb_cosx)
   call prof_stop("force_2e_deriv_cosx")
   if (RS_omega .gt. 0.0d0) then
      allocate(force_Ka_cosx_lr(natoms,3), force_Kb_cosx_lr(natoms,3))
      call prof_start("force_2e_deriv_cosx_lr")
      call exchange_force_cosx(nConts, natoms, Pa, Pb, RS_omega, force_Ka_cosx_lr, force_Kb_cosx_lr)
      call prof_stop("force_2e_deriv_cosx_lr")
   endif
endif

allocate(ecpForce(natoms,3))
call integrals_ecp_force(nConts, natoms, Pa+Pb, ecpForce)

allocate(W(nconts,nconts))
W = 0
do i = 1,nconts
do j = 1,nconts

   do k = 1,n_alpha
      W(i,j) =eLev_a(k)*C_a(i,k)*C_a(j,k) + W(i,j)
   enddo

   if  ( multi .eq. 1) then
      W(i,j) = W(i,j)*2
   else
      do k = 1,n_beta
         W(i,j) =eLev_b(k)*C_b(i,k)*C_b(j,k) + W(i,j)
      enddo
   endif

enddo
enddo

allocate(fxc_a(3,natoms))
allocate(fxc_b(3,natoms))
fxc_a = 0
fxc_b = 0

call prof_start("force_xc_grid")
mgga_active = xc_uses_tau()
vv10_active = xc_uses_vv10()
maxGauss = 1
do i = 1,natoms
   do j = 1,size(atoms(i)%shell)
      maxGauss = max(maxGauss, atoms(i)%shell(j)%nGauss)
      if (atoms(i)%shell(j)%angMoment+2 .gt. MAX_LPOW) then
         print *, 'ERROR: force grid Hessian angular momentum too high:', &
                  atoms(i)%shell(j)%angMoment
         stop
      endif
   enddo
enddo
SIG_CUTOFF = XC_SIG_TIGHT
block
  character(32) :: sigenv
  call get_environment_variable("ENGINE_FORCE_XC_SIG", sigenv)
  if (len_trim(sigenv) .gt. 0) read(sigenv,*) SIG_CUTOFF
end block
if (HF_exchange_frac .lt. 1.0d0) then

l = 0
do i = 1,natoms
   do j = 1,size(atoms(i)%shell)
      if (engine_puream .and. atoms(i)%shell(j)%angMoment .ge. 5) then
         print *, 'ERROR: puream force grid Hessian for Lt=', &
                  atoms(i)%shell(j)%angMoment, ' not implemented'
         stop 1
      else if (engine_puream .and. atoms(i)%shell(j)%angMoment .eq. 4) then
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_term_a(1,l)=3;bf_term_b(1,l)=1;bf_term_c(1,l)=0; bf_term_idim(1,l)=2; bf_term_coef(1,l)=sqrt(5.0d0)/2.0d0
         bf_term_a(2,l)=1;bf_term_b(2,l)=3;bf_term_c(2,l)=0; bf_term_idim(2,l)=7; bf_term_coef(2,l)=-sqrt(5.0d0)/2.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_term_a(1,l)=2;bf_term_b(1,l)=1;bf_term_c(1,l)=1; bf_term_idim(1,l)=5; bf_term_coef(1,l)=3.0d0*sqrt(2.0d0)/4.0d0
         bf_term_a(2,l)=0;bf_term_b(2,l)=3;bf_term_c(2,l)=1; bf_term_idim(2,l)=12; bf_term_coef(2,l)=-sqrt(10.0d0)/4.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_term_a(1,l)=1;bf_term_b(1,l)=1;bf_term_c(1,l)=2; bf_term_idim(1,l)=9; bf_term_coef(1,l)=3.0d0*sqrt(7.0d0)/7.0d0
         bf_term_a(2,l)=3;bf_term_b(2,l)=1;bf_term_c(2,l)=0; bf_term_idim(2,l)=2; bf_term_coef(2,l)=-sqrt(35.0d0)/14.0d0
         bf_term_a(3,l)=1;bf_term_b(3,l)=3;bf_term_c(3,l)=0; bf_term_idim(3,l)=7; bf_term_coef(3,l)=-sqrt(35.0d0)/14.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_term_a(1,l)=0;bf_term_b(1,l)=1;bf_term_c(1,l)=3; bf_term_idim(1,l)=14; bf_term_coef(1,l)=sqrt(70.0d0)/7.0d0
         bf_term_a(2,l)=2;bf_term_b(2,l)=1;bf_term_c(2,l)=1; bf_term_idim(2,l)=5; bf_term_coef(2,l)=-3.0d0*sqrt(14.0d0)/28.0d0
         bf_term_a(3,l)=0;bf_term_b(3,l)=3;bf_term_c(3,l)=1; bf_term_idim(3,l)=12; bf_term_coef(3,l)=-3.0d0*sqrt(70.0d0)/28.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=6
         bf_term_a(1,l)=0;bf_term_b(1,l)=0;bf_term_c(1,l)=4; bf_term_idim(1,l)=15; bf_term_coef(1,l)=1.0d0
         bf_term_a(2,l)=2;bf_term_b(2,l)=0;bf_term_c(2,l)=2; bf_term_idim(2,l)=6; bf_term_coef(2,l)=-3.0d0*sqrt(105.0d0)/35.0d0
         bf_term_a(3,l)=0;bf_term_b(3,l)=2;bf_term_c(3,l)=2; bf_term_idim(3,l)=13; bf_term_coef(3,l)=-3.0d0*sqrt(105.0d0)/35.0d0
         bf_term_a(4,l)=4;bf_term_b(4,l)=0;bf_term_c(4,l)=0; bf_term_idim(4,l)=1; bf_term_coef(4,l)=3.0d0/8.0d0
         bf_term_a(5,l)=2;bf_term_b(5,l)=2;bf_term_c(5,l)=0; bf_term_idim(5,l)=4; bf_term_coef(5,l)=3.0d0*sqrt(105.0d0)/140.0d0
         bf_term_a(6,l)=0;bf_term_b(6,l)=4;bf_term_c(6,l)=0; bf_term_idim(6,l)=11; bf_term_coef(6,l)=3.0d0/8.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_term_a(1,l)=1;bf_term_b(1,l)=0;bf_term_c(1,l)=3; bf_term_idim(1,l)=10; bf_term_coef(1,l)=sqrt(70.0d0)/7.0d0
         bf_term_a(2,l)=3;bf_term_b(2,l)=0;bf_term_c(2,l)=1; bf_term_idim(2,l)=3; bf_term_coef(2,l)=-3.0d0*sqrt(70.0d0)/28.0d0
         bf_term_a(3,l)=1;bf_term_b(3,l)=2;bf_term_c(3,l)=1; bf_term_idim(3,l)=8; bf_term_coef(3,l)=-3.0d0*sqrt(14.0d0)/28.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=4
         bf_term_a(1,l)=2;bf_term_b(1,l)=0;bf_term_c(1,l)=2; bf_term_idim(1,l)=6; bf_term_coef(1,l)=3.0d0*sqrt(21.0d0)/14.0d0
         bf_term_a(2,l)=0;bf_term_b(2,l)=2;bf_term_c(2,l)=2; bf_term_idim(2,l)=13; bf_term_coef(2,l)=-3.0d0*sqrt(21.0d0)/14.0d0
         bf_term_a(3,l)=4;bf_term_b(3,l)=0;bf_term_c(3,l)=0; bf_term_idim(3,l)=1; bf_term_coef(3,l)=-sqrt(5.0d0)/4.0d0
         bf_term_a(4,l)=0;bf_term_b(4,l)=4;bf_term_c(4,l)=0; bf_term_idim(4,l)=11; bf_term_coef(4,l)=sqrt(5.0d0)/4.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_term_a(1,l)=3;bf_term_b(1,l)=0;bf_term_c(1,l)=1; bf_term_idim(1,l)=3; bf_term_coef(1,l)=sqrt(10.0d0)/4.0d0
         bf_term_a(2,l)=1;bf_term_b(2,l)=2;bf_term_c(2,l)=1; bf_term_idim(2,l)=8; bf_term_coef(2,l)=-3.0d0*sqrt(2.0d0)/4.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_term_a(1,l)=4;bf_term_b(1,l)=0;bf_term_c(1,l)=0; bf_term_idim(1,l)=1; bf_term_coef(1,l)=sqrt(35.0d0)/8.0d0
         bf_term_a(2,l)=2;bf_term_b(2,l)=2;bf_term_c(2,l)=0; bf_term_idim(2,l)=4; bf_term_coef(2,l)=-3.0d0*sqrt(3.0d0)/4.0d0
         bf_term_a(3,l)=0;bf_term_b(3,l)=4;bf_term_c(3,l)=0; bf_term_idim(3,l)=11; bf_term_coef(3,l)=sqrt(35.0d0)/8.0d0
      else if (engine_puream .and. atoms(i)%shell(j)%angMoment .eq. 3) then
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_term_a(1,l)=2;bf_term_b(1,l)=1;bf_term_c(1,l)=0; bf_term_idim(1,l)=2; bf_term_coef(1,l)=3.0d0*sqrt(2.0d0)/4.0d0
         bf_term_a(2,l)=0;bf_term_b(2,l)=3;bf_term_c(2,l)=0; bf_term_idim(2,l)=7; bf_term_coef(2,l)=-sqrt(10.0d0)/4.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
         bf_term_a(1,l)=1;bf_term_b(1,l)=1;bf_term_c(1,l)=1; bf_term_idim(1,l)=5; bf_term_coef(1,l)=1.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_term_a(1,l)=0;bf_term_b(1,l)=1;bf_term_c(1,l)=2; bf_term_idim(1,l)=9; bf_term_coef(1,l)=sqrt(30.0d0)/5.0d0
         bf_term_a(2,l)=2;bf_term_b(2,l)=1;bf_term_c(2,l)=0; bf_term_idim(2,l)=2; bf_term_coef(2,l)=-sqrt(30.0d0)/20.0d0
         bf_term_a(3,l)=0;bf_term_b(3,l)=3;bf_term_c(3,l)=0; bf_term_idim(3,l)=7; bf_term_coef(3,l)=-sqrt(6.0d0)/4.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_term_a(1,l)=0;bf_term_b(1,l)=0;bf_term_c(1,l)=3; bf_term_idim(1,l)=10; bf_term_coef(1,l)=1.0d0
         bf_term_a(2,l)=2;bf_term_b(2,l)=0;bf_term_c(2,l)=1; bf_term_idim(2,l)=3; bf_term_coef(2,l)=-3.0d0*sqrt(5.0d0)/10.0d0
         bf_term_a(3,l)=0;bf_term_b(3,l)=2;bf_term_c(3,l)=1; bf_term_idim(3,l)=8; bf_term_coef(3,l)=-3.0d0*sqrt(5.0d0)/10.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_term_a(1,l)=1;bf_term_b(1,l)=0;bf_term_c(1,l)=2; bf_term_idim(1,l)=6; bf_term_coef(1,l)=sqrt(30.0d0)/5.0d0
         bf_term_a(2,l)=3;bf_term_b(2,l)=0;bf_term_c(2,l)=0; bf_term_idim(2,l)=1; bf_term_coef(2,l)=-sqrt(6.0d0)/4.0d0
         bf_term_a(3,l)=1;bf_term_b(3,l)=2;bf_term_c(3,l)=0; bf_term_idim(3,l)=4; bf_term_coef(3,l)=-sqrt(30.0d0)/20.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_term_a(1,l)=2;bf_term_b(1,l)=0;bf_term_c(1,l)=1; bf_term_idim(1,l)=3; bf_term_coef(1,l)=sqrt(3.0d0)/2.0d0
         bf_term_a(2,l)=0;bf_term_b(2,l)=2;bf_term_c(2,l)=1; bf_term_idim(2,l)=8; bf_term_coef(2,l)=-sqrt(3.0d0)/2.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_term_a(1,l)=3;bf_term_b(1,l)=0;bf_term_c(1,l)=0; bf_term_idim(1,l)=1; bf_term_coef(1,l)=sqrt(10.0d0)/4.0d0
         bf_term_a(2,l)=1;bf_term_b(2,l)=2;bf_term_c(2,l)=0; bf_term_idim(2,l)=4; bf_term_coef(2,l)=-3.0d0*sqrt(2.0d0)/4.0d0
      else if (engine_puream .and. atoms(i)%shell(j)%angMoment .eq. 2) then
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
         bf_term_a(1,l)=1; bf_term_b(1,l)=1; bf_term_c(1,l)=0
         bf_term_idim(1,l)=2; bf_term_coef(1,l)=1.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
         bf_term_a(1,l)=0; bf_term_b(1,l)=1; bf_term_c(1,l)=1
         bf_term_idim(1,l)=5; bf_term_coef(1,l)=1.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_term_a(1,l)=0; bf_term_b(1,l)=0; bf_term_c(1,l)=2
         bf_term_idim(1,l)=6; bf_term_coef(1,l)=1.0d0
         bf_term_a(2,l)=2; bf_term_b(2,l)=0; bf_term_c(2,l)=0
         bf_term_idim(2,l)=1; bf_term_coef(2,l)=-0.5d0
         bf_term_a(3,l)=0; bf_term_b(3,l)=2; bf_term_c(3,l)=0
         bf_term_idim(3,l)=4; bf_term_coef(3,l)=-0.5d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
         bf_term_a(1,l)=1; bf_term_b(1,l)=0; bf_term_c(1,l)=1
         bf_term_idim(1,l)=3; bf_term_coef(1,l)=1.0d0
         l = l+1
         bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_term_a(1,l)=2; bf_term_b(1,l)=0; bf_term_c(1,l)=0
         bf_term_idim(1,l)=1; bf_term_coef(1,l)=sqrt(3.0d0)/2.0d0
         bf_term_a(2,l)=0; bf_term_b(2,l)=2; bf_term_c(2,l)=0
         bf_term_idim(2,l)=4; bf_term_coef(2,l)=-sqrt(3.0d0)/2.0d0
      else
         idim0 = 0
         do a = atoms(i)%shell(j)%angMoment,0,-1
            do b = atoms(i)%shell(j)%angMoment-a,0,-1
               c = atoms(i)%shell(j)%angMoment-a-b
               l = l+1
               idim0 = idim0+1
               bf_atom(l) = i
               bf_shell(l) = j
               bf_nterm(l) = 1
               bf_term_a(1,l) = a
               bf_term_b(1,l) = b
               bf_term_c(1,l) = c
               bf_term_idim(1,l) = idim0
               bf_term_coef(1,l) = 1.0d0
            enddo
         enddo
      endif
   enddo
enddo

nbatches = (ngrids + BATCH - 1) / BATCH

wide_ld = max(max_batch_nsig, max_dft_batch_nsig)

saved_geo_level = xc_geo_level
xc_geo_level = XC_NLEVEL

!$omp parallel &
!$omp&   private(ib,batch_start,nb,ib_cache,local_start,igrid,ip,ii,kk,ll,n_sig,n_sig_wide,sig_idx,narrow_idx) &
!$omp&   private(l,i,j,a,b,c,k,dd,kk2) &
!$omp&   private(coor_batch,COORR_batch,eval_batch) &
!$omp&   private(HessPacked,val1_sig,Pbatch,Pbatch_b,Xop,Xop_b,ABC,ABC_b,nrhs) &
!$omp&   private(hpx,hpy,hpz,hex,cur_atom,cur_shell) &
!$omp&   private(val0_wide,val1_wide,maxval0_bf,maxval1_bf,wide_idx) &
!$omp&   private(i2) &
!$omp&   private(fxc_a_local,fxc_b_local,jj) &
!$omp&   private(Pgat,Pgat_b,accf,accf_b,tauw)
allocate(HessPacked(nconts,6,BATCH))
nrhs = merge(5,2,mgga_active)
allocate(Xop(nconts,nrhs*BATCH),ABC(nconts,nrhs*BATCH))
allocate(val1_sig(nconts,3,BATCH))
allocate(val0_wide(wide_ld,BATCH),val1_wide(wide_ld,3,BATCH))
allocate(maxval0_bf(wide_ld),maxval1_bf(wide_ld))
allocate(narrow_idx(wide_ld),wide_idx(wide_ld))
allocate(hpx(0:MAX_LPOW,BATCH),hpy(0:MAX_LPOW,BATCH),hpz(0:MAX_LPOW,BATCH))
allocate(hex(BATCH,maxGauss))
allocate(Pbatch(nconts,3,BATCH))
if (multi .ne. 1) then
   allocate(Pbatch_b(nconts,3,BATCH))
   allocate(Xop_b(nconts,nrhs*BATCH),ABC_b(nconts,nrhs*BATCH))
endif
allocate(Pgat(nconts,nconts))
allocate(accf(nconts,3),tauw(BATCH))
allocate(fxc_a_local(3,natoms))
fxc_a_local = 0
allocate(fxc_b_local(3,natoms))
fxc_b_local = 0
if (multi .ne. 1) then
   allocate(Pgat_b(nconts,nconts))
   allocate(accf_b(nconts,3))
endif

!$omp do schedule(dynamic)
do ib = 0,nbatches-1
   batch_start = ib*BATCH + 1
   nb = min(BATCH, ngrids-batch_start+1)

   do ip = 1,nb
      coor_batch(:,ip) = grids(batch_start+ip-1)%coor
   enddo

   ib_cache = (batch_start-1)/CACHE_BATCH + 1

   if (grid_cache_mode .and. batch_start+nb-1 .le. grid_cache_capacity) then
      n_sig_wide = batch_nsig(ib_cache)
      wide_idx(1:n_sig_wide) = batch_sig_idx(1:n_sig_wide,ib_cache)
      local_start = batch_start - (ib_cache-1)*CACHE_BATCH
      val0_wide(1:n_sig_wide,1:nb) = val_blocks(ib_cache)%val0(1:n_sig_wide,local_start:local_start+nb-1)
      val1_wide(1:n_sig_wide,1:3,1:nb) = val_blocks(ib_cache)%val1(1:n_sig_wide,1:3,local_start:local_start+nb-1)
   else
      n_sig_wide = dft_batch_shells(ib+1)%n_sig
      wide_idx(1:n_sig_wide) = dft_batch_shells(ib+1)%sig_idx(1:n_sig_wide)
      call GTOeval_batch_compact(ib+1,batch_start,nb,wide_ld,val0_wide,val1_wide)
   endif

   maxval0_bf(1:n_sig_wide) = 0.0d0
   maxval1_bf(1:n_sig_wide) = 0.0d0
   do ip = 1,nb
      maxval0_bf(1:n_sig_wide) = max(maxval0_bf(1:n_sig_wide), abs(val0_wide(1:n_sig_wide,ip)))
      maxval1_bf(1:n_sig_wide) = max(maxval1_bf(1:n_sig_wide), &
                                      abs(val1_wide(1:n_sig_wide,1,ip)), &
                                      abs(val1_wide(1:n_sig_wide,2,ip)), &
                                      abs(val1_wide(1:n_sig_wide,3,ip)))
   enddo
   n_sig = 0
   do ii = 1,n_sig_wide
      if (maxval0_bf(ii) .gt. SIG_CUTOFF .or. maxval1_bf(ii) .gt. SIG_CUTOFF) then
         n_sig = n_sig + 1
         sig_idx(n_sig) = wide_idx(ii)
         narrow_idx(n_sig) = ii
      endif
   enddo
   Xop(1:n_sig,1:nb) = val0_wide(narrow_idx(1:n_sig),1:nb)
   val1_sig(1:n_sig,1:3,1:nb) = val1_wide(narrow_idx(1:n_sig),1:3,1:nb)
   if (mgga_active) then
      do dd = 1,3
         Xop(1:n_sig,(1+dd)*nb+1:(2+dd)*nb) = val1_sig(1:n_sig,dd,1:nb)
      enddo
   endif

   cur_atom  = 0
   cur_shell = 0
   do ii = 1,n_sig
      l = sig_idx(ii)
      i = bf_atom(l)
      j = bf_shell(l)
      if (i .ne. cur_atom .or. j .ne. cur_shell) then
         do ip = 1,nb
            COORR_batch(:,ip) = coor_batch(:,ip)-atoms(i)%coor*ans2bohr
         enddo
         call GTOHess_build_tables(COORR_batch(:,1:nb),nb,atoms(i)%shell(j)%angMoment, &
                                   atoms(i)%shell(j)%nGauss,atoms(i)%shell(j)%exponents, &
                                   MAX_LPOW,BATCH,hpx,hpy,hpz,hex)
         cur_atom  = i
         cur_shell = j
      endif
      HessPacked(ii,:,1:nb) = 0.0d0
      do kk2 = 1,bf_nterm(l)
         call GTOHess_batch_tab(nb,bf_term_a(kk2,l),bf_term_b(kk2,l),bf_term_c(kk2,l), &
                             atoms(i)%shell(j)%nGauss, &
                             atoms(i)%shell(j)%cnVal(:,bf_term_idim(kk2,l)), &
                             atoms(i)%shell(j)%exponents, &
                             MAX_LPOW,BATCH,hpx,hpy,hpz,hex,eval_batch)
         do ip = 1,nb
            do kk = 1,6
               HessPacked(ii,kk,ip) = HessPacked(ii,kk,ip) + bf_term_coef(kk2,l)*eval_batch(kk,ip)
            enddo
         enddo
      enddo
   enddo

   do ip = 1,nb
      igrid = batch_start+ip-1
      Pbatch(1:n_sig,1,ip) = ( val1_sig(1:n_sig,1,ip)*TempD_all(1,igrid) &
                             + TempD_all(2,igrid)*HessPacked(1:n_sig,1,ip) &
                             + TempD_all(3,igrid)*HessPacked(1:n_sig,2,ip) &
                             + TempD_all(4,igrid)*HessPacked(1:n_sig,4,ip) ) * grids(igrid)%weight
      Pbatch(1:n_sig,2,ip) = ( val1_sig(1:n_sig,2,ip)*TempD_all(1,igrid) &
                             + TempD_all(2,igrid)*HessPacked(1:n_sig,2,ip) &
                             + TempD_all(3,igrid)*HessPacked(1:n_sig,3,ip) &
                             + TempD_all(4,igrid)*HessPacked(1:n_sig,5,ip) ) * grids(igrid)%weight
      Pbatch(1:n_sig,3,ip) = ( val1_sig(1:n_sig,3,ip)*TempD_all(1,igrid) &
                             + TempD_all(2,igrid)*HessPacked(1:n_sig,4,ip) &
                             + TempD_all(3,igrid)*HessPacked(1:n_sig,5,ip) &
                             + TempD_all(4,igrid)*HessPacked(1:n_sig,6,ip) ) * grids(igrid)%weight
      Xop(1:n_sig,nb+ip) = ( TempD_all(2,igrid)*val1_sig(1:n_sig,1,ip) &
                           + TempD_all(3,igrid)*val1_sig(1:n_sig,2,ip) &
                           + TempD_all(4,igrid)*val1_sig(1:n_sig,3,ip) ) * grids(igrid)%weight
   enddo
   if (multi .ne. 1) then
      do ip = 1,nb
         igrid = batch_start+ip-1
         Pbatch_b(1:n_sig,1,ip) = ( val1_sig(1:n_sig,1,ip)*TempD_all(5,igrid) &
                                  + TempD_all(6,igrid)*HessPacked(1:n_sig,1,ip) &
                                  + TempD_all(7,igrid)*HessPacked(1:n_sig,2,ip) &
                                  + TempD_all(8,igrid)*HessPacked(1:n_sig,4,ip) ) * grids(igrid)%weight
         Pbatch_b(1:n_sig,2,ip) = ( val1_sig(1:n_sig,2,ip)*TempD_all(5,igrid) &
                                  + TempD_all(6,igrid)*HessPacked(1:n_sig,2,ip) &
                                  + TempD_all(7,igrid)*HessPacked(1:n_sig,3,ip) &
                                  + TempD_all(8,igrid)*HessPacked(1:n_sig,5,ip) ) * grids(igrid)%weight
         Pbatch_b(1:n_sig,3,ip) = ( val1_sig(1:n_sig,3,ip)*TempD_all(5,igrid) &
                                  + TempD_all(6,igrid)*HessPacked(1:n_sig,4,ip) &
                                  + TempD_all(7,igrid)*HessPacked(1:n_sig,5,ip) &
                                  + TempD_all(8,igrid)*HessPacked(1:n_sig,6,ip) ) * grids(igrid)%weight
         Xop_b(1:n_sig,nb+ip) = ( TempD_all(6,igrid)*val1_sig(1:n_sig,1,ip) &
                                + TempD_all(7,igrid)*val1_sig(1:n_sig,2,ip) &
                                + TempD_all(8,igrid)*val1_sig(1:n_sig,3,ip) ) * grids(igrid)%weight
      enddo
   endif

   if (multi .eq. 1) then
      do jj = 1,n_sig
         Pgat(1:n_sig,jj) = Pa(sig_idx(1:n_sig),sig_idx(jj)) + Pb(sig_idx(1:n_sig),sig_idx(jj))
      enddo
   else
      do jj = 1,n_sig
         Pgat(1:n_sig,jj)   = Pa(sig_idx(1:n_sig),sig_idx(jj))
         Pgat_b(1:n_sig,jj) = Pb(sig_idx(1:n_sig),sig_idx(jj))
      enddo
   endif

   call dgemm('N','N',n_sig,nrhs*nb,n_sig,1.0d0,Pgat,nconts,Xop,nconts,0.0d0,ABC,nconts)
   if (mgga_active) then
      do ip = 1,nb
         igrid = batch_start+ip-1
         tauw(ip) = 0.5d0*TempD_all(9,igrid)*grids(igrid)%weight
      enddo
   endif

   accf(1:n_sig,1:3) = 0.0d0
   do k = 1,3
      do ip = 1,nb
         accf(1:n_sig,k) = accf(1:n_sig,k) &
                         + Pbatch(1:n_sig,k,ip)*ABC(1:n_sig,ip) &
                         + val1_sig(1:n_sig,k,ip)*ABC(1:n_sig,nb+ip)
      enddo
      if (mgga_active) then
         do dd = 1,3
            do ip = 1,nb
               accf(1:n_sig,k) = accf(1:n_sig,k) &
                  + tauw(ip)*HessPacked(1:n_sig,HK(dd,k),ip)*ABC(1:n_sig,(1+dd)*nb+ip)
            enddo
         enddo
      endif
   enddo
   do ii = 1,n_sig
      i2 = bf_atom(sig_idx(ii))
      fxc_a_local(1,i2) = fxc_a_local(1,i2) + accf(ii,1)
      fxc_a_local(2,i2) = fxc_a_local(2,i2) + accf(ii,2)
      fxc_a_local(3,i2) = fxc_a_local(3,i2) + accf(ii,3)
   enddo

   if (multi .ne. 1) then
      Xop_b(1:n_sig,1:nb) = Xop(1:n_sig,1:nb)
      if (mgga_active) Xop_b(1:n_sig,2*nb+1:5*nb) = Xop(1:n_sig,2*nb+1:5*nb)
      call dgemm('N','N',n_sig,nrhs*nb,n_sig,1.0d0,Pgat_b,nconts,Xop_b,nconts,0.0d0,ABC_b,nconts)
      if (mgga_active) then
         do ip = 1,nb
            igrid = batch_start+ip-1
            tauw(ip) = 0.5d0*TempD_all(10,igrid)*grids(igrid)%weight
         enddo
      endif
      accf_b(1:n_sig,1:3) = 0.0d0
      do k = 1,3
         do ip = 1,nb
            accf_b(1:n_sig,k) = accf_b(1:n_sig,k) &
                              + Pbatch_b(1:n_sig,k,ip)*ABC_b(1:n_sig,ip) &
                              + val1_sig(1:n_sig,k,ip)*ABC_b(1:n_sig,nb+ip)
         enddo
         if (mgga_active) then
            do dd = 1,3
               do ip = 1,nb
                  accf_b(1:n_sig,k) = accf_b(1:n_sig,k) &
                     + tauw(ip)*HessPacked(1:n_sig,HK(dd,k),ip)*ABC_b(1:n_sig,(1+dd)*nb+ip)
               enddo
            enddo
         endif
      enddo
      do ii = 1,n_sig
         i2 = bf_atom(sig_idx(ii))
         fxc_b_local(1,i2) = fxc_b_local(1,i2) + accf_b(ii,1)
         fxc_b_local(2,i2) = fxc_b_local(2,i2) + accf_b(ii,2)
         fxc_b_local(3,i2) = fxc_b_local(3,i2) + accf_b(ii,3)
      enddo
   endif
enddo
!$omp end do

deallocate(HessPacked,val1_sig,Pbatch,Xop,ABC)
deallocate(hpx,hpy,hpz,hex)
deallocate(Pgat,accf,tauw)
deallocate(val0_wide,val1_wide,maxval0_bf,maxval1_bf,narrow_idx,wide_idx)
if (allocated(Pbatch_b)) deallocate(Pbatch_b,Pgat_b,accf_b,Xop_b,ABC_b)

!$omp critical
fxc_a = fxc_a + fxc_a_local
fxc_b = fxc_b + fxc_b_local
!$omp end critical
deallocate(fxc_a_local,fxc_b_local)
!$omp end parallel

xc_geo_level = saved_geo_level
endif
call prof_stop("force_xc_grid")

l = 0

allocate( nucderiv(natoms,3) )
call prof_start("force_nucderiv")
call integrals_nuc_deriv_allatoms(nConts, natoms, Pa+Pb, nucderiv)
call prof_stop("force_nucderiv")

call prof_start("force_assemble")
block
   real(8) :: dbg_hcore(3), dbg_s(3), dbg_j(3), dbg_k(3)
   character(32) :: dbgenv
   integer :: dbgstat
   dbg_hcore=0; dbg_s=0; dbg_j=0; dbg_k=0
   call get_environment_variable("ENGINE_FORCE_DEBUG_PARTS", dbgenv, status=dbgstat)
do i = 1,natoms
   atoms(i)%atmForce = 0
   atoms(i)%atmForce = atoms(i)%atmForce - 2.0D0*fxc_a(:,i) - 2.0D0*fxc_b(:,i)
   do j = 1,atoms(i)%nconts
      l = 1+l
      do k = 1,nconts
          atoms(i)%atmForce= atoms(i)%atmForce -&
           2.0D0*(Pa(l,k)+Pb(l,k))*dHcore(l,k,:)+&
           2.0D0*W(l,k)*dS(l,k,:)
          if (dbgstat==0 .and. i==1) then
             dbg_hcore = dbg_hcore - 2.0D0*(Pa(l,k)+Pb(l,k))*dHcore(l,k,:)
             dbg_s = dbg_s + 2.0D0*W(l,k)*dS(l,k,:)
          endif
          if (.not. engine_use_df_j) then
             atoms(i)%atmForce = atoms(i)%atmForce -&
              2.0D0*(Pa(l,k)+Pb(l,k))*dJi(l,k,:)
             if (dbgstat==0 .and. i==1) then
                dbg_j = dbg_j - 2.0D0*(Pa(l,k)+Pb(l,k))*dJi(l,k,:)
             endif
          endif
          if ((.not. cosx_enabled) .and. (.not. engine_use_df_k)) then
             atoms(i)%atmForce = atoms(i)%atmForce +&
              2.0D0*HF_exchange_frac*(Pa(l,k)*dKa(l,k,:)+Pb(l,k)*dKb(l,k,:))
          if (dbgstat==0 .and. i==1) then
             dbg_k = dbg_k + 2.0D0*HF_exchange_frac*(Pa(l,k)*dKa(l,k,:)+Pb(l,k)*dKb(l,k,:))
          endif
             if (RS_omega .gt. 0.0d0) then
                atoms(i)%atmForce = atoms(i)%atmForce + &
                   2.0D0*RS_beta*(Pa(l,k)*dKa_lr(l,k,:)+Pb(l,k)*dKb_lr(l,k,:))
             endif
          endif
      enddo
   enddo
   if (dbgstat==0 .and. i==1) then
      print *, 'FORCE_DEBUG_PARTS atom1: dbg_hcore=', dbg_hcore
      print *, 'FORCE_DEBUG_PARTS atom1: dbg_s     =', dbg_s
      print *, 'FORCE_DEBUG_PARTS atom1: dbg_j     =', dbg_j
      print *, 'FORCE_DEBUG_PARTS atom1: dbg_k     =', dbg_k
   endif
   if (engine_use_df_j) then
      atoms(i)%atmForce = atoms(i)%atmForce - force_J_df(i,:)
   endif
   if (engine_use_df_k) then
      atoms(i)%atmForce = atoms(i)%atmForce + &
         HF_exchange_frac*(force_Ka_df(i,:)+force_Kb_df(i,:))
      if (RS_omega .gt. 0.0d0) then
         atoms(i)%atmForce = atoms(i)%atmForce + &
            0.5d0*RS_beta*(force_Ka_lr_df(i,:)+force_Kb_lr_df(i,:))
      endif
   endif
   if (cosx_enabled .and. HF_exchange_frac .gt. 0.0d0) then
      atoms(i)%atmForce = atoms(i)%atmForce + &
         2.0D0*HF_exchange_frac*(force_Ka_cosx(i,:)+force_Kb_cosx(i,:))
      if (RS_omega .gt. 0.0d0) then
         atoms(i)%atmForce = atoms(i)%atmForce + &
            2.0D0*RS_beta*(force_Ka_cosx_lr(i,:)+force_Kb_cosx_lr(i,:))
      endif
   endif
   atoms(i)%atmForce = atoms(i)%atmForce + ecpForce(i,:)
   Nucforce = 0
   do j = 1,natoms
      if (i .ne. j) then
         NucVec = atoms(i)%coor - atoms(j)%coor
         NucVec = NucVec * ans2bohr
         Nucforce = -(atoms(i)%charge-atoms(i)%ecpCoreElec) * (atoms(j)%charge-atoms(j)%ecpCoreElec) *NucVec /&
                    (NucVec(1)**2+NucVec(2)**2+NucVec(3)**2) **(3.0D0/2) + &
                    Nucforce
      endif
   enddo
   do j = 1,nPointCharges
      NucVec = atoms(i)%coor - pointcharge_coor(:,j)
      NucVec = NucVec * ans2bohr
      Nucforce = -(atoms(i)%charge-atoms(i)%ecpCoreElec) * pointcharge_q(j) *NucVec /&
                 (NucVec(1)**2+NucVec(2)**2+NucVec(3)**2) **(3.0D0/2) + &
                 Nucforce
   enddo

   Drforce = nucderiv(i,:)

   atoms(i)%atmForce = atoms(i)%atmForce + Nucforce - 2.0D0*Drforce*(atoms(i)%charge-atoms(i)%ecpCoreElec)
enddo
end block

if (cosmo_enabled) then
   allocate(coor_ang_cosmo(3,natoms), force_cosmo(3,natoms), ao_atom_cosmo(nConts))
   l = 0
   do i = 1,natoms
      coor_ang_cosmo(:,i) = atoms(i)%coor
      do j = 1,atoms(i)%nconts
         l = l+1
         ao_atom_cosmo(l) = i
      enddo
   enddo
   call cosmo_force_step(natoms, coor_ang_cosmo, ao_atom_cosmo, nConts, Pa+Pb, force_cosmo)
   do i = 1,natoms
      atoms(i)%atmForce = atoms(i)%atmForce + force_cosmo(:,i)
   enddo
    deallocate(coor_ang_cosmo, force_cosmo, ao_atom_cosmo)
 endif

if (vv10_active) then
   call prof_start("force_vv10")
   call calc_force_vv10()
   call prof_stop("force_vv10")
endif

call prof_stop("force_assemble")

contains

subroutine calc_force_vv10()
    implicit none
    interface
       subroutine gridgen_nlc(nrad, nsph, npts, coor_out, weight_out, coor_in)
       integer,intent(in) :: nrad, nsph
       integer,intent(out) :: npts
       real(8),allocatable,intent(out) :: coor_out(:,:), weight_out(:)
       real(8),optional,intent(in) :: coor_in(3,*)
       end subroutine gridgen_nlc
       subroutine gridgen_nlc_subset(nrad, nsph, atom_list, n_list, npts, coor_out, weight_out, atom_of_out, coor_in)
       integer,intent(in) :: nrad, nsph, n_list
       integer,intent(in) :: atom_list(n_list)
       integer,intent(out) :: npts
       real(8),allocatable,intent(out) :: coor_out(:,:), weight_out(:)
       integer,allocatable,intent(out) :: atom_of_out(:)
       real(8),optional,intent(in) :: coor_in(3,*)
       end subroutine gridgen_nlc_subset
       subroutine vv10_nlc_weight_gradient(np, coor, atom_of, nper, rho, eps, wquad, beta_v, gv_wt)
       integer,intent(in) :: np, nper
       real(8),intent(in) :: coor(3,np)
       integer,intent(in) :: atom_of(np)
       real(8),intent(in) :: rho(np), eps(np), wquad(np), beta_v
       real(8),intent(out) :: gv_wt(*)
       end subroutine vv10_nlc_weight_gradient
    end interface
    integer :: np, i, mu, dd, ee, iat, A, nrad, nsph, n
    integer,parameter :: PK(3,3) = reshape([1,2,4, 2,3,5, 4,5,6],[3,3])
    real(8),allocatable :: coor(:,:), weight(:)
    real(8),allocatable :: rho(:), grad(:,:), sig(:), eps(:), Fn(:), Fg(:)
    real(8),allocatable :: gv(:,:), atom_of(:)
    real(8),allocatable :: gv_gbf(:,:), gv_grid(:,:), gv_wt(:,:)
    real(8),allocatable :: coor_p0(:,:)
    real(8) :: Enl
    real(8) :: Pdens(nconts,nconts)
    integer :: atom_bf_start(natoms+1), atom_bf_list(nconts), atom_bf_cursor(natoms)
    real(8) :: beta_v, wdel
    integer :: nper, np2
    real(8) :: atom_rcut2(natoms), dist2_A
    integer :: jsh
    real(8),parameter :: VV10_GRAD_SIG_CUTOFF = 1.0d-6
    real(8) :: bf_rcut2_grad(nconts), dist2_bf
    integer :: kk, nsig_grad, sig_idx_grad(nconts)
    real(8),parameter :: VV10_HESS_SIG_CUTOFF = 1.0d-7
    real(8) :: bf_rcut2_hess(nconts)
    integer :: nsig_hess, sig_idx_hess(nconts)

    call vv10_nlc_grid_build(np, coor, weight)
    nper = np/natoms
    allocate(rho(np), grad(3,np), sig(np), eps(np), Fn(np), Fg(np))
    allocate(gv(natoms,3), atom_of(np))
    allocate(gv_gbf(natoms,3), gv_grid(natoms,3), gv_wt(natoms,3))
    allocate(coor_p0(3,natoms))
    Pdens = Pa + Pb
    do i = 1,natoms
       coor_p0(:,i) = atoms(i)%coor*ans2bohr
    enddo
    do i = 1,np
       atom_of(i) = (i-1)/nper + 1
    enddo
    gv = 0.0d0
    gv_gbf = 0.0d0
    gv_grid = 0.0d0
    gv_wt = 0.0d0

    atom_bf_start = 0
    do mu = 1,nconts
       atom_bf_start(bf_atom(mu)+1) = atom_bf_start(bf_atom(mu)+1) + 1
    enddo
    do A = 1,natoms
       atom_bf_start(A+1) = atom_bf_start(A+1) + atom_bf_start(A)
    enddo
    atom_bf_cursor = atom_bf_start(1:natoms)
    do mu = 1,nconts
       A = bf_atom(mu)
       atom_bf_cursor(A) = atom_bf_cursor(A) + 1
       atom_bf_list(atom_bf_cursor(A)) = mu
    enddo

    do mu = 1,nconts
       bf_rcut2_grad(mu) = (sqrt(-log(VV10_GRAD_SIG_CUTOFF)/ &
                             minval(atoms(bf_atom(mu))%shell(bf_shell(mu))%exponents)) + 2.0d0)**2
       bf_rcut2_hess(mu) = (sqrt(-log(VV10_HESS_SIG_CUTOFF)/ &
                             minval(atoms(bf_atom(mu))%shell(bf_shell(mu))%exponents)) + 2.0d0)**2
    enddo

    call prof_start("force_vv10_pass1")
    !$omp parallel do schedule(static) private(i,mu,dist2_bf,kk,nsig_grad,sig_idx_grad)
    do i = 1,np
       block
          real(8) :: val0(nconts), val1(nconts,3), wtot(nconts)
          integer :: n
          nsig_grad = 0
          do mu = 1,nconts
             dist2_bf = sum((coor(:,i)-coor_p0(:,bf_atom(mu)))**2)
             if (dist2_bf .le. bf_rcut2_grad(mu)) then
                nsig_grad = nsig_grad+1
                sig_idx_grad(nsig_grad) = mu
             endif
          enddo
          call GTOeval_point(coor(:,i), val0, val1)
          wtot = 0.0d0
          do kk = 1,nsig_grad
             n = sig_idx_grad(kk)
             wtot(:) = wtot(:) + Pdens(:,n)*val0(n)
          enddo
          rho(i) = dot_product(val0, wtot)
          grad(1,i) = 2.0d0*dot_product(val1(:,1), wtot)
          grad(2,i) = 2.0d0*dot_product(val1(:,2), wtot)
          grad(3,i) = 2.0d0*dot_product(val1(:,3), wtot)
          sig(i) = grad(1,i)**2 + grad(2,i)**2 + grad(3,i)**2
       end block
    enddo
    !$omp end parallel do
    call prof_stop("force_vv10_pass1")

    call prof_start("force_vv10_pass2")
    call vv10_evaluate(np, rho, sig, weight, coor, Enl, eps, Fn, Fg)
    call prof_stop("force_vv10_pass2")

    do A = 1,natoms
       atom_rcut2(A) = 0.0d0
       do jsh = 1,atoms(A)%nShell
          atom_rcut2(A) = max(atom_rcut2(A), rcut2_shared(jsh,A))
       enddo
    enddo
    call prof_start("force_vv10_pass3")
    !$omp parallel do schedule(dynamic,64) private(i,mu,dd,ee,A,iat,n,dist2_A,kk,nsig_grad,dist2_bf,sig_idx_grad) &
    !$omp&   private(nsig_hess,sig_idx_hess) &
    !$omp&   reduction(+:gv_gbf)
    do i = 1,np
       block
          real(8) :: val0(nconts), val1(nconts,3), Hess(nconts,6)
          real(8) :: wtot(nconts), C(nconts,3), TCV(nconts), hdt(nconts,3), C_A(nconts,3)
          real(8) :: T1, T2, T3, sn, sh, sc
          nsig_grad = 0
          nsig_hess = 0
          do mu = 1,nconts
             dist2_bf = sum((coor(:,i)-coor_p0(:,bf_atom(mu)))**2)
             if (dist2_bf .le. bf_rcut2_grad(mu)) then
                nsig_grad = nsig_grad+1
                sig_idx_grad(nsig_grad) = mu
             endif
             if (dist2_bf .le. bf_rcut2_hess(mu)) then
                nsig_hess = nsig_hess+1
                sig_idx_hess(nsig_hess) = mu
             endif
          enddo
          call GTOeval_point(coor(:,i), val0, val1)
          call Hessian_at_point_subset(coor(:,i), sig_idx_hess, nsig_hess, Hess)
          wtot = 0.0d0
          C = 0.0d0
          do kk = 1,nsig_grad
             n = sig_idx_grad(kk)
             block
                real(8) :: v0n, v1n1, v1n2, v1n3, pdmn
                v0n = val0(n); v1n1 = val1(n,1); v1n2 = val1(n,2); v1n3 = val1(n,3)
                do mu = 1,nconts
                   pdmn = Pdens(mu,n)
                   wtot(mu) = wtot(mu) + pdmn*v0n
                   C(mu,1) = C(mu,1) + pdmn*v1n1
                   C(mu,2) = C(mu,2) + pdmn*v1n2
                   C(mu,3) = C(mu,3) + pdmn*v1n3
                enddo
             end block
          enddo
          T1 = 2.0d0*Fg(i)*grad(1,i)
          T2 = 2.0d0*Fg(i)*grad(2,i)
          T3 = 2.0d0*Fg(i)*grad(3,i)
          do mu = 1,nconts
             TCV(mu) = T1*val1(mu,1) + T2*val1(mu,2) + T3*val1(mu,3)
          enddo
          hdt = 0.0d0
          do kk = 1,nsig_hess
             mu = sig_idx_hess(kk)
             do dd = 1,3
                hdt(mu,dd) = T1*Hess(mu,PK(1,dd)) + T2*Hess(mu,PK(2,dd)) + &
                             T3*Hess(mu,PK(3,dd))
             enddo
          enddo
          iat = atom_of(i)
          do A = 1,natoms
             if (A .ne. iat) then
                dist2_A = sum((coor(:,i)-coor_p0(:,A))**2)
                if (dist2_A .gt. atom_rcut2(A)) cycle
             endif
             do kk = 1,nsig_grad
                mu = sig_idx_grad(kk)
                C_A(mu,1) = 0.0d0
                C_A(mu,2) = 0.0d0
                C_A(mu,3) = 0.0d0
             enddo
             do n = atom_bf_start(A)+1, atom_bf_start(A+1)
                block
                   integer :: nn
                   real(8) :: v1n1, v1n2, v1n3, pdmn2
                   nn = atom_bf_list(n)
                   v1n1 = val1(nn,1); v1n2 = val1(nn,2); v1n3 = val1(nn,3)
                   do kk = 1,nsig_grad
                      mu = sig_idx_grad(kk)
                      pdmn2 = Pdens(mu,nn)
                      C_A(mu,1) = C_A(mu,1) + pdmn2*v1n1
                      C_A(mu,2) = C_A(mu,2) + pdmn2*v1n2
                      C_A(mu,3) = C_A(mu,3) + pdmn2*v1n3
                   enddo
                end block
             enddo
             do dd = 1,3
                sn = 0.0d0; sh = 0.0d0; sc = 0.0d0
                if (A .eq. iat) then
                   do mu = 1,nconts
                      if (bf_atom(mu) .eq. iat) cycle
                      sn = sn + Fn(i)*wtot(mu)*val1(mu,dd)
                      sh = sh + wtot(mu)*hdt(mu,dd)
                   enddo
                   do kk = 1,nsig_grad
                      mu = sig_idx_grad(kk)
                      sc = sc + (C(mu,dd) - C_A(mu,dd))*TCV(mu)
                   enddo
                   gv_gbf(A,dd) = gv_gbf(A,dd) + 2.0d0*weight(i)*(sn + sh + sc)
                else
                   do n = atom_bf_start(A)+1, atom_bf_start(A+1)
                      mu = atom_bf_list(n)
                      sn = sn + Fn(i)*wtot(mu)*val1(mu,dd)
                      sh = sh + wtot(mu)*hdt(mu,dd)
                   enddo
                   do kk = 1,nsig_grad
                      mu = sig_idx_grad(kk)
                      sc = sc + C_A(mu,dd)*TCV(mu)
                   enddo
                   gv_gbf(A,dd) = gv_gbf(A,dd) - 2.0d0*weight(i)*(sn + sh + sc)
                endif
             enddo
          enddo
       end block
    enddo
    !$omp end parallel do
    call prof_stop("force_vv10_pass3")

    call prof_start("force_vv10_grid_kernel")
    call vv10_grid_force(np, rho, sig, weight, coor, nper, natoms, gv_grid)
    call prof_stop("force_vv10_grid_kernel")

    beta_v = vv10_beta_val()
    call prof_start("force_vv10_pass5")
    call vv10_nlc_weight_gradient(np, coor, nint(atom_of), nper, rho, eps, weight, beta_v, gv_wt)
    call prof_stop("force_vv10_pass5")

    gv = gv_gbf + gv_grid + gv_wt

    do iat = 1,natoms
       atoms(iat)%atmForce(1) = atoms(iat)%atmForce(1) + gv(iat,1)
       atoms(iat)%atmForce(2) = atoms(iat)%atmForce(2) + gv(iat,2)
       atoms(iat)%atmForce(3) = atoms(iat)%atmForce(3) + gv(iat,3)
    enddo

    deallocate(coor, weight, rho, grad, sig, eps, Fn, Fg, gv, atom_of)
    deallocate(gv_gbf, gv_grid, gv_wt)
    deallocate(coor_p0)
end subroutine calc_force_vv10

subroutine Hessian_at_point_subset(coor_p, idx_list, n_idx, Hess)
    use GRID_info, only: rcut2_shared
    implicit none
    real(8),intent(in) :: coor_p(3)
    integer,intent(in) :: n_idx, idx_list(n_idx)
    real(8),intent(out) :: Hess(nconts,6)
    real(8) :: coor_rel(3), val6(6)
    integer :: kk, l, i, j, kk2
    Hess = 0.0d0
    do kk = 1,n_idx
       l = idx_list(kk)
       i = bf_atom(l)
       j = bf_shell(l)
       coor_rel(1) = coor_p(1) - atoms(i)%coor(1)*ans2bohr
       coor_rel(2) = coor_p(2) - atoms(i)%coor(2)*ans2bohr
       coor_rel(3) = coor_p(3) - atoms(i)%coor(3)*ans2bohr
       if (sum(coor_rel**2) .gt. rcut2_shared(j,i)) cycle
       do kk2 = 1,bf_nterm(l)
          call GTOHess(coor_rel, bf_term_a(kk2,l), bf_term_b(kk2,l), bf_term_c(kk2,l), &
                       atoms(i)%shell(j)%nGauss, atoms(i)%shell(j)%cnVal(:,bf_term_idim(kk2,l)), &
                       atoms(i)%shell(j)%exponents, val6)
          Hess(l,:) = Hess(l,:) + bf_term_coef(kk2,l)*val6
       enddo
    enddo
end subroutine Hessian_at_point_subset

end subroutine

