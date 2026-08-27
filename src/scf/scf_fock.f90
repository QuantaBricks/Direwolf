! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! scf_fock: builds this iteration's Fock matrices (Fa/Fb, MOL_info module


subroutine scf_build_fock(Pa_n, Pb_n, diis_debug, iter, Prms, Exc, E_n, dft_dt, need_deriv)
use MOL_info
use mod_integrals, only: integrals_build_coulomb, integrals_set_accuracy, &
                        integrals_build_jk_fused, jk_fusion_active, &
                          engine_use_df_j, engine_use_df_k
use GRID_info, only: xc_sig_now, XC_SIG_TIGHT, XC_SIG_LOOSE
use GRID_info, only: xc_geo_level, XC_NLEVEL, XC_LEVEL_PRMS
use GRID_info, only: xcgrid_dynamic, xcgrid_refined, xcgrid_switch_prms
use mod_vv10, only: vv10_dynamic, vv10_active_now, vv10_just_activated, vv10_switch_prms
use mod_exchange, only: HF_exchange_frac, RS_omega, RS_beta, exchange_build, exchange_build_lr, &
                         exchange_build_cosx_lr, exchange_build_cosx_sr, &
                         exchange_build_df, exchange_build_df_lr, exchange_build_cosx, cosx_enabled
use mod_density_fitting, only: df_build_coulomb
use mod_profile, only: prof_start, prof_stop
use mod_cosmo, only: cosmo_enabled, cosmo_scf_step
use mod_scf_history, only: scf_hist_set_ecomponents
implicit none
integer,intent(in) :: iter
real(8),intent(in) :: Pa_n(nconts,nconts),Pb_n(nconts,nconts)
logical,intent(in) :: diis_debug
real(8),intent(in) :: Prms
real(8),intent(out) :: Exc,E_n
real(8),intent(out) :: dft_dt
logical,intent(in) :: need_deriv

integer :: i,j,info,jlv
integer(8) :: wc3,wc4,wc_rate
real(8),allocatable :: Fxc_a(:,:),Fxc_b(:,:)
real(8),allocatable :: Fra(:,:),Frb(:,:)
real(8) :: ecoul_acc, exact_exchange_acc
real(8),allocatable :: Ka(:,:),Kb(:,:)
logical :: jk_fused_this_iter
real(8),allocatable :: Ptot(:,:),Jmat(:,:)
real(8),allocatable :: Fock_cosmo(:,:)
real(8) :: E_cosmo_nuc

allocate(Fxc_a(nconts,nconts))
allocate(Fxc_b(nconts,nconts))
allocate(Fra(nconts,nconts))
allocate(Frb(nconts,nconts))
if (xcgrid_dynamic .and. .not. xcgrid_refined) then
   if (Prms .gt. 0.0d0 .and. Prms_round8(Prms) .lt. xcgrid_switch_prms) call xcgrid_refine()
endif
if (vv10_dynamic .and. .not. vv10_active_now) then
   if (Prms .gt. 0.0d0 .and. Prms_round8(Prms) .lt. vv10_switch_prms) then
      vv10_active_now = .true.
      vv10_just_activated = .true.
      print *,"VV10: activating self-consistent NLC term (dP <",Prms,")"
      call flush(6)
   endif
endif
call system_clock(count_rate=wc_rate)
call system_clock(wc3)
if (diis_debug) print *, "RTDBG entering DFT_calc iter=",iter
if (diis_debug) call flush(6)
call prof_start("dft_xc")
if (HF_exchange_frac .ge. 1.0d0) then
   Exc = 0.0d0
   Fxc_a = 0.0d0
   Fxc_b = 0.0d0
else
   call  DFT_calc(Exc,Fxc_a,Fxc_b,info,need_deriv)
endif
call prof_stop("dft_xc")
call system_clock(wc4)
dft_dt = real(wc4-wc3,8)/wc_rate
if (diis_debug) print *, "RTDBG DFT_calc done, dt=",dft_dt
if (diis_debug) call flush(6)
allocate(Ptot(nconts,nconts))
allocate(Jmat(nconts,nconts))
Ptot = Pa_n + Pb_n
call integrals_set_accuracy(Prms)
if (Prms .gt. 0.0d0) then
   xc_sig_now = min(XC_SIG_LOOSE, max(XC_SIG_TIGHT, Prms*1.0d-3))
else
   xc_sig_now = XC_SIG_TIGHT
endif
xc_geo_level = XC_NLEVEL
do jlv = 1,XC_NLEVEL-1
   if (Prms .gt. XC_LEVEL_PRMS(jlv)) then
      xc_geo_level = jlv
      exit
   endif
enddo
call prof_start("coulomb_build")
jk_fused_this_iter = .false.
if (engine_use_df_j) then
   call df_build_coulomb(nconts, Ptot, Jmat)
else if (jk_fusion_active() .and. .not. cosx_enabled &
         .and. HF_exchange_frac .gt. 0.0d0) then
   if (.not. allocated(Ka)) allocate(Ka(nconts,nconts),Kb(nconts,nconts))
   call integrals_build_jk_fused(nconts, Ptot, Pa_n, Pb_n, Jmat, Ka, Kb)
   jk_fused_this_iter = .true.
else
   call integrals_build_coulomb(nconts, Ptot, Jmat)
endif
call prof_stop("coulomb_build")
if (diis_debug) print *, "RTDBG coulomb_build done"
if (diis_debug) call flush(6)
E_cosmo_nuc = 0.0d0
if (cosmo_enabled) then
   allocate(Fock_cosmo(nconts,nconts))
   call prof_start("cosmo_scf_step")
   call cosmo_scf_step(nconts, Ptot, Fock_cosmo, E_cosmo_nuc, iter)
   call prof_stop("cosmo_scf_step")
endif
exact_exchange_acc = 0.0d0
if (HF_exchange_frac .gt. 0.0d0) then
   if (.not. allocated(Ka)) allocate(Ka(nconts,nconts),Kb(nconts,nconts))
   call prof_start("exchange_build")
   if (cosx_enabled .and. iter .eq. 0) then
      Ka = 0.0d0
      Kb = 0.0d0
   else if (cosx_enabled) then
      call exchange_build_cosx(nconts, Pa_n, Pb_n, Ka, Kb, need_deriv)
   else if (engine_use_df_k) then
      call exchange_build_df(nconts, C_a, n_alpha, C_b, n_beta, Ka, Kb)
   else if (.not. jk_fused_this_iter) then
      call exchange_build(nconts, Pa_n, Pb_n, Ka, Kb)
   endif
   call prof_stop("exchange_build")
   Fra = Hcore + Jmat - HF_exchange_frac*Ka
   Frb = Hcore + Jmat - HF_exchange_frac*Kb
   if (RS_omega .gt. 0.0d0) then
      block
         real(8),allocatable :: Ka_lr(:,:), Kb_lr(:,:)
         allocate(Ka_lr(nconts,nconts),Kb_lr(nconts,nconts))
         call prof_start("exchange_build_lr")
         if (cosx_enabled .and. iter .eq. 0) then
            Ka_lr = 0.0d0
            Kb_lr = 0.0d0
         else if (cosx_enabled) then
            call exchange_build_cosx_lr(nconts, Pa_n, Pb_n, Ka_lr, Kb_lr)
         else if (engine_use_df_k) then
            call exchange_build_df_lr(nconts, C_a, n_alpha, C_b, n_beta, Ka_lr, Kb_lr)
         else
            call exchange_build_lr(nconts, Pa_n, Pb_n, Ka_lr, Kb_lr)
         endif
         call prof_stop("exchange_build_lr")
         exact_exchange_acc = exact_exchange_acc &
            - 0.5d0*RS_beta*(sum(Pa_n*Ka_lr)+sum(Pb_n*Kb_lr))
         Fra = Fra - RS_beta*Ka_lr
         Frb = Frb - RS_beta*Kb_lr
         deallocate(Ka_lr,Kb_lr)
      end block
   endif
   ecoul_acc = 0.5d0*sum((Pa_n+Pb_n)*Jmat)
   exact_exchange_acc = exact_exchange_acc &
      - 0.5d0*HF_exchange_frac*(sum(Pa_n*Ka)+sum(Pb_n*Kb))
   call scf_hist_set_ecomponents(ecoul_acc, exact_exchange_acc, &
      sum((Pa_n+Pb_n)*Hcore), E_rep)
   deallocate(Ka,Kb)
else
   ecoul_acc = 0.5d0*sum((Pa_n+Pb_n)*Jmat)
   call scf_hist_set_ecomponents(ecoul_acc, 0.0d0, &
      sum((Pa_n+Pb_n)*Hcore), E_rep)
   Fra = Hcore + Jmat
   Frb = Hcore + Jmat
endif
if (cosmo_enabled) then
   Fra = Fra + Fock_cosmo
   Frb = Frb + Fock_cosmo
   deallocate(Fock_cosmo)
endif
deallocate(Ptot)
deallocate(Jmat)
Fa = Fra + Fxc_a
Fb = Frb + Fxc_b
deallocate(Fxc_a)
deallocate(Fxc_b)

E_n = 0
do i =1,nconts
   do j =1,nconts
      E_n = E_n + (Pa_n(i,j)+Pb_n(i,j))*Hcore(i,j)+&
                  Pa_n(i,j)*Fra(i,j)+Pb_n(i,j)*Frb(i,j)
   enddo
enddo
deallocate(Fra)
deallocate(Frb)
E_n = E_n*0.5 + Exc + E_cosmo_nuc

contains

real(8) function Prms_round8(x) result(r)
implicit none
real(8),intent(in) :: x
real(8) :: scale
scale = 1.0d7/10.0d0**floor(log10(x))
r = anint(x*scale)/scale
end function Prms_round8

end subroutine scf_build_fock
