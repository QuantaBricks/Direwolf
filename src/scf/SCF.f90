! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Main SCF iteration driver (Fock build, DIIS, convergence check) - SCFcycle.


subroutine SCFcycle(info,Emax,Pmax,itmax,iconv_out,econv_out,do_force)
use MOL_info
use mod_integrals, only: integrals_build_coulomb, integrals_set_accuracy, cosx_force_hi_grid, &
                         cosx_promote_md, near_singular_overlap, cosx_set_force_kscreen, &
                         incr_since_full, incr_full_period, cosx_all_hi, cosx_all_hi_force
use mod_exchange, only: HF_exchange_frac, exchange_build, cosx_enabled
use mod_density_fitting, only: df_build_coulomb
use mod_profile, only: prof_start, prof_stop
use mod_output_format, only: set_buildmode, print_buildmode
use mod_scf_history, only: scf_hist_reset, scf_hist_record
use GRID_info, only: xcgrid_just_refined
use mod_xc, only: xc_uses_vv10
use mod_vv10, only: engine_vv10_nonself, vv10_just_activated, vv10_flush_on_switch
use mod_soscf, only: soscf_reset, soscf_step, newton_cg_step
    implicit none
INCLUDE 'parameter.h'
    real(8)    :: Emax,Pmax
    integer    :: itmax
    integer,intent(out) :: iconv_out
    real(8),intent(out) :: econv_out
    logical,intent(in) :: do_force
    integer    :: iter,iconv
    integer    :: i,j,k,n,m,info,info_sol
    real(8)    :: E_n
    real(8)    :: Exc
    real(8)    :: Prms
    real(8),parameter :: COSX_MD_PRMS = 1.0d-4
    real(8),parameter :: COSX_HI_PRMS = 2.0d-5
    real(8),allocatable    :: Pa_n(:,:), Pb_n(:,:),Pc(:,:)
    real(8),allocatable    :: Pa_old(:,:), Pb_old(:,:)
    real(8) :: dft_dt
    real(8) :: Enl_postscf

    integer    :: DIIS_MAX,diis_n
    logical    :: hi_frozen, hi_rebase, fz_phase, hi_stage_done
    real(8)    :: Emax_eff, Pmax_eff, cosx_lo_pmax
    logical    :: hi_reconv
    logical,parameter :: COSX_ALL_HI_RUN = .false., COSX_HI_RECONVERGE = .false.
    logical,parameter :: COSX_HI_FLUSH = .false., COSX_HI_FROZEN = .false.
    real(8)    :: dE_fz, E_lo_conv
    real(8), allocatable :: dFa_fz(:,:), dFb_fz(:,:)
    real(8),allocatable :: dFa_hi(:,:),dFb_hi(:,:)
    logical    :: diis_active
    logical    :: diis_debug
    character(len=8) :: diis_dbg_env
    character(len=32) :: diis_env
    real(8),allocatable :: diis_Fa(:,:,:),diis_Fb(:,:,:)
    real(4),allocatable :: diis_Pa(:,:,:),diis_Pb(:,:,:)
    real(4),allocatable :: diis_ea(:,:,:),diis_eb(:,:,:)

    integer(8) :: wc1,wc2,wc6,wc_rate
    logical :: soscf_enable, soscf_active, newton_enable
    integer :: newton_start_iter
    logical :: newton_auto
    integer :: stab_mode, stab_rounds
    logical :: stab_pending, stab_took
    real(8) :: stab_lam
    character(len=16) :: stab_env
    real(8), allocatable :: prms_h(:)
    character(len=32) :: newton_iter_env
    character(len=32) :: newton_env, newton_start_env
    real(8) :: newton_grad_max
    character(len=8) :: soscf_env
    character(len=32) :: soscf_start_env
    real(8) :: soscf_grad_max
    iter = 1
    iconv = 0
    Prms = 0.0d0
    call system_clock(count_rate=wc_rate)
    call scf_hist_reset()

    soscf_env = ""
    call get_environment_variable("ENGINE_SOSCF", soscf_env)
    soscf_enable = (trim(soscf_env) .eq. "1")
    soscf_grad_max = 1.0d-2
    newton_env = ""
    call get_environment_variable("ENGINE_NEWTON", newton_env)
    newton_enable = (trim(newton_env) .eq. "1")
    newton_grad_max = 5.0d-3
    newton_start_iter = 64
    newton_auto = .false.
    stab_mode = 0; stab_rounds = 0; stab_pending = .false.
    call get_environment_variable("ENGINE_STABILITY", stab_env)
    if (len_trim(stab_env) .gt. 0) read(stab_env,*) stab_mode
    allocate(prms_h(0:max(itmax,1)+1)); prms_h = 0.0d0
    call get_environment_variable("ENGINE_NEWTON_ITER", newton_iter_env)
    if (len_trim(newton_iter_env) .gt. 0) read(newton_iter_env,*) newton_start_iter
    if (newton_enable) newton_start_iter = 3
    call get_environment_variable("ENGINE_NEWTON_START", newton_start_env)
    if (len_trim(newton_start_env) .gt. 0) read(newton_start_env,*) newton_grad_max
    call get_environment_variable("ENGINE_SOSCF_START", soscf_start_env)
    if (len_trim(soscf_start_env) .gt. 0) read(soscf_start_env,*) soscf_grad_max
    call soscf_reset()

    diis_dbg_env = ""
    call get_environment_variable("ENGINE_DIIS_DEBUG", diis_dbg_env)
    diis_debug = (trim(diis_dbg_env) .eq. "1")
    DIIS_MAX = 10
    do i = 1, Natoms
       associate (z => atoms(i)%charge)
       if ((z.ge.21 .and. z.le.30) .or. (z.ge.39 .and. z.le.48) .or. &
           (z.ge.57 .and. z.le.80) .or. (z.ge.89 .and. z.le.112)) DIIS_MAX = 15
       end associate
    enddo
    call get_environment_variable("ENGINE_DIIS_WINDOW", diis_env)
    if (len_trim(diis_env) .gt. 0) read(diis_env,*) DIIS_MAX
    if (DIIS_MAX .lt. 1) DIIS_MAX = 1
    if (.not. engine_quiet) print *,"DIIS window size:",DIIS_MAX
    diis_n = 0
    hi_rebase = .false.
    fz_phase = .false.
    hi_stage_done = .false.
    dE_fz = 0
    cosx_lo_pmax = 1.0d-5
    cosx_all_hi = cosx_enabled .and. COSX_ALL_HI_RUN
    if (cosx_all_hi) then
       cosx_force_hi_grid = .true.
       hi_stage_done = .true.
       cosx_all_hi_force = do_force
       if (do_force) call cosx_set_force_kscreen(.true.)
    endif
    hi_reconv = COSX_HI_RECONVERGE
    hi_frozen = COSX_HI_FROZEN
    allocate(diis_Fa(nconts,nconts,DIIS_MAX))
    allocate(diis_ea(nconts,nconts,DIIS_MAX))
    allocate(diis_Pa(nconts,nconts,DIIS_MAX))
    allocate(diis_Fb(nconts,nconts,DIIS_MAX))
    allocate(diis_eb(nconts,nconts,DIIS_MAX))
    allocate(diis_Pb(nconts,nconts,DIIS_MAX))

    if (.not. engine_quiet) &
       print *,"cycle     Density_change     E_change(Hartree)  Total Energy(Hartree)   Exc(Hartree)  Diag.(s)  DFT(s)  Total(s)"
       call engine_startup_mark("SCF header")
    do while (iter .le. itmax)
        call system_clock(wc1)
        diis_active = .false.
        if (xcgrid_just_refined .or. (vv10_just_activated .and. vv10_flush_on_switch)) then
           diis_n = 0
           if (diis_debug) print *,"RTDBG DIIS history flushed after XC grid promotion / VV10 activation"
        endif
        xcgrid_just_refined = .false.
        vv10_just_activated = .false.
        soscf_active = .false.
        if (stab_pending) then
           soscf_active = .true.
           stab_pending = .false.
           diis_n = 0
        endif
        if (.not. soscf_active .and. newton_start_iter .gt. 0 .and. iter .ge. newton_start_iter) then
           call prof_start("diag")
           call newton_cg_step(C_a, C_b, eLev_a, eLev_b, n_alpha, n_beta, Prms, newton_grad_max, 0.3d0, 12, soscf_active)
           call prof_stop("diag")
           if (soscf_active) then
              diis_n = 0
              if (.not. engine_quiet) print *, "ENGINE_NEWTON: Newton-CG orbital step"
           endif
        endif
        if (soscf_enable .and. iter .gt. 2 .and. .not. soscf_active) then
           call prof_start("diag")
           call soscf_step(nconts, n_alpha, Fa, S, C_a, eLev_a, soscf_grad_max, soscf_active, 1)
           if (soscf_active .and. Multi .ne. 1) then
              call soscf_step(nconts, n_beta, Fb, S, C_b, eLev_b, 1.0d30, soscf_active, 2)
           endif
           call prof_stop("diag")
           if (soscf_active) then
              if (Multi .eq. 1) then
                 C_b = C_a
                 eLev_b = eLev_a
              endif
              if (diis_n .ne. 0) then
                 if (.not. engine_quiet) print *, &
                    "ENGINE_SOSCF: orbital gradient small enough, switching from DIIS to SOSCF orbital rotation"
                 diis_n = 0
              endif
           endif
        endif
        if (.not. soscf_active) then
           if (iter .gt. 1) then
              call prof_start("diis_extrapolate")
              call diis_extrapolate(nconts, multi, iter, Fa, Fb, Pa, Pb, S, X, &
                                     DIIS_MAX, diis_n, diis_Fa, diis_Fb, diis_Pa, diis_Pb, &
                                     diis_ea, diis_eb, diis_debug, diis_active)
              call prof_stop("diis_extrapolate")
           endif
           call prof_start("diag")
           call solHFR_KS(info_sol, Prms, near_singular_overlap)
           call prof_stop("diag")
        endif
        call system_clock(wc2)

        call prof_start("density_build")
        allocate(Pa_n(nconts,nconts))
        allocate(Pb_n(nconts,nconts))
        Pa_n = 0
        Pb_n = 0
        do i = 1,nconts
           do j = i,nconts
              do k = 1,n_alpha
                 Pa_n(i,j) = C_a(i,k)*C_a(j,k) + Pa_n(i,j)
              enddo
                 Pa_n(j,i) = Pa_n(i,j)
           enddo
        enddo
        if (multi .eq.1) then
           Pb_n = Pa_n
        else
        do i = 1,nconts
           do j = i,nconts
              do k = 1,n_beta
                 Pb_n(i,j) = C_b(i,k)*C_b(j,k) + Pb_n(i,j)
              enddo
                 Pb_n(j,i) = Pb_n(i,j)
           enddo
        enddo

        endif
        call prof_stop("density_build")
    if (.not. diis_active .and. .not. soscf_active .and. iter .gt. 1) then
       Pa_n = Pa_n*0.5 + Pa*0.5
       Pb_n = Pb_n*0.5 + Pb*0.5
    endif
    allocate(Pa_old(nconts,nconts),Pb_old(nconts,nconts))
    Pa_old = Pa
    Pb_old = Pb
    Pa = Pa_n
    Pb = Pb_n
    Fa = 0
    Fb = 0
        call scf_build_fock(Pa_n, Pb_n, diis_debug, iter, Prms, Exc, E_n, dft_dt, do_force)
        if (fz_phase) then
           Fa = Fa + dFa_fz
           if (multi .ne. 1) Fb = Fb + dFb_fz
           E_n = E_n + dE_fz
        endif
        if (iter .eq. 1) then
           block
           use mod_integrals, only: engine_use_df_j, engine_use_df_k, direct_mode, df_direct_mode
           use GRID_info, only: grid_cache_mode
           character(len=32) :: j_mode, k_mode, vxc_mode
           if (engine_use_df_j) then
              j_mode = "DF, "//trim(merge("DIRECT","STORE ",df_direct_mode))
           else
              j_mode = "exact, "//trim(merge("DIRECT","STORE ",direct_mode))
           endif
           if (cosx_enabled) then
              k_mode = "COSX (grid, always recomputed)"
           else if (HF_exchange_frac .le. 0.0d0) then
              k_mode = "N/A (no exact exchange)"
           else if (engine_use_df_k) then
              k_mode = "DF, "//trim(merge("DIRECT","STORE ",df_direct_mode))
           else
              k_mode = "exact, "//trim(merge("DIRECT","STORE ",direct_mode))
           endif
           vxc_mode = trim(merge("STORE (grid cached)      ","DIRECT (recomputed/batch)",grid_cache_mode))
           if (.not. engine_quiet) call set_buildmode(j_mode, k_mode, vxc_mode)
           end block
        endif
        allocate(Pc(nconts,nconts))
        Pc = (Pa_n+Pb_n)-(Pa_old+Pb_old)
        deallocate(Pa_old,Pb_old)

        Prms = 0
        do i = 1,nconts
           do j = i,nconts
              Prms = Pc(i,j)**2 + Prms
           enddo
        enddo

        Prms = (Prms/((nconts+1)*nconts/2))**0.5
        if (iter .le. itmax) prms_h(iter) = Prms
        call system_clock(wc6)
        if (.not. engine_quiet) then
           write(*,"(I4,F16.9,4X,F20.9,4X,F16.9,F16.9,2X,F8.3,F8.3,F8.3)")  &
                  iter,Prms,(E_n+E_rep)-E,E_n+E_rep,Exc, &
                  real(wc2-wc1,8)/wc_rate, dft_dt, real(wc6-wc1,8)/wc_rate
           call flush(6)
        endif
        call scf_hist_record(iter, Prms, (E_n+E_rep)-E, E_n+E_rep, Exc, &
                              real(wc2-wc1,8)/wc_rate, dft_dt, real(wc6-wc1,8)/wc_rate)

        deallocate(Pc)
        econv_out = abs(E_n+E_rep-E)
        if (soscf_active .and. (E_n+E_rep-E) .gt. 1.0d-6) call soscf_reset()
        Emax_eff = Emax
        Pmax_eff = Pmax
        if (cosx_enabled .and. do_force .and. hi_reconv .and. .not. hi_stage_done .and. &
            .not. cosx_force_hi_grid .and. .not. fz_phase .and. cosx_lo_pmax .gt. 0.0d0) then
           Pmax_eff = max(Pmax, cosx_lo_pmax)
           Emax_eff = max(Emax, 1.0d-2*cosx_lo_pmax)
        endif
        if (econv_out.le.Emax_eff .and. Prms.le.Pmax_eff) then
             if (.not. engine_quiet) &
                print '("Convergence: dE=",ES10.3," Hartree, dP=",ES10.3)',E_n+E_rep-E,Prms
            if (cosx_enabled .and. fz_phase) then
               fz_phase = .false.
               hi_stage_done = .true.
               cosx_force_hi_grid = .true.
               call cosx_set_force_kscreen(.true.)
               allocate(dFa_hi(nconts,nconts),dFb_hi(nconts,nconts))
               dFa_hi = -Fa
               dFb_hi = -Fb
               incr_since_full = incr_full_period
               call scf_build_fock(Pa_n, Pb_n, diis_debug, iter, Prms, Exc, E_n, dft_dt, do_force)
               dFa_hi = dFa_hi + Fa
               dFb_hi = dFb_hi + Fb
               call diis_shift_history(nconts, multi, dFa_hi, dFb_hi, S, X, diis_n, DIIS_MAX, &
                                       diis_Fa, diis_Fb, diis_Pa, diis_Pb, diis_ea, diis_eb)
               deallocate(dFa_hi,dFb_hi,dFa_fz,dFb_fz)
               E = E_n+E_rep
               Pa = Pa_n
               Pb = Pb_n
               deallocate(Pa_n)
               deallocate(Pb_n)
               iter = iter + 1
               cycle
            endif
            if (cosx_enabled .and. .not. hi_stage_done .and. .not. cosx_force_hi_grid) then
               cosx_force_hi_grid = .true.
               if (do_force) call cosx_set_force_kscreen(.true.)
               hi_rebase = do_force .and. hi_reconv .and. diis_n .gt. 0 .and. .not. COSX_HI_FLUSH
               if (hi_rebase) then
                  allocate(dFa_hi(nconts,nconts),dFb_hi(nconts,nconts))
                  dFa_hi = -Fa
                  dFb_hi = -Fb
               endif
               E_lo_conv = E_n
               call scf_build_fock(Pa_n, Pb_n, diis_debug, iter, Prms, Exc, E_n, dft_dt, do_force)
               if (hi_rebase) then
                  dFa_hi = dFa_hi + Fa
                  dFb_hi = dFb_hi + Fb
               endif
               if (.not. hi_reconv) goto 7711
               if (do_force) then
                  if (hi_rebase) then
                     call diis_shift_history(nconts, multi, dFa_hi, dFb_hi, S, X, diis_n, DIIS_MAX, &
                                             diis_Fa, diis_Fb, diis_Pa, diis_Pb, diis_ea, diis_eb)
                     if (hi_frozen) then
                        allocate(dFa_fz(nconts,nconts),dFb_fz(nconts,nconts))
                        dFa_fz = dFa_hi
                        dFb_fz = dFb_hi
                        dE_fz = E_n - E_lo_conv
                        fz_phase = .true.
                        cosx_force_hi_grid = .false.
                        call cosx_set_force_kscreen(.false.)
                     endif
                     deallocate(dFa_hi,dFb_hi)
                  else
                     diis_n = 0
                  endif
                  E = E_n+E_rep
                  Pa = Pa_n
                  Pb = Pb_n
                  deallocate(Pa_n)
                  deallocate(Pb_n)
                  iter = iter + 1
                  cycle
               endif
            endif
7711        continue
            if (stab_mode .gt. 0 .and. stab_rounds .lt. 3) then
               stab_rounds = stab_rounds + 1
               Pa = Pa_n
               Pb = Pb_n
               call newton_cg_step(C_a, C_b, eLev_a, eLev_b, n_alpha, n_beta, Prms, 1.0d30, 0.5d0, 12, &
                                   stab_took, stab_mode, stab_lam)
               if (stab_took) then
                  print *, "STABILITY: rotated orbitals along the unstable mode, re-converging SCF"
                  stab_pending = .true.
                  diis_n = 0
                  E = E_n+E_rep
                  deallocate(Pa_n)
                  deallocate(Pb_n)
                  iter = iter + 1
                  cycle
               endif
            endif
            iconv = 1
            E = E_n+E_rep
            Pa = Pa_n
            Pb = Pb_n
            deallocate(Pa_n)
            deallocate(Pb_n)
            if (xc_uses_vv10() .and. engine_vv10_nonself) then
               call DFT_vv10_postscf_aca(Enl_postscf)
               E = E + Enl_postscf
            endif
            exit
        endif

        E = E_n+E_rep
        Pa = Pa_n
        Pb = Pb_n
        deallocate(Pa_n)
        deallocate(Pb_n)
        iter = iter +1
    enddo

    deallocate(diis_Fa,diis_ea,diis_Pa)
    deallocate(diis_Fb,diis_eb,diis_Pb)

    if (iconv .eq. 0) then
       print *,"SCF FAILED"
    else if (.not. engine_quiet) then
       print *,"SCF converged"
    endif
    if (.not. engine_quiet) then
       print *, '[SCFEND]'
       print *
       call print_buildmode()
    endif

    iconv_out = iconv

end subroutine
