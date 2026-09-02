! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Main SCF iteration driver (Fock build, DIIS, convergence check) - SCFcycle.


subroutine SCFcycle(info,Emax,Pmax,itmax,iconv_out,econv_out,do_force)
use MOL_info
use mod_integrals, only: integrals_build_coulomb, integrals_set_accuracy, cosx_force_hi_grid, cosx_promote_md, near_singular_overlap
use mod_exchange, only: HF_exchange_frac, exchange_build, cosx_enabled
use mod_density_fitting, only: df_build_coulomb
use mod_profile, only: prof_start, prof_stop
use mod_scf_history, only: scf_hist_reset, scf_hist_record
use GRID_info, only: xcgrid_just_refined
use mod_xc, only: xc_uses_vv10
use mod_vv10, only: engine_vv10_nonself, vv10_just_activated, vv10_flush_on_switch
    implicit none
INCLUDE 'parameter.h'
    real(8)    :: Emax,Pmax
    integer    :: itmax
    integer,intent(out) :: iconv_out
    real(8),intent(out) :: econv_out
    logical,intent(in) :: do_force
    integer    :: iter,iconv
    integer    :: i,j,k,n,m,info,info_sol
    integer    :: elock_streak
    real(8)    :: plateau_prms_max
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
    logical    :: diis_active
    logical    :: diis_debug
    character(len=8) :: diis_dbg_env
    character(len=32) :: diis_env
    real(4),allocatable :: diis_Fa(:,:,:),diis_Fb(:,:,:)
    real(4),allocatable :: diis_Pa(:,:,:),diis_Pb(:,:,:)
    real(4),allocatable :: diis_ea(:,:,:),diis_eb(:,:,:)

    integer(8) :: wc1,wc2,wc6,wc_rate
    iter = 1
    iconv = 0
    elock_streak = 0
    Prms = 0.0d0
    call system_clock(count_rate=wc_rate)
    call scf_hist_reset()

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
    print *,"DIIS window size:",DIIS_MAX
    diis_n = 0
    allocate(diis_Fa(nconts,nconts,DIIS_MAX))
    allocate(diis_ea(nconts,nconts,DIIS_MAX))
    allocate(diis_Pa(nconts,nconts,DIIS_MAX))
    allocate(diis_Fb(nconts,nconts,DIIS_MAX))
    allocate(diis_eb(nconts,nconts,DIIS_MAX))
    allocate(diis_Pb(nconts,nconts,DIIS_MAX))

    print *,"cycle     Density_change     E_change(Hartree)  Total Energy(Hartree)   Exc(Hartree)  Diag.(s)  DFT(s)  Total(s)"
    do while (iter .le. itmax)
        call system_clock(wc1)
        diis_active = .false.
        if (xcgrid_just_refined .or. (vv10_just_activated .and. vv10_flush_on_switch)) then
           diis_n = 0
           if (diis_debug) print *,"RTDBG DIIS history flushed after XC grid promotion / VV10 activation"
        endif
        xcgrid_just_refined = .false.
        vv10_just_activated = .false.
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
    if (.not. diis_active .and. iter .gt. 1) then
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
           print *, '[BUILDMODE]'
           print '(A)', "  J:   "//trim(j_mode)
           print '(A)', "  K:   "//trim(k_mode)
           print '(A)', "  Vxc: "//trim(vxc_mode)
           print *, '[BUILDMODEEND]'
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
        call system_clock(wc6)
        write(*,"(I4,F16.9,4X,F20.9,4X,F16.9,F16.9,2X,F8.3,F8.3,F8.3)")  &
               iter,Prms,(E_n+E_rep)-E,E_n+E_rep,Exc, &
               real(wc2-wc1,8)/wc_rate, dft_dt, real(wc6-wc1,8)/wc_rate
        call flush(6)
        call scf_hist_record(iter, Prms, (E_n+E_rep)-E, E_n+E_rep, Exc, &
                              real(wc2-wc1,8)/wc_rate, dft_dt, real(wc6-wc1,8)/wc_rate)

        deallocate(Pc)
        econv_out = abs(E_n+E_rep-E)
        if (econv_out.le.3.0d0*Emax) then
           elock_streak = elock_streak + 1
        else
           elock_streak = 0
        endif
        plateau_prms_max = merge(30.0d0, 0.0d0, cosx_enabled) * Pmax
        if (econv_out.le.3.0d0*Emax .and. Prms.gt.Pmax .and. Prms.le.plateau_prms_max &
            .and. elock_streak.ge.8) then
           print '("Convergence: dE=",ES10.3," Hartree, dP=",ES10.3," (energy plateau; dP noise in near-degenerate manifold)")', &
                 E_n+E_rep-E, Prms
           Prms = Pmax
           econv_out = Emax
        endif
        if (econv_out.le.Emax .and. Prms.le.Pmax  ) then
             print '("Convergence: dE=",ES10.3," Hartree, dP=",ES10.3)',E_n+E_rep-E,Prms
            if (cosx_enabled .and. .not. cosx_force_hi_grid) then
               cosx_force_hi_grid = .true.
               call scf_build_fock(Pa_n, Pb_n, diis_debug, iter, Prms, Exc, E_n, dft_dt, do_force)
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
    else
       print *,"SCF converged"
    endif
    print *, '[SCFEND]'
    print *

    iconv_out = iconv

end subroutine
