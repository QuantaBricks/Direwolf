! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_xc: the sole boundary between the rest of Engine and the

module mod_xc
use xc_f03_lib_m
use iso_c_binding, only: c_int, c_size_t, c_double
use mod_vv10, only: vv10_set
implicit none
private
public :: xc_eval
public :: xc_select_functional
public :: xc_uses_tau
public :: xc_uses_vv10

logical :: initialized = .false.
logical :: use_mgga = .false.
logical :: use_lda = .false.
logical :: use_c = .false.
logical :: use_vv10 = .false.
logical :: use_hybrid = .false.
logical :: is_rsh = .false.
logical :: is_hf = .false.
integer(c_int) :: func_x_id = XC_GGA_X_PBE
integer(c_int) :: func_c_id = XC_GGA_C_PBE
type(xc_f03_func_t) :: func_x_upol, func_c_upol, func_x_pol, func_c_pol

contains

subroutine xc_select_functional(name, hf_frac, rs_omega, rs_beta)
    implicit none
    character(len=*),intent(in) :: name
    real(8),intent(out) :: hf_frac
    real(8),intent(out) :: rs_omega, rs_beta
    integer(c_int) :: family, number, flags
    real(c_double) :: omega_c, alpha_c, beta_c
    type(xc_f03_func_info_t) :: info_x
    character(len=len(name)) :: name_upper
    integer :: ic, jc

    name_upper = trim(name)
    do ic = 1, len_trim(name)
       jc = iachar(name_upper(ic:ic))
       if (jc .ge. iachar('a') .and. jc .le. iachar('z')) &
          name_upper(ic:ic) = achar(jc - (iachar('a')-iachar('A')))
    enddo

    is_hf = (trim(name_upper) .eq. "HF")
    use_vv10 = .false.

    select case (trim(name_upper))
    case ("B3LYP","B3LYP_HYB")
       func_x_id = XC_HYB_GGA_XC_B3LYP
       func_c_id = XC_GGA_C_PBE
    case ("CAM-B3LYP")
       func_x_id = XC_HYB_GGA_XC_CAM_B3LYP
       func_c_id = XC_GGA_C_PBE
    case ("TPSS")
       func_x_id = XC_MGGA_X_TPSS
       func_c_id = XC_MGGA_C_TPSS
    case ("R2SCAN")
       func_x_id = XC_MGGA_X_R2SCAN
       func_c_id = XC_MGGA_C_R2SCAN
    case ("R2SCAN0")
       func_x_id = XC_HYB_MGGA_XC_R2SCAN0
       func_c_id = XC_MGGA_C_R2SCAN
    case ("M06-2X","M06-2X_HYB")
       func_x_id = XC_HYB_MGGA_X_M06_2X
       func_c_id = XC_MGGA_C_M06_2X
    case ("WB97X-D","WB97X_D")
       func_x_id = XC_HYB_GGA_XC_WB97X_D
       func_c_id = XC_GGA_C_PBE
    case ("WB97X")
       func_x_id = XC_HYB_GGA_XC_WB97X
       func_c_id = XC_GGA_C_PBE
    case ("PBE0")
       func_x_id = XC_HYB_GGA_XC_PBEH
       func_c_id = XC_GGA_C_PBE
    case ("B97-3C","B97_3C")
       func_x_id = XC_GGA_XC_B97_3C
       func_c_id = XC_GGA_C_PBE
    case ("PBE","PBE_PBE","GGA_X_PBE,GGA_C_PBE")
       func_x_id = XC_GGA_X_PBE
       func_c_id = XC_GGA_C_PBE
    case ("BLYP")
       func_x_id = XC_GGA_X_B88
       func_c_id = XC_GGA_C_LYP
    case ("BP86")
       func_x_id = XC_GGA_X_B88
       func_c_id = XC_GGA_C_P86
       block
         character(len=8) :: p86env
         p86env = ""
         call get_environment_variable("ENGINE_P86_VWN", p86env)
         if (trim(p86env) .eq. "1") func_c_id = XC_GGA_C_P86VWN
       end block
    case ("M06","M06_HYB")
       func_x_id = XC_HYB_MGGA_X_M06
       func_c_id = XC_MGGA_C_M06
    case ("M05-2X","M05_2X")
       func_x_id = XC_HYB_MGGA_X_M05_2X
       func_c_id = XC_MGGA_C_M05_2X
    case ("M06L","M06-L","M06_L")
       func_x_id = XC_MGGA_X_M06_L
       func_c_id = XC_MGGA_C_M06_L
    case ("MN15")
       func_x_id = XC_HYB_MGGA_X_MN15
       func_c_id = XC_MGGA_C_MN15
    case ("MN15L","MN15-L","MN15_L")
       func_x_id = XC_MGGA_X_MN15_L
       func_c_id = XC_MGGA_C_MN15_L
    case ("LDA")
       func_x_id = XC_LDA_X
       func_c_id = XC_LDA_C_VWN_RPA
    case ("VV10")
       func_x_id = XC_GGA_X_RPW86
       func_c_id = XC_GGA_C_PBE
       use_vv10 = .true.
       call vv10_set(5.9d0, 0.0093d0)
     case ("RPW86PBE","RPW86","GGA_X_RPW86,GGA_C_PBE")
        func_x_id = XC_GGA_X_RPW86
        func_c_id = XC_GGA_C_PBE
     case ("WB97M-V","WB97M_V")
        func_x_id = XC_HYB_MGGA_XC_WB97M_V
        func_c_id = XC_GGA_C_PBE
        use_vv10 = .true.
        call vv10_set(6.0d0, 0.01d0)
    case ("WB97X-V","WB97X_V")
       func_x_id = XC_HYB_GGA_XC_WB97X_V
       func_c_id = XC_GGA_C_PBE
       use_vv10 = .true.
       call vv10_set(6.0d0, 0.01d0)
    case ("WB97X-3C","WB97X_3C")
       func_x_id = XC_HYB_GGA_XC_WB97X_V
       func_c_id = XC_GGA_C_PBE
    case ("B97M-V","B97M_V")
       func_x_id = XC_MGGA_XC_B97M_V
       func_c_id = XC_GGA_C_PBE
       use_vv10 = .true.
       call vv10_set(6.0d0, 0.01d0)
    case ("HF")
       func_x_id = XC_GGA_X_PBE
       func_c_id = XC_GGA_C_PBE
    case default
       print *, "Error: unrecognized functional '"//trim(name)//"' - check spelling " &
                //"(e.g. WB97M-V, not WB97MV)."
       stop 1
    end select

    family = xc_f03_family_from_id(func_x_id, number=number)
    use_mgga   = (family .eq. XC_FAMILY_MGGA) .or. (family .eq. XC_FAMILY_HYB_MGGA)
    use_hybrid = (family .eq. XC_FAMILY_HYB_GGA) .or. (family .eq. XC_FAMILY_HYB_MGGA) &
            .or. (family .eq. XC_FAMILY_HYB_LDA)
    use_lda    = (family .eq. XC_FAMILY_LDA) .or. (family .eq. XC_FAMILY_HYB_LDA)

    initialized = .false.
    call ensure_functionals_initialized()

    rs_omega = 0.0d0
    rs_beta  = 0.0d0
    if (is_hf) then
       hf_frac = 1.0d0
    else if (.not. use_hybrid) then
       hf_frac = 0.0d0
    else
       info_x = xc_f03_func_get_info(func_x_upol)
       flags  = xc_f03_func_info_get_flags(info_x)
       is_rsh = (iand(flags, XC_FLAGS_HYB_CAM)  .ne. 0) .or. (iand(flags, XC_FLAGS_HYB_CAMY) .ne. 0) &
           .or. (iand(flags, XC_FLAGS_HYB_LC)   .ne. 0) .or. (iand(flags, XC_FLAGS_HYB_LCY)  .ne. 0)
       if (is_rsh) then
          call xc_f03_hyb_cam_coef(func_x_upol, omega_c, alpha_c, beta_c)
          hf_frac  = alpha_c + beta_c
          rs_omega = omega_c
          rs_beta  = -beta_c
          block
             character(len=8) :: rshdbg
             call get_environment_variable('ENGINE_RSH_DEBUG', rshdbg)
             if (trim(rshdbg) .eq. '1') then
                write(6,'(a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6)') &
                   ' RSHDEBUG omega_c=',omega_c,' alpha_c=',alpha_c,' beta_c=',beta_c, &
                   ' hf_frac=',hf_frac,' rs_beta=',rs_beta
             endif
          end block
       else
          hf_frac = xc_f03_hyb_exx_coef(func_x_upol)
       endif
    endif
end subroutine xc_select_functional

subroutine ensure_functionals_initialized()
    implicit none
    type(xc_f03_func_info_t) :: info_x
    if (initialized) return
    call xc_f03_func_init(func_x_upol, func_x_id, XC_UNPOLARIZED)
    call xc_f03_func_init(func_x_pol,  func_x_id, XC_POLARIZED)
    info_x = xc_f03_func_get_info(func_x_upol)
    use_c = (xc_f03_func_info_get_kind(info_x) .ne. XC_EXCHANGE_CORRELATION)
    if (use_c) then
       call xc_f03_func_init(func_c_upol, func_c_id, XC_UNPOLARIZED)
       call xc_f03_func_init(func_c_pol,  func_c_id, XC_POLARIZED)
    endif
    initialized = .true.
end subroutine ensure_functionals_initialized

logical function xc_uses_tau()
    implicit none
    xc_uses_tau = use_mgga
end function xc_uses_tau

logical function xc_uses_vv10()
    implicit none
    xc_uses_vv10 = use_vv10
end function xc_uses_vv10

subroutine xc_eval(imult, idrv, np, rho_a, rho_b, sigma_aa, sigma_bb, sigma_ab, sigma_t, &
                    ex, ec, d1rho, d1sig, vrhoa, vrhob, vsigmaaa, vsigmabb, vsigmaab, &
                    tau_a, tau_b, vtaua, vtaub)
    implicit none
    integer,intent(in)  :: imult, idrv, np
    real(8),intent(in)  :: rho_a(np), rho_b(np), sigma_aa(np), sigma_bb(np), sigma_ab(np), sigma_t(np)
    real(8),intent(out) :: ex(np), ec(np)
    real(8),intent(out) :: d1rho(np), d1sig(np)
    real(8),intent(out) :: vrhoa(np), vrhob(np), vsigmaaa(np), vsigmabb(np), vsigmaab(np)
    real(8),intent(in),optional  :: tau_a(np), tau_b(np)
    real(8),intent(out),optional :: vtaua(np), vtaub(np)

    real(c_double) :: rho_u(np), sigma_u(np), exc_u(np), vrho_u(np), vsigma_u(np)
    real(c_double) :: rho_p(2*np), sigma_p(3*np)
    real(c_double) :: exc_x(np), exc_c(np)
    real(c_double) :: vrho_x(2*np), vrho_c(2*np), vsigma_x(3*np), vsigma_c(3*np)
    real(c_double) :: lapl_u(np), tau_u(np), vlapl_u(np), vtau_u(np)
    real(c_double) :: lapl_p(2*np), tau_p(2*np), vlapl_p(2*np), vtau_p(2*np)
    real(c_double) :: exc_mx(np), exc_mc(np)
    real(c_double) :: vrho_mx(2*np), vrho_mc(2*np), vsigma_mx(3*np), vsigma_mc(3*np)
    real(c_double) :: vtau_mx(2*np), vtau_mc(2*np)
    integer :: ip

    ex = 0; ec = 0
    d1rho = 0; d1sig = 0
    vrhoa = 0; vrhob = 0; vsigmaaa = 0; vsigmabb = 0; vsigmaab = 0
    if (present(vtaua)) vtaua = 0
    if (present(vtaub)) vtaub = 0

    call ensure_functionals_initialized()

    if (imult .eq. 1) then
        rho_u = rho_a + rho_b
        sigma_u = sigma_t
        if (use_mgga) then
           lapl_u = 0.0d0
           tau_u  = tau_a + tau_b
           call xc_f03_mgga_exc_vxc(func_x_upol, int(np,c_size_t), rho_u, sigma_u, lapl_u, tau_u, &
                                     exc_mx, vrho_u, vsigma_u, vlapl_u, vtau_u)
           ex = exc_mx*rho_u
           d1rho = vrho_u
           d1sig = vsigma_u
           vtaua = vtau_u
           if (use_c) then
              call xc_f03_mgga_exc_vxc(func_c_upol, int(np,c_size_t), rho_u, sigma_u, lapl_u, tau_u, &
                                        exc_mc, vrho_u, vsigma_u, vlapl_u, vtau_u)
              ec = exc_mc*rho_u
              d1rho = d1rho + vrho_u
              d1sig = d1sig + vsigma_u
              vtaua = vtaua + vtau_u
           endif
           vtaub = vtaua
        else if (use_lda) then
           call xc_f03_lda_exc_vxc(func_x_upol, int(np,c_size_t), rho_u, exc_u, vrho_u)
           ex = exc_u*rho_u
           d1rho = vrho_u
           if (use_c) then
              call xc_f03_lda_exc_vxc(func_c_upol, int(np,c_size_t), rho_u, exc_u, vrho_u)
              ec = exc_u*rho_u
              d1rho = d1rho + vrho_u
           endif
        else
           call xc_f03_gga_exc_vxc(func_x_upol, int(np,c_size_t), rho_u, sigma_u, exc_u, vrho_u, vsigma_u)
           ex = exc_u*rho_u
           d1rho = vrho_u
           d1sig = vsigma_u
           if (use_c) then
              call xc_f03_gga_exc_vxc(func_c_upol, int(np,c_size_t), rho_u, sigma_u, exc_u, vrho_u, vsigma_u)
              ec = exc_u*rho_u
              d1rho = d1rho + vrho_u
              d1sig = d1sig + vsigma_u
           endif
        endif
    else
        do ip = 1,np
           rho_p(2*ip-1) = rho_a(ip)
           rho_p(2*ip)   = rho_b(ip)
           sigma_p(3*ip-2) = sigma_aa(ip)
           sigma_p(3*ip-1) = sigma_ab(ip)
           sigma_p(3*ip)   = sigma_bb(ip)
        enddo
        if (use_mgga) then
           do ip = 1,np
              lapl_p(2*ip-1) = 0.0d0
              lapl_p(2*ip)   = 0.0d0
              tau_p(2*ip-1)  = tau_a(ip)
              tau_p(2*ip)    = tau_b(ip)
           enddo
           call xc_f03_mgga_exc_vxc(func_x_pol, int(np,c_size_t), rho_p, sigma_p, lapl_p, tau_p, &
                                     exc_mx, vrho_mx, vsigma_mx, vlapl_p, vtau_mx)
           do ip = 1,np
              ex(ip) = exc_mx(ip)*(rho_a(ip)+rho_b(ip))
              vrhoa(ip) = vrho_mx(2*ip-1)
              vrhob(ip) = vrho_mx(2*ip)
              vsigmaaa(ip) = vsigma_mx(3*ip-2)
              vsigmaab(ip) = vsigma_mx(3*ip-1)
              vsigmabb(ip) = vsigma_mx(3*ip)
              vtaua(ip) = vtau_mx(2*ip-1)
              vtaub(ip) = vtau_mx(2*ip)
           enddo
           if (use_c) then
              call xc_f03_mgga_exc_vxc(func_c_pol, int(np,c_size_t), rho_p, sigma_p, lapl_p, tau_p, &
                                        exc_mc, vrho_mc, vsigma_mc, vlapl_p, vtau_mc)
              do ip = 1,np
                 ec(ip) = exc_mc(ip)*(rho_a(ip)+rho_b(ip))
                 vrhoa(ip) = vrhoa(ip) + vrho_mc(2*ip-1)
                 vrhob(ip) = vrhob(ip) + vrho_mc(2*ip)
                 vsigmaaa(ip) = vsigmaaa(ip) + vsigma_mc(3*ip-2)
                 vsigmaab(ip) = vsigmaab(ip) + vsigma_mc(3*ip-1)
                 vsigmabb(ip) = vsigmabb(ip) + vsigma_mc(3*ip)
                 vtaua(ip) = vtaua(ip) + vtau_mc(2*ip-1)
                 vtaub(ip) = vtaub(ip) + vtau_mc(2*ip)
              enddo
           endif
        else if (use_lda) then
           call xc_f03_lda_exc_vxc(func_x_pol, int(np,c_size_t), rho_p, exc_x, vrho_x)
           do ip = 1,np
              ex(ip) = exc_x(ip)*(rho_a(ip)+rho_b(ip))
              vrhoa(ip) = vrho_x(2*ip-1)
              vrhob(ip) = vrho_x(2*ip)
           enddo
           if (use_c) then
              call xc_f03_lda_exc_vxc(func_c_pol, int(np,c_size_t), rho_p, exc_c, vrho_c)
              do ip = 1,np
                 ec(ip) = exc_c(ip)*(rho_a(ip)+rho_b(ip))
                 vrhoa(ip) = vrhoa(ip) + vrho_c(2*ip-1)
                 vrhob(ip) = vrhob(ip) + vrho_c(2*ip)
              enddo
           endif
        else
           call xc_f03_gga_exc_vxc(func_x_pol, int(np,c_size_t), rho_p, sigma_p, exc_x, vrho_x, vsigma_x)
           do ip = 1,np
              ex(ip) = exc_x(ip)*(rho_a(ip)+rho_b(ip))
              vrhoa(ip) = vrho_x(2*ip-1)
              vrhob(ip) = vrho_x(2*ip)
              vsigmaaa(ip) = vsigma_x(3*ip-2)
              vsigmaab(ip) = vsigma_x(3*ip-1)
              vsigmabb(ip) = vsigma_x(3*ip)
           enddo
           if (use_c) then
              call xc_f03_gga_exc_vxc(func_c_pol, int(np,c_size_t), rho_p, sigma_p, exc_c, vrho_c, vsigma_c)
              do ip = 1,np
                 ec(ip) = exc_c(ip)*(rho_a(ip)+rho_b(ip))
                 vrhoa(ip) = vrhoa(ip) + vrho_c(2*ip-1)
                 vrhob(ip) = vrhob(ip) + vrho_c(2*ip)
                 vsigmaaa(ip) = vsigmaaa(ip) + vsigma_c(3*ip-2)
                 vsigmaab(ip) = vsigmaab(ip) + vsigma_c(3*ip-1)
                 vsigmabb(ip) = vsigmabb(ip) + vsigma_c(3*ip)
              enddo
           endif
        endif
    endif
end subroutine xc_eval

end module mod_xc
