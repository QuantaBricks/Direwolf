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
public :: xc_eval_fxc
public :: xc_select_functional
public :: xc_uses_tau
public :: is_rsh
public :: xc_eval_fxc_mgga
public :: xc_eval_fxc_mgga_pol
public :: xc_uses_vv10
public :: xc_write_references

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
    character(len=len(name)) :: name_canon
    integer :: ic, jc

    name_upper = trim(name)
    do ic = 1, len_trim(name)
       jc = iachar(name_upper(ic:ic))
       if (jc .ge. iachar('a') .and. jc .le. iachar('z')) &
          name_upper(ic:ic) = achar(jc - (iachar('a')-iachar('A')))
    enddo
    name_canon = ''
    jc = 0
    do ic = 1, len_trim(name_upper)
       if (name_upper(ic:ic) == '-' .or. name_upper(ic:ic) == '_' .or. &
           name_upper(ic:ic) == ' ') cycle
       jc = jc + 1
       name_canon(jc:jc) = name_upper(ic:ic)
    enddo

    is_hf = (trim(name_canon) .eq. "HF")
    use_vv10 = .false.

    select case (trim(name_canon))
    case ("B3LYP","B3LYPHYB")
       func_x_id = XC_HYB_GGA_XC_B3LYP
       func_c_id = XC_GGA_C_PBE
    case ("CAMB3LYP")
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
    case ("M062X","M062XHYB")
       func_x_id = XC_HYB_MGGA_X_M06_2X
       func_c_id = XC_MGGA_C_M06_2X
    case ("WB97XD")
       func_x_id = XC_HYB_GGA_XC_WB97X_D
       func_c_id = XC_GGA_C_PBE
    case ("WB97X")
       func_x_id = XC_HYB_GGA_XC_WB97X
       func_c_id = XC_GGA_C_PBE
    case ("PBE0")
       func_x_id = XC_HYB_GGA_XC_PBEH
       func_c_id = XC_GGA_C_PBE
    case ("B973C")
       func_x_id = XC_GGA_XC_B97_3C
       func_c_id = XC_GGA_C_PBE
    case ("PBE","PBEPBE","GGAXPBE,GGACPBE")
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
    case ("M06","M06HYB")
       func_x_id = XC_HYB_MGGA_X_M06
       func_c_id = XC_MGGA_C_M06
    case ("M052X")
       func_x_id = XC_HYB_MGGA_X_M05_2X
       func_c_id = XC_MGGA_C_M05_2X
    case ("M06L")
       func_x_id = XC_MGGA_X_M06_L
       func_c_id = XC_MGGA_C_M06_L
    case ("MN15")
       func_x_id = XC_HYB_MGGA_X_MN15
       func_c_id = XC_MGGA_C_MN15
    case ("MN15L")
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
     case ("RPW86PBE","RPW86","GGAXRPW86,GGACPBE")
        func_x_id = XC_GGA_X_RPW86
        func_c_id = XC_GGA_C_PBE
     case ("WB97MV")
        func_x_id = XC_HYB_MGGA_XC_WB97M_V
        func_c_id = XC_GGA_C_PBE
        use_vv10 = .true.
        call vv10_set(6.0d0, 0.01d0)
    case ("WB97XV")
       func_x_id = XC_HYB_GGA_XC_WB97X_V
       func_c_id = XC_GGA_C_PBE
       use_vv10 = .true.
       call vv10_set(6.0d0, 0.01d0)
    case ("WB97X3C")
       func_x_id = XC_HYB_GGA_XC_WB97X_V
       func_c_id = XC_GGA_C_PBE
    case ("B97MV")
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

subroutine xc_write_references(utxt, ubib)
integer, intent(in) :: utxt, ubib
character(len=128) :: vs, lref, ldoi, fname, lkey
character(len=1024) :: sref, sdoi, sbib
type(xc_f03_func_t) :: ftmp
type(xc_f03_func_info_t) :: finfo
type(xc_f03_func_reference_t) :: fr
integer(c_int) :: ids(2), iref, k, nid, major, minor, micro
integer :: j, nseen
logical :: dup
character(len=1100) :: seen(16)
if (is_hf) return
nseen = 0
call xc_f03_version(major, minor, micro)
call xc_f03_version_string(vs)
call xc_f03_reference(lref)
call xc_f03_reference_doi(ldoi)
write(utxt,'(A)') '  Libxc '//trim(vs)
write(utxt,'(A)') '    '//trim(lref)
write(utxt,'(A)') '    doi:'//trim(ldoi)
if (ubib /= 0) then
   call xc_f03_reference_key(lkey)
   write(ubib,'(A)') '@article{'//trim(lkey)//','
   write(ubib,'(A)') '  author = {Lehtola, Susi and Steigemann, Conrad and Oliveira, Micael J. T. and Marques, Miguel A. L.},'
   write(ubib,'(A)') '  title = {{Recent developments in libxc - A comprehensive library of functionals for density '// &
        'functional theory}},'
   write(ubib,'(A)') '  journaltitle = {SoftwareX},'
   write(ubib,'(A)') '  date = {2018},'
   write(ubib,'(A)') '  volume = {7},'
   write(ubib,'(A)') '  pages = {1--5},'
   write(ubib,'(A)') '  doi = {'//trim(ldoi)//'},'
   write(ubib,'(A)') '  note = {Libxc '//trim(vs)//'},'
   write(ubib,'(A)') '}'
   write(ubib,'(A)') ''
endif
nid = 1 ; ids(1) = func_x_id
if (use_c) then
   nid = 2 ; ids(2) = func_c_id
endif
do k = 1, nid
   call xc_f03_func_init(ftmp, ids(k), XC_UNPOLARIZED)
   finfo = xc_f03_func_get_info(ftmp)
   fname = xc_f03_func_info_get_name(finfo)
   write(utxt,'(A)') '  Functional '//trim(fname)
   iref = 0
   do while (iref >= 0)
      fr = xc_f03_func_info_get_references(finfo, iref)
      sref = xc_f03_func_reference_get_ref(fr)
      sdoi = xc_f03_func_reference_get_doi(fr)
      dup = .false.
      do j = 1, nseen
         if (trim(seen(j)) == trim(sref)//'|'//trim(sdoi)) dup = .true.
      enddo
      if (.not. dup) then
         if (nseen < size(seen)) then
            nseen = nseen + 1
            seen(nseen) = trim(sref)//'|'//trim(sdoi)
         endif
         write(utxt,'(A)') '    '//trim(sref)
         if (len_trim(sdoi) > 0) write(utxt,'(A)') '    doi:'//trim(sdoi)
         if (ubib /= 0) then
            sbib = xc_f03_func_reference_get_bibtex(fr)
            if (len_trim(sbib) > 0) write(ubib,'(A)') trim(sbib)
         endif
      endif
   enddo
   call xc_f03_func_end(ftmp)
enddo
end subroutine xc_write_references

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

    if (use_mgga .and. .not. (present(tau_a) .and. present(vtaua))) &
       error stop 'xc_eval: meta-GGA functional called without tau_a/vtaua'

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

subroutine xc_eval_fxc(imult, np, rho_a, rho_b, sigma_aa, sigma_bb, sigma_ab, sigma_t, &
                        v2rho2, v2rhosigma, v2sigma2)
    implicit none
    integer,intent(in)  :: imult, np
    real(8),intent(in)  :: rho_a(np), rho_b(np), sigma_aa(np), sigma_bb(np), sigma_ab(np), sigma_t(np)
    real(8),intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)

    real(c_double) :: rho_u(np), sigma_u(np)
    real(c_double) :: rho_p(2*np), sigma_p(3*np)
    real(c_double) :: t2rho(np), t2rhosig(3*np), t2sig(6*np)
    real(c_double) :: t2rho_p(3*np), t2rhosig_p(6*np), t2sig_p(6*np)
    real(c_double) :: c2rho_p(3*np), c2rhosig_p(6*np), c2sig_p(6*np)
    integer :: ip, k

    call ensure_functionals_initialized()
    if (use_mgga) then
       print *, 'ERROR: xc_eval_fxc: meta-GGA fxc not implemented yet.'
       stop 1
    endif
    do k = 1,merge(6,3,imult.ne.1)*np
       v2rhosigma(k) = 0.0d0
    enddo
    do k = 1,6*np
       v2sigma2(k) = 0.0d0
    enddo

    if (imult .eq. 1) then
        rho_u = rho_a + rho_b
        sigma_u = sigma_t
        if (use_lda) then
           call xc_f03_lda_fxc(func_x_upol, int(np,c_size_t), rho_u, t2rho)
           if (use_c) then
              call xc_f03_lda_fxc(func_c_upol, int(np,c_size_t), rho_u, t2rho_p(1:np))
              t2rho = t2rho + t2rho_p(1:np)
           endif
           do k = 1,np
              v2rho2(k) = t2rho(k)
           enddo
        else
           call xc_f03_gga_fxc(func_x_upol, int(np,c_size_t), rho_u, sigma_u, &
                                t2rho, t2rhosig, t2sig)
           if (use_c) then
              call xc_f03_gga_fxc(func_c_upol, int(np,c_size_t), rho_u, sigma_u, &
                                   t2rho_p(1:np), t2rhosig_p(1:3*np), t2sig_p(1:6*np))
              t2rho = t2rho + t2rho_p(1:np)
              do k = 1,3*np
                 t2rhosig(k) = t2rhosig(k) + t2rhosig_p(k)
              enddo
              do k = 1,6*np
                 t2sig(k) = t2sig(k) + t2sig_p(k)
              enddo
           endif
           do k = 1,np
              v2rho2(k) = t2rho(k)
           enddo
           do k = 1,3*np
              v2rhosigma(k) = t2rhosig(k)
           enddo
           do k = 1,6*np
              v2sigma2(k) = t2sig(k)
           enddo
        endif
    else
        do ip = 1,np
           rho_p(2*ip-1) = rho_a(ip)
           rho_p(2*ip)   = rho_b(ip)
           sigma_p(3*ip-2) = sigma_aa(ip)
           sigma_p(3*ip-1) = sigma_ab(ip)
           sigma_p(3*ip)   = sigma_bb(ip)
        enddo
        if (use_lda) then
           call xc_f03_lda_fxc(func_x_pol, int(np,c_size_t), rho_p, t2rho_p)
           if (use_c) then
              call xc_f03_lda_fxc(func_c_pol, int(np,c_size_t), rho_p, c2rho_p)
              t2rho_p = t2rho_p + c2rho_p
           endif
           do k = 1,3*np
              v2rho2(k) = t2rho_p(k)
           enddo
        else
           call xc_f03_gga_fxc(func_x_pol, int(np,c_size_t), rho_p, sigma_p, &
                                t2rho_p, t2rhosig_p, t2sig_p)
           if (use_c) then
              call xc_f03_gga_fxc(func_c_pol, int(np,c_size_t), rho_p, sigma_p, &
                                   c2rho_p, c2rhosig_p, c2sig_p)
              t2rho_p = t2rho_p + c2rho_p
              t2rhosig_p = t2rhosig_p + c2rhosig_p
              t2sig_p = t2sig_p + c2sig_p
           endif
           do k = 1,3*np
              v2rho2(k) = t2rho_p(k)
           enddo
           do k = 1,6*np
              v2rhosigma(k) = t2rhosig_p(k)
           enddo
           do k = 1,6*np
              v2sigma2(k) = t2sig_p(k)
           enddo
        endif
    endif
end subroutine xc_eval_fxc

subroutine xc_eval_fxc_mgga(np, rho, sigma, tau, v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2sigmatau, v2tau2)
    implicit none
    integer,intent(in)  :: np
    real(8),intent(in)  :: rho(np), sigma(np), tau(np)
    real(8),intent(out) :: v2rho2(np), v2rhosigma(np), v2sigma2(np), v2rhotau(np), v2sigmatau(np), v2tau2(np)
    real(c_double) :: lapl(np), z1(np), z2(np), z3(np), z4(np), z5(np)
    real(c_double) :: a1(np), a2(np), a3(np), a4(np), a5(np), a6(np)
    lapl = 0.0d0
    call ensure_functionals_initialized()
    if (.not. use_mgga) then
       print *, 'ERROR: xc_eval_fxc_mgga: not a meta-GGA.'
       stop 1
    endif
    call xc_f03_mgga_fxc(func_x_upol, int(np,c_size_t), rho, sigma, lapl, tau, &
         a1, a2, z1, a3, a4, z2, a5, z3, z4, a6)
    v2rho2 = a1 ; v2rhosigma = a2 ; v2rhotau = a3
    v2sigma2 = a4 ; v2sigmatau = a5 ; v2tau2 = a6
    if (use_c) then
       call xc_f03_mgga_fxc(func_c_upol, int(np,c_size_t), rho, sigma, lapl, tau, &
            a1, a2, z1, a3, a4, z2, a5, z3, z4, a6)
       v2rho2 = v2rho2 + a1 ; v2rhosigma = v2rhosigma + a2 ; v2rhotau = v2rhotau + a3
       v2sigma2 = v2sigma2 + a4 ; v2sigmatau = v2sigmatau + a5 ; v2tau2 = v2tau2 + a6
    endif
end subroutine xc_eval_fxc_mgga

subroutine xc_eval_fxc_mgga_pol(np, ra, rb, saa, sab, sbb, ta, tb, K)
    implicit none
    integer,intent(in)  :: np
    real(8),intent(in)  :: ra(np), rb(np), saa(np), sab(np), sbb(np), ta(np), tb(np)
    real(8),intent(out) :: K(7,7,np)
    real(c_double) :: rho_p(2*np), sig_p(3*np), tau_p(2*np), lapl_p(2*np)
    real(c_double) :: f_rr(3*np), f_rs(6*np), f_rl(6*np), f_rt(4*np), f_ss(6*np), f_sl(9*np)
    real(c_double) :: f_st(6*np), f_ll(3*np), f_lt(4*np), f_tt(3*np)
    integer :: ip, pass
    call ensure_functionals_initialized()
    if (.not. use_mgga) then
       print *, 'ERROR: xc_eval_fxc_mgga_pol: not a meta-GGA.'
       stop 1
    endif
    do ip = 1,np
       rho_p(2*ip-1) = ra(ip) ; rho_p(2*ip) = rb(ip)
       sig_p(3*ip-2) = saa(ip) ; sig_p(3*ip-1) = sab(ip) ; sig_p(3*ip) = sbb(ip)
       tau_p(2*ip-1) = ta(ip) ; tau_p(2*ip) = tb(ip)
       lapl_p(2*ip-1) = 0.0d0 ; lapl_p(2*ip) = 0.0d0
    enddo
    K = 0.0d0
    do pass = 1,2
       if (pass .eq. 2 .and. .not. use_c) exit
       if (pass .eq. 1) then
          call xc_f03_mgga_fxc(func_x_pol, int(np,c_size_t), rho_p, sig_p, lapl_p, tau_p, &
               f_rr, f_rs, f_rl, f_rt, f_ss, f_sl, f_st, f_ll, f_lt, f_tt)
       else
          call xc_f03_mgga_fxc(func_c_pol, int(np,c_size_t), rho_p, sig_p, lapl_p, tau_p, &
               f_rr, f_rs, f_rl, f_rt, f_ss, f_sl, f_st, f_ll, f_lt, f_tt)
       endif
       do ip = 1,np
          K(1,1,ip) = K(1,1,ip) + f_rr(3*ip-2) ; K(1,2,ip) = K(1,2,ip) + f_rr(3*ip-1) ; K(2,2,ip) = K(2,2,ip) + f_rr(3*ip)
          K(1,3:5,ip) = K(1,3:5,ip) + f_rs(6*ip-5:6*ip-3)
          K(2,3:5,ip) = K(2,3:5,ip) + f_rs(6*ip-2:6*ip)
          K(1,6:7,ip) = K(1,6:7,ip) + f_rt(4*ip-3:4*ip-2)
          K(2,6:7,ip) = K(2,6:7,ip) + f_rt(4*ip-1:4*ip)
          K(3,3,ip) = K(3,3,ip) + f_ss(6*ip-5) ; K(3,4,ip) = K(3,4,ip) + f_ss(6*ip-4) ; K(3,5,ip) = K(3,5,ip) + f_ss(6*ip-3)
          K(4,4,ip) = K(4,4,ip) + f_ss(6*ip-2) ; K(4,5,ip) = K(4,5,ip) + f_ss(6*ip-1) ; K(5,5,ip) = K(5,5,ip) + f_ss(6*ip)
          K(3,6:7,ip) = K(3,6:7,ip) + f_st(6*ip-5:6*ip-4)
          K(4,6:7,ip) = K(4,6:7,ip) + f_st(6*ip-3:6*ip-2)
          K(5,6:7,ip) = K(5,6:7,ip) + f_st(6*ip-1:6*ip)
          K(6,6,ip) = K(6,6,ip) + f_tt(3*ip-2) ; K(6,7,ip) = K(6,7,ip) + f_tt(3*ip-1) ; K(7,7,ip) = K(7,7,ip) + f_tt(3*ip)
       enddo
    enddo
    do ip = 1,np
       do pass = 1,7
          K(pass+1:7,pass,ip) = K(pass,pass+1:7,ip)
       enddo
    enddo
end subroutine xc_eval_fxc_mgga_pol

end module mod_xc
