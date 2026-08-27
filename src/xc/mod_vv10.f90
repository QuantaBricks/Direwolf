! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_vv10: the Vydrov-Van Voorhis 2010 nonlocal correlation functional

module mod_vv10
implicit none
private
public :: vv10_set
public :: vv10_evaluate
public :: vv10_report
public :: vv10_nlc_grid_reset
public :: vv10_nlc_grid_build
public :: vv10_nlc_grid_retrieve
public :: vv10_nlc_basis_build
public :: NLC_CHUNK
public :: vv10_nlc_chunk_info
public :: vv10_nlc_chunk_val
public :: vv10_omega_kappa
public :: vv10_grid_force
public :: vv10_beta_val
public :: vv10_rcut_val
public :: engine_vv10_nonself
public :: vv10_dynamic
public :: vv10_active_now
public :: vv10_just_activated
public :: vv10_switch_prms
public :: vv10_flush_on_switch
public :: VV10_DYNGRID_NCONTS
public :: vv10_report_set
public :: vv10_nlc_radial
public :: vv10_nlc_angular

interface
   subroutine gridgen_nlc(nrad, nsph, npts, coor_out, weight_out, coor_in, per_atom_period_scale, intacc_eps)
   integer,intent(in) :: nrad, nsph
   integer,intent(out) :: npts
   real(8),allocatable,intent(out) :: coor_out(:,:), weight_out(:)
   real(8),optional,intent(in) :: coor_in(3,*)
   logical,optional,intent(in) :: per_atom_period_scale
   real(8),optional,intent(in) :: intacc_eps
   end subroutine gridgen_nlc
end interface

logical,save :: vv10_initialized = .false.
logical,save :: engine_vv10_nonself = .false.
logical,save :: vv10_dynamic = .false.
logical,save :: vv10_active_now = .true.
logical,save :: vv10_just_activated = .false.
real(8),save :: vv10_switch_prms = 3.0d-4
logical,save :: vv10_flush_on_switch = .true.
integer,parameter :: VV10_DYNGRID_NCONTS = 350
real(8),save :: vv10_b = 5.9d0
real(8),save :: vv10_C = 0.0093d0
real(8),save :: vv10_beta = 0.0d0
real(8),save :: vv10_rcut = 17.0d0
real(8),save :: vv10_rhocut = 1.0d-9
logical,save :: vv10_computed = .false.
real(8),save :: vv10_enl_last = 0.0d0
logical,save :: nlc_grid_built = .false.
integer,save :: nlc_npts = 0
real(8),allocatable,save :: nlc_coor(:,:), nlc_weight(:)
integer,save :: nlc_nrad = 20
integer,save :: nlc_nsph = 110
integer,parameter :: NLC_CHUNK = 128
type :: NlcValBlock
   real(4),allocatable :: val0(:,:)
   real(4),allocatable :: val1(:,:,:)
end type NlcValBlock
logical,save :: nlc_basis_built = .false.
type(NlcValBlock),allocatable,save,target :: nlc_val_blocks(:)
integer,allocatable,save :: nlc_chunk_nsig(:)
integer,allocatable,save :: nlc_chunk_sig_idx(:,:)
integer,save :: nlc_n_chunks = 0

contains

subroutine vv10_set(b_in, C_in)
    implicit none
    real(8),intent(in) :: b_in, C_in
    character(len=64) :: buf
    vv10_b = b_in
    vv10_C = C_in
    vv10_beta = (1.0d0/32.0d0) * (3.0d0/(vv10_b*vv10_b))**0.75d0
    buf = ""
    call get_environment_variable("ENGINE_VV10_RCUT", buf)
    if (len_trim(buf) .gt. 0) read(buf,*) vv10_rcut
    buf = ""
    call get_environment_variable("ENGINE_VV10_RHOCUT", buf)
    if (len_trim(buf) .gt. 0) read(buf,*) vv10_rhocut
    vv10_initialized = .true.
end subroutine vv10_set

subroutine vv10_evaluate(np, rho, gamma, weight, coor, Enl, eps, Fn, Fg)
    implicit none
    integer,intent(in) :: np
    real(8),intent(in)  :: rho(np), gamma(np), weight(np), coor(3,np)
    real(8),intent(out) :: Enl
    real(8),intent(out) :: eps(np), Fn(np), Fg(np)
    real(8),allocatable :: omega0(:), kappa(:), domega_dn(:), domega_dg(:), dkappa_dn(:)
    real(8),allocatable :: x(:), y(:), z(:)
    integer,allocatable :: act(:)
    real(8) :: pi
    real(4) :: dx4, dy4, dz4, R24, gi4, gj4, phi4, invgi4, invgj4, invgij4, omi4, ki4
    real(8) :: A, U, W, n_i
    integer :: i, nj, ni, n_active, ia, ja
    real(4),allocatable :: xs(:), ys(:), zs(:), wr(:), os(:), ks(:)
    real(8),allocatable :: rs(:)
    integer,parameter :: cellrad = 4
    integer :: nx, ny, nz, ncell, ci, b_s, b_e, ddx, ddy, ddz, ccx, ccy, ccz, nc
    integer,allocatable :: cellidx(:), cellcnt(:), cellstart(:), fillpos(:), actsorted(:)
    integer,allocatable :: pcx(:), pcy(:), pcz(:)
    real(8) :: xmin, xmax, ymin, ymax, zmin, zmax, cellsz
    logical :: noscreen
    if (.not. vv10_initialized) call vv10_set(vv10_b, vv10_C)
    pi = 4.0d0*atan(1.0d0)

    allocate(omega0(np), kappa(np), domega_dn(np), domega_dg(np), dkappa_dn(np))
    allocate(x(np), y(np), z(np))
    allocate(act(np))
    do i = 1,np
       x(i) = coor(1,i)
       y(i) = coor(2,i)
       z(i) = coor(3,i)
    enddo

    eps = vv10_beta
    Fn  = vv10_beta
    Fg  = 0.0d0
    n_active = 0
    do i = 1,np
       if (rho(i) .le. vv10_rhocut) cycle
       n_active = n_active + 1
       act(n_active) = i
       omega0(i) = sqrt(vv10_C*gamma(i)*gamma(i)/rho(i)**4 + 4.0d0*pi*rho(i)/3.0d0)
       kappa(i)  = vv10_b*(3.0d0*pi/2.0d0)*(rho(i)/(9.0d0*pi))**(1.0d0/6.0d0)
       domega_dn(i) = (4.0d0*pi/3.0d0 - 4.0d0*vv10_C*gamma(i)*gamma(i)/rho(i)**5)/(2.0d0*omega0(i))
       domega_dg(i) = vv10_C*gamma(i)/(omega0(i)*rho(i)**4)
       dkappa_dn(i) = kappa(i)/(6.0d0*rho(i))
    enddo

    noscreen = .true.
    if (n_active .gt. 0) then
       xmin = 1.0d30; xmax = -1.0d30
       ymin = 1.0d30; ymax = -1.0d30
       zmin = 1.0d30; zmax = -1.0d30
       do i = 1,n_active
          ia = act(i)
          if (x(ia) .lt. xmin) xmin = x(ia)
          if (x(ia) .gt. xmax) xmax = x(ia)
          if (y(ia) .lt. ymin) ymin = y(ia)
          if (y(ia) .gt. ymax) ymax = y(ia)
          if (z(ia) .lt. zmin) zmin = z(ia)
          if (z(ia) .gt. zmax) zmax = z(ia)
       enddo
       cellsz = vv10_rcut/dble(cellrad)
       nx = int((xmax-xmin)/cellsz) + 1
       ny = int((ymax-ymin)/cellsz) + 1
       nz = int((zmax-zmin)/cellsz) + 1
       noscreen = (nx .le. 2*cellrad+1) .and. (ny .le. 2*cellrad+1) .and. (nz .le. 2*cellrad+1)
    endif
    if (n_active .gt. 0 .and. .not. noscreen) then
       ncell = nx*ny*nz
       allocate(cellidx(n_active))
       do i = 1,n_active
          ia = act(i)
          cellidx(i) = min(nx-1,int((x(ia)-xmin)/cellsz)) &
                     + nx*(min(ny-1,int((y(ia)-ymin)/cellsz)) &
                     + ny*min(nz-1,int((z(ia)-zmin)/cellsz)))
       enddo
       allocate(cellcnt(0:ncell-1), cellstart(0:ncell))
       cellcnt = 0
       do i = 1,n_active
          cellcnt(cellidx(i)) = cellcnt(cellidx(i)) + 1
       enddo
       cellstart(0) = 1
       do ci = 0,ncell-1
          cellstart(ci+1) = cellstart(ci) + cellcnt(ci)
       enddo
       allocate(fillpos(0:ncell-1), actsorted(n_active))
       fillpos(0:ncell-1) = cellstart(0:ncell-1)
       do i = 1,n_active
          ci = cellidx(i)
          actsorted(fillpos(ci)) = act(i)
          fillpos(ci) = fillpos(ci) + 1
       enddo
       act(1:n_active) = actsorted(1:n_active)
       deallocate(actsorted, fillpos, cellidx, cellcnt)
    endif

    allocate(xs(n_active), ys(n_active), zs(n_active), wr(n_active), &
             rs(n_active), os(n_active), ks(n_active))
    do i = 1,n_active
       ia = act(i)
       xs(i) = real(x(ia),4)
       ys(i) = real(y(ia),4)
       zs(i) = real(z(ia),4)
       rs(i) = rho(ia)
       wr(i) = real(weight(ia)*rho(ia),4)
       os(i) = real(omega0(ia),4)
       ks(i) = real(kappa(ia),4)
    enddo

    if (.not. noscreen) then
       allocate(pcx(n_active), pcy(n_active), pcz(n_active))
       do i = 1,n_active
          ia = act(i)
          pcx(i) = min(nx-1, int((x(ia)-xmin)/cellsz))
          pcy(i) = min(ny-1, int((y(ia)-ymin)/cellsz))
          pcz(i) = min(nz-1, int((z(ia)-zmin)/cellsz))
       enddo
    endif

    if (noscreen) then
       !$omp parallel do schedule(dynamic,64) &
       !$omp& private(nj,dx4,dy4,dz4,R24,gi4,gj4,phi4,invgi4,invgj4,invgij4,A,U,W,n_i,omi4,ki4)
       do ni = 1,n_active
          n_i = rs(ni)
          omi4 = os(ni)
          ki4 = ks(ni)
          A = 0.0d0
          U = 0.0d0
          W = 0.0d0
          do nj = 1,n_active
             dx4 = xs(ni) - xs(nj)
             dy4 = ys(ni) - ys(nj)
             dz4 = zs(ni) - zs(nj)
             R24 = dx4*dx4 + dy4*dy4 + dz4*dz4
             gi4 = omi4*R24 + ki4
             gj4 = os(nj)*R24 + ks(nj)
             invgi4  = 1.0_4/gi4
             invgj4  = 1.0_4/gj4
             invgij4 = 1.0_4/(gi4+gj4)
             phi4 = -1.5_4*invgi4*invgj4*invgij4
             A = A + dble(wr(nj))*dble(phi4)
             U = U - dble(wr(nj))*dble(phi4)*dble(invgi4 + invgij4)
             W = W - dble(wr(nj))*dble(phi4)*dble(R24)*dble(invgi4 + invgij4)
          enddo
          eps(act(ni)) = vv10_beta + 0.5d0*A
          Fn(act(ni))  = vv10_beta + A + n_i*(dkappa_dn(act(ni))*U + domega_dn(act(ni))*W)
          Fg(act(ni))  = n_i*domega_dg(act(ni))*W
       enddo
       !$omp end parallel do
    else
       !$omp parallel do schedule(dynamic,64) &
       !$omp& private(nj,dx4,dy4,dz4,R24,gi4,gj4,phi4,invgi4,invgj4,invgij4,A,U,W,n_i,omi4,ki4, &
       !$omp& ddx,ddy,ddz,ccx,ccy,ccz,nc,b_s,b_e)
       do ni = 1,n_active
          n_i = rs(ni)
          omi4 = os(ni)
          ki4 = ks(ni)
          A = 0.0d0
          U = 0.0d0
          W = 0.0d0
          do ddz = -cellrad,cellrad
             ccz = pcz(ni)+ddz
             if (ccz .lt. 0 .or. ccz .ge. nz) cycle
             do ddy = -cellrad,cellrad
                ccy = pcy(ni)+ddy
                if (ccy .lt. 0 .or. ccy .ge. ny) cycle
                do ddx = -cellrad,cellrad
                   ccx = pcx(ni)+ddx
                   if (ccx .lt. 0 .or. ccx .ge. nx) cycle
                   nc = ccx + nx*(ccy + ny*ccz)
                   b_s = cellstart(nc)
                   b_e = cellstart(nc+1)-1
                   do nj = b_s,b_e
                      dx4 = xs(ni) - xs(nj)
                      dy4 = ys(ni) - ys(nj)
                      dz4 = zs(ni) - zs(nj)
                      R24 = dx4*dx4 + dy4*dy4 + dz4*dz4
                      gi4 = omi4*R24 + ki4
                      gj4 = os(nj)*R24 + ks(nj)
                      invgi4  = 1.0_4/gi4
                      invgj4  = 1.0_4/gj4
                      invgij4 = 1.0_4/(gi4+gj4)
                      phi4 = -1.5_4*invgi4*invgj4*invgij4
                      A = A + dble(wr(nj))*dble(phi4)
                      U = U - dble(wr(nj))*dble(phi4)*dble(invgi4 + invgij4)
                      W = W - dble(wr(nj))*dble(phi4)*dble(R24)*dble(invgi4 + invgij4)
                   enddo
                enddo
             enddo
          enddo
          eps(act(ni)) = vv10_beta + 0.5d0*A
          Fn(act(ni))  = vv10_beta + A + n_i*(dkappa_dn(act(ni))*U + domega_dn(act(ni))*W)
          Fg(act(ni))  = n_i*domega_dg(act(ni))*W
       enddo
       !$omp end parallel do
    endif
    deallocate(xs, ys, zs, wr, rs, os, ks)
    if (.not. noscreen) deallocate(pcx, pcy, pcz, cellstart)

    Enl = 0.0d0
    do i = 1,np
       Enl = Enl + weight(i)*rho(i)*eps(i)
    enddo
    vv10_enl_last = Enl
    vv10_computed = .true.

    deallocate(omega0, kappa, domega_dn, domega_dg, dkappa_dn)
    deallocate(x, y, z)
    deallocate(act)
end subroutine vv10_evaluate

integer function vv10_nlc_radial()
    implicit none
    vv10_nlc_radial = nlc_nrad
end function vv10_nlc_radial

integer function vv10_nlc_angular()
    implicit none
    vv10_nlc_angular = nlc_nsph
end function vv10_nlc_angular

real(8) function vv10_beta_val()
    implicit none
    vv10_beta_val = vv10_beta
end function vv10_beta_val

real(8) function vv10_rcut_val()
    implicit none
    vv10_rcut_val = vv10_rcut
end function vv10_rcut_val

subroutine vv10_omega_kappa(np, rho, sig, omega0, kappa, n_active, act)
    implicit none
    integer,intent(in) :: np
    real(8),intent(in)  :: rho(np), sig(np)
    real(8),intent(out) :: omega0(np), kappa(np)
    integer,intent(out) :: n_active
    integer,intent(out) :: act(np)
    real(8) :: pi
    integer :: i
    if (.not. vv10_initialized) call vv10_set(vv10_b, vv10_C)
    pi = 4.0d0*atan(1.0d0)
    n_active = 0
    do i = 1,np
       if (rho(i) .le. vv10_rhocut) cycle
       n_active = n_active + 1
       act(n_active) = i
       omega0(i) = sqrt(vv10_C*sig(i)*sig(i)/rho(i)**4 + 4.0d0*pi*rho(i)/3.0d0)
       kappa(i)  = vv10_b*(3.0d0*pi/2.0d0)*(rho(i)/(9.0d0*pi))**(1.0d0/6.0d0)
    enddo
end subroutine vv10_omega_kappa

subroutine vv10_grid_force(np, rho, sig, weight, coor, nper, natoms, gv)
    implicit none
    integer,intent(in) :: np, nper, natoms
    real(8),intent(in)  :: rho(np), sig(np), weight(np), coor(3,np)
    real(8),intent(inout) :: gv(natoms,3)
    real(8),allocatable :: omega0(:), kappa(:), x(:), y(:), z(:)
    real(8),allocatable :: gv_local(:,:)
    integer,allocatable :: act(:), atid(:)
    integer :: n_active, ia, ja, i, ati, atj
    real(8) :: rcut2, pref, gx, gy, gz
    real(4) :: dx4, dy4, dz4, r24, gi4, gj4, phi4, invgi4, invgj4, invgij4, Q4
    real(4),allocatable :: xs(:), ys(:), zs(:), wr(:), os(:), ks(:)
    integer,parameter :: cellrad = 4
    integer :: nx, ny, nz, ncell, ci, b_s, b_e, ddx, ddy, ddz, ccx, ccy, ccz, nc
    integer,allocatable :: cellidx(:), cellcnt(:), cellstart(:), fillpos(:), actsorted(:)
    integer,allocatable :: pcx(:), pcy(:), pcz(:)
    real(8) :: xmin, xmax, ymin, ymax, zmin, zmax, cellsz
    logical :: noscreen

    if (.not. vv10_initialized) call vv10_set(vv10_b, vv10_C)
    allocate(omega0(np), kappa(np), act(np))
    call vv10_omega_kappa(np, rho, sig, omega0, kappa, n_active, act)
    rcut2 = vv10_rcut*vv10_rcut
    if (n_active .eq. 0) then
       deallocate(omega0, kappa, act)
       return
    endif

    allocate(x(np), y(np), z(np))
    do i = 1,np
       x(i) = coor(1,i)
       y(i) = coor(2,i)
       z(i) = coor(3,i)
    enddo

    xmin = 1.0d30; xmax = -1.0d30
    ymin = 1.0d30; ymax = -1.0d30
    zmin = 1.0d30; zmax = -1.0d30
    do i = 1,n_active
       ia = act(i)
       if (x(ia) .lt. xmin) xmin = x(ia)
       if (x(ia) .gt. xmax) xmax = x(ia)
       if (y(ia) .lt. ymin) ymin = y(ia)
       if (y(ia) .gt. ymax) ymax = y(ia)
       if (z(ia) .lt. zmin) zmin = z(ia)
       if (z(ia) .gt. zmax) zmax = z(ia)
    enddo
    cellsz = vv10_rcut/dble(cellrad)
    nx = int((xmax-xmin)/cellsz) + 1
    ny = int((ymax-ymin)/cellsz) + 1
    nz = int((zmax-zmin)/cellsz) + 1
    noscreen = (nx .le. 2*cellrad+1) .and. (ny .le. 2*cellrad+1) .and. (nz .le. 2*cellrad+1)

    if (.not. noscreen) then
       ncell = nx*ny*nz
       allocate(cellidx(n_active))
       do i = 1,n_active
          ia = act(i)
          cellidx(i) = min(nx-1,int((x(ia)-xmin)/cellsz)) &
                     + nx*(min(ny-1,int((y(ia)-ymin)/cellsz)) &
                     + ny*min(nz-1,int((z(ia)-zmin)/cellsz)))
       enddo
       allocate(cellcnt(0:ncell-1), cellstart(0:ncell))
       cellcnt = 0
       do i = 1,n_active
          cellcnt(cellidx(i)) = cellcnt(cellidx(i)) + 1
       enddo
       cellstart(0) = 1
       do ci = 0,ncell-1
          cellstart(ci+1) = cellstart(ci) + cellcnt(ci)
       enddo
       allocate(fillpos(0:ncell-1), actsorted(n_active))
       fillpos(0:ncell-1) = cellstart(0:ncell-1)
       do i = 1,n_active
          ci = cellidx(i)
          actsorted(fillpos(ci)) = act(i)
          fillpos(ci) = fillpos(ci) + 1
       enddo
       act(1:n_active) = actsorted(1:n_active)
       deallocate(actsorted, fillpos, cellidx, cellcnt)
    endif

    allocate(xs(n_active), ys(n_active), zs(n_active), wr(n_active), &
             os(n_active), ks(n_active), atid(n_active))
    do i = 1,n_active
       ia = act(i)
       xs(i) = real(x(ia),4)
       ys(i) = real(y(ia),4)
       zs(i) = real(z(ia),4)
       wr(i) = real(weight(ia)*rho(ia),4)
       os(i) = real(omega0(ia),4)
       ks(i) = real(kappa(ia),4)
       atid(i) = (ia-1)/nper + 1
    enddo

    if (.not. noscreen) then
       allocate(pcx(n_active), pcy(n_active), pcz(n_active))
       do i = 1,n_active
          ia = act(i)
          pcx(i) = min(nx-1, int((x(ia)-xmin)/cellsz))
          pcy(i) = min(ny-1, int((y(ia)-ymin)/cellsz))
          pcz(i) = min(nz-1, int((z(ia)-zmin)/cellsz))
       enddo
    endif

    allocate(gv_local(natoms,3))
    gv_local = 0.0d0

    if (noscreen) then
       !$omp parallel do schedule(dynamic,64) &
       !$omp& private(ja,dx4,dy4,dz4,r24,gi4,gj4,phi4,invgi4,invgj4,invgij4,Q4,pref,ati,atj) &
       !$omp& reduction(+:gv_local) private(gx,gy,gz)
       do ia = 1,n_active
          ati = atid(ia)
          gx = 0.0d0; gy = 0.0d0; gz = 0.0d0
          do ja = 1,n_active
             atj = atid(ja)
             if (atj .eq. ati) cycle
             dx4 = xs(ia) - xs(ja)
             dy4 = ys(ia) - ys(ja)
             dz4 = zs(ia) - zs(ja)
             r24 = dx4*dx4 + dy4*dy4 + dz4*dz4
             if (dble(r24) .gt. rcut2) cycle
             gi4 = os(ia)*r24 + ks(ia)
             gj4 = os(ja)*r24 + ks(ja)
             invgi4  = 1.0_4/gi4
             invgj4  = 1.0_4/gj4
             invgij4 = 1.0_4/(gi4+gj4)
             phi4 = -1.5_4*invgi4*invgj4*invgij4
             Q4 = -2.0_4*phi4*(os(ia)*invgi4 + os(ja)*invgj4 + (os(ia)+os(ja))*invgij4)
             pref = dble(wr(ia))*dble(wr(ja))*dble(Q4)
             gx = gx + pref*dble(dx4)
             gy = gy + pref*dble(dy4)
             gz = gz + pref*dble(dz4)
          enddo
          gv_local(ati,1) = gv_local(ati,1) + gx
          gv_local(ati,2) = gv_local(ati,2) + gy
          gv_local(ati,3) = gv_local(ati,3) + gz
       enddo
       !$omp end parallel do
    else
       !$omp parallel do schedule(dynamic,64) &
       !$omp& private(ja,dx4,dy4,dz4,r24,gi4,gj4,phi4,invgi4,invgj4,invgij4,Q4,pref,ati,atj, &
       !$omp& ddx,ddy,ddz,ccx,ccy,ccz,nc,b_s,b_e) &
       !$omp& reduction(+:gv_local) private(gx,gy,gz)
       do ia = 1,n_active
          ati = atid(ia)
          gx = 0.0d0; gy = 0.0d0; gz = 0.0d0
          do ddz = -cellrad,cellrad
             ccz = pcz(ia)+ddz
             if (ccz .lt. 0 .or. ccz .ge. nz) cycle
             do ddy = -cellrad,cellrad
                ccy = pcy(ia)+ddy
                if (ccy .lt. 0 .or. ccy .ge. ny) cycle
                do ddx = -cellrad,cellrad
                   ccx = pcx(ia)+ddx
                   if (ccx .lt. 0 .or. ccx .ge. nx) cycle
                   nc = ccx + nx*(ccy + ny*ccz)
                   b_s = cellstart(nc)
                   b_e = cellstart(nc+1)-1
                   do ja = b_s,b_e
                      atj = atid(ja)
                      if (atj .eq. ati) cycle
                      dx4 = xs(ia) - xs(ja)
                      dy4 = ys(ia) - ys(ja)
                      dz4 = zs(ia) - zs(ja)
                      r24 = dx4*dx4 + dy4*dy4 + dz4*dz4
                      if (dble(r24) .gt. rcut2) cycle
                      gi4 = os(ia)*r24 + ks(ia)
                      gj4 = os(ja)*r24 + ks(ja)
                      invgi4  = 1.0_4/gi4
                      invgj4  = 1.0_4/gj4
                      invgij4 = 1.0_4/(gi4+gj4)
                      phi4 = -1.5_4*invgi4*invgj4*invgij4
                      Q4 = -2.0_4*phi4*(os(ia)*invgi4 + os(ja)*invgj4 + (os(ia)+os(ja))*invgij4)
                      pref = dble(wr(ia))*dble(wr(ja))*dble(Q4)
                      gx = gx + pref*dble(dx4)
                      gy = gy + pref*dble(dy4)
                      gz = gz + pref*dble(dz4)
                   enddo
                enddo
             enddo
          enddo
          gv_local(ati,1) = gv_local(ati,1) + gx
          gv_local(ati,2) = gv_local(ati,2) + gy
          gv_local(ati,3) = gv_local(ati,3) + gz
       enddo
       !$omp end parallel do
    endif

    gv = gv + gv_local

    deallocate(omega0, kappa, x, y, z, act)
    deallocate(xs, ys, zs, wr, os, ks, atid)
    deallocate(gv_local)
    if (.not. noscreen) deallocate(pcx, pcy, pcz, cellstart)
end subroutine vv10_grid_force

subroutine vv10_report(Enl, used)
    implicit none
    real(8),intent(out) :: Enl
    logical,intent(out) :: used
    Enl = vv10_enl_last
    used = vv10_computed
end subroutine vv10_report

subroutine vv10_report_set(Enl)
    implicit none
    real(8),intent(in) :: Enl
    vv10_enl_last = Enl
    vv10_computed = .true.
end subroutine vv10_report_set

subroutine vv10_nlc_grid_reset()
    implicit none
    integer :: ib
    if (allocated(nlc_coor))   deallocate(nlc_coor)
    if (allocated(nlc_weight)) deallocate(nlc_weight)
    if (allocated(nlc_val_blocks)) then
       do ib = 1,size(nlc_val_blocks)
          if (allocated(nlc_val_blocks(ib)%val0)) deallocate(nlc_val_blocks(ib)%val0)
          if (allocated(nlc_val_blocks(ib)%val1)) deallocate(nlc_val_blocks(ib)%val1)
       enddo
       deallocate(nlc_val_blocks)
    endif
    if (allocated(nlc_chunk_nsig))    deallocate(nlc_chunk_nsig)
    if (allocated(nlc_chunk_sig_idx)) deallocate(nlc_chunk_sig_idx)
    nlc_n_chunks = 0
    nlc_grid_built = .false.
    nlc_basis_built = .false.
    nlc_npts = 0
end subroutine vv10_nlc_grid_reset

subroutine vv10_nlc_grid_build(npts, coor, weight)
    implicit none
    integer,intent(out) :: npts
    real(8),allocatable,intent(out) :: coor(:,:), weight(:)
    character(len=32) :: buf
    if (nlc_grid_built) then
       npts = nlc_npts
       allocate(coor(3,npts), weight(npts))
       coor = nlc_coor
       weight = nlc_weight
       return
    endif
    buf = ""
    call get_environment_variable("ENGINE_VV10_RADIAL", buf)
    if (len_trim(buf) .gt. 0) read(buf,*) nlc_nrad
    buf = ""
    call get_environment_variable("ENGINE_VV10_ANGULAR", buf)
    if (len_trim(buf) .gt. 0) read(buf,*) nlc_nsph
    call gridgen_nlc(nlc_nrad, nlc_nsph, nlc_npts, nlc_coor, nlc_weight)
    nlc_grid_built = .true.
    npts = nlc_npts
    allocate(coor(3,npts), weight(npts))
    coor = nlc_coor
    weight = nlc_weight
    print '(A,I0,A,I0,A,I0,A)', "  VV10 NLC grid points: ",nlc_npts, &
         " (",nlc_nrad," radial x ",nlc_nsph," angular/atom)"
end subroutine vv10_nlc_grid_build

subroutine vv10_nlc_grid_retrieve(npts, coor, weight, built)
    implicit none
    integer,intent(out) :: npts
    real(8),allocatable,intent(out) :: coor(:,:), weight(:)
    logical,intent(out) :: built
    built = nlc_grid_built
    npts = nlc_npts
    if (built) then
       allocate(coor(3,npts), weight(npts))
       coor = nlc_coor
       weight = nlc_weight
    endif
end subroutine vv10_nlc_grid_retrieve

subroutine vv10_nlc_basis_build(npts)
    use MOL_info, only: nconts, engine_verbose
    implicit none
    integer,intent(in) :: npts
    real(8),parameter :: SIG_CUTOFF = 1.0d-10
    integer :: ib, cstart, cn, i, k, nsig
    if (.not. nlc_basis_built) then
       nlc_n_chunks = (npts+NLC_CHUNK-1)/NLC_CHUNK
       allocate(nlc_val_blocks(nlc_n_chunks))
       allocate(nlc_chunk_nsig(nlc_n_chunks))
       allocate(nlc_chunk_sig_idx(nconts,nlc_n_chunks))
       !$omp parallel do schedule(dynamic) private(ib,cstart,cn,i,k,nsig)
       do ib = 1,nlc_n_chunks
          block
             real(8),allocatable :: val0_full(:,:), val1_full(:,:,:)
             logical,allocatable :: sig(:)
             cstart = (ib-1)*NLC_CHUNK + 1
             cn = min(NLC_CHUNK, npts-cstart+1)
             allocate(val0_full(nconts,cn), val1_full(nconts,3,cn), sig(nconts))
             do i = 1,cn
                block
                  real(8) :: val0(nconts), val1(nconts,3)
                  call GTOeval_point(nlc_coor(:,cstart+i-1), val0, val1)
                  val0_full(:,i) = val0
                  val1_full(:,:,i) = val1
                end block
             enddo
             sig = (maxval(abs(val0_full),dim=2) .gt. SIG_CUTOFF)
             nsig = count(sig)
             nlc_chunk_nsig(ib) = nsig
             k = 0
             do i = 1,nconts
                if (sig(i)) then
                   k = k+1
                   nlc_chunk_sig_idx(k,ib) = i
                endif
             enddo
             allocate(nlc_val_blocks(ib)%val0(nsig,cn), nlc_val_blocks(ib)%val1(nsig,3,cn))
             nlc_val_blocks(ib)%val0(1:nsig,1:cn) = real(val0_full(nlc_chunk_sig_idx(1:nsig,ib),1:cn),4)
             nlc_val_blocks(ib)%val1(1:nsig,:,1:cn) = real(val1_full(nlc_chunk_sig_idx(1:nsig,ib),:,1:cn),4)
             deallocate(val0_full, val1_full, sig)
          end block
       enddo
       !$omp end parallel do
       if (engine_verbose .ge. 2) then
          block
            real(8) :: avg_nsig
            avg_nsig = real(sum(nlc_chunk_nsig),8)/real(nlc_n_chunks,8)
            print '("NLC compact cache: avg nsig=",F7.1," of nconts=",I0," (",I0," chunks of ",I0," points)")', &
                  avg_nsig, nconts, nlc_n_chunks, NLC_CHUNK
            call flush(6)
          end block
       endif
       nlc_basis_built = .true.
    endif
end subroutine vv10_nlc_basis_build

subroutine vv10_nlc_chunk_info(ib, cstart, cn, nsig, sig_idx)
    implicit none
    integer,intent(in) :: ib
    integer,intent(out) :: cstart, cn, nsig
    integer,allocatable,intent(out) :: sig_idx(:)
    cstart = (ib-1)*NLC_CHUNK + 1
    cn = min(NLC_CHUNK, nlc_npts-cstart+1)
    nsig = nlc_chunk_nsig(ib)
    allocate(sig_idx(nsig))
    sig_idx = nlc_chunk_sig_idx(1:nsig,ib)
end subroutine vv10_nlc_chunk_info

subroutine vv10_nlc_chunk_val(ib, val0_ptr, val1_ptr)
    implicit none
    integer,intent(in) :: ib
    real(4),pointer,intent(out) :: val0_ptr(:,:), val1_ptr(:,:,:)
    val0_ptr => nlc_val_blocks(ib)%val0
    val1_ptr => nlc_val_blocks(ib)%val1
end subroutine vv10_nlc_chunk_val

end module mod_vv10
