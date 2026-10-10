! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_vv10_aca: ACA (adaptive cross approximation) / H-matrix acceleration

module mod_vv10_aca
implicit none
private
public :: vv10_evaluate_aca
public :: vv10_plan, vv10_aca_plan_build, vv10_aca_plan_evaluate, vv10_aca_plan_free
public :: vv10_aca_plan_evaluate_full
public :: dbg_rank_sum, dbg_rank_cnt

type :: vv10_plan
   integer,allocatable :: act(:)
   integer,allocatable :: perm(:)
   integer :: n = 0
   integer,allocatable :: near_loA(:), near_hiA(:), near_loB(:), near_hiB(:)
   logical,allocatable :: near_is_self(:)
   integer :: n_near = 0
   integer,allocatable :: far_loA(:), far_hiA(:), far_loB(:), far_hiB(:)
   integer :: n_far = 0
end type vv10_plan

integer(8),save :: dbg_rank_sum(3) = 0_8
integer(8),save :: dbg_rank_cnt(3) = 0_8

integer,parameter :: LEAF_SIZE_DEFAULT = 256
real(8),parameter :: ETA_DEFAULT = 1.5d0
real(8),parameter :: ACA_TOL_DEFAULT = 1.0d-7
integer,parameter :: ACA_MAXRANK_DEFAULT = 60
real(8),parameter :: VV10_ACA_RCUT_DEFAULT = 20.0d0
integer(8),parameter :: TASK_CUTOFF = 2000000_8
real(8),parameter :: FWEIGHT_SKIP = 1.0d-12

type :: vv10_node
   integer :: lo = 0, hi = 0
   real(8) :: lox(3) = 0.0d0, hix(3) = 0.0d0
   real(8) :: diam = 0.0d0
   logical :: is_leaf = .true.
   type(vv10_node), pointer :: left => null()
   type(vv10_node), pointer :: right => null()
end type vv10_node

contains

subroutine vv10_evaluate_aca(np, rho, gamma, weight, coor, omega0, kappa, act, n_active, &
                              beta, Enl, A_out, leaf_size_in, eta_in, tol_in, maxrank_in, &
                              rcut_in, n_far, n_near, n_self, kernel_touched)
    implicit none
    integer,intent(in) :: np, n_active
    real(8),intent(in) :: rho(np), gamma(np), weight(np), coor(3,np)
    real(8),intent(in) :: omega0(np), kappa(np), beta
    integer,intent(in) :: act(np)
    real(8),intent(out) :: Enl
    real(8),intent(out) :: A_out(np)
    integer,intent(in),optional :: leaf_size_in, maxrank_in
    real(8),intent(in),optional :: eta_in, tol_in, rcut_in
    integer,intent(out),optional :: n_far, n_near, n_self
    integer(8),intent(out),optional :: kernel_touched

    real(8),allocatable :: xs(:), ys(:), zs(:), fs(:), os(:), ks(:)
    integer,allocatable :: perm(:)
    type(vv10_node),pointer :: root
    integer :: leaf_size, maxrank, i, ia
    real(8) :: eta, tol, rcut
    real(8),allocatable :: Aacc(:)
    integer :: cnt_far, cnt_near, cnt_self
    integer(8) :: touched

    leaf_size = LEAF_SIZE_DEFAULT;  if (present(leaf_size_in)) leaf_size = leaf_size_in
    eta       = ETA_DEFAULT;        if (present(eta_in))       eta       = eta_in
    tol       = ACA_TOL_DEFAULT;    if (present(tol_in))       tol       = tol_in
    maxrank   = ACA_MAXRANK_DEFAULT;if (present(maxrank_in))   maxrank   = maxrank_in
    rcut      = VV10_ACA_RCUT_DEFAULT; if (present(rcut_in))   rcut      = rcut_in

    allocate(xs(n_active), ys(n_active), zs(n_active), fs(n_active), os(n_active), ks(n_active))
    allocate(perm(n_active), Aacc(n_active))
    do i = 1,n_active
       ia = act(i)
       xs(i) = coor(1,ia); ys(i) = coor(2,ia); zs(i) = coor(3,ia)
       fs(i) = weight(ia)*rho(ia)
       os(i) = omega0(ia); ks(i) = kappa(ia)
       perm(i) = i
    enddo
    Aacc = 0.0d0

    nullify(root)
    call build_node(root, perm, 1, n_active, xs, ys, zs, leaf_size)

    cnt_far = 0; cnt_near = 0; cnt_self = 0; touched = 0_8
    !$omp parallel
    !$omp single
    call dual_self(root, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, &
                    cnt_far, cnt_near, cnt_self, touched)
    !$omp end single
    !$omp end parallel

    Enl = beta*sum(fs)
    A_out = 0.0d0
    do i = 1,n_active
       A_out(act(i)) = Aacc(i)
       Enl = Enl + 0.5d0*fs(i)*Aacc(i)
    enddo

    if (present(n_far))  n_far  = cnt_far
    if (present(n_near)) n_near = cnt_near
    if (present(n_self)) n_self = cnt_self
    if (present(kernel_touched)) kernel_touched = touched

    call free_tree(root)
    deallocate(xs, ys, zs, fs, os, ks, perm, Aacc)
end subroutine vv10_evaluate_aca

recursive subroutine build_node(node, perm, lo, hi, xs, ys, zs, leaf_size)
    implicit none
    type(vv10_node),pointer,intent(inout) :: node
    integer,intent(inout) :: perm(:)
    integer,intent(in) :: lo, hi, leaf_size
    real(8),intent(in) :: xs(:), ys(:), zs(:)
    integer :: axis, mid
    real(8) :: span(3)

    allocate(node)
    node%lo = lo; node%hi = hi
    call bbox(perm, lo, hi, xs, ys, zs, node%lox, node%hix)
    node%diam = sqrt(sum((node%hix-node%lox)**2))

    if (hi-lo+1 .le. leaf_size) then
       node%is_leaf = .true.
       return
    endif
    node%is_leaf = .false.

    span = node%hix - node%lox
    axis = maxloc(span, dim=1)
    mid = (lo+hi)/2
    call quickselect(perm, lo, hi, mid, axis, xs, ys, zs)

    call build_node(node%left,  perm, lo,    mid, xs, ys, zs, leaf_size)
    call build_node(node%right, perm, mid+1, hi,  xs, ys, zs, leaf_size)
end subroutine build_node

subroutine bbox(perm, lo, hi, xs, ys, zs, lox, hix)
    implicit none
    integer,intent(in) :: perm(:), lo, hi
    real(8),intent(in) :: xs(:), ys(:), zs(:)
    real(8),intent(out) :: lox(3), hix(3)
    integer :: k, p
    lox = (/ xs(perm(lo)), ys(perm(lo)), zs(perm(lo)) /)
    hix = lox
    do k = lo,hi
       p = perm(k)
       lox(1) = min(lox(1), xs(p)); hix(1) = max(hix(1), xs(p))
       lox(2) = min(lox(2), ys(p)); hix(2) = max(hix(2), ys(p))
       lox(3) = min(lox(3), zs(p)); hix(3) = max(hix(3), zs(p))
    enddo
end subroutine bbox

real(8) function keyval(p, axis, xs, ys, zs)
    implicit none
    integer,intent(in) :: p, axis
    real(8),intent(in) :: xs(:), ys(:), zs(:)
    select case(axis)
    case(1); keyval = xs(p)
    case(2); keyval = ys(p)
    case default; keyval = zs(p)
    end select
end function keyval

recursive subroutine quickselect(perm, lo, hi, mid, axis, xs, ys, zs)
    implicit none
    integer,intent(inout) :: perm(:)
    integer,intent(in) :: lo, hi, mid, axis
    real(8),intent(in) :: xs(:), ys(:), zs(:)
    integer :: p, i, j, tmp
    real(8) :: pivotval
    if (lo .ge. hi) return
    p = perm((lo+hi)/2)
    pivotval = keyval(p, axis, xs, ys, zs)
    i = lo; j = hi
    do
       do while (keyval(perm(i), axis, xs, ys, zs) .lt. pivotval)
          i = i+1
       enddo
       do while (keyval(perm(j), axis, xs, ys, zs) .gt. pivotval)
          j = j-1
       enddo
       if (i .le. j) then
          tmp = perm(i); perm(i) = perm(j); perm(j) = tmp
          i = i+1; j = j-1
       endif
       if (i .gt. j) exit
    enddo
    if (mid .le. j) then
       call quickselect(perm, lo, j, mid, axis, xs, ys, zs)
    else if (mid .ge. i) then
       call quickselect(perm, i, hi, mid, axis, xs, ys, zs)
    endif
end subroutine quickselect

recursive subroutine free_tree(node)
    implicit none
    type(vv10_node),pointer,intent(inout) :: node
    if (.not. associated(node)) return
    if (.not. node%is_leaf) then
       call free_tree(node%left)
       call free_tree(node%right)
    endif
    deallocate(node)
    nullify(node)
end subroutine free_tree

logical function admissible(A, B, eta)
    implicit none
    type(vv10_node),intent(in) :: A, B
    real(8),intent(in) :: eta
    real(8) :: sep, d(3)
    d = max(0.0d0, max(A%lox-B%hix, B%lox-A%hix))
    sep = sqrt(sum(d*d))
    admissible = sep .ge. eta*max(A%diam, B%diam)
end function admissible

logical function beyond_rcut(A, B, rcut)
    implicit none
    type(vv10_node),intent(in) :: A, B
    real(8),intent(in) :: rcut
    real(8) :: sep, d(3)
    d = max(0.0d0, max(A%lox-B%hix, B%lox-A%hix))
    sep = sqrt(sum(d*d))
    beyond_rcut = sep .gt. rcut
end function beyond_rcut

real(4) function phi4(i, xA, yA, zA, oA, kA, j, xB, yB, zB, oB, kB)
    implicit none
    integer,intent(in) :: i, j
    real(4),intent(in) :: xA(:), yA(:), zA(:), oA(:), kA(:)
    real(4),intent(in) :: xB(:), yB(:), zB(:), oB(:), kB(:)
    real(4) :: dx, dy, dz, R2, gi, gj
    dx = xA(i)-xB(j); dy = yA(i)-yB(j); dz = zA(i)-zB(j)
    R2 = dx*dx+dy*dy+dz*dz
    gi = oA(i)*R2+kA(i)
    gj = oB(j)*R2+kB(j)
    phi4 = -1.5_4/(gi*gj*(gi+gj))
end function phi4

real(4) function kern_eval(kern_id, i, xA, yA, zA, oA, kA, j, xB, yB, zB, oB, kB)
    implicit none
    integer,intent(in) :: kern_id, i, j
    real(4),intent(in) :: xA(:), yA(:), zA(:), oA(:), kA(:)
    real(4),intent(in) :: xB(:), yB(:), zB(:), oB(:), kB(:)
    real(4) :: dx, dy, dz, R2, gi, gj, phi, uwfac
    dx = xA(i)-xB(j); dy = yA(i)-yB(j); dz = zA(i)-zB(j)
    R2 = dx*dx+dy*dy+dz*dz
    gi = oA(i)*R2+kA(i)
    gj = oB(j)*R2+kB(j)
    phi = -1.5_4/(gi*gj*(gi+gj))
    select case(kern_id)
    case(1)
       kern_eval = phi
    case(2)
       uwfac = 1.0_4/gi + 1.0_4/(gi+gj)
       kern_eval = phi*uwfac
    case default
       uwfac = 1.0_4/gi + 1.0_4/(gi+gj)
       kern_eval = phi*uwfac*R2
    end select
end function kern_eval

subroutine kern_pair_auw(i, xA, yA, zA, oA, kA, j, xB, yB, zB, oB, kB, val1, val2, val3)
    implicit none
    integer,intent(in) :: i, j
    real(4),intent(in) :: xA(:), yA(:), zA(:), oA(:), kA(:)
    real(4),intent(in) :: xB(:), yB(:), zB(:), oB(:), kB(:)
    real(8),intent(out) :: val1, val2, val3
    real(8) :: dx, dy, dz, R2, gi, gj, phi, uwfac
    dx = real(xA(i),8)-real(xB(j),8); dy = real(yA(i),8)-real(yB(j),8); dz = real(zA(i),8)-real(zB(j),8)
    R2 = dx*dx+dy*dy+dz*dz
    gi = real(oA(i),8)*R2+real(kA(i),8)
    gj = real(oB(j),8)*R2+real(kB(j),8)
    phi = -1.5d0/(gi*gj*(gi+gj))
    uwfac = 1.0d0/gi + 1.0d0/(gi+gj)
    val1 = phi
    val2 = phi*uwfac
    val3 = val2*R2
end subroutine kern_pair_auw

subroutine kern_pair_uw(i, xA, yA, zA, oA, kA, j, xB, yB, zB, oB, kB, val2, val3)
    implicit none
    integer,intent(in) :: i, j
    real(4),intent(in) :: xA(:), yA(:), zA(:), oA(:), kA(:)
    real(4),intent(in) :: xB(:), yB(:), zB(:), oB(:), kB(:)
    real(8),intent(out) :: val2, val3
    real(8) :: dx, dy, dz, R2, gi, gj, phi, uwfac
    dx = real(xA(i),8)-real(xB(j),8); dy = real(yA(i),8)-real(yB(j),8); dz = real(zA(i),8)-real(zB(j),8)
    R2 = dx*dx+dy*dy+dz*dz
    gi = real(oA(i),8)*R2+real(kA(i),8)
    gj = real(oB(j),8)*R2+real(kB(j),8)
    phi = -1.5d0/(gi*gj*(gi+gj))
    uwfac = 1.0d0/gi + 1.0d0/(gi+gj)
    val2 = phi*uwfac
    val3 = val2*R2
end subroutine kern_pair_uw

subroutine gather_block(idx, n, xs, ys, zs, os, ks, fs, xg, yg, zg, og, kg, fg)
    implicit none
    integer,intent(in) :: idx(n), n
    real(8),intent(in) :: xs(:), ys(:), zs(:), os(:), ks(:), fs(:)
    real(4),allocatable,intent(out) :: xg(:), yg(:), zg(:), og(:), kg(:), fg(:)
    integer :: k, p
    allocate(xg(n), yg(n), zg(n), og(n), kg(n), fg(n))
    do k = 1,n
       p = idx(k)
       xg(k) = real(xs(p),4); yg(k) = real(ys(p),4); zg(k) = real(zs(p),4)
       og(k) = real(os(p),4); kg(k) = real(ks(p),4); fg(k) = real(fs(p),4)
    enddo
end subroutine gather_block

subroutine aca_compress(nA, xA, yA, zA, oA, kA, nB, xB, yB, zB, oB, kB, tol, maxrank, U, V, rank)
    implicit none
    integer,intent(in) :: nA, nB, maxrank
    real(4),intent(in) :: xA(:), yA(:), zA(:), oA(:), kA(:)
    real(4),intent(in) :: xB(:), yB(:), zB(:), oB(:), kB(:)
    real(8),intent(in) :: tol
    real(8),allocatable,intent(out) :: U(:,:), V(:,:)
    integer,intent(out) :: rank
    real(8),allocatable :: row(:), col(:), Utmp(:,:), Vtmp(:,:)
    logical,allocatable :: used(:)
    integer :: k, i_local, j_local, m, r, pidx
    real(8) :: piv, nrmu, nrmv, fro2, term, cross

    allocate(Utmp(nA,maxrank), Vtmp(nB,maxrank))
    allocate(row(nB), col(nA), used(nA))
    used = .false.
    i_local = 1
    r = 0
    fro2 = 0.0d0
    do k = 1,min(maxrank,nA,nB)
       do j_local = 1,nB
          row(j_local) = real(phi4(i_local, xA, yA, zA, oA, kA, j_local, xB, yB, zB, oB, kB),8)
       enddo
       do m = 1,r
          row = row - Utmp(i_local,m)*Vtmp(:,m)
       enddo
       j_local = maxloc(abs(row), dim=1)
       piv = row(j_local)
       if (abs(piv) .lt. 1.0d-300) then
          used(i_local) = .true.
          i_local = 0
          do m = 1,nA
             if (.not. used(m)) then; i_local = m; exit; endif
          enddo
          if (i_local .eq. 0) exit
          cycle
       endif
       r = r+1
       Vtmp(:,r) = row/piv
       do pidx = 1,nA
          col(pidx) = real(phi4(pidx, xA, yA, zA, oA, kA, j_local, xB, yB, zB, oB, kB),8)
       enddo
       do m = 1,r-1
          col = col - Vtmp(j_local,m)*Utmp(:,m)
       enddo
       Utmp(:,r) = col
       nrmu = sqrt(sum(Utmp(:,r)**2)); nrmv = sqrt(sum(Vtmp(:,r)**2))
       term = nrmu*nrmv
       cross = 0.0d0
       do m = 1,r-1
          cross = cross + 2.0d0*sum(Utmp(:,r)*Utmp(:,m))*sum(Vtmp(:,r)*Vtmp(:,m))
       enddo
       fro2 = fro2 + term*term + cross
       used(i_local) = .true.
       if (fro2 .gt. 0.0d0 .and. term .le. tol*sqrt(max(fro2,term*term))) exit
       i_local = 0
       piv = -1.0d0
       do m = 1,nA
          if (.not. used(m) .and. abs(Utmp(m,r)) .gt. piv) then
             piv = abs(Utmp(m,r)); i_local = m
          endif
       enddo
       if (i_local .eq. 0) exit
    enddo
    rank = r
    allocate(U(nA,max(rank,1)), V(nB,max(rank,1)))
    if (rank .gt. 0) then
       U(:,1:rank) = Utmp(:,1:rank)
       V(:,1:rank) = Vtmp(:,1:rank)
    endif
    deallocate(Utmp, Vtmp, row, col, used)
end subroutine aca_compress

subroutine aca_compress_k(kern_id, nA, xA, yA, zA, oA, kA, nB, xB, yB, zB, oB, kB, tol, maxrank, U, V, rank)
    implicit none
    integer,intent(in) :: kern_id, nA, nB, maxrank
    real(4),intent(in) :: xA(:), yA(:), zA(:), oA(:), kA(:)
    real(4),intent(in) :: xB(:), yB(:), zB(:), oB(:), kB(:)
    real(8),intent(in) :: tol
    real(8),allocatable,intent(out) :: U(:,:), V(:,:)
    integer,intent(out) :: rank
    real(8),allocatable :: row(:), col(:), Utmp(:,:), Vtmp(:,:)
    logical,allocatable :: used(:)
    integer :: k, i_local, j_local, m, r, pidx
    real(8) :: piv, nrmu, nrmv, fro2, term, cross

    allocate(Utmp(nA,maxrank), Vtmp(nB,maxrank))
    allocate(row(nB), col(nA), used(nA))
    used = .false.
    i_local = 1
    r = 0
    fro2 = 0.0d0
    do k = 1,min(maxrank,nA,nB)
       do j_local = 1,nB
          row(j_local) = real(kern_eval(kern_id, i_local, xA, yA, zA, oA, kA, j_local, xB, yB, zB, oB, kB),8)
       enddo
       do m = 1,r
          row = row - Utmp(i_local,m)*Vtmp(:,m)
       enddo
       j_local = maxloc(abs(row), dim=1)
       piv = row(j_local)
       if (abs(piv) .lt. 1.0d-300) then
          used(i_local) = .true.
          i_local = 0
          do m = 1,nA
             if (.not. used(m)) then; i_local = m; exit; endif
          enddo
          if (i_local .eq. 0) exit
          cycle
       endif
       r = r+1
       Vtmp(:,r) = row/piv
       do pidx = 1,nA
          col(pidx) = real(kern_eval(kern_id, pidx, xA, yA, zA, oA, kA, j_local, xB, yB, zB, oB, kB),8)
       enddo
       do m = 1,r-1
          col = col - Vtmp(j_local,m)*Utmp(:,m)
       enddo
       Utmp(:,r) = col
       nrmu = sqrt(sum(Utmp(:,r)**2)); nrmv = sqrt(sum(Vtmp(:,r)**2))
       term = nrmu*nrmv
       cross = 0.0d0
       do m = 1,r-1
          cross = cross + 2.0d0*sum(Utmp(:,r)*Utmp(:,m))*sum(Vtmp(:,r)*Vtmp(:,m))
       enddo
       fro2 = fro2 + term*term + cross
       used(i_local) = .true.
       if (fro2 .gt. 0.0d0 .and. term .le. tol*sqrt(max(fro2,term*term))) exit
       i_local = 0
       piv = -1.0d0
       do m = 1,nA
          if (.not. used(m) .and. abs(Utmp(m,r)) .gt. piv) then
             piv = abs(Utmp(m,r)); i_local = m
          endif
       enddo
       if (i_local .eq. 0) exit
    enddo
    rank = r
    allocate(U(nA,max(rank,1)), V(nB,max(rank,1)))
    if (rank .gt. 0) then
       U(:,1:rank) = Utmp(:,1:rank)
       V(:,1:rank) = Vtmp(:,1:rank)
    endif
    deallocate(Utmp, Vtmp, row, col, used)
end subroutine aca_compress_k

subroutine aca_compress_uw(nA, xA, yA, zA, oA, kA, nB, xB, yB, zB, oB, kB, tol, maxrank, &
                            U2, V2, rank2, U3, V3, rank3)
    implicit none
    integer,intent(in) :: nA, nB, maxrank
    real(4),intent(in) :: xA(:), yA(:), zA(:), oA(:), kA(:)
    real(4),intent(in) :: xB(:), yB(:), zB(:), oB(:), kB(:)
    real(8),intent(in) :: tol
    real(8),allocatable,intent(out) :: U2(:,:), V2(:,:), U3(:,:), V3(:,:)
    integer,intent(out) :: rank2, rank3
    real(8),allocatable :: row2(:), col2(:), Utmp2(:,:), Vtmp2(:,:)
    real(8),allocatable :: row3(:), col3(:), Utmp3(:,:), Vtmp3(:,:)
    logical,allocatable :: used2(:), used3(:)
    real(8),allocatable :: rowcache2(:,:), rowcache3(:,:), colcache2(:,:), colcache3(:,:)
    integer,allocatable :: row_slot(:), col_slot(:)
    integer :: cache_cap, n_row_slot, n_col_slot, slot
    integer :: k, i2, j2, i3, j3, m, r2, r3, pidx, active
    real(8) :: piv2, piv3, nrmu, nrmv, fro2_2, fro2_3, term2, term3, cross
    logical :: conv2, conv3
    real(8) :: v2tmp, v3tmp

    allocate(Utmp2(nA,maxrank), Vtmp2(nB,maxrank), Utmp3(nA,maxrank), Vtmp3(nB,maxrank))
    allocate(row2(nB), col2(nA), row3(nB), col3(nA), used2(nA), used3(nA))
    cache_cap = 2*maxrank + 2
    allocate(rowcache2(nB,cache_cap), rowcache3(nB,cache_cap))
    allocate(colcache2(nA,cache_cap), colcache3(nA,cache_cap))
    allocate(row_slot(nA), col_slot(nB))
    row_slot = 0; col_slot = 0
    n_row_slot = 0; n_col_slot = 0
    used2 = .false.; used3 = .false.
    i2 = 1; i3 = 1
    r2 = 0; r3 = 0
    fro2_2 = 0.0d0; fro2_3 = 0.0d0
    conv2 = .false.; conv3 = .false.
    active = min(maxrank,nA,nB)
    do k = 1,active
       if (conv2 .and. conv3) exit
       if (.not. conv2) then
          if (row_slot(i2) .eq. 0) then
             n_row_slot = n_row_slot + 1
             row_slot(i2) = n_row_slot
             do j2 = 1,nB
                call kern_pair_uw(i2, xA, yA, zA, oA, kA, j2, xB, yB, zB, oB, kB, v2tmp, v3tmp)
                rowcache2(j2,n_row_slot) = v2tmp; rowcache3(j2,n_row_slot) = v3tmp
             enddo
          endif
          row2 = rowcache2(:,row_slot(i2))
          do m = 1,r2
             row2 = row2 - Utmp2(i2,m)*Vtmp2(:,m)
          enddo
          j2 = maxloc(abs(row2), dim=1)
          piv2 = row2(j2)
          if (abs(piv2) .lt. 1.0d-300) then
             used2(i2) = .true.
             i2 = 0
             do m = 1,nA
                if (.not. used2(m)) then; i2 = m; exit; endif
             enddo
             if (i2 .eq. 0) then; conv2 = .true.; endif
          else
             r2 = r2+1
             Vtmp2(:,r2) = row2/piv2
             if (col_slot(j2) .eq. 0) then
                n_col_slot = n_col_slot + 1
                col_slot(j2) = n_col_slot
                do pidx = 1,nA
                   call kern_pair_uw(pidx, xA, yA, zA, oA, kA, j2, xB, yB, zB, oB, kB, v2tmp, v3tmp)
                   colcache2(pidx,n_col_slot) = v2tmp; colcache3(pidx,n_col_slot) = v3tmp
                enddo
             endif
             col2 = colcache2(:,col_slot(j2))
             do m = 1,r2-1
                col2 = col2 - Vtmp2(j2,m)*Utmp2(:,m)
             enddo
             Utmp2(:,r2) = col2
             nrmu = sqrt(sum(Utmp2(:,r2)**2)); nrmv = sqrt(sum(Vtmp2(:,r2)**2))
             term2 = nrmu*nrmv
             cross = 0.0d0
             do m = 1,r2-1
                cross = cross + 2.0d0*sum(Utmp2(:,r2)*Utmp2(:,m))*sum(Vtmp2(:,r2)*Vtmp2(:,m))
             enddo
             fro2_2 = fro2_2 + term2*term2 + cross
             used2(i2) = .true.
             if (fro2_2 .gt. 0.0d0 .and. term2 .le. tol*sqrt(max(fro2_2,term2*term2))) then
                conv2 = .true.
             else
                i2 = 0; piv2 = -1.0d0
                do m = 1,nA
                   if (.not. used2(m) .and. abs(Utmp2(m,r2)) .gt. piv2) then
                      piv2 = abs(Utmp2(m,r2)); i2 = m
                   endif
                enddo
                if (i2 .eq. 0) conv2 = .true.
             endif
          endif
       endif
       if (r2 .ge. maxrank) conv2 = .true.
       if (.not. conv3) then
          if (row_slot(i3) .eq. 0) then
             n_row_slot = n_row_slot + 1
             row_slot(i3) = n_row_slot
             do j3 = 1,nB
                call kern_pair_uw(i3, xA, yA, zA, oA, kA, j3, xB, yB, zB, oB, kB, v2tmp, v3tmp)
                rowcache2(j3,n_row_slot) = v2tmp; rowcache3(j3,n_row_slot) = v3tmp
             enddo
          endif
          row3 = rowcache3(:,row_slot(i3))
          do m = 1,r3
             row3 = row3 - Utmp3(i3,m)*Vtmp3(:,m)
          enddo
          j3 = maxloc(abs(row3), dim=1)
          piv3 = row3(j3)
          if (abs(piv3) .lt. 1.0d-300) then
             used3(i3) = .true.
             i3 = 0
             do m = 1,nA
                if (.not. used3(m)) then; i3 = m; exit; endif
             enddo
             if (i3 .eq. 0) then; conv3 = .true.; endif
          else
             r3 = r3+1
             Vtmp3(:,r3) = row3/piv3
             if (col_slot(j3) .eq. 0) then
                n_col_slot = n_col_slot + 1
                col_slot(j3) = n_col_slot
                do pidx = 1,nA
                   call kern_pair_uw(pidx, xA, yA, zA, oA, kA, j3, xB, yB, zB, oB, kB, v2tmp, v3tmp)
                   colcache2(pidx,n_col_slot) = v2tmp; colcache3(pidx,n_col_slot) = v3tmp
                enddo
             endif
             col3 = colcache3(:,col_slot(j3))
             do m = 1,r3-1
                col3 = col3 - Vtmp3(j3,m)*Utmp3(:,m)
             enddo
             Utmp3(:,r3) = col3
             nrmu = sqrt(sum(Utmp3(:,r3)**2)); nrmv = sqrt(sum(Vtmp3(:,r3)**2))
             term3 = nrmu*nrmv
             cross = 0.0d0
             do m = 1,r3-1
                cross = cross + 2.0d0*sum(Utmp3(:,r3)*Utmp3(:,m))*sum(Vtmp3(:,r3)*Vtmp3(:,m))
             enddo
             fro2_3 = fro2_3 + term3*term3 + cross
             used3(i3) = .true.
             if (fro2_3 .gt. 0.0d0 .and. term3 .le. tol*sqrt(max(fro2_3,term3*term3))) then
                conv3 = .true.
             else
                i3 = 0; piv3 = -1.0d0
                do m = 1,nA
                   if (.not. used3(m) .and. abs(Utmp3(m,r3)) .gt. piv3) then
                      piv3 = abs(Utmp3(m,r3)); i3 = m
                   endif
                enddo
                if (i3 .eq. 0) conv3 = .true.
             endif
          endif
       endif
       if (r3 .ge. maxrank) conv3 = .true.
    enddo
    rank2 = r2
    rank3 = r3
    allocate(U2(nA,max(rank2,1)), V2(nB,max(rank2,1)))
    allocate(U3(nA,max(rank3,1)), V3(nB,max(rank3,1)))
    if (rank2 .gt. 0) then
       U2(:,1:rank2) = Utmp2(:,1:rank2)
       V2(:,1:rank2) = Vtmp2(:,1:rank2)
    endif
    if (rank3 .gt. 0) then
       U3(:,1:rank3) = Utmp3(:,1:rank3)
       V3(:,1:rank3) = Vtmp3(:,1:rank3)
    endif
    deallocate(Utmp2, Vtmp2, Utmp3, Vtmp3, row2, col2, row3, col3, used2, used3)
    deallocate(rowcache2, rowcache3, colcache2, colcache3, row_slot, col_slot)
end subroutine aca_compress_uw

subroutine aca_compress_auw(nA, xA, yA, zA, oA, kA, nB, xB, yB, zB, oB, kB, tol, maxrank, &
                             U1, V1, rank1, U2, V2, rank2, U3, V3, rank3)
    implicit none
    integer,intent(in) :: nA, nB, maxrank
    real(4),intent(in) :: xA(:), yA(:), zA(:), oA(:), kA(:)
    real(4),intent(in) :: xB(:), yB(:), zB(:), oB(:), kB(:)
    real(8),intent(in) :: tol
    real(8),allocatable,intent(out) :: U1(:,:), V1(:,:), U2(:,:), V2(:,:), U3(:,:), V3(:,:)
    integer,intent(out) :: rank1, rank2, rank3
    real(8),allocatable :: row1(:), col1(:), Utmp1(:,:), Vtmp1(:,:)
    real(8),allocatable :: row2(:), col2(:), Utmp2(:,:), Vtmp2(:,:)
    real(8),allocatable :: row3(:), col3(:), Utmp3(:,:), Vtmp3(:,:)
    logical,allocatable :: used1(:), used2(:), used3(:)
    real(8),allocatable :: rowcache1(:,:), rowcache2(:,:), rowcache3(:,:)
    real(8),allocatable :: colcache1(:,:), colcache2(:,:), colcache3(:,:)
    integer,allocatable :: row_slot(:), col_slot(:)
    integer :: cache_cap, n_row_slot, n_col_slot
    integer :: k, i1, j1, i2, j2, i3, j3, m, r1, r2, r3, pidx, active
    real(8) :: piv1, piv2, piv3, nrmu, nrmv, fro2_1, fro2_2, fro2_3, term1, term2, term3, cross
    logical :: conv1, conv2, conv3
    real(8) :: v1tmp, v2tmp, v3tmp

    allocate(Utmp1(nA,maxrank), Vtmp1(nB,maxrank))
    allocate(Utmp2(nA,maxrank), Vtmp2(nB,maxrank), Utmp3(nA,maxrank), Vtmp3(nB,maxrank))
    allocate(row1(nB), col1(nA), row2(nB), col2(nA), row3(nB), col3(nA))
    allocate(used1(nA), used2(nA), used3(nA))
    cache_cap = 3*maxrank + 2
    allocate(rowcache1(nB,cache_cap), rowcache2(nB,cache_cap), rowcache3(nB,cache_cap))
    allocate(colcache1(nA,cache_cap), colcache2(nA,cache_cap), colcache3(nA,cache_cap))
    allocate(row_slot(nA), col_slot(nB))
    row_slot = 0; col_slot = 0
    n_row_slot = 0; n_col_slot = 0
    used1 = .false.; used2 = .false.; used3 = .false.
    i1 = 1; i2 = 1; i3 = 1
    r1 = 0; r2 = 0; r3 = 0
    fro2_1 = 0.0d0; fro2_2 = 0.0d0; fro2_3 = 0.0d0
    conv1 = .false.; conv2 = .false.; conv3 = .false.
    active = min(maxrank,nA,nB)
    do k = 1,active
       if (conv1 .and. conv2 .and. conv3) exit
       if (.not. conv1) then
          if (row_slot(i1) .eq. 0) then
             n_row_slot = n_row_slot + 1
             row_slot(i1) = n_row_slot
             do j1 = 1,nB
                call kern_pair_auw(i1, xA, yA, zA, oA, kA, j1, xB, yB, zB, oB, kB, v1tmp, v2tmp, v3tmp)
                rowcache1(j1,n_row_slot) = v1tmp; rowcache2(j1,n_row_slot) = v2tmp; rowcache3(j1,n_row_slot) = v3tmp
             enddo
          endif
          row1 = rowcache1(:,row_slot(i1))
          do m = 1,r1
             row1 = row1 - Utmp1(i1,m)*Vtmp1(:,m)
          enddo
          j1 = maxloc(abs(row1), dim=1)
          piv1 = row1(j1)
          if (abs(piv1) .lt. 1.0d-300) then
             used1(i1) = .true.
             i1 = 0
             do m = 1,nA
                if (.not. used1(m)) then; i1 = m; exit; endif
             enddo
             if (i1 .eq. 0) then; conv1 = .true.; endif
          else
             r1 = r1+1
             Vtmp1(:,r1) = row1/piv1
             if (col_slot(j1) .eq. 0) then
                n_col_slot = n_col_slot + 1
                col_slot(j1) = n_col_slot
                do pidx = 1,nA
                   call kern_pair_auw(pidx, xA, yA, zA, oA, kA, j1, xB, yB, zB, oB, kB, v1tmp, v2tmp, v3tmp)
                   colcache1(pidx,n_col_slot) = v1tmp; colcache2(pidx,n_col_slot) = v2tmp; colcache3(pidx,n_col_slot) = v3tmp
                enddo
             endif
             col1 = colcache1(:,col_slot(j1))
             do m = 1,r1-1
                col1 = col1 - Vtmp1(j1,m)*Utmp1(:,m)
             enddo
             Utmp1(:,r1) = col1
             nrmu = sqrt(sum(Utmp1(:,r1)**2)); nrmv = sqrt(sum(Vtmp1(:,r1)**2))
             term1 = nrmu*nrmv
             cross = 0.0d0
             do m = 1,r1-1
                cross = cross + 2.0d0*sum(Utmp1(:,r1)*Utmp1(:,m))*sum(Vtmp1(:,r1)*Vtmp1(:,m))
             enddo
             fro2_1 = fro2_1 + term1*term1 + cross
             used1(i1) = .true.
             if (fro2_1 .gt. 0.0d0 .and. term1 .le. tol*sqrt(max(fro2_1,term1*term1))) then
                conv1 = .true.
             else
                i1 = 0; piv1 = -1.0d0
                do m = 1,nA
                   if (.not. used1(m) .and. abs(Utmp1(m,r1)) .gt. piv1) then
                      piv1 = abs(Utmp1(m,r1)); i1 = m
                   endif
                enddo
                if (i1 .eq. 0) conv1 = .true.
             endif
          endif
       endif
       if (r1 .ge. maxrank) conv1 = .true.
       if (.not. conv2) then
          if (row_slot(i2) .eq. 0) then
             n_row_slot = n_row_slot + 1
             row_slot(i2) = n_row_slot
             do j2 = 1,nB
                call kern_pair_auw(i2, xA, yA, zA, oA, kA, j2, xB, yB, zB, oB, kB, v1tmp, v2tmp, v3tmp)
                rowcache1(j2,n_row_slot) = v1tmp; rowcache2(j2,n_row_slot) = v2tmp; rowcache3(j2,n_row_slot) = v3tmp
             enddo
          endif
          row2 = rowcache2(:,row_slot(i2))
          do m = 1,r2
             row2 = row2 - Utmp2(i2,m)*Vtmp2(:,m)
          enddo
          j2 = maxloc(abs(row2), dim=1)
          piv2 = row2(j2)
          if (abs(piv2) .lt. 1.0d-300) then
             used2(i2) = .true.
             i2 = 0
             do m = 1,nA
                if (.not. used2(m)) then; i2 = m; exit; endif
             enddo
             if (i2 .eq. 0) then; conv2 = .true.; endif
          else
             r2 = r2+1
             Vtmp2(:,r2) = row2/piv2
             if (col_slot(j2) .eq. 0) then
                n_col_slot = n_col_slot + 1
                col_slot(j2) = n_col_slot
                do pidx = 1,nA
                   call kern_pair_auw(pidx, xA, yA, zA, oA, kA, j2, xB, yB, zB, oB, kB, v1tmp, v2tmp, v3tmp)
                   colcache1(pidx,n_col_slot) = v1tmp; colcache2(pidx,n_col_slot) = v2tmp; colcache3(pidx,n_col_slot) = v3tmp
                enddo
             endif
             col2 = colcache2(:,col_slot(j2))
             do m = 1,r2-1
                col2 = col2 - Vtmp2(j2,m)*Utmp2(:,m)
             enddo
             Utmp2(:,r2) = col2
             nrmu = sqrt(sum(Utmp2(:,r2)**2)); nrmv = sqrt(sum(Vtmp2(:,r2)**2))
             term2 = nrmu*nrmv
             cross = 0.0d0
             do m = 1,r2-1
                cross = cross + 2.0d0*sum(Utmp2(:,r2)*Utmp2(:,m))*sum(Vtmp2(:,r2)*Vtmp2(:,m))
             enddo
             fro2_2 = fro2_2 + term2*term2 + cross
             used2(i2) = .true.
             if (fro2_2 .gt. 0.0d0 .and. term2 .le. tol*sqrt(max(fro2_2,term2*term2))) then
                conv2 = .true.
             else
                i2 = 0; piv2 = -1.0d0
                do m = 1,nA
                   if (.not. used2(m) .and. abs(Utmp2(m,r2)) .gt. piv2) then
                      piv2 = abs(Utmp2(m,r2)); i2 = m
                   endif
                enddo
                if (i2 .eq. 0) conv2 = .true.
             endif
          endif
       endif
       if (r2 .ge. maxrank) conv2 = .true.
       if (.not. conv3) then
          if (row_slot(i3) .eq. 0) then
             n_row_slot = n_row_slot + 1
             row_slot(i3) = n_row_slot
             do j3 = 1,nB
                call kern_pair_auw(i3, xA, yA, zA, oA, kA, j3, xB, yB, zB, oB, kB, v1tmp, v2tmp, v3tmp)
                rowcache1(j3,n_row_slot) = v1tmp; rowcache2(j3,n_row_slot) = v2tmp; rowcache3(j3,n_row_slot) = v3tmp
             enddo
          endif
          row3 = rowcache3(:,row_slot(i3))
          do m = 1,r3
             row3 = row3 - Utmp3(i3,m)*Vtmp3(:,m)
          enddo
          j3 = maxloc(abs(row3), dim=1)
          piv3 = row3(j3)
          if (abs(piv3) .lt. 1.0d-300) then
             used3(i3) = .true.
             i3 = 0
             do m = 1,nA
                if (.not. used3(m)) then; i3 = m; exit; endif
             enddo
             if (i3 .eq. 0) then; conv3 = .true.; endif
          else
             r3 = r3+1
             Vtmp3(:,r3) = row3/piv3
             if (col_slot(j3) .eq. 0) then
                n_col_slot = n_col_slot + 1
                col_slot(j3) = n_col_slot
                do pidx = 1,nA
                   call kern_pair_auw(pidx, xA, yA, zA, oA, kA, j3, xB, yB, zB, oB, kB, v1tmp, v2tmp, v3tmp)
                   colcache1(pidx,n_col_slot) = v1tmp; colcache2(pidx,n_col_slot) = v2tmp; colcache3(pidx,n_col_slot) = v3tmp
                enddo
             endif
             col3 = colcache3(:,col_slot(j3))
             do m = 1,r3-1
                col3 = col3 - Vtmp3(j3,m)*Utmp3(:,m)
             enddo
             Utmp3(:,r3) = col3
             nrmu = sqrt(sum(Utmp3(:,r3)**2)); nrmv = sqrt(sum(Vtmp3(:,r3)**2))
             term3 = nrmu*nrmv
             cross = 0.0d0
             do m = 1,r3-1
                cross = cross + 2.0d0*sum(Utmp3(:,r3)*Utmp3(:,m))*sum(Vtmp3(:,r3)*Vtmp3(:,m))
             enddo
             fro2_3 = fro2_3 + term3*term3 + cross
             used3(i3) = .true.
             if (fro2_3 .gt. 0.0d0 .and. term3 .le. tol*sqrt(max(fro2_3,term3*term3))) then
                conv3 = .true.
             else
                i3 = 0; piv3 = -1.0d0
                do m = 1,nA
                   if (.not. used3(m) .and. abs(Utmp3(m,r3)) .gt. piv3) then
                      piv3 = abs(Utmp3(m,r3)); i3 = m
                   endif
                enddo
                if (i3 .eq. 0) conv3 = .true.
             endif
          endif
       endif
       if (r3 .ge. maxrank) conv3 = .true.
    enddo
    rank1 = r1; rank2 = r2; rank3 = r3
    allocate(U1(nA,max(rank1,1)), V1(nB,max(rank1,1)))
    allocate(U2(nA,max(rank2,1)), V2(nB,max(rank2,1)))
    allocate(U3(nA,max(rank3,1)), V3(nB,max(rank3,1)))
    if (rank1 .gt. 0) then
       U1(:,1:rank1) = Utmp1(:,1:rank1)
       V1(:,1:rank1) = Vtmp1(:,1:rank1)
    endif
    if (rank2 .gt. 0) then
       U2(:,1:rank2) = Utmp2(:,1:rank2)
       V2(:,1:rank2) = Vtmp2(:,1:rank2)
    endif
    if (rank3 .gt. 0) then
       U3(:,1:rank3) = Utmp3(:,1:rank3)
       V3(:,1:rank3) = Vtmp3(:,1:rank3)
    endif
    deallocate(Utmp1, Vtmp1, Utmp2, Vtmp2, Utmp3, Vtmp3, row1, col1, row2, col2, row3, col3)
    deallocate(used1, used2, used3)
    deallocate(rowcache1, rowcache2, rowcache3, colcache1, colcache2, colcache3, row_slot, col_slot)
end subroutine aca_compress_auw

subroutine far_block(node_perm, iAloc, nA, iBloc, nB, xs, ys, zs, os, ks, fs, Aacc, &
                      tol, maxrank, touched)
    implicit none
    integer,intent(in) :: node_perm(:)
    integer,intent(in) :: iAloc(nA), iBloc(nB), nA, nB, maxrank
    real(8),intent(in) :: xs(:), ys(:), zs(:), os(:), ks(:), fs(:), tol
    real(8),intent(inout) :: Aacc(:)
    integer(8),intent(inout) :: touched
    real(8),allocatable :: U(:,:), V(:,:), UtfB(:), VtfA(:)
    real(4),allocatable :: xA(:), yA(:), zA(:), oA(:), kA(:), fA4(:)
    real(4),allocatable :: xB(:), yB(:), zB(:), oB(:), kB(:), fB4(:)
    real(8),allocatable :: fa8(:), fb8(:)
    integer :: rank, k
    call gather_block(iAloc, nA, xs, ys, zs, os, ks, fs, xA, yA, zA, oA, kA, fA4)
    call gather_block(iBloc, nB, xs, ys, zs, os, ks, fs, xB, yB, zB, oB, kB, fB4)
    if (maxval(abs(fA4)) .lt. FWEIGHT_SKIP .or. maxval(abs(fB4)) .lt. FWEIGHT_SKIP) then
       deallocate(xA,yA,zA,oA,kA,fA4,xB,yB,zB,oB,kB,fB4)
       return
    endif
    call aca_compress(nA, xA, yA, zA, oA, kA, nB, xB, yB, zB, oB, kB, tol, maxrank, U, V, rank)
    if (rank .eq. 0) then
       deallocate(xA,yA,zA,oA,kA,fA4,xB,yB,zB,oB,kB,fB4)
       return
    endif
    allocate(fb8(nB), fa8(nA), UtfB(rank), VtfA(rank))
    fb8 = real(fB4,8); fa8 = real(fA4,8)
    UtfB = matmul(transpose(V), fb8)
    do k = 1,nA
       !$omp atomic update
       Aacc(iAloc(k)) = Aacc(iAloc(k)) + sum(U(k,1:rank)*UtfB(1:rank))
    enddo
    VtfA = matmul(transpose(U), fa8)
    do k = 1,nB
       !$omp atomic update
       Aacc(iBloc(k)) = Aacc(iBloc(k)) + sum(V(k,1:rank)*VtfA(1:rank))
    enddo
    !$omp atomic update
    touched = touched + int(rank,8)*int(nA+nB,8)
    deallocate(U, V, fb8, fa8, UtfB, VtfA, xA,yA,zA,oA,kA,fA4,xB,yB,zB,oB,kB,fB4)
end subroutine far_block

subroutine far_block_full(iAloc, nA, iBloc, nB, xs, ys, zs, os, ks, fs, Aacc, Uacc, Wacc, &
                           tol, maxrank, touched)
    implicit none
    integer,intent(in) :: iAloc(nA), iBloc(nB), nA, nB, maxrank
    real(8),intent(in) :: xs(:), ys(:), zs(:), os(:), ks(:), fs(:), tol
    real(8),intent(inout) :: Aacc(:), Uacc(:), Wacc(:)
    integer(8),intent(inout) :: touched
    real(8),allocatable :: UU(:,:), VV(:,:), UtfB(:), VtfA(:)
    real(4),allocatable :: xA(:), yA(:), zA(:), oA(:), kA(:), fA4(:)
    real(4),allocatable :: xB(:), yB(:), zB(:), oB(:), kB(:), fB4(:)
    real(8),allocatable :: fa8(:), fb8(:)
    integer :: rank, k
    call gather_block(iAloc, nA, xs, ys, zs, os, ks, fs, xA, yA, zA, oA, kA, fA4)
    call gather_block(iBloc, nB, xs, ys, zs, os, ks, fs, xB, yB, zB, oB, kB, fB4)
    if (maxval(abs(fA4)) .lt. FWEIGHT_SKIP .or. maxval(abs(fB4)) .lt. FWEIGHT_SKIP) then
       deallocate(xA,yA,zA,oA,kA,fA4,xB,yB,zB,oB,kB,fB4)
       return
    endif
    allocate(fb8(nB), fa8(nA))
    fb8 = real(fB4,8); fa8 = real(fA4,8)

    block
      real(8),allocatable :: U3A(:,:), V3A(:,:), U3(:,:), V3(:,:), UtfB3(:), VtfA3(:)
      real(8),allocatable :: VtfC2(:), VtfC3(:)
      integer :: rankA1, rankA2, rankA3
      call aca_compress_auw(nA, xA,yA,zA,oA,kA, nB, xB,yB,zB,oB,kB, tol, maxrank, &
                             UU, VV, rankA1, U3A, V3A, rankA2, U3, V3, rankA3)
      !$omp atomic update
      dbg_rank_sum(1) = dbg_rank_sum(1) + int(rankA1,8)
      !$omp atomic update
      dbg_rank_cnt(1) = dbg_rank_cnt(1) + 1_8
      !$omp atomic update
      dbg_rank_sum(2) = dbg_rank_sum(2) + int(rankA2,8)
      !$omp atomic update
      dbg_rank_cnt(2) = dbg_rank_cnt(2) + 1_8
      !$omp atomic update
      dbg_rank_sum(3) = dbg_rank_sum(3) + int(rankA3,8)
      !$omp atomic update
      dbg_rank_cnt(3) = dbg_rank_cnt(3) + 1_8
      if (rankA1 .gt. 0) then
         allocate(UtfB3(rankA1), VtfA3(rankA1))
         UtfB3 = matmul(transpose(VV), fb8)
         do k = 1,nA
            !$omp atomic update
            Aacc(iAloc(k)) = Aacc(iAloc(k)) + sum(UU(k,1:rankA1)*UtfB3(1:rankA1))
         enddo
         VtfA3 = matmul(transpose(UU), fa8)
         do k = 1,nB
            !$omp atomic update
            Aacc(iBloc(k)) = Aacc(iBloc(k)) + sum(VV(k,1:rankA1)*VtfA3(1:rankA1))
         enddo
         deallocate(UtfB3, VtfA3)
      endif
      if (rankA2 .gt. 0) then
         allocate(VtfC2(rankA2))
         VtfC2 = matmul(transpose(V3A), fb8)
         do k = 1,nA
            !$omp atomic update
            Uacc(iAloc(k)) = Uacc(iAloc(k)) + sum(U3A(k,1:rankA2)*VtfC2(1:rankA2))
         enddo
         deallocate(VtfC2)
      endif
      if (rankA3 .gt. 0) then
         allocate(VtfC3(rankA3))
         VtfC3 = matmul(transpose(V3), fb8)
         do k = 1,nA
            !$omp atomic update
            Wacc(iAloc(k)) = Wacc(iAloc(k)) + sum(U3(k,1:rankA3)*VtfC3(1:rankA3))
         enddo
         deallocate(VtfC3)
      endif
      !$omp atomic update
      touched = touched + int(rankA1+rankA2+rankA3,8)*int(nA+nB,8)
      deallocate(UU, VV, U3A, V3A, U3, V3)
    end block

    call far_block_directional_uw(iBloc, nB, xB,yB,zB,oB,kB,fb8, iAloc, nA, xA,yA,zA,oA,kA,fa8, &
                                   Uacc, Wacc, tol, maxrank, touched)

    deallocate(fa8, fb8, xA,yA,zA,oA,kA,fA4,xB,yB,zB,oB,kB,fB4)
end subroutine far_block_full

subroutine far_block_directional(kern_id, iRloc, nR, xR,yR,zR,oR,kR,fR8, iCloc, nC, xC,yC,zC,oC,kC,fC8, &
                                  Racc, tol, maxrank, touched)
    implicit none
    integer,intent(in) :: kern_id, iRloc(nR), nR, iCloc(nC), nC, maxrank
    real(4),intent(in) :: xR(:),yR(:),zR(:),oR(:),kR(:), xC(:),yC(:),zC(:),oC(:),kC(:)
    real(8),intent(in) :: fR8(:), fC8(:), tol
    real(8),intent(inout) :: Racc(:)
    integer(8),intent(inout) :: touched
    real(8),allocatable :: UU(:,:), VV(:,:), VtfC(:)
    integer :: rank, k
    call aca_compress_k(kern_id, nR, xR,yR,zR,oR,kR, nC, xC,yC,zC,oC,kC, tol, maxrank, UU, VV, rank)
    !$omp atomic update
    dbg_rank_sum(kern_id) = dbg_rank_sum(kern_id) + int(rank,8)
    !$omp atomic update
    dbg_rank_cnt(kern_id) = dbg_rank_cnt(kern_id) + 1_8
    if (rank .gt. 0) then
       allocate(VtfC(rank))
       VtfC = matmul(transpose(VV), fC8(1:nC))
       do k = 1,nR
          !$omp atomic update
          Racc(iRloc(k)) = Racc(iRloc(k)) + sum(UU(k,1:rank)*VtfC(1:rank))
       enddo
       !$omp atomic update
       touched = touched + int(rank,8)*int(nR+nC,8)
       deallocate(VtfC)
    endif
    deallocate(UU, VV)
end subroutine far_block_directional

subroutine far_block_directional_uw(iRloc, nR, xR,yR,zR,oR,kR,fR8, iCloc, nC, xC,yC,zC,oC,kC,fC8, &
                                     Uacc, Wacc, tol, maxrank, touched)
    implicit none
    integer,intent(in) :: iRloc(nR), nR, iCloc(nC), nC, maxrank
    real(4),intent(in) :: xR(:),yR(:),zR(:),oR(:),kR(:), xC(:),yC(:),zC(:),oC(:),kC(:)
    real(8),intent(in) :: fR8(:), fC8(:), tol
    real(8),intent(inout) :: Uacc(:), Wacc(:)
    integer(8),intent(inout) :: touched
    real(8),allocatable :: U2(:,:), V2(:,:), U3(:,:), V3(:,:), VtfC2(:), VtfC3(:)
    integer :: rank2, rank3, k
    call aca_compress_uw(nR, xR,yR,zR,oR,kR, nC, xC,yC,zC,oC,kC, tol, maxrank, U2,V2,rank2, U3,V3,rank3)
    !$omp atomic update
    dbg_rank_sum(2) = dbg_rank_sum(2) + int(rank2,8)
    !$omp atomic update
    dbg_rank_cnt(2) = dbg_rank_cnt(2) + 1_8
    !$omp atomic update
    dbg_rank_sum(3) = dbg_rank_sum(3) + int(rank3,8)
    !$omp atomic update
    dbg_rank_cnt(3) = dbg_rank_cnt(3) + 1_8
    if (rank2 .gt. 0) then
       allocate(VtfC2(rank2))
       VtfC2 = matmul(transpose(V2), fC8(1:nC))
       do k = 1,nR
          !$omp atomic update
          Uacc(iRloc(k)) = Uacc(iRloc(k)) + sum(U2(k,1:rank2)*VtfC2(1:rank2))
       enddo
       deallocate(VtfC2)
    endif
    if (rank3 .gt. 0) then
       allocate(VtfC3(rank3))
       VtfC3 = matmul(transpose(V3), fC8(1:nC))
       do k = 1,nR
          !$omp atomic update
          Wacc(iRloc(k)) = Wacc(iRloc(k)) + sum(U3(k,1:rank3)*VtfC3(1:rank3))
       enddo
       deallocate(VtfC3)
    endif
    !$omp atomic update
    touched = touched + int(rank2,8)*int(nR+nC,8)
    deallocate(U2, V2, U3, V3)
end subroutine far_block_directional_uw

subroutine near_block(iAloc, nA, iBloc, nB, xs, ys, zs, os, ks, fs, Aacc, touched, is_self)
    implicit none
    integer,intent(in) :: iAloc(nA), iBloc(nB), nA, nB
    real(8),intent(in) :: xs(:), ys(:), zs(:), os(:), ks(:), fs(:)
    real(8),intent(inout) :: Aacc(:)
    integer(8),intent(inout) :: touched
    logical,intent(in) :: is_self
    integer :: k, m
    real(8) :: pd
    real(8),allocatable :: sRowA(:), sRowB(:)
    real(4),allocatable :: xA(:), yA(:), zA(:), oA(:), kA(:), fA4(:)
    real(4),allocatable :: xB(:), yB(:), zB(:), oB(:), kB(:), fB4(:)

    call gather_block(iAloc, nA, xs, ys, zs, os, ks, fs, xA, yA, zA, oA, kA, fA4)
    if (is_self) then
       allocate(sRowA(nA)); sRowA = 0.0d0
       do k = 1,nA
          pd = real(phi4(k, xA,yA,zA,oA,kA, k, xA,yA,zA,oA,kA),8)
          sRowA(k) = sRowA(k) + real(fA4(k),8)*pd
          do m = k+1,nA
             pd = real(phi4(k, xA,yA,zA,oA,kA, m, xA,yA,zA,oA,kA),8)
             sRowA(k) = sRowA(k) + real(fA4(m),8)*pd
             sRowA(m) = sRowA(m) + real(fA4(k),8)*pd
          enddo
       enddo
       do k = 1,nA
          !$omp atomic update
          Aacc(iAloc(k)) = Aacc(iAloc(k)) + sRowA(k)
       enddo
       !$omp atomic update
       touched = touched + (int(nA,8)*int(nA+1,8))/2_8
       deallocate(xA,yA,zA,oA,kA,fA4,sRowA)
    else
       call gather_block(iBloc, nB, xs, ys, zs, os, ks, fs, xB, yB, zB, oB, kB, fB4)
       allocate(sRowA(nA), sRowB(nB)); sRowA = 0.0d0; sRowB = 0.0d0
       do k = 1,nA
          do m = 1,nB
             pd = real(phi4(k, xA,yA,zA,oA,kA, m, xB,yB,zB,oB,kB),8)
             sRowA(k) = sRowA(k) + real(fB4(m),8)*pd
             sRowB(m) = sRowB(m) + real(fA4(k),8)*pd
          enddo
       enddo
       do k = 1,nA
          !$omp atomic update
          Aacc(iAloc(k)) = Aacc(iAloc(k)) + sRowA(k)
       enddo
       do k = 1,nB
          !$omp atomic update
          Aacc(iBloc(k)) = Aacc(iBloc(k)) + sRowB(k)
       enddo
       !$omp atomic update
       touched = touched + int(nA,8)*int(nB,8)
       deallocate(xA,yA,zA,oA,kA,fA4,xB,yB,zB,oB,kB,fB4,sRowA,sRowB)
    endif
end subroutine near_block

subroutine near_block_full(iAloc, nA, iBloc, nB, xs, ys, zs, os, ks, fs, Aacc, Uacc, Wacc, touched, is_self)
    implicit none
    integer,intent(in) :: iAloc(nA), iBloc(nB), nA, nB
    real(8),intent(in) :: xs(:), ys(:), zs(:), os(:), ks(:), fs(:)
    real(8),intent(inout) :: Aacc(:), Uacc(:), Wacc(:)
    integer(8),intent(inout) :: touched
    logical,intent(in) :: is_self
    integer :: k, m
    real(8) :: phid_k, phid_m, uwfac_k, uwfac_m
    real(4) :: dx,dy,dz,R2,gi,gj,phi
    real(4),allocatable :: xA(:), yA(:), zA(:), oA(:), kA(:), fA4(:)
    real(4),allocatable :: xB(:), yB(:), zB(:), oB(:), kB(:), fB4(:)
    real(8),allocatable :: sRowA(:), sRowU(:), sRowW(:)
    real(8),allocatable :: sRowAb(:), sRowUb(:), sRowWb(:)

    call gather_block(iAloc, nA, xs, ys, zs, os, ks, fs, xA, yA, zA, oA, kA, fA4)
    if (is_self) then
       allocate(sRowA(nA), sRowU(nA), sRowW(nA))
       sRowA = 0.0d0; sRowU = 0.0d0; sRowW = 0.0d0
       do k = 1,nA
          gi = kA(k); gj = kA(k)
          phi = -1.5_4/(gi*gj*(gi+gj))
          phid_k = real(fA4(k),8)*real(phi,8)
          uwfac_k = real(1.0_4/gi + 1.0_4/(gi+gj),8)
          sRowA(k) = sRowA(k) + phid_k
          sRowU(k) = sRowU(k) + phid_k*uwfac_k
          do m = k+1,nA
             dx=xA(k)-xA(m); dy=yA(k)-yA(m); dz=zA(k)-zA(m)
             R2=dx*dx+dy*dy+dz*dz
             gi=oA(k)*R2+kA(k); gj=oA(m)*R2+kA(m)
             phi = -1.5_4/(gi*gj*(gi+gj))
             phid_k = real(fA4(m),8)*real(phi,8)
             phid_m = real(fA4(k),8)*real(phi,8)
             uwfac_k = real(1.0_4/gi + 1.0_4/(gi+gj),8)
             uwfac_m = real(1.0_4/gj + 1.0_4/(gi+gj),8)
             sRowA(k) = sRowA(k) + phid_k
             sRowA(m) = sRowA(m) + phid_m
             sRowU(k) = sRowU(k) + phid_k*uwfac_k
             sRowU(m) = sRowU(m) + phid_m*uwfac_m
             sRowW(k) = sRowW(k) + phid_k*uwfac_k*real(R2,8)
             sRowW(m) = sRowW(m) + phid_m*uwfac_m*real(R2,8)
          enddo
       enddo
       do k = 1,nA
          !$omp atomic update
          Aacc(iAloc(k)) = Aacc(iAloc(k)) + sRowA(k)
          !$omp atomic update
          Uacc(iAloc(k)) = Uacc(iAloc(k)) + sRowU(k)
          !$omp atomic update
          Wacc(iAloc(k)) = Wacc(iAloc(k)) + sRowW(k)
       enddo
       !$omp atomic update
       touched = touched + (int(nA,8)*int(nA+1,8))/2_8
       deallocate(xA,yA,zA,oA,kA,fA4,sRowA,sRowU,sRowW)
    else
       call gather_block(iBloc, nB, xs, ys, zs, os, ks, fs, xB, yB, zB, oB, kB, fB4)
       allocate(sRowA(nA), sRowU(nA), sRowW(nA), sRowAb(nB), sRowUb(nB), sRowWb(nB))
       sRowA = 0.0d0; sRowU = 0.0d0; sRowW = 0.0d0
       sRowAb = 0.0d0; sRowUb = 0.0d0; sRowWb = 0.0d0
       do k = 1,nA
          do m = 1,nB
             dx=xA(k)-xB(m); dy=yA(k)-yB(m); dz=zA(k)-zB(m)
             R2=dx*dx+dy*dy+dz*dz
             gi=oA(k)*R2+kA(k); gj=oB(m)*R2+kB(m)
             phi = -1.5_4/(gi*gj*(gi+gj))
             phid_k = real(fB4(m),8)*real(phi,8)
             phid_m = real(fA4(k),8)*real(phi,8)
             uwfac_k = real(1.0_4/gi + 1.0_4/(gi+gj),8)
             uwfac_m = real(1.0_4/gj + 1.0_4/(gi+gj),8)
             sRowA(k) = sRowA(k) + phid_k
             sRowAb(m) = sRowAb(m) + phid_m
             sRowU(k) = sRowU(k) + phid_k*uwfac_k
             sRowUb(m) = sRowUb(m) + phid_m*uwfac_m
             sRowW(k) = sRowW(k) + phid_k*uwfac_k*real(R2,8)
             sRowWb(m) = sRowWb(m) + phid_m*uwfac_m*real(R2,8)
          enddo
       enddo
       do k = 1,nA
          !$omp atomic update
          Aacc(iAloc(k)) = Aacc(iAloc(k)) + sRowA(k)
          !$omp atomic update
          Uacc(iAloc(k)) = Uacc(iAloc(k)) + sRowU(k)
          !$omp atomic update
          Wacc(iAloc(k)) = Wacc(iAloc(k)) + sRowW(k)
       enddo
       do k = 1,nB
          !$omp atomic update
          Aacc(iBloc(k)) = Aacc(iBloc(k)) + sRowAb(k)
          !$omp atomic update
          Uacc(iBloc(k)) = Uacc(iBloc(k)) + sRowUb(k)
          !$omp atomic update
          Wacc(iBloc(k)) = Wacc(iBloc(k)) + sRowWb(k)
       enddo
       !$omp atomic update
       touched = touched + int(nA,8)*int(nB,8)
       deallocate(xA,yA,zA,oA,kA,fA4,xB,yB,zB,oB,kB,fB4,sRowA,sRowU,sRowW,sRowAb,sRowUb,sRowWb)
    endif
end subroutine near_block_full

recursive subroutine dual_cross(A, B, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, &
                                 n_far, n_near, touched)
    implicit none
    type(vv10_node),intent(in) :: A, B
    integer,intent(in) :: perm(:)
    real(8),intent(in) :: xs(:), ys(:), zs(:), os(:), ks(:), fs(:), eta, tol, rcut
    integer,intent(in) :: maxrank
    real(8),intent(inout) :: Aacc(:)
    integer,intent(inout) :: n_far, n_near
    integer(8),intent(inout) :: touched

    logical :: big

    if (beyond_rcut(A, B, rcut)) return
    if (admissible(A, B, eta)) then
       call far_block(perm, perm(A%lo:A%hi), A%hi-A%lo+1, perm(B%lo:B%hi), B%hi-B%lo+1, &
                       xs, ys, zs, os, ks, fs, Aacc, tol, maxrank, touched)
       !$omp atomic update
       n_far = n_far+1
       return
    endif
    if (A%is_leaf .and. B%is_leaf) then
       call near_block(perm(A%lo:A%hi), A%hi-A%lo+1, perm(B%lo:B%hi), B%hi-B%lo+1, &
                        xs, ys, zs, os, ks, fs, Aacc, touched, .false.)
       !$omp atomic update
       n_near = n_near+1
       return
    endif
    big = int(A%hi-A%lo+1,8)*int(B%hi-B%lo+1,8) .gt. TASK_CUTOFF
    if (A%is_leaf) then
       !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
       call dual_cross(A, B%left,  perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
       !$omp end task
       !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
       call dual_cross(A, B%right, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
       !$omp end task
       !$omp taskwait
    else if (B%is_leaf) then
       !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
       call dual_cross(A%left,  B, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
       !$omp end task
       !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
       call dual_cross(A%right, B, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
       !$omp end task
       !$omp taskwait
    else if (A%diam .ge. B%diam) then
       !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
       call dual_cross(A%left,  B, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
       !$omp end task
       !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
       call dual_cross(A%right, B, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
       !$omp end task
       !$omp taskwait
    else
       !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
       call dual_cross(A, B%left,  perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
       !$omp end task
       !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
       call dual_cross(A, B%right, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
       !$omp end task
       !$omp taskwait
    endif
end subroutine dual_cross

recursive subroutine dual_self(A, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, &
                                n_far, n_near, n_self, touched)
    implicit none
    type(vv10_node),intent(in) :: A
    integer,intent(in) :: perm(:)
    real(8),intent(in) :: xs(:), ys(:), zs(:), os(:), ks(:), fs(:), eta, tol, rcut
    integer,intent(in) :: maxrank
    real(8),intent(inout) :: Aacc(:)
    integer,intent(inout) :: n_far, n_near, n_self
    integer(8),intent(inout) :: touched

    logical :: big

    if (A%is_leaf) then
       call near_block(perm(A%lo:A%hi), A%hi-A%lo+1, perm(A%lo:A%hi), A%hi-A%lo+1, &
                        xs, ys, zs, os, ks, fs, Aacc, touched, .true.)
       !$omp atomic update
       n_self = n_self+1
       return
    endif
    big = int(A%hi-A%lo+1,8)*int(A%hi-A%lo+1,8) .gt. TASK_CUTOFF
    !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,n_self,touched) firstprivate(eta,tol,maxrank,rcut)
    call dual_self(A%left,  perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, n_self, touched)
    !$omp end task
    !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,n_self,touched) firstprivate(eta,tol,maxrank,rcut)
    call dual_self(A%right, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, n_self, touched)
    !$omp end task
    !$omp task if(big) shared(perm,xs,ys,zs,os,ks,fs,Aacc,n_far,n_near,touched) firstprivate(eta,tol,maxrank,rcut)
    call dual_cross(A%left, A%right, perm, xs, ys, zs, os, ks, fs, Aacc, eta, tol, maxrank, rcut, n_far, n_near, touched)
    !$omp end task
    !$omp taskwait
end subroutine dual_self

subroutine push_near(plan, loA, hiA, loB, hiB, is_self)
    implicit none
    type(vv10_plan),intent(inout) :: plan
    integer,intent(in) :: loA, hiA, loB, hiB
    logical,intent(in) :: is_self
    integer,allocatable :: tmpI(:)
    logical,allocatable :: tmpL(:)
    integer :: cap, newcap
    cap = size(plan%near_loA)
    if (plan%n_near .eq. cap) then
       newcap = cap*2
       allocate(tmpI(newcap)); tmpI(1:cap) = plan%near_loA; call move_alloc(tmpI, plan%near_loA)
       allocate(tmpI(newcap)); tmpI(1:cap) = plan%near_hiA; call move_alloc(tmpI, plan%near_hiA)
       allocate(tmpI(newcap)); tmpI(1:cap) = plan%near_loB; call move_alloc(tmpI, plan%near_loB)
       allocate(tmpI(newcap)); tmpI(1:cap) = plan%near_hiB; call move_alloc(tmpI, plan%near_hiB)
       allocate(tmpL(newcap)); tmpL(1:cap) = plan%near_is_self; call move_alloc(tmpL, plan%near_is_self)
    endif
    plan%n_near = plan%n_near + 1
    plan%near_loA(plan%n_near) = loA; plan%near_hiA(plan%n_near) = hiA
    plan%near_loB(plan%n_near) = loB; plan%near_hiB(plan%n_near) = hiB
    plan%near_is_self(plan%n_near) = is_self
end subroutine push_near

subroutine push_far(plan, loA, hiA, loB, hiB)
    implicit none
    type(vv10_plan),intent(inout) :: plan
    integer,intent(in) :: loA, hiA, loB, hiB
    integer,allocatable :: tmpI(:)
    integer :: cap, newcap
    cap = size(plan%far_loA)
    if (plan%n_far .eq. cap) then
       newcap = cap*2
       allocate(tmpI(newcap)); tmpI(1:cap) = plan%far_loA; call move_alloc(tmpI, plan%far_loA)
       allocate(tmpI(newcap)); tmpI(1:cap) = plan%far_hiA; call move_alloc(tmpI, plan%far_hiA)
       allocate(tmpI(newcap)); tmpI(1:cap) = plan%far_loB; call move_alloc(tmpI, plan%far_loB)
       allocate(tmpI(newcap)); tmpI(1:cap) = plan%far_hiB; call move_alloc(tmpI, plan%far_hiB)
    endif
    plan%n_far = plan%n_far + 1
    plan%far_loA(plan%n_far) = loA; plan%far_hiA(plan%n_far) = hiA
    plan%far_loB(plan%n_far) = loB; plan%far_hiB(plan%n_far) = hiB
end subroutine push_far

recursive subroutine record_cross(A, B, eta, plan)
    implicit none
    type(vv10_node),intent(in) :: A, B
    real(8),intent(in) :: eta
    type(vv10_plan),intent(inout) :: plan
    if (admissible(A, B, eta)) then
       call push_far(plan, A%lo, A%hi, B%lo, B%hi)
       return
    endif
    if (A%is_leaf .and. B%is_leaf) then
       call push_near(plan, A%lo, A%hi, B%lo, B%hi, .false.)
       return
    endif
    if (A%is_leaf) then
       call record_cross(A, B%left,  eta, plan)
       call record_cross(A, B%right, eta, plan)
    else if (B%is_leaf) then
       call record_cross(A%left,  B, eta, plan)
       call record_cross(A%right, B, eta, plan)
    else if (A%diam .ge. B%diam) then
       call record_cross(A%left,  B, eta, plan)
       call record_cross(A%right, B, eta, plan)
    else
       call record_cross(A, B%left,  eta, plan)
       call record_cross(A, B%right, eta, plan)
    endif
end subroutine record_cross

recursive subroutine record_self(A, eta, plan)
    implicit none
    type(vv10_node),intent(in) :: A
    real(8),intent(in) :: eta
    type(vv10_plan),intent(inout) :: plan
    if (A%is_leaf) then
       call push_near(plan, A%lo, A%hi, A%lo, A%hi, .true.)
       return
    endif
    call record_self(A%left,  eta, plan)
    call record_self(A%right, eta, plan)
    call record_cross(A%left, A%right, eta, plan)
end subroutine record_self

subroutine vv10_aca_plan_build(np, coor, act, n_active, plan, leaf_size_in, eta_in)
    implicit none
    integer,intent(in) :: np, n_active, act(np)
    real(8),intent(in) :: coor(3,np)
    type(vv10_plan),intent(out) :: plan
    integer,intent(in),optional :: leaf_size_in
    real(8),intent(in),optional :: eta_in
    real(8),allocatable :: xs(:), ys(:), zs(:)
    type(vv10_node),pointer :: root
    integer :: leaf_size, i, ia
    real(8) :: eta

    leaf_size = LEAF_SIZE_DEFAULT; if (present(leaf_size_in)) leaf_size = leaf_size_in
    eta       = ETA_DEFAULT;       if (present(eta_in))       eta       = eta_in

    plan%n = n_active
    allocate(plan%act(n_active))
    allocate(plan%perm(n_active))
    allocate(xs(n_active), ys(n_active), zs(n_active))
    do i = 1,n_active
       ia = act(i)
       plan%act(i) = ia
       xs(i) = coor(1,ia); ys(i) = coor(2,ia); zs(i) = coor(3,ia)
       plan%perm(i) = i
    enddo

    nullify(root)
    call build_node(root, plan%perm, 1, n_active, xs, ys, zs, leaf_size)

    allocate(plan%near_loA(16), plan%near_hiA(16), plan%near_loB(16), plan%near_hiB(16), plan%near_is_self(16))
    allocate(plan%far_loA(16), plan%far_hiA(16), plan%far_loB(16), plan%far_hiB(16))
    plan%n_near = 0; plan%n_far = 0

    call record_self(root, eta, plan)

    call free_tree(root)
    deallocate(xs, ys, zs)
end subroutine vv10_aca_plan_build

subroutine vv10_aca_plan_evaluate(plan, np, rho, weight, coor, omega0, kappa, beta, Enl, A_out, &
                                   tol_in, maxrank_in, kernel_touched)
    implicit none
    type(vv10_plan),intent(in) :: plan
    integer,intent(in) :: np
    real(8),intent(in) :: rho(np), weight(np), coor(3,np), omega0(np), kappa(np), beta
    real(8),intent(out) :: Enl, A_out(np)
    real(8),intent(in),optional :: tol_in
    integer,intent(in),optional :: maxrank_in
    integer(8),intent(out),optional :: kernel_touched

    real(8),allocatable :: xs(:), ys(:), zs(:), fs(:), os(:), ks(:), Aacc(:)
    integer :: i, ia, k, maxrank
    integer(8) :: touched
    real(8) :: tol

    tol     = ACA_TOL_DEFAULT;     if (present(tol_in))     tol     = tol_in
    maxrank = ACA_MAXRANK_DEFAULT; if (present(maxrank_in)) maxrank = maxrank_in

    allocate(xs(plan%n), ys(plan%n), zs(plan%n), fs(plan%n), os(plan%n), ks(plan%n), Aacc(plan%n))
    do i = 1,plan%n
       ia = plan%act(i)
       xs(i) = coor(1,ia); ys(i) = coor(2,ia); zs(i) = coor(3,ia)
       fs(i) = weight(ia)*rho(ia)
       os(i) = omega0(ia); ks(i) = kappa(ia)
    enddo
    Aacc = 0.0d0
    touched = 0_8

    !$omp parallel do schedule(dynamic) reduction(+:Aacc,touched)
    do k = 1,plan%n_near
       call near_block(plan%perm(plan%near_loA(k):plan%near_hiA(k)), plan%near_hiA(k)-plan%near_loA(k)+1, &
                        plan%perm(plan%near_loB(k):plan%near_hiB(k)), plan%near_hiB(k)-plan%near_loB(k)+1, &
                        xs, ys, zs, os, ks, fs, Aacc, touched, plan%near_is_self(k))
    enddo
    !$omp end parallel do

    !$omp parallel do schedule(dynamic) reduction(+:Aacc,touched)
    do k = 1,plan%n_far
       call far_block(plan%perm, plan%perm(plan%far_loA(k):plan%far_hiA(k)), plan%far_hiA(k)-plan%far_loA(k)+1, &
                       plan%perm(plan%far_loB(k):plan%far_hiB(k)), plan%far_hiB(k)-plan%far_loB(k)+1, &
                       xs, ys, zs, os, ks, fs, Aacc, tol, maxrank, touched)
    enddo
    !$omp end parallel do

    Enl = beta*sum(fs)
    A_out = 0.0d0
    do i = 1,plan%n
       A_out(plan%act(i)) = Aacc(i)
       Enl = Enl + 0.5d0*fs(i)*Aacc(i)
    enddo
    if (present(kernel_touched)) kernel_touched = touched
    deallocate(xs, ys, zs, fs, os, ks, Aacc)
end subroutine vv10_aca_plan_evaluate

subroutine vv10_aca_plan_evaluate_full(plan, np, rho, gamma, weight, coor, omega0, kappa, beta, vv10_C, &
                                        Enl, Fn, Fg, tol_in, maxrank_in, kernel_touched)
    implicit none
    type(vv10_plan),intent(in) :: plan
    integer,intent(in) :: np
    real(8),intent(in) :: rho(np), gamma(np), weight(np), coor(3,np), omega0(np), kappa(np)
    real(8),intent(in) :: beta, vv10_C
    real(8),intent(out) :: Enl, Fn(np), Fg(np)
    real(8),intent(in),optional :: tol_in
    integer,intent(in),optional :: maxrank_in
    integer(8),intent(out),optional :: kernel_touched

    real(8),allocatable :: xs(:), ys(:), zs(:), fs(:), os(:), ks(:)
    real(8),allocatable :: Aacc(:), Uacc(:), Wacc(:)
    real(8),allocatable :: domega_dn(:), domega_dg(:), dkappa_dn(:)
    integer :: i, ia, k, maxrank
    integer(8) :: touched
    real(8) :: tol, pi

    tol     = ACA_TOL_DEFAULT;     if (present(tol_in))     tol     = tol_in
    maxrank = ACA_MAXRANK_DEFAULT; if (present(maxrank_in)) maxrank = maxrank_in
    pi = 4.0d0*atan(1.0d0)

    allocate(xs(plan%n), ys(plan%n), zs(plan%n), fs(plan%n), os(plan%n), ks(plan%n))
    allocate(Aacc(plan%n), Uacc(plan%n), Wacc(plan%n))
    allocate(domega_dn(plan%n), domega_dg(plan%n), dkappa_dn(plan%n))
    do i = 1,plan%n
       ia = plan%act(i)
       xs(i) = coor(1,ia); ys(i) = coor(2,ia); zs(i) = coor(3,ia)
       fs(i) = weight(ia)*rho(ia)
       os(i) = omega0(ia); ks(i) = kappa(ia)
       domega_dn(i) = (4.0d0*pi/3.0d0 - 4.0d0*vv10_C*gamma(ia)*gamma(ia)/rho(ia)**5)/(2.0d0*omega0(ia))
       domega_dg(i) = vv10_C*gamma(ia)/(omega0(ia)*rho(ia)**4)
       dkappa_dn(i) = kappa(ia)/(6.0d0*rho(ia))
    enddo
    Aacc = 0.0d0; Uacc = 0.0d0; Wacc = 0.0d0
    touched = 0_8

    block
      integer(8) :: cn0, cn1, cn2, cnrate
      call system_clock(cn0, cnrate)
    !$omp parallel do schedule(dynamic) reduction(+:Aacc,Uacc,Wacc,touched)
    do k = 1,plan%n_near
       call near_block_full(plan%perm(plan%near_loA(k):plan%near_hiA(k)), plan%near_hiA(k)-plan%near_loA(k)+1, &
                             plan%perm(plan%near_loB(k):plan%near_hiB(k)), plan%near_hiB(k)-plan%near_loB(k)+1, &
                             xs, ys, zs, os, ks, fs, Aacc, Uacc, Wacc, touched, plan%near_is_self(k))
    enddo
    !$omp end parallel do
      call system_clock(cn1)
      print '("  [near_block_full: ",F8.2,"s]")', real(cn1-cn0,8)/real(cnrate,8)
    end block

    block
      integer(8) :: cf0, cf1, cfrate
      call system_clock(cf0, cfrate)
    !$omp parallel do schedule(dynamic) reduction(+:Aacc,Uacc,Wacc,touched)
    do k = 1,plan%n_far
       call far_block_full(plan%perm(plan%far_loA(k):plan%far_hiA(k)), plan%far_hiA(k)-plan%far_loA(k)+1, &
                            plan%perm(plan%far_loB(k):plan%far_hiB(k)), plan%far_hiB(k)-plan%far_loB(k)+1, &
                            xs, ys, zs, os, ks, fs, Aacc, Uacc, Wacc, tol, maxrank, touched)
    enddo
    !$omp end parallel do
      call system_clock(cf1)
      print '("  [far_block_full: ",F8.2,"s]")', real(cf1-cf0,8)/real(cfrate,8)
    end block

    Enl = beta*sum(fs)
    Fn = 0.0d0; Fg = 0.0d0
    do i = 1,plan%n
       ia = plan%act(i)
       Enl = Enl + 0.5d0*fs(i)*Aacc(i)
       Fn(ia) = beta + Aacc(i) + rho(ia)*(dkappa_dn(i)*(-Uacc(i)) + domega_dn(i)*(-Wacc(i)))
       Fg(ia) = rho(ia)*domega_dg(i)*(-Wacc(i))
    enddo
    if (present(kernel_touched)) kernel_touched = touched
    deallocate(xs, ys, zs, fs, os, ks, Aacc, Uacc, Wacc, domega_dn, domega_dg, dkappa_dn)
end subroutine vv10_aca_plan_evaluate_full

subroutine vv10_aca_plan_free(plan)
    implicit none
    type(vv10_plan),intent(inout) :: plan
    if (allocated(plan%act)) deallocate(plan%act)
    if (allocated(plan%perm)) deallocate(plan%perm)
    if (allocated(plan%near_loA)) deallocate(plan%near_loA, plan%near_hiA, plan%near_loB, plan%near_hiB, plan%near_is_self)
    if (allocated(plan%far_loA)) deallocate(plan%far_loA, plan%far_hiA, plan%far_loB, plan%far_hiB)
    plan%n = 0; plan%n_near = 0; plan%n_far = 0
end subroutine vv10_aca_plan_free

end module mod_vv10_aca
