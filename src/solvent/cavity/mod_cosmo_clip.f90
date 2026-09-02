! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_clip: GePol sequential multi-sphere clipping of one

module mod_cosmo_clip
use mod_cosmo_constants, only: MAXV, COSMO_EDGE_TOL, COSMO_INTER_TOL, COSMO_INTER_MAXIT
use mod_cosmo_dual, only: dscal, dvec3, cosmo_cross
implicit none

contains

function cosmo_two_sphere_circle_center(ci, Ri, cj, Rj) result(cc)
real(8),intent(in) :: ci(3), Ri, cj(3), Rj
real(8) :: cc(3), de2
de2 = sum((cj-ci)**2)
cc = ci + (cj-ci)*(Ri*Ri - Rj*Rj + de2)/(2.0d0*de2)
end function cosmo_two_sphere_circle_center

subroutine cosmo_two_sphere_circle_center_deriv(ci, Ri, cj, Rj, cc, dcc_dci, dcc_dcj)
real(8),intent(in) :: ci(3), Ri, cj(3), Rj
real(8),intent(out) :: cc(3), dcc_dci(3,3), dcc_dcj(3,3)
real(8) :: d(3), s, t, k
integer :: a,b
d = cj - ci
s = sum(d*d)
t = (Ri*Ri - Rj*Rj + s)/(2.0d0*s)
cc = ci + d*t
k = (Ri*Ri - Rj*Rj)/(s*s)
dcc_dci = 0.0d0
dcc_dcj = 0.0d0
do a = 1,3
   dcc_dci(a,a) = (1.0d0-t)
   dcc_dcj(a,a) = t
   do b = 1,3
      dcc_dci(a,b) = dcc_dci(a,b) + k*d(a)*d(b)
      dcc_dcj(a,b) = dcc_dcj(a,b) - k*d(a)*d(b)
   enddo
enddo
end subroutine cosmo_two_sphere_circle_center_deriv

subroutine cosmo_solve3x3(J, rhs, x)
real(8),intent(in) :: J(3,3), rhs(3)
real(8),intent(out) :: x(3)
real(8) :: detJ, Jc(3,3)
detJ = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) &
     - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
     + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
Jc = J; Jc(:,1) = rhs
x(1) = cosmo_det3(Jc)/detJ
Jc = J; Jc(:,2) = rhs
x(2) = cosmo_det3(Jc)/detJ
Jc = J; Jc(:,3) = rhs
x(3) = cosmo_det3(Jc)/detJ
end subroutine cosmo_solve3x3

function cosmo_det3(J) result(d)
real(8),intent(in) :: J(3,3)
real(8) :: d
d = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) &
  - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
  + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
end function cosmo_det3

subroutine cosmo_vertex_jac(v, Ci, Ri, has_jprev, Cjprev, Rjprev, OVa, OVb, Cnew, Rnew, &
                             dv_dCi, dv_dCjprev, dv_dCnew)
real(8),intent(in) :: v(3), Ci(3), Ri
logical,intent(in) :: has_jprev
real(8),intent(in) :: Cjprev(3), Rjprev
real(8),intent(in) :: OVa(3), OVb(3)
real(8),intent(in) :: Cnew(3), Rnew
real(8),intent(out) :: dv_dCi(3,3), dv_dCjprev(3,3), dv_dCnew(3,3)

real(8) :: row0dv(3), row1dv(3), row2dv(3)
real(8) :: row0dCi(3), row1dCi(3), row2dCi(3), row1dCjprev(3), row2dCnew(3)
real(8) :: Jmat(3,3), rhs(3), n_edge(3)
integer :: c

row0dv = 2.0d0*(v-Ci)
row0dCi = -2.0d0*(v-Ci)

if (has_jprev) then
   row1dv = 2.0d0*(Cjprev-Ci)
   row1dCi = -2.0d0*(v-Ci)
   row1dCjprev = 2.0d0*(v-Cjprev)
else
   call cosmo_cross(OVa-Ci, OVb-Ci, n_edge)
   row1dv = n_edge
   row1dCi = -n_edge
   row1dCjprev = 0.0d0
endif

row2dv = 2.0d0*(Cnew-Ci)
row2dCi = -2.0d0*(v-Ci)
row2dCnew = 2.0d0*(v-Cnew)

Jmat(1,:) = row0dv
Jmat(2,:) = row1dv
Jmat(3,:) = row2dv

do c = 1,3
   rhs = (/ -row0dCi(c), -row1dCi(c), -row2dCi(c) /)
   call cosmo_solve3x3(Jmat, rhs, dv_dCi(:,c))
enddo

dv_dCjprev = 0.0d0
if (has_jprev) then
   do c = 1,3
      rhs = (/ 0.0d0, -row1dCjprev(c), 0.0d0 /)
      call cosmo_solve3x3(Jmat, rhs, dv_dCjprev(:,c))
   enddo
endif

do c = 1,3
   rhs = (/ 0.0d0, 0.0d0, -row2dCnew(c) /)
   call cosmo_solve3x3(Jmat, rhs, dv_dCnew(:,c))
enddo
end subroutine cosmo_vertex_jac

function cosmo_inter(P1, P2, P3, R, cj, Rj, intcas) result(P4)
real(8),intent(in) :: P1(3), P2(3), P3(3), R, cj(3), Rj
integer,intent(in) :: intcas
real(8) :: P4(3)
real(8) :: alphat, delta, diff, dnorm
integer :: m

if (abs(sqrt(sum((P1-cj)**2)) - Rj) < 1.0d-10 .and. &
    abs(sqrt(sum((P2-cj)**2)) - Rj) >= 1.0d-10) then
   P4 = P1
   return
else if (abs(sqrt(sum((P2-cj)**2)) - Rj) < 1.0d-10 .and. &
         abs(sqrt(sum((P1-cj)**2)) - Rj) >= 1.0d-10) then
   P4 = P2
   return
endif

alphat = 0.5d0
delta = 0.0d0
m = 1
do
   alphat = alphat + delta
   P4 = P1 + alphat*(P2-P1) - P3
   dnorm = sqrt(sum(P4*P4))
   P4 = P4*R/dnorm + P3
   diff = sqrt(sum((P4-cj)**2)) - Rj
   if (abs(diff) < COSMO_INTER_TOL) return
   if (intcas == 0) then
      if (diff > 0.0d0) then
         delta = 1.0d0/(2.0d0**(m+1))
      else
         delta = -1.0d0/(2.0d0**(m+1))
      endif
   else
      if (diff > 0.0d0) then
         delta = -1.0d0/(2.0d0**(m+1))
      else
         delta = 1.0d0/(2.0d0**(m+1))
      endif
   endif
   m = m+1
   if (m > COSMO_INTER_MAXIT) return
enddo
end function cosmo_inter

subroutine cosmo_clip_tessera(tri_verts, i_atom, natoms, Rvdw, coor_bohr, &
                               verts, ccc, intsph, nv, discarded)
real(8),intent(in) :: tri_verts(3,3)
integer,intent(in) :: i_atom, natoms
real(8),intent(in) :: Rvdw(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(out) :: verts(3,MAXV), ccc(3,MAXV)
integer,intent(out) :: intsph(MAXV)
integer,intent(out) :: nv
logical,intent(out) :: discarded

real(8) :: verts_new(3,MAXV), ccc_new(3,MAXV)
integer :: intsph_new(MAXV)
integer :: nv_new
integer :: ind(MAXV), ltyp(MAXV)
real(8) :: pointl(3,MAXV)
integer :: ja, l, iv1, iv2, n, ncov, ncut, ii
real(8) :: p1(3), p2(3), p3(3), p4(3), point(3), rc
logical :: hit

discarded = .false.
nv = 3
verts(:,1:3) = tri_verts
do l = 1,3
   ccc(:,l) = coor_bohr(:,i_atom)
   intsph(l) = i_atom
enddo

do ja = 1,natoms
   if (ja == i_atom) cycle

   ncov = 0
   do n = 1,nv
      ind(n) = 0
      if (sqrt(sum((verts(:,n)-coor_bohr(:,ja))**2)) - Rvdw(ja) < 0.0d0) then
         ind(n) = 1
         ncov = ncov+1
      endif
   enddo
   if (ncov == nv) then
      discarded = .true.
      return
   endif
   if (ncov == 0) cycle

   do l = 1,nv
      iv1 = l
      iv2 = l+1; if (l == nv) iv2 = 1
      if (ind(iv1) == 1 .and. ind(iv2) == 1) then
         ltyp(l) = 0
      else if (ind(iv1) == 0 .and. ind(iv2) == 1) then
         ltyp(l) = 1
      else if (ind(iv1) == 1 .and. ind(iv2) == 0) then
         ltyp(l) = 2
      else
         ltyp(l) = 4
         rc = sqrt(sum((ccc(:,l)-verts(:,iv1))**2))
         hit = .false.
         do ii = 1,11
            point = verts(:,iv1) + real(ii,8)*(verts(:,iv2)-verts(:,iv1))/11.0d0
            point = point - ccc(:,l)
            point = point*rc/sqrt(sum(point*point)) + ccc(:,l)
            if (sqrt(sum((point-coor_bohr(:,ja))**2)) - Rvdw(ja) < -1.0d-10) then
               ltyp(l) = 3
               pointl(:,l) = point
               hit = .true.
               exit
            endif
         enddo
      endif
   enddo

   ncut = 0
   do l = 1,nv
      if (ltyp(l) == 1 .or. ltyp(l) == 2) ncut = ncut+1
      if (ltyp(l) == 3) ncut = ncut+2
   enddo
   if (ncut/2 > 1) then
      discarded = .true.
      return
   endif

   n = 0
   do l = 1,nv
      iv1 = l
      iv2 = l+1; if (l == nv) iv2 = 1
      select case (ltyp(l))
      case (0)
      case (1)
         n = n+1
         verts_new(:,n) = verts(:,iv1); ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
         p1 = verts(:,iv1); p2 = verts(:,iv2); p3 = ccc(:,iv1)
         p4 = cosmo_inter(p1, p2, p3, sqrt(sum((p1-p3)**2)), coor_bohr(:,ja), Rvdw(ja), 0)
         n = n+1
         verts_new(:,n) = p4
         ccc_new(:,n) = cosmo_two_sphere_circle_center(coor_bohr(:,i_atom), Rvdw(i_atom), &
                                                        coor_bohr(:,ja), Rvdw(ja))
         intsph_new(n) = ja
      case (2)
         p1 = verts(:,iv1); p2 = verts(:,iv2); p3 = ccc(:,iv1)
         p4 = cosmo_inter(p1, p2, p3, sqrt(sum((p1-p3)**2)), coor_bohr(:,ja), Rvdw(ja), 1)
         n = n+1
         verts_new(:,n) = p4; ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
      case (3)
         n = n+1
         verts_new(:,n) = verts(:,iv1); ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
         p1 = verts(:,iv1); p2 = pointl(:,l); p3 = ccc(:,iv1)
         p4 = cosmo_inter(p1, p2, p3, sqrt(sum((p1-p3)**2)), coor_bohr(:,ja), Rvdw(ja), 0)
         n = n+1
         verts_new(:,n) = p4
         ccc_new(:,n) = cosmo_two_sphere_circle_center(coor_bohr(:,i_atom), Rvdw(i_atom), &
                                                        coor_bohr(:,ja), Rvdw(ja))
         intsph_new(n) = ja
         p1 = pointl(:,l); p2 = verts(:,iv2); p3 = ccc(:,iv1)
         p4 = cosmo_inter(p1, p2, p3, sqrt(sum((p1-p3)**2)), coor_bohr(:,ja), Rvdw(ja), 1)
         n = n+1
         verts_new(:,n) = p4; ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
      case (4)
         n = n+1
         verts_new(:,n) = verts(:,iv1); ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
      end select
      if (n > MAXV-3) then
         discarded = .true.
         return
      endif
   enddo

   if (n < 3) then
      discarded = .true.
      return
   endif
   nv_new = n
   verts(:,1:nv_new) = verts_new(:,1:nv_new)
   ccc(:,1:nv_new) = ccc_new(:,1:nv_new)
   intsph(1:nv_new) = intsph_new(1:nv_new)
   nv = nv_new
enddo

block
   real(8) :: verts_clean(3,MAXV), ccc_clean(3,MAXV)
   integer :: intsph_clean(MAXV), nv_clean, i2
   nv_clean = 0
   do l = 1,nv
      ii = l+1; if (l == nv) ii = 1
      if (sqrt(sum((verts(:,l)-verts(:,ii))**2)) >= COSMO_EDGE_TOL) then
         nv_clean = nv_clean+1
         verts_clean(:,nv_clean) = verts(:,l)
         ccc_clean(:,nv_clean) = ccc(:,l)
         intsph_clean(nv_clean) = intsph(l)
      endif
   enddo
   if (nv_clean < 3) then
      discarded = .true.
      return
   endif
   nv = nv_clean
   verts(:,1:nv) = verts_clean(:,1:nv)
   ccc(:,1:nv) = ccc_clean(:,1:nv)
   intsph(1:nv) = intsph_clean(1:nv)
end block
end subroutine cosmo_clip_tessera

subroutine cosmo_clip_tessera_deriv(tri_verts, i_atom, natoms, Rvdw, coor_bohr, &
                                     verts, ccc, intsph, nv, discarded, &
                                     vJacI, vAtomA, vJacA, vAtomB, vJacB, cJacI, cJacGov)
real(8),intent(in) :: tri_verts(3,3)
integer,intent(in) :: i_atom, natoms
real(8),intent(in) :: Rvdw(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(out) :: verts(3,MAXV), ccc(3,MAXV)
integer,intent(out) :: intsph(MAXV)
integer,intent(out) :: nv
logical,intent(out) :: discarded
real(8),intent(out) :: vJacI(3,3,MAXV)
integer,intent(out) :: vAtomA(MAXV)
real(8),intent(out) :: vJacA(3,3,MAXV)
integer,intent(out) :: vAtomB(MAXV)
real(8),intent(out) :: vJacB(3,3,MAXV)
real(8),intent(out) :: cJacI(3,3,MAXV)
real(8),intent(out) :: cJacGov(3,3,MAXV)

real(8) :: verts_new(3,MAXV), ccc_new(3,MAXV)
integer :: intsph_new(MAXV)
real(8) :: vJacI_new(3,3,MAXV), vJacA_new(3,3,MAXV), vJacB_new(3,3,MAXV)
integer :: vAtomA_new(MAXV), vAtomB_new(MAXV)
real(8) :: cJacI_new(3,3,MAXV), cJacGov_new(3,3,MAXV)
integer :: nv_new
integer :: ind(MAXV), ltyp(MAXV)
real(8) :: pointl(3,MAXV)
integer :: ja, l, iv1, iv2, n, ncov, ncut, ii
real(8) :: p1(3), p2(3), p3(3), p4(3), point(3), rc
logical :: hit
logical :: has_jprev
real(8) :: Cjprev(3), Rjprev, dvi(3,3), dvprev(3,3), dvnew(3,3)
real(8) :: ccval(3), dccI(3,3), dccJa(3,3)
real(8) :: Idn(3,3)
integer :: aa,bb

Idn = 0.0d0
do aa = 1,3
   Idn(aa,aa) = 1.0d0
enddo

discarded = .false.
nv = 3
verts(:,1:3) = tri_verts
do l = 1,3
   ccc(:,l) = coor_bohr(:,i_atom)
   intsph(l) = i_atom
   vJacI(:,:,l) = Idn
   vAtomA(l) = 0; vJacA(:,:,l) = 0.0d0
   vAtomB(l) = 0; vJacB(:,:,l) = 0.0d0
   cJacI(:,:,l) = Idn
   cJacGov(:,:,l) = 0.0d0
enddo

do ja = 1,natoms
   if (ja == i_atom) cycle

   ncov = 0
   do n = 1,nv
      ind(n) = 0
      if (sqrt(sum((verts(:,n)-coor_bohr(:,ja))**2)) - Rvdw(ja) < 0.0d0) then
         ind(n) = 1
         ncov = ncov+1
      endif
   enddo
   if (ncov == nv) then
      discarded = .true.
      return
   endif
   if (ncov == 0) cycle

   do l = 1,nv
      iv1 = l
      iv2 = l+1; if (l == nv) iv2 = 1
      if (ind(iv1) == 1 .and. ind(iv2) == 1) then
         ltyp(l) = 0
      else if (ind(iv1) == 0 .and. ind(iv2) == 1) then
         ltyp(l) = 1
      else if (ind(iv1) == 1 .and. ind(iv2) == 0) then
         ltyp(l) = 2
      else
         ltyp(l) = 4
         rc = sqrt(sum((ccc(:,l)-verts(:,iv1))**2))
         hit = .false.
         do ii = 1,11
            point = verts(:,iv1) + real(ii,8)*(verts(:,iv2)-verts(:,iv1))/11.0d0
            point = point - ccc(:,l)
            point = point*rc/sqrt(sum(point*point)) + ccc(:,l)
            if (sqrt(sum((point-coor_bohr(:,ja))**2)) - Rvdw(ja) < -1.0d-10) then
               ltyp(l) = 3
               pointl(:,l) = point
               hit = .true.
               exit
            endif
         enddo
      endif
   enddo

   ncut = 0
   do l = 1,nv
      if (ltyp(l) == 1 .or. ltyp(l) == 2) ncut = ncut+1
      if (ltyp(l) == 3) ncut = ncut+2
   enddo
   if (ncut/2 > 1) then
      discarded = .true.
      return
   endif

   n = 0
   do l = 1,nv
      iv1 = l
      iv2 = l+1; if (l == nv) iv2 = 1

      has_jprev = (intsph(iv1) /= i_atom)
      if (has_jprev) then
         Cjprev = coor_bohr(:,intsph(iv1)); Rjprev = Rvdw(intsph(iv1))
      endif

      select case (ltyp(l))
      case (0)
      case (1)
         n = n+1
         verts_new(:,n) = verts(:,iv1); ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
         vJacI_new(:,:,n) = vJacI(:,:,iv1)
         vAtomA_new(n) = vAtomA(iv1); vJacA_new(:,:,n) = vJacA(:,:,iv1)
         vAtomB_new(n) = vAtomB(iv1); vJacB_new(:,:,n) = vJacB(:,:,iv1)
         cJacI_new(:,:,n) = cJacI(:,:,iv1); cJacGov_new(:,:,n) = cJacGov(:,:,iv1)

         p1 = verts(:,iv1); p2 = verts(:,iv2); p3 = ccc(:,iv1)
         p4 = cosmo_inter(p1, p2, p3, sqrt(sum((p1-p3)**2)), coor_bohr(:,ja), Rvdw(ja), 0)
         call cosmo_vertex_jac(p4, coor_bohr(:,i_atom), Rvdw(i_atom), has_jprev, Cjprev, Rjprev, &
                                verts(:,iv1), verts(:,iv2), coor_bohr(:,ja), Rvdw(ja), &
                                dvi, dvprev, dvnew)
         n = n+1
         verts_new(:,n) = p4
         call cosmo_two_sphere_circle_center_deriv(coor_bohr(:,i_atom), Rvdw(i_atom), &
                                                     coor_bohr(:,ja), Rvdw(ja), ccval, dccI, dccJa)
         ccc_new(:,n) = ccval
         intsph_new(n) = ja
         vJacI_new(:,:,n) = dvi
         if (has_jprev) then
            vAtomA_new(n) = intsph(iv1); vJacA_new(:,:,n) = dvprev
         else
            vAtomA_new(n) = 0; vJacA_new(:,:,n) = 0.0d0
         endif
         vAtomB_new(n) = ja; vJacB_new(:,:,n) = dvnew
         cJacI_new(:,:,n) = dccI; cJacGov_new(:,:,n) = dccJa
      case (2)
         p1 = verts(:,iv1); p2 = verts(:,iv2); p3 = ccc(:,iv1)
         p4 = cosmo_inter(p1, p2, p3, sqrt(sum((p1-p3)**2)), coor_bohr(:,ja), Rvdw(ja), 1)
         call cosmo_vertex_jac(p4, coor_bohr(:,i_atom), Rvdw(i_atom), has_jprev, Cjprev, Rjprev, &
                                verts(:,iv1), verts(:,iv2), coor_bohr(:,ja), Rvdw(ja), &
                                dvi, dvprev, dvnew)
         n = n+1
         verts_new(:,n) = p4; ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
         vJacI_new(:,:,n) = dvi
         if (has_jprev) then
            vAtomA_new(n) = intsph(iv1); vJacA_new(:,:,n) = dvprev
         else
            vAtomA_new(n) = 0; vJacA_new(:,:,n) = 0.0d0
         endif
         vAtomB_new(n) = ja; vJacB_new(:,:,n) = dvnew
         cJacI_new(:,:,n) = cJacI(:,:,iv1); cJacGov_new(:,:,n) = cJacGov(:,:,iv1)
      case (3)
         n = n+1
         verts_new(:,n) = verts(:,iv1); ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
         vJacI_new(:,:,n) = vJacI(:,:,iv1)
         vAtomA_new(n) = vAtomA(iv1); vJacA_new(:,:,n) = vJacA(:,:,iv1)
         vAtomB_new(n) = vAtomB(iv1); vJacB_new(:,:,n) = vJacB(:,:,iv1)
         cJacI_new(:,:,n) = cJacI(:,:,iv1); cJacGov_new(:,:,n) = cJacGov(:,:,iv1)

         p1 = verts(:,iv1); p2 = pointl(:,l); p3 = ccc(:,iv1)
         p4 = cosmo_inter(p1, p2, p3, sqrt(sum((p1-p3)**2)), coor_bohr(:,ja), Rvdw(ja), 0)
         call cosmo_vertex_jac(p4, coor_bohr(:,i_atom), Rvdw(i_atom), has_jprev, Cjprev, Rjprev, &
                                verts(:,iv1), verts(:,iv2), coor_bohr(:,ja), Rvdw(ja), &
                                dvi, dvprev, dvnew)
         n = n+1
         verts_new(:,n) = p4
         call cosmo_two_sphere_circle_center_deriv(coor_bohr(:,i_atom), Rvdw(i_atom), &
                                                     coor_bohr(:,ja), Rvdw(ja), ccval, dccI, dccJa)
         ccc_new(:,n) = ccval
         intsph_new(n) = ja
         vJacI_new(:,:,n) = dvi
         if (has_jprev) then
            vAtomA_new(n) = intsph(iv1); vJacA_new(:,:,n) = dvprev
         else
            vAtomA_new(n) = 0; vJacA_new(:,:,n) = 0.0d0
         endif
         vAtomB_new(n) = ja; vJacB_new(:,:,n) = dvnew
         cJacI_new(:,:,n) = dccI; cJacGov_new(:,:,n) = dccJa

         p1 = pointl(:,l); p2 = verts(:,iv2); p3 = ccc(:,iv1)
         p4 = cosmo_inter(p1, p2, p3, sqrt(sum((p1-p3)**2)), coor_bohr(:,ja), Rvdw(ja), 1)
         call cosmo_vertex_jac(p4, coor_bohr(:,i_atom), Rvdw(i_atom), has_jprev, Cjprev, Rjprev, &
                                verts(:,iv1), verts(:,iv2), coor_bohr(:,ja), Rvdw(ja), &
                                dvi, dvprev, dvnew)
         n = n+1
         verts_new(:,n) = p4; ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
         vJacI_new(:,:,n) = dvi
         if (has_jprev) then
            vAtomA_new(n) = intsph(iv1); vJacA_new(:,:,n) = dvprev
         else
            vAtomA_new(n) = 0; vJacA_new(:,:,n) = 0.0d0
         endif
         vAtomB_new(n) = ja; vJacB_new(:,:,n) = dvnew
         cJacI_new(:,:,n) = cJacI(:,:,iv1); cJacGov_new(:,:,n) = cJacGov(:,:,iv1)
      case (4)
         n = n+1
         verts_new(:,n) = verts(:,iv1); ccc_new(:,n) = ccc(:,iv1); intsph_new(n) = intsph(iv1)
         vJacI_new(:,:,n) = vJacI(:,:,iv1)
         vAtomA_new(n) = vAtomA(iv1); vJacA_new(:,:,n) = vJacA(:,:,iv1)
         vAtomB_new(n) = vAtomB(iv1); vJacB_new(:,:,n) = vJacB(:,:,iv1)
         cJacI_new(:,:,n) = cJacI(:,:,iv1); cJacGov_new(:,:,n) = cJacGov(:,:,iv1)
      end select
      if (n > MAXV-3) then
         discarded = .true.
         return
      endif
   enddo

   if (n < 3) then
      discarded = .true.
      return
   endif
   nv_new = n
   verts(:,1:nv_new) = verts_new(:,1:nv_new)
   ccc(:,1:nv_new) = ccc_new(:,1:nv_new)
   intsph(1:nv_new) = intsph_new(1:nv_new)
   vJacI(:,:,1:nv_new) = vJacI_new(:,:,1:nv_new)
   vAtomA(1:nv_new) = vAtomA_new(1:nv_new); vJacA(:,:,1:nv_new) = vJacA_new(:,:,1:nv_new)
   vAtomB(1:nv_new) = vAtomB_new(1:nv_new); vJacB(:,:,1:nv_new) = vJacB_new(:,:,1:nv_new)
   cJacI(:,:,1:nv_new) = cJacI_new(:,:,1:nv_new)
   cJacGov(:,:,1:nv_new) = cJacGov_new(:,:,1:nv_new)
   nv = nv_new
enddo

block
   real(8) :: verts_clean(3,MAXV), ccc_clean(3,MAXV)
   integer :: intsph_clean(MAXV), nv_clean, i2
   real(8) :: vJacI_c(3,3,MAXV), vJacA_c(3,3,MAXV), vJacB_c(3,3,MAXV)
   integer :: vAtomA_c(MAXV), vAtomB_c(MAXV)
   real(8) :: cJacI_c(3,3,MAXV), cJacGov_c(3,3,MAXV)
   nv_clean = 0
   do l = 1,nv
      ii = l+1; if (l == nv) ii = 1
      if (sqrt(sum((verts(:,l)-verts(:,ii))**2)) >= COSMO_EDGE_TOL) then
         nv_clean = nv_clean+1
         verts_clean(:,nv_clean) = verts(:,l)
         ccc_clean(:,nv_clean) = ccc(:,l)
         intsph_clean(nv_clean) = intsph(l)
         vJacI_c(:,:,nv_clean) = vJacI(:,:,l)
         vAtomA_c(nv_clean) = vAtomA(l); vJacA_c(:,:,nv_clean) = vJacA(:,:,l)
         vAtomB_c(nv_clean) = vAtomB(l); vJacB_c(:,:,nv_clean) = vJacB(:,:,l)
         cJacI_c(:,:,nv_clean) = cJacI(:,:,l)
         cJacGov_c(:,:,nv_clean) = cJacGov(:,:,l)
      endif
   enddo
   if (nv_clean < 3) then
      discarded = .true.
      return
   endif
   nv = nv_clean
   verts(:,1:nv) = verts_clean(:,1:nv)
   ccc(:,1:nv) = ccc_clean(:,1:nv)
   intsph(1:nv) = intsph_clean(1:nv)
   vJacI(:,:,1:nv) = vJacI_c(:,:,1:nv)
   vAtomA(1:nv) = vAtomA_c(1:nv); vJacA(:,:,1:nv) = vJacA_c(:,:,1:nv)
   vAtomB(1:nv) = vAtomB_c(1:nv); vJacB(:,:,1:nv) = vJacB_c(:,:,1:nv)
   cJacI(:,:,1:nv) = cJacI_c(:,:,1:nv)
   cJacGov(:,:,1:nv) = cJacGov_c(:,:,1:nv)
end block
end subroutine cosmo_clip_tessera_deriv

end module mod_cosmo_clip
