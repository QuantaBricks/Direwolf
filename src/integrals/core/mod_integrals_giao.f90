! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! GIAO integrals for NMR shieldings: first-order (magnetic field) one-electron matrices, the nuclear-moment operators and the two-electron GIAO derivative contractions.

submodule (mod_integrals) giao_impl
implicit none
contains

subroutine giao1e_engine(which, buf, shls, atm, natm, bas, nbas, env)
implicit none
integer,intent(in) :: which
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),parameter :: opt = 0_8
if (engine_puream) then
   select case (which)
   case (1); call cint1e_igovlp_sph(buf, shls, atm, natm, bas, nbas, env, opt)
   case (2); call cint1e_giao_irjxp_sph(buf, shls, atm, natm, bas, nbas, env, opt)
   case (3); call cint1e_igkin_sph(buf, shls, atm, natm, bas, nbas, env, opt)
   case (4); call cint1e_ignuc_sph(buf, shls, atm, natm, bas, nbas, env, opt)
   case (5); call cint1e_prinvxp_sph(buf, shls, atm, natm, bas, nbas, env, opt)
   case (6); call cint1e_giao_a11part_sph(buf, shls, atm, natm, bas, nbas, env, opt)
   case (7); call cint1e_a01gp_sph(buf, shls, atm, natm, bas, nbas, env, opt)
   end select
else
   select case (which)
   case (1); call cint1e_igovlp_cart(buf, shls, atm, natm, bas, nbas, env, opt)
   case (2); call cint1e_giao_irjxp_cart(buf, shls, atm, natm, bas, nbas, env, opt)
   case (3); call cint1e_igkin_cart(buf, shls, atm, natm, bas, nbas, env, opt)
   case (4); call cint1e_ignuc_cart(buf, shls, atm, natm, bas, nbas, env, opt)
   case (5); call cint1e_prinvxp_cart(buf, shls, atm, natm, bas, nbas, env, opt)
   case (6); call cint1e_giao_a11part_cart(buf, shls, atm, natm, bas, nbas, env, opt)
   case (7); call cint1e_a01gp_cart(buf, shls, atm, natm, bas, nbas, env, opt)
   end select
endif
end subroutine giao1e_engine

subroutine ig1_engine(buf, shls, atm, natm, bas, nbas, env)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),parameter :: opt = 0_8
if (engine_puream) then
   call cint2e_ig1_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2e_ig1_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine ig1_engine

module subroutine integrals_giao_1e(nbf, S1, H1)
implicit none
integer,intent(in) :: nbf
real(8),intent(out) :: S1(nbf,nbf,3), H1(nbf,nbf,3)
integer,allocatable :: ao_offset(:)
real(8),allocatable :: buf(:)
integer :: i0, j0, di, dj, n, num, e1, e2, c, shls(4), w, ia, ja, kb
real(8) :: fac
allocate(ao_offset(0:nBases-1))
num = 0 ; n = 0
do i0 = 0,nBases-1
   ao_offset(i0) = num
   num = num + cgto_engine(i0,bas)
   n = max(n, cgto_engine(i0,bas))
enddo
S1 = 0.0d0 ; H1 = 0.0d0
!$omp parallel private(buf, i0, j0, di, dj, e1, e2, c, shls, w, fac, ia, ja, kb)
allocate(buf(9*n*n))
!$omp do schedule(dynamic)
do i0 = 0,nBases-1
   do j0 = 0,nBases-1
      shls(1) = i0 ; shls(2) = j0 ; shls(3) = 0 ; shls(4) = 0
      di = cgto_engine(i0,bas) ; dj = cgto_engine(j0,bas)
      do w = 1,4
         buf = 0.0d0
         if (w == 4) then
            call giao1e_engine(4, buf, shls, atm_nuc, nAtoms_nuc, bas, nBases, env_nuc)
         else
            call giao1e_engine(w, buf, shls, atm, size(atm,2), bas, nBases, env)
         endif
         do c = 1,3
            do e2 = 1,dj
               do e1 = 1,di
                  ia = ao_offset(i0) + e1 ; ja = ao_offset(j0) + e2
                  fac = NorVEC(ia)*NorVEC(ja)
                  kb = e1 + di*((e2-1) + dj*(c-1))
                  select case (w)
                  case (1); S1(ia,ja,c) = buf(kb)*fac
                  case (2); H1(ia,ja,c) = H1(ia,ja,c) + 0.5d0*buf(kb)*fac
                  case default; H1(ia,ja,c) = H1(ia,ja,c) + buf(kb)*fac
                  end select
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf)
!$omp end parallel
if (necpbas > 0) call giao_ecp_add(nbf, ao_offset, H1)
end subroutine integrals_giao_1e

module subroutine integrals_nmr_nucleus(nbf, rn, PSO, A11, A01)
implicit none
integer,intent(in) :: nbf
real(8),intent(in) :: rn(3)
real(8),intent(out) :: PSO(nbf,nbf,3), A11(nbf,nbf,9), A01(nbf,nbf,9)
integer,allocatable :: ao_offset(:)
real(8),allocatable :: buf(:), env_loc(:)
integer,parameter :: PTR_RINV_ORIG = 4
integer :: i0, j0, di, dj, n, num, e1, e2, c, shls(4), kb
real(8) :: fac, ai, aj
real(8),allocatable :: amin(:), ctr(:,:)
allocate(ao_offset(0:nBases-1), amin(0:nBases-1), ctr(3,0:nBases-1))
num = 0 ; n = 0
do i0 = 0,nBases-1
   ao_offset(i0) = num
   num = num + cgto_engine(i0,bas)
   n = max(n, cgto_engine(i0,bas))
   amin(i0) = minval(env(bas(6,i0+1)+1:bas(6,i0+1)+bas(3,i0+1)))
   ctr(:,i0) = env(atm(2,bas(1,i0+1)+1)+1:atm(2,bas(1,i0+1)+1)+3)
enddo
allocate(env_loc(size(env)))
env_loc = env
env_loc(PTR_RINV_ORIG+1:PTR_RINV_ORIG+3) = rn
PSO = 0.0d0 ; A11 = 0.0d0 ; A01 = 0.0d0
!$omp parallel private(buf, i0, j0, di, dj, e1, e2, c, shls, fac, kb, ai, aj)
allocate(buf(9*n*n))
!$omp do schedule(dynamic)
do i0 = 0,nBases-1
   do j0 = 0,nBases-1
      ai = amin(i0) ; aj = amin(j0)
      if (ai*aj/(ai + aj)*sum((ctr(:,i0) - ctr(:,j0))**2) >= 30.0d0) cycle
      shls(1) = i0 ; shls(2) = j0 ; shls(3) = 0 ; shls(4) = 0
      di = cgto_engine(i0,bas) ; dj = cgto_engine(j0,bas)
      buf = 0.0d0
      call giao1e_engine(5, buf, shls, atm, size(atm,2), bas, nBases, env_loc)
      do c = 1,3
         do e2 = 1,dj
            do e1 = 1,di
               fac = NorVEC(ao_offset(i0)+e1)*NorVEC(ao_offset(j0)+e2)
               kb = e1 + di*((e2-1) + dj*(c-1))
               PSO(ao_offset(i0)+e1,ao_offset(j0)+e2,c) = buf(kb)*fac
            enddo
         enddo
      enddo
      buf = 0.0d0
      call giao1e_engine(6, buf, shls, atm, size(atm,2), bas, nBases, env_loc)
      do c = 1,9
         do e2 = 1,dj
            do e1 = 1,di
               fac = NorVEC(ao_offset(i0)+e1)*NorVEC(ao_offset(j0)+e2)
               kb = e1 + di*((e2-1) + dj*(c-1))
               A11(ao_offset(i0)+e1,ao_offset(j0)+e2,c) = buf(kb)*fac
            enddo
         enddo
      enddo
      buf = 0.0d0
      call giao1e_engine(7, buf, shls, atm, size(atm,2), bas, nBases, env_loc)
      do c = 1,9
         do e2 = 1,dj
            do e1 = 1,di
               fac = NorVEC(ao_offset(i0)+e1)*NorVEC(ao_offset(j0)+e2)
               kb = e1 + di*((e2-1) + dj*(c-1))
               A01(ao_offset(i0)+e1,ao_offset(j0)+e2,c) = buf(kb)*fac
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf)
!$omp end parallel
end subroutine integrals_nmr_nucleus

module subroutine integrals_giao_jk(nbf, D, Jx, Kx, do_k, omega)
implicit none
integer,intent(in) :: nbf
real(8),intent(in) :: D(nbf,nbf)
real(8),intent(out) :: Jx(nbf,nbf,3), Kx(nbf,nbf,3)
logical,intent(in) :: do_k
real(8),intent(in) :: omega
real(8),allocatable :: envw(:)
integer,allocatable :: ao_offset(:)
real(8),allocatable :: buf(:), Jt(:,:,:), Kt(:,:,:), Q(:,:)
integer :: I, J, K, L, di, dj, dk, dl, n, num, a, b, c, d4, x, shls(4), ia, ja, ka, la, nthr
real(8) :: v, qmax, thr
real(8),allocatable :: Q2(:,:)
logical :: first
allocate(ao_offset(0:nBases-1))
num = 0 ; n = 0
do I = 0,nBases-1
   ao_offset(I) = num
   num = num + cgto_engine(I,bas)
   n = max(n, cgto_engine(I,bas))
enddo
allocate(Q(0:nBases-1,0:nBases-1))
Q = 0.0d0
!$omp parallel private(buf, I, J, di, dj, a, b, shls, v)
allocate(buf(n**4))
!$omp do schedule(dynamic)
do I = 0,nBases-1
   do J = 0,I
      shls = [I, J, I, J]
      di = cgto_engine(I,bas) ; dj = cgto_engine(J,bas)
      buf = 0.0d0
      call twoe_engine(buf, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      v = 0.0d0
      do b = 1,dj
         do a = 1,di
            v = max(v, abs(buf(a + di*((b-1) + dj*((a-1) + di*(b-1))))))
         enddo
      enddo
      Q(I,J) = sqrt(v) ; Q(J,I) = Q(I,J)
   enddo
enddo
!$omp end do
deallocate(buf)
!$omp end parallel
qmax = maxval(Q)
thr = 1.0d-12
allocate(envw(size(env)))
envw = env ; envw(9) = omega
Jx = 0.0d0 ; Kx = 0.0d0
!$omp parallel private(buf, Jt, Kt, I, J, K, L, di, dj, dk, dl, a, b, c, d4, x, shls, ia, ja, ka, la, v)
allocate(buf(3*n**4), Jt(nbf,nbf,3), Kt(nbf,nbf,3))
Jt = 0.0d0 ; Kt = 0.0d0
!$omp do schedule(dynamic)
do I = 0,nBases-1
   do J = 0,nBases-1
      if (Q(I,J)*qmax*50.0d0 < thr) cycle
      di = cgto_engine(I,bas) ; dj = cgto_engine(J,bas)
      do K = 0,nBases-1
         do L = 0,nBases-1
            if (Q(I,J)*Q(K,L)*50.0d0 < thr) cycle
            dk = cgto_engine(K,bas) ; dl = cgto_engine(L,bas)
            shls = [I, J, K, L]
            buf(1:3*di*dj*dk*dl) = 0.0d0
            call ig1_engine(buf, shls, atm, size(atm,2), bas, nBases, envw)
            do x = 1,3
               do d4 = 1,dl
                  la = ao_offset(L) + d4
                  do c = 1,dk
                     ka = ao_offset(K) + c
                     do b = 1,dj
                        ja = ao_offset(J) + b
                        do a = 1,di
                           ia = ao_offset(I) + a
                           v = buf(a + di*((b-1) + dj*((c-1) + dk*((d4-1) + dl*(x-1))))) &
                               *NorVEC(ia)*NorVEC(ja)*NorVEC(ka)*NorVEC(la)
                           Jt(ia,ja,x) = Jt(ia,ja,x) + v*D(la,ka)
                           if (do_k) then
                              Kt(ia,la,x) = Kt(ia,la,x) + v*D(ja,ka)
                              Kt(ka,ja,x) = Kt(ka,ja,x) + v*D(la,ia)
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
!$omp critical
Jx = Jx + Jt
Kx = Kx + Kt
!$omp end critical
deallocate(buf, Jt, Kt)
!$omp end parallel
end subroutine integrals_giao_jk

module subroutine integrals_k_general(nbf, nd, dm, K, omega)
implicit none
integer,intent(in) :: nbf, nd
real(8),intent(in) :: dm(nbf,nbf,nd)
real(8),intent(out) :: K(nbf,nbf,nd)
real(8),intent(in) :: omega
real(8),allocatable :: envw(:)
integer,allocatable :: ao_offset(:)
real(8),allocatable :: buf(:), Kt(:,:,:), Q(:,:)
integer :: I, J, KK, L, di, dj, dk, dl, n, num, a, b, c, d4, id, shls(4), ia, ja, ka, la
integer :: ij, kl
real(8) :: v, f, vmax, qmax, thr
allocate(ao_offset(0:nBases-1))
num = 0 ; n = 0
do I = 0,nBases-1
   ao_offset(I) = num
   num = num + cgto_engine(I,bas)
   n = max(n, cgto_engine(I,bas))
enddo
allocate(Q(0:nBases-1,0:nBases-1))
Q = 0.0d0
!$omp parallel private(buf, I, J, di, dj, a, b, shls, v)
allocate(buf(n**4))
!$omp do schedule(dynamic)
do I = 0,nBases-1
   do J = 0,I
      shls = [I, J, I, J]
      di = cgto_engine(I,bas) ; dj = cgto_engine(J,bas)
      buf = 0.0d0
      call twoe_engine(buf, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      v = 0.0d0
      do b = 1,dj
         do a = 1,di
            v = max(v, abs(buf(a + di*((b-1) + dj*((a-1) + di*(b-1))))))
         enddo
      enddo
      Q(I,J) = sqrt(v) ; Q(J,I) = Q(I,J)
   enddo
enddo
!$omp end do
deallocate(buf)
!$omp end parallel
thr = 1.0d-12
allocate(envw(size(env)))
envw = env ; envw(9) = omega
K = 0.0d0
!$omp parallel private(buf, Kt, I, J, KK, L, di, dj, dk, dl, a, b, c, d4, id, shls, ia, ja, ka, la, v, f, ij, kl)
allocate(buf(n**4), Kt(nbf,nbf,nd))
Kt = 0.0d0
!$omp do schedule(dynamic)
do I = 0,nBases-1
   do J = 0,I
      ij = I*(I+1)/2 + J
      di = cgto_engine(I,bas) ; dj = cgto_engine(J,bas)
      do KK = 0,I
         do L = 0,KK
            kl = KK*(KK+1)/2 + L
            if (kl > ij) cycle
            if (Q(I,J)*Q(KK,L) < thr) cycle
            dk = cgto_engine(KK,bas) ; dl = cgto_engine(L,bas)
            shls = [I, J, KK, L]
            buf(1:di*dj*dk*dl) = 0.0d0
            call twoe_engine(buf, shls, atm, size(atm,2), bas, nBases, envw, 0_8)
            f = 1.0d0
            if (I == J) f = f*0.5d0
            if (KK == L) f = f*0.5d0
            if (ij == kl) f = f*0.5d0
            do d4 = 1,dl
               la = ao_offset(L) + d4
               do c = 1,dk
                  ka = ao_offset(KK) + c
                  do b = 1,dj
                     ja = ao_offset(J) + b
                     do a = 1,di
                        ia = ao_offset(I) + a
                        v = f*buf(a + di*((b-1) + dj*((c-1) + dk*(d4-1))))*NorVEC(ia)*NorVEC(ja)*NorVEC(ka)*NorVEC(la)
                        do id = 1,nd
                           Kt(ia,ka,id) = Kt(ia,ka,id) + v*dm(ja,la,id)
                           Kt(ja,ka,id) = Kt(ja,ka,id) + v*dm(ia,la,id)
                           Kt(ia,la,id) = Kt(ia,la,id) + v*dm(ja,ka,id)
                           Kt(ja,la,id) = Kt(ja,la,id) + v*dm(ia,ka,id)
                           Kt(ka,ia,id) = Kt(ka,ia,id) + v*dm(la,ja,id)
                           Kt(la,ia,id) = Kt(la,ia,id) + v*dm(ka,ja,id)
                           Kt(ka,ja,id) = Kt(ka,ja,id) + v*dm(la,ia,id)
                           Kt(la,ja,id) = Kt(la,ja,id) + v*dm(ka,ia,id)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
!$omp critical
K = K + Kt
!$omp end critical
deallocate(buf, Kt)
!$omp end parallel
end subroutine integrals_k_general

module subroutine integrals_ao_atom(nbf, ao_atom)
implicit none
integer,intent(in) :: nbf
integer,intent(out) :: ao_atom(nbf)
integer :: i0, k, num
integer,parameter :: ATOM_OF = 1
num = 0
do i0 = 0,nBases-1
   do k = 1,cgto_engine(i0,bas)
      num = num + 1
      ao_atom(num) = bas(ATOM_OF,i0+1) + 1
   enddo
enddo
end subroutine integrals_ao_atom

module subroutine integrals_giao_ri_j(nbf, D, Jx)
implicit none
integer,intent(in) :: nbf
real(8),intent(in) :: D(nbf,nbf)
real(8),intent(out) :: Jx(nbf,nbf,3)
integer,allocatable :: ao_offset(:), aux_offset(:)
real(8),allocatable :: buf(:), bt(:), bvec(:), cvec(:), tmp(:), Jt(:,:,:)
integer :: I, J, K, di, dj, dk, n, nk, num, a, b, c, x, shls(4), ia, ja, kp, natm_loc
real(8) :: v, f, r2, ai, aj
logical,allocatable :: keep(:,:)
real(8),allocatable :: amin(:), ctr(:,:)
call ensure_df_built()
allocate(ao_offset(0:nBases-1), aux_offset(0:nBasesAux-1))
allocate(amin(0:nBases-1), ctr(3,0:nBases-1), keep(0:nBases-1,0:nBases-1))
do I = 0,nBases-1
   amin(I) = minval(env(bas(6,I+1)+1:bas(6,I+1)+bas(3,I+1)))
   ctr(:,I) = env(atm(2,bas(1,I+1)+1)+1:atm(2,bas(1,I+1)+1)+3)
enddo
do I = 0,nBases-1
   do J = 0,nBases-1
      ai = amin(I) ; aj = amin(J)
      r2 = sum((ctr(:,I) - ctr(:,J))**2)
      keep(J,I) = (ai*aj/(ai + aj)*r2 < 30.0d0)
   enddo
enddo
num = 0 ; n = 0
do I = 0,nBases-1
   ao_offset(I) = num
   num = num + cgto_engine(I,bas)
   n = max(n, cgto_engine(I,bas))
enddo
num = 0 ; nk = 0
do K = 0,nBasesAux-1
   aux_offset(K) = num
   num = num + cgto_engine(nBases+K,basDF)
   nk = max(nk, cgto_engine(nBases+K,basDF))
enddo
natm_loc = size(atm,2)
allocate(bvec(nContsAux), cvec(nContsAux), tmp(nContsAux))
bvec = 0.0d0
!$omp parallel private(buf, bt, I, J, K, di, dj, dk, a, b, c, shls, ia, ja, kp, v, f)
allocate(buf(n*n*nk), bt(nContsAux))
bt = 0.0d0
!$omp do schedule(dynamic)
do I = 0,nBases-1
   do J = 0,I
      if (.not. keep(J,I)) cycle
      di = cgto_engine(I,bas) ; dj = cgto_engine(J,bas)
      f = 1.0d0
      if (I /= J) f = 2.0d0
      do K = 0,nBasesAux-1
         dk = cgto_engine(nBases+K,basDF)
         shls = [I, J, nBases+K, 0]
         buf(1:di*dj*dk) = 0.0d0
         call threec2e_engine(buf, shls, atm, natm_loc, basDF, nBases+nBasesAux, envDF, 0_8)
         do c = 1,dk
            kp = aux_offset(K) + c
            do b = 1,dj
               ja = ao_offset(J) + b
               do a = 1,di
                  ia = ao_offset(I) + a
                  bt(kp) = bt(kp) + f*buf(a + di*((b-1) + dj*(c-1)))*NorVEC(ia)*NorVEC(ja)*NorVECAux(kp)*D(ia,ja)
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
!$omp critical
bvec = bvec + bt
!$omp end critical
deallocate(buf, bt)
!$omp end parallel
call dgemv('T', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, bvec, 1, 0.0d0, tmp, 1)
tmp = tmp*df_evalinv
call dgemv('N', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, tmp, 1, 0.0d0, cvec, 1)
Jx = 0.0d0
!$omp parallel private(buf, Jt, I, J, K, di, dj, dk, a, b, c, x, shls, ia, ja, kp, v)
allocate(buf(3*n*n*nk), Jt(nbf,nbf,3))
Jt = 0.0d0
!$omp do schedule(dynamic)
do I = 0,nBases-1
   do J = 0,I-1
      if (bas(1,I+1) == bas(1,J+1)) cycle
      if (.not. keep(J,I)) cycle
      di = cgto_engine(I,bas) ; dj = cgto_engine(J,bas)
      do K = 0,nBasesAux-1
         dk = cgto_engine(nBases+K,basDF)
         shls = [I, J, nBases+K, 0]
         buf(1:3*di*dj*dk) = 0.0d0
         call giao3c_engine(buf, shls, atm, natm_loc, basDF, nBases+nBasesAux, envDF)
         do x = 1,3
            do c = 1,dk
               kp = aux_offset(K) + c
               do b = 1,dj
                  ja = ao_offset(J) + b
                  do a = 1,di
                     ia = ao_offset(I) + a
                     v = buf(a + di*((b-1) + dj*((c-1) + dk*(x-1))))*NorVEC(ia)*NorVEC(ja)*NorVECAux(kp)*cvec(kp)
                     Jt(ia,ja,x) = Jt(ia,ja,x) + v
                     Jt(ja,ia,x) = Jt(ja,ia,x) - v
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
!$omp critical
Jx = Jx + Jt
!$omp end critical
deallocate(buf, Jt)
!$omp end parallel
end subroutine integrals_giao_ri_j

subroutine giao3c_engine(buf, shls, atm, natm, bas_, nbas, env_)
use iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_long, c_loc
implicit none
real(8),target :: buf(*)
integer,target :: shls(4)
integer,target :: atm(*)
integer,intent(in) :: natm
integer,target :: bas_(*)
integer,intent(in) :: nbas
real(8),target :: env_(*)
interface
   function c_ig1_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(C, name='int3c2e_ig1_sph') result(r)
   import :: c_int, c_double, c_ptr, c_long
   type(c_ptr),value :: out, dims, shls, atm, bas, env, opt, cache
   integer(c_int),value :: natm, nbas
   integer(c_long) :: r
   end function c_ig1_sph
   function c_ig1_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(C, name='int3c2e_ig1_cart') result(r)
   import :: c_int, c_double, c_ptr, c_long
   type(c_ptr),value :: out, dims, shls, atm, bas, env, opt, cache
   integer(c_int),value :: natm, nbas
   integer(c_long) :: r
   end function c_ig1_cart
end interface
integer(c_long) :: rr
integer(c_int),target :: shls_c(4)
shls_c = int(shls, c_int)
if (engine_puream) then
   rr = c_ig1_sph(c_loc(buf), c_null_ptr, c_loc(shls_c), c_loc(atm), int(natm,c_int), &
                  c_loc(bas_), int(nbas,c_int), c_loc(env_), c_null_ptr, c_null_ptr)
else
   rr = c_ig1_cart(c_loc(buf), c_null_ptr, c_loc(shls_c), c_loc(atm), int(natm,c_int), &
                   c_loc(bas_), int(nbas,c_int), c_loc(env_), c_null_ptr, c_null_ptr)
endif
end subroutine giao3c_engine

module subroutine integrals_giao_cosx_k(nbf, D, nat, coor_bohr, Kx, nx, dx, Kdx)
implicit none
integer,intent(in) :: nbf, nat, nx
real(8),intent(in) :: D(nbf,nbf), coor_bohr(3,nat), dx(nbf,nbf,nx)
real(8),intent(out) :: Kx(nbf,nbf,3), Kdx(nbf,nbf,nx)
real(8),allocatable :: G(:,:,:), T(:,:,:), K1(:,:)
integer,allocatable :: ao_atom(:), kind(:)
integer :: m, a, b, c, mu, nu, sg
real(8) :: Mac(nbf), eps3(3,3,3)
eps3 = 0.0d0
eps3(1,2,3) = 1.0d0 ; eps3(2,3,1) = 1.0d0 ; eps3(3,1,2) = 1.0d0
eps3(1,3,2) = -1.0d0 ; eps3(3,2,1) = -1.0d0 ; eps3(2,1,3) = -1.0d0
allocate(ao_atom(nbf))
call integrals_ao_atom(nbf, ao_atom)
allocate(G(nbf,nbf,9+nx), T(nbf,nbf,9+nx), kind(9+nx), K1(nbf,nbf))
do c = 1,3
   G(:,:,c) = D
   kind(c) = c
enddo
m = 3
do a = 1,3
   do c = 1,3
      if (a == c) cycle
      m = m + 1
      b = 6 - a - c
      do nu = 1,nbf
         G(:,nu,m) = D(:,nu)*0.5d0*eps3(a,b,c)*coor_bohr(b,ao_atom(nu))
      enddo
      kind(m) = c
   enddo
enddo
do m = 1,nx
   G(:,:,9+m) = dx(:,:,m)
   kind(9+m) = 0
enddo
call cosx_set_bra_kind(kind, 9+nx)
call integrals_build_exchange_cosx_batch_raw(nbf, 9+nx, G, T)
call cosx_set_bra_kind(kind, 0)
do m = 1,nx
   Kdx(:,:,m) = 0.5d0*(transpose(T(:,:,9+m)) - T(:,:,9+m))
enddo
Kx = 0.0d0
do a = 1,3
   K1 = 0.0d0
   m = 3
   do sg = 1,3
      do c = 1,3
         if (sg == c) cycle
         m = m + 1
         if (sg /= a) cycle
         K1 = K1 - T(:,:,m)
         b = 6 - a - c
         do nu = 1,nbf
            do mu = 1,nbf
               K1(mu,nu) = K1(mu,nu) + 0.5d0*eps3(a,b,c)*coor_bohr(b,ao_atom(mu))*T(mu,nu,c)
            enddo
         enddo
      enddo
   enddo
   do nu = 1,nbf
      do mu = 1,nbf
         Kx(mu,nu,a) = K1(mu,nu) - K1(nu,mu)
      enddo
   enddo
enddo
end subroutine integrals_giao_cosx_k

subroutine giao_ecp_add(nbf, ao_offset, H1)
implicit none
integer,intent(in) :: nbf, ao_offset(0:nBases-1)
real(8),intent(inout) :: H1(nbf,nbf,3)
real(8),allocatable :: V(:,:,:), envC(:), buf(:), Rao(:,:)
real(8) :: eps(3,3,3), Cc(3), g
integer :: s0, s1, ia_c, i0, j0, di, dj, n, e1, e2, c, d, e, ia, ja, a1, k
eps = 0.0d0
eps(1,2,3) = 1.0d0 ; eps(2,3,1) = 1.0d0 ; eps(3,1,2) = 1.0d0
eps(1,3,2) = -1.0d0 ; eps(3,2,1) = -1.0d0 ; eps(2,1,3) = -1.0d0
n = 0
do i0 = 0,nBases-1
   n = max(n, cgto_engine(i0,bas))
enddo
allocate(V(nbf,nbf,3), envC(size(envECP)), Rao(3,nbf))
do i0 = 0,nBases-1
   a1 = bas(1,i0+1) + 1
   do k = 1,cgto_engine(i0,bas)
      Rao(:,ao_offset(i0)+k) = env(atm(2,a1)+1:atm(2,a1)+3)
   enddo
enddo
s0 = nBases + 1
do while (s0 <= nBases + necpbas)
   ia_c = basECP(1,s0)
   s1 = s0
   do while (s1 + 1 <= nBases + necpbas)
      if (basECP(1,s1+1) /= ia_c) exit
      s1 = s1 + 1
   enddo
   envC = envECP
   envC(19) = dble(s0 - 1)
   envC(20) = dble(s1 - s0 + 1)
   Cc = env(atm(2,ia_c+1)+1:atm(2,ia_c+1)+3)
   V = 0.0d0
   !$omp parallel private(buf, i0, j0, di, dj, e1, e2, c, ia, ja)
   allocate(buf(3*n*n))
   !$omp do schedule(dynamic)
   do i0 = 0,nBases-1
      do j0 = 0,nBases-1
         di = cgto_engine(i0,bas) ; dj = cgto_engine(j0,bas)
         buf = 0.0d0
         call giao_ecp_rbra(buf, di, dj, [i0, j0, 0, 0], envC)
         do c = 1,3
            do e2 = 1,dj
               do e1 = 1,di
                  ia = ao_offset(i0) + e1 ; ja = ao_offset(j0) + e2
                  V(ia,ja,c) = buf(e1 + di*((e2-1) + dj*(c-1)))*NorVEC(ia)*NorVEC(ja)
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end do
   deallocate(buf)
   !$omp end parallel
   do c = 1,3
      do ja = 1,nbf
         do ia = 1,nbf
            g = 0.0d0
            do d = 1,3
               do e = 1,3
                  if (eps(c,d,e) == 0.0d0) cycle
                  g = g + eps(c,d,e)*((Rao(d,ia) - Cc(d))*V(ia,ja,e) - (Rao(d,ja) - Cc(d))*V(ja,ia,e))
               enddo
            enddo
            H1(ia,ja,c) = H1(ia,ja,c) - 0.5d0*g
         enddo
      enddo
   enddo
   s0 = s1 + 1
enddo
end subroutine giao_ecp_add

subroutine giao_ecp_rbra(buf, di, dj, shls, envC)
use iso_c_binding, only: c_int, c_ptr, c_null_ptr, c_loc
implicit none
real(8),target :: buf(*)
integer,intent(in) :: di, dj, shls(4)
real(8),intent(in) :: envC(:)
interface
   function ecp_ig_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(C, name='ECPscalar_rbra_sph') result(r)
   import :: c_int, c_ptr
   type(c_ptr),value :: out, dims, shls, atm, bas, env, opt, cache
   integer(c_int),value :: natm, nbas
   integer(c_int) :: r
   end function ecp_ig_sph
   function ecp_ig_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(C, name='ECPscalar_rbra_cart') result(r)
   import :: c_int, c_ptr
   type(c_ptr),value :: out, dims, shls, atm, bas, env, opt, cache
   integer(c_int),value :: natm, nbas
   integer(c_int) :: r
   end function ecp_ig_cart
end interface
integer(c_int),target :: dims_c(2), shls_c(2)
integer(c_int),allocatable,target :: atm_c(:,:), bas_c(:,:)
real(8),allocatable,target :: env_c(:)
integer(c_int) :: rr
dims_c = [int(di,c_int), int(dj,c_int)]
shls_c = [int(shls(1),c_int), int(shls(2),c_int)]
allocate(atm_c(size(atm,1),size(atm,2)), bas_c(size(basECP,1),size(basECP,2)), env_c(size(envC)))
atm_c = int(atm,c_int) ; bas_c = int(basECP,c_int) ; env_c = envC
if (engine_puream) then
   rr = ecp_ig_sph(c_loc(buf), c_loc(dims_c), c_loc(shls_c), c_loc(atm_c), int(size(atm,2),c_int), c_loc(bas_c), &
                   int(nBases+necpbas,c_int), c_loc(env_c), c_null_ptr, c_null_ptr)
else
   rr = ecp_ig_cart(c_loc(buf), c_loc(dims_c), c_loc(shls_c), c_loc(atm_c), int(size(atm,2),c_int), c_loc(bas_c), &
                    int(nBases+necpbas,c_int), c_loc(env_c), c_null_ptr, c_null_ptr)
endif
end subroutine giao_ecp_rbra

module subroutine integrals_giao_charges(nbf, n, pos, q, H1x)
implicit none
integer,intent(in) :: nbf, n
real(8),intent(in) :: pos(3,n), q(n)
real(8),intent(out) :: H1x(nbf,nbf,3)
integer,allocatable :: ao_offset(:), atm_x(:,:)
real(8),allocatable :: buf(:), env_x(:)
integer :: i0, j0, di, dj, nmax, num, e1, e2, c, shls(4), ia, ja, kb, k, off, na
real(8) :: fac
na = size(atm,2)
allocate(atm_x(6,na+n), env_x(size(env)+4*n))
atm_x = 0
atm_x(:,1:na) = atm(:,1:na)
atm_x(1,1:na) = 0
env_x = 0.0d0
env_x(1:size(env)) = env
off = size(env)
do k = 1,n
   atm_x(2,na+k) = off
   env_x(off+1:off+3) = pos(:,k)
   off = off + 3
   atm_x(3,na+k) = 3
   atm_x(5,na+k) = off
   env_x(off+1) = q(k)
   off = off + 1
enddo
allocate(ao_offset(0:nBases-1))
num = 0 ; nmax = 0
do i0 = 0,nBases-1
   ao_offset(i0) = num
   num = num + cgto_engine(i0,bas)
   nmax = max(nmax, cgto_engine(i0,bas))
enddo
H1x = 0.0d0
!$omp parallel private(buf, i0, j0, di, dj, e1, e2, c, shls, fac, ia, ja, kb)
allocate(buf(3*nmax*nmax))
!$omp do schedule(dynamic)
do i0 = 0,nBases-1
   do j0 = 0,nBases-1
      shls = [i0, j0, 0, 0]
      di = cgto_engine(i0,bas) ; dj = cgto_engine(j0,bas)
      buf = 0.0d0
      call giao1e_engine(4, buf, shls, atm_x, na+n, bas, nBases, env_x)
      do c = 1,3
         do e2 = 1,dj
            do e1 = 1,di
               ia = ao_offset(i0) + e1 ; ja = ao_offset(j0) + e2
               fac = NorVEC(ia)*NorVEC(ja)
               kb = e1 + di*((e2-1) + dj*(c-1))
               H1x(ia,ja,c) = buf(kb)*fac
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf)
!$omp end parallel
end subroutine integrals_giao_charges

end submodule giao_impl
