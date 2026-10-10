! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! One-electron piece of the analytic Hessian (S + kinetic + nuclear

subroutine calc_hessian_1e(nConts, natoms, ao_atom, Ptot, Wtot, znum, ecpce, coor_bohr, d2E_1e)
use omp_lib
use mod_integrals, only: integrals_hessian_1e_shell, integrals_hessian_nuc_atom, integrals_hessian_nuc_atom_raw, &
                         integrals_hessian_ecp_atom_raw
implicit none
integer,intent(in) :: nConts, natoms
integer,intent(in) :: ao_atom(nConts)
real(8),intent(in) :: Ptot(nConts,nConts), Wtot(nConts,nConts)
integer,intent(in) :: znum(natoms), ecpce(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(out) :: d2E_1e(3,natoms,3,natoms)

real(8),allocatable :: s1aa(:,:,:), s1ab(:,:,:), h1aa(:,:,:), h1ab(:,:,:)
real(8),allocatable :: CC(:,:,:), CB(:,:,:), CK(:,:,:)
real(8),allocatable :: dp(:,:,:,:,:)
integer :: mu,nu,ia,ja,ic,nthr,tid,a,b,k,t
real(8) :: Zeff, pm, wm, x, y
real(8) :: tw0, t_shell, t_accA, t_nucint, t_nucacc, t_rawint, t_diag, t1
logical :: timing
character(len=4) :: tenv
call get_environment_variable("ENGINE_HESS1E_TIMING", tenv)
timing = (trim(tenv) .eq. "1")
t_shell = 0 ; t_accA = 0 ; t_nucint = 0 ; t_nucacc = 0 ; t_rawint = 0 ; t_diag = 0

nthr = omp_get_max_threads()
allocate(dp(3,natoms,3,natoms,0:nthr-1))

allocate(s1aa(nConts,nConts,9), s1ab(nConts,nConts,9))
allocate(h1aa(nConts,nConts,9), h1ab(nConts,nConts,9))
tw0 = omp_get_wtime()
call integrals_hessian_1e_shell(nConts, s1aa, s1ab, h1aa, h1ab)
t_shell = omp_get_wtime() - tw0
tw0 = omp_get_wtime()

d2E_1e = 0.0d0
dp = 0.0d0

!$omp parallel default(shared) private(mu,nu,ia,ja,tid,a,b,k,pm,wm,x,y)
tid = omp_get_thread_num()
!$omp do schedule(static)
do mu = 1,nConts
   ia = ao_atom(mu)
   do nu = 1,nConts
      ja = ao_atom(nu)
      pm = Ptot(mu,nu) ; wm = Wtot(mu,nu)
      do b = 1,3
         do a = 1,3
            k = b + 3*(a-1)
            x = h1aa(mu,nu,k)*pm - s1aa(mu,nu,k)*wm
            dp(a,ia,b,ia,tid) = dp(a,ia,b,ia,tid) + x
            y = h1ab(mu,nu,k)*pm - s1ab(mu,nu,k)*wm
            dp(a,ia,b,ja,tid) = dp(a,ia,b,ja,tid) + y
            dp(b,ja,a,ia,tid) = dp(b,ja,a,ia,tid) + y
         enddo
      enddo
   enddo
enddo
!$omp end do
!$omp end parallel
call add_partials()
t_accA = omp_get_wtime() - tw0

block
   real(8),allocatable :: AA(:,:,:), AB(:,:,:), AAe(:,:,:), ABe(:,:,:)
   real(8) :: dh_l(3,3), ds_l(3,3), d2diag(3,3)
   real(8),allocatable :: dh(:,:,:), ds(:,:,:)
   real(8) :: ccv, cbv, ckv, zb, hblk(9), sblk(9)
   logical :: hasecp
   integer :: p,q,kt
   allocate(AA(nConts,nConts,9), AB(nConts,nConts,9), AAe(nConts,nConts,9), ABe(nConts,nConts,9))
   allocate(dh(3,3,0:nthr-1), ds(3,3,0:nthr-1))
   do ic = 1,natoms
      Zeff = real(znum(ic)-ecpce(ic),8)
      t1 = omp_get_wtime()
      call integrals_hessian_nuc_atom_raw(coor_bohr(:,ic), nConts, AA, AB)
      call integrals_hessian_ecp_atom_raw(ic, nConts, AAe, ABe, hasecp)
      t_nucint = t_nucint + omp_get_wtime() - t1
      t1 = omp_get_wtime()
      zb = Zeff
      if (hasecp) zb = 1.0d0
      dp = 0.0d0
      !$omp parallel default(shared) private(mu,nu,ia,ja,tid,a,b,k,kt,pm,x,y,ccv,cbv,ckv)
      tid = omp_get_thread_num()
      !$omp do schedule(static)
      do mu = 1,nConts
         ia = ao_atom(mu)
         do nu = 1,nConts
            ja = ao_atom(nu)
            pm = Ptot(mu,nu)
            do b = 1,3
               do a = 1,3
                  k = b + 3*(a-1)
                  ccv = AA(mu,nu,k) + 2.0d0*AB(mu,nu,k) + AA(nu,mu,k)
                  if (hasecp) ccv = Zeff*ccv - (AAe(mu,nu,k) + 2.0d0*ABe(mu,nu,k) + AAe(nu,mu,k))
                  dp(a,ic,b,ic,tid) = dp(a,ic,b,ic,tid) - zb*ccv*pm
                  kt = a + 3*(b-1)
                  cbv = -(AA(mu,nu,kt) + AB(mu,nu,kt))
                  if (hasecp) cbv = Zeff*cbv + (AAe(mu,nu,kt) + ABe(mu,nu,kt))
                  x = -zb*cbv*pm
                  dp(a,ic,b,ia,tid) = dp(a,ic,b,ia,tid) + x
                  dp(b,ia,a,ic,tid) = dp(b,ia,a,ic,tid) + x
                  ckv = -(AB(mu,nu,k) + AA(nu,mu,k))
                  if (hasecp) ckv = Zeff*ckv + (ABe(mu,nu,k) + AAe(nu,mu,k))
                  y = -zb*ckv*pm
                  dp(a,ic,b,ja,tid) = dp(a,ic,b,ja,tid) + y
                  dp(b,ja,a,ic,tid) = dp(b,ja,a,ic,tid) + y
               enddo
            enddo
         enddo
      enddo
      !$omp end do
      !$omp end parallel
      call add_partials()
      t_nucacc = t_nucacc + omp_get_wtime() - t1
      t1 = omp_get_wtime()
      if (hasecp) then
         !$omp parallel do collapse(2) default(shared) private(kt,q,p) schedule(static)
         do kt = 1,9
            do q = 1,nConts
               do p = 1,nConts
                  AA(p,q,kt) = Zeff*AA(p,q,kt) - AAe(p,q,kt)
                  AB(p,q,kt) = Zeff*AB(p,q,kt) - ABe(p,q,kt)
               enddo
            enddo
         enddo
         !$omp end parallel do
         Zeff = 1.0d0
      endif
      dh = 0.0d0 ; ds = 0.0d0
      !$omp parallel default(shared) private(p,q,tid,a,b,hblk,sblk,dh_l,ds_l)
      tid = omp_get_thread_num()
      dh_l = 0.0d0 ; ds_l = 0.0d0
      !$omp do schedule(static)
      do p = 1,nConts
         do q = 1,nConts
            hblk = -Zeff*AA(p,q,:) - Zeff*AB(p,q,:)
            sblk = 0.0d0
            if (ao_atom(p) == ic) then
               hblk = hblk + h1aa(p,q,:) + Zeff*AA(p,q,:) + Zeff*AB(p,q,:)
               sblk = sblk + s1aa(p,q,:)
            endif
            if (ao_atom(q) == ic) then
               hblk = hblk + Zeff*AA(q,p,:) + Zeff*AB(p,q,:)
            endif
            if (ao_atom(p) == ic .and. ao_atom(q) == ic) then
               hblk = hblk + h1ab(p,q,:)
               sblk = sblk + s1ab(p,q,:)
            endif
            do b = 1,3
               do a = 1,3
                  dh_l(a,b) = dh_l(a,b) + hblk(b+3*(a-1))*Ptot(p,q)
                  ds_l(a,b) = ds_l(a,b) + sblk(b+3*(a-1))*Wtot(p,q)
               enddo
            enddo
         enddo
      enddo
      !$omp end do
      dh(:,:,tid) = dh_l ; ds(:,:,tid) = ds_l
      !$omp end parallel
      d2diag = 0.0d0
      do t = 0,nthr-1
         d2diag = d2diag + 2.0d0*dh(:,:,t) - 2.0d0*ds(:,:,t)
      enddo
      d2E_1e(:,ic,:,ic) = d2diag
      t_diag = t_diag + omp_get_wtime() - t1
   enddo
   deallocate(AA,AB,AAe,ABe,dh,ds)
end block
deallocate(s1aa, s1ab, h1aa, h1ab, dp)
if (timing) print '(A,5F8.3)', ' HESS1E shell,accA,ints,nuc_acc,diag:', &
   t_shell, t_accA, t_nucint, t_nucacc, t_diag

contains
   subroutine add_partials()
   integer :: tt
   do tt = 0,nthr-1
      d2E_1e = d2E_1e + dp(:,:,:,:,tt)
   enddo
   end subroutine add_partials

end subroutine calc_hessian_1e
