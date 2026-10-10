! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_hessian: analytic COSMO (electrostatic) contribution to the

module mod_cosmo_hessian
use mod_cosmo_cavity, only: cosmo_build_cavity_deriv, bohr2ang, PI_cosmo
use mod_cosmo_state
implicit none

logical, save :: cosmo_hess_active = .false.
real(8), allocatable, save :: cosmo_hess_Fx(:,:,:,:)

contains

subroutine cosmo_hess_release()
cosmo_hess_active = .false.
if (allocated(cosmo_hess_Fx)) deallocate(cosmo_hess_Fx)
end subroutine cosmo_hess_release

subroutine cosmo_hess_cavity(natoms, coor_bohr, nt, tc, ta, Js, Ja)
use mod_cosmo_cavity, only: cosmo_smd_radius_bohr
implicit none
integer,intent(in) :: natoms
real(8),intent(in) :: coor_bohr(3,natoms)
integer,intent(out) :: nt
real(8),allocatable,intent(out) :: tc(:,:), ta(:), Js(:,:,:), Ja(:,:)
real(8),allocatable :: d_tc(:,:), d_ta(:), Rvdw(:), d_dA(:,:,:), d_dC(:,:,:,:)
integer,allocatable :: d_tat(:), d_nslot(:), d_aos(:,:)
integer :: d_n, k, kk, s, ia, j, ia2

if (use_smd_radii) then
   allocate(Rvdw(natoms))
   do ia = 1,natoms
      Rvdw(ia) = cosmo_smd_radius_bohr(saved_Z(ia), saved_smd_alpha)
   enddo
   call cosmo_build_cavity_deriv(natoms, saved_Z, coor_bohr, saved_radii_scale, saved_avg_area_ang2, &
                                  d_tc, d_ta, d_tat, d_n, d_nslot, d_aos, d_dA, d_dC, Rvdw)
   deallocate(Rvdw)
else
   call cosmo_build_cavity_deriv(natoms, saved_Z, coor_bohr, saved_radii_scale, saved_avg_area_ang2, &
                                  d_tc, d_ta, d_tat, d_n, d_nslot, d_aos, d_dA, d_dC)
endif
nt = 0
do k = 1,d_n
   if (d_ta(k)*bohr2ang*bohr2ang >= 1.0d-4) nt = nt+1
enddo
allocate(tc(3,nt), ta(nt), Js(3,3*natoms,nt), Ja(3*natoms,nt))
Js = 0.0d0 ; Ja = 0.0d0
kk = 0
do k = 1,d_n
   if (d_ta(k)*bohr2ang*bohr2ang < 1.0d-4) cycle
   kk = kk+1
   tc(:,kk) = d_tc(:,k)
   ta(kk) = d_ta(k)
   do s = 1,d_nslot(k)
      ia2 = d_aos(s,k)
      do j = 1,3
         Js(:,3*(ia2-1)+j,kk) = Js(:,3*(ia2-1)+j,kk) + d_dC(:,j,s,k)
         Ja(3*(ia2-1)+j,kk) = Ja(3*(ia2-1)+j,kk) + d_dA(j,s,k)
      enddo
   enddo
enddo
end subroutine cosmo_hess_cavity

subroutine cosmo_hess_response(nConts, nb, D, dF)
use mod_integrals, only: cosmo_build_one_tess_matrix, esp_at_grid_batch, cosmo_fock_batch
implicit none
integer,intent(in) :: nConts, nb
real(8),intent(in) :: D(nConts,nConts,nb)
real(8),intent(out) :: dF(nConts,nConts,nb)
real(8),allocatable :: V(:,:), Bk(:,:)
integer :: k, b, info
real(8) :: lam

allocate(V(n_tess,nb))
if (store_mode) then
   call dgemm('T','N', n_tess, nb, nConts*nConts, 1.0d0, B_store, nConts*nConts, D, nConts*nConts, &
              0.0d0, V, n_tess)
else if (COSMO_DIRECT_BATCH) then
   block
      real(8), allocatable :: Ds(:,:)
      allocate(Ds(nConts,nConts))
      do b = 1,nb
         Ds = 0.5d0*(D(:,:,b) + transpose(D(:,:,b)))
         call esp_at_grid_batch(n_tess, tess_coor, Ds, V(:,b))
      enddo
   end block
else
   allocate(Bk(nConts,nConts))
   do k = 1,n_tess
      call cosmo_build_one_tess_matrix(tess_coor(:,k), Bk)
      do b = 1,nb
         V(k,b) = sum(D(:,:,b)*Bk)
      enddo
   enddo
   deallocate(Bk)
endif
do b = 1,nb
   call dsytrs('U', n_tess, 1, A_fac, n_tess, A_ipiv, V(:,b), n_tess, info)
enddo
do b = 1,nb
   lam = sum(V(:,b))/saved_sum_Ainv1
   V(:,b) = -f_eps*(V(:,b) - lam*saved_Ainv1)
enddo
if (store_mode) then
   call dgemm('N','N', nConts*nConts, nb, n_tess, 1.0d0, B_store, nConts*nConts, V, n_tess, &
              0.0d0, dF, nConts*nConts)
else if (COSMO_DIRECT_BATCH) then
   dF = 0.0d0
   call cosmo_fock_batch(n_tess, tess_coor, nb, V, nConts, dF)
else
   allocate(Bk(nConts,nConts))
   dF = 0.0d0
   do k = 1,n_tess
      call cosmo_build_one_tess_matrix(tess_coor(:,k), Bk)
      do b = 1,nb
         dF(:,:,b) = dF(:,:,b) + V(k,b)*Bk
      enddo
   enddo
   deallocate(Bk)
endif
deallocate(V)
end subroutine cosmo_hess_response

subroutine cosmo_hess_add(nConts, nb, D, Jbuf)
implicit none
integer,intent(in) :: nConts, nb
real(8),intent(in) :: D(nConts,nConts,nb)
real(8),intent(inout) :: Jbuf(nConts,nConts,nb)
real(8),allocatable :: dF(:,:,:)
allocate(dF(nConts,nConts,nb))
call cosmo_hess_response(nConts, nb, D, dF)
Jbuf = Jbuf + dF
deallocate(dF)
end subroutine cosmo_hess_add

subroutine cosmo_hess_add_spin(nConts, nb, Dtot, Vout)
implicit none
integer,intent(in) :: nConts, nb
real(8),intent(in) :: Dtot(nConts,nConts,nb)
real(8),intent(inout) :: Vout(nConts,nConts,nb,2)
real(8),allocatable :: dF(:,:,:)
allocate(dF(nConts,nConts,nb))
call cosmo_hess_response(nConts, nb, Dtot, dF)
Vout(:,:,:,1) = Vout(:,:,:,1) + dF
Vout(:,:,:,2) = Vout(:,:,:,2) + dF
deallocate(dF)
end subroutine cosmo_hess_add_spin

subroutine cosmo_hess_vind_c(nConts, nb, nocc, Co, Lf, VC)
implicit none
integer,intent(in) :: nConts, nb, nocc
real(8),intent(in) :: Co(nConts,nocc), Lf(nConts,nocc,nb)
real(8),intent(inout) :: VC(nConts,nocc,nb)
real(8),allocatable :: Dch(:,:,:), Fch(:,:,:)
integer :: c0, nch, b
allocate(Dch(nConts,nConts,16), Fch(nConts,nConts,16))
do c0 = 0, nb-1, 16
   nch = min(16, nb-c0)
   do b = 1,nch
      call dgemm('N','T',nConts,nConts,nocc,1.0d0,Lf(1,1,c0+b),nConts,Co,nConts,0.0d0,Dch(1,1,b),nConts)
      Dch(:,:,b) = Dch(:,:,b) + transpose(Dch(:,:,b))
   enddo
   call cosmo_hess_response(nConts, nch, Dch(:,:,1:nch), Fch(:,:,1:nch))
   do b = 1,nch
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Fch(1,1,b),nConts,Co,nConts,1.0d0,VC(1,1,c0+b),nConts)
   enddo
enddo
deallocate(Dch, Fch)
end subroutine cosmo_hess_vind_c

subroutine cosmo_hess_vind_spin_c(nConts, nb, noa, nob, Ca, Cb, La, Lb, VCa, VCb)
implicit none
integer,intent(in) :: nConts, nb, noa, nob
real(8),intent(in) :: Ca(nConts,noa), Cb(nConts,nob), La(nConts,noa,nb), Lb(nConts,nob,nb)
real(8),intent(inout) :: VCa(nConts,noa,nb), VCb(nConts,nob,nb)
real(8),allocatable :: Dt(:,:,:), Dtmp(:,:), Fch(:,:,:)
integer :: c0, nch, b
allocate(Dt(nConts,nConts,16), Dtmp(nConts,nConts), Fch(nConts,nConts,16))
do c0 = 0, nb-1, 16
   nch = min(16, nb-c0)
   do b = 1,nch
      call dgemm('N','T',nConts,nConts,noa,1.0d0,La(1,1,c0+b),nConts,Ca,nConts,0.0d0,Dt(1,1,b),nConts)
      call dgemm('N','T',nConts,nConts,nob,1.0d0,Lb(1,1,c0+b),nConts,Cb,nConts,0.0d0,Dtmp,nConts)
      Dt(:,:,b) = Dt(:,:,b) + Dtmp
      Dt(:,:,b) = Dt(:,:,b) + transpose(Dt(:,:,b))
   enddo
   call cosmo_hess_response(nConts, nch, Dt(:,:,1:nch), Fch(:,:,1:nch))
   do b = 1,nch
      call dgemm('N','N',nConts,noa,nConts,1.0d0,Fch(1,1,b),nConts,Ca,nConts,1.0d0,VCa(1,1,c0+b),nConts)
      call dgemm('N','N',nConts,nob,nConts,1.0d0,Fch(1,1,b),nConts,Cb,nConts,1.0d0,VCb(1,1,c0+b),nConts)
   enddo
enddo
deallocate(Dt, Dtmp, Fch)
end subroutine cosmo_hess_vind_spin_c

subroutine cosmo_hess_prepare(nConts, natoms, ao_atom, coor_bohr, Ptot, d2E)
use omp_lib, only: omp_get_wtime
use mod_integrals, only: cosmo_build_one_tess_matrix, cosmo_build_one_tess_matrix_shellderiv, cosmo_fock_batch, &
                         integrals_hessian_nuc_atom_raw, cosmo_hess_elec_batch
use mod_cosmo_scf, only: cosmo_scf_step
use mod_cosmo_force, only: cosmo_force_step
implicit none
integer,intent(in) :: nConts, natoms
integer,intent(in) :: ao_atom(nConts)
real(8),intent(in) :: coor_bohr(3,natoms), Ptot(nConts,nConts)
real(8),intent(out) :: d2E(3,natoms,3,natoms)

integer :: ntot, nt, k, l, x, y, i, j, a, b, A1, B1, mu, nu, r, info, i9, nt2
real(8),allocatable :: tc(:,:), ta(:), Js(:,:,:), Ja(:,:)
real(8),allocatable :: tc2(:,:), ta2(:), Jsp(:,:,:), Jap(:,:), Jsm(:,:,:), Jam(:,:)
real(8),allocatable :: gA(:,:), gs(:,:), G(:,:), w(:,:), dq(:,:), Us(:,:)
real(8),allocatable :: Dk(:,:,:), AA(:,:,:), AB(:,:,:), Fdum(:,:), Sym(:,:,:)
real(8),allocatable :: Es(:,:), Ea(:), Laa(:), pref(:), hk2(:,:), hAs(:,:), cdis(:,:)
real(8),allocatable :: Zk(:,:), Mkk(:,:), Tm(:,:), coor2(:,:)
real(8) :: gk(3,natoms), hk(3,natoms,3,natoms), hss(3,3), Z, d(3), dist, u(3), p, Enuc, lam, step
real(8) :: cc(3), q_k
character(len=32) :: envv
logical :: dbg
real(8) :: tw0, tw(6)
logical,parameter :: COSMO_HESS_BATCH = .true.
logical,parameter :: COSMO_HESS_CHECK = .false.
real(8),allocatable :: Cb(:,:,:), gshb(:,:,:), hrowb(:,:,:,:), Gdb(:,:,:,:), Wb(:,:,:)
real(8),allocatable :: Gnew(:,:), gAnew(:,:), Fxnew(:,:,:,:)
real(8) :: hkn(3,3)

ntot = 3*natoms
call get_environment_variable("ENGINE_COSMO_HESS_DEBUG", envv)
dbg = (trim(envv) .eq. "1")
step = 2.0d-4
call get_environment_variable("ENGINE_COSMO_HESS_STEP", envv)
if (len_trim(envv) > 0) read(envv,*) step

allocate(Fdum(nConts,nConts))
call cosmo_scf_step(nConts, Ptot, Fdum, Enuc, 0)
deallocate(Fdum)

call cosmo_hess_cavity(natoms, coor_bohr, nt, tc, ta, Js, Ja)
if (nt /= n_tess) then
   print '("COSMO Hessian error: tessera count ",I0," differs from the SCF cavity ",I0)', nt, n_tess
   stop 1
endif
if (maxval(abs(tc - tess_coor)) > 1.0d-6) then
   print '("COSMO Hessian error: cavity differs from the SCF one, max |ds| = ",ES10.3)', maxval(abs(tc - tess_coor))
   stop 1
endif

tw0 = omp_get_wtime()
allocate(gA(ntot,nt), gs(3,nt), G(ntot,ntot), hk2(ntot,ntot), hAs(ntot,3))
allocate(Dk(nConts,nConts,3), AA(nConts,nConts,9), AB(nConts,nConts,9))
G = 0.0d0

if (COSMO_HESS_BATCH) then
   allocate(Cb(nt,3,3+ntot), gshb(3,natoms,nt), hrowb(3,3,natoms,nt), Gdb(3,natoms,3,natoms), &
            Wb(nConts,nConts,3+ntot))
   Cb = 0.0d0
   do k = 1,nt
      do i = 1,3
         Cb(k,i,i) = q_last(k)
         Cb(k,i,4:3+ntot) = q_last(k)*Js(i,:,k)
      enddo
   enddo
   call cosmo_hess_elec_batch(nt, tc, q_last(1:nt), Ptot, ao_atom, natoms, 3+ntot, Cb, &
                              gshb, hrowb, Gdb, Wb)
   G = reshape(Gdb, (/ntot,ntot/))
   do k = 1,nt
      q_k = q_last(k)
      gk = 2.0d0*gshb(:,:,k)
      hAs = 0.0d0
      do A1 = 1,natoms
         Z = real(saved_Z(A1)-saved_ecpCoreElec(A1),8)
         d = tc(:,k) - coor_bohr(:,A1)
         dist = sqrt(sum(d*d))
         gk(:,A1) = gk(:,A1) + Z*d/dist**3
         do i = 1,3
            do j = 1,3
               hkn(i,j) = Z*(3.0d0*d(i)*d(j)/dist**5)
            enddo
            hkn(i,i) = hkn(i,i) - Z/dist**3
         enddo
         do i = 1,3
            do j = 1,3
               G(3*(A1-1)+i,3*(A1-1)+j) = G(3*(A1-1)+i,3*(A1-1)+j) + q_k*hkn(i,j)
               hAs(3*(A1-1)+i,j) = -(hrowb(i,j,A1,k) + hkn(i,j))
            enddo
         enddo
      enddo
      gA(:,k) = reshape(gk,(/ntot/))
      do i = 1,3
         gs(i,k) = -sum(gk(i,:))
      enddo
      do i = 1,3
         do j = 1,3
            hss(i,j) = 0.0d0
            do A1 = 1,natoms
               hss(i,j) = hss(i,j) - hAs(3*(A1-1)+i,j)
            enddo
         enddo
      enddo
      G = G + q_k*matmul(hAs, Js(:,:,k))
      G = G + q_k*transpose(matmul(hAs, Js(:,:,k)))
      G = G + q_k*matmul(transpose(Js(:,:,k)), matmul(hss, Js(:,:,k)))
   enddo
   deallocate(Cb, gshb, hrowb, Gdb)
   if (COSMO_HESS_CHECK) then
      allocate(Gnew(ntot,ntot), gAnew(ntot,nt))
      Gnew = G ; gAnew = gA
      G = 0.0d0
   endif
endif

do k = 1,nt
   if (COSMO_HESS_BATCH .and. .not. COSMO_HESS_CHECK) exit
   q_k = q_last(k)
   call cosmo_build_one_tess_matrix_shellderiv(tc(:,k), Dk)
   call integrals_hessian_nuc_atom_raw(tc(:,k), nConts, AA, AB)
   gk = 0.0d0 ; hk = 0.0d0
   do A1 = 1,natoms
      Z = real(saved_Z(A1)-saved_ecpCoreElec(A1),8)
      d = tc(:,k) - coor_bohr(:,A1)
      dist = sqrt(sum(d*d))
      gk(:,A1) = gk(:,A1) + Z*d/dist**3
      do i = 1,3
         do j = 1,3
            hk(i,A1,j,A1) = hk(i,A1,j,A1) + Z*(3.0d0*d(i)*d(j)/dist**5)
         enddo
         hk(i,A1,i,A1) = hk(i,A1,i,A1) - Z/dist**3
      enddo
   enddo
   do nu = 1,nConts
      B1 = ao_atom(nu)
      do mu = 1,nConts
         p = Ptot(mu,nu)
         if (p == 0.0d0) cycle
         A1 = ao_atom(mu)
         gk(:,A1) = gk(:,A1) - p*Dk(mu,nu,:)
         gk(:,B1) = gk(:,B1) - p*Dk(nu,mu,:)
         do a = 1,3
            do b = 1,3
               i9 = 3*(a-1)+b
               hk(a,A1,b,A1) = hk(a,A1,b,A1) - p*AA(mu,nu,i9)
               hk(a,B1,b,B1) = hk(a,B1,b,B1) - p*AA(nu,mu,i9)
               hk(a,A1,b,B1) = hk(a,A1,b,B1) - p*AB(mu,nu,i9)
               hk(a,B1,b,A1) = hk(a,B1,b,A1) - p*AB(mu,nu,3*(b-1)+a)
            enddo
         enddo
      enddo
   enddo
   gA(:,k) = reshape(gk,(/ntot/))
   do i = 1,3
      gs(i,k) = -sum(gk(i,:))
   enddo
   hk2 = reshape(hk,(/ntot,ntot/))
   do x = 1,ntot
      do i = 1,3
         hAs(x,i) = 0.0d0
         do B1 = 1,natoms
            hAs(x,i) = hAs(x,i) - hk2(x,3*(B1-1)+i)
         enddo
      enddo
   enddo
   do i = 1,3
      do j = 1,3
         hss(i,j) = 0.0d0
         do A1 = 1,natoms
            hss(i,j) = hss(i,j) - hAs(3*(A1-1)+i,j)
         enddo
      enddo
   enddo
   G = G + q_k*hk2
   G = G + q_k*matmul(hAs, Js(:,:,k))
   G = G + q_k*transpose(matmul(hAs, Js(:,:,k)))
   G = G + q_k*matmul(transpose(Js(:,:,k)), matmul(hss, Js(:,:,k)))
enddo
deallocate(Dk, AA, AB, hk2, hAs)
if (COSMO_HESS_BATCH .and. COSMO_HESS_CHECK) then
   print '(A,2ES12.4)', 'COSMO HESS CHECK G  max|new-old|, max|old| ', maxval(abs(Gnew-G)), maxval(abs(G))
   print '(A,2ES12.4)', 'COSMO HESS CHECK gA max|new-old|, max|old| ', maxval(abs(gAnew-gA)), maxval(abs(gA))
   deallocate(Gnew, gAnew)
endif
tw(1) = omp_get_wtime()-tw0; tw0 = omp_get_wtime()

allocate(Es(3,nt), Ea(nt), Laa(nt), pref(nt), w(nt,ntot), Zk(3,ntot), Mkk(3,3), Tm(3,3))
do k = 1,nt
   pref(k) = -0.5d0*1.07d0*sqrt(4.0d0*PI_cosmo)*ta(k)**(-1.5d0)
   Ea(k) = q_last(k)**2/(2.0d0*f_eps)*pref(k)
   Laa(k) = q_last(k)**2/(2.0d0*f_eps)*0.75d0*1.07d0*sqrt(4.0d0*PI_cosmo)*ta(k)**(-2.5d0)
   Es(:,k) = q_last(k)*gs(:,k)
   do x = 1,ntot
      w(k,x) = gA(x,k) + dot_product(gs(:,k), Js(:,x,k)) + q_last(k)*pref(k)*Ja(x,k)/f_eps
   enddo
enddo
!$omp parallel do private(k,l,d,dist,u,Zk,Mkk,Tm,i,j,x) schedule(dynamic) reduction(+:G)
do k = 1,nt
   Zk = 0.0d0 ; Mkk = 0.0d0
   do l = 1,nt
      if (l == k) cycle
      d = tc(:,k) - tc(:,l)
      dist = sqrt(sum(d*d))
      u = d/dist**3
      Es(:,k) = Es(:,k) - q_last(k)*q_last(l)*u/f_eps
      do x = 1,ntot
         w(k,x) = w(k,x) + (q_last(l)/f_eps)*dot_product(u, Js(:,x,l) - Js(:,x,k))
      enddo
      do i = 1,3
         do j = 1,3
            Tm(i,j) = 3.0d0*d(i)*d(j)/dist**5
         enddo
         Tm(i,i) = Tm(i,i) - 1.0d0/dist**3
      enddo
      Mkk = Mkk + (q_last(k)*q_last(l)/f_eps)*Tm
      Zk = Zk - (q_last(k)*q_last(l)/f_eps)*matmul(Tm, Js(:,:,l))
   enddo
   Zk = Zk + matmul(Mkk, Js(:,:,k))
   G = G + matmul(transpose(Js(:,:,k)), Zk)
   do y = 1,ntot
      do x = 1,ntot
         G(x,y) = G(x,y) + Laa(k)*Ja(x,k)*Ja(y,k)
      enddo
   enddo
enddo
!$omp end parallel do

tw(2) = omp_get_wtime()-tw0; tw0 = omp_get_wtime()
if (dbg) then
   block
      use mod_cosmo_force, only: cosmo_force_step
      real(8) :: gmine(ntot), fref(3,natoms), cang(3,natoms)
      gmine = 0.0d0
      do k = 1,nt
         gmine = gmine + q_last(k)*gA(:,k)
         do x = 1,ntot
            gmine(x) = gmine(x) + dot_product(Es(:,k), Js(:,x,k)) + Ea(k)*Ja(x,k)
         enddo
      enddo
      cang = coor_bohr*bohr2ang
      call cosmo_force_step(natoms, cang, ao_atom, nConts, Ptot, fref)
      print '("COSMO HESS DEBUG: gradient max|mine - cosmo_force_step| = ",ES10.3,"  max|g| = ",ES10.3)', &
            maxval(abs(gmine - reshape(fref,(/ntot/)))), maxval(abs(gmine))
   end block
endif

allocate(dq(nt,ntot))
dq = w
do x = 1,ntot
   call dsytrs('U', nt, 1, A_fac, nt, A_ipiv, dq(:,x), nt, info)
enddo
do x = 1,ntot
   lam = sum(dq(:,x))/saved_sum_Ainv1
   dq(:,x) = -f_eps*(dq(:,x) - lam*saved_Ainv1)
enddo
G = G + matmul(transpose(w), dq)

tw(3) = omp_get_wtime()-tw0; tw0 = omp_get_wtime()
!$omp parallel do private(y,A1,r,coor2,nt2,tc2,ta2,Jsp,Jap,Jsm,Jam,k,x) schedule(dynamic)
do y = 1,ntot
   allocate(coor2(3,natoms))
   A1 = (y-1)/3+1 ; r = mod(y-1,3)+1
   coor2 = coor_bohr
   coor2(r,A1) = coor_bohr(r,A1) + step
   call cosmo_hess_cavity(natoms, coor2, nt2, tc2, ta2, Jsp, Jap)
   if (nt2 /= nt) then
      print '("COSMO Hessian error: tessera count changed under displacement ",I0," -> ",I0)', nt, nt2
      print *, 'try a smaller ENGINE_COSMO_HESS_STEP'
      stop 1
   endif
   coor2(r,A1) = coor_bohr(r,A1) - step
   call cosmo_hess_cavity(natoms, coor2, nt2, tc2, ta2, Jsm, Jam)
   if (nt2 /= nt) then
      print '("COSMO Hessian error: tessera count changed under displacement ",I0," -> ",I0)', nt, nt2
      print *, 'try a smaller ENGINE_COSMO_HESS_STEP'
      stop 1
   endif
   do k = 1,nt
      do x = 1,ntot
         G(x,y) = G(x,y) + (dot_product(Es(:,k), Jsp(:,x,k)-Jsm(:,x,k)) + Ea(k)*(Jap(x,k)-Jam(x,k)))/(2.0d0*step)
      enddo
   enddo
   deallocate(coor2)
enddo
!$omp end parallel do

tw(4) = omp_get_wtime()-tw0; tw0 = omp_get_wtime()
G = 0.5d0*(G + transpose(G))
d2E = reshape(G, (/3,natoms,3,natoms/))

if (allocated(cosmo_hess_Fx)) deallocate(cosmo_hess_Fx)
allocate(cosmo_hess_Fx(nConts,nConts,3,natoms))
cosmo_hess_Fx = 0.0d0
if (COSMO_HESS_BATCH) then
   do nu = 1,nConts
      B1 = ao_atom(nu)
      do mu = 1,nConts
         A1 = ao_atom(mu)
         cosmo_hess_Fx(mu,nu,:,A1) = cosmo_hess_Fx(mu,nu,:,A1) + Wb(mu,nu,1:3)
         cosmo_hess_Fx(mu,nu,:,B1) = cosmo_hess_Fx(mu,nu,:,B1) + Wb(nu,mu,1:3)
      enddo
   enddo
   do x = 1,ntot
      A1 = (x-1)/3+1 ; r = mod(x-1,3)+1
      cosmo_hess_Fx(:,:,r,A1) = cosmo_hess_Fx(:,:,r,A1) - Wb(:,:,3+x) - transpose(Wb(:,:,3+x))
   enddo
   deallocate(Wb)
   if (COSMO_HESS_CHECK) then
      allocate(Fxnew(nConts,nConts,3,natoms))
      Fxnew = cosmo_hess_Fx
      cosmo_hess_Fx = 0.0d0
   endif
endif
allocate(Dk(nConts,nConts,3), Sym(nConts,nConts,3))
do k = 1,nt
   if (COSMO_HESS_BATCH .and. .not. COSMO_HESS_CHECK) exit
   q_k = q_last(k)
   call cosmo_build_one_tess_matrix_shellderiv(tc(:,k), Dk)
   do nu = 1,nConts
      B1 = ao_atom(nu)
      do mu = 1,nConts
         A1 = ao_atom(mu)
         cosmo_hess_Fx(mu,nu,:,A1) = cosmo_hess_Fx(mu,nu,:,A1) - q_k*Dk(mu,nu,:)
         cosmo_hess_Fx(mu,nu,:,B1) = cosmo_hess_Fx(mu,nu,:,B1) - q_k*Dk(nu,mu,:)
      enddo
   enddo
   do i = 1,3
      Sym(:,:,i) = Dk(:,:,i) + transpose(Dk(:,:,i))
   enddo
   do x = 1,ntot
      cc = q_k*Js(:,x,k)
      if (all(cc == 0.0d0)) cycle
      A1 = (x-1)/3+1 ; r = mod(x-1,3)+1
      do i = 1,3
         cosmo_hess_Fx(:,:,r,A1) = cosmo_hess_Fx(:,:,r,A1) + cc(i)*Sym(:,:,i)
      enddo
   enddo
enddo
deallocate(Dk, Sym)
if (COSMO_HESS_BATCH .and. COSMO_HESS_CHECK) then
   print '(A,2ES12.4)', 'COSMO HESS CHECK Fx max|new-old|, max|old| ', maxval(abs(Fxnew-cosmo_hess_Fx)), &
         maxval(abs(cosmo_hess_Fx))
   cosmo_hess_Fx = Fxnew
   deallocate(Fxnew)
endif
tw(5) = omp_get_wtime()-tw0; tw0 = omp_get_wtime()
if (store_mode) then
   call dgemm('N','N', nConts*nConts, ntot, nt, 1.0d0, B_store, nConts*nConts, dq, nt, &
              1.0d0, cosmo_hess_Fx, nConts*nConts)
else if (COSMO_DIRECT_BATCH) then
   call cosmo_fock_batch(nt, tc, ntot, dq, nConts, cosmo_hess_Fx)
else
   allocate(Dk(nConts,nConts,1))
   do k = 1,nt
      call cosmo_build_one_tess_matrix(tc(:,k), Dk(:,:,1))
      do x = 1,ntot
         A1 = (x-1)/3+1 ; r = mod(x-1,3)+1
         cosmo_hess_Fx(:,:,r,A1) = cosmo_hess_Fx(:,:,r,A1) + dq(k,x)*Dk(:,:,1)
      enddo
   enddo
   deallocate(Dk)
endif
tw(6) = omp_get_wtime()-tw0
if (dbg) print '(A,6F9.3)', 'COSMO HESS TIMING k-loop/pairs/dq/cavFD/Fx/Fxdq ', tw
cosmo_hess_active = .true.
end subroutine cosmo_hess_prepare

end module mod_cosmo_hessian
