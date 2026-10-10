! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_soscf: EXPERIMENTAL, opt-in (ENGINE_SOSCF=1) second-order-SCF

module mod_soscf
implicit none
private
public :: soscf_reset, soscf_step, newton_cg_step

integer,parameter :: SOSCF_MAX = 10
type :: hist_t
   integer :: nvec = 0
   real(8),allocatable :: s_hist(:,:), y_hist(:,:), rho_hist(:)
   real(8),allocatable :: g_prev(:), s_prev(:)
   logical :: have_prev = .false.
end type hist_t
type(hist_t),save :: H(2)

contains

subroutine soscf_reset()
implicit none
integer :: k
do k = 1,2
   H(k)%nvec = 0
   H(k)%have_prev = .false.
   if (allocated(H(k)%s_hist)) deallocate(H(k)%s_hist)
   if (allocated(H(k)%y_hist)) deallocate(H(k)%y_hist)
   if (allocated(H(k)%rho_hist)) deallocate(H(k)%rho_hist)
   if (allocated(H(k)%g_prev)) deallocate(H(k)%g_prev)
   if (allocated(H(k)%s_prev)) deallocate(H(k)%s_prev)
enddo
end subroutine soscf_reset

subroutine soscf_step(nconts, nocc, Fa, S, C_a, eLev_a, grad_max, took_step, spin)
implicit none
integer,intent(in) :: nconts, nocc
real(8),intent(in) :: Fa(nconts,nconts), S(nconts,nconts)
real(8),intent(in) :: grad_max
real(8),intent(inout) :: C_a(nconts,nconts), eLev_a(nconts)
logical,intent(out) :: took_step
integer,intent(in) :: spin

integer :: nvir, nd, i, a, k
real(8),allocatable :: Tmp(:,:), Fmo_ov(:,:)
real(8),allocatable :: gflat(:), q(:), d(:), h0inv(:)
real(8),allocatable :: alpha_hist(:)
real(8),allocatable :: ynew(:)
real(8) :: sy, beta_k
real(8),allocatable :: kappa(:,:)

nvir = nconts - nocc
nd = nvir*nocc
took_step = .false.
if (nd .le. 0) return

allocate(Tmp(nconts,nocc))
call DGEMM('N','N',nconts,nocc,nconts,1.0d0,Fa,nconts,C_a(1,1),nconts,0.0d0,Tmp,nconts)
allocate(Fmo_ov(nvir,nocc))
call DGEMM('T','N',nvir,nocc,nconts,1.0d0,C_a(1,nocc+1),nconts,Tmp,nconts,0.0d0,Fmo_ov,nvir)
deallocate(Tmp)
allocate(gflat(nd))
gflat = reshape(Fmo_ov,(/nd/))
deallocate(Fmo_ov)

block
   real(8) :: grms
   grms = sqrt(sum(gflat**2)/real(nd,8))
   if (grms .gt. grad_max) then
      deallocate(gflat)
      return
   endif
end block

allocate(h0inv(nd))
do i = 1,nocc
   do a = 1,nvir
      h0inv((i-1)*nvir+a) = 1.0d0/max(eLev_a(nocc+a)-eLev_a(i), 1.0d-3)
   enddo
enddo

if (H(spin)%have_prev) then
   allocate(ynew(nd))
   ynew = gflat - H(spin)%g_prev
   sy = dot_product(H(spin)%s_prev, ynew)
   if (sy .gt. 1.0d-12) then
      if (H(spin)%nvec .eq. SOSCF_MAX) then
         H(spin)%s_hist(:,1:SOSCF_MAX-1) = H(spin)%s_hist(:,2:SOSCF_MAX)
         H(spin)%y_hist(:,1:SOSCF_MAX-1) = H(spin)%y_hist(:,2:SOSCF_MAX)
         H(spin)%rho_hist(1:SOSCF_MAX-1) = H(spin)%rho_hist(2:SOSCF_MAX)
         H(spin)%nvec = SOSCF_MAX - 1
      endif
      H(spin)%nvec = H(spin)%nvec + 1
      H(spin)%s_hist(:,H(spin)%nvec) = H(spin)%s_prev
      H(spin)%y_hist(:,H(spin)%nvec) = ynew
      H(spin)%rho_hist(H(spin)%nvec) = 1.0d0/sy
   endif
   deallocate(ynew)
endif

allocate(q(nd)); q = gflat
allocate(alpha_hist(max(H(spin)%nvec,1)))
do k = H(spin)%nvec,1,-1
   alpha_hist(k) = H(spin)%rho_hist(k) * dot_product(H(spin)%s_hist(:,k), q)
   q = q - alpha_hist(k)*H(spin)%y_hist(:,k)
enddo
allocate(d(nd))
d = h0inv * q
do k = 1,H(spin)%nvec
   beta_k = H(spin)%rho_hist(k) * dot_product(H(spin)%y_hist(:,k), d)
   d = d + H(spin)%s_hist(:,k)*(alpha_hist(k)-beta_k)
enddo
deallocate(q,alpha_hist,h0inv)

allocate(kappa(nvir,nocc))
kappa = -reshape(d,(/nvir,nocc/))
deallocate(d)

block
   real(8),parameter :: SOSCF_TRUST = 1.0d-1
   real(8) :: knorm
   knorm = sqrt(sum(kappa**2))
   if (knorm .gt. SOSCF_TRUST) kappa = kappa*(SOSCF_TRUST/knorm)
end block

if (.not. allocated(H(spin)%s_hist)) then
   allocate(H(spin)%s_hist(nd,SOSCF_MAX), H(spin)%y_hist(nd,SOSCF_MAX), H(spin)%rho_hist(SOSCF_MAX))
endif
if (allocated(H(spin)%g_prev)) deallocate(H(spin)%g_prev)
if (allocated(H(spin)%s_prev)) deallocate(H(spin)%s_prev)
allocate(H(spin)%g_prev(nd), H(spin)%s_prev(nd))
H(spin)%g_prev = gflat
H(spin)%s_prev = reshape(kappa,(/nd/))
H(spin)%have_prev = .true.
deallocate(gflat)

call soscf_rotate(nconts, nocc, kappa, Fa, C_a, eLev_a)
deallocate(kappa)
took_step = .true.

end subroutine soscf_step

subroutine soscf_rotate(nconts, nocc, kappa, Fa, C_a, eLev_a)
implicit none
integer,intent(in) :: nconts, nocc
real(8),intent(in) :: kappa(nconts-nocc,nocc), Fa(nconts,nconts)
real(8),intent(inout) :: C_a(nconts,nconts), eLev_a(nconts)

integer :: nvir, r, p, info, lwork
real(8),allocatable :: Cnew(:,:)
real(8),allocatable :: Umat(:,:), VTmat(:,:), Sigma(:), kcopy(:,:), work(:)
real(8),allocatable :: Dblk_occ(:,:), Dblk_vir(:,:), Tmp(:,:)
real(8),allocatable :: Sinv_occ(:,:), Sinv_vir(:,:)
real(8),allocatable :: FC(:,:)
real(8) :: wq(1)

nvir = nconts - nocc
r = min(nvir,nocc)

allocate(Cnew(nconts,nconts))
Cnew = C_a
call DGEMM('N','N',nconts,nocc,nvir,1.0d0,C_a(1,nocc+1),nconts,kappa,nvir,1.0d0,Cnew(1,1),nconts)
call DGEMM('N','T',nconts,nvir,nocc,-1.0d0,C_a(1,1),nconts,kappa,nvir,1.0d0,Cnew(1,nocc+1),nconts)

allocate(Umat(nvir,nvir), VTmat(nocc,nocc), Sigma(r), kcopy(nvir,nocc))
kcopy = kappa
call DGESVD('A','A',nvir,nocc,kcopy,nvir,Sigma,Umat,nvir,VTmat,nocc,wq,-1,info)
lwork = int(wq(1))
allocate(work(lwork))
call DGESVD('A','A',nvir,nocc,kcopy,nvir,Sigma,Umat,nvir,VTmat,nocc,work,lwork,info)
deallocate(work,kcopy)

allocate(Dblk_vir(nvir,nvir)); Dblk_vir = 0.0d0
do p = 1,nvir
   if (p .le. r) then
      Dblk_vir(p,p) = 1.0d0/sqrt(1.0d0+Sigma(p)**2)
   else
      Dblk_vir(p,p) = 1.0d0
   endif
enddo
allocate(Tmp(nvir,nvir))
call DGEMM('N','N',nvir,nvir,nvir,1.0d0,Umat,nvir,Dblk_vir,nvir,0.0d0,Tmp,nvir)
allocate(Sinv_vir(nvir,nvir))
call DGEMM('N','T',nvir,nvir,nvir,1.0d0,Tmp,nvir,Umat,nvir,0.0d0,Sinv_vir,nvir)
deallocate(Tmp,Dblk_vir,Umat)

allocate(Dblk_occ(nocc,nocc)); Dblk_occ = 0.0d0
do p = 1,nocc
   if (p .le. r) then
      Dblk_occ(p,p) = 1.0d0/sqrt(1.0d0+Sigma(p)**2)
   else
      Dblk_occ(p,p) = 1.0d0
   endif
enddo
allocate(Tmp(nocc,nocc))
call DGEMM('T','N',nocc,nocc,nocc,1.0d0,VTmat,nocc,Dblk_occ,nocc,0.0d0,Tmp,nocc)
allocate(Sinv_occ(nocc,nocc))
call DGEMM('N','N',nocc,nocc,nocc,1.0d0,Tmp,nocc,VTmat,nocc,0.0d0,Sinv_occ,nocc)
deallocate(Tmp,Dblk_occ,VTmat,Sigma)

call DGEMM('N','N',nconts,nocc,nocc,1.0d0,Cnew(1,1),nconts,Sinv_occ,nocc,0.0d0,C_a(1,1),nconts)
call DGEMM('N','N',nconts,nvir,nvir,1.0d0,Cnew(1,nocc+1),nconts,Sinv_vir,nvir,0.0d0,C_a(1,nocc+1),nconts)
deallocate(Cnew,Sinv_occ,Sinv_vir)

allocate(FC(nconts,nconts))
call DGEMM('N','N',nconts,nconts,nconts,1.0d0,Fa,nconts,C_a,nconts,0.0d0,FC,nconts)
do p = 1,nconts
   eLev_a(p) = dot_product(C_a(:,p),FC(:,p))
enddo
deallocate(FC)

end subroutine soscf_rotate

subroutine newton_cg_step(C_a, C_b, eLev_a, eLev_b, na, nb, Prms, grad_max, trust, ncgmax, took_step, stab_mode, lam_out)
use MOL_info, only: Pa, Pb, Fa, Fb, nconts, Multi
use mod_exchange, only: cosx_enabled
use mod_integrals, only: cosx_set_response_mode, cosx_incr_stash
implicit none
real(8),intent(inout) :: C_a(:,:), C_b(:,:), eLev_a(:), eLev_b(:)
integer,intent(in) :: na, nb, ncgmax
real(8),intent(in) :: Prms, grad_max, trust
logical,intent(out) :: took_step
integer,intent(in),optional :: stab_mode
real(8),intent(out),optional :: lam_out

integer :: nsp, s, no(2), nv(2), nd(2), ndt, off(2), it, a, i
real(8),allocatable :: F0(:,:,:), P0(:,:,:), Cs(:,:,:), Foo(:,:,:), Fvv(:,:,:)
real(8),allocatable :: F0r(:,:,:)
logical :: have_F0r
logical :: an_done, an_ok, an_check
integer :: an_ncheck
real(8),allocatable :: g(:), x(:), r(:), z(:), p(:), Ap(:), Mdiag(:)
real(8) :: gn, rz, rz_new, alpha, beta, pAp, xn, rn
real(8),allocatable :: Tmp(:,:)
real(8),allocatable :: kap(:,:)
real(8) :: cgtol, stab_lam
integer :: stab_nit
character(len=32) :: cgenv

cgtol = 0.1d0
call get_environment_variable("ENGINE_NEWTON_CGTOL", cgenv)
if (len_trim(cgenv) .gt. 0) read(cgenv,*) cgtol
took_step = .false.
have_F0r = .false.
an_done = .false.
nsp = 2
if (multi .eq. 1) nsp = 1
no(1) = na; no(2) = nb
nv = nconts - no
nd = nv*no
if (nd(1) .le. 0) return
if (nsp .eq. 2 .and. nd(2) .le. 0) nsp = 1
if (nsp .eq. 1 .and. multi .ne. 1) return
ndt = nd(1); off(1) = 0; off(2) = nd(1)
if (nsp .eq. 2) ndt = nd(1) + nd(2)

allocate(F0(nconts,nconts,2), P0(nconts,nconts,2), Cs(nconts,nconts,2))
F0(:,:,1) = Fa; F0(:,:,2) = Fb
Cs(:,:,1) = C_a; Cs(:,:,2) = C_b
if (nsp .eq. 1) then
   F0(:,:,2) = Fa; Cs(:,:,2) = C_a
endif
do s = 1,2
   call DGEMM('N','T',nconts,nconts,no(s),1.0d0,Cs(1,1,s),nconts,Cs(1,1,s),nconts,0.0d0,P0(1,1,s),nconts)
enddo

allocate(Foo(nconts,nconts,2), Fvv(nconts,nconts,2), g(ndt), Mdiag(ndt))
do s = 1,nsp
   allocate(Tmp(nconts,nconts))
   call DGEMM('N','N',nconts,nconts,nconts,1.0d0,F0(1,1,s),nconts,Cs(1,1,s),nconts,0.0d0,Tmp,nconts)
   call DGEMM('T','N',no(s),no(s),nconts,1.0d0,Cs(1,1,s),nconts,Tmp(1,1),nconts,0.0d0,Foo(1,1,s),nconts)
   call DGEMM('T','N',nv(s),nv(s),nconts,1.0d0,Cs(1,no(s)+1,s),nconts,Tmp(1,no(s)+1),nconts,0.0d0,Fvv(1,1,s),nconts)
   allocate(kap(nv(s),no(s)))
   call DGEMM('T','N',nv(s),no(s),nconts,1.0d0,Cs(1,no(s)+1,s),nconts,Tmp(1,1),nconts,0.0d0,kap,nv(s))
   g(off(s)+1:off(s)+nd(s)) = reshape(kap,(/nd(s)/))
   deallocate(kap, Tmp)
   do i = 1,no(s)
      do a = 1,nv(s)
         if (s .eq. 1) then
            Mdiag(off(s)+(i-1)*nv(s)+a) = max(eLev_a(no(s)+a)-eLev_a(i), 0.1d0)
         else
            Mdiag(off(s)+(i-1)*nv(s)+a) = max(eLev_b(no(s)+a)-eLev_b(i), 0.1d0)
         endif
      enddo
   enddo
enddo

if (present(stab_mode)) then
   allocate(x(ndt), r(ndt))
   call stab_lanczos(stab_lam, x, stab_nit)
   if (present(lam_out)) lam_out = stab_lam
   print '(A,ES11.3,A,I3,A)', ' STABILITY: lowest orbital-Hessian eigenvalue =', stab_lam, &
      '  (', stab_nit, ' Lanczos steps)'
   if (stab_lam .lt. -1.0d-3) then
      print '(A)', ' STABILITY: wavefunction is UNSTABLE (negative eigenvalue)'
   else
      print '(A)', ' STABILITY: wavefunction is stable'
   endif
   if (stab_mode .ge. 2 .and. stab_lam .lt. -1.0d-3) then
      x = x*(trust/max(1.0d-12,sqrt(dot_product(x,x))))
      do s = 1,nsp
         allocate(kap(nv(s),no(s)))
         kap = reshape(x(off(s)+1:off(s)+nd(s)),(/nv(s),no(s)/))
         if (s .eq. 1) then
            call soscf_rotate(nconts, no(s), kap, F0(1,1,1), C_a, eLev_a)
         else
            call soscf_rotate(nconts, no(s), kap, F0(1,1,2), C_b, eLev_b)
         endif
         deallocate(kap)
      enddo
      if (nsp .eq. 1 .and. multi .eq. 1) then
         C_b = C_a; eLev_b = eLev_a
      endif
      took_step = .true.
   endif
   deallocate(F0,P0,Cs,Foo,Fvv,g,Mdiag,x,r)
   return
endif

gn = sqrt(sum(g**2)/real(ndt,8))
if (gn .gt. grad_max .or. gn .lt. 1.0d-9) then
   deallocate(F0,P0,Cs,Foo,Fvv,g,Mdiag)
   return
endif

allocate(x(ndt), r(ndt), z(ndt), p(ndt), Ap(ndt))
x = 0.0d0
r = -g
z = r/Mdiag
p = z
rz = dot_product(r,z)
rn = sqrt(dot_product(r,r))
do it = 1,ncgmax
   call apply_A(p, Ap)
   pAp = dot_product(p,Ap)
   if (pAp .le. 1.0d-14*dot_product(p,p)) then
      if (it .eq. 1) x = x + p*(0.1d0/ max(1.0d-12,sqrt(dot_product(p,p))))
      exit
   endif
   alpha = rz/pAp
   x = x + alpha*p
   r = r - alpha*Ap
   if (sqrt(dot_product(r,r)) .lt. cgtol*rn) exit
   z = r/Mdiag
   rz_new = dot_product(r,z)
   beta = rz_new/rz
   rz = rz_new
   p = z + beta*p
enddo

xn = sqrt(dot_product(x,x))
print '(A,ES10.3,A,I3,A,ES10.3,A,ES10.3,A,ES10.3)', ' NEWTON gn=',gn,' ncg=',it, &
   ' res/rn=',sqrt(dot_product(r,r))/rn,' |x|=',xn,' -g.x=',-dot_product(g,x)
if (xn .gt. trust) x = x*(trust/xn)

do s = 1,nsp
   allocate(kap(nv(s),no(s)))
   kap = reshape(x(off(s)+1:off(s)+nd(s)),(/nv(s),no(s)/))
   if (s .eq. 1) then
      call soscf_rotate(nconts, no(s), kap, F0(1,1,1), C_a, eLev_a)
   else
      call soscf_rotate(nconts, no(s), kap, F0(1,1,2), C_b, eLev_b)
   endif
   deallocate(kap)
enddo
if (nsp .eq. 1 .and. multi .eq. 1) then
   C_b = C_a; eLev_b = eLev_a
endif
took_step = .true.
deallocate(F0,P0,Cs,Foo,Fvv,g,Mdiag,x,r,z,p,Ap)

contains

subroutine stab_lanczos(lam, vec, nit)
real(8),intent(out) :: lam, vec(ndt)
integer,intent(out) :: nit
integer,parameter :: mmax = 40
real(8),allocatable :: Q(:,:), al(:), be(:), w(:), Tm(:,:), ev(:), wk(:)
integer :: j, k, i2, info2, rep
real(8) :: nrm, resid
allocate(Q(ndt,mmax+1), al(mmax), be(mmax), w(ndt))
do i2 = 1,ndt
   Q(i2,1) = sin(1.234d0*i2 + 0.37d0*mod(i2,7)) + 0.3d0
enddo
Q(:,1) = Q(:,1)/sqrt(dot_product(Q(:,1),Q(:,1)))
lam = 0.0d0
nit = 0
do j = 1,mmax
   call apply_A(Q(:,j), w)
   al(j) = dot_product(Q(:,j),w)
   w = w - al(j)*Q(:,j)
   if (j .gt. 1) w = w - be(j-1)*Q(:,j-1)
   do rep = 1,2
      do k = 1,j
         w = w - dot_product(Q(:,k),w)*Q(:,k)
      enddo
   enddo
   nrm = sqrt(dot_product(w,w))
   be(j) = nrm
   nit = j
   if (j .ge. 4 .and. (mod(j,3) .eq. 0 .or. nrm .lt. 1.0d-8 .or. j .eq. mmax)) then
      allocate(Tm(j,j), ev(j), wk(3*j))
      Tm = 0.0d0
      do k = 1,j
         Tm(k,k) = al(k)
         if (k .lt. j) then
            Tm(k,k+1) = be(k); Tm(k+1,k) = be(k)
         endif
      enddo
      call DSYEV('V','U',j,Tm,j,ev,wk,3*j,info2)
      lam = ev(1)
      resid = abs(nrm*Tm(j,1))
      if (resid .lt. 1.0d-4 .or. nrm .lt. 1.0d-8 .or. j .eq. mmax) then
         vec = matmul(Q(:,1:j), Tm(:,1))
         deallocate(Tm,ev,wk)
         exit
      endif
      deallocate(Tm,ev,wk)
   endif
   Q(:,j+1) = w/nrm
enddo
deallocate(Q,al,be,w)
end subroutine stab_lanczos

subroutine apply_A(v, y)
real(8),intent(in) :: v(ndt)
real(8),intent(out) :: y(ndt)
integer :: s2
real(8),allocatable :: dP(:,:,:), K2(:,:), T2(:,:), G(:,:), Gall(:,:,:), Gfd(:,:,:)
real(8) :: dmax
allocate(dP(nconts,nconts,2), K2(nconts,nconts), T2(nconts,nconts))
dP = 0.0d0
do s2 = 1,nsp
   allocate(kap(nv(s2),no(s2)))
   kap = reshape(v(off(s2)+1:off(s2)+nd(s2)),(/nv(s2),no(s2)/))
   call DGEMM('N','N',nconts,no(s2),nv(s2),1.0d0,Cs(1,no(s2)+1,s2),nconts,kap,nv(s2),0.0d0,T2,nconts)
   call DGEMM('N','T',nconts,nconts,no(s2),1.0d0,T2,nconts,Cs(1,1,s2),nconts,0.0d0,K2,nconts)
   dP(:,:,s2) = K2 + transpose(K2)
   deallocate(kap)
enddo
if (nsp .eq. 1) dP(:,:,2) = dP(:,:,1)
dmax = maxval(abs(dP))
if (dmax .lt. 1.0d-300) then
   y = 0.0d0
   return
endif
allocate(Gall(nconts,nconts,2))
call an_setup()
if (an_ok) then
   call an_G(dP, Gall)
   if (an_check .and. an_ncheck .lt. 3) then
      allocate(Gfd(nconts,nconts,2))
      call fd_G(dP, dmax, Gfd)
      an_ncheck = an_ncheck + 1
      print '(A,I2,A,ES10.3,A,ES10.3)', ' NEWTON check ', an_ncheck, ': |G_analytic - G_fd| / |G_fd| =', &
         sqrt(sum((Gall(:,:,1:nsp)-Gfd(:,:,1:nsp))**2)/max(sum(Gfd(:,:,1:nsp)**2),1.0d-300)), &
         '   |G_fd| =', sqrt(sum(Gfd(:,:,1:nsp)**2))
      deallocate(Gfd)
   endif
else
   call fd_G(dP, dmax, Gall)
endif
allocate(G(nconts,nconts))
do s2 = 1,nsp
   G = Gall(:,:,s2)
   allocate(kap(nv(s2),no(s2)))
   kap = reshape(v(off(s2)+1:off(s2)+nd(s2)),(/nv(s2),no(s2)/))
   call DGEMM('N','N',nconts,no(s2),nconts,1.0d0,G,nconts,Cs(1,1,s2),nconts,0.0d0,T2,nconts)
   call DGEMM('T','N',nv(s2),no(s2),nconts,1.0d0,Cs(1,no(s2)+1,s2),nconts,T2,nconts,0.0d0,K2,nconts)
   y(off(s2)+1:off(s2)+nd(s2)) = reshape(K2(1:nv(s2),1:no(s2)),(/nd(s2)/))
   call DGEMM('N','N',nv(s2),no(s2),nv(s2),1.0d0,Fvv(1,1,s2),nconts,kap,nv(s2),0.0d0,K2,nconts)
   y(off(s2)+1:off(s2)+nd(s2)) = y(off(s2)+1:off(s2)+nd(s2)) + reshape(K2(1:nv(s2),1:no(s2)),(/nd(s2)/))
   call DGEMM('N','N',nv(s2),no(s2),no(s2),-1.0d0,kap,nv(s2),Foo(1,1,s2),nconts,0.0d0,K2,nconts)
   y(off(s2)+1:off(s2)+nd(s2)) = y(off(s2)+1:off(s2)+nd(s2)) + reshape(K2(1:nv(s2),1:no(s2)),(/nd(s2)/))
   deallocate(kap)
enddo
deallocate(dP,K2,T2,G,Gall)
end subroutine apply_A

subroutine cpu_time_wall(t)
use omp_lib, only: omp_get_wtime
real(8),intent(out) :: t
t = omp_get_wtime()
end subroutine cpu_time_wall

subroutine an_setup()
use mod_exchange, only: HF_exchange_frac, RS_omega, cosx_enabled
use mod_integrals, only: engine_use_df_k
use mod_xc, only: xc_uses_vv10
use mod_cosmo, only: cosmo_enabled
character(len=8) :: e
if (an_done) return
an_done = .true.
an_ok = .false.; an_check = .false.; an_ncheck = 0
call get_environment_variable("ENGINE_NEWTON_ANALYTIC", e)
if (trim(e) .eq. "1") an_ok = .true.
call get_environment_variable("ENGINE_NEWTON_CHECK", e)
if (trim(e) .eq. "1") an_check = .true.
if (RS_omega .gt. 0.0d0) an_ok = .false.
if (xc_uses_vv10()) an_ok = .false.
if (cosmo_enabled) an_ok = .false.
if (HF_exchange_frac .gt. 0.0d0 .and. .not. cosx_enabled .and. engine_use_df_k) an_ok = .false.
end subroutine an_setup

subroutine an_G(dP, Gall)
use mod_exchange, only: HF_exchange_frac, cosx_enabled
use mod_integrals, only: integrals_build_coulomb_batch, integrals_build_coulomb_df_batch, &
                         integrals_build_exchange_batch, integrals_build_exchange_cosx_batch, engine_use_df_j, &
                         cosx_set_resp_follow_scf, integrals_set_accuracy
use mod_xc_hess, only: xc_fxc_apply_u
real(8),intent(in) :: dP(nconts,nconts,2)
real(8),intent(out) :: Gall(nconts,nconts,2)
real(8),allocatable :: Dt(:,:,:), Vj(:,:,:), Kx(:,:,:)
real(8) :: tj0, tj1, tk1, tx1
allocate(Dt(nconts,nconts,1), Vj(nconts,nconts,1), Kx(nconts,nconts,2))
call cpu_time_wall(tj0)
Dt(:,:,1) = dP(:,:,1) + dP(:,:,2)
call integrals_set_accuracy(sqrt(sum(dP**2)/(2.0d0*dble(nconts)**2)))
if (engine_use_df_j) then
   call integrals_build_coulomb_df_batch(nconts, 1, Dt, Vj)
else
   call integrals_build_coulomb_batch(nconts, 1, Dt, Vj)
endif
Gall(:,:,1) = Vj(:,:,1)
Gall(:,:,2) = Vj(:,:,1)
call cpu_time_wall(tj1)
if (HF_exchange_frac .gt. 0.0d0) then
   if (cosx_enabled) then
      call cosx_incr_stash(.true.)
      call cosx_set_resp_follow_scf(.true.)
      if (nsp .eq. 1) then
         call integrals_build_exchange_cosx_batch(nconts, 1, dP(:,:,1:1), Kx(:,:,1:1))
         Kx(:,:,2) = Kx(:,:,1)
      else
         call integrals_build_exchange_cosx_batch(nconts, 2, dP, Kx)
      endif
      call cosx_set_resp_follow_scf(.false.)
      call cosx_incr_stash(.false.)
   else
      call integrals_build_exchange_batch(nconts, 1, dP(:,:,1:1), dP(:,:,2:2), Kx(:,:,1:1), Kx(:,:,2:2))
   endif
   Gall = Gall - HF_exchange_frac*Kx
endif
call cpu_time_wall(tk1)
if (HF_exchange_frac .lt. 1.0d0) &
   call xc_fxc_apply_u(nconts, 1, P0(:,:,1), P0(:,:,2), dP(:,:,1:1), dP(:,:,2:2), Gall(:,:,1:1), Gall(:,:,2:2))
call cpu_time_wall(tx1)
if (an_check) then
   print '(A,3F8.3,A)', ' NEWTON product time J, K, XC:', tj1-tj0, tk1-tj1, tx1-tk1, ' s'
   call flush(6)
   block
      character(len=8) :: es
      integer,save :: nstop = -1, nseen = 0
      if (nstop .lt. 0) then
         call get_environment_variable("ENGINE_NEWTON_STOP", es)
         nstop = 0
         if (len_trim(es) .gt. 0) read(es,*) nstop
      endif
      nseen = nseen + 1
      if (nstop .gt. 0 .and. nseen .ge. nstop) stop
   end block
endif
deallocate(Dt,Vj,Kx)
end subroutine an_G

subroutine fd_G(dP, dmax, Gall)
use mod_exchange, only: cosx_enabled
use mod_integrals, only: cosx_set_response_mode, cosx_incr_stash
real(8),intent(in) :: dP(nconts,nconts,2), dmax
real(8),intent(out) :: Gall(nconts,nconts,2)
integer :: s2
real(8),allocatable :: Pa_sv(:,:), Pb_sv(:,:), Pa_p(:,:), Pb_p(:,:)
real(8) :: h, Exc_l, En_l, dt_l
logical :: fd_inc
character(len=8) :: einc
call get_environment_variable("ENGINE_NEWTON_FDINC", einc)
fd_inc = (trim(einc) .ne. "0")
h = 0.01d0/dmax
block
   character(len=32) :: he
   call get_environment_variable("ENGINE_NEWTON_H", he)
   if (len_trim(he) .gt. 0) read(he,*) h
   if (len_trim(he) .gt. 0) h = h/dmax
end block
allocate(Pa_sv(nconts,nconts), Pb_sv(nconts,nconts), Pa_p(nconts,nconts), Pb_p(nconts,nconts))
Pa_sv = Pa; Pb_sv = Pb
Pa_p = P0(:,:,1) + h*dP(:,:,1)
Pb_p = P0(:,:,2) + h*dP(:,:,2)
Pa = Pa_p; Pb = Pb_p
if (cosx_enabled .and. fd_inc) then
   call cosx_incr_stash(.true.)
else if (cosx_enabled) then
   call cosx_incr_stash(.true.)
   call cosx_set_response_mode(.true.)
   if (.not. have_F0r) then
      if (.not. allocated(F0r)) allocate(F0r(nconts,nconts,2))
      Pa = P0(:,:,1); Pb = P0(:,:,2)
      call scf_build_fock(P0(:,:,1), P0(:,:,2), .false., 99, Prms, Exc_l, En_l, dt_l, .false.)
      F0r(:,:,1) = Fa; F0r(:,:,2) = Fb
      have_F0r = .true.
      Pa = Pa_p; Pb = Pb_p
   endif
endif
call scf_build_fock(Pa_p, Pb_p, .false., 99, Prms, Exc_l, En_l, dt_l, .false.)
if (cosx_enabled .and. fd_inc) then
   call cosx_incr_stash(.false.)
else if (cosx_enabled) then
   call cosx_set_response_mode(.false.)
   call cosx_incr_stash(.false.)
endif
do s2 = 1,nsp
   if (cosx_enabled .and. .not. fd_inc) then
      if (s2 .eq. 1) then
         Gall(:,:,1) = (Fa - F0r(:,:,1))/h
      else
         Gall(:,:,2) = (Fb - F0r(:,:,2))/h
      endif
   else if (s2 .eq. 1) then
      Gall(:,:,1) = (Fa - F0(:,:,1))/h
   else
      Gall(:,:,2) = (Fb - F0(:,:,2))/h
   endif
enddo
Pa = Pa_sv; Pb = Pb_sv
Fa = F0(:,:,1); Fb = F0(:,:,2)
if (nsp .eq. 1) Fb = F0(:,:,1)
deallocate(Pa_sv,Pb_sv,Pa_p,Pb_p)
end subroutine fd_G

end subroutine newton_cg_step

end module mod_soscf
