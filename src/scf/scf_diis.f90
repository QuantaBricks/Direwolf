! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! scf_diis: the DIIS/ADIIS Fock-matrix extrapolation step extracted from


subroutine diis_extrapolate(nconts, multi, iter, Fa, Fb, Pa, Pb, S, X, &
                             DIIS_MAX, diis_n, diis_Fa, diis_Fb, diis_Pa, diis_Pb, &
                             diis_ea, diis_eb, diis_debug, diis_active)
use omp_lib, only: omp_get_max_threads
implicit none
interface
   subroutine openblas_set_num_threads(num_threads) bind(c, name="openblas_set_num_threads")
      use iso_c_binding, only: c_int
      integer(c_int), value :: num_threads
   end subroutine openblas_set_num_threads
end interface
integer,intent(in) :: nconts,multi,iter,DIIS_MAX
real(8),intent(inout) :: Fa(nconts,nconts),Fb(nconts,nconts)
real(8),intent(in) :: Pa(nconts,nconts),Pb(nconts,nconts)
real(8),intent(in) :: S(nconts,nconts),X(nconts,nconts)
integer,intent(inout) :: diis_n
real(8),intent(inout) :: diis_Fa(nconts*(nconts+1)/2,DIIS_MAX),diis_Fb(nconts*(nconts+1)/2,DIIS_MAX)
real(4),intent(inout) :: diis_Pa(nconts*(nconts+1)/2,DIIS_MAX),diis_Pb(nconts*(nconts+1)/2,DIIS_MAX)
real(4),intent(inout) :: diis_ea(nconts*(nconts+1)/2,DIIS_MAX),diis_eb(nconts*(nconts+1)/2,DIIS_MAX)
logical,intent(in) :: diis_debug
logical,intent(out) :: diis_active

integer :: p,q,nsz,info_diis
character(len=8) :: noadiis
character(len=16) :: epse
real(8) :: eps_eff
real(8) :: raw_diag_n,diis_error,aediis_coeff
real(8),allocatable :: e_ortho(:,:)
real(8),allocatable :: Bmat(:,:),rhs(:),diagscale(:),Bsvd(:),Fsum(:)
integer :: npk
real(8),external :: diis_wdot44
real(8),allocatable :: work_lapack(:)
integer :: lwork_lapack,rank_lapack
real(8),allocatable :: adiis_lin(:),adiis_quad(:,:),adiis_w(:),blend_w(:)
real(8) :: DIIS_COEF_MAX, diis_rcond
character(len=32) :: envtmp
real(8),parameter :: DIIS_EPSILON = 0.05d0, DIIS_THRESHOLD = 1.0d-5
real(8),parameter :: DIIS_REF_NCONTS = 108.0d0

diis_active = .false.
DIIS_COEF_MAX = 1.0d6
diis_rcond = 1.0d-14
call get_environment_variable("ENGINE_DIIS_COEFMAX", envtmp)
if (len_trim(envtmp) .gt. 0) read(envtmp,*) DIIS_COEF_MAX
call get_environment_variable("ENGINE_DIIS_RCOND", envtmp)
if (len_trim(envtmp) .gt. 0) read(envtmp,*) diis_rcond

allocate(e_ortho(nconts,nconts))
call openblas_set_num_threads(omp_get_max_threads())
call diis_error_matrix(nconts, Fa, Pa, S, X, e_ortho)
if (diis_n .eq. DIIS_MAX) then
   diis_Fa(:,1:DIIS_MAX-1) = diis_Fa(:,2:DIIS_MAX)
   diis_ea(:,1:DIIS_MAX-1) = diis_ea(:,2:DIIS_MAX)
   diis_Pa(:,1:DIIS_MAX-1) = diis_Pa(:,2:DIIS_MAX)
   if (multi .ne. 1) then
      diis_Fb(:,1:DIIS_MAX-1) = diis_Fb(:,2:DIIS_MAX)
      diis_eb(:,1:DIIS_MAX-1) = diis_eb(:,2:DIIS_MAX)
      diis_Pb(:,1:DIIS_MAX-1) = diis_Pb(:,2:DIIS_MAX)
   endif
   diis_n = DIIS_MAX - 1
endif
diis_n = diis_n + 1
npk = nconts*(nconts+1)/2
call diis_pack8(nconts, Fa, diis_Fa(:,diis_n))
call diis_pack4(nconts, e_ortho, diis_ea(:,diis_n))
call diis_pack4(nconts, Pa, diis_Pa(:,diis_n))
if (multi .ne. 1) then
   call diis_error_matrix(nconts, Fb, Pb, S, X, e_ortho)
   call diis_pack8(nconts, Fb, diis_Fb(:,diis_n))
   call diis_pack4(nconts, e_ortho, diis_eb(:,diis_n))
   call diis_pack4(nconts, Pb, diis_Pb(:,diis_n))
endif
call openblas_set_num_threads(1)
deallocate(e_ortho)

if (diis_n .ge. 2) then
   nsz = diis_n+1
   allocate(Bmat(nsz,nsz),rhs(nsz),diagscale(nsz),Bsvd(nsz))
   !$omp parallel do schedule(dynamic,1) default(shared) private(p,q)
   do p = 1,diis_n
      do q = p,diis_n
         Bmat(p,q) = diis_wdot44(nconts, diis_ea(:,p), diis_ea(:,q))
         if (multi .ne. 1) Bmat(p,q) = 0.5d0*(Bmat(p,q) + diis_wdot44(nconts, diis_eb(:,p), diis_eb(:,q)))
         Bmat(q,p) = Bmat(p,q)
      enddo
   enddo
   !$omp end parallel do
   raw_diag_n = Bmat(diis_n,diis_n)
   diis_error = sqrt(max(raw_diag_n,0.0d0)) * (DIIS_REF_NCONTS/real(nconts,8))
   Bmat(1:diis_n,nsz) = -1.0d0
   Bmat(nsz,1:diis_n) = -1.0d0
   Bmat(nsz,nsz) = 0.0d0
   rhs = 0.0d0
   rhs(nsz) = -1.0d0
   diagscale = 1.0d0
   if (minval((/(Bmat(p,p),p=1,diis_n)/)) .gt. 0.0d0) then
      do p = 1,diis_n
         diagscale(p) = 1.0d0/sqrt(Bmat(p,p))
      enddo
   endif
   do p = 1,nsz
      Bmat(p,:) = Bmat(p,:)*diagscale(p)
      Bmat(:,p) = Bmat(:,p)*diagscale(p)
   enddo
   lwork_lapack = max(1, 3*nsz + max(2*nsz,nsz,1), 5*nsz)
   allocate(work_lapack(lwork_lapack))
   if (diis_debug) then
      print *,"RTDBG entering DGELSS iter=",iter,"diis_n=",diis_n
      call flush(6)
   endif
   call DGELSS(nsz,nsz,1,Bmat,nsz,rhs,nsz,Bsvd,diis_rcond,rank_lapack, &
               work_lapack,lwork_lapack,info_diis)
   if (diis_debug) then
      print *,"RTDBG DGELSS done"
      call flush(6)
   endif
   deallocate(work_lapack)
   rhs = rhs*diagscale
   deallocate(diagscale,Bsvd)
   deallocate(Bmat)

   eps_eff = 1.0d0
   call get_environment_variable("ENGINE_DIIS_EPS", epse)
   if (len_trim(epse) .gt. 0) read(epse,*) eps_eff
   if (diis_error .lt. DIIS_THRESHOLD) then
      aediis_coeff = 0.0d0
   else if (diis_error .lt. eps_eff) then
      aediis_coeff = (diis_error-DIIS_THRESHOLD)/(eps_eff-DIIS_THRESHOLD)
   else
      aediis_coeff = 1.0d0
   endif

   call get_environment_variable("ENGINE_NO_ADIIS", noadiis)
   if (trim(noadiis) .eq. "1") aediis_coeff = 0.0d0
   allocate(blend_w(diis_n))
   blend_w = rhs(1:diis_n)
   if (diis_debug) then
      print '("  DIIS dbg iter",I3," diis_n",I3," Dnorm(rms)=",ES12.5," aediis_coeff=",F8.5)', &
            iter,diis_n,diis_error,aediis_coeff
      print '("    DIIS  weights: ",20F10.4)', rhs(1:diis_n)
   endif
   if (aediis_coeff .gt. 0.0d0) then
      allocate(adiis_lin(diis_n),adiis_quad(diis_n,diis_n),adiis_w(diis_n))
      if (diis_debug) then
         print *,"RTDBG entering build_adiis_terms"
         call flush(6)
      endif
      call build_adiis_terms(diis_n,nconts,diis_Pa,diis_Fa, &
           diis_Pb,diis_Fb,multi.ne.1,adiis_lin,adiis_quad)
      if (diis_debug) then
         print *,"RTDBG build_adiis_terms done, entering adiis_powell"
         call flush(6)
      endif
      call adiis_powell(diis_n,adiis_quad,adiis_lin,adiis_w)
      if (diis_debug) then
         print *,"RTDBG adiis_powell done"
         call flush(6)
      endif
      if (diis_debug) print '("    ADIIS weights: ",20F10.4)', adiis_w
      blend_w = aediis_coeff*adiis_w + (1.0d0-aediis_coeff)*rhs(1:diis_n)
      deallocate(adiis_lin,adiis_quad,adiis_w)
   endif
   if (diis_debug) print '("    Blend weights: ",20F10.4)', blend_w

   if (info_diis .eq. 0 .and. maxval(abs(blend_w)) .le. DIIS_COEF_MAX) then
      allocate(Fsum(npk))
      Fsum = 0
      do p = 1,diis_n
         Fsum = Fsum + blend_w(p)*diis_Fa(:,p)
      enddo
      call diis_unpack8(nconts, Fsum, Fa)
      if (multi .ne. 1) then
         Fsum = 0
         do p = 1,diis_n
            Fsum = Fsum + blend_w(p)*diis_Fb(:,p)
         enddo
         call diis_unpack8(nconts, Fsum, Fb)
      endif
      deallocate(Fsum)
      diis_active = .true.
   endif
   deallocate(blend_w)
   deallocate(rhs)
endif
end subroutine diis_extrapolate

subroutine diis_error_matrix(nconts, F, P, S, X, e_ortho)
use mod_gemm_big, only: gemm_big
implicit none
integer,intent(in) :: nconts
real(8),intent(in) :: F(nconts,nconts),P(nconts,nconts)
real(8),intent(in) :: S(nconts,nconts),X(nconts,nconts)
real(8),intent(out) :: e_ortho(nconts,nconts)
real(8),allocatable :: t1(:,:),t2(:,:)
integer :: i, j
real(8) :: bij

allocate(t1(nconts,nconts),t2(nconts,nconts))

call gemm_big('N','N',nconts,P,S,t1)
call gemm_big('N','N',nconts,F,t1,t2)
call gemm_big('N','N',nconts,t2,X,t1)
call gemm_big('T','N',nconts,X,t1,t2)
!$omp parallel do private(i, j, bij) schedule(static)
do j = 1, nconts
   do i = 1, j - 1
      bij = t2(i,j) - t2(j,i)
      e_ortho(i,j) = bij
      e_ortho(j,i) = -bij
   enddo
   e_ortho(j,j) = 0.0d0
enddo
!$omp end parallel do

deallocate(t1,t2)
end subroutine diis_error_matrix

subroutine diis_shift_history(nconts, multi, dFa, dFb, S, X, diis_n, DIIS_MAX, &
                              diis_Fa, diis_Fb, diis_Pa, diis_Pb, diis_ea, diis_eb)
implicit none
integer,intent(in) :: nconts, multi, diis_n, DIIS_MAX
real(8),intent(in) :: dFa(nconts,nconts), dFb(nconts,nconts)
real(8),intent(in) :: S(nconts,nconts), X(nconts,nconts)
real(8),intent(inout) :: diis_Fa(nconts*(nconts+1)/2,DIIS_MAX), diis_Fb(nconts*(nconts+1)/2,DIIS_MAX)
real(4),intent(in) :: diis_Pa(nconts*(nconts+1)/2,DIIS_MAX), diis_Pb(nconts*(nconts+1)/2,DIIS_MAX)
real(4),intent(inout) :: diis_ea(nconts*(nconts+1)/2,DIIS_MAX), diis_eb(nconts*(nconts+1)/2,DIIS_MAX)
real(8),allocatable :: e_ortho(:,:), P8(:,:), Ff(:,:)
integer :: p
allocate(e_ortho(nconts,nconts), P8(nconts,nconts), Ff(nconts,nconts))
do p = 1,diis_n
   call diis_unpack8(nconts, diis_Fa(:,p), Ff)
   Ff = Ff + dFa
   call diis_pack8(nconts, Ff, diis_Fa(:,p))
   call diis_unpack4(nconts, diis_Pa(:,p), P8)
   call diis_error_matrix(nconts, Ff, P8, S, X, e_ortho)
   call diis_pack4(nconts, e_ortho, diis_ea(:,p))
   if (multi .ne. 1) then
      call diis_unpack8(nconts, diis_Fb(:,p), Ff)
      Ff = Ff + dFb
      call diis_pack8(nconts, Ff, diis_Fb(:,p))
      call diis_unpack4(nconts, diis_Pb(:,p), P8)
      call diis_error_matrix(nconts, Ff, P8, S, X, e_ortho)
      call diis_pack4(nconts, e_ortho, diis_eb(:,p))
   endif
enddo
deallocate(e_ortho, P8, Ff)
end subroutine diis_shift_history

subroutine diis_pack8(n, A, Ap)
implicit none
integer,intent(in) :: n
real(8),intent(in) :: A(n,n)
real(8),intent(out) :: Ap(n*(n+1)/2)
integer :: i, j, off
!$omp parallel do schedule(static) default(shared) private(i, off)
do j = 1, n
   off = j*(j-1)/2
   do i = 1, j
      Ap(off+i) = A(i,j)
   enddo
enddo
!$omp end parallel do
end subroutine diis_pack8

subroutine diis_pack4(n, A, Ap)
implicit none
integer,intent(in) :: n
real(8),intent(in) :: A(n,n)
real(4),intent(out) :: Ap(n*(n+1)/2)
integer :: i, j, off
!$omp parallel do schedule(static) default(shared) private(i, off)
do j = 1, n
   off = j*(j-1)/2
   do i = 1, j
      Ap(off+i) = real(A(i,j),4)
   enddo
enddo
!$omp end parallel do
end subroutine diis_pack4

subroutine diis_unpack8(n, Ap, A)
implicit none
integer,intent(in) :: n
real(8),intent(in) :: Ap(n*(n+1)/2)
real(8),intent(out) :: A(n,n)
integer :: i, j, off
!$omp parallel do schedule(static) default(shared) private(i, off)
do j = 1, n
   off = j*(j-1)/2
   do i = 1, j
      A(i,j) = Ap(off+i)
      A(j,i) = Ap(off+i)
   enddo
enddo
!$omp end parallel do
end subroutine diis_unpack8

subroutine diis_unpack4(n, Ap, A)
implicit none
integer,intent(in) :: n
real(4),intent(in) :: Ap(n*(n+1)/2)
real(8),intent(out) :: A(n,n)
integer :: i, j, off
!$omp parallel do schedule(static) default(shared) private(i, off)
do j = 1, n
   off = j*(j-1)/2
   do i = 1, j
      A(i,j) = real(Ap(off+i),8)
      A(j,i) = real(Ap(off+i),8)
   enddo
enddo
!$omp end parallel do
end subroutine diis_unpack4

function diis_wdot44(n, a, b) result(r)
implicit none
integer,intent(in) :: n
real(4),intent(in) :: a(n*(n+1)/2), b(n*(n+1)/2)
real(8) :: r
real(8) :: tri, dg
integer :: k, j
tri = 0.0d0
do k = 1, n*(n+1)/2
   tri = tri + real(a(k),8)*real(b(k),8)
enddo
dg = 0.0d0
do j = 1, n
   dg = dg + real(a(j*(j+1)/2),8)*real(b(j*(j+1)/2),8)
enddo
r = 2.0d0*tri - dg
end function diis_wdot44
