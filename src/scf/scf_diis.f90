! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! scf_diis: the DIIS/ADIIS Fock-matrix extrapolation step extracted from


subroutine diis_extrapolate(nconts, multi, iter, Fa, Fb, Pa, Pb, S, X, &
                             DIIS_MAX, diis_n, diis_Fa, diis_Fb, diis_Pa, diis_Pb, &
                             diis_ea, diis_eb, diis_debug, diis_active)
implicit none
integer,intent(in) :: nconts,multi,iter,DIIS_MAX
real(8),intent(inout) :: Fa(nconts,nconts),Fb(nconts,nconts)
real(8),intent(in) :: Pa(nconts,nconts),Pb(nconts,nconts)
real(8),intent(in) :: S(nconts,nconts),X(nconts,nconts)
integer,intent(inout) :: diis_n
real(8),intent(inout) :: diis_Fa(nconts,nconts,DIIS_MAX),diis_Fb(nconts,nconts,DIIS_MAX)
real(4),intent(inout) :: diis_Pa(nconts,nconts,DIIS_MAX),diis_Pb(nconts,nconts,DIIS_MAX)
real(4),intent(inout) :: diis_ea(nconts,nconts,DIIS_MAX),diis_eb(nconts,nconts,DIIS_MAX)
logical,intent(in) :: diis_debug
logical,intent(out) :: diis_active

integer :: p,q,nsz,info_diis
real(8) :: raw_diag_n,diis_error,aediis_coeff
real(8),allocatable :: e_ortho(:,:)
real(8),allocatable :: Bmat(:,:),rhs(:),diagscale(:),Bsvd(:)
real(8),allocatable :: work_lapack(:)
integer :: lwork_lapack,rank_lapack
real(8),allocatable :: adiis_lin(:),adiis_quad(:,:),adiis_w(:),blend_w(:)
real(8),parameter :: DIIS_COEF_MAX = 10.0d0
real(8),parameter :: DIIS_EPSILON = 0.05d0, DIIS_THRESHOLD = 1.0d-5
real(8),parameter :: DIIS_REF_NCONTS = 108.0d0

diis_active = .false.

allocate(e_ortho(nconts,nconts))
call diis_error_matrix(nconts, Fa, Pa, S, X, e_ortho)
if (diis_n .eq. DIIS_MAX) then
   diis_Fa(:,:,1:DIIS_MAX-1) = diis_Fa(:,:,2:DIIS_MAX)
   diis_ea(:,:,1:DIIS_MAX-1) = diis_ea(:,:,2:DIIS_MAX)
   diis_Pa(:,:,1:DIIS_MAX-1) = diis_Pa(:,:,2:DIIS_MAX)
   if (multi .ne. 1) then
      diis_Fb(:,:,1:DIIS_MAX-1) = diis_Fb(:,:,2:DIIS_MAX)
      diis_eb(:,:,1:DIIS_MAX-1) = diis_eb(:,:,2:DIIS_MAX)
      diis_Pb(:,:,1:DIIS_MAX-1) = diis_Pb(:,:,2:DIIS_MAX)
   endif
   diis_n = DIIS_MAX - 1
endif
diis_n = diis_n + 1
diis_Fa(:,:,diis_n) = Fa
diis_ea(:,:,diis_n) = e_ortho
diis_Pa(:,:,diis_n) = Pa
if (multi .ne. 1) then
   call diis_error_matrix(nconts, Fb, Pb, S, X, e_ortho)
   diis_Fb(:,:,diis_n) = Fb
   diis_eb(:,:,diis_n) = e_ortho
   diis_Pb(:,:,diis_n) = Pb
endif
deallocate(e_ortho)

if (diis_n .ge. 2) then
   nsz = diis_n+1
   allocate(Bmat(nsz,nsz),rhs(nsz),diagscale(nsz),Bsvd(nsz))
   do p = 1,diis_n
      do q = 1,diis_n
         Bmat(p,q) = sum(real(diis_ea(:,:,p),8)*real(diis_ea(:,:,q),8))
         if (multi .ne. 1) Bmat(p,q) = 0.5d0*(Bmat(p,q) + sum(real(diis_eb(:,:,p),8)*real(diis_eb(:,:,q),8)))
      enddo
   enddo
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
   call DGELSS(nsz,nsz,1,Bmat,nsz,rhs,nsz,Bsvd,1.0d-5,rank_lapack, &
               work_lapack,lwork_lapack,info_diis)
   if (diis_debug) then
      print *,"RTDBG DGELSS done"
      call flush(6)
   endif
   deallocate(work_lapack)
   rhs = rhs*diagscale
   deallocate(diagscale,Bsvd)
   deallocate(Bmat)

   if (diis_error .lt. DIIS_THRESHOLD) then
      aediis_coeff = 0.0d0
   else if (diis_error .lt. DIIS_EPSILON) then
      aediis_coeff = (diis_error-DIIS_THRESHOLD)/(DIIS_EPSILON-DIIS_THRESHOLD)
   else
      aediis_coeff = 1.0d0
   endif

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
      Fa = 0
      if (multi .ne. 1) Fb = 0
      do p = 1,diis_n
         Fa = Fa + blend_w(p)*diis_Fa(:,:,p)
         if (multi .ne. 1) Fb = Fb + blend_w(p)*diis_Fb(:,:,p)
      enddo
      diis_active = .true.
   endif
   deallocate(blend_w)
   deallocate(rhs)
endif
end subroutine diis_extrapolate

subroutine diis_error_matrix(nconts, F, P, S, X, e_ortho)
implicit none
integer,intent(in) :: nconts
real(8),intent(in) :: F(nconts,nconts),P(nconts,nconts)
real(8),intent(in) :: S(nconts,nconts),X(nconts,nconts)
real(8),intent(out) :: e_ortho(nconts,nconts)
real(8),allocatable :: t1(:,:),e_ao(:,:),t3(:,:)

allocate(t1(nconts,nconts),e_ao(nconts,nconts),t3(nconts,nconts))

call dgemm('N','N',nconts,nconts,nconts,1.0d0,P,nconts,S,nconts,0.0d0,t1,nconts)
call dgemm('N','N',nconts,nconts,nconts,1.0d0,F,nconts,t1,nconts,0.0d0,e_ao,nconts)
call dgemm('N','N',nconts,nconts,nconts,1.0d0,P,nconts,F,nconts,0.0d0,t1,nconts)
call dgemm('N','N',nconts,nconts,nconts,-1.0d0,S,nconts,t1,nconts,1.0d0,e_ao,nconts)

call dgemm('N','N',nconts,nconts,nconts,1.0d0,e_ao,nconts,X,nconts,0.0d0,t3,nconts)
call dgemm('T','N',nconts,nconts,nconts,1.0d0,X,nconts,t3,nconts,0.0d0,e_ortho,nconts)

deallocate(t1,e_ao,t3)
end subroutine diis_error_matrix
