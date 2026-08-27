! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Top-level Hartree-Fock/Kohn-Sham solver entry point (solHFR_KS).

subroutine solHFR_KS(infor, Prms_prev, near_singular_S)
use MOL_info
use omp_lib, only: omp_get_max_threads
    implicit none
INCLUDE 'parameter.h'
    interface
       subroutine openblas_set_num_threads(num_threads) bind(c, name="openblas_set_num_threads")
       use iso_c_binding, only: c_int
       integer(c_int),value :: num_threads
       end subroutine openblas_set_num_threads
    end interface
    integer    :: i,j,k,l,n,m,start,info,infor
    integer    :: ij,kl
    integer    :: LWORK,LIWORK
    real(8),allocatable   ::  WORK(:)
    integer,allocatable   ::  IWORK(:)
    real(8),allocatable   ::  e1(:),e2(:)
    real(8),allocatable   ::  Tmp(:,:)
    real(8),intent(in) :: Prms_prev
    logical,intent(in) :: near_singular_S
    real(8),save :: level_shift = 0.0d0
    real(8),parameter :: LEVEL_SHIFT_OFF_PRMS = 1.0d-4
    logical,save :: level_shift_env_read = .false.
    character(len=32) :: ls_envchar
    integer :: ls_envstat, nvirt_a, nvirt_b
    real(8),allocatable :: SCv(:,:)
    real(8),allocatable :: C_a_prev(:,:), C_b_prev(:,:)
    logical :: apply_shift

    LWORK = 1 + 6*nconts + 2*nconts**2
    LIWORK = 3 + 5*nconts
    allocate (WORK(LWORK))
    allocate (IWORK(LIWORK))
    allocate (Tmp(nconts,nconts))

    call openblas_set_num_threads(omp_get_max_threads())

    if (.not. level_shift_env_read) then
       call get_environment_variable("ENGINE_LEVEL_SHIFT", ls_envchar, status=ls_envstat)
       if (ls_envstat .eq. 0) then
          read(ls_envchar,*) level_shift
       else if (near_singular_S) then
          level_shift = 0.3d0
          print *,"near-singular overlap matrix detected: auto-enabling", &
                  " level_shift=0.3 Ha to stabilize DIIS (override with ENGINE_LEVEL_SHIFT)"
       endif
       level_shift_env_read = .true.
    endif
    apply_shift = (level_shift .ne. 0.0d0) .and. (Prms_prev .ge. LEVEL_SHIFT_OFF_PRMS)
    if (apply_shift) then
       allocate(C_a_prev(nconts,nconts), C_b_prev(nconts,nconts))
       C_a_prev = C_a
       C_b_prev = C_b
    endif

    if (apply_shift) then
       nvirt_a = nconts - n_alpha
       if (nvirt_a .gt. 0) then
          allocate(SCv(nconts,nvirt_a))
          call DGEMM('N','N',nconts,nvirt_a,nconts,1.0d0,S,nconts, &
                     C_a_prev(1,n_alpha+1),nconts,0.0d0,SCv,nconts)
          call DGEMM('N','T',nconts,nconts,nvirt_a,level_shift,SCv,nconts, &
                     SCv,nconts,1.0d0,Fa,nconts)
          deallocate(SCv)
       endif
    endif
    call DGEMM('N','N',nconts,nconts,nconts,1.0d0,Fa,nconts,X,nconts,0.0d0,Tmp,nconts)
    call DGEMM('T','N',nconts,nconts,nconts,1.0d0,X,nconts,Tmp,nconts,0.0d0,Fa,nconts)
    C_a = 0
    C_b = 0

    allocate(e1(nconts))
    call DSYEVD('V','U',nconts,Fa,nconts,e1,WORK,LWORK,IWORK,LIWORK,INFO)
    eLev_a = e1
    deallocate(e1)

    call DGEMM('N','N',nconts,nconts,nconts,1.0d0,X,nconts,Fa,nconts,0.0d0,C_a,nconts)

    if(multi .eq. 1) then
        C_b = C_a
        eLev_b = eLev_a
    else
        allocate(e2(nconts))
        if (apply_shift) then
           nvirt_b = nconts - n_beta
           if (nvirt_b .gt. 0) then
              allocate(SCv(nconts,nvirt_b))
              call DGEMM('N','N',nconts,nvirt_b,nconts,1.0d0,S,nconts, &
                         C_b_prev(1,n_beta+1),nconts,0.0d0,SCv,nconts)
              call DGEMM('N','T',nconts,nconts,nvirt_b,level_shift,SCv,nconts, &
                         SCv,nconts,1.0d0,Fb,nconts)
              deallocate(SCv)
           endif
        endif
        call DGEMM('N','N',nconts,nconts,nconts,1.0d0,Fb,nconts,X,nconts,0.0d0,Tmp,nconts)
        call DGEMM('T','N',nconts,nconts,nconts,1.0d0,X,nconts,Tmp,nconts,0.0d0,Fb,nconts)
        call DSYEVD('V','U',nconts,Fb,nconts,e2,WORK,LWORK,IWORK,LIWORK,INFO)
        eLev_b = e2
        call DGEMM('N','N',nconts,nconts,nconts,1.0d0,X,nconts,Fb,nconts,0.0d0,C_b,nconts)
        deallocate(e2)
    endif
    if (apply_shift) deallocate(C_a_prev, C_b_prev)

    call openblas_set_num_threads(1)

    deallocate(Tmp)
    deallocate(WORK)
    deallocate(IWORK)

end subroutine
