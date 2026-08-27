! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! test_coulomb: correctness gate for Phase 3's symmetrized Coulomb

program test_coulomb
use MOL_info
use mod_integrals, only: integrals_init, eri_get, integrals_build_coulomb, &
                          integrals_finalize
implicit none
include "parameter.h"

integer :: info, i, j, k, l
real(8),allocatable :: Ptot(:,:), J_brute(:,:), J_store(:,:), J_direct(:,:)
real(8) :: maxdiff_A, maxdiff_B
integer :: n

print *,"=== Test A: build_coulomb_store vs brute-force reference (STORE mode) ==="
call setup_water()
call integrals_init(info)
n = nConts
allocate(Ptot(n,n), J_brute(n,n), J_store(n,n))
Ptot = S

J_brute = 0
do i = 1,n
   do j = 1,n
      do k = 1,n
         do l = 1,n
            J_brute(i,j) = J_brute(i,j) + eri_get(i,j,k,l)*Ptot(k,l)
         enddo
      enddo
   enddo
enddo

call integrals_build_coulomb(n, Ptot, J_store)

maxdiff_A = maxval(abs(J_brute - J_store))
print *,"max|J_brute - J_store| =",maxdiff_A
if (maxdiff_A < 1.0d-10) then
   print *,"TEST A: PASS"
else
   print *,"TEST A: FAIL"
   stop 1
endif

print *,"=== Test B: build_coulomb_direct vs build_coulomb_store (same molecule) ==="
call setenv_max_gb("1000")
call reset_engine_state()
call setup_water()
call integrals_init(info)
call integrals_build_coulomb(nConts, Ptot, J_store)

call setenv_max_gb("0.0000001")
call reset_engine_state()
call setup_water()
call integrals_init(info)
allocate(J_direct(n,n))
call integrals_build_coulomb(nConts, Ptot, J_direct)

maxdiff_B = maxval(abs(J_store - J_direct))
print *,"max|J_store - J_direct| =",maxdiff_B
if (maxdiff_B < 1.0d-10) then
   print *,"TEST B: PASS"
else
   print *,"TEST B: FAIL"
   stop 1
endif

print *,"=== Test C: build_coulomb_store vs brute-force reference (REAL density) ==="
block
real(8),allocatable :: Fprime(:,:), Cprime(:,:), Cmo(:,:), eval(:), WORK(:)
real(8),allocatable :: Ptot_real(:,:), J_brute_real(:,:), J_store_real(:,:)
integer :: LWORK, INFO2, iocc
real(8) :: maxdiff_C

call setenv_max_gb("1000")
call reset_engine_state()
call setup_water()
call integrals_init(info)
n = nConts
allocate(Fprime(n,n), Cprime(n,n), Cmo(n,n), eval(n))
LWORK = 1 + 6*n + n**2
allocate(WORK(LWORK))
Fprime = matmul(transpose(X), matmul(Hcore, X))
Cprime = Fprime
call DSYEV('V','U',n,Cprime,n,eval,WORK,LWORK,INFO2)
Cmo = matmul(X,Cprime)

allocate(Ptot_real(n,n), J_brute_real(n,n), J_store_real(n,n))
Ptot_real = 0
do i = 1,n
   do j = 1,n
      do iocc = 1,5
         Ptot_real(i,j) = Ptot_real(i,j) + 2.0d0*Cmo(i,iocc)*Cmo(j,iocc)
      enddo
   enddo
enddo

J_brute_real = 0
do i = 1,n
   do j = 1,n
      do k = 1,n
         do l = 1,n
            J_brute_real(i,j) = J_brute_real(i,j) + eri_get(i,j,k,l)*Ptot_real(k,l)
         enddo
      enddo
   enddo
enddo

call integrals_build_coulomb(n, Ptot_real, J_store_real)

maxdiff_C = maxval(abs(J_brute_real - J_store_real))
print *,"max|J_brute_real - J_store_real| =",maxdiff_C
if (maxdiff_C < 1.0d-10) then
   print *,"TEST C: PASS"
else
   print *,"TEST C: FAIL"
   stop 1
endif
end block

print *,"ALL TESTS PASSED"

contains

subroutine setup_water()
    implicit none
    integer :: ii
    integer :: ncenters
    integer :: atomchg(3)
    real(8) :: coord(3,3)
    character(30) :: baselabel

    ncenters = 3
    atomchg = (/8,1,1/)
    coord(1,:) = (/0.0d0, 0.0d0, 0.0d0/)
    coord(2,:) = (/0.7584d0, 0.0d0, 0.5861d0/)
    coord(3,:) = (/-0.7584d0, 0.0d0, 0.5861d0/)
    baselabel = '6-31g'

    Natoms = ncenters
    Charge = 0
    Multi = 1
    allocate(atoms(Natoms))
    do ii = 1,Natoms
       atoms(ii)%coor = coord(ii,:)
       atoms(ii)%charge = atomchg(ii)
       write (atoms(ii)%base,"(I0.2,A1,A27)") atomchg(ii),'-',baselabel
    enddo
end subroutine setup_water

subroutine setenv_max_gb(s)
    use iso_c_binding, only: c_char, c_int
    implicit none
    character(*),intent(in) :: s
    interface
       function c_setenv(name,value,overwrite) bind(C,name="setenv") result(r)
          import :: c_char, c_int
          character(kind=c_char),intent(in) :: name(*), value(*)
          integer(c_int),value :: overwrite
          integer(c_int) :: r
       end function c_setenv
    end interface
    integer(c_int) :: r
    r = c_setenv("ENGINE_MAX_TWOEI_GB"//char(0), trim(s)//char(0), 1_c_int)
end subroutine setenv_max_gb

end program test_coulomb
