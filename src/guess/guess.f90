! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! read basis set and integral

subroutine guess(info,Gtype,Pa_chk,Pb_chk)

use MOL_info
use GRID_info, only: xcgrid_dynamic
use mod_vv10, only: vv10_active_now, vv10_dynamic
use mod_integrals, only: near_singular_overlap
    implicit none
INCLUDE 'parameter.h'
    integer    :: info,Gtype
    real(8)    :: Pa_chk(nconts,nconts),Pb_chk(nconts,nconts)
    integer    :: i,j,k,l,m,n
    real(8)    :: Exc_dummy,E_dummy,dft_dt_dummy

    allocate(C_a(nconts,nconts))
    allocate(C_b(nconts,nconts))

    allocate(Fa(nConts,nConts))
    allocate(Fb(nConts,nConts))

    allocate(Pa(nConts,nConts))
    allocate(Pb(nConts,nConts))
    allocate(eLev_a(nConts),eLev_b(nConts))
    Pa = 0
    Pb = 0
    if (Gtype .eq. 1) then
      print *,"pseudo Huckel initial guess"
      do i = 1,nconts
         do j = i,nconts
            if (i .eq. j)  then
                Fa(i,j) = Hcore(i,j)
            else
                Fa(i,j) = 0.5*1.75*(Hcore(i,i)+Hcore(j,j)) *S(i,j)
                Fa(j,i) = Fa(i,j)
            endif
         enddo
      enddo
      Fb = Fa
      return
    endif
    if (Gtype .eq. 2) then
        print *,"Core Hamiltonian initial guess"
        Fa=Hcore
        Fb=Hcore
    endif
    if (Gtype .eq. 3) then
       if (engine_verbose .ge. 2) print *,"SAD-like atomic-density initial guess"
       call guess_sad(info)
    endif
    if (Gtype .eq. 4) then
       print *,"Checkpoint-file warm-start initial guess"
       Pa = Pa_chk
       Pb = Pb_chk
       Fa = 0.0d0
       Fb = 0.0d0
       C_a = 0.0d0
       C_b = 0.0d0
       if (xcgrid_dynamic) call xcgrid_refine()
       if (vv10_dynamic) vv10_active_now = .true.
       call scf_build_fock(Pa, Pb, .false., 0, 0.0d0, Exc_dummy, E_dummy, dft_dt_dummy, .false.)
       call solHFR_KS(info, 0.0d0, near_singular_overlap)
       call scf_build_fock(Pa, Pb, .false., 1, 0.0d0, Exc_dummy, E_dummy, dft_dt_dummy, .false.)
    endif
end subroutine

subroutine guess_sad(info)
use MOL_info
use mod_integrals, only: compute_atomic_density, engine_harris_guess
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: info
real(8),allocatable :: Fa_gwh(:,:)
real(8),allocatable :: Ablk(:,:),Sblk(:,:),Wblk(:),WORKblk(:)
real(8),allocatable :: Pa_sad(:,:),Pb_sad(:,:)
integer :: ia,off,cnt,ii,jj,p,q,LWORKblk,INFOblk
integer :: Z_atom
real(8) :: Exc_dummy,E_dummy,dft_dt_dummy
real(8) :: remaining_elec,group_elec,occ_per_orbital,occ_a_per_orbital,occ_b_per_orbital
integer :: gstart,gend,k
real(8),parameter :: DEGEN_TOL = 1.0d-4

type :: atomcache_t
   integer :: Z
   character(len=30) :: base
   integer :: cnt
   real(8),allocatable :: Pa(:,:), Pb(:,:)
end type atomcache_t
type(atomcache_t) :: acache(natoms)
integer :: nacache, ci
real(8),allocatable :: Pa_atom_tmp(:,:), Pb_atom_tmp(:,:)
logical :: got_real
logical :: try_harris
character(len=8) :: envval

try_harris = engine_harris_guess
call get_environment_variable("ENGINE_HARRIS_GUESS", envval)
if (trim(envval) .eq. "0") try_harris = .false.
if (trim(envval) .eq. "1") try_harris = .true.

nacache = 0

allocate(Fa_gwh(nconts,nconts))
do ii = 1,nconts
   do jj = ii,nconts
      if (ii .eq. jj) then
         Fa_gwh(ii,jj) = Hcore(ii,jj)
      else
         Fa_gwh(ii,jj) = 0.5d0*1.75d0*(Hcore(ii,ii)+Hcore(jj,jj))*S(ii,jj)
         Fa_gwh(jj,ii) = Fa_gwh(ii,jj)
      endif
   enddo
enddo

allocate(Pa_sad(nconts,nconts),Pb_sad(nconts,nconts))
Pa_sad = 0.0d0
Pb_sad = 0.0d0

off = 0
do ia = 1,natoms
   cnt = atoms(ia)%nconts
   Z_atom = atoms(ia)%charge

   got_real = .false.
   if (try_harris) then
   do ci = 1,nacache
      if (acache(ci)%Z .eq. Z_atom .and. trim(acache(ci)%base) .eq. trim(atoms(ia)%base) &
          .and. acache(ci)%cnt .eq. cnt) then
         Pa_sad(off+1:off+cnt,off+1:off+cnt) = acache(ci)%Pa
         Pb_sad(off+1:off+cnt,off+1:off+cnt) = acache(ci)%Pb
         got_real = .true.
         exit
      endif
   enddo
   if (.not. got_real) then
      call compute_atomic_density(ia, (multi .eq. 1), Pa_atom_tmp, Pb_atom_tmp, got_real)
      if (got_real) then
         Pa_sad(off+1:off+cnt,off+1:off+cnt) = Pa_atom_tmp
         Pb_sad(off+1:off+cnt,off+1:off+cnt) = Pb_atom_tmp
         nacache = nacache + 1
         acache(nacache)%Z = Z_atom
         acache(nacache)%base = atoms(ia)%base
         acache(nacache)%cnt = cnt
         call move_alloc(Pa_atom_tmp, acache(nacache)%Pa)
         call move_alloc(Pb_atom_tmp, acache(nacache)%Pb)
      endif
   endif
   endif

   if (.not. got_real) then
   allocate(Ablk(cnt,cnt),Sblk(cnt,cnt),Wblk(cnt))
   Ablk = Fa_gwh(off+1:off+cnt, off+1:off+cnt)
   Sblk = S(off+1:off+cnt, off+1:off+cnt)
   LWORKblk = 1 + 6*cnt + cnt**2
   allocate(WORKblk(LWORKblk))
   call DSYGV(1,'V','U',cnt,Ablk,cnt,Sblk,cnt,Wblk,WORKblk,LWORKblk,INFOblk)
   remaining_elec = real(Z_atom,8)
   gstart = 1
   do while (gstart .le. cnt .and. remaining_elec .gt. 1.0d-12)
      gend = gstart
      do while (gend .lt. cnt)
         if (abs(Wblk(gend+1)-Wblk(gstart)) .gt. DEGEN_TOL*max(1.0d0,abs(Wblk(gstart)))) exit
         gend = gend + 1
      enddo
      group_elec = min(remaining_elec, real(2*(gend-gstart+1),8))
      if (multi .eq. 1) then
         occ_per_orbital = 0.5d0*group_elec/real(gend-gstart+1,8)
         occ_a_per_orbital = occ_per_orbital
         occ_b_per_orbital = occ_per_orbital
      else
         occ_a_per_orbital = min(1.0d0, group_elec/real(gend-gstart+1,8))
         occ_b_per_orbital = max(0.0d0, group_elec/real(gend-gstart+1,8) - 1.0d0)
      endif
      do k = gstart,gend
         do p = 1,cnt
            do q = 1,cnt
               Pa_sad(off+p,off+q) = Pa_sad(off+p,off+q) + occ_a_per_orbital*Ablk(p,k)*Ablk(q,k)
               Pb_sad(off+p,off+q) = Pb_sad(off+p,off+q) + occ_b_per_orbital*Ablk(p,k)*Ablk(q,k)
            enddo
         enddo
      enddo
      remaining_elec = remaining_elec - group_elec
      gstart = gend + 1
   enddo
   deallocate(Ablk,Sblk,Wblk,WORKblk)
   endif
   off = off + cnt
enddo
deallocate(Fa_gwh)

Pa = Pa_sad
Pb = Pb_sad
Fa = 0.0d0
Fb = 0.0d0
C_a = 0.0d0
C_b = 0.0d0
call scf_build_fock(Pa_sad, Pb_sad, .false., 0, 0.0d0, Exc_dummy, E_dummy, dft_dt_dummy, .false.)
deallocate(Pa_sad,Pb_sad)

end subroutine guess_sad
