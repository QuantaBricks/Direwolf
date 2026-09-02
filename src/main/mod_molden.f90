! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Writes Molden-format orbital/geometry output files.

module mod_molden
implicit none
private
public :: write_molden, read_molden

contains

function build_ao_perm(is_sph) result(perm)
use MOL_info, only: Natoms, atoms, nConts
implicit none
logical,intent(in) :: is_sph
integer,allocatable :: perm(:)
integer :: ia, ish, l, ao, nl

allocate(perm(nConts))
ao = 0
do ia = 1,Natoms
   do ish = 1,atoms(ia)%nShell
      l = atoms(ia)%shell(ish)%angMoment
      call shell_perm(ao, l, is_sph, perm)
      nl = merge(2*l+1, (l+1)*(l+2)/2, is_sph)
      ao = ao + nl
   enddo
enddo
end function build_ao_perm

subroutine write_molden(filename)
use MOL_info
use mod_integrals, only: engine_puream
implicit none
INCLUDE 'parameter.h'
character(len=*),intent(in) :: filename
character(len=1),parameter :: shell_letter(0:4) = ['s','p','d','f','g']
integer,allocatable :: perm(:)
integer :: unit_md, ia, ish, ip, l, mo

perm = build_ao_perm(engine_puream)

open(newunit=unit_md, file=trim(filename), status='replace', action='write')

write(unit_md,'(A)') '[Molden Format]'

write(unit_md,'(A)') '[Atoms] AU'
do ia = 1,Natoms
   write(unit_md,'(A2,2I5,3F20.10)') elem_symbol(atoms(ia)%charge), ia, atoms(ia)%charge, &
        atoms(ia)%coor(1)*ans2bohr, atoms(ia)%coor(2)*ans2bohr, atoms(ia)%coor(3)*ans2bohr
enddo

write(unit_md,'(A)') '[GTO]'
do ia = 1,Natoms
   write(unit_md,'(I5,I5)') ia, 0
   do ish = 1,atoms(ia)%nShell
      l = atoms(ia)%shell(ish)%angMoment
      write(unit_md,'(A1,I4,A)') shell_letter(l), atoms(ia)%shell(ish)%nGauss, '  1.00'
      do ip = 1,atoms(ia)%shell(ish)%nGauss
         write(unit_md,'(2E20.10)') atoms(ia)%shell(ish)%exponents(ip), &
                                     atoms(ia)%shell(ish)%contrCoeff(ip)
      enddo
   enddo
   write(unit_md,*)
enddo

if (engine_puream) then
   write(unit_md,'(A)') '[5D]'
   write(unit_md,'(A)') '[7F]'
   write(unit_md,'(A)') '[9G]'
else
   write(unit_md,'(A)') '[6D]'
   write(unit_md,'(A)') '[10F]'
   write(unit_md,'(A)') '[15G]'
endif

write(unit_md,'(A)') '[MO]'
if (Multi .eq. 1) then
   do mo = 1,nConts
      call write_mo_block(unit_md, mo, eLev_a(mo), 'Alpha', &
                           merge(2.0d0,0.0d0,mo.le.n_alpha), C_a(:,mo), perm)
   enddo
else
   do mo = 1,nConts
      call write_mo_block(unit_md, mo, eLev_a(mo), 'Alpha', &
                           merge(1.0d0,0.0d0,mo.le.n_alpha), C_a(:,mo), perm)
   enddo
   do mo = 1,nConts
      call write_mo_block(unit_md, mo, eLev_b(mo), 'Beta ', &
                           merge(1.0d0,0.0d0,mo.le.n_beta), C_b(:,mo), perm)
   enddo
endif

close(unit_md)

end subroutine write_molden

subroutine shell_perm(ao_base, l, is_sph, perm)
implicit none
integer,intent(in) :: ao_base, l
logical,intent(in) :: is_sph
integer,intent(inout) :: perm(:)
integer,parameter :: d2(6)  = [1,4,6,2,3,5]
integer,parameter :: f3(10) = [1,7,10,4,2,3,6,9,8,5]
integer,parameter :: g4(15) = [1,11,15,2,3,7,12,10,14,4,6,13,5,8,9]
integer :: p, k

if (l .lt. 2) then
   do k = 1,2*l+1
      perm(ao_base+k) = ao_base + k
   enddo
   return
endif

if (is_sph) then
   perm(ao_base+1) = ao_base + l+1
   do p = 1,l
      perm(ao_base+2*p)   = ao_base + l+1+p
      perm(ao_base+2*p+1) = ao_base + l+1-p
   enddo
else
   select case (l)
   case (2)
      do k = 1,6
         perm(ao_base+k) = ao_base + d2(k)
      enddo
   case (3)
      do k = 1,10
         perm(ao_base+k) = ao_base + f3(k)
      enddo
   case (4)
      do k = 1,15
         perm(ao_base+k) = ao_base + g4(k)
      enddo
   end select
endif

end subroutine shell_perm

subroutine write_mo_block(unit_md, mo, ene, spin, occ, coeff, perm)
implicit none
integer,intent(in) :: unit_md, mo, perm(:)
real(8),intent(in) :: ene, occ, coeff(:)
character(len=*),intent(in) :: spin
integer :: k

write(unit_md,'(A,I0)')     'Sym= ', mo
write(unit_md,'(A,F16.8)')  'Ene= ', ene
write(unit_md,'(A,A)')      'Spin= ', trim(spin)
write(unit_md,'(A,F10.6)')  'Occup= ', occ
do k = 1,size(perm)
   write(unit_md,'(I5,F20.10)') k, coeff(perm(k))
enddo

end subroutine write_mo_block

subroutine read_molden(filename, Pa_out, Pb_out, success)
use MOL_info, only: Natoms, atoms, nConts
use mod_integrals, only: engine_puream
implicit none
character(len=*),intent(in) :: filename
real(8),intent(out) :: Pa_out(nConts,nConts), Pb_out(nConts,nConts)
logical,intent(out) :: success

integer,allocatable :: perm(:)
integer :: unit_md, ios, ia, ish, ip, l_file, ng_file
integer :: idummy
real(8) :: rdummy1, rdummy2
character(len=4096) :: line, upline
logical :: file_exists, is_sph_file, tag_seen, has_beta
real(8),allocatable :: mo_coeff(:), mo_engine(:)
real(8) :: ene, occ
character(len=16) :: spin, letter_buf
logical :: block_ok

success = .false.
Pa_out = 0.0d0
Pb_out = 0.0d0

inquire(file=trim(filename), exist=file_exists)
if (.not. file_exists) return

open(newunit=unit_md, file=trim(filename), status='old', action='read', iostat=ios)
if (ios .ne. 0) return

if (.not. find_line(unit_md, '[GTO]')) then
   close(unit_md); return
endif

do ia = 1,Natoms
   if (.not. next_nonblank(unit_md, line)) then
      close(unit_md); return
   endif
   read(line, *, iostat=ios) idummy, idummy
   if (ios .ne. 0) then
      close(unit_md); return
   endif
   do ish = 1,atoms(ia)%nShell
      if (.not. next_nonblank(unit_md, line)) then
         close(unit_md); return
      endif
      read(line, *, iostat=ios) letter_buf, ng_file, rdummy1
      if (ios .ne. 0) then
         close(unit_md); return
      endif
      l_file = letter_to_l(letter_buf)
      if (l_file .lt. 0 .or. l_file .ne. atoms(ia)%shell(ish)%angMoment &
          .or. ng_file .ne. atoms(ia)%shell(ish)%nGauss) then
         close(unit_md); return
      endif
      do ip = 1,ng_file
         read(unit_md, '(A)', iostat=ios) line
         if (ios .ne. 0) then
            close(unit_md); return
         endif
         read(line, *, iostat=ios) rdummy1, rdummy2
         if (ios .ne. 0) then
            close(unit_md); return
         endif
      enddo
   enddo
enddo

is_sph_file = .false.
tag_seen = .false.
do
   read(unit_md, '(A)', iostat=ios) line
   if (ios .ne. 0) then
      close(unit_md); return
   endif
   upline = to_upper(line)
   if (index(upline, '[MO]') .gt. 0) exit
   if (index(upline, '5D') .gt. 0) then
      is_sph_file = .true.; tag_seen = .true.
   else if (index(upline, '6D') .gt. 0) then
      is_sph_file = .false.; tag_seen = .true.
   endif
enddo
if (.not. tag_seen .or. is_sph_file .neqv. engine_puream) then
   close(unit_md); return
endif

perm = build_ao_perm(engine_puream)
allocate(mo_coeff(nConts), mo_engine(nConts))
has_beta = .false.

do
   if (.not. next_nonblank(unit_md, line)) exit
   upline = to_upper(line)
   if (index(upline, 'SYM=') .gt. 0) then
      if (.not. next_nonblank(unit_md, line)) exit
      upline = to_upper(line)
   endif
   if (index(upline, 'ENE=') .eq. 0) then
      close(unit_md); success = .false.; return
   endif
   ene = parse_real_after_eq(line)
   if (.not. next_nonblank(unit_md, line)) then
      close(unit_md); return
   endif
   upline = to_upper(line)
   if (index(upline, 'SPIN=') .eq. 0) then
      close(unit_md); return
   endif
   spin = adjustl(line(index(upline,'SPIN=')+5:))
   if (.not. next_nonblank(unit_md, line)) then
      close(unit_md); return
   endif
   upline = to_upper(line)
   if (index(upline, 'OCCUP=') .eq. 0) then
      close(unit_md); return
   endif
   occ = parse_real_after_eq(line)

   block_ok = .true.
   do ip = 1,nConts
      read(unit_md, '(A)', iostat=ios) line
      if (ios .ne. 0) then
         block_ok = .false.; exit
      endif
      read(line, *, iostat=ios) idummy, mo_coeff(ip)
      if (ios .ne. 0) then
         block_ok = .false.; exit
      endif
   enddo
   if (.not. block_ok) then
      close(unit_md); return
   endif

   do ip = 1,nConts
      mo_engine(perm(ip)) = mo_coeff(ip)
   enddo

   upline = to_upper(spin)
   if (index(upline, 'BETA') .gt. 0) then
      has_beta = .true.
      if (occ .gt. 0.0d0) Pb_out = Pb_out + occ * outer(mo_engine)
   else
      if (occ .gt. 0.0d0) Pa_out = Pa_out + occ * outer(mo_engine)
   endif
enddo

close(unit_md)

if (.not. has_beta) then
   Pa_out = 0.5d0 * Pa_out
   Pb_out = Pa_out
endif

success = .true.

contains

   function outer(v) result(m)
   implicit none
   real(8),intent(in) :: v(:)
   real(8) :: m(size(v),size(v))
   m = spread(v,2,size(v)) * spread(v,1,size(v))
   end function outer

end subroutine read_molden

integer function letter_to_l(word)
implicit none
character(len=*),intent(in) :: word
character(len=1) :: c
c = to_upper(trim(adjustl(word)))
select case (c)
case ('S'); letter_to_l = 0
case ('P'); letter_to_l = 1
case ('D'); letter_to_l = 2
case ('F'); letter_to_l = 3
case ('G'); letter_to_l = 4
case default; letter_to_l = -1
end select
end function letter_to_l

logical function next_nonblank(unit_md, line)
implicit none
integer,intent(in) :: unit_md
character(len=*),intent(out) :: line
integer :: ios
do
   read(unit_md, '(A)', iostat=ios) line
   if (ios .ne. 0) then
      next_nonblank = .false.; return
   endif
   if (len_trim(line) .gt. 0) then
      next_nonblank = .true.; return
   endif
enddo
end function next_nonblank

logical function find_line(unit_md, tag)
implicit none
integer,intent(in) :: unit_md
character(len=*),intent(in) :: tag
character(len=4096) :: line
integer :: ios
do
   read(unit_md, '(A)', iostat=ios) line
   if (ios .ne. 0) then
      find_line = .false.; return
   endif
   if (trim(to_upper(adjustl(line))) .eq. trim(to_upper(tag))) then
      find_line = .true.; return
   endif
enddo
end function find_line

real(8) function parse_real_after_eq(line)
implicit none
character(len=*),intent(in) :: line
integer :: ieq, ios
ieq = index(line, '=')
parse_real_after_eq = 0.0d0
if (ieq .eq. 0) return
read(line(ieq+1:), *, iostat=ios) parse_real_after_eq
if (ios .ne. 0) parse_real_after_eq = 0.0d0
end function parse_real_after_eq

function to_upper(s) result(u)
implicit none
character(len=*),intent(in) :: s
character(len=len(s)) :: u
integer :: i, ic
u = s
do i = 1,len(s)
   ic = iachar(s(i:i))
   if (ic .ge. iachar('a') .and. ic .le. iachar('z')) u(i:i) = achar(ic - 32)
enddo
end function to_upper

character(len=2) function elem_symbol(z)
implicit none
integer,intent(in) :: z
character(len=2),parameter :: syms(103) = [character(len=2) :: &
   'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
   'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
   'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
   'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
   'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
   'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
   'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
   'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
   'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
   'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
   'Md','No','Lr']
if (z .ge. 1 .and. z .le. 103) then
   elem_symbol = syms(z)
else
   elem_symbol = '? '
endif
end function elem_symbol

end module mod_molden
