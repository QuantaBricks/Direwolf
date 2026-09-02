! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_basis_files: maps a "<Z>-<label>" data/bases block's label to the

module mod_basis_files
implicit none
private
public :: basis_file_for_label, set_basis_dir_override

integer,parameter :: NFAMILY = 11
character(len=20),parameter :: FAMILY_FILES(NFAMILY) = [character(len=20) :: &
   "631g.bas","6311g.bas","321g.bas","ccpv.bas","ccpv_jkfit.bas", &
   "def2.bas","mtzvp.bas","vdzp.bas","lanl2dz.bas","lanl08.bas","sdd.bas"]

character(len=512), save :: basis_dir = ''
logical, save :: basis_dir_ready = .false.

interface
   function engine_get_exe_path(buf, bufsize) bind(C, name="engine_get_exe_path") result(n)
   use, intrinsic :: iso_c_binding, only: c_char, c_long
   character(kind=c_char), intent(out) :: buf(*)
   integer(c_long), value :: bufsize
   integer(c_long) :: n
   end function
end interface

contains

function basis_file_for_label(label) result(path)
character(len=*),intent(in) :: label
character(len=200) :: path
character(len=30) :: token
integer :: ifile, unit_try, ios
logical :: found

if (.not. basis_dir_ready) call resolve_basis_dir()

path = ""
do ifile = 1, NFAMILY
   open(newunit=unit_try, file=trim(basis_dir)//trim(FAMILY_FILES(ifile)), &
        status="old", action="read", iostat=ios)
   if (ios /= 0) cycle

   found = .false.
   do
      read(unit_try,*,iostat=ios) token
      if (ios /= 0) exit
      if (trim(token) == trim(label)) then
         found = .true.
         exit
      endif
   enddo
   close(unit_try)

   if (found) then
      path = trim(basis_dir)//trim(FAMILY_FILES(ifile))
      return
   endif
enddo
end function basis_file_for_label

subroutine set_basis_dir_override(dir)
character(len=*), intent(in) :: dir
if (len_trim(dir) == 0) return
basis_dir = trim(dir)//"/"
basis_dir_ready = .true.
end subroutine set_basis_dir_override

subroutine resolve_basis_dir()
use, intrinsic :: iso_c_binding, only: c_long
character(len=512) :: exe_path, candidate, env_dir
integer(c_long) :: n
integer :: slash, ios, env_len

basis_dir_ready = .true.
basis_dir = "data/bases/"

call get_environment_variable("ENGINE_DATA_DIR", env_dir, env_len, ios)
if (ios == 0 .and. env_len > 0) then
   candidate = trim(env_dir)//"/bases/"
   if (verify_candidate(candidate)) then
      basis_dir = candidate
      return
   endif
endif

n = engine_get_exe_path(exe_path, int(len(exe_path), c_long))
if (n > 0) then
   slash = index(exe_path(1:int(n)), '/', back=.true.)
   if (slash > 0) then
      candidate = exe_path(1:slash)//"data/bases/"
      if (verify_candidate(candidate)) basis_dir = candidate
   endif
endif
end subroutine resolve_basis_dir

logical function verify_candidate(candidate)
character(len=*), intent(in) :: candidate
integer :: ios
inquire(file=trim(candidate)//"631g.bas", exist=verify_candidate, iostat=ios)
if (ios /= 0) verify_candidate = .false.
end function verify_candidate

end module mod_basis_files
