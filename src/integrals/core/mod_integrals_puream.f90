! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! libcint-backed 1-electron integral wrappers (overlap/kinetic/nuclear) plus spherical-harmonic (puream) transform helpers.

submodule (mod_integrals) puream_impl
implicit none
contains

module function cgto_engine(ish, bas) result(n)
implicit none
integer,intent(in) :: ish
integer,intent(in) :: bas(*)
integer :: n
integer,external :: CINTcgto_cart, CINTcgto_spheric
if (engine_puream) then
   n = CINTcgto_spheric(ish, bas)
else
   n = CINTcgto_cart(ish, bas)
endif
end function cgto_engine

module subroutine ovlp1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine ovlp1e_engine

module subroutine kin1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_kin_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine kin1e_engine

module subroutine nuc1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine nuc1e_engine

module subroutine r1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_r_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_r_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine r1e_engine

module subroutine ecp1e_engine(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
implicit none
real(8) :: buf(*)
integer :: dims(2)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
integer(8),intent(in) :: cache
external :: ECPscalar_cart, ECPscalar_sph
if (engine_puream) then
   call ECPscalar_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
else
   call ECPscalar_cart(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
endif
end subroutine ecp1e_engine

module subroutine ovlp1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipovlp_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipovlp_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine ovlp1e_ip_engine

module subroutine kin1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipkin_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipkin_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine kin1e_ip_engine

module subroutine nuc1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipnuc_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipnuc_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine nuc1e_ip_engine

module subroutine rinv1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_iprinv_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_iprinv_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine rinv1e_ip_engine

module subroutine grids1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_grids_ip_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_grids_ip_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine grids1e_ip_engine

module subroutine grids1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
use iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
interface
   function cint1e_grids_sph_c(out,shls,atm,natm,bas,nbas,env,copt) &
            bind(c,name="cint1e_grids_sph") result(r)
   import :: c_int, c_double, c_ptr
   integer(c_int) :: r
   real(c_double) :: out(*)
   integer(c_int) :: shls(*)
   integer(c_int) :: atm(*)
   integer(c_int),value :: natm
   integer(c_int) :: bas(*)
   integer(c_int),value :: nbas
   real(c_double) :: env(*)
   type(c_ptr),value :: copt
   end function cint1e_grids_sph_c
   function cint1e_grids_cart_c(out,shls,atm,natm,bas,nbas,env,copt) &
            bind(c,name="cint1e_grids_cart") result(r)
   import :: c_int, c_double, c_ptr
   integer(c_int) :: r
   real(c_double) :: out(*)
   integer(c_int) :: shls(*)
   integer(c_int) :: atm(*)
   integer(c_int),value :: natm
   integer(c_int) :: bas(*)
   integer(c_int),value :: nbas
   real(c_double) :: env(*)
   type(c_ptr),value :: copt
   end function cint1e_grids_cart_c
end interface
integer(c_int) :: rc
if (engine_puream) then
   rc = cint1e_grids_sph_c(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
else
   rc = cint1e_grids_cart_c(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
endif
end subroutine grids1e_engine

module subroutine grids1e_engine_cached(buf, shls, atm, natm, bas, nbas, env, opt, cache)
use iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
real(8) :: cache(*)
interface
   function int1e_grids_sph_cached_c(out,dims,shls,atm,natm,bas,nbas,env,copt,ccache) &
            bind(c,name="int1e_grids_sph") result(r)
   import :: c_int, c_double, c_ptr
   integer(c_int) :: r
   real(c_double) :: out(*)
   type(c_ptr),value :: dims
   integer(c_int) :: shls(*)
   integer(c_int) :: atm(*)
   integer(c_int),value :: natm
   integer(c_int) :: bas(*)
   integer(c_int),value :: nbas
   real(c_double) :: env(*)
   type(c_ptr),value :: copt
   real(c_double) :: ccache(*)
   end function int1e_grids_sph_cached_c
   function int1e_grids_cart_cached_c(out,dims,shls,atm,natm,bas,nbas,env,copt,ccache) &
            bind(c,name="int1e_grids_cart") result(r)
   import :: c_int, c_double, c_ptr
   integer(c_int) :: r
   real(c_double) :: out(*)
   type(c_ptr),value :: dims
   integer(c_int) :: shls(*)
   integer(c_int) :: atm(*)
   integer(c_int),value :: natm
   integer(c_int) :: bas(*)
   integer(c_int),value :: nbas
   real(c_double) :: env(*)
   type(c_ptr),value :: copt
   real(c_double) :: ccache(*)
   end function int1e_grids_cart_cached_c
end interface
integer(c_int) :: rc
if (engine_puream) then
   rc = int1e_grids_sph_cached_c(buf, c_null_ptr, shls, atm, natm, bas, nbas, env, c_null_ptr, cache)
else
   rc = int1e_grids_cart_cached_c(buf, c_null_ptr, shls, atm, natm, bas, nbas, env, c_null_ptr, cache)
endif
end subroutine grids1e_engine_cached

module subroutine twoe_ip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine twoe_ip1_engine

module subroutine twoe_ip1_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine twoe_ip1_optimizer_engine

module subroutine ovlp1e_ipip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipipovlp_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipipovlp_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine ovlp1e_ipip_engine

module subroutine ovlp1e_ipovlpip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipovlpip_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipovlpip_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine ovlp1e_ipovlpip_engine

module subroutine kin1e_ipip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipipkin_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipipkin_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine kin1e_ipip_engine

module subroutine kin1e_ipkinip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipkinip_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipkinip_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine kin1e_ipkinip_engine

module subroutine nuc1e_ipip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipipnuc_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipipnuc_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine nuc1e_ipip_engine

module subroutine nuc1e_ipnucip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipnucip_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipnucip_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine nuc1e_ipnucip_engine

module subroutine twoe_ipip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint2e_ipip1_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2e_ipip1_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine twoe_ipip1_engine

module subroutine twoe_ip1ip2_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint2e_ip1ip2_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2e_ip1ip2_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine twoe_ip1ip2_engine

module subroutine twoe_ipvip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint2e_ipvip1_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2e_ipvip1_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine twoe_ipvip1_engine

module subroutine twoe_ipip1_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint2e_ipip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint2e_ipip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine twoe_ipip1_optimizer_engine

module subroutine rinv1e_ipiprinv_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_ipiprinv_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_ipiprinv_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine rinv1e_ipiprinv_engine

module subroutine rinv1e_iprinvip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint1e_iprinvip_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint1e_iprinvip_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine rinv1e_iprinvip_engine

module subroutine threec2e_ip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint3c2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint3c2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine threec2e_ip1_engine

module subroutine threec2e_ip2_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint3c2e_ip2_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint3c2e_ip2_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine threec2e_ip2_engine

module subroutine threec2e_ip1_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint3c2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint3c2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine threec2e_ip1_optimizer_engine

module subroutine threec2e_ip2_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint3c2e_ip2_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint3c2e_ip2_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine threec2e_ip2_optimizer_engine

module subroutine twoc2e_ip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint2c2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2c2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine twoc2e_ip1_engine

module subroutine twoc2e_ip2_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint2c2e_ip2_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2c2e_ip2_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine twoc2e_ip2_engine

module subroutine twoc2e_ip1_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint2c2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint2c2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine twoc2e_ip1_optimizer_engine

module subroutine twoc2e_ip2_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint2c2e_ip2_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint2c2e_ip2_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine twoc2e_ip2_optimizer_engine

module subroutine ecp1e_ipnuc_engine(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
implicit none
real(8) :: buf(*)
integer :: dims(2)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
integer(8),intent(in) :: cache
external :: ECPscalar_ipnuc_cart, ECPscalar_ipnuc_sph
if (engine_puream) then
   call ECPscalar_ipnuc_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
else
   call ECPscalar_ipnuc_cart(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
endif
end subroutine ecp1e_ipnuc_engine

module subroutine ecp1e_iprinv_engine(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
implicit none
real(8) :: buf(*)
integer :: dims(2)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
integer(8),intent(in) :: cache
external :: ECPscalar_iprinv_cart, ECPscalar_iprinv_sph
if (engine_puream) then
   call ECPscalar_iprinv_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
else
   call ECPscalar_iprinv_cart(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
endif
end subroutine ecp1e_iprinv_engine

module subroutine twoe_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2e_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine twoe_engine

module subroutine threec2e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint3c2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint3c2e_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine threec2e_engine

module subroutine twoc2e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
if (engine_puream) then
   call cint2c2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
else
   call cint2c2e_cart(buf, shls, atm, natm, bas, nbas, env, opt)
endif
end subroutine twoc2e_engine

module subroutine twoe_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine twoe_optimizer_engine

module subroutine threec2e_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint3c2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint3c2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine threec2e_optimizer_engine

module subroutine twoc2e_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
if (engine_puream) then
   call cint2c2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
else
   call cint2c2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
endif
end subroutine twoc2e_optimizer_engine

end submodule puream_impl
