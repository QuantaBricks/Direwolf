! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! NMR magnetic shielding tensors with gauge-including atomic orbitals (GIAO), Hartree-Fock / Kohn-Sham, restricted or unrestricted.

module mod_nmr
implicit none
private
public :: nmr_shielding, nmr_isotropy_anisotropy, nmr_xc_giao, nmr_xc_resp
real(8),allocatable,save :: wv_store(:,:)
real(8),save :: t_tick = -1.0d0

contains

subroutine nmr_tick(label)
use omp_lib, only: omp_get_wtime
use MOL_info, only: engine_verbose
implicit none
character(len=*),intent(in) :: label
real(8) :: t
t = omp_get_wtime()
if (engine_verbose >= 2 .and. t_tick >= 0.0d0) print '(A,A,A,F9.2,A)', '[nmr] ', label, '  ', t - t_tick, ' s'
t_tick = t
end subroutine nmr_tick

subroutine nmr_shielding(nConts, nocc, C_a, eLev_a, Pa, Pb, cx, nat, coor_bohr, sdia, spara, niter, dft, ri_j, ub, &
                         nocc_b, C_b, eLev_b, ub_b, nprobe, probe_bohr, sdia_p, spara_p)
use mod_integrals, only: integrals_giao_1e, integrals_nmr_nucleus, integrals_giao_jk, integrals_k_general
use mod_integrals, only: integrals_giao_ri_j, integrals_giao_cosx_k, integrals_build_exchange_cosx_batch_raw
use mod_integrals, only: integrals_build_exchange_cosx_batch_c2g
use mod_integrals, only: cosx_set_resp_grid, cosx_set_resp_relscreen, cosx_set_force_kscreen
use mod_exchange, only: RS_omega, RS_beta, cosx_enabled
use mod_xc, only: xc_uses_tau
use mod_cosmo_state, only: cosmo_enabled, n_tess, tess_coor, q_last
use mod_integrals, only: integrals_giao_charges
implicit none
integer,intent(in) :: nConts, nocc, nat
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts), Pa(nConts,nConts), Pb(nConts,nConts), cx, coor_bohr(3,nat)
real(8),intent(out) :: sdia(3,3,nat), spara(3,3,nat)
integer,intent(out) :: niter
real(8),intent(out),optional :: ub(nConts,nocc,3)
logical,intent(in) :: ri_j
logical,intent(in) :: dft
integer,intent(in),optional :: nocc_b
real(8),intent(in),optional :: C_b(nConts,nConts), eLev_b(nConts)
real(8),intent(out),optional :: ub_b(nConts,*)
integer,intent(in),optional :: nprobe
real(8),intent(in),optional :: probe_bohr(3,*)
real(8),intent(out),optional :: sdia_p(3,3,*), spara_p(3,3,*)
real(8),parameter :: ALPHA = 1.0d0/137.035999084d0, SC = ALPHA*ALPHA*1.0d6, CONV0 = 1.0d-5
integer,parameter :: MAXIT = 100
integer :: nv, x, ia, a, b, k, it, i, nsp, isp, no
real(8),allocatable :: D(:,:), S1(:,:,:), H1(:,:,:), Jx(:,:,:), Kx(:,:,:), Koo(:,:,:), doo(:,:,:,:), Deff(:,:,:)
real(8),allocatable :: Co(:,:), Cv(:,:), rhs(:,:,:), u(:,:,:), r(:,:,:), z(:,:,:), p(:,:,:), Ap(:,:,:), de(:,:), tmp(:,:)
real(8),allocatable :: PSO(:,:,:), A11(:,:,:), A01(:,:,:), dtot(:,:,:), fa(:,:,:), Vxc(:,:,:,:), Cs(:,:), es(:)
real(8),allocatable :: Jlr(:,:,:), Klr(:,:,:), Kc(:,:,:), Kdoo(:,:,:,:), Ksp(:,:,:,:), Jt(:,:,:)
logical :: has_k, meta, use_cosx, kb_cosx, kd_cosx, opn
logical :: kproj_done(2)
real(8),allocatable :: BKs(:,:), BKds(:,:)
character(len=4) :: kbgrid
real(8) :: CONV, wsp
real(8) :: rz(3), rzn(3), alp, bet, T11(3,3), T01(3,3), rmax, sd, sp
integer :: npr

opn = present(C_b)
nsp = merge(2, 1, opn)
allocate(D(nConts,nConts), Deff(nConts,nConts,nsp))
D = Pa + Pb
if (opn) then
   Deff(:,:,1) = 2.0d0*Pa ; Deff(:,:,2) = 2.0d0*Pb
else
   Deff(:,:,1) = D
endif
allocate(S1(nConts,nConts,3), H1(nConts,nConts,3), Jx(nConts,nConts,3), Ksp(nConts,nConts,3,nsp))
call nmr_tick('start')
call integrals_giao_1e(nConts, S1, H1)
call nmr_tick('GIAO one-electron integrals')
if (cosmo_enabled) then
   block
      real(8),allocatable :: H1x(:,:,:)
      allocate(H1x(nConts,nConts,3))
      call integrals_giao_charges(nConts, n_tess, tess_coor, q_last, H1x)
      H1 = H1 + H1x
   end block
endif
meta = dft .and. xc_uses_tau()
has_k = (cx /= 0.0d0 .or. RS_beta /= 0.0d0 .or. meta)
allocate(doo(nConts,nConts,3,nsp), tmp(nConts,nConts))
do isp = 1,nsp
   do x = 1,3
      call dgemm('N','N',nConts,nConts,nConts,1.0d0,Deff(1,1,isp),nConts,S1(1,1,x),nConts,0.0d0,tmp,nConts)
      call dgemm('N','N',nConts,nConts,nConts,-0.5d0,tmp,nConts,Deff(1,1,isp),nConts,0.0d0,doo(1,1,x,isp),nConts)
   enddo
enddo

use_cosx = cosx_enabled .and. (cx /= 0.0d0 .or. RS_beta /= 0.0d0)
kb_cosx = use_cosx
kd_cosx = use_cosx
if (use_cosx) then
   kbgrid = 'hi'
   call cosx_set_resp_grid(kbgrid)
   call cosx_set_resp_relscreen(.true.)
endif
Jx = 0.0d0 ; Ksp = 0.0d0
allocate(Kx(nConts,nConts,3))
if (.not. ri_j .or. (.not. kb_cosx .and. cx /= 0.0d0)) then
   allocate(Jt(nConts,nConts,3))
   do isp = 1,nsp
      call integrals_giao_jk(nConts, Deff(1,1,isp), Jt, Kx, cx /= 0.0d0 .and. .not. kb_cosx, 0.0d0)
      if (.not. ri_j) Jx = Jx + Jt/real(nsp,8)
      if (cx /= 0.0d0 .and. .not. kb_cosx) Ksp(:,:,:,isp) = Kx
   enddo
   deallocate(Jt)
endif
if (ri_j) call integrals_giao_ri_j(nConts, D, Jx)
call nmr_tick('J^B (RI-J or exact J/K)')
if (kb_cosx) allocate(Kdoo(nConts,nConts,3,nsp))
kproj_done = .false.
if (kb_cosx) allocate(BKs(3*nConts*nConts,nsp), BKds(3*nConts*nConts,nsp))
do isp = 1,nsp
   if (kb_cosx) then
   endif
   if (kproj_done(isp)) then
      Ksp(:,:,:,isp) = 0.0d0
      Kdoo(:,:,:,isp) = 0.0d0
   else if (kb_cosx) then
      allocate(Kc(nConts,nConts,3))
      call integrals_giao_cosx_k(nConts, Deff(1,1,isp), nat, coor_bohr, Kc, 3, doo(1,1,1,isp), Kdoo(1,1,1,isp))
      Ksp(:,:,:,isp) = -merge(cx, RS_beta, cx /= 0.0d0)*Kc
      deallocate(Kc)
   else
      Ksp(:,:,:,isp) = cx*Ksp(:,:,:,isp)
      if (RS_omega > 0.0d0 .and. RS_beta /= 0.0d0) then
         allocate(Jlr(nConts,nConts,3), Klr(nConts,nConts,3))
         call integrals_giao_jk(nConts, Deff(1,1,isp), Jlr, Klr, .true., RS_omega)
         Ksp(:,:,:,isp) = Ksp(:,:,:,isp) + RS_beta*Klr
         deallocate(Jlr, Klr)
      endif
   endif
enddo
deallocate(Kx)
call nmr_tick('K^B (COSX or exact)')
if (dft) then
   allocate(Vxc(nConts,nConts,3,nsp))
   if (opn) then
      block
         real(8),allocatable :: Ps(:,:,:)
         allocate(Ps(nConts,nConts,2))
         Ps(:,:,1) = Pa ; Ps(:,:,2) = Pb
         call nmr_xc_giao(nConts, 2, Ps, nat, coor_bohr, Vxc)
      end block
   else
      call nmr_xc_giao(nConts, 1, D, nat, coor_bohr, Vxc)
   endif
endif

call nmr_tick('XC GIAO term')
niter = 0
sdia = 0.0d0 ; spara = 0.0d0
allocate(dtot(nConts,nConts,3), Cs(nConts,nConts), es(nConts))
dtot = 0.0d0
wsp = merge(0.5d0, 1.0d0, opn)
do isp = 1,nsp
   if (isp == 1) then
      no = nocc ; Cs = C_a ; es = eLev_a
   else
      no = nocc_b ; Cs = C_b ; es = eLev_b
   endif
   nv = nConts - no
   if (no == 0 .or. nv == 0) cycle
   allocate(Co(nConts,no), Cv(nConts,nv))
   Co = Cs(:,1:no) ; Cv = Cs(:,no+1:nConts)

   allocate(fa(nConts,nConts,3))
   fa = H1 + Jx - 0.5d0*Ksp(:,:,:,isp)
   if (dft) fa = fa + Vxc(:,:,:,isp)
   if (has_k) then
      allocate(Koo(nConts,nConts,3))
      if (kd_cosx) then
         call k_eff(doo(:,:,:,isp), Koo, Kdoo(:,:,:,isp))
      else
         call k_eff(doo(:,:,:,isp), Koo)
      endif
      fa = fa - 0.5d0*Koo
      deallocate(Koo)
   endif

   allocate(rhs(nv,no,3), de(nv,no), u(nv,no,3), r(nv,no,3), z(nv,no,3), p(nv,no,3), Ap(nv,no,3))
   do i = 1,no
      do a = 1,nv
         de(a,i) = es(no+a) - es(i)
      enddo
   enddo
   block
   real(8),allocatable :: t1(:,:), t2(:,:)
   allocate(t1(nConts,no), t2(nv,no))
   do x = 1,3
      call dgemm('N','N',nConts,no,nConts,1.0d0,fa(1,1,x),nConts,Co,nConts,0.0d0,t1,nConts)
      call dgemm('T','N',nv,no,nConts,1.0d0,Cv,nConts,t1,nConts,0.0d0,rhs(1,1,x),nv)
      call dgemm('N','N',nConts,no,nConts,1.0d0,S1(1,1,x),nConts,Co,nConts,0.0d0,t1,nConts)
      call dgemm('T','N',nv,no,nConts,1.0d0,Cv,nConts,t1,nConts,0.0d0,t2,nv)
      do i = 1,no
         rhs(:,i,x) = -(rhs(:,i,x) - t2(:,i)*es(i))
      enddo
   enddo
   if (kproj_done(isp)) then
      rhs = rhs - 0.5d0*merge(cx, RS_beta, cx /= 0.0d0) &
          *(reshape(BKs(1:nv*no*3,isp), [nv,no,3]) - reshape(BKds(1:nv*no*3,isp), [nv,no,3]))
   endif
   deallocate(t1, t2)
   end block
   deallocate(fa)

   CONV = CONV0
   if (use_cosx) call cosx_set_force_kscreen(.true.)
   if (.not. has_k) then
      do x = 1,3
         u(:,:,x) = rhs(:,:,x)/de
      enddo
   else
      u = 0.0d0
      r = rhs
      do x = 1,3
         z(:,:,x) = r(:,:,x)/de
      enddo
      p = z
      do x = 1,3
         rz(x) = sum(r(:,:,x)*z(:,:,x))
      enddo
      do it = 1,MAXIT
         niter = max(niter, it)
         rmax = 0.0d0
         do x = 1,3
            rmax = max(rmax, maxval(abs(z(:,:,x))))
         enddo
         if (rmax < CONV) exit
         call apply_a(p, Ap)
         do x = 1,3
            alp = rz(x)/sum(p(:,:,x)*Ap(:,:,x))
            u(:,:,x) = u(:,:,x) + alp*p(:,:,x)
            r(:,:,x) = r(:,:,x) - alp*Ap(:,:,x)
            z(:,:,x) = r(:,:,x)/de
            rzn(x) = sum(r(:,:,x)*z(:,:,x))
            bet = rzn(x)/rz(x)
            p(:,:,x) = z(:,:,x) + bet*p(:,:,x)
            rz(x) = rzn(x)
         enddo
      enddo
   endif
   if (use_cosx) call cosx_set_force_kscreen(.false.)

   if (present(ub)) then
      block
         real(8),allocatable :: t1(:,:), uo(:,:,:)
         allocate(t1(nConts,no), uo(nConts,no,3))
         do x = 1,3
            uo(no+1:nConts,:,x) = u(:,:,x)
            call dgemm('N','N',nConts,no,nConts,1.0d0,S1(1,1,x),nConts,Co,nConts,0.0d0,t1,nConts)
            call dgemm('T','N',no,no,nConts,-0.5d0,Co,nConts,t1,nConts,0.0d0,uo(1,1,x),nConts)
         enddo
         if (isp == 1) then
            ub = uo
         else
            do x = 1,3
               do i = 1,no
                  ub_b(:,i+(x-1)*no) = uo(:,i,x)
               enddo
            enddo
         endif
      end block
   else
      block
      real(8),allocatable :: t1(:,:)
      allocate(t1(nConts,no))
      do x = 1,3
         call dgemm('N','N',nConts,no,nv,1.0d0,Cv,nConts,u(1,1,x),nv,0.0d0,t1,nConts)
         call dgemm('N','T',nConts,nConts,no,1.0d0,t1,nConts,Co,nConts,0.0d0,tmp,nConts)
         do b = 1,nConts
            do a = 1,nConts
               dtot(a,b,x) = dtot(a,b,x) + wsp*(2.0d0*(tmp(a,b) - tmp(b,a)) + doo(a,b,x,isp))
            enddo
         enddo
      enddo
      deallocate(t1)
      end block
   endif
   deallocate(Co, Cv, rhs, de, u, r, z, p, Ap)
enddo

if (.not. present(ub)) then
   do isp = 1,nsp
      if (isp == 1) then
         no = nocc
      else
         no = nocc_b
      endif
      if (no > 0 .and. no < nConts) cycle
      dtot = dtot + wsp*doo(:,:,:,isp)
   enddo
   call nmr_tick('response equations (CG, three directions)')
   allocate(PSO(nConts,nConts,3), A11(nConts,nConts,9), A01(nConts,nConts,9))
   npr = 0
   if (present(nprobe)) npr = nprobe
   do ia = 1,nat + npr
      if (ia <= nat) then
         call integrals_nmr_nucleus(nConts, coor_bohr(:,ia), PSO, A11, A01)
      else
         call integrals_nmr_nucleus(nConts, probe_bohr(:,ia-nat), PSO, A11, A01)
      endif
      do k = 1,9
         a = (k-1)/3 + 1 ; b = mod(k-1,3) + 1
         T11(a,b) = sum(A11(:,:,k)*D)
         T01(a,b) = sum(A01(:,:,k)*D)
      enddo
      do b = 1,3
         do a = 1,3
            sd = SC*(T11(a,b) + T01(a,b))
            if (a == b) sd = sd - SC*(T11(1,1) + T11(2,2) + T11(3,3))
            sp = SC*sum(dtot(:,:,a)*PSO(:,:,b))
            if (ia <= nat) then
               sdia(a,b,ia) = sd ; spara(a,b,ia) = sp
            else
               sdia_p(a,b,ia-nat) = sd ; spara_p(a,b,ia-nat) = sp
            endif
         enddo
      enddo
   enddo
endif

call nmr_tick('PSO operators and the shielding sums')
if (use_cosx) then
   call cosx_set_resp_grid('')
   call cosx_set_resp_relscreen(.false.)
endif

contains

   subroutine k_eff(dm, kk, kpre, no_cosx)
   real(8),intent(in) :: dm(nConts,nConts,3)
   real(8),intent(out) :: kk(nConts,nConts,3)
   real(8),intent(in),optional :: kpre(nConts,nConts,3)
   logical,intent(in),optional :: no_cosx
   real(8),allocatable :: kl(:,:,:)
   integer :: kk_i
   kk = 0.0d0
   if (meta) then
      allocate(kl(nConts,nConts,3))
      call nmr_xc_resp(nConts, 3, dm, kl, isp)
      kk = 2.0d0*kl
      deallocate(kl)
   endif
   if (kd_cosx .and. present(no_cosx)) then
      continue
   else if (kd_cosx .and. present(kpre)) then
      kk = kk + merge(cx, RS_beta, cx /= 0.0d0)*kpre
   else if (kd_cosx) then
      allocate(kl(nConts,nConts,3))
      call cosx_set_resp_grid('lo')
      call integrals_build_exchange_cosx_batch_raw(nConts, 3, dm, kl)
      call cosx_set_resp_grid(kbgrid)
      do kk_i = 1,3
         kk(:,:,kk_i) = kk(:,:,kk_i) + merge(cx, RS_beta, cx /= 0.0d0)*0.5d0*(transpose(kl(:,:,kk_i)) - kl(:,:,kk_i))
      enddo
      deallocate(kl)
   else if (cx /= 0.0d0) then
      allocate(kl(nConts,nConts,3))
      call integrals_k_general(nConts, 3, dm, kl, 0.0d0)
      kk = kk + cx*kl
      deallocate(kl)
   endif
   if (.not. kd_cosx .and. RS_omega > 0.0d0 .and. RS_beta /= 0.0d0) then
      allocate(kl(nConts,nConts,3))
      call integrals_k_general(nConts, 3, dm, kl, RS_omega)
      kk = kk + RS_beta*kl
      deallocate(kl)
   endif
   end subroutine k_eff

   subroutine apply_a(pp, app)
   real(8),intent(in) :: pp(nv,no,3)
   real(8),intent(out) :: app(nv,no,3)
   real(8),allocatable :: dd(:,:,:), kk(:,:,:), t1(:,:), t2(:,:), L1(:,:,:), L2(:,:,:), KC(:,:,:), KC2(:,:,:)
   integer :: xx, aa, bb
   allocate(dd(nConts,nConts,3), kk(nConts,nConts,3), t1(nConts,no), t2(nConts,nConts))
   do xx = 1,3
      call dgemm('N','N',nConts,no,nv,1.0d0,Cv,nConts,pp(1,1,xx),nv,0.0d0,t1,nConts)
      call dgemm('N','T',nConts,nConts,no,1.0d0,t1,nConts,Co,nConts,0.0d0,t2,nConts)
      do bb = 1,nConts
         do aa = 1,nConts
            dd(aa,bb,xx) = 2.0d0*(t2(aa,bb) - t2(bb,aa))
         enddo
      enddo
   enddo
   if (kd_cosx) then
      call k_eff(dd, kk, no_cosx=.true.)
   else
      call k_eff(dd, kk)
   endif
   if (kd_cosx) then
      allocate(L1(nConts,no,3), L2(nConts,no,3), KC(nConts,no,3), KC2(nConts,no,3))
      do xx = 1,3
         call dgemm('N','N',nConts,no,nv,2.0d0,Cv,nConts,pp(1,1,xx),nv,0.0d0,L1(1,1,xx),nConts)
         L2(:,:,xx) = -L1(:,:,xx)
      enddo
      call cosx_set_resp_grid('lo')
      call integrals_build_exchange_cosx_batch_c2g(nConts, 3, L1, L2, no, Cs, KC, KC2)
      call cosx_set_resp_grid(kbgrid)
   endif
   do xx = 1,3
      call dgemm('N','N',nConts,no,nConts,1.0d0,kk(1,1,xx),nConts,Co,nConts,0.0d0,t1,nConts)
      if (kd_cosx) t1 = t1 + merge(cx, RS_beta, cx /= 0.0d0)*0.5d0*(KC2(:,:,xx) - KC(:,:,xx))
      call dgemm('T','N',nv,no,nConts,1.0d0,Cv,nConts,t1,nConts,0.0d0,app(1,1,xx),nv)
      app(:,:,xx) = de*pp(:,:,xx) - 0.5d0*app(:,:,xx)
   enddo
   if (kd_cosx) deallocate(L1, L2, KC, KC2)
   deallocate(dd, kk, t1, t2)
   end subroutine apply_a

end subroutine nmr_shielding
subroutine nmr_xc_giao(nConts, nsp, Ds, nat, coor_bohr, Vxc)
use GRID_info, only: Grids, ngrids
use mod_xc, only: xc_eval, xc_uses_tau, xc_uses_vv10
use mod_integrals, only: integrals_ao_atom
use mod_cube, only: cube_screen_begin, cube_screen_end
use omp_lib, only: omp_get_num_threads, omp_get_thread_num
implicit none
integer,intent(in) :: nConts, nat, nsp
real(8),intent(in) :: Ds(nConts,nConts,nsp), coor_bohr(3,nat)
real(8),intent(out) :: Vxc(nConts,nConts,3,nsp)
interface
   subroutine GTOeval_point(coor_grid, val0, val1)
   use MOL_info, only: nconts
   real(8),intent(in) :: coor_grid(3)
   real(8),intent(out) :: val0(nconts), val1(nconts,3)
   end subroutine GTOeval_point
end interface
integer,parameter :: NB = 1024, NBT = 128
real(8),allocatable :: Phi(:,:), dPhi(:,:,:), T(:,:), X(:,:), Q(:,:,:,:), v0(:), v1(:,:), Ql(:,:,:,:), Dm(:,:,:), D(:,:)
real(8),allocatable :: Hcl(:,:,:,:), Aal(:,:,:,:), Bjl(:,:,:,:)
integer :: nblk, ib, itd, isp, ns, j
real(8),allocatable :: Pk(:,:), dPk(:,:,:), Dsub(:,:), Qs(:,:)
integer,allocatable :: idx(:)
real(8),parameter :: ATHR = 1.0d-11
real(8),allocatable :: rho(:,:), gr(:,:,:), sig(:), wt(:), rr(:,:), vr(:), zz(:), ex(:), ec(:), d1r(:), d1s(:)
real(8),allocatable :: vra(:), vrb(:), vsa(:), vsb(:), vsab(:), saa(:), sbb(:), sab(:)
real(8),allocatable :: vrs(:,:), Gs(:,:,:), vts(:,:), reff(:,:)
integer,allocatable :: ao_atom(:)
integer :: ip0, nbn, k, g, a, b, mu, nu
logical :: own, meta
real(8),allocatable :: Ym(:,:), Hc(:,:,:,:), Aa(:,:,:,:), Bj(:,:,:,:), Rao(:,:), tauv(:,:), vtau(:), tb(:), vtb(:)
real(8) :: dd(3), eps3(3,3,3), v1s

meta = xc_uses_tau()
if (.not. allocated(Grids) .or. ngrids <= 0) error stop 'nmr: no xc grid (the functional has to be a density functional)'

allocate(Dm(nConts,nConts,nsp), D(nConts,nConts))
if (nsp == 1) then
   Dm = Ds ; D = Ds(:,:,1)
else
   Dm = 2.0d0*Ds ; D = Ds(:,:,1) + Ds(:,:,2)
endif
allocate(Q(nConts,nConts,3,nsp))
Q = 0.0d0
allocate(ao_atom(nConts))
call integrals_ao_atom(nConts, ao_atom)
if (meta) then
   allocate(Hc(nConts,nConts,3,nsp), Aa(nConts,nConts,3,nsp), Bj(nConts,nConts,3,nsp), Rao(3,nConts))
   Hc = 0.0d0 ; Aa = 0.0d0 ; Bj = 0.0d0
   do nu = 1,nConts
      Rao(:,nu) = coor_bohr(:,ao_atom(nu))
   enddo
   if (allocated(wv_store)) deallocate(wv_store)
   allocate(wv_store(ngrids,nsp))
   wv_store = 0.0d0
endif
call cube_screen_begin(own)
nblk = (ngrids + NBT - 1)/NBT
!$omp parallel private(Phi, dPhi, T, X, rho, gr, sig, wt, rr, vr, zz, ex, ec, d1r, d1s, vra, vrb, vsa, vsb, vsab, &
!$omp&                 saa, sbb, sab, vrs, Gs, vts, reff, isp, &
!$omp&                 tauv, vtau, tb, vtb, v0, v1, Ql, Hcl, Aal, Bjl, ib, ip0, nbn, k, g, nu, itd, &
!$omp&                 Pk, dPk, Dsub, Qs, idx, ns, j)
allocate(Phi(NBT,nConts), dPhi(NBT,nConts,3), T(NBT,nConts), X(NBT,nConts), Ql(nConts,nConts,3,nsp))
allocate(Pk(NBT,nConts), dPk(NBT,nConts,3), Dsub(nConts,nConts), Qs(nConts,nConts), idx(nConts))
allocate(rho(NBT,nsp), gr(NBT,3,nsp), sig(NBT), wt(NBT), rr(NBT,3), vr(NBT), zz(NBT), ex(NBT), ec(NBT), d1r(NBT), d1s(NBT))
allocate(vra(NBT), vrb(NBT), vsa(NBT), vsb(NBT), vsab(NBT), v0(nConts), v1(nConts,3))
allocate(saa(NBT), sbb(NBT), sab(NBT), vrs(NBT,nsp), Gs(NBT,3,nsp), vts(NBT,nsp), reff(NBT,nsp))
Ql = 0.0d0
zz = 0.0d0
vts = 0.0d0
allocate(tauv(NBT,nsp), vtau(NBT), tb(NBT), vtb(NBT))
if (meta) then
   allocate(Hcl(nConts,nConts,3,nsp), Aal(nConts,nConts,3,nsp), Bjl(nConts,nConts,3,nsp))
   Hcl = 0.0d0 ; Aal = 0.0d0 ; Bjl = 0.0d0
endif
!$omp do schedule(static,1)
do ib = 1,nblk
   ip0 = 1 + (ib-1)*NBT
   nbn = min(NBT, ngrids - ip0 + 1)
   do k = 1,nbn
      call GTOeval_point(Grids(ip0+k-1)%coor, v0, v1)
      Phi(k,:) = v0
      dPhi(k,:,1) = v1(:,1) ; dPhi(k,:,2) = v1(:,2) ; dPhi(k,:,3) = v1(:,3)
      wt(k) = Grids(ip0+k-1)%weight
      rr(k,:) = Grids(ip0+k-1)%coor
   enddo
   ns = 0
   do nu = 1,nConts
      if (meta .or. maxval(abs(Phi(1:nbn,nu))) > ATHR .or. maxval(abs(dPhi(1:nbn,nu,1))) > ATHR &
          .or. maxval(abs(dPhi(1:nbn,nu,2))) > ATHR .or. maxval(abs(dPhi(1:nbn,nu,3))) > ATHR) then
         ns = ns + 1
         idx(ns) = nu
      endif
   enddo
   do j = 1,ns
      Pk(:,j) = Phi(:,idx(j))
      dPk(:,j,1) = dPhi(:,idx(j),1) ; dPk(:,j,2) = dPhi(:,idx(j),2) ; dPk(:,j,3) = dPhi(:,idx(j),3)
   enddo
   do isp = 1,nsp
      do j = 1,ns
         Dsub(1:ns,j) = Ds(idx(1:ns),idx(j),isp)
      enddo
      call dgemm('N','N',nbn,ns,ns,1.0d0,Pk,NBT,Dsub,nConts,0.0d0,T,NBT)
      do k = 1,nbn
         rho(k,isp) = max(sum(T(k,1:ns)*Pk(k,1:ns)), 0.0d0)
         do g = 1,3
            gr(k,g,isp) = 2.0d0*sum(T(k,1:ns)*dPk(k,1:ns,g))
         enddo
      enddo
      if (meta) then
         tauv(1:nbn,isp) = 0.0d0
         do g = 1,3
            call dgemm('N','N',nbn,nConts,nConts,1.0d0,dPhi(1,1,g),NBT,Ds(1,1,isp),nConts,0.0d0,X,NBT)
            do k = 1,nbn
               tauv(k,isp) = tauv(k,isp) + 0.5d0*sum(X(k,:)*dPhi(k,:,g))
            enddo
         enddo
      endif
   enddo
   if (nsp == 1) then
      do k = 1,nbn
         sig(k) = gr(k,1,1)**2 + gr(k,2,1)**2 + gr(k,3,1)**2
         vr(k) = 0.5d0*rho(k,1)
      enddo
      if (meta) then
         tb(1:nbn) = 0.5d0*tauv(1:nbn,1)
         call xc_eval(1, 1, nbn, vr(1:nbn), vr(1:nbn), zz(1:nbn), zz(1:nbn), zz(1:nbn), sig(1:nbn), &
                      ex(1:nbn), ec(1:nbn), d1r(1:nbn), d1s(1:nbn), vra(1:nbn), vrb(1:nbn), &
                      vsa(1:nbn), vsb(1:nbn), vsab(1:nbn), tb(1:nbn), tb(1:nbn), vtau(1:nbn), vtb(1:nbn))
         vts(1:nbn,1) = vtau(1:nbn)
      else
         call xc_eval(1, 1, nbn, vr(1:nbn), vr(1:nbn), zz(1:nbn), zz(1:nbn), zz(1:nbn), sig(1:nbn), &
                      ex(1:nbn), ec(1:nbn), d1r(1:nbn), d1s(1:nbn), vra(1:nbn), vrb(1:nbn), &
                      vsa(1:nbn), vsb(1:nbn), vsab(1:nbn))
      endif
      do k = 1,nbn
         vrs(k,1) = d1r(k)
         Gs(k,:,1) = 2.0d0*d1s(k)*gr(k,:,1)
         reff(k,1) = rho(k,1)
      enddo
   else
      do k = 1,nbn
         saa(k) = gr(k,1,1)**2 + gr(k,2,1)**2 + gr(k,3,1)**2
         sbb(k) = gr(k,1,2)**2 + gr(k,2,2)**2 + gr(k,3,2)**2
         sab(k) = gr(k,1,1)*gr(k,1,2) + gr(k,2,1)*gr(k,2,2) + gr(k,3,1)*gr(k,3,2)
         sig(k) = saa(k) + sbb(k) + 2.0d0*sab(k)
      enddo
      if (meta) then
         call xc_eval(2, 1, nbn, rho(1:nbn,1), rho(1:nbn,2), saa(1:nbn), sbb(1:nbn), sab(1:nbn), sig(1:nbn), &
                      ex(1:nbn), ec(1:nbn), d1r(1:nbn), d1s(1:nbn), vra(1:nbn), vrb(1:nbn), &
                      vsa(1:nbn), vsb(1:nbn), vsab(1:nbn), tauv(1:nbn,1), tauv(1:nbn,2), vtau(1:nbn), vtb(1:nbn))
         vts(1:nbn,1) = vtau(1:nbn) ; vts(1:nbn,2) = vtb(1:nbn)
      else
         call xc_eval(2, 1, nbn, rho(1:nbn,1), rho(1:nbn,2), saa(1:nbn), sbb(1:nbn), sab(1:nbn), sig(1:nbn), &
                      ex(1:nbn), ec(1:nbn), d1r(1:nbn), d1s(1:nbn), vra(1:nbn), vrb(1:nbn), &
                      vsa(1:nbn), vsb(1:nbn), vsab(1:nbn))
      endif
      do k = 1,nbn
         vrs(k,1) = vra(k) ; vrs(k,2) = vrb(k)
         Gs(k,:,1) = 2.0d0*vsa(k)*gr(k,:,1) + vsab(k)*gr(k,:,2)
         Gs(k,:,2) = 2.0d0*vsb(k)*gr(k,:,2) + vsab(k)*gr(k,:,1)
         reff(k,1) = 2.0d0*rho(k,1) ; reff(k,2) = 2.0d0*rho(k,2)
      enddo
   endif
   do isp = 1,nsp
      if (meta) then
         do k = 1,nbn
            if (reff(k,isp) > 1.0d-12) wv_store(ip0+k-1,isp) = wt(k)*vts(k,isp)/reff(k,isp)
         enddo
         call nmr_tau_block(nbn, nConts, NBT, Phi, dPhi, rr, wt, vts(1,isp), reff(1,isp), Dm(1,1,isp), Rao, &
                            Hcl(1,1,1,isp), Aal(1,1,1,isp), Bjl(1,1,1,isp), wv_store(ip0,isp))
      endif
      do g = 1,3
         do j = 1,ns
            do k = 1,nbn
               X(k,j) = wt(k)*(0.5d0*vrs(k,isp)*rr(k,g)*Pk(k,j) &
                         + rr(k,g)*(Gs(k,1,isp)*dPk(k,j,1) + Gs(k,2,isp)*dPk(k,j,2) + Gs(k,3,isp)*dPk(k,j,3)) &
                         + 0.5d0*Gs(k,g,isp)*Pk(k,j))
            enddo
         enddo
         call dgemm('T','N',ns,ns,nbn,1.0d0,Pk,NBT,X,NBT,0.0d0,Qs,nConts)
         do j = 1,ns
            Ql(idx(1:ns),idx(j),g,isp) = Ql(idx(1:ns),idx(j),g,isp) + Qs(1:ns,j)
         enddo
      enddo
   enddo
enddo
!$omp end do
do itd = 0,omp_get_num_threads()-1
   if (itd == omp_get_thread_num()) then
      Q = Q + Ql
      if (meta) then
         Hc = Hc + Hcl ; Aa = Aa + Aal ; Bj = Bj + Bjl
      endif
   endif
   !$omp barrier
enddo
deallocate(Pk, dPk, Dsub, Qs, idx)
deallocate(Phi, dPhi, T, X, Ql, rho, gr, sig, wt, rr, vr, zz, ex, ec, d1r, d1s, vra, vrb, vsa, vsb, vsab, v0, v1)
deallocate(saa, sbb, sab, vrs, Gs, vts, reff, tauv, vtau, tb, vtb)
if (meta) deallocate(Hcl, Aal, Bjl)
!$omp end parallel
allocate(Phi(NB,nConts), dPhi(NB,nConts,3), T(NB,nConts), X(NB,nConts))
if (xc_uses_vv10()) then
   block
      use mod_vv10, only: vv10_nlc_grid_retrieve, vv10_evaluate
      real(8),allocatable :: nc(:,:), nw(:), nrho(:), ngam(:), ngr(:,:), neps(:), nFn(:), nFg(:)
      integer :: nnp, kk0
      logical :: built
      real(8) :: Enl
      call vv10_nlc_grid_retrieve(nnp, nc, nw, built)
      if (.not. built) error stop 'nmr: the VV10 grid of the SCF is not available'
      allocate(nrho(nnp), ngam(nnp), ngr(3,nnp), neps(nnp), nFn(nnp), nFg(nnp))
      do ip0 = 1,nnp,NB
         nbn = min(NB, nnp - ip0 + 1)
         !$omp parallel private(v0, v1, k)
         allocate(v0(nConts), v1(nConts,3))
         !$omp do schedule(dynamic,16)
         do k = 1,nbn
            call GTOeval_point(nc(:,ip0+k-1), v0, v1)
            Phi(k,:) = v0
            dPhi(k,:,1) = v1(:,1) ; dPhi(k,:,2) = v1(:,2) ; dPhi(k,:,3) = v1(:,3)
         enddo
         !$omp end do
         deallocate(v0, v1)
         !$omp end parallel
         call dgemm('N','N',nbn,nConts,nConts,1.0d0,Phi,NB,D,nConts,0.0d0,T,NB)
         do k = 1,nbn
            kk0 = ip0+k-1
            nrho(kk0) = max(sum(T(k,:)*Phi(k,:)), 0.0d0)
            do g = 1,3
               ngr(g,kk0) = 2.0d0*sum(T(k,:)*dPhi(k,:,g))
            enddo
            ngam(kk0) = ngr(1,kk0)**2 + ngr(2,kk0)**2 + ngr(3,kk0)**2
         enddo
      enddo
      call vv10_evaluate(nnp, nrho, ngam, nw, nc, Enl, neps, nFn, nFg)
      do ip0 = 1,nnp,NB
         nbn = min(NB, nnp - ip0 + 1)
         !$omp parallel private(v0, v1, k)
         allocate(v0(nConts), v1(nConts,3))
         !$omp do schedule(dynamic,16)
         do k = 1,nbn
            call GTOeval_point(nc(:,ip0+k-1), v0, v1)
            Phi(k,:) = v0
            dPhi(k,:,1) = v1(:,1) ; dPhi(k,:,2) = v1(:,2) ; dPhi(k,:,3) = v1(:,3)
         enddo
         !$omp end do
         deallocate(v0, v1)
         !$omp end parallel
         do g = 1,3
            do nu = 1,nConts
               do k = 1,nbn
                  kk0 = ip0+k-1
                  X(k,nu) = nw(kk0)*(0.5d0*nFn(kk0)*nc(g,kk0)*Phi(k,nu) &
                            + 2.0d0*nFg(kk0)*nc(g,kk0)*(ngr(1,kk0)*dPhi(k,nu,1) + ngr(2,kk0)*dPhi(k,nu,2) &
                            + ngr(3,kk0)*dPhi(k,nu,3)) + nFg(kk0)*ngr(g,kk0)*Phi(k,nu))
               enddo
            enddo
            do isp = 1,nsp
               call dgemm('T','N',nConts,nConts,nbn,1.0d0,Phi,NB,X,NB,1.0d0,Q(1,1,g,isp),nConts)
            enddo
         enddo
      enddo
   end block
endif
call cube_screen_end(own)
do isp = 1,nsp
   do g = 1,3
      do nu = 1,nConts
         do mu = 1,nu
            T(1,1) = Q(mu,nu,g,isp) + Q(nu,mu,g,isp)
            Q(mu,nu,g,isp) = T(1,1) ; Q(nu,mu,g,isp) = T(1,1)
         enddo
      enddo
   enddo
enddo
eps3 = 0.0d0
eps3(1,2,3) = 1.0d0 ; eps3(2,3,1) = 1.0d0 ; eps3(3,1,2) = 1.0d0
eps3(1,3,2) = -1.0d0 ; eps3(3,2,1) = -1.0d0 ; eps3(2,1,3) = -1.0d0
Vxc = 0.0d0
do isp = 1,nsp
   do nu = 1,nConts
      do mu = 1,nConts
         dd = coor_bohr(:,ao_atom(mu)) - coor_bohr(:,ao_atom(nu))
         do a = 1,3
            do b = 1,3
               do g = 1,3
                  if (eps3(a,b,g) /= 0.0d0) Vxc(mu,nu,a,isp) = Vxc(mu,nu,a,isp) - 0.5d0*eps3(a,b,g)*dd(b)*Q(mu,nu,g,isp)
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
if (meta) then
   allocate(Ym(nConts,nConts))
   do isp = 1,nsp
      do a = 1,3
         do nu = 1,nConts
            do mu = 1,nConts
               v1s = 0.0d0
               do b = 1,3
                  do g = 1,3
                     if (eps3(g,a,b) /= 0.0d0) then
                        v1s = v1s + 0.5d0*eps3(g,a,b)*(Rao(b,mu)*Hc(mu,nu,g,isp) - Rao(b,nu)*Hc(nu,mu,g,isp))
                     endif
                  enddo
               enddo
               Ym(mu,nu) = (Aa(mu,nu,a,isp) - Aa(nu,mu,a,isp)) + 0.5d0*v1s + 0.5d0*(Bj(mu,nu,a,isp) - Bj(nu,mu,a,isp))
            enddo
         enddo
         Vxc(:,:,a,isp) = Vxc(:,:,a,isp) - Ym
      enddo
   enddo
   deallocate(Ym, Hc, Aa, Bj, Rao)
endif
deallocate(Phi, dPhi, T, X, Q, Dm, D)
end subroutine nmr_xc_giao

subroutine nmr_tau_block(nbn, nConts, NB, Phi, dPhi, rr, wt, vtau, rho, D, Rao, Hc, Aa, Bj, wv)
integer,intent(in) :: nbn, nConts, NB
real(8),intent(in) :: Phi(NB,nConts), dPhi(NB,nConts,3), rr(NB,3), wt(NB), vtau(NB), rho(NB), D(nConts,nConts), Rao(3,nConts)
real(8),intent(inout) :: Hc(nConts,nConts,3), Aa(nConts,nConts,3), Bj(nConts,nConts,3)
real(8),intent(in) :: wv(*)
real(8),allocatable :: DT(:,:), DG(:,:,:), W1(:,:), W2(:,:), PH(:,:), jx(:,:)
real(8) :: ac(3), e(3)
integer :: k, n, a, c, c1, c2
allocate(DT(NB,nConts), DG(NB,nConts,3), W1(NB,nConts), W2(NB,nConts), PH(NB,nConts), jx(NB,3))
call dgemm('N','N',nbn,nConts,nConts,1.0d0,Phi,NB,D,nConts,0.0d0,DT,NB)
do c = 1,3
   call dgemm('N','N',nbn,nConts,nConts,1.0d0,dPhi(1,1,c),NB,D,nConts,0.0d0,DG(1,1,c),NB)
   do n = 1,nConts
      do k = 1,nbn
         W1(k,n) = wt(k)*vtau(k)*dPhi(k,n,c)
      enddo
   enddo
   call dgemm('T','N',nConts,nConts,nbn,1.0d0,Phi,NB,W1,NB,1.0d0,Hc(1,1,c),nConts)
enddo
do a = 1,3
   e = 0.0d0 ; e(a) = 1.0d0
   do n = 1,nConts
      do k = 1,nbn
         c1 = mod(a,3) + 1 ; c2 = mod(a+1,3) + 1
         PH(k,n) = 0.5d0*(Rao(c1,n)*rr(k,c2) - Rao(c2,n)*rr(k,c1))
      enddo
   enddo
   do c = 1,3
      do n = 1,nConts
         do k = 1,nbn
            W1(k,n) = 0.5d0*wt(k)*vtau(k)*PH(k,n)*dPhi(k,n,c)
         enddo
      enddo
      call dgemm('T','N',nConts,nConts,nbn,1.0d0,W1,NB,dPhi(1,1,c),NB,1.0d0,Aa(1,1,a),nConts)
   enddo
   jx = 0.0d0
   do c = 1,3
      do n = 1,nConts
         c1 = mod(c,3) + 1 ; c2 = mod(c+1,3) + 1
         ac(c) = 0.5d0*(e(c1)*Rao(c2,n) - e(c2)*Rao(c1,n))
         do k = 1,nbn
            jx(k,c) = jx(k,c) + PH(k,n)*Phi(k,n)*DG(k,n,c) - PH(k,n)*DT(k,n)*dPhi(k,n,c) - ac(c)*DT(k,n)*Phi(k,n)
         enddo
      enddo
   enddo
   do n = 1,nConts
      do k = 1,nbn
         W2(k,n) = wv(k)*(jx(k,1)*dPhi(k,n,1) + jx(k,2)*dPhi(k,n,2) + jx(k,3)*dPhi(k,n,3))
      enddo
   enddo
   call dgemm('T','N',nConts,nConts,nbn,1.0d0,Phi,NB,W2,NB,1.0d0,Bj(1,1,a),nConts)
enddo
deallocate(DT, DG, W1, W2, PH, jx)
end subroutine nmr_tau_block

subroutine nmr_xc_resp(nConts, nd, dm, Fr, isp)
use GRID_info, only: Grids, ngrids
use mod_cube, only: cube_screen_begin, cube_screen_end
implicit none
integer,intent(in) :: nConts, nd, isp
real(8),intent(in) :: dm(nConts,nConts,nd)
real(8),intent(out) :: Fr(nConts,nConts,nd)
interface
   subroutine GTOeval_point(coor_grid, val0, val1)
   use MOL_info, only: nconts
   real(8),intent(in) :: coor_grid(3)
   real(8),intent(out) :: val0(nconts), val1(nconts,3)
   end subroutine GTOeval_point
end interface
integer,parameter :: NB = 1024
real(8),allocatable :: Phi(:,:), dPhi(:,:,:), Tm(:,:), W2(:,:), v0(:), v1(:,:), jx(:,:), B(:,:)
integer :: ip0, nbn, k, c, n, id, mu, nu
logical :: own
allocate(Phi(NB,nConts), dPhi(NB,nConts,3), Tm(NB,nConts), W2(NB,nConts), jx(NB,3), B(nConts,nConts))
Fr = 0.0d0
call cube_screen_begin(own)
do ip0 = 1,ngrids,NB
   nbn = min(NB, ngrids - ip0 + 1)
   !$omp parallel private(v0, v1, k)
   allocate(v0(nConts), v1(nConts,3))
   !$omp do schedule(dynamic,16)
   do k = 1,nbn
      call GTOeval_point(Grids(ip0+k-1)%coor, v0, v1)
      Phi(k,:) = v0
      dPhi(k,:,1) = v1(:,1) ; dPhi(k,:,2) = v1(:,2) ; dPhi(k,:,3) = v1(:,3)
   enddo
   !$omp end do
   deallocate(v0, v1)
   !$omp end parallel
   do id = 1,nd
      B = 0.0d0
      do c = 1,3
         call dgemm('N','T',nbn,nConts,nConts,1.0d0,dPhi(1,1,c),NB,dm(1,1,id),nConts,0.0d0,Tm,NB)
         do k = 1,nbn
            jx(k,c) = sum(Phi(k,:)*Tm(k,:))
         enddo
      enddo
      do n = 1,nConts
         do k = 1,nbn
            W2(k,n) = wv_store(ip0+k-1,isp)*(jx(k,1)*dPhi(k,n,1) + jx(k,2)*dPhi(k,n,2) + jx(k,3)*dPhi(k,n,3))
         enddo
      enddo
      call dgemm('T','N',nConts,nConts,nbn,1.0d0,Phi,NB,W2,NB,1.0d0,Fr(1,1,id),nConts)
   enddo
enddo
call cube_screen_end(own)
do id = 1,nd
   B = Fr(:,:,id)
   do nu = 1,nConts
      do mu = 1,nConts
         Fr(mu,nu,id) = 0.5d0*(B(mu,nu) - B(nu,mu))
      enddo
   enddo
enddo
end subroutine nmr_xc_resp

subroutine nmr_isotropy_anisotropy(t, iso, ani, pv)
real(8),intent(in) :: t(3,3)
real(8),intent(out) :: iso, ani, pv(3)
real(8) :: s(3,3), w(3), work(20), tmp
integer :: info, i, j
s = 0.5d0*(t + transpose(t))
call dsyev('N', 'U', 3, s, 3, w, work, 20, info)
iso = (t(1,1) + t(2,2) + t(3,3))/3.0d0
do i = 1,2
   do j = i+1,3
      if (abs(w(j)) < abs(w(i))) then
         tmp = w(i) ; w(i) = w(j) ; w(j) = tmp
      endif
   enddo
enddo
ani = w(3) - 0.5d0*(w(1) + w(2))
pv = w
end subroutine nmr_isotropy_anisotropy

end module mod_nmr
