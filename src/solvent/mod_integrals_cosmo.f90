! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_integrals_cosmo: the ONE integral routine COSMO needs from the

submodule (mod_integrals) cosmo_integrals_impl
implicit none
contains

module subroutine cosmo_build_one_tess_matrix(tess_pos_bohr, Bk)
use MOL_info, only: Natoms
implicit none
real(8),intent(in) :: tess_pos_bohr(3)
real(8),intent(out) :: Bk(:,:)

integer,allocatable :: atm_k(:,:)
real(8),allocatable :: env_k(:)
real(8),allocatable :: buf1eV(:,:)
integer :: nAtoms_k, offpoint_k
integer :: i,j,di,dj,n,shls(4)

nAtoms_k = Natoms + 1
allocate(atm_k(6,nAtoms_k))
atm_k = 0
atm_k(:,1:Natoms) = atm(:,1:Natoms)
atm_k(1,1:Natoms) = 0

allocate(env_k(size(env)+4))
env_k = 0.0d0
env_k(1:size(env)) = env(1:size(env))
offpoint_k = size(env)
atm_k(2,nAtoms_k) = offpoint_k
env_k(offpoint_k+1) = tess_pos_bohr(1)
env_k(offpoint_k+2) = tess_pos_bohr(2)
env_k(offpoint_k+3) = tess_pos_bohr(3)
offpoint_k = offpoint_k + 3
atm_k(3,nAtoms_k) = 3
atm_k(5,nAtoms_k) = offpoint_k
env_k(offpoint_k+1) = 1.0d0
offpoint_k = offpoint_k + 1

n = 0
do i = 0,nBases-1
   n = max(n, cgto_engine(i,bas))
enddo
allocate(buf1eV(n,n))
Bk = 0.0d0
do i = 0,nBases-1
   do j = i,nBases-1
      shls(1) = i
      shls(2) = j
      di = cgto_engine(i,bas)
      dj = cgto_engine(j,bas)
      call nuc1e_engine(buf1eV(1:di,1:dj), shls, atm_k, nAtoms_k, bas, nBases, env_k, 0_8)
      call store1e(shls, di, dj, bas, nBases, buf1eV(1:di,1:dj), NorVEC, size(Bk,1), Bk)
   enddo
enddo

deallocate(buf1eV, atm_k, env_k)
end subroutine cosmo_build_one_tess_matrix

module subroutine esp_at_grid_batch(npts, grid_bohr, Ptot, V_elec)
use MOL_info, only: Natoms
implicit none
integer,intent(in) :: npts
real(8),intent(in) :: grid_bohr(3,npts)
real(8),intent(in) :: Ptot(:,:)
real(8),intent(out) :: V_elec(npts)

integer,parameter :: CHUNK = 256
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),allocatable :: env_g(:), buf(:), int_cache(:)
integer,allocatable :: ao_offset(:)
integer :: nenv, off_grids, ip, ish, jsh, di, dj, is, js, nmax
integer :: g0, g1, ng, g, p, q, e1, e2, shls(4)
real(8) :: w, nij

nenv = size(env)
off_grids = nenv
allocate(env_g(nenv + 3*npts))
env_g = 0.0d0
env_g(1:nenv) = env(1:nenv)
do ip = 1,npts
   env_g(off_grids + 3*(ip-1) + 1) = grid_bohr(1,ip)
   env_g(off_grids + 3*(ip-1) + 2) = grid_bohr(2,ip)
   env_g(off_grids + 3*(ip-1) + 3) = grid_bohr(3,ip)
enddo
env_g(12) = dble(npts)
env_g(13) = dble(off_grids)

allocate(ao_offset(0:nBases-1))
ip = 0
nmax = 0
do ish = 0,nBases-1
   ao_offset(ish) = ip
   ip = ip + cgto_engine(ish,bas)
   nmax = max(nmax, cgto_engine(ish,bas))
enddo

!$omp parallel private(buf,int_cache,g0,g1,ng,ish,jsh,di,dj,is,js,shls,g,p,q,e1,e2,w,nij)
allocate(buf(nmax*nmax*CHUNK), int_cache(INT_CACHE_SIZE))
!$omp do schedule(dynamic)
do g0 = 1,npts,CHUNK
   g1 = min(g0+CHUNK-1, npts)
   ng = g1 - g0 + 1
   V_elec(g0:g1) = 0.0d0
   do ish = 0,nBases-1
      di = cgto_engine(ish,bas)
      is = ao_offset(ish)
      do jsh = ish,nBases-1
         dj = cgto_engine(jsh,bas)
         js = ao_offset(jsh)
         shls(1) = ish
         shls(2) = jsh
         shls(3) = g0-1
         shls(4) = g1
         call grids1e_engine_cached(buf, shls, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         do q = 1,dj
            e2 = js + q
            do p = 1,di
               e1 = is + p
               if (ish .eq. jsh .and. e2 .lt. e1) cycle
               nij = NorVEC(e1)*NorVEC(e2)
               if (nij .eq. 0.0d0) cycle
               w = Ptot(e1,e2)*nij
               if (e2 .gt. e1) w = 2.0d0*w
               do g = 1,ng
                  V_elec(g0+g-1) = V_elec(g0+g-1) - w*buf(g + ng*((p-1) + di*(q-1)))
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf, int_cache)
!$omp end parallel

deallocate(env_g, ao_offset)
end subroutine esp_at_grid_batch

module subroutine cosmo_fock_batch(npts, grid_bohr, nq, Q, nao, F)
implicit none
integer,intent(in) :: npts, nq, nao
real(8),intent(in) :: grid_bohr(3,npts), Q(npts,nq)
real(8),intent(inout) :: F(nao,nao,*)

integer,parameter :: CHUNK = 256
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),allocatable :: env_g(:), buf(:), int_cache(:), blk(:,:)
integer,allocatable :: ao_offset(:)
integer :: nenv, off_grids, ip, ish, jsh, di, dj, is, js, nmax, iq
integer :: g0, g1, ng, p, qq, e1, e2, shls(4)
real(8) :: nij

nenv = size(env)
off_grids = nenv
allocate(env_g(nenv + 3*npts))
env_g = 0.0d0
env_g(1:nenv) = env(1:nenv)
do ip = 1,npts
   env_g(off_grids + 3*(ip-1) + 1) = grid_bohr(1,ip)
   env_g(off_grids + 3*(ip-1) + 2) = grid_bohr(2,ip)
   env_g(off_grids + 3*(ip-1) + 3) = grid_bohr(3,ip)
enddo
env_g(12) = dble(npts)
env_g(13) = dble(off_grids)

allocate(ao_offset(0:nBases-1))
ip = 0
nmax = 0
do ish = 0,nBases-1
   ao_offset(ish) = ip
   ip = ip + cgto_engine(ish,bas)
   nmax = max(nmax, cgto_engine(ish,bas))
enddo

!$omp parallel private(buf,int_cache,blk,g0,g1,ng,jsh,di,dj,is,js,shls,p,qq,e1,e2,nij,iq)
allocate(buf(nmax*nmax*CHUNK), int_cache(INT_CACHE_SIZE), blk(nmax*nmax,nq))
!$omp do schedule(dynamic)
do ish = 0,nBases-1
   di = cgto_engine(ish,bas)
   is = ao_offset(ish)
   do jsh = ish,nBases-1
      dj = cgto_engine(jsh,bas)
      js = ao_offset(jsh)
      blk(1:di*dj,1:nq) = 0.0d0
      do g0 = 1,npts,CHUNK
         g1 = min(g0+CHUNK-1, npts)
         ng = g1 - g0 + 1
         shls(1) = ish
         shls(2) = jsh
         shls(3) = g0-1
         shls(4) = g1
         call grids1e_engine_cached(buf, shls, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         call dgemm('T','N', di*dj, nq, ng, 1.0d0, buf, ng, Q(g0,1), npts, 1.0d0, blk, nmax*nmax)
      enddo
      do qq = 1,dj
         e2 = js + qq
         do p = 1,di
            e1 = is + p
            nij = NorVEC(e1)*NorVEC(e2)
            do iq = 1,nq
               F(e1,e2,iq) = F(e1,e2,iq) - nij*blk(p + di*(qq-1), iq)
               if (jsh .ne. ish) F(e2,e1,iq) = F(e2,e1,iq) - nij*blk(p + di*(qq-1), iq)
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf, int_cache, blk)
!$omp end parallel

deallocate(env_g, ao_offset)
end subroutine cosmo_fock_batch

module subroutine cosmo_force_elec_batch(npts, grid_bohr, qk, Ptot, ao_atom, natoms, F_shell, field)
implicit none
integer,intent(in) :: npts, natoms
real(8),intent(in) :: grid_bohr(3,npts), qk(npts)
real(8),intent(in) :: Ptot(:,:)
integer,intent(in) :: ao_atom(:)
real(8),intent(out) :: F_shell(3,natoms), field(3,npts)

integer,parameter :: CHUNK = 64
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),allocatable :: env_g(:), buf(:), int_cache(:), Fch(:,:,:)
integer,allocatable :: ao_offset(:)
integer :: nenv, off_grids, ip, ish, jsh, di, dj, is, js, nmax, nch, ich
integer :: g0, g1, ng, g, p, q, e1, e2, x, ia, shls(4), base
real(8) :: w, sx

nenv = size(env)
off_grids = nenv
allocate(env_g(nenv + 3*npts))
env_g = 0.0d0
env_g(1:nenv) = env(1:nenv)
do ip = 1,npts
   env_g(off_grids + 3*(ip-1) + 1) = grid_bohr(1,ip)
   env_g(off_grids + 3*(ip-1) + 2) = grid_bohr(2,ip)
   env_g(off_grids + 3*(ip-1) + 3) = grid_bohr(3,ip)
enddo
env_g(12) = dble(npts)
env_g(13) = dble(off_grids)

allocate(ao_offset(0:nBases-1))
ip = 0
nmax = 0
do ish = 0,nBases-1
   ao_offset(ish) = ip
   ip = ip + cgto_engine(ish,bas)
   nmax = max(nmax, cgto_engine(ish,bas))
enddo

nch = (npts + CHUNK - 1) / CHUNK
allocate(Fch(3,natoms,nch))

!$omp parallel private(buf,int_cache,ich,g0,g1,ng,ish,jsh,di,dj,is,js,shls,g,p,q,e1,e2,x,ia,w,sx,base)
allocate(buf(nmax*nmax*CHUNK*3), int_cache(INT_CACHE_SIZE))
!$omp do schedule(dynamic)
do ich = 1,nch
   g0 = (ich-1)*CHUNK + 1
   g1 = min(g0+CHUNK-1, npts)
   ng = g1 - g0 + 1
   field(:,g0:g1) = 0.0d0
   Fch(:,:,ich) = 0.0d0
   do ish = 0,nBases-1
      di = cgto_engine(ish,bas)
      is = ao_offset(ish)
      do jsh = 0,nBases-1
         dj = cgto_engine(jsh,bas)
         js = ao_offset(jsh)
         shls(1) = ish
         shls(2) = jsh
         shls(3) = g0-1
         shls(4) = g1
         call grids1e_ip_engine_cached(buf, shls, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         do q = 1,dj
            e2 = js + q
            do p = 1,di
               e1 = is + p
               w = Ptot(e1,e2)*NorVEC(e1)*NorVEC(e2)
               if (w .eq. 0.0d0) cycle
               ia = ao_atom(e1)
               do x = 1,3
                  base = ng*((p-1) + di*((q-1) + dj*(x-1)))
                  sx = 0.0d0
                  do g = 1,ng
                     field(x,g0+g-1) = field(x,g0+g-1) + w*buf(base+g)
                     sx = sx + qk(g0+g-1)*buf(base+g)
                  enddo
                  Fch(x,ia,ich) = Fch(x,ia,ich) + w*sx
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf, int_cache)
!$omp end parallel

F_shell = 0.0d0
do ich = 1,nch
   F_shell = F_shell + Fch(:,:,ich)
enddo
deallocate(env_g, ao_offset, Fch)
end subroutine cosmo_force_elec_batch

module subroutine cosmo_hess_elec_batch(npts, grid_bohr, qk, Ptot, ao_atom, natoms, ncol, C, &
                                        gsh, hrow, Gd, W)
use omp_lib, only: omp_get_max_threads, omp_get_thread_num
implicit none
integer,intent(in) :: npts, natoms, ncol
real(8),intent(in) :: grid_bohr(3,npts), qk(npts), Ptot(:,:), C(npts,3,ncol)
integer,intent(in) :: ao_atom(:)
real(8),intent(out) :: gsh(3,natoms,npts), hrow(3,3,natoms,npts), Gd(3,natoms,3,natoms)
real(8),intent(out) :: W(:,:,:)

integer,parameter :: CHUNK = 64
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),allocatable :: env_g(:), b1(:), baa(:), bab(:), int_cache(:), Wb(:,:)
real(8),allocatable :: gsh_t(:,:,:,:), hrow_t(:,:,:,:,:), Gd_t(:,:,:,:,:)
integer,allocatable :: ao_offset(:), pi_(:), pj_(:)
integer :: nenv, off_grids, ip, ish, jsh, di, dj, is, js, nmax, npair, ipair, nth, it
integer :: g0, g1, ng, g, p, q, e1, e2, a, b, iA, iB, shls(4), dd, i, icol
real(8) :: pw, v, nij

nenv = size(env)
off_grids = nenv
allocate(env_g(nenv + 3*npts))
env_g = 0.0d0
env_g(1:nenv) = env(1:nenv)
do ip = 1,npts
   env_g(off_grids + 3*(ip-1) + 1:off_grids + 3*ip) = grid_bohr(:,ip)
enddo
env_g(12) = dble(npts)
env_g(13) = dble(off_grids)

allocate(ao_offset(0:nBases-1))
ip = 0
nmax = 0
do ish = 0,nBases-1
   ao_offset(ish) = ip
   ip = ip + cgto_engine(ish,bas)
   nmax = max(nmax, cgto_engine(ish,bas))
enddo
npair = nBases*nBases
allocate(pi_(npair), pj_(npair))
ipair = 0
do ish = 0,nBases-1
   do jsh = 0,nBases-1
      ipair = ipair+1
      pi_(ipair) = ish ; pj_(ipair) = jsh
   enddo
enddo

nth = omp_get_max_threads()
allocate(gsh_t(3,natoms,npts,nth), hrow_t(3,3,natoms,npts,nth), Gd_t(3,natoms,3,natoms,nth))
gsh_t = 0.0d0 ; hrow_t = 0.0d0 ; Gd_t = 0.0d0
W = 0.0d0

!$omp parallel private(b1,baa,bab,int_cache,Wb,it,ipair,ish,jsh,di,dj,is,js,dd,g0,g1,ng,shls, &
!$omp&   g,p,q,e1,e2,a,b,iA,iB,pw,v,nij,i,icol)
allocate(b1(nmax*nmax*CHUNK*3), baa(nmax*nmax*CHUNK*9), bab(nmax*nmax*CHUNK*9), &
         int_cache(INT_CACHE_SIZE), Wb(nmax*nmax,ncol))
it = omp_get_thread_num() + 1
!$omp do schedule(static)
do ipair = 1,npair
   ish = pi_(ipair) ; jsh = pj_(ipair)
   di = cgto_engine(ish,bas) ; is = ao_offset(ish)
   dj = cgto_engine(jsh,bas) ; js = ao_offset(jsh)
   dd = di*dj
   Wb(1:dd,:) = 0.0d0
   do g0 = 1,npts,CHUNK
      g1 = min(g0+CHUNK-1, npts)
      ng = g1 - g0 + 1
      shls(1) = ish ; shls(2) = jsh ; shls(3) = g0-1 ; shls(4) = g1
      call grids1e_ip_engine_cached(b1, shls, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
      call grids1e_ipip_engine_cached(baa, shls, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
      call grids1e_ipvip_engine_cached(bab, shls, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
      do i = 1,3
         call dgemm('T','N', dd, ncol, ng, 1.0d0, b1(1+ng*dd*(i-1)), ng, C(g0,i,1), npts*3, &
                    1.0d0, Wb, nmax*nmax)
      enddo
      do q = 1,dj
         e2 = js + q
         iB = ao_atom(e2)
         do p = 1,di
            e1 = is + p
            nij = NorVEC(e1)*NorVEC(e2)
            pw = Ptot(e1,e2)*nij
            if (pw .eq. 0.0d0) cycle
            iA = ao_atom(e1)
            do g = 1,ng
               do a = 1,3
                  gsh_t(a,iA,g0+g-1,it) = gsh_t(a,iA,g0+g-1,it) &
                       + pw*b1(g + ng*((p-1) + di*((q-1) + dj*(a-1))))
               enddo
               do b = 1,3
                  do a = 1,3
                     v = -2.0d0*pw*baa(g + ng*((p-1) + di*((q-1) + dj*(3*(a-1)+b-1))))
                     hrow_t(a,b,iA,g0+g-1,it) = hrow_t(a,b,iA,g0+g-1,it) + v
                     Gd_t(a,iA,b,iA,it) = Gd_t(a,iA,b,iA,it) + qk(g0+g-1)*v
                     v = -pw*bab(g + ng*((p-1) + di*((q-1) + dj*(3*(a-1)+b-1))))
                     hrow_t(a,b,iA,g0+g-1,it) = hrow_t(a,b,iA,g0+g-1,it) + v
                     hrow_t(b,a,iB,g0+g-1,it) = hrow_t(b,a,iB,g0+g-1,it) + v
                     Gd_t(a,iA,b,iB,it) = Gd_t(a,iA,b,iB,it) + qk(g0+g-1)*v
                     Gd_t(b,iB,a,iA,it) = Gd_t(b,iB,a,iA,it) + qk(g0+g-1)*v
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   do q = 1,dj
      e2 = js + q
      do p = 1,di
         e1 = is + p
         nij = NorVEC(e1)*NorVEC(e2)
         do icol = 1,ncol
            W(e1,e2,icol) = Wb(p + di*(q-1), icol)*nij
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(b1, baa, bab, int_cache, Wb)
!$omp end parallel

gsh = 0.0d0 ; hrow = 0.0d0 ; Gd = 0.0d0
do it = 1,nth
   gsh = gsh + gsh_t(:,:,:,it)
   hrow = hrow + hrow_t(:,:,:,:,it)
   Gd = Gd + Gd_t(:,:,:,:,it)
enddo
deallocate(env_g, ao_offset, pi_, pj_, gsh_t, hrow_t, Gd_t)
end subroutine cosmo_hess_elec_batch

module subroutine esp_at_grid_batch_df(npts, grid_bohr, Ptot, V_elec, ok)
use MOL_info, only: Natoms
implicit none
integer,intent(in) :: npts
real(8),intent(in) :: grid_bohr(3,npts)
real(8),intent(in) :: Ptot(:,:)
real(8),intent(out) :: V_elec(npts)
logical,intent(out) :: ok

integer,parameter :: CHUNK = 256
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),parameter :: UNIT_S_COEFF = 2.0d0*sqrt(acos(-1.0d0))
real(8),allocatable :: envB(:), buf(:), int_cache(:)
integer,allocatable :: basB(:,:)
integer :: nenvB, off_grids, off_dummy, ip, ksh, dk, offP, nmax
integer :: g0, g1, ng, g, r, shls(4), nshB, idummy

ok = .false.
if (.not. df_built) return
if (.not. df_direct_cvec_cache_valid) return
if (.not. allocated(df_direct_cvec_cache)) return
if (.not. allocated(df_direct_cvec_cache_Ptot)) return
if (size(df_direct_cvec_cache_Ptot,1) .ne. size(Ptot,1)) return
if (any(df_direct_cvec_cache_Ptot .ne. Ptot)) return
ok = .true.

nshB = nBases + nBasesAux + 1
idummy = nBases + nBasesAux
nenvB = size(envDF)
off_dummy = nenvB
off_grids = nenvB + 2
allocate(basB(8,nshB), envB(nenvB + 2 + 3*npts))
basB = 0
basB(1:8,1:nBases+nBasesAux) = basDF(1:8,1:nBases+nBasesAux)
envB = 0.0d0
envB(1:nenvB) = envDF(1:nenvB)

basB(1,nshB) = 0
basB(2,nshB) = 0
basB(3,nshB) = 1
basB(4,nshB) = 1
basB(6,nshB) = off_dummy
basB(7,nshB) = off_dummy + 1
envB(off_dummy+1) = 0.0d0
envB(off_dummy+2) = UNIT_S_COEFF

do ip = 1,npts
   envB(off_grids + 3*(ip-1) + 1) = grid_bohr(1,ip)
   envB(off_grids + 3*(ip-1) + 2) = grid_bohr(2,ip)
   envB(off_grids + 3*(ip-1) + 3) = grid_bohr(3,ip)
enddo
envB(12) = dble(npts)
envB(13) = dble(off_grids)

nmax = 0
do ksh = nBases,nBases+nBasesAux-1
   nmax = max(nmax, df_shell_ncgto(ksh))
enddo

!$omp parallel private(buf,int_cache,g0,g1,ng,ksh,dk,offP,shls,g,r)
allocate(buf(nmax*CHUNK), int_cache(INT_CACHE_SIZE))
!$omp do schedule(dynamic)
do g0 = 1,npts,CHUNK
   g1 = min(g0+CHUNK-1, npts)
   ng = g1 - g0 + 1
   V_elec(g0:g1) = 0.0d0
   do ksh = nBases,nBases+nBasesAux-1
      dk = df_shell_ncgto(ksh)
      offP = df_shell_offset(ksh)
      shls(1) = ksh
      shls(2) = idummy
      shls(3) = g0-1
      shls(4) = g1
      call grids1e_engine_cached(buf, shls, atm, size(atm,2), basB, nshB, envB, 0_8, int_cache)
      do r = 1,dk
         do g = 1,ng
            V_elec(g0+g-1) = V_elec(g0+g-1) &
               - df_direct_cvec_cache(offP+r)*NorVECAux(offP+r)*buf(g + ng*(r-1))
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf, int_cache)
!$omp end parallel

deallocate(basB, envB)
end subroutine esp_at_grid_batch_df

module subroutine cosmo_build_one_tess_matrix_shellderiv(tess_pos_bohr, dBk)
use MOL_info, only: Natoms
implicit none
real(8),intent(in) :: tess_pos_bohr(3)
real(8),intent(out) :: dBk(:,:,:)

integer,allocatable :: atm_k(:,:)
real(8),allocatable :: env_k(:)
real(8),allocatable :: buf1eV(:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: nAtoms_k, offpoint_k
integer :: i,j,di,dj,shls(2)
integer :: si_off, num_off

nAtoms_k = Natoms + 1
allocate(atm_k(6,nAtoms_k))
atm_k = 0
atm_k(:,1:Natoms) = atm(:,1:Natoms)
atm_k(1,1:Natoms) = 0

allocate(env_k(size(env)+4))
env_k = 0.0d0
env_k(1:size(env)) = env(1:size(env))
offpoint_k = size(env)
atm_k(2,nAtoms_k) = offpoint_k
env_k(offpoint_k+1) = tess_pos_bohr(1)
env_k(offpoint_k+2) = tess_pos_bohr(2)
env_k(offpoint_k+3) = tess_pos_bohr(3)
offpoint_k = offpoint_k + 3
atm_k(3,nAtoms_k) = 3
atm_k(5,nAtoms_k) = offpoint_k
env_k(offpoint_k+1) = 1.0d0
offpoint_k = offpoint_k + 1

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo

dBk = 0.0d0
do i = 0,nBases-1
   do j = 0,nBases-1
      shls(1) = i
      shls(2) = j
      di = shell_dim(i)
      dj = shell_dim(j)
      allocate(buf1eV(di,dj,3))
      call nuc1e_ip_engine(buf1eV, shls, atm_k, nAtoms_k, bas, nBases, env_k, 0_8)
      call store1edrv(shls, di, dj, ao_offset, nBases, buf1eV, NorVEC, size(dBk,1), dBk)
      deallocate(buf1eV)
   enddo
enddo

deallocate(atm_k, env_k, ao_offset, shell_dim)
end subroutine cosmo_build_one_tess_matrix_shellderiv

end submodule cosmo_integrals_impl
