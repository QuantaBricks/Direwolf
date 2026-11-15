! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_tddft: singlet excitation energies of a closed-shell RKS reference

module mod_tddft
use mod_integrals, only: cosx_set_resp_grid
implicit none
private
public :: tddft_rks, tddft_uks, tddft_main_transitions, tddft_polar_hybrid, tddft_stability
logical,save,public :: tddft_tda = .false.
character(len=240),save,public :: tddft_note = ''
logical,save :: td_unstable = .false.
real(8),save :: td_unstable_val = 0.0d0
real(8),save :: td_neg_min = 0.0d0, td_pos_nth = 0.0d0
logical,save :: td_uks = .false.
integer,save :: td_nob = 0
real(8),allocatable,save :: td_Cb(:,:), td_eb(:)

contains

subroutine tddft_rks(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                     w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
implicit none
integer,intent(in) :: nConts, nocc, nroots
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: w(nroots), osc(nroots), tdip(3,nroots)
real(8),intent(out) :: XpY(td_ntot(nConts,nocc),nroots), XmY(td_ntot(nConts,nocc),nroots)
integer,intent(out) :: napply
real(8),intent(in) :: conv_thr
logical,intent(in) :: robust
td_unstable = .false.
if (tddft_tda) then
   call tddft_tda_core(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                       w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
   tddft_note = 'Tamm-Dancoff approximation (tddft_tda)'
   return
endif
call tddft_full(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
if (td_unstable) then
   print '(A)', 'WARNING: TDDFT: the full problem has an imaginary excitation energy (w^2 < 0), the reference is'
   print '(A,ES10.3,A)', '         unstable toward a lower solution (lowest w^2 found: ', td_unstable_val, ' a.u.).'
   print '(A)', '         The states below are Tamm-Dancoff (real by construction); a stability analysis or an'
   print '(A)', '         unstable-to-unrestricted re-optimisation of the ground state is advisable.'
   call tddft_tda_core(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                       w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
   tddft_note = 'Tamm-Dancoff approximation: the full problem has an imaginary excitation energy (unstable reference)'
else
   tddft_note = 'full TDDFT'
endif
end subroutine tddft_rks

subroutine tddft_stability(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, triplet, lam, nit)
use mod_cpscf, only: cpscf_apb_rhf, cpscf_apb_amb_rhf, cpscf_response_begin, cpscf_response_end, cpscf_triplet
use mod_exchange, only: RS_omega
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
logical,intent(in) :: triplet
real(8),intent(out) :: lam
integer,intent(out) :: nit
integer :: nvir, nt, nb0, maxsub, m, it, k, j, q, info, lw, nnew, imin, jj
real(8),allocatable :: dl(:), V(:,:), AV(:,:), Hm(:,:), ev(:), wk(:), x(:), r(:), T(:,:), Xr(:,:,:), Pr(:,:,:), Mr(:,:,:)
logical,allocatable :: used(:)
logical :: exx
real(8) :: rn, den, nrm, lamold
td_uks = .false.
nvir = nConts - nocc
nt = nvir*nocc
exx = (HF_exchange_frac /= 0.0d0 .or. RS_omega > 0.0d0)
allocate(dl(nt))
call td_diag(nConts, nocc, eLev_a, dl)
nb0 = min(nt, 3)
maxsub = min(nt, 30)
allocate(V(nt,maxsub), AV(nt,maxsub), x(nt), r(nt), T(nt,nb0), used(nt))
m = nb0
V(:,1:m) = 0.0d0
used = .false.
do k = 1,m
   imin = minloc(dl, 1, mask=.not. used)
   used(imin) = .true.
   V(imin,k) = 1.0d0
enddo
deallocate(used)
cpscf_triplet = triplet
call cpscf_response_begin(nConts, Pa, Pb, HF_exchange_frac, scf_grid=.true.)
call apply_apb(1, m)
lam = 0.0d0 ; lamold = 1.0d30 ; nit = 0
do it = 1,20
   nit = it
   allocate(Hm(m,m), ev(m), wk(1))
   call dgemm('T','N',m,m,nt,1.0d0,V,nt,AV,nt,0.0d0,Hm,m)
   Hm = 0.5d0*(Hm + transpose(Hm))
   call dsyev('V','U',m,Hm,m,ev,wk,-1,info)
   lw = int(wk(1)) ; deallocate(wk) ; allocate(wk(lw))
   call dsyev('V','U',m,Hm,m,ev,wk,lw,info)
   if (info /= 0) error stop 'TDDFT stability: subspace diagonalization failed'
   lam = ev(1)
   nnew = 0
   do q = 1,min(2, m)
      call dgemv('N',nt,m,1.0d0,V,nt,Hm(1,q),1,0.0d0,x,1)
      call dgemv('N',nt,m,1.0d0,AV,nt,Hm(1,q),1,0.0d0,r,1)
      r = r - ev(q)*x
      rn = sqrt(sum(r*r))
      if (rn < 2.0d-4 .and. q == 1) exit
      if (rn < 2.0d-4) cycle
      nnew = nnew + 1
      do k = 1,nt
         den = ev(q) - dl(k)
         if (abs(den) < 1.0d-3) den = sign(1.0d-3, den)
         T(k,nnew) = r(k)/den
      enddo
   enddo
   deallocate(Hm, ev, wk)
   if (nnew == 0) exit
   if (m + nnew > maxsub) exit
   jj = m
   do q = 1,nnew
      nrm = sqrt(sum(T(:,q)**2))
      if (nrm == 0.0d0) cycle
      do k = 1,2
         do j = 1,jj
            T(:,q) = T(:,q) - sum(V(:,j)*T(:,q))*V(:,j)
         enddo
      enddo
      nrm = sqrt(sum(T(:,q)**2))
      if (nrm < 1.0d-6) cycle
      jj = jj + 1
      V(:,jj) = T(:,q)/nrm
   enddo
   if (jj == m) exit
   call apply_apb(m+1, jj)
   m = jj
enddo
call cpscf_response_end()
cpscf_triplet = .false.
deallocate(dl, V, AV, x, r, T)

contains

subroutine apply_apb(j1, j2)
integer,intent(in) :: j1, j2
integer :: nb, jx
nb = j2 - j1 + 1
allocate(Xr(nvir,nocc,nb), Pr(nvir,nocc,nb), Mr(nvir,nocc,1))
do jx = 1,nb
   Xr(:,:,jx) = reshape(V(:,j1+jx-1), (/nvir,nocc/))
enddo
if (exx) then
   call cpscf_apb_amb_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, nb, Xr, Pr, 0, Mr, Mr)
else
   call cpscf_apb_rhf(nConts, nocc, nvir, nb, C_a, eLev_a, HF_exchange_frac, Xr, Pr)
endif
do jx = 1,nb
   AV(:,j1+jx-1) = reshape(Pr(:,:,jx), (/nt/))
enddo
deallocate(Xr, Pr, Mr)
end subroutine apply_apb

end subroutine tddft_stability

subroutine tddft_tda_core(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                          w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
use mod_cpscf, only: cpscf_apb_rhf, cpscf_apb_amb_rhf, cpscf_apb_amb_uhf, cpscf_response_end
use mod_cpscf, only: cphf_t_vj, cphf_t_vk, cphf_t_vxc
use mod_exchange, only: RS_omega
use MOL_info, only: engine_verbose
use omp_lib, only: omp_get_wtime
use mod_xc_hess, only: xc_fxc_ngc
implicit none
integer,intent(in) :: nConts, nocc, nroots
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: w(nroots), osc(nroots), tdip(3,nroots)
real(8),intent(out) :: XpY(td_ntot(nConts,nocc),nroots), XmY(td_ntot(nConts,nocc),nroots)
integer,intent(out) :: napply
real(8),intent(in) :: conv_thr
logical,intent(in) :: robust
integer :: nvir, nt, nwork, nkeep, ngs, maxsub, m, it, k, r, q, qq, j, info, lw, nnew, imin, nreal, nk
real(8),allocatable :: dl(:), V(:,:), AV(:,:), Hm(:,:), Hc(:,:), evr(:), evi(:), vl(:,:), vr(:,:), wk(:), &
                       Rv(:), Ax(:), Tn(:,:), wr(:), Cs(:,:), S(:,:), Wt(:,:)
real(8),allocatable :: scl(:), Cq(:,:), Rb(:,:), Ab(:,:)
integer,allocatable :: ord(:)
logical,allocatable :: conv(:), used(:)
logical :: exx
real(8) :: rn, den, nrm, nrm0, thr, rmax, t0, tapp, cn

nvir = nConts - nocc
nt = td_ntot(nConts, nocc)
if (nroots < 1 .or. nroots > nt) error stop 'TDDFT: tddft_nroots must be between 1 and nocc*nvir'
exx = (HF_exchange_frac /= 0.0d0 .or. RS_omega > 0.0d0)
allocate(dl(nt))
call td_diag(nConts, nocc, eLev_a, dl)
if (minval(dl) <= 0.0d0) error stop 'TDDFT: a virtual orbital lies below an occupied one'
nwork = nroots
if (robust) nwork = min(nt, nroots + max(nroots, 6))
nkeep = min(nt, 2*nwork)
ngs = min(nt, max(3*nroots, 2*nwork))
maxsub = min(nt, max(20*nwork, ngs + 2*nwork, 60))
allocate(V(nt,maxsub), AV(nt,maxsub), conv(nwork), used(nt), Rv(nt), Ax(nt), Tn(nt,nwork), wr(nwork))
allocate(scl(maxsub))
scl = 1.0d0
m = min(ngs, maxsub)
V(:,1:m) = 0.0d0
used = .false.
do k = 1,m
   imin = minloc(dl, 1, mask=.not. used)
   used(imin) = .true.
   V(imin,k) = 1.0d0
enddo
deallocate(used)
call td_begin(nConts, Pa, Pb, HF_exchange_frac)
xc_fxc_ngc = 512
napply = 0
t0 = omp_get_wtime()
tapp = 0.0d0
call apply_a(1, m)
conv = .false.
do it = 1,200
   allocate(Hm(m,m), Hc(m,m), evr(m), evi(m), vl(1,1), vr(m,m), wk(1), ord(m))
   call dgemm('T','N',m,m,nt,1.0d0,V,nt,AV,nt,0.0d0,Hm,m)
   Hc = Hm
   call dgeev('N','V',m,Hc,m,evr,evi,vl,1,vr,m,wk,-1,info)
   lw = int(wk(1)) ; deallocate(wk) ; allocate(wk(lw))
   call dgeev('N','V',m,Hc,m,evr,evi,vl,1,vr,m,wk,lw,info)
   if (info /= 0) error stop 'TDDFT (TDA): subspace diagonalization failed'
   nreal = 0
   do q = 1,m
      if (abs(evi(q)) <= 1.0d-8*max(1.0d0, abs(evr(q)))) then
         nreal = nreal + 1 ; ord(nreal) = q
      endif
   enddo
   do q = 2,nreal
      qq = ord(q) ; j = q - 1
      do while (j >= 1)
         if (evr(ord(j)) <= evr(qq)) exit
         ord(j+1) = ord(j) ; j = j - 1
      enddo
      ord(j+1) = qq
   enddo
   if (nreal < min(nwork, m)) error stop 'TDDFT (TDA): too few real roots in the subspace'
   nnew = 0 ; rmax = 0.0d0
   block
      integer :: nr_, rr_
      nr_ = min(nwork, m)
      allocate(Cq(m,nr_), Rb(nt,nr_), Ab(nt,nr_))
      do rr_ = 1,nr_
         Cq(:,rr_) = vr(1:m,ord(rr_))
      enddo
      call dgemm('N','N',nt,nr_,m,1.0d0,V,nt,Cq,m,0.0d0,Rb,nt)
      call dgemm('N','N',nt,nr_,m,1.0d0,AV,nt,Cq,m,0.0d0,Ab,nt)
   end block
   do r = 1,min(nwork, m)
      qq = ord(r)
      wr(r) = evr(qq)
      Rv = Rb(:,r) ; Ax = Ab(:,r)
      nrm = sqrt(sum(Rv*Rv))
      Rv = Rv/nrm ; Ax = Ax/nrm
      Ax = Ax - wr(r)*Rv
      rn = sqrt(sum(Ax*Ax))
      thr = conv_thr/(2.0d0*max(abs(wr(r)), 0.05d0))
      conv(r) = rn < thr
      if (r <= nroots) rmax = max(rmax, rn/thr*conv_thr)
      if (conv(r)) cycle
      nnew = nnew + 1
      do k = 1,nt
         den = wr(r) - dl(k)
         if (abs(den) < 1.0d-4) den = sign(1.0d-4, den)
         Tn(k,nnew) = Ax(k)/den
      enddo
   enddo
   deallocate(Cq, Rb, Ab)
   if (engine_verbose >= 2) print '(A,I4,A,I5,A,I4,A,ES9.2,A,F9.2,A)', '[tddft] it', it, '  space', m, &
      '  conv', count(conv(1:nroots)), '  max res', rmax, '  t', omp_get_wtime() - t0, ' s'
   if (all(conv(1:nroots))) then
      do r = 1,nroots
         qq = ord(r)
         call dgemv('N',nt,m,1.0d0,V,nt,vr(1,qq),1,0.0d0,Rv,1)
         Rv = Rv/sqrt(sum(Rv*Rv))
         w(r) = evr(qq)
         XpY(:,r) = Rv ; XmY(:,r) = Rv
      enddo
      deallocate(Hm, Hc, evr, evi, vl, vr, wk, ord)
      exit
   endif
   if (m + nnew > maxsub) then
      allocate(Cs(m,nkeep), Wt(nt,nkeep))
      nk = 0
      do r = 1,min(nkeep, nreal)
         qq = ord(r)
         Cs(:,nk+1) = vr(1:m,qq)
         do q = 1,2
            do j = 1,nk
               Cs(:,nk+1) = Cs(:,nk+1) - sum(Cs(:,j)*Cs(:,nk+1))*Cs(:,j)
            enddo
         enddo
         cn = sqrt(sum(Cs(:,nk+1)**2))
         if (cn < 1.0d-8) cycle
         nk = nk + 1
         Cs(:,nk) = Cs(:,nk)/cn
      enddo
      call dgemm('N','N',nt,nk,m,1.0d0,V,nt,Cs,m,0.0d0,Wt,nt)
      V(:,1:nk) = Wt(:,1:nk)
      call dgemm('N','N',nt,nk,m,1.0d0,AV,nt,Cs,m,0.0d0,Wt,nt)
      AV(:,1:nk) = Wt(:,1:nk)
      m = nk
      deallocate(Cs, Wt)
   endif
   deallocate(Hm, Hc, evr, evi, vl, vr, wk, ord)
   k = m
   do r = 1,nnew
      nrm0 = sqrt(sum(Tn(:,r)**2))
      if (nrm0 == 0.0d0) cycle
      do q = 1,2
         do j = 1,k
            Tn(:,r) = Tn(:,r) - sum(V(:,j)*Tn(:,r))*V(:,j)
         enddo
      enddo
      nrm = sqrt(sum(Tn(:,r)**2))
      if (nrm < 1.0d-6*nrm0 .or. k >= maxsub) cycle
      k = k + 1
      V(:,k) = Tn(:,r)/nrm
      scl(k) = max(1.0d-4, min(1.0d0, nrm))
   enddo
   if (k == m) then
      print '(A)', 'WARNING: TDDFT (TDA) Davidson made no new directions; stopping with the current roots'
      exit
   endif
   call apply_a(m+1, k)
   m = k
enddo
if (.not. all(conv(1:nroots))) then
   print '(A)', 'WARNING: TDDFT (TDA) Davidson did not converge every root'
   if (allocated(Hm)) then
      do r = 1,nroots
         qq = ord(r)
         call dgemv('N',nt,m,1.0d0,V,nt,vr(1,qq),1,0.0d0,Rv,1)
         Rv = Rv/sqrt(sum(Rv*Rv))
         w(r) = evr(qq)
         XpY(:,r) = Rv ; XmY(:,r) = Rv
      enddo
      deallocate(Hm, Hc, evr, evi, vl, vr, wk, ord)
   endif
endif
if (engine_verbose >= 2) print '(A,F9.2,A,F9.2,A,F9.2,A,F9.2,A,F9.2)', '[tddft] total', omp_get_wtime() - t0, &
   '  apply', tapp, '  J', cphf_t_vj, '  K', cphf_t_vk, '  XC', cphf_t_vxc
call td_end()
xc_fxc_ngc = 128
call td_props(nConts, nocc, C_a, nroots, w, XpY, XmY, tdip, osc)
deallocate(dl, V, AV, conv, Rv, Ax, Tn, wr)

contains

subroutine apply_a(j1, j2)
integer,intent(in) :: j1, j2
integer :: nb, jj
real(8) :: ta
real(8),allocatable :: Xf(:,:), Qp(:,:), Qm(:,:), Xr(:,:,:), Pr(:,:,:), Mr(:,:,:)
ta = omp_get_wtime()
nb = j2 - j1 + 1
allocate(Xf(nt,nb), Qp(nt,nb), Qm(nt,nb))
do jj = 1,nb
   Xf(:,jj) = scl(j1+jj-1)*V(:,j1+jj-1)
enddo
if (td_uks) then
   call cpscf_apb_amb_uhf(nConts, nocc, td_nob, C_a, eLev_a, td_Cb, td_eb, HF_exchange_frac, nb, Xf, Qp, nb, Xf, Qm)
   AV(:,j1:j2) = 0.5d0*(Qp + Qm)
else
   allocate(Xr(nvir,nocc,nb), Pr(nvir,nocc,nb), Mr(nvir,nocc,nb))
   do jj = 1,nb
      Xr(:,:,jj) = reshape(Xf(:,jj), (/nvir,nocc/))
   enddo
   if (exx) then
      call cpscf_apb_amb_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, nb, Xr, Pr, nb, Xr, Mr)
      do jj = 1,nb
         AV(:,j1+jj-1) = 0.5d0*(reshape(Pr(:,:,jj), (/nt/)) + reshape(Mr(:,:,jj), (/nt/)))
      enddo
   else
      call cpscf_apb_rhf(nConts, nocc, nvir, nb, C_a, eLev_a, HF_exchange_frac, Xr, Pr)
      do jj = 1,nb
         AV(:,j1+jj-1) = 0.5d0*(reshape(Pr(:,:,jj), (/nt/)) + dl*Xf(:,jj))
      enddo
   endif
   deallocate(Xr, Pr, Mr)
endif
do jj = 1,nb
   AV(:,j1+jj-1) = AV(:,j1+jj-1)/scl(j1+jj-1)
enddo
napply = napply + nb
tapp = tapp + (omp_get_wtime() - ta)
deallocate(Xf, Qp, Qm)
end subroutine apply_a

end subroutine tddft_tda_core

subroutine tddft_full(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                      w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
use mod_cpscf, only: cpscf_apb_rhf, cpscf_apb_amb_uhf, cpscf_response_end
use mod_exchange, only: RS_omega
use mod_cpscf, only: cphf_t_vj, cphf_t_vk, cphf_t_vxc
use MOL_info, only: engine_verbose
use omp_lib, only: omp_get_wtime
use mod_xc_hess, only: xc_fxc_ngc
implicit none
integer,intent(in) :: nConts, nocc, nroots
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: w(nroots), osc(nroots), tdip(3,nroots)
real(8),intent(out) :: XpY(td_ntot(nConts,nocc),nroots), XmY(td_ntot(nConts,nocc),nroots)
integer,intent(out) :: napply
real(8),intent(in) :: conv_thr
logical,intent(in) :: robust
integer :: nvir, nt, maxsub, m, nnew, it, k, r, i, a, info, lwork, nkeep, x
real(8),allocatable :: sd(:), d2(:), V(:,:), W2(:,:), G(:,:), ev(:), work(:), S(:,:)
real(8),allocatable :: T(:,:), Rv(:), Xb(:,:,:), Qb(:,:,:), Rb(:,:)
logical,allocatable :: conv(:), used(:)
real(8) :: rn, den, nrm, nrm0
integer :: imin, nwork, ngs
real(8) :: t0, tapp, rmax

if (HF_exchange_frac /= 0.0d0 .or. RS_omega > 0.0d0) then
   block
      character(8) :: hm
      call get_environment_variable('ENGINE_TDDFT_SOLVER', hm)
      if (trim(hm) == 'single' .and. .not. td_uks) then
         call tddft_hybrid1(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                            w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
      else
         call tddft_hybrid(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                           w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
      endif
   end block
   return
endif
nvir = nConts - nocc
nt = td_ntot(nConts, nocc)
if (nroots < 1 .or. nroots > nt) error stop 'TDDFT: tddft_nroots must be between 1 and nocc*nvir'

allocate(sd(nt), d2(nt))
call td_diag(nConts, nocc, eLev_a, d2)
if (minval(d2) <= 0.0d0) error stop 'TDDFT: a virtual orbital lies below an occupied one'
sd = sqrt(d2)
d2 = d2*d2

nwork = min(nt, nroots + 2)
if (robust) nwork = min(nt, nroots + max(nroots, 6))
ngs = max(3*nroots, 2*nwork)
maxsub = min(nt, max(20*nwork, ngs + 2*nwork, 60))
allocate(V(nt,maxsub), W2(nt,maxsub), conv(nwork), used(nt), T(nt,nwork), Rv(nt))

m = min(nt, ngs, maxsub)
m = min(m, maxsub)
V(:,1:m) = 0.0d0
used = .false.
do k = 1,m
   imin = minloc(d2, 1, mask=.not. used)
   used(imin) = .true.
   V(imin,k) = 1.0d0
enddo
deallocate(used)

call td_begin(nConts, Pa, Pb, HF_exchange_frac)
xc_fxc_ngc = 512
napply = 0
t0 = omp_get_wtime()
tapp = 0.0d0
call apply_omega(1, m)
conv = .false.
do it = 1,200
   allocate(G(m,m), ev(m), work(1))
   call dgemm('T','N',m,m,nt,1.0d0,V,nt,W2,nt,0.0d0,G,m)
   G = 0.5d0*(G + transpose(G))
   call dsyev('V','U',m,G,m,ev,work,-1,info)
   lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
   call dsyev('V','U',m,G,m,ev,work,lwork,info)
   if (info /= 0) error stop 'TDDFT: subspace diagonalization failed'
   nnew = 0
   rmax = 0.0d0
   block
      integer :: nr_, rr_
      real(8),allocatable :: Gs_(:,:)
      nr_ = min(nwork, m)
      allocate(Rb(nt,nr_), Gs_(m,nr_))
      do rr_ = 1,nr_
         Gs_(:,rr_) = -ev(rr_)*G(1:m,rr_)
      enddo
      call dgemm('N','N',nt,nr_,m,1.0d0,W2,nt,G,m,0.0d0,Rb,nt)
      call dgemm('N','N',nt,nr_,m,1.0d0,V,nt,Gs_,m,1.0d0,Rb,nt)
   end block
   do r = 1,min(nwork, m)
      Rv = Rb(:,r)
      rn = sqrt(sum(Rv*Rv))
      conv(r) = rn < conv_thr
      if (r <= nroots) rmax = max(rmax, rn)
      if (conv(r)) cycle
      nnew = nnew + 1
      do k = 1,nt
         den = ev(r) - d2(k)
         if (abs(den) < 1.0d-4) den = sign(1.0d-4, den)
         T(k,nnew) = Rv(k)/den
      enddo
   enddo
   deallocate(Rb)
   if (engine_verbose >= 2) print '(A,I4,A,I5,A,I4,A,ES9.2,A,F9.2,A)', '[tddft] it', it, '  space', m, &
      '  conv', count(conv(1:nroots)), '  max res', rmax, '  t', omp_get_wtime() - t0, ' s'
   if (all(conv(1:nroots))) exit
   if (m + nnew > maxsub) then
      nkeep = min(m, 2*nwork)
      allocate(S(nt,nkeep))
      call dgemm('N','N',nt,nkeep,m,1.0d0,V,nt,G,m,0.0d0,S,nt)
      V(:,1:nkeep) = S
      call dgemm('N','N',nt,nkeep,m,1.0d0,W2,nt,G,m,0.0d0,S,nt)
      W2(:,1:nkeep) = S
      deallocate(S)
      m = nkeep
   endif
   deallocate(G, ev, work)
   k = m
   do r = 1,nnew
      nrm0 = sqrt(sum(T(:,r)**2))
      if (nrm0 == 0.0d0) cycle
      do x = 1,2
         call orth_against(k, T(:,r))
      enddo
      nrm = sqrt(sum(T(:,r)**2))
      if (nrm < 1.0d-6*nrm0 .or. k >= maxsub) cycle
      k = k + 1
      V(:,k) = T(:,r)/nrm
   enddo
   if (k == m) then
      print '(A)', 'WARNING: TDDFT Davidson made no new directions; stopping with the current roots'
      exit
   endif
   call apply_omega(m+1, k)
   m = k
enddo
if (.not. all(conv(1:nroots))) print '(A)', 'WARNING: TDDFT Davidson did not converge every root'
if (.not. allocated(G)) then
   allocate(G(m,m), ev(m), work(1))
   call dgemm('T','N',m,m,nt,1.0d0,V,nt,W2,nt,0.0d0,G,m)
   G = 0.5d0*(G + transpose(G))
   call dsyev('V','U',m,G,m,ev,work,-1,info)
   lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
   call dsyev('V','U',m,G,m,ev,work,lwork,info)
endif
if (engine_verbose >= 2) print '(A,F9.2,A,F9.2,A,F9.2,A,F9.2,A,F9.2)', '[tddft] total', omp_get_wtime() - t0, &
   '  apply', tapp, '  J', cphf_t_vj, '  K', cphf_t_vk, '  XC', cphf_t_vxc
call td_end()
xc_fxc_ngc = 128

if (ev(1) < -1.0d-8) then
   td_unstable = .true. ; td_unstable_val = ev(1)
   deallocate(sd, d2, V, W2, G, ev, work, conv, T, Rv)
   return
endif
do r = 1,nroots
   w(r) = sqrt(max(ev(r), 0.0d0))
   call dgemv('N',nt,m,1.0d0,V,nt,G(1,r),1,0.0d0,Rv,1)
   Rv = Rv/sqrt(sum(Rv*Rv))
   XpY(:,r) = sd*Rv/sqrt(w(r))
   XmY(:,r) = sqrt(w(r))*Rv/sd
enddo
call td_props(nConts, nocc, C_a, nroots, w, XpY, XmY, tdip, osc)
deallocate(sd, d2, V, W2, G, ev, work, conv, T, Rv)

contains

subroutine apply_omega(j1, j2)
integer,intent(in) :: j1, j2
integer :: nb, j
real(8) :: ta
ta = omp_get_wtime()
nb = j2 - j1 + 1
if (td_uks) then
   block
      real(8),allocatable :: Xf(:,:), Qf(:,:), Zf(:,:)
      allocate(Xf(nt,nb), Qf(nt,nb), Zf(nt,1))
      do j = 1,nb
         Xf(:,j) = sd*V(:,j1+j-1)
      enddo
      call cpscf_apb_amb_uhf(nConts, nocc, td_nob, C_a, eLev_a, td_Cb, td_eb, HF_exchange_frac, &
                             nb, Xf, Qf, 0, Zf, Zf)
      do j = 1,nb
         W2(:,j1+j-1) = sd*Qf(:,j)
      enddo
      deallocate(Xf, Qf, Zf)
   end block
   napply = napply + nb
   tapp = tapp + (omp_get_wtime() - ta)
   return
endif
allocate(Xb(nvir,nocc,nb), Qb(nvir,nocc,nb))
do j = 1,nb
   Xb(:,:,j) = reshape(sd*V(:,j1+j-1), (/nvir,nocc/))
enddo
call cpscf_apb_rhf(nConts, nocc, nvir, nb, C_a, eLev_a, HF_exchange_frac, Xb, Qb)
do j = 1,nb
   W2(:,j1+j-1) = sd*reshape(Qb(:,:,j), (/nt/))
enddo
napply = napply + nb
tapp = tapp + (omp_get_wtime() - ta)
deallocate(Xb, Qb)
end subroutine apply_omega

subroutine orth_against(n, t)
integer,intent(in) :: n
real(8),intent(inout) :: t(nt)
real(8),allocatable :: c(:)
if (n < 1) return
allocate(c(n))
call dgemv('T',nt,n,1.0d0,V,nt,t,1,0.0d0,c,1)
call dgemv('N',nt,n,-1.0d0,V,nt,c,1,1.0d0,t,1)
deallocate(c)
end subroutine orth_against

end subroutine tddft_full

subroutine tddft_props(nConts, nocc, C_a, nroots, w, XpY, XmY, tdip, osc)
use mod_integrals, only: integrals_dipole_ao
implicit none
integer,intent(in) :: nConts, nocc, nroots
real(8),intent(in) :: C_a(nConts,nConts), w(nroots)
real(8),intent(inout) :: XpY(nConts-nocc,nocc,nroots), XmY(nConts-nocc,nocc,nroots)
real(8),intent(out) :: tdip(3,nroots), osc(nroots)
real(8),allocatable :: Rao(:,:,:), Tm(:,:), Vmo(:,:,:)
integer :: nvir, x, r
integer :: imx(2)
nvir = nConts - nocc
do r = 1,nroots
   imx = maxloc(abs(XpY(:,:,r)))
   if (XpY(imx(1),imx(2),r) < 0.0d0) then
      XpY(:,:,r) = -XpY(:,:,r)
      XmY(:,:,r) = -XmY(:,:,r)
   endif
enddo
allocate(Rao(nConts,nConts,3), Tm(nConts,nocc), Vmo(nvir,nocc,3))
call integrals_dipole_ao(nConts, Rao)
do x = 1,3
   call dgemm('N','N',nConts,nocc,nConts,1.0d0,Rao(1,1,x),nConts,C_a,nConts,0.0d0,Tm,nConts)
   call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,Tm,nConts,0.0d0,Vmo(1,1,x),nvir)
enddo
do r = 1,nroots
   do x = 1,3
      tdip(x,r) = sqrt(2.0d0)*sum(Vmo(:,:,x)*XpY(:,:,r))
   enddo
   osc(r) = 2.0d0/3.0d0*w(r)*sum(tdip(:,r)**2)
enddo
deallocate(Rao, Tm, Vmo)
end subroutine tddft_props

subroutine tddft_hybrid(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                        w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
use mod_cpscf, only: cpscf_apb_amb_rhf, cpscf_apb_amb_uhf, cpscf_response_end
use mod_cpscf, only: cphf_t_vj, cphf_t_vk, cphf_t_vxc
use MOL_info, only: engine_verbose
use omp_lib, only: omp_get_wtime
use mod_xc_hess, only: xc_fxc_ngc
implicit none
integer,intent(in) :: nConts, nocc, nroots
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: w(nroots), osc(nroots), tdip(3,nroots)
real(8),intent(out) :: XpY(td_ntot(nConts,nocc),nroots), XmY(td_ntot(nConts,nocc),nroots)
integer,intent(out) :: napply
real(8),intent(in) :: conv_thr
logical,intent(in) :: robust
integer :: nvir, nt, nwork, nkeep, ngs, maxsub, mp, mm, it, k, r, i, a, info, np, nm, kp, km, imin
real(8),allocatable :: dl(:), Vp(:,:), Vm(:,:), P(:,:), Mv(:,:), Cc(:,:), Dd(:,:), wr(:), R1(:), R2(:), Tp(:,:), Tm(:,:)
real(8),allocatable :: sclp(:), sclm(:), Rb1(:,:), Rb2(:,:), Cw(:,:), Dw(:,:)
logical,allocatable :: conv(:), used(:)
real(8) :: om, den, rmax, t0, tapp, rn

nvir = nConts - nocc
nt = td_ntot(nConts, nocc)
allocate(dl(nt))
call td_diag(nConts, nocc, eLev_a, dl)
if (minval(dl) <= 0.0d0) error stop 'TDDFT: a virtual orbital lies below an occupied one'
nwork = nroots
if (robust) nwork = min(nt, nroots + max(nroots, 6))
nkeep = min(nt, 2*nwork)
ngs = min(nt, max(3*nroots, 2*nwork))
maxsub = min(nt, max(ngs + 2*nwork, nkeep + 2*nwork, 20*nwork, 60))
allocate(sclp(maxsub), sclm(maxsub))
sclp = 1.0d0 ; sclm = 1.0d0
allocate(Vp(nt,maxsub), Vm(nt,maxsub), P(nt,maxsub), Mv(nt,maxsub), conv(nwork), wr(nkeep), &
         Cc(maxsub,nkeep), Dd(maxsub,nkeep), R1(nt), R2(nt), Tp(nt,nwork), Tm(nt,nwork), used(nt))
allocate(Rb1(nt,nwork), Rb2(nt,nwork), Cw(maxsub,nkeep), Dw(maxsub,nkeep))
mp = ngs
Vp(:,1:mp) = 0.0d0
used = .false.
do k = 1,mp
   imin = minloc(dl, 1, mask=.not. used)
   used(imin) = .true.
   Vp(imin,k) = 1.0d0
enddo
deallocate(used)
mm = mp
Vm(:,1:mm) = Vp(:,1:mp)

call td_begin(nConts, Pa, Pb, HF_exchange_frac)
xc_fxc_ngc = 512
napply = 0
t0 = omp_get_wtime()
tapp = 0.0d0
call apply_pm(1, mp, 1, mm)
conv = .false.
do it = 1,200
   call subspace()
   if (td_unstable) exit
   np = 0 ; nm = 0
   rmax = 0.0d0
   do r = 1,nwork
      Cw(1:mp,r) = wr(r)*Cc(1:mp,r)
      Dw(1:mm,r) = wr(r)*Dd(1:mm,r)
   enddo
   call dgemm('N','N',nt,nwork,mp,1.0d0,P,nt,Cc,maxsub,0.0d0,Rb1,nt)
   call dgemm('N','N',nt,nwork,mm,-1.0d0,Vm,nt,Dw,maxsub,1.0d0,Rb1,nt)
   call dgemm('N','N',nt,nwork,mm,1.0d0,Mv,nt,Dd,maxsub,0.0d0,Rb2,nt)
   call dgemm('N','N',nt,nwork,mp,-1.0d0,Vp,nt,Cw,maxsub,1.0d0,Rb2,nt)
   do r = 1,nwork
      om = wr(r)
      R1 = Rb1(:,r) ; R2 = Rb2(:,r)
      rn = max(sqrt(sum(R1*R1)), sqrt(sum(R2*R2)))
      conv(r) = rn < conv_thr
      if (r <= nroots) rmax = max(rmax, rn)
      if (conv(r)) cycle
      np = np + 1 ; nm = nm + 1
      do k = 1,nt
         den = dl(k)*dl(k) - om*om
         if (abs(den) < 1.0d-8) den = sign(1.0d-8, den)
         Tp(k,np) = -(dl(k)*R1(k) + om*R2(k))/den
         Tm(k,nm) = -(dl(k)*R2(k) + om*R1(k))/den
      enddo
   enddo
   if (engine_verbose >= 2) print '(A,I4,A,2I5,A,I4,A,ES9.2,A,F9.2,A)', '[tddft] it', it, '  spaces', mp, mm, &
      '  conv', count(conv(1:nroots)), '  max res', rmax, '  t', omp_get_wtime() - t0, ' s'
   if (all(conv(1:nroots))) exit
   if (mp + np > maxsub .or. mm + nm > maxsub) call collapse()
   kp = mp
   call add_dirs(Vp, kp, np, Tp, sclp)
   km = mm
   call add_dirs(Vm, km, nm, Tm, sclm)
   if (kp == mp .and. km == mm) then
      print '(A)', 'WARNING: TDDFT Davidson made no new directions; stopping with the current roots'
      exit
   endif
   call apply_pm(mp+1, kp, mm+1, km)
   mp = kp ; mm = km
enddo
if (.not. td_unstable .and. td_neg_min < 0.0d0 .and. abs(td_neg_min) > td_pos_nth) then
   td_unstable = .true. ; td_unstable_val = 1.0d0/td_neg_min
endif
if (td_unstable) then
   call td_end()
   xc_fxc_ngc = 128
   deallocate(dl, Vp, Vm, P, Mv, conv, wr, Cc, Dd, R1, R2, Tp, Tm)
   return
endif
if (.not. all(conv(1:nroots))) print '(A)', 'WARNING: TDDFT Davidson did not converge every root'
if (engine_verbose >= 2) print '(A,F9.2,A,F9.2,A,F9.2,A,F9.2,A,F9.2)', '[tddft] total', omp_get_wtime() - t0, &
   '  apply', tapp, '  J', cphf_t_vj, '  K', cphf_t_vk, '  XC', cphf_t_vxc
call td_end()
xc_fxc_ngc = 128
do r = 1,nroots
   w(r) = wr(r)
   call dgemv('N',nt,mp,1.0d0,Vp,nt,Cc(1,r),1,0.0d0,XpY(1,r),1)
   call dgemv('N',nt,mm,1.0d0,Vm,nt,Dd(1,r),1,0.0d0,XmY(1,r),1)
enddo
call td_props(nConts, nocc, C_a, nroots, w, XpY, XmY, tdip, osc)
deallocate(dl, Vp, Vm, P, Mv, conv, wr, Cc, Dd, R1, R2, Tp, Tm)

contains

real(8) function sca(j, which)
integer,intent(in) :: j, which
logical,save :: init = .false., on = .false.
character(len=8) :: e
if (.not. init) then
   call get_environment_variable("ENGINE_TDDFT_SCALE", e)
   on = (trim(e) .eq. "1")
   init = .true.
endif
if (.not. on) then
   sca = 1.0d0
else if (which == 1) then
   sca = sclp(j)
else
   sca = sclm(j)
endif
end function sca

subroutine apply_pm(jp1, jp2, jm1, jm2)
integer,intent(in) :: jp1, jp2, jm1, jm2
real(8),allocatable :: Xa(:,:,:), Xb(:,:,:), Qa(:,:,:), Qb(:,:,:)
integer :: na, nb, jj
real(8) :: ta
ta = omp_get_wtime()
na = max(0, jp2 - jp1 + 1) ; nb = max(0, jm2 - jm1 + 1)
if (td_uks) then
   block
      real(8),allocatable :: Xf(:,:), Qf(:,:), Yf(:,:), Rf(:,:)
      allocate(Xf(nt,max(na,1)), Qf(nt,max(na,1)), Yf(nt,max(nb,1)), Rf(nt,max(nb,1)))
      do jj = 1,na
         Xf(:,jj) = sca(jp1+jj-1,1)*Vp(:,jp1+jj-1)
      enddo
      do jj = 1,nb
         Yf(:,jj) = sca(jm1+jj-1,2)*Vm(:,jm1+jj-1)
      enddo
      call cpscf_apb_amb_uhf(nConts, nocc, td_nob, C_a, eLev_a, td_Cb, td_eb, HF_exchange_frac, &
                             na, Xf, Qf, nb, Yf, Rf)
      do jj = 1,na
         P(:,jp1+jj-1) = Qf(:,jj)/sca(jp1+jj-1,1)
      enddo
      do jj = 1,nb
         Mv(:,jm1+jj-1) = Rf(:,jj)/sca(jm1+jj-1,2)
      enddo
      deallocate(Xf, Qf, Yf, Rf)
   end block
   napply = napply + na + nb
   tapp = tapp + (omp_get_wtime() - ta)
   return
endif
allocate(Xa(nvir,nocc,max(na,1)), Qa(nvir,nocc,max(na,1)), Xb(nvir,nocc,max(nb,1)), Qb(nvir,nocc,max(nb,1)))
do jj = 1,na
   Xa(:,:,jj) = reshape(sca(jp1+jj-1,1)*Vp(:,jp1+jj-1), (/nvir,nocc/))
enddo
do jj = 1,nb
   Xb(:,:,jj) = reshape(sca(jm1+jj-1,2)*Vm(:,jm1+jj-1), (/nvir,nocc/))
enddo
call cpscf_apb_amb_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, na, Xa, Qa, nb, Xb, Qb)
do jj = 1,na
   P(:,jp1+jj-1) = reshape(Qa(:,:,jj), (/nt/))/sca(jp1+jj-1,1)
enddo
do jj = 1,nb
   Mv(:,jm1+jj-1) = reshape(Qb(:,:,jj), (/nt/))/sca(jm1+jj-1,2)
enddo
napply = napply + na + nb
tapp = tapp + (omp_get_wtime() - ta)
deallocate(Xa, Xb, Qa, Qb)
end subroutine apply_pm

subroutine subspace()
real(8),allocatable :: Ap(:,:), Am(:,:), Sv(:,:), X1(:,:), G(:,:), U(:,:), ev(:), wi(:), vl(:,:), wk(:), tv(:), Af(:,:)
integer,allocatable :: ord(:), ipv(:)
integer :: lw, q, qq, j, ng
real(8) :: cn
allocate(Ap(mp,mp), Am(mm,mm), Sv(mp,mm), X1(mp,mm), G(mm,mm), U(mm,mm), ev(mm), wi(mm), vl(1,1), wk(1), &
         tv(mp), ord(mm), ipv(max(mp,mm)), Af(mp,mp))
call dgemm('T','N',mp,mp,nt,1.0d0,Vp,nt,P,nt,0.0d0,Ap,mp)
call dgemm('T','N',mm,mm,nt,1.0d0,Vm,nt,Mv,nt,0.0d0,Am,mm)
call dgemm('T','N',mp,mm,nt,1.0d0,Vp,nt,Vm,nt,0.0d0,Sv,mp)
Af = Ap
X1 = Sv
call dgesv(mp,mm,Af,mp,ipv,X1,mp,info)
if (info /= 0) error stop 'TDDFT: singular (A+B) in the subspace'
call dgemm('T','N',mm,mm,mp,1.0d0,Sv,mp,X1,mp,0.0d0,G,mm)
U = Am
call dgesv(mm,mm,U,mm,ipv,G,mm,info)
if (info /= 0) error stop 'TDDFT: singular (A-B) in the subspace'
call dgeev('N','V',mm,G,mm,ev,wi,vl,1,U,mm,wk,-1,info)
lw = int(wk(1)) ; deallocate(wk) ; allocate(wk(lw))
call dgeev('N','V',mm,G,mm,ev,wi,vl,1,U,mm,wk,lw,info)
if (info /= 0) error stop 'TDDFT: subspace diagonalization failed'
ng = 0
do q = 1,mm
   if (ev(q) > 0.0d0 .and. abs(wi(q)) <= 1.0d-8*ev(q)) then
      ng = ng + 1
      ord(ng) = q
   endif
enddo
do q = 2,ng
   qq = ord(q)
   j = q - 1
   do while (j >= 1)
      if (ev(ord(j)) >= ev(qq)) exit
      ord(j+1) = ord(j)
      j = j - 1
   enddo
   ord(j+1) = qq
enddo
if (ng < min(nwork, mm)) then
   td_unstable = .true.
   td_unstable_val = 0.0d0
   do q = 1,mm
      if (abs(wi(q)) <= 1.0d-8*abs(ev(q))) td_unstable_val = min(td_unstable_val, ev(q))
   enddo
   deallocate(Ap, Am, Sv, X1, G, U, ev, wi, vl, wk, tv, ord, ipv, Af)
   return
endif
td_neg_min = 0.0d0
do q = 1,mm
   if (abs(wi(q)) <= 1.0d-8*abs(ev(q)) .and. ev(q) < td_neg_min) td_neg_min = ev(q)
enddo
td_pos_nth = ev(ord(max(1, min(nroots, ng))))
Cc = 0.0d0 ; Dd = 0.0d0
wr = huge(1.0d0)
do q = 1,min(nkeep, ng)
   qq = ord(q)
   wr(q) = 1.0d0/sqrt(ev(qq))
   Dd(1:mm,q) = U(:,qq)
   call dgemv('N',mp,mm,wr(q),X1,mp,Dd(1,q),1,0.0d0,Cc(1,q),1)
   call dgemv('N',mp,mm,1.0d0,Sv,mp,Dd(1,q),1,0.0d0,tv,1)
   cn = abs(sum(Cc(1:mp,q)*tv))
   if (cn == 0.0d0) error stop 'TDDFT: (X+Y)^T(X-Y) vanishes in the subspace'
   Cc(1:mp,q) = Cc(1:mp,q)/sqrt(cn)
   Dd(1:mm,q) = Dd(1:mm,q)/sqrt(cn)
enddo
deallocate(Ap, Am, Sv, X1, G, U, ev, wi, vl, wk, tv, ord, ipv, Af)
end subroutine subspace

subroutine collapse()
call shrink(Vp, P, mp, Cc)
call shrink(Vm, Mv, mm, Dd)
end subroutine collapse

subroutine shrink(Vs, Ws, ms, Co)
real(8),intent(inout) :: Vs(nt,maxsub), Ws(nt,maxsub)
integer,intent(inout) :: ms
real(8),intent(in) :: Co(maxsub,nkeep)
real(8),allocatable :: S(:,:), Wk(:,:)
integer :: q, qq, j, nk
real(8) :: cn
allocate(S(ms,nkeep))
nk = 0
do q = 1,min(nkeep, ms)
   S(:,nk+1) = Co(1:ms,q)
   do qq = 1,2
      do j = 1,nk
         S(:,nk+1) = S(:,nk+1) - sum(S(:,j)*S(:,nk+1))*S(:,j)
      enddo
   enddo
   cn = sqrt(sum(S(:,nk+1)**2))
   if (cn < 1.0d-8) cycle
   nk = nk + 1
   S(:,nk) = S(:,nk)/cn
enddo
allocate(Wk(nt,nk))
call dgemm('N','N',nt,nk,ms,1.0d0,Vs,nt,S,ms,0.0d0,Wk,nt)
Vs(:,1:nk) = Wk
call dgemm('N','N',nt,nk,ms,1.0d0,Ws,nt,S,ms,0.0d0,Wk,nt)
Ws(:,1:nk) = Wk
deallocate(Wk, S)
ms = nk
end subroutine shrink

subroutine add_dirs(Vs, k, n, Tn, scl)
real(8),intent(inout) :: Vs(nt,maxsub), Tn(nt,*), scl(maxsub)
integer,intent(inout) :: k
integer,intent(in) :: n
real(8),allocatable :: c(:), Cb(:,:), nr0(:)
real(8) :: nrm
integer :: q, pass, k0, j
if (n < 1) return
k0 = k
allocate(c(maxsub), nr0(n))
do q = 1,n
   nr0(q) = sqrt(sum(Tn(:,q)**2))
enddo
if (k0 > 0) then
   allocate(Cb(k0,n))
   do pass = 1,2
      call dgemm('T','N',k0,n,nt,1.0d0,Vs,nt,Tn,nt,0.0d0,Cb,k0)
      call dgemm('N','N',nt,n,k0,-1.0d0,Vs,nt,Cb,k0,1.0d0,Tn,nt)
   enddo
   deallocate(Cb)
endif
do q = 1,n
   if (nr0(q) == 0.0d0 .or. k >= maxsub) cycle
   if (k > k0) then
      do pass = 1,2
         call dgemv('T',nt,k-k0,1.0d0,Vs(1,k0+1),nt,Tn(1,q),1,0.0d0,c,1)
         call dgemv('N',nt,k-k0,-1.0d0,Vs(1,k0+1),nt,c,1,1.0d0,Tn(1,q),1)
      enddo
   endif
   nrm = sqrt(sum(Tn(:,q)**2))
   if (nrm < 1.0d-6*nr0(q)) cycle
   k = k + 1
   Vs(:,k) = Tn(:,q)/nrm
   scl(k) = max(1.0d-4, min(1.0d0, nrm))
enddo
deallocate(c, nr0)
end subroutine add_dirs

end subroutine tddft_hybrid

subroutine tddft_hybrid1(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nroots, &
                        w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
use mod_cpscf, only: cpscf_apb_amb1_rhf, cpscf_response_begin, cpscf_response_end
use mod_cpscf, only: cphf_t_vj, cphf_t_vk, cphf_t_vxc
use MOL_info, only: engine_verbose
use omp_lib, only: omp_get_wtime
use mod_xc_hess, only: xc_fxc_ngc
implicit none
integer,intent(in) :: nConts, nocc, nroots
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: w(nroots), osc(nroots), tdip(3,nroots)
real(8),intent(out) :: XpY(nConts-nocc,nocc,nroots), XmY(nConts-nocc,nocc,nroots)
integer,intent(out) :: napply
real(8),intent(in) :: conv_thr
logical,intent(in) :: robust
integer :: nvir, nt, nwork, ngs, maxsub, m, it, k, r, i, a, info, nnew, kk, imin, j
real(8),allocatable :: dl(:), V(:,:), P(:,:), Mv(:,:), Ap(:,:), Am(:,:), ev(:), U(:,:)
real(8),allocatable :: Cc(:,:), Dd(:,:), wr(:), R1(:), R2(:), T(:,:), S(:,:), Wk(:,:)
logical,allocatable :: conv(:), used(:)
real(8) :: om, den, nrm, nrm0, rmax, t0, tapp, rn

nvir = nConts - nocc
nt = nvir*nocc
allocate(dl(nt))
do i = 1,nocc
   do a = 1,nvir
      dl(a+(i-1)*nvir) = eLev_a(nocc+a) - eLev_a(i)
   enddo
enddo
if (minval(dl) <= 0.0d0) error stop 'TDDFT: a virtual orbital lies below an occupied one'
nwork = min(nt, nroots + 2)
if (robust) nwork = min(nt, nroots + max(nroots, 6))
ngs = min(nt, max(3*nroots, 2*nwork))
maxsub = min(nt, max(ngs + 4*nwork, 8*nwork, 60))
allocate(V(nt,maxsub), P(nt,maxsub), Mv(nt,maxsub), conv(nwork), wr(nwork), Cc(maxsub,nwork), Dd(maxsub,nwork))
allocate(R1(nt), R2(nt), T(nt,2*nwork), used(nt))
m = ngs
V(:,1:m) = 0.0d0
used = .false.
do k = 1,m
   imin = minloc(dl, 1, mask=.not. used)
   used(imin) = .true.
   V(imin,k) = 1.0d0
enddo
deallocate(used)

block
   use mod_xc_hess, only: xc_resp_td_light
   use mod_cpscf, only: cpscf_triplet
   call xc_resp_td_light(.true., cpscf_triplet)
end block
call cpscf_response_begin(nConts, Pa, Pb, HF_exchange_frac)
call cosx_set_resp_grid('td')
xc_fxc_ngc = 512
napply = 0
t0 = omp_get_wtime()
tapp = 0.0d0
call apply_pm(1, m)
conv = .false.
do it = 1,200
   call subspace()
   nnew = 0
   rmax = 0.0d0
   do r = 1,nwork
      om = wr(r)
      call dgemv('N',nt,m,1.0d0,P,nt,Cc(1,r),1,0.0d0,R1,1)
      call dgemv('N',nt,m,-om,V,nt,Dd(1,r),1,1.0d0,R1,1)
      call dgemv('N',nt,m,1.0d0,Mv,nt,Dd(1,r),1,0.0d0,R2,1)
      call dgemv('N',nt,m,-om,V,nt,Cc(1,r),1,1.0d0,R2,1)
      rn = max(sqrt(sum(R1*R1)), sqrt(sum(R2*R2)))
      conv(r) = rn < conv_thr
      if (r <= nroots) rmax = max(rmax, rn)
      if (conv(r)) cycle
      do k = 1,nt
         den = dl(k)*dl(k) - om*om
         if (abs(den) < 1.0d-8) den = sign(1.0d-8, den)
         T(k,nnew+1) = -(dl(k)*R1(k) + om*R2(k))/den
         T(k,nnew+2) = -(dl(k)*R2(k) + om*R1(k))/den
      enddo
      nnew = nnew + 2
   enddo
   if (engine_verbose >= 2) print '(A,I4,A,I5,A,I4,A,ES9.2,A,F9.2,A)', '[tddft] it', it, '  space', m, &
      '  conv', count(conv(1:nroots)), '  max res', rmax, '  t', omp_get_wtime() - t0, ' s'
   if (all(conv(1:nroots))) exit
   if (m + nnew > maxsub) call collapse()
   k = m
   do r = 1,nnew
      nrm0 = sqrt(sum(T(:,r)**2))
      if (nrm0 == 0.0d0) cycle
      do kk = 1,2
         call orth_against(k, T(:,r))
      enddo
      nrm = sqrt(sum(T(:,r)**2))
      if (nrm < 1.0d-6*nrm0 .or. k >= maxsub) cycle
      k = k + 1
      V(:,k) = T(:,r)/nrm
   enddo
   if (k == m) then
      print '(A)', 'WARNING: TDDFT Davidson made no new directions; stopping with the current roots'
      exit
   endif
   call apply_pm(m+1, k)
   m = k
enddo
if (.not. all(conv(1:nroots))) print '(A)', 'WARNING: TDDFT Davidson did not converge every root'
if (engine_verbose >= 2) print '(A,F9.2,A,F9.2,A,F9.2,A,F9.2,A,F9.2)', '[tddft] total', omp_get_wtime() - t0, &
   '  apply', tapp, '  J', cphf_t_vj, '  K', cphf_t_vk, '  XC', cphf_t_vxc
call td_end()
xc_fxc_ngc = 128
do r = 1,nroots
   w(r) = wr(r)
   call dgemv('N',nt,m,1.0d0,V,nt,Cc(1,r),1,0.0d0,R1,1)
   XpY(:,:,r) = reshape(R1, (/nvir,nocc/))
   call dgemv('N',nt,m,1.0d0,V,nt,Dd(1,r),1,0.0d0,R1,1)
   XmY(:,:,r) = reshape(R1, (/nvir,nocc/))
enddo
call tddft_props(nConts, nocc, C_a, nroots, w, XpY, XmY, tdip, osc)
deallocate(dl, V, P, Mv, conv, wr, Cc, Dd, R1, R2, T)

contains

subroutine apply_pm(j1, j2)
integer,intent(in) :: j1, j2
real(8),allocatable :: Xb(:,:,:), Qp(:,:,:), Qm(:,:,:)
integer :: nb, jj
real(8) :: ta
ta = omp_get_wtime()
nb = j2 - j1 + 1
allocate(Xb(nvir,nocc,nb), Qp(nvir,nocc,nb), Qm(nvir,nocc,nb))
do jj = 1,nb
   Xb(:,:,jj) = reshape(V(:,j1+jj-1), (/nvir,nocc/))
enddo
call cpscf_apb_amb1_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, nb, Xb, Qp, Qm)
do jj = 1,nb
   P(:,j1+jj-1) = reshape(Qp(:,:,jj), (/nt/))
   Mv(:,j1+jj-1) = reshape(Qm(:,:,jj), (/nt/))
enddo
napply = napply + nb
tapp = tapp + (omp_get_wtime() - ta)
deallocate(Xb, Qp, Qm)
end subroutine apply_pm

subroutine subspace()
real(8),allocatable :: wk(:), wi(:), vl(:,:), Z(:,:)
integer :: lw, q, qq
integer,allocatable :: ord(:)
real(8) :: cn
allocate(Ap(m,m), Am(m,m), ev(m), U(m,m), wk(1), wi(m), vl(1,1), Z(m,m), ord(m))
call dgemm('T','N',m,m,nt,1.0d0,V,nt,P,nt,0.0d0,Ap,m)
call dgemm('T','N',m,m,nt,1.0d0,V,nt,Mv,nt,0.0d0,Am,m)
call dgemm('N','N',m,m,m,1.0d0,Am,m,Ap,m,0.0d0,Z,m)
call dgeev('N','V',m,Z,m,ev,wi,vl,1,U,m,wk,-1,info)
lw = int(wk(1)) ; deallocate(wk) ; allocate(wk(lw))
call dgeev('N','V',m,Z,m,ev,wi,vl,1,U,m,wk,lw,info)
if (info /= 0) error stop 'TDDFT: subspace diagonalization failed'
do q = 1,m
   ord(q) = q
enddo
do q = 2,m
   qq = ord(q)
   j = q - 1
   do while (j >= 1)
      if (ev(ord(j)) <= ev(qq)) exit
      ord(j+1) = ord(j)
      j = j - 1
   enddo
   ord(j+1) = qq
enddo
if (ev(ord(1)) <= 0.0d0) error stop 'TDDFT: non-positive w^2 in the subspace (unstable reference?)'
do q = 1,min(nwork, m)
   qq = ord(q)
   wr(q) = sqrt(ev(qq))
   Cc(1:m,q) = U(:,qq)
   call dgemv('N',m,m,1.0d0,Ap,m,Cc(1,q),1,0.0d0,Dd(1,q),1)
   cn = sum(Cc(1:m,q)*Dd(1:m,q))/wr(q)
   if (cn <= 0.0d0) error stop 'TDDFT: (X+Y)^T(X-Y) not positive in the subspace'
   Cc(1:m,q) = Cc(1:m,q)/sqrt(cn)
   Dd(1:m,q) = Dd(1:m,q)/(sqrt(cn)*wr(q))
enddo
deallocate(Ap, Am, ev, U, wk, wi, vl, Z, ord)
end subroutine subspace

subroutine collapse()
integer :: q, nk, qq
real(8) :: cn
allocate(S(m,2*nwork))
nk = 0
do q = 1,2*nwork
   if (q <= nwork) then
      S(:,nk+1) = Cc(1:m,q)
   else
      S(:,nk+1) = Dd(1:m,q-nwork)
   endif
   do qq = 1,2
      do j = 1,nk
         S(:,nk+1) = S(:,nk+1) - sum(S(:,j)*S(:,nk+1))*S(:,j)
      enddo
   enddo
   cn = sqrt(sum(S(:,nk+1)**2))
   if (cn < 1.0d-8) cycle
   nk = nk + 1
   S(:,nk) = S(:,nk)/cn
enddo
allocate(Wk(nt,nk))
call dgemm('N','N',nt,nk,m,1.0d0,V,nt,S,m,0.0d0,Wk,nt)
V(:,1:nk) = Wk
call dgemm('N','N',nt,nk,m,1.0d0,P,nt,S,m,0.0d0,Wk,nt)
P(:,1:nk) = Wk
call dgemm('N','N',nt,nk,m,1.0d0,Mv,nt,S,m,0.0d0,Wk,nt)
Mv(:,1:nk) = Wk
deallocate(Wk, S)
m = nk
end subroutine collapse

subroutine orth_against(n, t)
integer,intent(in) :: n
real(8),intent(inout) :: t(nt)
real(8),allocatable :: c(:)
if (n < 1) return
allocate(c(n))
call dgemv('T',nt,n,1.0d0,V,nt,t,1,0.0d0,c,1)
call dgemv('N',nt,n,-1.0d0,V,nt,c,1,1.0d0,t,1)
deallocate(c)
end subroutine orth_against

end subroutine tddft_hybrid1

subroutine tddft_polar_hybrid(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, omega, alpha, napply)
use mod_cpscf, only: cpscf_apb_amb_rhf, cpscf_response_begin, cpscf_response_end
use mod_integrals, only: integrals_dipole_ao
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac, omega
real(8),intent(out) :: alpha(3,3)
integer,intent(out) :: napply
integer :: nvir, nt, maxsub, mp, mm, it, k, r, i, a, x, info, np, nm, kp, km, nb
real(8),allocatable :: dl(:), v(:,:), Vp(:,:), Vm(:,:), P(:,:), Mv(:,:), U(:,:), W(:,:), R1(:), R2(:), &
                       Tp(:,:), Tm(:,:), Rao(:,:,:), Tq(:,:), Bm(:,:), rhs(:,:), Cs(:,:)
integer,allocatable :: ipv(:)
logical :: conv(3)
real(8) :: den, rn

nvir = nConts - nocc
nt = nvir*nocc
allocate(dl(nt), v(nt,3), Rao(nConts,nConts,3), Tq(nConts,nocc))
do i = 1,nocc
   do a = 1,nvir
      dl(a+(i-1)*nvir) = eLev_a(nocc+a) - eLev_a(i)
   enddo
enddo
if (omega >= minval(dl)) error stop 'dynamic polarizability: omega is not below the HOMO-LUMO gap (resonance region)'
call integrals_dipole_ao(nConts, Rao)
do x = 1,3
   call dgemm('N','N',nConts,nocc,nConts,1.0d0,Rao(1,1,x),nConts,C_a,nConts,0.0d0,Tq,nConts)
   call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,Tq,nConts,0.0d0,v(1,x),nvir)
enddo
deallocate(Rao, Tq)
maxsub = min(nt, 120)
allocate(Vp(nt,maxsub), Vm(nt,maxsub), P(nt,maxsub), Mv(nt,maxsub), U(nt,3), W(nt,3), R1(nt), R2(nt), &
         Tp(nt,3), Tm(nt,3))
mp = 0 ; mm = 0
do x = 1,3
   Tp(:,x) = -v(:,x)*dl/(dl*dl - omega*omega)
   Tm(:,x) = -v(:,x)*omega/(dl*dl - omega*omega)
enddo
block
   character(8) :: eg
   call get_environment_variable("ENGINE_POLAR_GRID", eg)
   call cpscf_response_begin(nConts, Pa, Pb, HF_exchange_frac, scf_grid=(trim(eg) == "scf"))
end block
napply = 0
kp = 0 ; call add_dirs(Vp, kp, 3, Tp)
km = 0 ; call add_dirs(Vm, km, 3, Tm)
call apply_pm(1, kp, 1, km)
mp = kp ; mm = km
conv = .false.
do it = 1,100
   call solve_reduced()
   np = 0 ; nm = 0
   do x = 1,3
      call dgemv('N',nt,mp,1.0d0,P,nt,Cs(1,x),1,0.0d0,R1,1)
      if (mm > 0) call dgemv('N',nt,mm,-omega,Vm,nt,Cs(mp+1,x),1,1.0d0,R1,1)
      R1 = R1 + v(:,x)
      R2 = 0.0d0
      if (mm > 0) call dgemv('N',nt,mm,1.0d0,Mv,nt,Cs(mp+1,x),1,0.0d0,R2,1)
      call dgemv('N',nt,mp,-omega,Vp,nt,Cs(1,x),1,1.0d0,R2,1)
      rn = max(sqrt(sum(R1*R1)), sqrt(sum(R2*R2)))
      conv(x) = rn < 1.0d-7
      if (conv(x)) cycle
      np = np + 1 ; nm = nm + 1
      do k = 1,nt
         den = dl(k)*dl(k) - omega*omega
         Tp(k,np) = -(dl(k)*R1(k) + omega*R2(k))/den
         Tm(k,nm) = -(dl(k)*R2(k) + omega*R1(k))/den
      enddo
   enddo
   if (all(conv)) exit
   if (mp + np > maxsub .or. mm + nm > maxsub) error stop 'dynamic polarizability: trial space exhausted'
   kp = mp ; call add_dirs(Vp, kp, np, Tp)
   km = mm ; call add_dirs(Vm, km, nm, Tm)
   if (kp == mp .and. km == mm) exit
   call apply_pm(mp+1, kp, mm+1, km)
   mp = kp ; mm = km
enddo
if (.not. all(conv)) print '(A)', 'WARNING: dynamic polarizability response not converged'
if (omega > 0.0d0) call check_below_first_pole()
call cpscf_response_end()
call dgemm('N','N',nt,3,mp,1.0d0,Vp,nt,Cs,mp+mm,0.0d0,U,nt)
do x = 1,3
   do r = 1,3
      alpha(x,r) = -4.0d0*sum(v(:,x)*U(:,r))
   enddo
enddo
alpha = 0.5d0*(alpha + transpose(alpha))
deallocate(dl, v, Vp, Vm, P, Mv, U, W, R1, R2, Tp, Tm)
if (allocated(Cs)) deallocate(Cs)

contains

subroutine check_below_first_pole()
real(8),allocatable :: Ap(:,:), Am(:,:), Sv(:,:), X1(:,:), G(:,:), ev(:), wi(:), vl(:,:), wk(:), U(:,:)
integer,allocatable :: ipv2(:)
integer :: q, lw, inf
real(8) :: wlow, best
allocate(Ap(mp,mp), Am(mm,mm), Sv(mp,mm), X1(mp,mm), G(mm,mm), U(mm,mm), ev(mm), wi(mm), vl(1,1), wk(1), ipv2(max(mp,mm)))
call dgemm('T','N',mp,mp,nt,1.0d0,Vp,nt,P,nt,0.0d0,Ap,mp)
call dgemm('T','N',mm,mm,nt,1.0d0,Vm,nt,Mv,nt,0.0d0,Am,mm)
call dgemm('T','N',mp,mm,nt,1.0d0,Vp,nt,Vm,nt,0.0d0,Sv,mp)
X1 = Sv
call dgesv(mp,mm,Ap,mp,ipv2,X1,mp,inf)
if (inf /= 0) return
call dgemm('T','N',mm,mm,mp,1.0d0,Sv,mp,X1,mp,0.0d0,G,mm)
U = Am
call dgesv(mm,mm,U,mm,ipv2,G,mm,inf)
if (inf /= 0) return
call dgeev('N','N',mm,G,mm,ev,wi,vl,1,vl,1,wk,-1,inf)
lw = int(wk(1)) ; deallocate(wk) ; allocate(wk(lw))
call dgeev('N','N',mm,G,mm,ev,wi,vl,1,vl,1,wk,lw,inf)
if (inf /= 0) return
best = 0.0d0
do q = 1,mm
   if (ev(q) > best .and. abs(wi(q)) <= 1.0d-8*ev(q)) best = ev(q)
enddo
if (best <= 0.0d0) return
wlow = 1.0d0/sqrt(best)
if (omega >= wlow) then
   print '(A,F10.6,A,F10.6,A)', 'ERROR: dynamic polarizability: omega = ', omega, &
         ' a.u. is not below the lowest dipole-active excitation (', wlow, ' a.u.)'
   error stop 'dynamic polarizability: omega beyond the first excitation (resonance region)'
endif
end subroutine check_below_first_pole

subroutine solve_reduced()
integer :: nn
nn = mp + mm
if (allocated(Cs)) deallocate(Cs)
allocate(Bm(nn,nn), rhs(nn,3), ipv(nn), Cs(nn,3))
call dgemm('T','N',mp,mp,nt,1.0d0,Vp,nt,P,nt,0.0d0,Bm(1,1),nn)
if (mm > 0) then
   call dgemm('T','N',mm,mm,nt,1.0d0,Vm,nt,Mv,nt,0.0d0,Bm(mp+1,mp+1),nn)
   call dgemm('T','N',mp,mm,nt,-omega,Vp,nt,Vm,nt,0.0d0,Bm(1,mp+1),nn)
   Bm(mp+1:nn,1:mp) = transpose(Bm(1:mp,mp+1:nn))
endif
call dgemm('T','N',mp,3,nt,-1.0d0,Vp,nt,v,nt,0.0d0,rhs,nn)
rhs(mp+1:nn,:) = 0.0d0
call dgesv(nn,3,Bm,nn,ipv,rhs,nn,info)
if (info /= 0) error stop 'dynamic polarizability: singular reduced system'
Cs = rhs
deallocate(Bm, rhs, ipv)
end subroutine solve_reduced

subroutine apply_pm(jp1, jp2, jm1, jm2)
integer,intent(in) :: jp1, jp2, jm1, jm2
real(8),allocatable :: Xa(:,:,:), Xb(:,:,:), Qa(:,:,:), Qb(:,:,:)
integer :: na, nbb, jj
na = max(0, jp2 - jp1 + 1) ; nbb = max(0, jm2 - jm1 + 1)
allocate(Xa(nvir,nocc,max(na,1)), Qa(nvir,nocc,max(na,1)), Xb(nvir,nocc,max(nbb,1)), Qb(nvir,nocc,max(nbb,1)))
do jj = 1,na
   Xa(:,:,jj) = reshape(Vp(:,jp1+jj-1), (/nvir,nocc/))
enddo
do jj = 1,nbb
   Xb(:,:,jj) = reshape(Vm(:,jm1+jj-1), (/nvir,nocc/))
enddo
call cpscf_apb_amb_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, na, Xa, Qa, nbb, Xb, Qb)
do jj = 1,na
   P(:,jp1+jj-1) = reshape(Qa(:,:,jj), (/nt/))
enddo
do jj = 1,nbb
   Mv(:,jm1+jj-1) = reshape(Qb(:,:,jj), (/nt/))
enddo
napply = napply + na + nbb
deallocate(Xa, Xb, Qa, Qb)
end subroutine apply_pm

subroutine add_dirs(Vs, k, n, Tn)
real(8),intent(inout) :: Vs(nt,maxsub), Tn(nt,*)
integer,intent(inout) :: k
integer,intent(in) :: n
real(8),allocatable :: c(:)
real(8) :: nrm, nrm0
integer :: q, pass
allocate(c(maxsub))
do q = 1,n
   nrm0 = sqrt(sum(Tn(:,q)**2))
   if (nrm0 == 0.0d0 .or. k >= maxsub) cycle
   do pass = 1,2
      if (k > 0) then
         call dgemv('T',nt,k,1.0d0,Vs,nt,Tn(1,q),1,0.0d0,c,1)
         call dgemv('N',nt,k,-1.0d0,Vs,nt,c,1,1.0d0,Tn(1,q),1)
      endif
   enddo
   nrm = sqrt(sum(Tn(:,q)**2))
   if (nrm < 1.0d-10*nrm0) cycle
   k = k + 1
   Vs(:,k) = Tn(:,q)/nrm
enddo
deallocate(c)
end subroutine add_dirs

end subroutine tddft_polar_hybrid

subroutine tddft_main_transitions(nvir, nocc, XpY, XmY, thr, n, ifrom, ito, wt)
implicit none
integer,intent(in) :: nvir, nocc
real(8),intent(in) :: XpY(nvir,nocc), XmY(nvir,nocc), thr
integer,intent(out) :: n, ifrom(nvir*nocc), ito(nvir*nocc)
real(8),intent(out) :: wt(nvir*nocc)
integer :: i, a, k, j, it
real(8) :: x, t
n = 0
do i = 1,nocc
   do a = 1,nvir
      x = XpY(a,i)*XmY(a,i)
      if (x >= thr) then
         n = n + 1
         ifrom(n) = i ; ito(n) = nocc + a ; wt(n) = x
      endif
   enddo
enddo
do k = 2,n
   t = wt(k) ; i = ifrom(k) ; it = ito(k)
   j = k - 1
   do while (j >= 1)
      if (wt(j) >= t) exit
      wt(j+1) = wt(j) ; ifrom(j+1) = ifrom(j) ; ito(j+1) = ito(j)
      j = j - 1
   enddo
   wt(j+1) = t ; ifrom(j+1) = i ; ito(j+1) = it
enddo
end subroutine tddft_main_transitions

subroutine tddft_uks(nConts, noa, nob, Ca, ea, Cb, eb, Pa, Pb, HF_exchange_frac, nroots, &
                     w, osc, tdip, XpY, XmY, napply, conv_thr, robust)
implicit none
integer,intent(in) :: nConts, noa, nob, nroots
real(8),intent(in) :: Ca(nConts,nConts), ea(nConts), Cb(nConts,nConts), eb(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts), HF_exchange_frac, conv_thr
real(8),intent(out) :: w(nroots), osc(nroots), tdip(3,nroots)
real(8),intent(out) :: XpY((nConts-noa)*noa+(nConts-nob)*nob,nroots), XmY((nConts-noa)*noa+(nConts-nob)*nob,nroots)
integer,intent(out) :: napply
logical,intent(in) :: robust
td_uks = .true.
td_nob = nob
if (allocated(td_Cb)) deallocate(td_Cb, td_eb)
allocate(td_Cb(nConts,nConts), td_eb(nConts))
td_Cb = Cb ; td_eb = eb
call tddft_rks(nConts, noa, Ca, ea, Pa, Pb, HF_exchange_frac, nroots, w, osc, tdip, XpY, XmY, napply, &
               conv_thr, robust)
td_uks = .false.
td_nob = 0
deallocate(td_Cb, td_eb)
end subroutine tddft_uks

pure integer function td_ntot(nConts, nocc)
integer,intent(in) :: nConts, nocc
td_ntot = (nConts - nocc)*nocc
if (td_uks) td_ntot = td_ntot + (nConts - td_nob)*td_nob
end function td_ntot

subroutine td_diag(nConts, nocc, ea, d)
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: ea(nConts)
real(8),intent(out) :: d(*)
integer :: i, a, nv, na
nv = nConts - nocc
do i = 1,nocc
   do a = 1,nv
      d(a+(i-1)*nv) = ea(nocc+a) - ea(i)
   enddo
enddo
if (.not. td_uks) return
na = nv*nocc
nv = nConts - td_nob
do i = 1,td_nob
   do a = 1,nv
      d(na+a+(i-1)*nv) = td_eb(td_nob+a) - td_eb(i)
   enddo
enddo
end subroutine td_diag

subroutine td_begin(nConts, Pa, Pb, HF_exchange_frac)
use mod_cpscf, only: cpscf_response_begin, cpscf_response_begin_uks
use mod_xc_hess, only: xc_resp_td_light
use mod_cpscf, only: cpscf_triplet
integer,intent(in) :: nConts
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts), HF_exchange_frac
call xc_resp_td_light(.true., cpscf_triplet)
if (td_uks) then
   call cpscf_response_begin_uks(nConts, Pa, Pb, HF_exchange_frac)
else
   call cpscf_response_begin(nConts, Pa, Pb, HF_exchange_frac)
endif
call cosx_set_resp_grid('td')
end subroutine td_begin

subroutine td_end()
use mod_cpscf, only: cpscf_response_end
use mod_xc_hess, only: xc_resp_td_light
call cpscf_response_end()
call cosx_set_resp_grid('')
call xc_resp_td_light(.false., .false.)
end subroutine td_end

subroutine td_props(nConts, nocc, C_a, nroots, w, XpY, XmY, tdip, osc)
integer,intent(in) :: nConts, nocc, nroots
real(8),intent(in) :: C_a(nConts,nConts), w(nroots)
real(8),intent(inout) :: XpY(td_ntot(nConts,nocc),nroots), XmY(td_ntot(nConts,nocc),nroots)
real(8),intent(out) :: tdip(3,nroots), osc(nroots)
if (td_uks) then
   call td_props_uks(nConts, nocc, C_a, nroots, w, XpY, XmY, tdip, osc)
else
   call tddft_props(nConts, nocc, C_a, nroots, w, XpY, XmY, tdip, osc)
endif
end subroutine td_props

subroutine td_props_uks(nConts, noa, Ca, nroots, w, XpY, XmY, tdip, osc)
use mod_integrals, only: integrals_dipole_ao
integer,intent(in) :: nConts, noa, nroots
real(8),intent(in) :: Ca(nConts,nConts), w(nroots)
real(8),intent(inout) :: XpY(td_ntot(nConts,noa),nroots), XmY(td_ntot(nConts,noa),nroots)
real(8),intent(out) :: tdip(3,nroots), osc(nroots)
real(8),allocatable :: Rao(:,:,:), Tm(:,:), Vmo(:,:)
integer :: nva, nvb, na, nt, x, r, imx
nva = nConts - noa ; nvb = nConts - td_nob
na = nva*noa ; nt = na + nvb*td_nob
do r = 1,nroots
   imx = maxloc(abs(XpY(:,r)), 1)
   if (XpY(imx,r) < 0.0d0) then
      XpY(:,r) = -XpY(:,r) ; XmY(:,r) = -XmY(:,r)
   endif
enddo
allocate(Rao(nConts,nConts,3), Tm(nConts,max(noa,td_nob)), Vmo(nt,3))
call integrals_dipole_ao(nConts, Rao)
do x = 1,3
   call dgemm('N','N',nConts,noa,nConts,1.0d0,Rao(1,1,x),nConts,Ca,nConts,0.0d0,Tm,nConts)
   call dgemm('T','N',nva,noa,nConts,1.0d0,Ca(1,noa+1),nConts,Tm,nConts,0.0d0,Vmo(1,x),nva)
   if (td_nob > 0) then
      call dgemm('N','N',nConts,td_nob,nConts,1.0d0,Rao(1,1,x),nConts,td_Cb,nConts,0.0d0,Tm,nConts)
      call dgemm('T','N',nvb,td_nob,nConts,1.0d0,td_Cb(1,td_nob+1),nConts,Tm,nConts,0.0d0,Vmo(na+1,x),nvb)
   endif
enddo
do r = 1,nroots
   do x = 1,3
      tdip(x,r) = sum(Vmo(:,x)*XpY(:,r))
   enddo
   osc(r) = 2.0d0/3.0d0*w(r)*sum(tdip(:,r)**2)
enddo
deallocate(Rao, Tm, Vmo)
end subroutine td_props_uks

end module mod_tddft
