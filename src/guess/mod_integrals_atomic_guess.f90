! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Builds the atomic (SAD) initial density guess used to start SCF.

submodule (mod_integrals) atomic_guess_impl
implicit none
contains

module subroutine compute_atomic_density(ia, closed_pairing, Pa_atom, Pb_atom, ok)
use MOL_info
use mod_basis_files, only: basis_file_for_label
implicit none
integer,intent(in) :: ia
logical,intent(in) :: closed_pairing
real(8),allocatable,intent(out) :: Pa_atom(:,:), Pb_atom(:,:)
logical,intent(out) :: ok

integer,parameter :: ECP_MAX_CHANNELS_L = 6, ECP_MAX_TERMS_L = 20
integer :: ecpL_l(ECP_MAX_CHANNELS_L), ecpL_nterms(ECP_MAX_CHANNELS_L)
integer :: ecpL_power(ECP_MAX_CHANNELS_L,ECP_MAX_TERMS_L)
real(8) :: ecpL_exp(ECP_MAX_CHANNELS_L,ECP_MAX_TERMS_L), ecpL_coef(ECP_MAX_CHANNELS_L,ECP_MAX_TERMS_L)
integer :: ecpL_nchan
integer,allocatable :: basECP_l(:,:)
real(8),allocatable :: envECP_l(:)
integer :: necpbas_l, Z_eff
logical :: has_ecp_l
character(len=200) :: bpath_l
character(len=40) :: ecp_label_l, tok_l
character(len=2) :: letter_l
integer :: ios_l, nt_l, ch_l, t_l, run_start_l, shidx_l, offpoint_ecp, k_l

integer :: numShell, cnt, Z_atom
integer,allocatable :: atm_l(:,:), bas_l(:,:)
real(8),allocatable :: env_l(:)
integer :: offpoint, j, k, envsize
integer :: shls(4)
real(8),allocatable :: NorVEC_l(:)
real(8),allocatable :: S_l(:,:), Hcore_l(:,:)
real(8),allocatable :: eri(:,:,:,:)
integer :: i0,j0,k0,l0,di,dj,dk,dl,n
integer,allocatable :: shell_off(:)
integer :: p,q,e1,e2,e3,e4
real(8),allocatable :: Feff(:,:), Jmat(:,:), Ka(:,:), Kb(:,:), Ptot(:,:)
real(8),allocatable :: Ablk(:,:), Sblk(:,:), Wblk(:), WORKblk(:)
real(8),allocatable :: Pa_new(:,:), Pb_new(:,:)
integer :: LWORKblk, INFOblk
integer :: iter
real(8) :: rms, remaining_elec, group_elec, occ_a_per_orbital, occ_b_per_orbital
integer :: gstart, gend
real(8),parameter :: DEGEN_TOL = 1.0d-4
integer,parameter :: MAXITER = 40
integer,external :: CINTcgto_cart
real(8),external :: CINTgto_norm

ok = .false.

numShell = atoms(ia)%nshell
cnt = atoms(ia)%nconts
Z_atom = atoms(ia)%charge
if (numShell .le. 0 .or. cnt .le. 0) return

necpbas_l = 0
has_ecp_l = .false.
Z_eff = Z_atom
if (len_trim(atoms(ia)%ecpbase) .gt. 0) then
   ecp_label_l = trim(atoms(ia)%ecpbase)//".ecp"
   bpath_l = basis_file_for_label(ecp_label_l)
   if (len_trim(bpath_l) .gt. 0) then
      open(unit=203, file=trim(bpath_l))
      do
         read(203,*,iostat=ios_l) tok_l
         if (ios_l /= 0) exit
         if (trim(tok_l) .eq. trim(ecp_label_l)) then
            has_ecp_l = .true.
            exit
         endif
      enddo
      if (has_ecp_l) then
         block
            integer :: dum_core
            read(203,*) dum_core
            Z_eff = Z_atom - dum_core
         end block
         ecpL_nchan = 0
         do
            read(203,*,iostat=ios_l) letter_l, nt_l
            if (ios_l /= 0) exit
            ecpL_nchan = ecpL_nchan + 1
            select case (trim(letter_l))
            case ("ul"); ecpL_l(ecpL_nchan) = -1
            case ("s");  ecpL_l(ecpL_nchan) = 0
            case ("p");  ecpL_l(ecpL_nchan) = 1
            case ("d");  ecpL_l(ecpL_nchan) = 2
            case ("f");  ecpL_l(ecpL_nchan) = 3
            case ("g");  ecpL_l(ecpL_nchan) = 4
            end select
            ecpL_nterms(ecpL_nchan) = nt_l
            do t_l = 1,nt_l
               read(203,*) ecpL_power(ecpL_nchan,t_l), ecpL_exp(ecpL_nchan,t_l), ecpL_coef(ecpL_nchan,t_l)
            enddo
         enddo
      endif
      close(203)
   endif
endif
if (has_ecp_l) then
   if (.not. (atoms(ia)%ecpCoreElec .eq. 2  .or. atoms(ia)%ecpCoreElec .eq. 10 .or. &
              atoms(ia)%ecpCoreElec .eq. 18 .or. atoms(ia)%ecpCoreElec .eq. 36 .or. &
              atoms(ia)%ecpCoreElec .eq. 54 .or. atoms(ia)%ecpCoreElec .eq. 86)) then
      return
   endif
endif
if (has_ecp_l) then
   do ch_l = 1,ecpL_nchan
      nt_l = ecpL_nterms(ch_l)
      necpbas_l = necpbas_l + 1
      do t_l = 2,nt_l
         if (ecpL_power(ch_l,t_l) /= ecpL_power(ch_l,t_l-1)) necpbas_l = necpbas_l + 1
      enddo
   enddo
endif

allocate(atm_l(6,1))
atm_l = 0
atm_l(1,1) = Z_eff
atm_l(2,1) = 20

envsize = 24
do j = 1,numShell
   envsize = envsize + 2*atoms(ia)%shell(j)%nGauss
enddo
allocate(env_l(envsize))
env_l = 0.0d0

allocate(bas_l(8,numShell))
bas_l = 0
offpoint = 24
do j = 1,numShell
   bas_l(1,j) = 0
   bas_l(2,j) = atoms(ia)%shell(j)%angMoment
   bas_l(3,j) = atoms(ia)%shell(j)%nGauss
   bas_l(4,j) = 1
   bas_l(6,j) = offpoint
   do k = 1,atoms(ia)%shell(j)%nGauss
      env_l(offpoint+k) = atoms(ia)%shell(j)%exponents(k)
   enddo
   offpoint = offpoint + atoms(ia)%shell(j)%nGauss
   bas_l(7,j) = offpoint
   do k = 1,atoms(ia)%shell(j)%nGauss
      env_l(offpoint+k) = atoms(ia)%shell(j)%contrCoeff(k) &
                           *CINTgto_norm(atoms(ia)%shell(j)%angMoment, atoms(ia)%shell(j)%exponents(k))
   enddo
   offpoint = offpoint + atoms(ia)%shell(j)%nGauss
enddo

if (has_ecp_l) then
   allocate(basECP_l(8, numShell+necpbas_l))
   allocate(envECP_l(size(env_l) + 2*ECP_MAX_TERMS_L*max(necpbas_l,1)))
   basECP_l = 0
   envECP_l = 0.0d0
   basECP_l(1:8,1:numShell) = bas_l(1:8,1:numShell)
   envECP_l(1:size(env_l)) = env_l(1:size(env_l))
   offpoint_ecp = size(env_l)
   shidx_l = numShell
   do ch_l = 1,ecpL_nchan
      nt_l = ecpL_nterms(ch_l)
      t_l = 1
      do while (t_l <= nt_l)
         run_start_l = t_l
         do while (t_l < nt_l)
            if (ecpL_power(ch_l,t_l+1) /= ecpL_power(ch_l,run_start_l)) exit
            t_l = t_l + 1
         enddo
         shidx_l = shidx_l + 1
         basECP_l(1,shidx_l) = 0
         basECP_l(2,shidx_l) = ecpL_l(ch_l)
         basECP_l(3,shidx_l) = t_l - run_start_l + 1
         basECP_l(4,shidx_l) = ecpL_power(ch_l,run_start_l)
         basECP_l(5,shidx_l) = 0
         basECP_l(6,shidx_l) = offpoint_ecp
         do k_l = run_start_l,t_l
            envECP_l(offpoint_ecp+1+(k_l-run_start_l)) = ecpL_exp(ch_l,k_l)
         enddo
         offpoint_ecp = offpoint_ecp + (t_l-run_start_l+1)
         basECP_l(7,shidx_l) = offpoint_ecp
         do k_l = run_start_l,t_l
            envECP_l(offpoint_ecp+1+(k_l-run_start_l)) = ecpL_coef(ch_l,k_l)
         enddo
         offpoint_ecp = offpoint_ecp + (t_l-run_start_l+1)
         t_l = t_l + 1
      enddo
   enddo
   envECP_l(18) = 0.0d0
   envECP_l(19) = dble(numShell)
   envECP_l(20) = dble(necpbas_l)
endif

allocate(shell_off(0:numShell-1))
shell_off(0) = 0
do j0 = 1,numShell-1
   shell_off(j0) = shell_off(j0-1) + cgto_engine(j0-1,bas_l)
enddo

allocate(NorVEC_l(cnt))
if (engine_puream) then
   block
      real(8),allocatable :: buf1e(:,:)
      real(8) :: nfac
      integer :: aoidx, dj2
      do j = 0,numShell-1
         shls(1) = j; shls(2) = j
         dj2 = cgto_engine(j, bas_l)
         allocate(buf1e(dj2,dj2))
         call ovlp1e_engine(buf1e, shls, atm_l, 1, bas_l, numShell, env_l, 0_8)
         nfac = 1.0d0/sqrt(buf1e(1,1))
         do aoidx = 1,dj2
            NorVEC_l(shell_off(j)+aoidx) = nfac
         enddo
         deallocate(buf1e)
      enddo
   end block
else
   block
      real(8),allocatable :: buf1e(:,:)
      do j = 0,numShell-1
         shls(1) = j; shls(2) = j
         di = CINTcgto_cart(j, bas_l)
         allocate(buf1e(di,di))
         call cint1e_ovlp_cart(buf1e, shls, atm_l, 1, bas_l, numShell, env_l, 0_8)
         call Normal(j, di, bas_l, numShell, buf1e, NorVEC_l, cnt)
         deallocate(buf1e)
      enddo
   end block
endif

allocate(S_l(cnt,cnt), Hcore_l(cnt,cnt))
S_l = 0.0d0; Hcore_l = 0.0d0
block
   real(8),allocatable :: bufS(:,:), bufT(:,:), bufV(:,:), bufECP(:,:)
   integer :: dims_arr_l(2)
   n = 0
   do j = 0,numShell-1
      n = max(n, cgto_engine(j,bas_l))
   enddo
   allocate(bufS(n,n),bufT(n,n),bufV(n,n))
   if (has_ecp_l) allocate(bufECP(n,n))
   do i0 = 0,numShell-1
      do j0 = i0,numShell-1
         shls(1)=i0; shls(2)=j0
         di = cgto_engine(i0,bas_l); dj = cgto_engine(j0,bas_l)
         call ovlp1e_engine(bufS(1:di,1:dj), shls, atm_l, 1, bas_l, numShell, env_l, 0_8)
         call store1e(shls,di,dj,bas_l,numShell,bufS(1:di,1:dj),NorVEC_l,cnt,S_l)
         call kin1e_engine(bufT(1:di,1:dj), shls, atm_l, 1, bas_l, numShell, env_l, 0_8)
         call nuc1e_engine(bufV(1:di,1:dj), shls, atm_l, 1, bas_l, numShell, env_l, 0_8)
         if (has_ecp_l) then
            dims_arr_l(1) = di
            dims_arr_l(2) = dj
            call ecp1e_engine(bufECP(1:di,1:dj), dims_arr_l, shls, atm_l, &
                               1, basECP_l, numShell+necpbas_l, envECP_l, 0_8, 0_8)
            call store1e(shls,di,dj,bas_l,numShell, &
                         bufT(1:di,1:dj)+bufV(1:di,1:dj)+bufECP(1:di,1:dj),NorVEC_l,cnt,Hcore_l)
         else
            call store1e(shls,di,dj,bas_l,numShell,bufT(1:di,1:dj)+bufV(1:di,1:dj),NorVEC_l,cnt,Hcore_l)
         endif
      enddo
   enddo
   deallocate(bufS,bufT,bufV)
   if (allocated(bufECP)) deallocate(bufECP)
end block

if (.not. engine_puream) then
block
   real(8),parameter :: TRACE_PENALTY = 1.0d4
   real(8) :: cT(cnt)
   integer :: xx_idx, yy_idx, zz_idx
   do j0 = 0,numShell-1
      if (atoms(ia)%shell(j0+1)%angMoment .eq. 2) then
         cT = 0.0d0
         xx_idx = shell_off(j0)+1
         yy_idx = shell_off(j0)+4
         zz_idx = shell_off(j0)+6
         cT(xx_idx) = 1.0d0/sqrt(5.0d0)
         cT(yy_idx) = 1.0d0/sqrt(5.0d0)
         cT(zz_idx) = 1.0d0/sqrt(5.0d0)
         Hcore_l = Hcore_l + TRACE_PENALTY*spread(cT,2,cnt)*spread(cT,1,cnt)
      endif
   enddo
end block
endif

allocate(eri(cnt,cnt,cnt,cnt))
eri = 0.0d0
block
   real(8),allocatable :: buf2e(:,:,:,:)
   allocate(buf2e(n,n,n,n))
   do i0=0,numShell-1
    do j0=0,numShell-1
     do k0=0,numShell-1
      do l0=0,numShell-1
        shls = (/i0,j0,k0,l0/)
        di=cgto_engine(i0,bas_l); dj=cgto_engine(j0,bas_l)
        dk=cgto_engine(k0,bas_l); dl=cgto_engine(l0,bas_l)
        call twoe_engine(buf2e(1:di,1:dj,1:dk,1:dl), shls, atm_l, 1, bas_l, numShell, env_l, 0_8)
        do p = 1,di
         do q = 1,dj
          do k = 1,dk
           do j = 1,dl
              e1 = shell_off(i0)+p; e2 = shell_off(j0)+q
              e3 = shell_off(k0)+k; e4 = shell_off(l0)+j
              eri(e1,e2,e3,e4) = buf2e(p,q,k,j)*NorVEC_l(e1)*NorVEC_l(e2)*NorVEC_l(e3)*NorVEC_l(e4)
           enddo
          enddo
         enddo
        enddo
      enddo
     enddo
    enddo
   enddo
   deallocate(buf2e)
end block

allocate(Ablk(cnt,cnt), Sblk(cnt,cnt), Wblk(cnt))
LWORKblk = 1 + 6*cnt + cnt**2
allocate(WORKblk(LWORKblk))
allocate(Pa_atom(cnt,cnt), Pb_atom(cnt,cnt))
allocate(Pa_new(cnt,cnt), Pb_new(cnt,cnt))
allocate(Feff(cnt,cnt), Jmat(cnt,cnt), Ka(cnt,cnt), Kb(cnt,cnt), Ptot(cnt,cnt))
Pa_atom = 0.0d0; Pb_atom = 0.0d0
Feff = Hcore_l

block
   real(8),allocatable :: X_l(:,:)
   real(4),allocatable :: diis_Fa(:,:,:), diis_Fb(:,:,:)
   real(4),allocatable :: diis_Pa(:,:,:), diis_Pb(:,:,:)
   real(4),allocatable :: diis_ea(:,:,:), diis_eb(:,:,:)
   integer :: DIIS_MAXl, diis_nl, multi_l
   logical :: diis_active_l

   block
      real(8),allocatable :: S_e(:), S_temp(:,:), S_X(:,:), WORKs(:)
      integer :: LWORKs, INFOs, ii
      real(8),parameter :: X_EVAL_FLOOR = 1.0d-6
      allocate(S_e(cnt), S_temp(cnt,cnt), S_X(cnt,cnt))
      S_temp = S_l
      LWORKs = 1 + 6*cnt + cnt**2
      allocate(WORKs(LWORKs))
      call DSYEV('V','U',cnt,S_temp,cnt,S_e,WORKs,LWORKs,INFOs)
      S_X = 0.0d0
      do ii = 1,cnt
         if (S_e(ii) .gt. X_EVAL_FLOOR*S_e(cnt)) S_X(ii,ii) = S_e(ii)**(-0.5d0)
      enddo
      allocate(X_l(cnt,cnt))
      X_l = matmul(S_temp,matmul(S_X,transpose(S_temp)))
      deallocate(S_e,S_temp,S_X,WORKs)
   end block

   DIIS_MAXl = min(8,cnt)
   allocate(diis_Fa(cnt,cnt,DIIS_MAXl), diis_Fb(cnt,cnt,DIIS_MAXl))
   allocate(diis_Pa(cnt,cnt,DIIS_MAXl), diis_Pb(cnt,cnt,DIIS_MAXl))
   allocate(diis_ea(cnt,cnt,DIIS_MAXl), diis_eb(cnt,cnt,DIIS_MAXl))
   diis_nl = 0
   multi_l = merge(1,3,closed_pairing)

   do iter = 1,MAXITER
      if (iter .gt. 1) then
         call diis_extrapolate(cnt, multi_l, iter, Feff, Feff, Pa_atom, Pb_atom, S_l, X_l, &
                                DIIS_MAXl, diis_nl, diis_Fa, diis_Fb, diis_Pa, diis_Pb, &
                                diis_ea, diis_eb, .false., diis_active_l)
      endif
      Pa_new = 0.0d0; Pb_new = 0.0d0
      if (engine_puream) then
         call puream_atomic_fill(cnt, numShell, shell_off, bas_l, Z_atom, atoms(ia)%ecpCoreElec, &
                                  closed_pairing, Feff, S_l, Pa_new, Pb_new, INFOblk)
         if (INFOblk .ne. 0) then
            deallocate(Pa_atom,Pb_atom)
            ok = .false.
            return
         endif
      else
      Ablk = Feff
      Sblk = S_l
      call DSYGV(1,'V','U',cnt,Ablk,cnt,Sblk,cnt,Wblk,WORKblk,LWORKblk,INFOblk)
      if (INFOblk .ne. 0) then
         deallocate(Pa_atom,Pb_atom)
         ok = .false.
         return
      endif

      remaining_elec = real(Z_eff,8)
      gstart = 1
      do while (gstart .le. cnt .and. remaining_elec .gt. 1.0d-12)
         gend = gstart
         do while (gend .lt. cnt)
            if (abs(Wblk(gend+1)-Wblk(gstart)) .gt. DEGEN_TOL*max(1.0d0,abs(Wblk(gstart)))) exit
            gend = gend + 1
         enddo
         group_elec = min(remaining_elec, real(2*(gend-gstart+1),8))
         if (closed_pairing) then
            occ_a_per_orbital = 0.5d0*group_elec/real(gend-gstart+1,8)
            occ_b_per_orbital = occ_a_per_orbital
         else
            occ_a_per_orbital = min(1.0d0, group_elec/real(gend-gstart+1,8))
            occ_b_per_orbital = max(0.0d0, group_elec/real(gend-gstart+1,8) - 1.0d0)
         endif
         do k = gstart,gend
            do p = 1,cnt
               do q = 1,cnt
                  Pa_new(p,q) = Pa_new(p,q) + occ_a_per_orbital*Ablk(p,k)*Ablk(q,k)
                  Pb_new(p,q) = Pb_new(p,q) + occ_b_per_orbital*Ablk(p,k)*Ablk(q,k)
               enddo
            enddo
         enddo
         remaining_elec = remaining_elec - group_elec
         gstart = gend + 1
      enddo
      endif

      rms = sqrt(sum((Pa_new-Pa_atom)**2 + (Pb_new-Pb_atom)**2)/real(cnt*cnt,8))
      Pa_atom = Pa_new
      Pb_atom = Pb_new

      if (iter .gt. 1 .and. rms .lt. 1.0d-8) exit

      Ptot = Pa_atom + Pb_atom
      do p = 1,cnt
         do q = 1,cnt
            Jmat(p,q) = sum(Ptot(:,:)*eri(p,q,:,:))
            Ka(p,q)   = sum(Pa_atom(:,:)*eri(p,:,q,:))
            Kb(p,q)   = sum(Pb_atom(:,:)*eri(p,:,q,:))
         enddo
      enddo
      Feff = Hcore_l + Jmat - 0.5d0*(Ka+Kb)
   enddo
end block

ok = .true.

end subroutine compute_atomic_density

subroutine puream_atomic_fill(cnt, numShell, shell_off, bas_l, Z_atom, ecpCoreElec, &
                               closed_pairing, Feff, S_l, Pa_new, Pb_new, INFO)
implicit none
integer,intent(in) :: cnt, numShell, Z_atom, ecpCoreElec
integer,intent(in) :: shell_off(0:numShell-1)
integer,intent(in) :: bas_l(8,numShell)
logical,intent(in) :: closed_pairing
real(8),intent(in) :: Feff(cnt,cnt), S_l(cnt,cnt)
real(8),intent(inout) :: Pa_new(cnt,cnt), Pb_new(cnt,cnt)
integer,intent(out) :: INFO

integer :: lval, degen, nsh_l, j0, k2, p, q, mm, rr, ndocc, ndocc_eff
integer,allocatable :: lshell(:)
real(8),allocatable :: Fl(:,:), Sl2(:,:), Wl(:), WORKl(:)
integer :: LWORKl, INFOl
real(8) :: frac, occ_tot, occ_a, occ_b, acc, accS

INFO = 0
do lval = 0,6
   degen = 2*lval+1
   nsh_l = 0
   do j0 = 1,numShell
      if (bas_l(2,j0) .eq. lval) nsh_l = nsh_l + 1
   enddo
   if (nsh_l .eq. 0) cycle

   allocate(lshell(nsh_l))
   k2 = 0
   do j0 = 1,numShell
      if (bas_l(2,j0) .eq. lval) then
         k2 = k2 + 1
         lshell(k2) = j0 - 1
      endif
   enddo

   allocate(Fl(nsh_l,nsh_l), Sl2(nsh_l,nsh_l), Wl(nsh_l))
   do p = 1,nsh_l
      do q = 1,nsh_l
         acc = 0.0d0; accS = 0.0d0
         do mm = 1,degen
            acc  = acc  + Feff(shell_off(lshell(p))+mm, shell_off(lshell(q))+mm)
            accS = accS + S_l (shell_off(lshell(p))+mm, shell_off(lshell(q))+mm)
         enddo
         Fl(p,q)  = acc /real(degen,8)
         Sl2(p,q) = accS/real(degen,8)
      enddo
   enddo
   LWORKl = 1 + 6*nsh_l + nsh_l**2
   allocate(WORKl(LWORKl))
   call DSYGV(1,'V','U',nsh_l,Fl,nsh_l,Sl2,nsh_l,Wl,WORKl,LWORKl,INFOl)
   if (INFOl .ne. 0) then
      INFO = INFOl
      deallocate(lshell,Fl,Sl2,Wl,WORKl)
      return
   endif

   call atomic_l_occ(Z_atom, ecpCoreElec, lval, ndocc, frac)
   ndocc_eff = min(ndocc, nsh_l)

   do rr = 1,nsh_l
      if (rr .le. ndocc_eff) then
         occ_tot = 2.0d0
      else if (rr .eq. ndocc_eff+1 .and. frac .gt. 0.0d0) then
         occ_tot = frac
      else
         occ_tot = 0.0d0
      endif
      if (occ_tot .le. 0.0d0) cycle
      if (closed_pairing) then
         occ_a = 0.5d0*occ_tot
         occ_b = occ_a
      else
         occ_a = min(1.0d0, occ_tot)
         occ_b = max(0.0d0, occ_tot - 1.0d0)
      endif
      do p = 1,nsh_l
         do q = 1,nsh_l
            acc  = occ_a*Fl(p,rr)*Fl(q,rr)
            accS = occ_b*Fl(p,rr)*Fl(q,rr)
            do mm = 1,degen
               Pa_new(shell_off(lshell(p))+mm, shell_off(lshell(q))+mm) = &
                  Pa_new(shell_off(lshell(p))+mm, shell_off(lshell(q))+mm) + acc
               Pb_new(shell_off(lshell(p))+mm, shell_off(lshell(q))+mm) = &
                  Pb_new(shell_off(lshell(p))+mm, shell_off(lshell(q))+mm) + accS
            enddo
         enddo
      enddo
   enddo
   deallocate(lshell,Fl,Sl2,Wl,WORKl)
enddo
end subroutine puream_atomic_fill

subroutine atomic_l_occ(Z, ecpCoreElec, lval, ndocc, frac)
implicit none
integer,intent(in) :: Z, ecpCoreElec, lval
integer,intent(out) :: ndocc
real(8),intent(out) :: frac
integer,parameter :: NZ = 118
integer :: cfg(4,NZ)
integer :: ne, nd

ndocc = 0
frac = 0.0d0
if (lval .ge. 4 .or. Z .lt. 1 .or. Z .gt. NZ) return

cfg = reshape((/ &
     & 1, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 0,  &
     & 4, 1, 0, 0, 4, 2, 0, 0, 4, 3, 0, 0, 4, 4, 0, 0,  &
     & 4, 5, 0, 0, 4, 6, 0, 0, 5, 6, 0, 0, 6, 6, 0, 0,  &
     & 6, 7, 0, 0, 6, 8, 0, 0, 6, 9, 0, 0, 6, 10, 0, 0,  &
     & 6, 11, 0, 0, 6, 12, 0, 0, 7, 12, 0, 0, 8, 12, 0, 0,  &
     & 8, 13, 0, 0, 8, 12, 2, 0, 8, 12, 3, 0, 8, 12, 4, 0,  &
     & 6, 12, 7, 0, 6, 12, 8, 0, 6, 12, 9, 0, 6, 12, 10, 0,  &
     & 7, 12, 10, 0, 8, 12, 10, 0, 8, 13, 10, 0, 8, 14, 10, 0,  &
     & 8, 15, 10, 0, 8, 16, 10, 0, 8, 17, 10, 0, 8, 18, 10, 0,  &
     & 9, 18, 10, 0, 10, 18, 10, 0, 10, 19, 10, 0, 10, 18, 12, 0,  &
     & 10, 18, 13, 0, 8, 18, 16, 0, 8, 18, 17, 0, 8, 18, 18, 0,  &
     & 8, 18, 19, 0, 8, 18, 20, 0, 9, 18, 20, 0, 10, 18, 20, 0,  &
     & 10, 19, 20, 0, 10, 20, 20, 0, 10, 21, 20, 0, 10, 22, 20, 0,  &
     & 10, 23, 20, 0, 10, 24, 20, 0, 11, 24, 20, 0, 12, 24, 20, 0,  &
     & 12, 24, 21, 0, 12, 24, 22, 0, 12, 24, 21, 2, 12, 24, 20, 4,  &
     & 12, 24, 20, 5, 12, 24, 20, 6, 12, 24, 20, 7, 11, 24, 20, 9,  &
     & 10, 24, 20, 11, 10, 24, 20, 12, 10, 24, 20, 13, 10, 24, 20, 14,  &
     & 11, 24, 20, 14, 12, 24, 20, 14, 12, 25, 20, 14, 12, 24, 22, 14,  &
     & 12, 24, 23, 14, 10, 24, 26, 14, 10, 24, 27, 14, 10, 24, 28, 14,  &
     & 10, 24, 29, 14, 10, 24, 30, 14, 11, 24, 30, 14, 12, 24, 30, 14,  &
     & 12, 25, 30, 14, 12, 26, 30, 14, 12, 27, 30, 14, 12, 28, 30, 14,  &
     & 12, 29, 30, 14, 12, 30, 30, 14, 13, 30, 30, 14, 14, 30, 30, 14,  &
     & 14, 30, 31, 14, 14, 30, 32, 14, 14, 30, 30, 17, 14, 30, 30, 18,  &
     & 14, 30, 30, 19, 13, 30, 30, 21, 12, 30, 30, 23, 12, 30, 30, 24,  &
     & 12, 30, 30, 25, 12, 30, 30, 26, 12, 30, 30, 27, 12, 30, 30, 28,  &
     & 13, 30, 30, 28, 14, 30, 30, 28, 14, 30, 31, 28, 14, 30, 32, 28,  &
     & 14, 30, 33, 28, 12, 30, 36, 28, 12, 30, 37, 28, 12, 30, 38, 28,  &
     & 12, 30, 39, 28, 12, 30, 40, 28, 13, 30, 40, 28, 14, 30, 40, 28,  &
     & 14, 31, 40, 28, 14, 32, 40, 28, 14, 33, 40, 28, 14, 34, 40, 28,  &
     & 14, 35, 40, 28, 14, 36, 40, 28 /), (/4,NZ/))

ne = cfg(lval+1, Z)
if (ecpCoreElec .eq. 2 .or. ecpCoreElec .eq. 10 .or. ecpCoreElec .eq. 18 .or. &
    ecpCoreElec .eq. 36 .or. ecpCoreElec .eq. 54 .or. ecpCoreElec .eq. 86) &
   ne = ne - cfg(lval+1, ecpCoreElec)
nd = (2*lval+1)*2
ndocc = ne/nd
frac = (real(ne,8)/real(nd,8) - real(ndocc,8))*2.0d0
end subroutine atomic_l_occ

end submodule atomic_guess_impl
