! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

submodule (mod_integrals) core_impl
implicit none
contains

module subroutine integrals_init(info)
use MOL_info
use mod_meminfo, only: get_memory_budget_bytes, get_store_threshold_bytes
use mod_mem_predict, only: predict_core_overhead_bytes, MEM_SAFETY_MARGIN
use mod_basis_files, only: basis_file_for_label
use mod_exchange, only: cosx_enabled
use mod_profile, only: fmt_gb
    implicit none
INCLUDE 'parameter.h'
    integer    :: i,j,k,l,n,m,start,info
    integer    :: ij,kl,ijkl,kl_max
    integer,allocatable  ::  ishls(:),jshls(:)
    integer    :: numShell,numGauss
    character  :: base_temp*30,basis_char*30
    character(len=200) :: bpath
    integer    :: base_ios
    real(8)    :: sum_tmp
    integer :: shls(4)
    real(8),allocatable :: buf1e(:,:), buf2e(:,:,:,:),buf1eV(:,:),buf1eT(:,:)
    integer :: offpoint
    integer :: di, dj, dk, dl
    real(8),external :: CINTgto_norm
    integer,external :: CINTcgto_cart,CINTcgto_spheric
    real(8) :: a,b

    real(8) :: dis
    integer   :: LWORK,LIWORK
    real(8),allocatable   ::  WORK(:),IWORK(:)
    real(8),allocatable   ::  S_e(:),S_temp(:,:),S_X(:,:)
    real(8),parameter     :: X_EVAL_FLOOR = 1.0d-6

    integer(8) :: avail_bytes, need_bytes, overhead_bytes
    logical    :: skip_exact_2e

    real(8),allocatable :: ecp_atom_coor(:,:)
    integer :: n_ecp_atoms

    integer :: offpoint_pc, ipc

    nBases = 0
    do i = 1,Natoms
       base_temp =""
       atoms(i)%nsShell = 0
       atoms(i)%npShell = 0
       atoms(i)%ndShell = 0
       atoms(i)%nfShell = 0
       atoms(i)%ngShell = 0
       bpath = basis_file_for_label(atoms(i)%base)
       if (len_trim(bpath) == 0) then
          print *, 'Error: no data/bases block "', trim(atoms(i)%base), &
                   '" found (atom ', i, ') - check baselabel/basis for a typo.'
          stop 1
       endif
       open (unit = 201, file = trim(bpath))
       do while(atoms(i)%base /= base_temp)
            read(201,*,iostat=base_ios) base_temp
            if (base_ios /= 0) then
               print *, 'Error: no data/bases block "', trim(atoms(i)%base), &
                        '" found (atom ', i, ') - check baselabel/basis for a typo.'
               stop 1
            endif
       enddo
       do while(.True.)
           read(201,*), basis_char,basis_char
           if(basis_char == "END") exit
           if(basis_char == "s") atoms(i)%nsShell = atoms(i)%nsShell+1
           if(basis_char == "p") atoms(i)%npShell = atoms(i)%npShell+1
           if(basis_char == "d") atoms(i)%ndShell = atoms(i)%ndShell+1
           if(basis_char == "f") atoms(i)%nfShell = atoms(i)%nfShell+1
           if(basis_char == "g") atoms(i)%ngShell = atoms(i)%ngShell+1
       enddo
       rewind(201)
       base_temp =""

       do while(atoms(i)%base /= base_temp)
            read(201,*,iostat=base_ios) base_temp
            if (base_ios /= 0) then
               print *, 'Error: no data/bases block "', trim(atoms(i)%base), &
                        '" found (atom ', i, ') - check baselabel/basis for a typo.'
               stop 1
            endif
       enddo

       numShell = atoms(i)%nsShell + atoms(i)%npShell + atoms(i)%ndShell + atoms(i)%nfShell + atoms(i)%ngShell
       nBases = nBases + numShell
       if (engine_puream) then
          atoms(i)%nconts = atoms(i)%nsShell + 3*atoms(i)%npShell + 5*atoms(i)%ndShell + 7*atoms(i)%nfShell + 9*atoms(i)%ngShell
       else
          atoms(i)%nconts = atoms(i)%nsShell + 3*atoms(i)%npShell + 6*atoms(i)%ndShell + 10*atoms(i)%nfShell + 15*atoms(i)%ngShell
       endif
       atoms(i)%nshell = numShell
       allocate(atoms(i)%shell(numShell))
       do j = 1,numShell
          Read(201,*),numGauss,basis_char
          atoms(i)%shell(j)%nGauss = numGauss
          if (basis_char == "s") atoms(i)%shell(j)%angMoment = 0
          if (basis_char == "p") atoms(i)%shell(j)%angMoment = 1
          if (basis_char == "d") atoms(i)%shell(j)%angMoment = 2
          if (basis_char == "f") atoms(i)%shell(j)%angMoment = 3
          if (basis_char == "g") atoms(i)%shell(j)%angMoment = 4
          allocate(atoms(i)%shell(j)%exponents(numGauss))
          allocate(atoms(i)%shell(j)%contrCoeff(numGauss))
          do k = 1,numGauss
             Read(201,*),atoms(i)%shell(j)%exponents(k),atoms(i)%shell(j)%contrCoeff(k)
          enddo
       enddo
       close(201)
    enddo

    allocate(atm(6,nAtoms))
    allocate(bas(8,nBases))
    block
       integer :: env_need, ia, ish
       env_need = 20 + 3*nAtoms + 1
       do ia = 1,nAtoms
          do ish = 1,atoms(ia)%nshell
             env_need = env_need + 2*atoms(ia)%shell(ish)%nGauss
          enddo
       enddo
       allocate(env(env_need))
    end block
    atm =0
    bas =0
    env =0
    block
       character(len=64) :: expcut_env
       integer :: expcut_stat
       real(8) :: expcut_val
       expcut_val = 40.0d0
       call get_environment_variable('ENGINE_LIBCINT_EXPCUTOFF', expcut_env, status=expcut_stat)
       if (expcut_stat .eq. 0 .and. len_trim(expcut_env) .gt. 0) then
          read(expcut_env,*,iostat=expcut_stat) expcut_val
          if (expcut_stat .ne. 0) expcut_val = 40.0d0
       endif
       env(1) = expcut_val
    end block
    offpoint = 20
    do i = 1,nAtoms
       atm(1,i) = atoms(i)%charge
       atm(2,i) = offpoint
       env(offpoint+1) = atoms(i)%coor(1)*ans2bohr
       env(offpoint+2) = atoms(i)%coor(2)*ans2bohr
       env(offpoint+3) = atoms(i)%coor(3)*ans2bohr
       offpoint = offpoint + 3
    enddo
    offpoint = offpoint +1
    start = 0
    do i =1,nAtoms
        numShell = atoms(i)%nsShell + atoms(i)%npShell + atoms(i)%ndShell + atoms(i)%nfShell + atoms(i)%ngShell
        do j = 1,numShell
           bas(1,j+start) = i-1
           bas(2,j+start) = atoms(i)%shell(j)%angMoment
           bas(3,j+start) = atoms(i)%shell(j)%nGauss
           bas(4,j+start) = 1
           bas(6,j+start) = offpoint
           do k = 1,atoms(i)%shell(j)%nGauss
               env(offpoint+k) = atoms(i)%shell(j)%exponents(k)
           enddo
           offpoint = offpoint + atoms(i)%shell(j)%nGauss

           bas(7,j+start) = offpoint
           do k = 1,atoms(i)%shell(j)%nGauss
                env(offpoint+k) = atoms(i)%shell(j)%contrCoeff(k)*CINTgto_norm( &
                                  atoms(i)%shell(j)%angMoment, atoms(i)%shell(j)%exponents(k))
           enddo
           offpoint = offpoint + atoms(i)%shell(j)%nGauss

        enddo
        start  = start +numShell
    enddo

block
   integer,parameter :: ECP_MAX_CHANNELS = 6, ECP_MAX_TERMS = 20
   integer :: ecp_l(nAtoms,ECP_MAX_CHANNELS)
   integer :: ecp_nterms(nAtoms,ECP_MAX_CHANNELS)
   integer :: ecp_power(nAtoms,ECP_MAX_CHANNELS,ECP_MAX_TERMS)
   real(8) :: ecp_exp(nAtoms,ECP_MAX_CHANNELS,ECP_MAX_TERMS)
   real(8) :: ecp_coef(nAtoms,ECP_MAX_CHANNELS,ECP_MAX_TERMS)
   integer :: ecp_nchan(nAtoms)
   logical :: has_ecp(nAtoms)
   character(len=2) :: letter
   character(len=40) :: ecp_label, tok
   integer :: ios, nt, ch, t, run_start, shidx
   logical :: ecp_explicit(nAtoms)

   has_ecp = .false.
   ecp_nchan = 0
   do i = 1,nAtoms
      ecp_explicit(i) = (trim(atoms(i)%ecpbase) /= trim(atoms(i)%base))
      ecp_label = trim(atoms(i)%ecpbase)//".ecp"
      bpath = basis_file_for_label(ecp_label)
      if (len_trim(bpath) == 0) then
         if (ecp_explicit(i)) &
            print '(A,I0,A,A,A)', 'Note: atom ', i, ' - no data/bases block "', &
                                   trim(ecp_label), '" found; treated as all-electron.'
         cycle
      endif
      open(unit=202, file=trim(bpath))
      do
         read(202,*,iostat=ios) tok
         if (ios /= 0) exit
         if (trim(tok) == trim(ecp_label)) then
            has_ecp(i) = .true.
            exit
         endif
      enddo
      if (ecp_explicit(i) .and. .not. has_ecp(i)) &
         print '(A,I0,A,A,A)', 'Note: atom ', i, ' - no data/bases block "', &
                                trim(ecp_label), '" found; treated as all-electron.'
      if (.not. has_ecp(i)) then
         close(202)
         cycle
      endif
      read(202,*) atoms(i)%ecpCoreElec
      ch = 0
      do
         read(202,*,iostat=ios) letter, nt
         if (ios /= 0) exit
         ch = ch + 1
         ecp_nchan(i) = ch
         select case (trim(letter))
         case ("ul"); ecp_l(i,ch) = -1
         case ("s");  ecp_l(i,ch) = 0
         case ("p");  ecp_l(i,ch) = 1
         case ("d");  ecp_l(i,ch) = 2
         case ("f");  ecp_l(i,ch) = 3
         case ("g");  ecp_l(i,ch) = 4
         end select
         ecp_nterms(i,ch) = nt
         do t = 1,nt
            read(202,*) ecp_power(i,ch,t), ecp_exp(i,ch,t), ecp_coef(i,ch,t)
         enddo
      enddo
      close(202)
   enddo

   if (any(ecp_explicit) .and. .not. any(has_ecp .and. ecp_explicit)) then
      print *, 'WARNING: an ECP name was requested (ecplabel or a per-atom ecp ' // &
               'override) but matched NO atom in data/bases - likely a typo. ' // &
               'Every requested atom is running all-electron (no ECP applied).'
   endif

   n_ecp_atoms = count(has_ecp)
   allocate(ecp_atom_coor(3,max(n_ecp_atoms,1)))
   n_ecp_atoms = 0
   do i = 1,nAtoms
      if (.not. has_ecp(i)) cycle
      n_ecp_atoms = n_ecp_atoms + 1
      ecp_atom_coor(:,n_ecp_atoms) = atoms(i)%coor*ans2bohr
   enddo

   necpbas = 0
   do i = 1,nAtoms
      if (.not. has_ecp(i)) cycle
      do ch = 1,ecp_nchan(i)
         nt = ecp_nterms(i,ch)
         necpbas = necpbas + 1
         do t = 2,nt
            if (ecp_power(i,ch,t) /= ecp_power(i,ch,t-1)) necpbas = necpbas + 1
         enddo
      enddo
   enddo

   allocate(basECP(8, nBases+necpbas))
   allocate(envECP(size(env) + 2*ECP_MAX_TERMS*max(necpbas,1)))
   basECP = 0
   envECP = 0.0d0
   basECP(1:8,1:nBases) = bas(1:8,1:nBases)
   envECP(1:size(env)) = env(1:size(env))
   offpoint = size(env)
   shidx = nBases

   do i = 1,nAtoms
      if (.not. has_ecp(i)) cycle
      do ch = 1,ecp_nchan(i)
         nt = ecp_nterms(i,ch)
         t = 1
         do while (t <= nt)
            run_start = t
            do while (t < nt)
               if (ecp_power(i,ch,t+1) /= ecp_power(i,ch,run_start)) exit
               t = t + 1
            enddo
            shidx = shidx + 1
            basECP(1,shidx) = i-1
            basECP(2,shidx) = ecp_l(i,ch)
            basECP(3,shidx) = t - run_start + 1
            basECP(4,shidx) = ecp_power(i,ch,run_start)
            basECP(5,shidx) = 0
            basECP(6,shidx) = offpoint
            do k = run_start,t
               envECP(offpoint+1+(k-run_start)) = ecp_exp(i,ch,k)
            enddo
            offpoint = offpoint + (t-run_start+1)
            basECP(7,shidx) = offpoint
            do k = run_start,t
               envECP(offpoint+1+(k-run_start)) = ecp_coef(i,ch,k)
            enddo
            offpoint = offpoint + (t-run_start+1)
            t = t + 1
         enddo
      enddo
   enddo

   envECP(18) = 0.0d0
   envECP(19) = dble(nBases)
   envECP(20) = dble(necpbas)
end block

    nAtoms_nuc = nAtoms + nPointCharges
    allocate(atm_nuc(6,nAtoms_nuc))
    atm_nuc = 0
    atm_nuc(:,1:nAtoms) = atm(:,1:nAtoms)
    do i = 1,nAtoms
       atm_nuc(1,i) = atoms(i)%charge - atoms(i)%ecpCoreElec
    enddo
    allocate(env_nuc(size(env) + 4*nPointCharges))
    env_nuc = 0
    env_nuc(1:size(env)) = env(1:size(env))
    offpoint_pc = size(env)
    do ipc = 1,nPointCharges
       atm_nuc(2,nAtoms+ipc) = offpoint_pc
       env_nuc(offpoint_pc+1) = pointcharge_coor(1,ipc)*ans2bohr
       env_nuc(offpoint_pc+2) = pointcharge_coor(2,ipc)*ans2bohr
       env_nuc(offpoint_pc+3) = pointcharge_coor(3,ipc)*ans2bohr
       offpoint_pc = offpoint_pc + 3
       atm_nuc(3,nAtoms+ipc) = 3
       atm_nuc(5,nAtoms+ipc) = offpoint_pc
       env_nuc(offpoint_pc+1) = pointcharge_q(ipc)
       offpoint_pc = offpoint_pc + 1
    enddo

 do i = 1,nAtoms
    nConts = atoms(i)%nconts+nConts
 enddo

 do i = 1,nAtoms
    nContssph = nContssph + atoms(i)%nsShell + 3*atoms(i)%npShell &
              + 5*atoms(i)%ndShell + 7*atoms(i)%nfShell + 9*atoms(i)%ngShell
 enddo

nRec = INDEX_2E(nConts-1,nConts-1,nConts-1,nConts-1)+1

print '("nConts = ", I0, "   (# contracted basis functions)   nRec = ", I0, "   (# unique 2e-integral records)")', nConts, nRec

overhead_bytes = int(real(predict_core_overhead_bytes(nConts, nAtoms),8) * MEM_SAFETY_MARGIN, 8)
need_bytes = int(nRec,8) * 8_8
avail_bytes = get_memory_budget_bytes()
direct_mode = ((need_bytes + overhead_bytes) > get_store_threshold_bytes())
block
   character(len=8) :: fdir
   fdir = ""
   call get_environment_variable("ENGINE_FORCE_TWOEI_DIRECT", fdir)
   if (trim(fdir) .eq. "1") direct_mode = .true.
   if (trim(fdir) .eq. "0") direct_mode = .false.
end block
if (.not. engine_use_df) then
   print '(A)', "  TWOEI would need="//trim(fmt_gb(real(need_bytes,8)/1024**3))// &
           "GB  (plus other resident arrays)="//trim(fmt_gb(real(overhead_bytes,8)/1024**3))// &
           "GB  available="//trim(fmt_gb(real(avail_bytes,8)/1024**3))// &
           "GB  STORE threshold="//trim(fmt_gb(real(get_store_threshold_bytes(),8)/1024**3))//"GB"
   print '(A)', "  Build mode: "// &
           trim(merge("DIRECT (recomputed every SCF cycle)","STORE (held in memory)             ",direct_mode))
   call flush(6)
endif
print *, '[SYSTEMEND]'
print *

allocate(NorVEC(nConts))
allocate(NorVECsph(nBases))
NorVEC = 0
do i = 0,nBases-1
  shls(1)=i
  shls(2)=i
  if (.not. engine_puream) then
     di = CINTcgto_cart(i, bas)
     allocate (buf1e(di,di))
     call cint1e_ovlp_cart(buf1e, shls , atm, nAtoms,bas,nBases,env,0_8)
     call Normal(i,di,bas,nBases,buf1e,NorVEC,nConts)
     deallocate (buf1e)
  endif
  di = CINTcgto_spheric(i, bas)
  allocate (buf1e(di,di))
  call cint1e_ovlp_sph(buf1e,shls,atm,nAtoms,bas,nBases,env,0_8)
  NorVECsph(1+i) = dsqrt(1.0/buf1e(1,1))
  deallocate (buf1e)
enddo
if (engine_puream) NorVEC = 1.0d0
norvec_is_identity = all(NorVEC .eq. 1.0d0)

k = 1
do i =1,natoms
   do j = 1,atoms(i)%nshell
      atoms(i)%shell(j)%contrCoeff(:) =atoms(i)%shell(j)%contrCoeff(:)&
                                       *NorVECsph(k)
      k = k +1
   enddo
enddo

k = 1
start = 0
do i = 1,natoms
   do j = 1,atoms(i)%nshell
      do l = 1,atoms(i)%shell(j)%nGauss
         env(bas(7,j+start)+l) = atoms(i)%shell(j)%contrCoeff(l)*CINTgto_norm( &
                                  atoms(i)%shell(j)%angMoment, atoms(i)%shell(j)%exponents(l))
         if (allocated(envECP))  envECP(bas(7,j+start)+l)  = env(bas(7,j+start)+l)
         if (allocated(env_nuc)) env_nuc(bas(7,j+start)+l) = env(bas(7,j+start)+l)
      enddo
      k = k +1
   enddo
   start = start + atoms(i)%nshell
enddo

block
   integer :: ia,ish,ip,aa,bb,cc,idim,Lt,Lx1,Ly1,Lz1
   integer,external :: factorial2
   real(8) :: NorA
   do ia = 1,natoms
      do ish = 1,atoms(ia)%nshell
         Lt = atoms(ia)%shell(ish)%angMoment
         allocate(atoms(ia)%shell(ish)%cnVal(atoms(ia)%shell(ish)%nGauss, &
                  (Lt+1)*(Lt+2)/2))
         idim = 0
         do aa = Lt,0,-1
            do bb = Lt-aa,0,-1
               cc = Lt-aa-bb
               idim = idim+1
               Lx1 = factorial2(aa)
               Ly1 = factorial2(bb)
               Lz1 = factorial2(cc)
               do ip = 1,atoms(ia)%shell(ish)%nGauss
                  NorA = (2**(2*Lt+1.5)*atoms(ia)%shell(ish)%exponents(ip)**(Lt+1.5))**0.5 &
                         * (Lx1*Ly1*Lz1*PI**(1.5))**(-0.5)
                  atoms(ia)%shell(ish)%cnVal(ip,idim) = atoms(ia)%shell(ish)%contrCoeff(ip)*NorA
               enddo
            enddo
         enddo
      enddo
   enddo
end block

allocate( S (nConts,nConts) )
allocate( Hcore (nConts,nConts))
S = 0
Hcore = 0
n = 0
do i = 0,nBases-1
   n = max(n, cgto_engine(i,bas))
enddo
allocate(buf1e(n,n))
allocate(buf1eT(n,n))
allocate(buf1eV(n,n))
block
   real(8),allocatable :: buf1eECP(:,:)
   integer :: dims_arr(2)
   integer :: iatm_of, jatm_of, iecp
   real(8) :: r_i(3), r_j(3), r2_min
   logical :: ecp_pair_relevant
   real(8),parameter :: ECP_SKIP_R2 = 144.0d0
   if (necpbas > 0) allocate(buf1eECP(n,n))
   do i = 0,nBases-1
      do j = i,nBases-1
         shls(1)=i
         shls(2)=j
         di = cgto_engine(i, bas)
         dj = cgto_engine(j, bas)
         call ovlp1e_engine(buf1e(1:di,1:dj), shls , atm, nAtoms,bas,nBases,env,0_8)
         call store1e(shls,di,dj,bas,nBases,buf1e(1:di,1:dj),NorVEC,nConts,S)
         call nuc1e_engine(buf1eV(1:di,1:dj), shls , atm_nuc, nAtoms_nuc,bas,nBases,env_nuc,0_8)
         call kin1e_engine(buf1eT(1:di,1:dj), shls , atm, nAtoms,bas,nBases,env,0_8)
         if (necpbas > 0) then
            iatm_of = bas(1,i+1) + 1
            jatm_of = bas(1,j+1) + 1
            r_i = atoms(iatm_of)%coor*ans2bohr
            r_j = atoms(jatm_of)%coor*ans2bohr
            ecp_pair_relevant = .false.
            do iecp = 1,n_ecp_atoms
               r2_min = min(sum((r_i-ecp_atom_coor(:,iecp))**2), &
                             sum((r_j-ecp_atom_coor(:,iecp))**2))
               if (r2_min .le. ECP_SKIP_R2) then
                  ecp_pair_relevant = .true.
                  exit
               endif
            enddo
            if (.not. ecp_pair_relevant) then
               buf1eECP(1:di,1:dj) = 0.0d0
            else
               dims_arr(1) = di
               dims_arr(2) = dj
               call ecp1e_engine(buf1eECP(1:di,1:dj), dims_arr, shls, atm, &
                                  nAtoms, basECP, nBases+necpbas, envECP, 0_8, 0_8)
            endif
            call store1e(shls,di,dj,bas,nBases, &
                         buf1eT(1:di,1:dj)+buf1eV(1:di,1:dj)+buf1eECP(1:di,1:dj), &
                         NorVEC,nConts,Hcore)
         else
            call store1e(shls,di,dj,bas,nBases,buf1eT(1:di,1:dj)+buf1eV(1:di,1:dj),NorVEC,nConts,Hcore)
         endif
      enddo
   enddo
   if (allocated(buf1eECP)) deallocate(buf1eECP)
end block
deallocate(buf1e)
deallocate(buf1eT)
deallocate(buf1eV)

    LWORK = 1 + 6*nconts + nconts**2

    allocate (S_e(nconts))
    allocate (S_temp(nconts,nconts))
    allocate (WORK(LWORK))
    allocate (S_X(nconts,nconts))
    allocate (X(nconts,nconts))
    S_temp = S
    S_X = 0
    call DSYEV('V','U',nconts,S_temp,nconts,S_e,WORK,LWORK,INFO)
    do i = 1,nconts
        if (S_e(i) .gt. X_EVAL_FLOOR*S_e(nconts)) then
            S_X(i,i) = S_e(i)**(-(1.0)/2)
        else
            print *, "WARNING: near-linear-dependent basis function ", &
                     "detected, flooring overlap eigenvalue ", i, " = ", S_e(i)
            S_X(i,i) = 0.0d0
        endif
    enddo
    X = MATMUL(S_temp,MATMUL(S_X,TRANSPOSE(S_temp)))
    deallocate (S_e)
    deallocate (S_temp)
    deallocate (S_X)
    deallocate (WORK)

call compute_schwarz_bounds()

skip_exact_2e = (engine_use_df_j .and. (engine_use_df_k .or. cosx_enabled)) &
                .or. engine_estimate_only

if ((.not. direct_mode) .and. (.not. skip_exact_2e)) then
allocate(TWOEI(nRec))
TWOEI = 0

allocate(ishls(0:nBases*nBases-1))
allocate(jshls(0:nBases*nBases-1))

ij = 0

do i = 0,nBases-1
    do j=0,i
       ishls(ij) = i
       jshls(ij) = j
       ij = ij +1
    enddo
enddo

    n = 0
    do i = 0,nBases-1
       n = max(n, cgto_engine(i,bas))
    enddo
!$omp parallel private(i,j,di,dj,kl,k,l,dk,dl,shls,buf2e)
    allocate(buf2e(n,n,n,n))
!$omp do schedule(dynamic)
do ij = 0, nBases*(nBases+1)/2-1
    i = ishls(ij)
    j = jshls(ij)
    di = cgto_engine(i, bas)
    dj = cgto_engine(j, bas)
    do kl = 0, ij
             k = ishls(kl)
             l = jshls(kl)
             if (schwarz_bound(i,j)*schwarz_bound(k,l) .lt. SCHWARZ_CUTOFF) cycle
             dk = cgto_engine(k, bas)
             dl = cgto_engine(l, bas)
             shls(1) = i
             shls(2) = j
             shls(3) = k
             shls(4) = l
             call twoe_engine(buf2e(1:di,1:dj,1:dk,1:dl), shls, atm, nAtoms,bas,nBases,env,int2e_opt)
             call store2e(shls,bas,buf2e(1:di,1:dj,1:dk,1:dl),TWOEI,di,dj,dk,dl,nConts,nRec,nBases,NorVEC)
    enddo
enddo
!$omp end do
    deallocate(buf2e)
!$omp end parallel

block
   real(8),allocatable :: buf2e_chk(:,:,:,:)
   integer :: ij_chk,kl_chk,i_chk,j_chk,k_chk,l_chk,di_chk,dj_chk,dk_chk,dl_chk
   integer :: shls_chk(4)
   integer :: p_chk,q_chk,r_chk,s_chk,e1_chk,e2_chk,e3_chk,e4_chk
   real(8) :: stored_val,fresh_val,maxdiff,n1,n2,n3,n4
   integer :: nchecked
   maxdiff = 0.0d0
   nchecked = 0
   allocate(buf2e_chk(n,n,n,n))
   do ij_chk = 0, min(200_8, nBases*(nBases+1)/2-1)
      i_chk = ishls(ij_chk); j_chk = jshls(ij_chk)
      di_chk = cgto_engine(i_chk, bas); dj_chk = cgto_engine(j_chk, bas)
      do kl_chk = 0, ij_chk
         k_chk = ishls(kl_chk); l_chk = jshls(kl_chk)
         if (schwarz_bound(i_chk,j_chk)*schwarz_bound(k_chk,l_chk) .lt. SCHWARZ_CUTOFF) cycle
         dk_chk = cgto_engine(k_chk, bas); dl_chk = cgto_engine(l_chk, bas)
         shls_chk = (/i_chk,j_chk,k_chk,l_chk/)
         call twoe_engine(buf2e_chk(1:di_chk,1:dj_chk,1:dk_chk,1:dl_chk), shls_chk, atm, nAtoms,bas,nBases,env,int2e_opt)
         e1_chk = 0
         do p_chk = 0,i_chk-1
            e1_chk = e1_chk + cgto_engine(p_chk,bas)
         enddo
         e2_chk = 0
         do p_chk = 0,j_chk-1
            e2_chk = e2_chk + cgto_engine(p_chk,bas)
         enddo
         e3_chk = 0
         do p_chk = 0,k_chk-1
            e3_chk = e3_chk + cgto_engine(p_chk,bas)
         enddo
         e4_chk = 0
         do p_chk = 0,l_chk-1
            e4_chk = e4_chk + cgto_engine(p_chk,bas)
         enddo
         do p_chk = 1,di_chk
            n1 = NorVEC(e1_chk+p_chk)
            do q_chk = 1,dj_chk
               n2 = NorVEC(e2_chk+q_chk)
               do r_chk = 1,dk_chk
                  n3 = NorVEC(e3_chk+r_chk)
                  do s_chk = 1,dl_chk
                     n4 = NorVEC(e4_chk+s_chk)
                     fresh_val = buf2e_chk(p_chk,q_chk,r_chk,s_chk)*n1*n2*n3*n4
                     stored_val = eri_get(e1_chk+p_chk,e2_chk+q_chk,e3_chk+r_chk,e4_chk+s_chk)
                     maxdiff = max(maxdiff, abs(fresh_val-stored_val))
                     nchecked = nchecked + 1
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   deallocate(buf2e_chk)
   print *, 'DEBUG integral-recompute-vs-STORE check: n=',nchecked,' maxdiff=',maxdiff
end block

deallocate(ishls)
deallocate(jshls)
endif

     dis = 0
     E_rep = 0
     do i = 1,natoms
        do j = i+1,natoms
          dis = (atoms(i)%coor(1)-atoms(j)%coor(1))**2 +dis
          dis = (atoms(i)%coor(2)-atoms(j)%coor(2))**2 +dis
          dis = (atoms(i)%coor(3)-atoms(j)%coor(3))**2 +dis
          dis = dis**0.5*ans2bohr
          E_rep = (atoms(i)%charge-atoms(i)%ecpCoreElec) &
                * (atoms(j)%charge-atoms(j)%ecpCoreElec) /dis + E_rep
          dis = 0
        enddo
     enddo

     do i = 1,natoms
        do ipc = 1,nPointCharges
           dis = (atoms(i)%coor(1)-pointcharge_coor(1,ipc))**2 &
               + (atoms(i)%coor(2)-pointcharge_coor(2,ipc))**2 &
               + (atoms(i)%coor(3)-pointcharge_coor(3,ipc))**2
           dis = dis**0.5*ans2bohr
           E_rep = (atoms(i)%charge-atoms(i)%ecpCoreElec) * pointcharge_q(ipc) /dis + E_rep
        enddo
     enddo

end subroutine integrals_init

module subroutine ensure_int2e_optimizer()
    implicit none
    if (int2e_opt .eq. 0_8) then
       call twoe_optimizer_engine(int2e_opt, atm, size(atm,2), bas, nBases, env)
    endif
end subroutine ensure_int2e_optimizer

module subroutine ensure_int2e_optimizer_lr()
    use mod_exchange, only: RS_omega
    implicit none
    if (int2e_opt_lr .eq. 0_8) then
       env(9) = RS_omega
       call twoe_optimizer_engine(int2e_opt_lr, atm, size(atm,2), bas, nBases, env)
       env(9) = 0.0d0
    endif
end subroutine ensure_int2e_optimizer_lr

module subroutine integrals_set_accuracy(prms)
    implicit none
    real(8),intent(in) :: prms
    scf_prms_now = prms
    density_screen_cutoff = min(DENSITY_SCREEN_MAX, &
         max(DENSITY_SCREEN_INCR_FLOOR, prms*DENSITY_SCREEN_SCALE))
    density_screen_cutoff_incr = min(DENSITY_SCREEN_MAX, &
         max(DENSITY_SCREEN_INCR_FLOOR, prms*DENSITY_SCREEN_INCR_SCALE))
    density_screen_cutoff_incr_j = min(DENSITY_SCREEN_MAX, &
         max(DENSITY_SCREEN_INCR_FLOOR, prms*DENSITY_SCREEN_INCR_SCALE_J))
    block
      character(len=32) :: absenv
      real(8) :: scl
      absenv = ""
      call get_environment_variable("ENGINE_FULL_SCREEN_SCALE", absenv)
      if (len_trim(absenv) .gt. 0) then
         read(absenv,*) scl
         density_screen_cutoff = min(DENSITY_SCREEN_MAX, &
              max(DENSITY_SCREEN_INCR_FLOOR, prms*scl))
      endif
      absenv = ""
      call get_environment_variable("ENGINE_INCR_SCREEN_SCALE", absenv)
      if (len_trim(absenv) .gt. 0) then
         read(absenv,*) scl
         density_screen_cutoff_incr = min(DENSITY_SCREEN_MAX, &
              max(DENSITY_SCREEN_INCR_FLOOR, prms*scl))
      endif
      absenv = ""
      call get_environment_variable("ENGINE_INCR_SCREEN_SCALE_J", absenv)
      if (len_trim(absenv) .gt. 0) then
         read(absenv,*) scl
         density_screen_cutoff_incr_j = min(DENSITY_SCREEN_MAX, &
              max(DENSITY_SCREEN_INCR_FLOOR, prms*scl))
      endif
      absenv = ""
      call get_environment_variable("ENGINE_INCR_SCREEN_ABS", absenv)
      if (len_trim(absenv) .gt. 0) then
         read(absenv,*) density_screen_cutoff_incr
         density_screen_cutoff_incr_j = density_screen_cutoff_incr
      endif
    end block
end subroutine integrals_set_accuracy

module subroutine compute_schwarz_bounds()
    implicit none
    integer :: si,sj,di,dj,shls(4),p,q,offi,offj,ii,num
    real(8),allocatable :: buf2e(:,:,:,:)
    real(8) :: bound

    call ensure_int2e_optimizer()
    allocate(schwarz_bound(0:nBases-1,0:nBases-1))
    schwarz_bound = 0

    num = 0
    do si = 0,nBases-1
       num = num + cgto_engine(si,bas)
    enddo
    allocate(ao_shell(num))
    num = 0
    do si = 0,nBases-1
       di = cgto_engine(si,bas)
       do p = 1,di
          ao_shell(num+p) = si
       enddo
       num = num + di
    enddo
    do si = 0,nBases-1
       di = cgto_engine(si,bas)
       num = 0
       do ii = 0,si-1
          num = num + cgto_engine(ii,bas)
       enddo
       offi = num
       do sj = 0,si
          dj = cgto_engine(sj,bas)
          num = 0
          do ii = 0,sj-1
             num = num + cgto_engine(ii,bas)
          enddo
          offj = num
          shls(1) = si; shls(2) = sj; shls(3) = si; shls(4) = sj
          allocate(buf2e(di,dj,di,dj))
          call twoe_engine(buf2e, shls, atm, size(atm,2), bas, nBases, env, int2e_opt)
          bound = 0
          do p = 1,di
             do q = 1,dj
                bound = max(bound, abs(buf2e(p,q,p,q))*NorVEC(offi+p)**2*NorVEC(offj+q)**2)
             enddo
          enddo
          deallocate(buf2e)
          bound = sqrt(bound)
          schwarz_bound(si,sj) = bound
          schwarz_bound(sj,si) = bound
       enddo
    enddo
end subroutine compute_schwarz_bounds

module subroutine build_significant_pairs()
    implicit none
    integer :: si,sj,nshpair
    real(8) :: max_schwarz_bound

    if (sig_pairs_built) return

    max_schwarz_bound = maxval(schwarz_bound)
    nshpair = nBases*(nBases+1)/2

    n_sig_pairs = 0
    do si = 0,nBases-1
       do sj = 0,si
          if (schwarz_bound(si,sj)*max_schwarz_bound .ge. SCHWARZ_CUTOFF) then
             n_sig_pairs = n_sig_pairs + 1
          endif
       enddo
    enddo

    allocate(sig_pair_ishl(0:n_sig_pairs-1),sig_pair_jshl(0:n_sig_pairs-1))
    n_sig_pairs = 0
    do si = 0,nBases-1
       do sj = 0,si
          if (schwarz_bound(si,sj)*max_schwarz_bound .ge. SCHWARZ_CUTOFF) then
             sig_pair_ishl(n_sig_pairs) = si
             sig_pair_jshl(n_sig_pairs) = sj
             n_sig_pairs = n_sig_pairs + 1
          endif
       enddo
    enddo
    print *,"DIAG: build_significant_pairs: kept",n_sig_pairs," of",nshpair," shell pairs"
    sig_pairs_built = .true.
end subroutine build_significant_pairs

module real(8) function eri_get(i,j,k,l) result(val)
    implicit none
    integer,intent(in) :: i,j,k,l
    val = TWOEI(INDEX_2E(i-1,j-1,k-1,l-1)+1)
end function eri_get

module subroutine integrals_finalize()
    implicit none
    if (allocated(atm))       deallocate(atm)
    if (allocated(bas))       deallocate(bas)
    if (allocated(env))       deallocate(env)
    if (allocated(atm_nuc))   deallocate(atm_nuc)
    if (allocated(env_nuc))   deallocate(env_nuc)
    if (allocated(NorVEC))    deallocate(NorVEC)
    if (allocated(NorVECsph)) deallocate(NorVECsph)
    if (allocated(TWOEI))     deallocate(TWOEI)
    if (allocated(schwarz_bound)) deallocate(schwarz_bound)
    if (int2e_opt .ne. 0_8) then
       call cintdel_optimizer(int2e_opt)
       int2e_opt = 0_8
    endif
    if (int2e_opt_lr .ne. 0_8) then
       call cintdel_optimizer(int2e_opt_lr)
       int2e_opt_lr = 0_8
    endif
    if (allocated(ao_shell)) deallocate(ao_shell)
    if (allocated(basECP)) deallocate(basECP)
    if (allocated(envECP)) deallocate(envECP)
    if (allocated(cPK_pq)) deallocate(cPK_pq,cPK_rs,cPK_val)
    coulomb_list_built = .false.
    cPK_n = 0
    coulomb_npairs = 0
    nBases = 0
    nRec = 0
    if (allocated(incr_Ptot_prev)) deallocate(incr_Ptot_prev,incr_J_prev)
    incr_has_prev = .false.
    incr_since_full = 0
    incr_fock_env_checked = .false.
    incr_fock_enabled = .true.
    density_screen_cutoff = SCHWARZ_CUTOFF
    density_screen_cutoff_incr = SCHWARZ_CUTOFF
    density_screen_cutoff_incr_j = SCHWARZ_CUTOFF
    if (allocated(sig_pair_ishl)) deallocate(sig_pair_ishl,sig_pair_jshl)
    sig_pairs_built = .false.
    n_sig_pairs = 0
    if (allocated(cK_i)) deallocate(cK_i,cK_j,cK_k,cK_l,cK_val)
    exchange_list_built = .false.
    cK_n = 0
    if (allocated(basDF)) deallocate(basDF)
    if (allocated(envDF)) deallocate(envDF)
    if (allocated(NorVECAux)) deallocate(NorVECAux)
    if (allocated(dfB_compact)) deallocate(dfB_compact)
    if (allocated(dfB_compact_LR)) deallocate(dfB_compact_LR)
    if (allocated(df_evec)) deallocate(df_evec)
    if (allocated(df_evalinv)) deallocate(df_evalinv)
    if (allocated(aux_shell_bound)) deallocate(aux_shell_bound)
    if (int3c2e_opt .ne. 0_8) then
       call cintdel_optimizer(int3c2e_opt)
       int3c2e_opt = 0_8
    endif
    if (int3c2e_opt_lr .ne. 0_8) then
       call cintdel_optimizer(int3c2e_opt_lr)
       int3c2e_opt_lr = 0_8
    endif
    if (int2c2e_opt .ne. 0_8) then
       call cintdel_optimizer(int2c2e_opt)
       int2c2e_opt = 0_8
    endif
    df_built = .false.
    df_direct_mode = .false.
    nBasesAux = 0
    nContsAux = 0
    if (allocated(dfK)) deallocate(dfK)
    if (allocated(dfK_compact)) deallocate(dfK_compact)
    if (allocated(dfC_full_compact)) deallocate(dfC_full_compact)
    if (allocated(df_pair_row)) deallocate(df_pair_row,df_pair_col)
    if (allocated(df_pair_mirror)) deallocate(df_pair_mirror)
    if (allocated(df_row_start)) deallocate(df_row_start,df_row_count)
    df_total_pairs = 0
    df_max_row_count = 0
    dfK_built = .false.
    df_full_lr_metric_built = .false.
    if (allocated(df_Minvhalf)) deallocate(df_Minvhalf)
    df_Minvhalf_built = .false.
    if (allocated(df_Minv_full)) deallocate(df_Minv_full)
    df_Minv_full_built = .false.
    df_lr_optimizer_built = .false.
    if (allocated(df_triple_ish)) deallocate(df_triple_ish,df_triple_jsh,df_triple_ksh)
    if (allocated(df_shell_ncgto)) deallocate(df_shell_ncgto,df_shell_offset)
    df_n_triples = 0
end subroutine integrals_finalize

module subroutine Normal(iL,di,bas,nBasis,buf1e,Norfac,nConts)
    implicit none
    integer :: iL,di,nConts,nBasis
    real(8) :: buf1e(di,di)
    integer :: bas(8,nbasis)
    real(8) :: Norfac(nConts)

    real(8) :: Frac
    integer :: startnum
    integer :: i,j,k
    startnum = 0
    do i = 0,iL-1
       startnum = startnum + cgto_engine(i, bas)
    enddo

    do i = 1,di
       Frac = dsqrt(1.0/buf1e(i,i))
       Norfac(startnum+i) = Frac
    enddo
end subroutine Normal

module integer(8) function INDEX_2E(i,j,k,l)
    implicit none
    integer  ::   i,j,k,l
    integer(8)  ::   ij,kl

    if (i>j) then
        ij = int(i,8)*(i+1)/2 + j
    else
        ij = int(j,8)*(j+1)/2 + i
    endif

    if (k>l) then
        kl = int(k,8)*(k+1)/2 + l
    else
        kl = int(l,8)*(l+1)/2 + k
    endif

    if (ij > kl) then
        INDEX_2E = (ij+1)*ij/2 + kl
    else
        INDEX_2E = (kl+1)*kl/2 + ij
    endif
end function INDEX_2E

module subroutine  store2e(shls,bas,buf2e,TWOEI,di,dj,dk,dl,nConts,nRec,nBases,Norfac)
implicit none
integer :: di, dj, dk, dl, nConts, nBases
integer(8) :: nRec
integer :: shls(4)
integer :: bas(8,nBases)
real(8) :: buf2e(di,dj,dk,dl)
real(8) :: TWOEI(nRec)
real(8) :: Norfac(nConts)
real(8) :: intVal

integer :: is,js,ks,ls
integer :: i,j,k,l,num
integer :: e1,e2,e3,e4
integer :: ii

num = 0
do ii = 0,shls(1)-1
   num = num + cgto_engine(ii, bas)
enddo
is = num

num = 0
do ii = 0,shls(2)-1
   num = num + cgto_engine(ii, bas)
enddo
js = num

num = 0
do ii = 0,shls(3)-1
   num = num + cgto_engine(ii, bas)
enddo
ks = num

num = 0
do ii = 0,shls(4)-1
   num = num + cgto_engine(ii, bas)
enddo
ls = num

do i = 1,di
   do j = 1,dj
      do k = 1,dk
         do l = 1,dl
              e1 = is+i-1
              e2 = js+j-1
              e3 = ks+k-1
              e4 = ls+l-1
              intVal = buf2e(i,j,k,l)*Norfac(e1+1)*Norfac(e2+1)*Norfac(e3+1)*Norfac(e4+1)
              TWOEI(INDEX_2E(e1,e2,e3,e4)+1) = intVal
         enddo
      enddo
   enddo
enddo
end subroutine store2e

module subroutine  store1e(shls,di,dj,bas,nBases,buf1e,Norfac,nConts,S)
implicit none
integer :: di, dj, nConts,  nBases
integer :: shls(4)
integer :: bas(8,nBases)
real(8) :: S(nConts,nConts)
real(8) :: Norfac(nConts)
real(8) :: buf1e(di,dj)

real(8) :: intVal
integer :: is,js
integer :: i,j,num
integer :: e1,e2
integer :: ii

num = 0
do ii = 0,shls(1)-1
   num = num + cgto_engine(ii, bas)
enddo
is = num

num = 0
do ii = 0,shls(2)-1
   num = num + cgto_engine(ii, bas)
enddo
js = num

do i = 1,di
   do j = 1,dj
        e1 = is+i
        e2 = js+j
        if (e2 .ge. e1) then
          intVal = buf1e(i,j)*Norfac(e1)*Norfac(e2)
          S(e1,e2) = intVal
          S(e2,e1) = intVal
        endif
   enddo
enddo
end subroutine store1e

module subroutine  store1edrv(shls,di,dj,ao_offset,nBases,buf1e,Norfac,nConts,S)
implicit none
integer :: di, dj, nConts,  nBases
integer :: shls(2)
integer :: ao_offset(0:nBases-1)
real(8) :: S(nConts,nConts,3)
real(8) :: Norfac(nConts)
real(8) :: buf1e(di,dj,3)

real(8) :: intVal
integer :: is,js
integer :: i,j,k
integer :: e1,e2

is = ao_offset(shls(1))
js = ao_offset(shls(2))

do i = 1,di
   do j = 1,dj
      do k =1,3
         e1 = is+i
         e2 = js+j
         intVal = buf1e(i,j,k)*Norfac(e1)*Norfac(e2)
         S(e1,e2,k) = intVal
      enddo
   enddo
enddo
end subroutine store1edrv

module subroutine  store2edrv(shls,ao_offset,buf2e,dJi,dKa,dKb,Da,Db,di,dj,dk,dl,nConts,nBases,Norfac)
implicit none
integer :: di, dj, dk, dl, nConts, nBases
integer :: shls(4)
integer :: ao_offset(0:nBases-1)
real(8) :: buf2e(di,dj,dk,dl,3)
real(8) :: dJi(nConts,nConts,3),dKa(nConts,nConts,3),dKb(nConts,nConts,3)
real(8) :: Da(nConts,nConts),Db(nConts,nConts)
real(8) :: Norfac(nConts)
real(8) :: intVal(3)

integer :: is,js,ks,ls
integer :: i,j,k,l
integer :: e1,e2,e3,e4
logical :: ket_distinct
real(8) :: dJi_weight

is = ao_offset(shls(1))
js = ao_offset(shls(2))
ks = ao_offset(shls(3))
ls = ao_offset(shls(4))

ket_distinct = (shls(3) .ne. shls(4))
dJi_weight = merge(2.0d0, 1.0d0, ket_distinct)

if (norvec_is_identity) then
   do i = 1,di
      do j = 1,dj
         do k = 1,dk
            do l = 1,dl
                 e1 = is+i
                 e2 = js+j
                 e3 = ks+k
                 e4 = ls+l
                  intVal = buf2e(i,j,k,l,:)
                  dJi(e1,e2,:) = dJi(e1,e2,:) + dJi_weight*intVal*(Da(e3,e4) +Db(e3,e4))
                  dKa(e1,e3,:) = dKa(e1,e3,:) + intVal*Da(e2,e4)
                  dKb(e1,e3,:) = dKb(e1,e3,:) + intVal*Db(e2,e4)
                  if (ket_distinct) then
                     dKa(e1,e4,:) = dKa(e1,e4,:) + intVal*Da(e2,e3)
                     dKb(e1,e4,:) = dKb(e1,e4,:) + intVal*Db(e2,e3)
                  endif
            enddo
         enddo
      enddo
   enddo
else
do i = 1,di
   do j = 1,dj
      do k = 1,dk
         do l = 1,dl
              e1 = is+i
              e2 = js+j
              e3 = ks+k
              e4 = ls+l
               intVal = buf2e(i,j,k,l,:)*Norfac(e1)*Norfac(e2)*Norfac(e3)*Norfac(e4)
               dJi(e1,e2,:) = dJi(e1,e2,:) + dJi_weight*intVal*(Da(e3,e4) +Db(e3,e4))
               dKa(e1,e3,:) = dKa(e1,e3,:) + intVal*Da(e2,e4)
               dKb(e1,e3,:) = dKb(e1,e3,:) + intVal*Db(e2,e4)
               if (ket_distinct) then
                  dKa(e1,e4,:) = dKa(e1,e4,:) + intVal*Da(e2,e3)
                  dKb(e1,e4,:) = dKb(e1,e4,:) + intVal*Db(e2,e3)
               endif
         enddo
      enddo
   enddo
enddo
endif
end subroutine store2edrv

module subroutine compute_shell_density_bound(nConts, Ptot, dmax_shell)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: dmax_shell(0:nBases-1,0:nBases-1)
    integer :: p,q
    dmax_shell = 0.0d0
    do q = 1,nConts
       do p = 1,nConts
          dmax_shell(ao_shell(p),ao_shell(q)) = &
             max(dmax_shell(ao_shell(p),ao_shell(q)), abs(Ptot(p,q)))
       enddo
    enddo
end subroutine compute_shell_density_bound

end submodule core_impl
