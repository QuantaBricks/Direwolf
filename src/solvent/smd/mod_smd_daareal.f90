! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_smd_daareal: analytic accessible-solid-angle SASA (Liotard 1992 /

module mod_smd_daareal
use mod_smd_tables, only: PI_smd
implicit none

contains

subroutine smd_daareal(natoms, RAD, RIJ, URIJ, k, AREA0, NCROSS, NC, DAREA)
integer,intent(in) :: natoms, k
real(8),intent(in) :: RAD(natoms), RIJ(natoms,natoms), URIJ(3,natoms,natoms)
real(8),intent(out) :: AREA0
integer,intent(out) :: NCROSS
integer,intent(out) :: NC(0:natoms)
real(8),intent(out) :: DAREA(3,0:natoms)

real(8),parameter :: EPSI = 1.0d-11
integer :: MXSS
integer,allocatable :: LAB(:), NCNCT(:,:)
logical,allocatable :: CONECT(:,:)
real(8),allocatable :: CTHETA(:,:), STHETA(:), SIT(:)
real(8),allocatable :: DCSIT(:,:), DJCOSN(:,:,:,:), COSN(:,:,:), DSTETA(:,:)
real(8),allocatable :: DCTETA(:,:,:), DCOSN(:,:,:,:)
real(8),allocatable :: WORK(:), DIWORK(:,:), DJWORK(:,:), DKWORK(:,:), D0WORK(:,:)
real(8),allocatable :: DCAODD(:,:), DCAPLY(:,:), DICOSN(:,:,:,:), DCASLC(:,:)

integer :: i, j, l, ll, m, icor, jcor, lt, ia, ib, ibold, li, lj, lk
integer :: nclust, nfree, npoly, nphi
real(8) :: RK, RK2, RKINV, RIKINV, EPSK, X, Y, TWOPI, FOURPI
real(8) :: CISJ, SICJ, SISJ, TIJ, EPSIJ, ASLICE, APOLY, AEVEN, AODD
real(8) :: DRCTHT, COSNI, COSNIJ, SIN2IJ, AIJ, BIJ, CIJ, DCIJ, DSN2IJ
real(8) :: C2I, S1NPHI, PHI, PHI1, PHI2, CHEK, BUF
real(8) :: DX, DY, DIX, DIY, DJX, DJY, DKX, DKY, DIS2IJ, DJS2IJ
real(8) :: VN(3), DIVN(3,3), DJVN(3,3), VIJ(3), DIVIJ(3,3), DJVIJ(3,3)
real(8) :: CNIJ(3), DICNIJ(3,3), DJCNIJ(3,3), CNIK(3), DICNIK(3,3), DKCNIK(3,3)
real(8) :: DAPHI(3), DBPHI(3), DICOS(3,3), DJCOS(3,3), DIWIJ(3,3), DJWIJ(3,3)
real(8) :: DIAIJ(3), DJAIJ(3), DIBIJ(3), DJBIJ(3), DICIJ(3), DJCIJ(3)
real(8) :: buf3(3)
logical :: FREEIJ, FREEJI, LPOLY

MXSS = 2*natoms+1
allocate(LAB(natoms), NCNCT(MXSS,natoms), CONECT(natoms,natoms))
allocate(CTHETA(natoms,natoms), STHETA(natoms), SIT(natoms))
allocate(DCSIT(3,natoms), DJCOSN(3,3,natoms,natoms), COSN(3,natoms,natoms), DSTETA(3,natoms))
allocate(DCTETA(3,natoms,natoms), DCOSN(3,3,natoms,natoms))
allocate(WORK(MXSS), DIWORK(3,MXSS), DJWORK(3,MXSS), DKWORK(3,MXSS), D0WORK(3,MXSS))
allocate(DCAODD(3,0:natoms), DCAPLY(3,0:natoms), DICOSN(3,3,natoms,natoms), DCASLC(3,0:natoms))

LAB = 0; STHETA = 0.0d0; SIT = 0.0d0
DCSIT = 0.0d0; DSTETA = 0.0d0
NCNCT = 0
CONECT = .false.
CTHETA = 0.0d0
COSN = 0.0d0; DCTETA = 0.0d0
DJCOSN = 0.0d0; DCOSN = 0.0d0; DICOSN = 0.0d0
WORK = 0.0d0
DCAODD = 0.0d0; DCAPLY = 0.0d0; DCASLC = 0.0d0

TWOPI = PI_smd+PI_smd
FOURPI = TWOPI+TWOPI
RK = RAD(k)
NC(0) = k
NCROSS = 0
DAREA(:,0) = 0.0d0

if (RK <= 0.0d0) then
   AREA0 = 0.0d0
   deallocate(LAB,NCNCT,CONECT,CTHETA,STHETA,SIT,DCSIT,DJCOSN,COSN,DSTETA, &
              DCTETA,DCOSN,WORK,DIWORK,DJWORK,DKWORK,D0WORK,DCAODD,DCAPLY,DICOSN,DCASLC)
   return
endif
AREA0 = FOURPI

10 continue
NCROSS = 0
EPSK = EPSI*RK

do i = 1,natoms
   if (i == k .or. RAD(i) <= 0.0d0) cycle
   if (RK+RAD(i)-RIJ(i,k) < EPSK) cycle
   if (RIJ(i,k)-abs(RK-RAD(i)) < EPSK) then
      if (RK <= RAD(i)) then
         AREA0 = 0.0d0
         deallocate(LAB,NCNCT,CONECT,CTHETA,STHETA,SIT,DCSIT,DJCOSN,COSN,DSTETA, &
                    DCTETA,DCOSN,WORK,DIWORK,DJWORK,DKWORK,D0WORK,DCAODD,DCAPLY,DICOSN,DCASLC)
         return
      endif
   else
      NCROSS = NCROSS+1
      NC(NCROSS) = i
      NCNCT(NCROSS,1) = i
   endif
enddo
if (NCROSS == 0) then
   deallocate(LAB,NCNCT,CONECT,CTHETA,STHETA,SIT,DCSIT,DJCOSN,COSN,DSTETA, &
              DCTETA,DCOSN,WORK,DIWORK,DJWORK,DKWORK,D0WORK,DCAODD,DCAPLY,DICOSN,DCASLC)
   return
endif

RKINV = 0.5d0/RK
RK2 = RK**2
do i = 1,NCROSS
   li = NCNCT(i,1)
   RIKINV = 1.0d0/RIJ(li,k)
   CTHETA(i,i) = RKINV*(RIJ(li,k)+(RK2-RAD(li)**2)*RIKINV)
   STHETA(i) = sqrt(1.0d0-CTHETA(i,i)**2)
   CONECT(i,i) = .false.
   COSN(1,i,i) = URIJ(1,k,li); COSN(2,i,i) = URIJ(2,k,li); COSN(3,i,i) = URIJ(3,k,li)
   X = -CTHETA(i,i)/STHETA(i)
   DRCTHT = RKINV*(1.0d0-(RK2-RAD(li)**2)*RIKINV**2)
   do icor = 1,3
      DCTETA(icor,i,i) = DRCTHT*COSN(icor,i,i)
      DSTETA(icor,i) = X*DCTETA(icor,i,i)
      COSNI = COSN(icor,i,i)*RIKINV
      do jcor = 1,icor
         COSNIJ = -COSNI*COSN(jcor,i,i)
         DCOSN(icor,jcor,i,i) = COSNIJ
         DCOSN(jcor,icor,i,i) = COSNIJ
      enddo
      DCOSN(icor,icor,i,i) = DCOSN(icor,icor,i,i) + RIKINV
   enddo
enddo

do i = 2,NCROSS
   do j = 1,i-1
      if (CONECT(j,j)) cycle
      CISJ = CTHETA(i,i)*STHETA(j)
      SICJ = STHETA(i)*CTHETA(j,j)
      SISJ = STHETA(i)*STHETA(j)
      CTHETA(j,i) = dot_product(COSN(:,i,i),COSN(:,j,j))
      do icor = 1,3
         DCTETA(icor,i,j) = dot_product(DCOSN(:,icor,i,i),COSN(:,j,j))
         DCTETA(icor,j,i) = dot_product(COSN(:,i,i),DCOSN(:,icor,j,j))
      enddo
      TIJ = CTHETA(j,i)-CTHETA(i,i)*CTHETA(j,j)
      if (TIJ > SISJ-EPSI*abs(CISJ-SICJ)) then
         if (CTHETA(j,j) > CTHETA(i,i)) then
            CONECT(j,j) = .true.
         else
            CONECT(i,i) = .true.
            cycle
         endif
      else
         EPSIJ = EPSI*(SICJ+CISJ)
         if (SICJ+CISJ >= 0.0d0) then
            CONECT(j,i) = TIJ > EPSIJ-SISJ
         else
            if (TIJ <= -SISJ-EPSIJ) then
               NCROSS = 0
               AREA0 = 0.0d0
               deallocate(LAB,NCNCT,CONECT,CTHETA,STHETA,SIT,DCSIT,DJCOSN,COSN,DSTETA, &
                          DCTETA,DCOSN,WORK,DIWORK,DJWORK,DKWORK,D0WORK,DCAODD,DCAPLY,DICOSN,DCASLC)
               return
            else
               CONECT(j,i) = .true.
            endif
         endif
         CONECT(i,j) = CONECT(j,i)
      endif
   enddo
enddo

ASLICE = 0.0d0
DCASLC(:,0:NCROSS) = 0.0d0
NCLUST = 0
do i = 1,NCROSS
   if (CONECT(i,i)) cycle
   block
      logical :: found
      found = .false.
      do j = 1,NCROSS
         if (CONECT(j,j)) cycle
         if (CONECT(j,i)) then
            NCLUST = NCLUST+1
            LAB(NCLUST) = i
            WORK(NCLUST) = CTHETA(i,i)+EPSI*STHETA(i)
            WORK(NCLUST+NCROSS) = CTHETA(i,i)-EPSI*STHETA(i)
            found = .true.
            exit
         endif
      enddo
      if (found) cycle
   end block
   ASLICE = ASLICE + 1.0d0-CTHETA(i,i)
   do icor = 1,3
      DCASLC(icor,i) = DCASLC(icor,i) - DCTETA(icor,i,i)
      DCASLC(icor,0) = DCASLC(icor,0) + DCTETA(icor,i,i)
   enddo
enddo
ASLICE = TWOPI*ASLICE
DCASLC(:,0:NCROSS) = DCASLC(:,0:NCROSS)*TWOPI

if (NCLUST == 0) then
   AREA0 = FOURPI-ASLICE
   do i = 0,NCROSS
      DAREA(:,i) = -DCASLC(:,i)
   enddo
   deallocate(LAB,NCNCT,CONECT,CTHETA,STHETA,SIT,DCSIT,DJCOSN,COSN,DSTETA, &
              DCTETA,DCOSN,WORK,DIWORK,DJWORK,DKWORK,D0WORK,DCAODD,DCAPLY,DICOSN,DCASLC)
   return
endif

NFREE = 0
NCNCT(MXSS,LAB(1)) = 0
do i = 2,NCLUST
   li = LAB(i)
   NCNCT(MXSS,li) = 0
   do j = 1,i-1
      lj = LAB(j)
      if (.not. CONECT(lj,li)) cycle
      SIN2IJ = 1.0d0/(1.0d0-CTHETA(lj,li)**2)
      AIJ = (CTHETA(li,li)-CTHETA(lj,lj)*CTHETA(lj,li))*SIN2IJ
      BIJ = (CTHETA(lj,lj)-CTHETA(li,li)*CTHETA(lj,li))*SIN2IJ
      CIJ = sqrt((1.0d0-AIJ*CTHETA(li,li)-BIJ*CTHETA(lj,lj))*SIN2IJ)
      VN(1) = COSN(2,li,li)*COSN(3,lj,lj)-COSN(3,li,li)*COSN(2,lj,lj)
      VN(2) = COSN(3,li,li)*COSN(1,lj,lj)-COSN(1,li,li)*COSN(3,lj,lj)
      VN(3) = COSN(1,li,li)*COSN(2,lj,lj)-COSN(2,li,li)*COSN(1,lj,lj)
      DCIJ = 0.5d0/CIJ
      do icor = 1,3
         DIVN(1,icor) = DCOSN(2,icor,li,li)*COSN(3,lj,lj)-DCOSN(3,icor,li,li)*COSN(2,lj,lj)
         DIVN(2,icor) = DCOSN(3,icor,li,li)*COSN(1,lj,lj)-DCOSN(1,icor,li,li)*COSN(3,lj,lj)
         DIVN(3,icor) = DCOSN(1,icor,li,li)*COSN(2,lj,lj)-DCOSN(2,icor,li,li)*COSN(1,lj,lj)
         DJVN(1,icor) = COSN(2,li,li)*DCOSN(3,icor,lj,lj)-COSN(3,li,li)*DCOSN(2,icor,lj,lj)
         DJVN(2,icor) = COSN(3,li,li)*DCOSN(1,icor,lj,lj)-COSN(1,li,li)*DCOSN(3,icor,lj,lj)
         DJVN(3,icor) = COSN(1,li,li)*DCOSN(2,icor,lj,lj)-COSN(2,li,li)*DCOSN(1,icor,lj,lj)
      enddo
      DSN2IJ = 2.0d0*CTHETA(lj,li)*SIN2IJ

      do icor = 1,3
         COSN(icor,li,lj) = AIJ*COSN(icor,li,li)+BIJ*COSN(icor,lj,lj)+CIJ*VN(icor)
         COSN(icor,lj,li) = AIJ*COSN(icor,li,li)+BIJ*COSN(icor,lj,lj)-CIJ*VN(icor)
         DIS2IJ = DSN2IJ*DCTETA(icor,li,lj)
         DJS2IJ = DSN2IJ*DCTETA(icor,lj,li)
         DIAIJ(icor) = (DCTETA(icor,li,li)-CTHETA(lj,lj)*DCTETA(icor,li,lj))*SIN2IJ + AIJ*DIS2IJ
         DJAIJ(icor) = (-DCTETA(icor,lj,lj)*CTHETA(lj,li)-CTHETA(lj,lj)*DCTETA(icor,lj,li))*SIN2IJ + AIJ*DJS2IJ
         DJBIJ(icor) = (DCTETA(icor,lj,lj)-CTHETA(li,li)*DCTETA(icor,lj,li))*SIN2IJ + BIJ*DJS2IJ
         DIBIJ(icor) = (-DCTETA(icor,li,li)*CTHETA(lj,li)-CTHETA(li,li)*DCTETA(icor,li,lj))*SIN2IJ + BIJ*DIS2IJ
         DICIJ(icor) = -DCIJ*((DIAIJ(icor)*CTHETA(li,li)+AIJ*DCTETA(icor,li,li)+DIBIJ(icor)*CTHETA(lj,lj))*SIN2IJ) &
                       + 0.5d0*CIJ*DIS2IJ
         DJCIJ(icor) = -DCIJ*((DJBIJ(icor)*CTHETA(lj,lj)+BIJ*DCTETA(icor,lj,lj)+DJAIJ(icor)*CTHETA(li,li))*SIN2IJ) &
                       + 0.5d0*CIJ*DJS2IJ
      enddo

      do icor = 1,3
         do jcor = 1,3
            DICOS(icor,jcor) = DIAIJ(jcor)*COSN(icor,li,li) + AIJ*DCOSN(icor,jcor,li,li) + DIBIJ(jcor)*COSN(icor,lj,lj)
            DJCOS(icor,jcor) = DJAIJ(jcor)*COSN(icor,li,li) + DJBIJ(jcor)*COSN(icor,lj,lj) + BIJ*DCOSN(icor,jcor,lj,lj)
            DIWIJ(icor,jcor) = DICIJ(jcor)*VN(icor) + CIJ*DIVN(icor,jcor)
            DJWIJ(icor,jcor) = DJCIJ(jcor)*VN(icor) + CIJ*DJVN(icor,jcor)
            DICOSN(icor,jcor,li,lj) = DICOS(icor,jcor) + DIWIJ(icor,jcor)
            DJCOSN(icor,jcor,li,lj) = DJCOS(icor,jcor) + DJWIJ(icor,jcor)
            DICOSN(icor,jcor,lj,li) = DICOS(icor,jcor) - DIWIJ(icor,jcor)
            DJCOSN(icor,jcor,lj,li) = DJCOS(icor,jcor) - DJWIJ(icor,jcor)
         enddo
      enddo

      FREEIJ = .true.
      FREEJI = .true.

      block
         logical :: skip_rest
         skip_rest = .false.
         do l = 1,NCLUST
            if (l == i .or. l == j) cycle
            ll = LAB(l)
            if (.not.(CONECT(ll,li) .and. CONECT(ll,lj))) cycle
            if (FREEJI) then
               CHEK = dot_product(COSN(:,lj,li),COSN(:,ll,ll))
               if (CHEK > WORK(l)) then
                  FREEJI = .false.
               else
                  if (CHEK >= WORK(l+NCROSS)) then
                     RK = RK*(1.0d0+4.0d0*EPSI)
                     skip_rest = .true.
                     exit
                  endif
               endif
            endif
            if (FREEIJ) then
               CHEK = dot_product(COSN(:,li,lj),COSN(:,ll,ll))
               if (CHEK > WORK(l)) then
                  FREEIJ = .false.
               else if (CHEK >= WORK(l+NCROSS)) then
                  RK = RK*(1.0d0+4.0d0*EPSI)
                  skip_rest = .true.
                  exit
               endif
            endif
            if (.not.(FREEIJ .or. FREEJI)) exit
         enddo
         if (skip_rest) then
            deallocate(LAB,NCNCT,CONECT,CTHETA,STHETA,SIT,DCSIT,DJCOSN,COSN,DSTETA, &
                       DCTETA,DCOSN,WORK,DIWORK,DJWORK,DKWORK,D0WORK,DCAODD,DCAPLY,DICOSN,DCASLC)
            goto 10
         endif
      end block

      CTHETA(li,lj) = CTHETA(lj,li)
      if (FREEJI) then
         NFREE = NFREE+1
         m = NCNCT(MXSS,li)+1; NCNCT(m,li) = lj; NCNCT(MXSS,li) = m
         m = NCNCT(MXSS,lj)+1; NCNCT(m,lj) = -li; NCNCT(MXSS,lj) = m
      endif
      if (FREEIJ) then
         NFREE = NFREE+1
         m = NCNCT(MXSS,li)+1; NCNCT(m,li) = -lj; NCNCT(MXSS,li) = m
         m = NCNCT(MXSS,lj)+1; NCNCT(m,lj) = li; NCNCT(MXSS,lj) = m
      endif
   enddo
enddo
if (NFREE == 0) then
   NCROSS = 0
   AREA0 = 0.0d0
   deallocate(LAB,NCNCT,CONECT,CTHETA,STHETA,SIT,DCSIT,DJCOSN,COSN,DSTETA, &
              DCTETA,DCOSN,WORK,DIWORK,DJWORK,DKWORK,D0WORK,DCAODD,DCAPLY,DICOSN,DCASLC)
   return
endif

APOLY = 0.0d0
DCAPLY(:,0:NCROSS) = 0.0d0
do i = 1,NCLUST
   DCAODD(:,0:NCROSS) = 0.0d0
   li = LAB(i)
   nphi = NCNCT(MXSS,li)
   if (nphi == 0) cycle
   lj = NCNCT(1,li)
   if (lj > 0) then
      CNIJ = COSN(:,lj,li)
      if (li >= lj) then
         DICNIJ = DICOSN(:,:,lj,li); DJCNIJ = DJCOSN(:,:,lj,li)
      else
         DICNIJ = DJCOSN(:,:,lj,li); DJCNIJ = DICOSN(:,:,lj,li)
      endif
   else
      lj = -lj
      CNIJ = COSN(:,li,lj)
      if (li >= lj) then
         DICNIJ = DICOSN(:,:,li,lj); DJCNIJ = DJCOSN(:,:,li,lj)
      else
         DICNIJ = DJCOSN(:,:,li,lj); DJCNIJ = DICOSN(:,:,li,lj)
      endif
      NCNCT(1,li) = lj
   endif
   SIT(li) = 1.0d0/STHETA(li)
   C2I = CTHETA(li,li)**2
   VIJ(1) = COSN(2,li,li)*CNIJ(3)-COSN(3,li,li)*CNIJ(2)
   VIJ(2) = COSN(3,li,li)*CNIJ(1)-COSN(1,li,li)*CNIJ(3)
   VIJ(3) = COSN(1,li,li)*CNIJ(2)-COSN(2,li,li)*CNIJ(1)
   LPOLY = dot_product(VIJ,COSN(:,lj,lj)) > 0.0d0

   X = 2.0d0*CTHETA(li,li)
   Y = -SIT(li)**2
   do icor = 1,3
      DCSIT(icor,li) = Y*DSTETA(icor,li)
      DIVIJ(1,icor) = DCOSN(2,icor,li,li)*CNIJ(3)-DCOSN(3,icor,li,li)*CNIJ(2) &
                     + COSN(2,li,li)*DICNIJ(3,icor)-COSN(3,li,li)*DICNIJ(2,icor)
      DIVIJ(2,icor) = DCOSN(3,icor,li,li)*CNIJ(1)-DCOSN(1,icor,li,li)*CNIJ(3) &
                     + COSN(3,li,li)*DICNIJ(1,icor)-COSN(1,li,li)*DICNIJ(3,icor)
      DIVIJ(3,icor) = DCOSN(1,icor,li,li)*CNIJ(2)-DCOSN(2,icor,li,li)*CNIJ(1) &
                     + COSN(1,li,li)*DICNIJ(2,icor)-COSN(2,li,li)*DICNIJ(1,icor)
      DJVIJ(1,icor) = COSN(2,li,li)*DJCNIJ(3,icor)-COSN(3,li,li)*DJCNIJ(2,icor)
      DJVIJ(2,icor) = COSN(3,li,li)*DJCNIJ(1,icor)-COSN(1,li,li)*DJCNIJ(3,icor)
      DJVIJ(3,icor) = COSN(1,li,li)*DJCNIJ(2,icor)-COSN(2,li,li)*DJCNIJ(1,icor)
   enddo

   do j = 2,nphi
      lj = NCNCT(j,li)
      if (lj > 0) then
         CNIK = COSN(:,lj,li)
         if (li >= lj) then
            DICNIK = DICOSN(:,:,lj,li); DKCNIK = DJCOSN(:,:,lj,li)
         else
            DICNIK = DJCOSN(:,:,lj,li); DKCNIK = DICOSN(:,:,lj,li)
         endif
      else
         lj = -lj
         CNIK = COSN(:,li,lj)
         if (li >= lj) then
            DICNIK = DICOSN(:,:,li,lj); DKCNIK = DJCOSN(:,:,li,lj)
         else
            DICNIK = DJCOSN(:,:,li,lj); DKCNIK = DICOSN(:,:,li,lj)
         endif
         NCNCT(j,li) = lj
      endif

      X = dot_product(VIJ,CNIK)
      Y = dot_product(CNIK,CNIJ)-C2I
      WORK(j-1) = atan2(X,Y)
      if (WORK(j-1) <= 0.0d0) WORK(j-1) = WORK(j-1)+TWOPI
      DX = -X/(X**2+Y**2)
      DY = Y/(X**2+Y**2)
      do icor = 1,3
         DIX = dot_product(DIVIJ(:,icor),CNIK)+dot_product(VIJ,DICNIK(:,icor))
         DIY = dot_product(DICNIK(:,icor),CNIJ)+dot_product(CNIK,DICNIJ(:,icor))-2.0d0*CTHETA(li,li)*DCTETA(icor,li,li)
         DJX = dot_product(DJVIJ(:,icor),CNIK)
         DJY = dot_product(CNIK,DJCNIJ(:,icor))
         DKX = dot_product(VIJ,DKCNIK(:,icor))
         DKY = dot_product(DKCNIK(:,icor),CNIJ)
         DIWORK(icor,j-1) = DY*DIX+DX*DIY
         DJWORK(icor,j-1) = DY*DJX+DX*DJY
         DKWORK(icor,j-1) = DY*DKX+DX*DKY
         D0WORK(icor,j-1) = -DIWORK(icor,j-1)-DJWORK(icor,j-1)-DKWORK(icor,j-1)
      enddo
   enddo

   if (nphi == 2) then
      AODD = WORK(1)
      lj = NCNCT(1,li); lk = NCNCT(2,li)
      do icor = 1,3
         DCAODD(icor,li) = DCAODD(icor,li)+DIWORK(icor,1)
         DCAODD(icor,lj) = DCAODD(icor,lj)+DJWORK(icor,1)
         DCAODD(icor,lk) = DCAODD(icor,lk)+DKWORK(icor,1)
         DCAODD(icor,0) = DCAODD(icor,0)+D0WORK(icor,1)
      enddo
   else
      do j = 1,nphi-2
         do l = j+1,nphi-1
            if (WORK(j) > WORK(l)) then
               BUF = WORK(l); WORK(l) = WORK(j); WORK(j) = BUF
               m = NCNCT(l+1,li); NCNCT(l+1,li) = NCNCT(j+1,li); NCNCT(j+1,li) = m
               buf3 = DIWORK(:,l); DIWORK(:,l) = DIWORK(:,j); DIWORK(:,j) = buf3
               buf3 = DJWORK(:,l); DJWORK(:,l) = DJWORK(:,j); DJWORK(:,j) = buf3
               buf3 = DKWORK(:,l); DKWORK(:,l) = DKWORK(:,j); DKWORK(:,j) = buf3
               buf3 = D0WORK(:,l); D0WORK(:,l) = D0WORK(:,j); D0WORK(:,j) = buf3
            endif
         enddo
      enddo

      AODD = WORK(1)
      lj = NCNCT(1,li); lk = NCNCT(2,li)
      do icor = 1,3
         DCAODD(icor,li) = DCAODD(icor,li)+DIWORK(icor,1)
         DCAODD(icor,lj) = DCAODD(icor,lj)+DJWORK(icor,1)
         DCAODD(icor,lk) = DCAODD(icor,lk)+DKWORK(icor,1)
         DCAODD(icor,0) = DCAODD(icor,0)+D0WORK(icor,1)
      enddo
      do j = 3,nphi-1,2
         AODD = AODD+WORK(j)-WORK(j-1)
         lk = NCNCT(j,li); ll = NCNCT(j+1,li)
         do icor = 1,3
            DCAODD(icor,li) = DCAODD(icor,li)+DIWORK(icor,j)-DIWORK(icor,j-1)
            DCAODD(icor,lj) = DCAODD(icor,lj)+DJWORK(icor,j)-DJWORK(icor,j-1)
            DCAODD(icor,lk) = DCAODD(icor,lk)-DKWORK(icor,j-1)
            DCAODD(icor,ll) = DCAODD(icor,ll)+DKWORK(icor,j)
            DCAODD(icor,0) = DCAODD(icor,0)+D0WORK(icor,j)-D0WORK(icor,j-1)
         enddo
      enddo
   endif

   AEVEN = TWOPI-AODD
   X = 1.0d0-CTHETA(li,li)
   if (LPOLY) then
      APOLY = APOLY+AODD
      ASLICE = ASLICE+AEVEN*X
      do lt = 0,NCROSS
         DCAPLY(:,lt) = DCAPLY(:,lt)+DCAODD(:,lt)
         DCASLC(:,lt) = DCASLC(:,lt)-DCAODD(:,lt)*X
      enddo
      DCASLC(:,li) = DCASLC(:,li)-AEVEN*DCTETA(:,li,li)
      DCASLC(:,0) = DCASLC(:,0)+AEVEN*DCTETA(:,li,li)
   else
      APOLY = APOLY+AEVEN
      ASLICE = ASLICE+AODD*X
      do lt = 0,NCROSS
         DCAPLY(:,lt) = DCAPLY(:,lt)-DCAODD(:,lt)
         DCASLC(:,lt) = DCASLC(:,lt)+DCAODD(:,lt)*X
      enddo
      DCASLC(:,li) = DCASLC(:,li)-AODD*DCTETA(:,li,li)
      DCASLC(:,0) = DCASLC(:,0)+AODD*DCTETA(:,li,li)
      m = NCNCT(1,li)
      do j = 1,nphi-1
         NCNCT(j,li) = NCNCT(j+1,li)
      enddo
      NCNCT(nphi,li) = m
   endif
enddo

NPOLY = 0
do i = 1,NCLUST
   li = LAB(i)
   do j = 2,NCNCT(MXSS,li),2
      if (NCNCT(j,li) == 0) cycle
      ia = NCNCT(j-1,li)
      ib = li
      NCNCT(j,li) = 0
      PHI1 = CTHETA(ia,ib)-CTHETA(ia,ia)*CTHETA(ib,ib)
      PHI2 = SIT(ia)*SIT(ib)
      PHI = acos(PHI1*PHI2)
      APOLY = APOLY+PHI
      S1NPHI = 1.0d0/sin(PHI)
      PHI1 = S1NPHI*PHI1
      PHI2 = S1NPHI*PHI2
      do icor = 1,3
         DAPHI(icor) = (DCTETA(icor,ia,ib)-CTHETA(ib,ib)*DCTETA(icor,ia,ia))*PHI2+PHI1*DCSIT(icor,ia)*SIT(ib)
         DBPHI(icor) = (DCTETA(icor,ib,ia)-CTHETA(ia,ia)*DCTETA(icor,ib,ib))*PHI2+PHI1*SIT(ia)*DCSIT(icor,ib)
         DCAPLY(icor,ia) = DCAPLY(icor,ia)-DAPHI(icor)
         DCAPLY(icor,ib) = DCAPLY(icor,ib)-DBPHI(icor)
         DCAPLY(icor,0) = DCAPLY(icor,0)+DAPHI(icor)+DBPHI(icor)
      enddo

      do
         l = 2
         do
            if (l > NCNCT(MXSS,ia)) exit
            if (NCNCT(l,ia) == ib) then
               ibold = ib
               NCNCT(l,ia) = 0
               ib = ia
               ia = NCNCT(l-1,ib)
               if (ia /= ibold) then
                  PHI1 = CTHETA(ia,ib)-CTHETA(ia,ia)*CTHETA(ib,ib)
                  PHI2 = SIT(ia)*SIT(ib)
                  PHI = acos(PHI1*PHI2)
                  S1NPHI = 1.0d0/sin(PHI)
                  PHI1 = S1NPHI*PHI1
                  PHI2 = S1NPHI*PHI2
                  do icor = 1,3
                     DAPHI(icor) = (DCTETA(icor,ia,ib)-CTHETA(ib,ib)*DCTETA(icor,ia,ia))*PHI2+PHI1*DCSIT(icor,ia)*SIT(ib)
                     DBPHI(icor) = (DCTETA(icor,ib,ia)-CTHETA(ia,ia)*DCTETA(icor,ib,ib))*PHI2+PHI1*SIT(ia)*DCSIT(icor,ib)
                  enddo
               else
                  ib = li
                  ia = NCNCT(j-1,li)
               endif
               APOLY = APOLY+PHI
               do icor = 1,3
                  DCAPLY(icor,ia) = DCAPLY(icor,ia)-DAPHI(icor)
                  DCAPLY(icor,ib) = DCAPLY(icor,ib)-DBPHI(icor)
                  DCAPLY(icor,0) = DCAPLY(icor,0)+DAPHI(icor)+DBPHI(icor)
               enddo
               exit
            endif
            l = l+2
         enddo
         if (l > NCNCT(MXSS,ia)) then
            exit
         endif
      enddo
      NPOLY = NPOLY+1
   enddo
enddo

APOLY = APOLY+real(NPOLY-NFREE,8)*TWOPI
AREA0 = FOURPI-ASLICE-mod(APOLY,FOURPI)
do i = 0,NCROSS
   DAREA(:,i) = -DCASLC(:,i)-DCAPLY(:,i)
enddo

deallocate(LAB,NCNCT,CONECT,CTHETA,STHETA,SIT,DCSIT,DJCOSN,COSN,DSTETA, &
           DCTETA,DCOSN,WORK,DIWORK,DJWORK,DKWORK,D0WORK,DCAODD,DCAPLY,DICOSN,DCASLC)
end subroutine smd_daareal

end module mod_smd_daareal
