! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_smd_sts: per-atom surface-tension STS(i) and its gradient - the

module mod_smd_sts
use mod_smd_tables, only: NATCNV, smd_rkkval
implicit none

contains

subroutine smd_sts_and_grad(natoms, Z, RIJ, URIJ, SIGMA, HSIGMA, STS, DSTS)
integer,intent(in) :: natoms, Z(natoms)
real(8),intent(in) :: RIJ(natoms,natoms), URIJ(3,natoms,natoms)
real(8),intent(in) :: SIGMA(150), HSIGMA(150)
real(8),intent(out) :: STS(natoms), DSTS(3,natoms,natoms)

real(8) :: RKKVAL(15,15)
real(8) :: COT(natoms,natoms), DCOTDR(natoms,natoms)
integer :: i, j, k, l, itpc, jtpc, ktpc, ntp, ntpc, nktpc, kccc
real(8) :: RHLD, DELTAR, CUTOF1, CUTOF2, EXPONT, EXPVAL, BCCC
real(8) :: RTKK, RTKK2, RTKK3, RTKK5, DRTKK, RTKKS, RTKKS2, DHOLDER
real(8) :: SCOTC, SCOTO, COTJK, DCOTJK, RHLD2, HOLDRK

RKKVAL = smd_rkkval()
EXPVAL = exp(1.0d0)
BCCC = 1.0d0

COT = 0.0d0
DCOTDR = 0.0d0
STS = 0.0d0
DSTS = 0.0d0

do i = 1,natoms
   itpc = NATCNV(Z(i))
   do j = 1,natoms
      if (i == j) cycle
      jtpc = NATCNV(Z(j))
      if (itpc == 0 .or. jtpc == 0) cycle
      DELTAR = 0.30d0
      RHLD = RKKVAL(itpc,jtpc)
      if ((Z(i) == 8 .and. Z(j) == 6) .or. (Z(i) == 6 .and. Z(j) == 8)) then
         DELTAR = 0.10d0; RHLD = 1.330d0
      endif
      if (Z(i) == 8 .and. Z(j) == 8) then
         DELTAR = 0.30d0; RHLD = 1.80d0
      endif
      CUTOF1 = RHLD + DELTAR
      if (RIJ(i,j) < CUTOF1) then
         EXPONT = DELTAR/(RIJ(i,j)-CUTOF1)
         COT(i,j) = exp(EXPONT)
         DCOTDR(i,j) = -COT(i,j)*EXPONT/(RIJ(i,j)-CUTOF1)
      endif
   enddo
enddo

do i = 1,natoms
   STS(i) = SIGMA(Z(i))

   if (Z(i) == 1) then
      do j = 1,natoms
         ntp = Z(j)
         STS(i) = STS(i) + COT(i,j)*HSIGMA(ntp)
         do l = 1,3
            DSTS(l,i,i) = DSTS(l,i,i) + HSIGMA(ntp)*DCOTDR(i,j)*URIJ(l,j,i)
            DSTS(l,j,i) = DSTS(l,j,i) - HSIGMA(ntp)*DCOTDR(i,j)*URIJ(l,j,i)
         enddo
      enddo
   endif

   if (Z(i) == 8) then
      do j = 1,natoms
         ntp = Z(j)
         if (ntp == 6) then
            STS(i) = STS(i) + COT(i,j)*SIGMA(103)
            do l = 1,3
               DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(103)*DCOTDR(i,j)*URIJ(l,j,i)
               DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(103)*DCOTDR(i,j)*URIJ(l,j,i)
            enddo
         else if (ntp == 7) then
            STS(i) = STS(i) + COT(i,j)*SIGMA(106)
            do l = 1,3
               DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(106)*DCOTDR(i,j)*URIJ(l,j,i)
               DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(106)*DCOTDR(i,j)*URIJ(l,j,i)
            enddo
         else if (ntp == 15) then
            STS(i) = STS(i) + COT(i,j)*SIGMA(114)
            do l = 1,3
               DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(114)*DCOTDR(i,j)*URIJ(l,j,i)
               DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(114)*DCOTDR(i,j)*URIJ(l,j,i)
            enddo
         else if (ntp == 8 .and. i /= j) then
            STS(i) = STS(i) + COT(i,j)*SIGMA(104)
            do l = 1,3
               DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(104)*DCOTDR(i,j)*URIJ(l,j,i)
               DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(104)*DCOTDR(i,j)*URIJ(l,j,i)
            enddo
         endif
      enddo
   endif

   if (Z(i) == 7) then
      RTKKS = 0.0d0; RTKKS2 = 0.0d0
      do j = 1,natoms
         ntp = Z(j)
         if (ntp /= 6) cycle
         RTKK3 = 0.0d0; RTKK5 = 0.0d0
         do k = 1,natoms
            if (k == i .or. k == j) cycle
            nktpc = NATCNV(Z(k)); ntpc = NATCNV(Z(j))
            RHLD2 = RKKVAL(nktpc,ntpc)
            CUTOF2 = RHLD2 + 0.3d0
            if (RIJ(j,k) < CUTOF2) then
               RTKK2 = exp(1.0d0-(BCCC/(1.0d0-((RIJ(j,k)-RHLD2)/0.3d0))))/EXPVAL
               if (Z(k) == 8) RTKK5 = RTKK5 + RTKK2
               RTKK3 = RTKK3 + RTKK2
            endif
         enddo
         RTKKS = RTKKS + COT(i,j)*RTKK3**2
         RTKKS2 = RTKKS2 + COT(i,j)*RTKK5
      enddo
      DHOLDER = RTKKS**1.3d0
      STS(i) = STS(i) + DHOLDER*SIGMA(105)
      STS(i) = STS(i) + RTKKS2*SIGMA(111)

      do j = 1,natoms
         ntp = Z(j)
         if (ntp /= 6) cycle
         SCOTC = 0.0d0; SCOTO = 0.0d0
         do k = 1,natoms
            if (k == i .or. k == j) cycle
            nktpc = NATCNV(Z(k)); ntpc = NATCNV(Z(j))
            RHLD2 = RKKVAL(nktpc,ntpc)
            CUTOF2 = RHLD2 + 0.30d0
            COTJK = 0.0d0; DCOTJK = 0.0d0
            if (RIJ(j,k) < CUTOF2) then
               COTJK = exp(1.0d0-(BCCC/(1.0d0-((RIJ(j,k)-RHLD2)/0.30d0))))/EXPVAL
               DCOTJK = -COTJK*0.3d0/(RIJ(j,k)-CUTOF2)**2
            endif
            if (Z(k) == 8) then
               SCOTO = SCOTO + COTJK
               do l = 1,3
                  DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(111)*COTJK*DCOTDR(i,j)*URIJ(l,j,i) &
                               + SIGMA(111)*COT(i,j)*DCOTJK*URIJ(l,k,j)
                  DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(111)*COTJK*DCOTDR(i,j)*URIJ(l,j,i)
                  DSTS(l,k,i) = DSTS(l,k,i) - SIGMA(111)*COT(i,j)*DCOTJK*URIJ(l,k,j)
               enddo
            endif
            SCOTC = SCOTC + COTJK
         enddo
         do l = 1,3
            DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(105)*RTKKS**0.3d0*1.3d0*SCOTC**2*DCOTDR(i,j)*URIJ(l,j,i)
            DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(105)*RTKKS**0.3d0*1.3d0*SCOTC**2*DCOTDR(i,j)*URIJ(l,j,i)
         enddo
         do k = 1,natoms
            if (k == i .or. k == j) cycle
            nktpc = NATCNV(Z(k)); ntpc = NATCNV(Z(j))
            RHLD2 = RKKVAL(nktpc,ntpc)
            CUTOF2 = RHLD2 + 0.30d0
            COTJK = 0.0d0; DCOTJK = 0.0d0
            if (RIJ(j,k) < CUTOF2) then
               COTJK = exp(1.0d0-(1.0d0/(1.0d0-((RIJ(j,k)-RHLD2)/0.30d0))))/EXPVAL
               DCOTJK = -COTJK*0.3d0/(RIJ(j,k)-CUTOF2)**2
            endif
            do l = 1,3
               DSTS(l,j,i) = DSTS(l,j,i) + SIGMA(105)*RTKKS**0.3d0*1.3d0*2.0d0*SCOTC*COT(i,j)*DCOTJK*URIJ(l,k,j)
               DSTS(l,k,i) = DSTS(l,k,i) - SIGMA(105)*RTKKS**0.3d0*1.3d0*2.0d0*SCOTC*COT(i,j)*DCOTJK*URIJ(l,k,j)
            enddo
         enddo
      enddo

      RTKK = 0.0d0
      DELTAR = 0.065d0
      RHLD = 1.225d0
      CUTOF1 = RHLD+DELTAR
      do j = 1,natoms
         if (j == i) cycle
         if (Z(j) /= 6) cycle
         if (RIJ(i,j) < CUTOF1) then
            RTKK = RTKK + exp(1.0d0-(BCCC/(1.0d0-((RIJ(i,j)-RHLD)/DELTAR))))/EXPVAL
            EXPONT = DELTAR/(RIJ(i,j)-CUTOF1)
            DRTKK = -exp(EXPONT)*EXPONT/(RIJ(i,j)-CUTOF1)
            do l = 1,3
               DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(116)*DRTKK*URIJ(l,j,i)
               DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(116)*DRTKK*URIJ(l,j,i)
            enddo
         endif
      enddo
      STS(i) = STS(i) + RTKK*SIGMA(116)
   endif

   if (Z(i) == 6) then
      itpc = NATCNV(6)
      do kccc = 1,2
         RTKK = 0.0d0
         do j = 1,natoms
            if (j == i) cycle
            if (Z(j) /= 6) cycle
            ntpc = NATCNV(Z(j))
            if (kccc == 1) then
               RHLD = RKKVAL(itpc,ntpc); DELTAR = 0.30d0
            else
               RHLD = 1.27d0; DELTAR = 0.07d0
            endif
            CUTOF1 = RHLD+DELTAR
            if (RIJ(i,j) < CUTOF1) then
               RTKK = RTKK + exp(1.0d0-(BCCC/(1.0d0-((RIJ(i,j)-RHLD)/DELTAR))))/EXPVAL
               EXPONT = DELTAR/(RIJ(i,j)-CUTOF1)
               DRTKK = -exp(EXPONT)*EXPONT/(RIJ(i,j)-CUTOF1)
               if (kccc == 1) then
                  do l = 1,3
                     DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(101)*DRTKK*URIJ(l,j,i)
                     DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(101)*DRTKK*URIJ(l,j,i)
                  enddo
               else
                  do l = 1,3
                     DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(102)*DRTKK*URIJ(l,j,i)
                     DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(102)*DRTKK*URIJ(l,j,i)
                  enddo
               endif
            endif
         enddo
         if (kccc == 1) then
            STS(i) = STS(i) + RTKK*SIGMA(101)
         else
            STS(i) = STS(i) + RTKK*SIGMA(102)
         endif
      enddo

      RTKK = 0.0d0
      do j = 1,natoms
         if (Z(j) == 7) RTKK = RTKK + COT(i,j)
      enddo
      HOLDRK = RTKK*RTKK
      STS(i) = STS(i) + HOLDRK*SIGMA(110)
      do j = 1,natoms
         if (Z(j) /= 7) cycle
         do l = 1,3
            DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(110)*2.0d0*RTKK*DCOTDR(i,j)*URIJ(l,j,i)
            DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(110)*2.0d0*RTKK*DCOTDR(i,j)*URIJ(l,j,i)
         enddo
      enddo
   endif

   if (Z(i) == 16) then
      do j = 1,natoms
         if (Z(j) == 16 .and. j /= i) then
            STS(i) = STS(i) + COT(i,j)*SIGMA(107)
            do l = 1,3
               DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(107)*DCOTDR(i,j)*URIJ(l,j,i)
               DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(107)*DCOTDR(i,j)*URIJ(l,j,i)
            enddo
         endif
      enddo
      do j = 1,natoms
         if (Z(j) == 15) then
            STS(i) = STS(i) + COT(i,j)*SIGMA(115)
            do l = 1,3
               DSTS(l,i,i) = DSTS(l,i,i) + SIGMA(115)*DCOTDR(i,j)*URIJ(l,j,i)
               DSTS(l,j,i) = DSTS(l,j,i) - SIGMA(115)*DCOTDR(i,j)*URIJ(l,j,i)
            enddo
         endif
      enddo
   endif
enddo
end subroutine smd_sts_and_grad

end module mod_smd_sts
