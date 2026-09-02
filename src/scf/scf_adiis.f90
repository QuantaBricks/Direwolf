! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! scf_adiis: ADIIS (Hu and Yang, J. Chem. Phys. 132, 054109 (2010); Garza


subroutine build_adiis_terms(n,nconts,Pa_hist,Fa_hist,Pb_hist,Fb_hist,open_shell,lin,quad)
implicit none
integer,intent(in) :: n,nconts
real(4),intent(in) :: Pa_hist(nconts,nconts,n),Fa_hist(nconts,nconts,n)
real(4),intent(in) :: Pb_hist(nconts,nconts,n),Fb_hist(nconts,nconts,n)
logical,intent(in) :: open_shell
real(8),intent(out) :: lin(n),quad(n,n)
real(8),allocatable :: dPa(:,:,:),dPb(:,:,:)
integer :: i,j
real(8) :: qraw

allocate(dPa(nconts,nconts,n))
do i = 1,n
   dPa(:,:,i) = real(Pa_hist(:,:,i),8) - real(Pa_hist(:,:,n),8)
enddo
if (open_shell) then
   allocate(dPb(nconts,nconts,n))
   do i = 1,n
      dPb(:,:,i) = real(Pb_hist(:,:,i),8) - real(Pb_hist(:,:,n),8)
   enddo
endif

do i = 1,n
   lin(i) = 2.0d0*sum(dPa(:,:,i)*real(Fa_hist(:,:,n),8))
   if (open_shell) lin(i) = lin(i) + 2.0d0*sum(dPb(:,:,i)*real(Fb_hist(:,:,n),8))
enddo

do i = 1,n
   do j = 1,n
      qraw = sum(dPa(:,:,i)*(real(Fa_hist(:,:,j),8)-real(Fa_hist(:,:,n),8)))
      if (open_shell) qraw = qraw + sum(dPb(:,:,i)*(real(Fb_hist(:,:,j),8)-real(Fb_hist(:,:,n),8)))
      quad(i,j) = qraw
   enddo
enddo
quad = quad + transpose(quad)

deallocate(dPa)
if (open_shell) deallocate(dPb)
end subroutine build_adiis_terms

function adiis_fx(n,Amat,bvec,x) result(val)
implicit none
integer,intent(in) :: n
real(8),intent(in) :: Amat(n,n),bvec(n),x(n)
real(8) :: val
val = 0.5d0*dot_product(x,matmul(Amat,x)) + dot_product(bvec,x)
end function adiis_fx

subroutine adiis_powell(n,Amat,bvec,x)
implicit none
integer,intent(in) :: n
real(8),intent(in) :: Amat(n,n),bvec(n)
real(8),intent(out) :: x(n)
integer,parameter :: MAX_MACRO = 500
real(8),parameter :: DF_TOL = 1.0d-8
real(8) :: xguess(n,2*n+1),yguess(2*n+1)
real(8) :: cvec(n),direction(n),new_x(n),old_x(n),Ax(n)
real(8) :: cur_val,new_val,curval_before,dE,step,denom
integer :: i,imacro,best_idx,nguess
real(8),external :: adiis_fx

if (n .eq. 1) then
   x(1) = 1.0d0
   return
endif

nguess = 1
xguess(:,nguess) = 1.0d0/n
do i = 1,n
   nguess = nguess+1
   xguess(:,nguess) = 1.0d0/(n+2)
   xguess(i,nguess) = xguess(i,nguess)*3.0d0
enddo
do i = 1,n
   nguess = nguess+1
   xguess(:,nguess) = 0.0d0
   xguess(i,nguess) = 1.0d0
enddo

do i = 1,nguess
   yguess(i) = adiis_fx(n,Amat,bvec,xguess(:,i))
enddo
best_idx = 1
do i = 2,nguess
   if (yguess(i) .lt. yguess(best_idx)) best_idx = i
enddo
x = xguess(:,best_idx)

cur_val = adiis_fx(n,Amat,bvec,x)
old_x = x

do imacro = 1,MAX_MACRO
   curval_before = cur_val
   do i = 1,n
      cvec = 0.0d0
      cvec(i) = 1.0d0
      direction = cvec - x
      denom = dot_product(direction,matmul(Amat,direction))
      if (denom .eq. 0.0d0) cycle
      Ax = matmul(Amat,x)
      step = -(dot_product(direction,Ax) + dot_product(bvec,direction))/denom
      if (step .ne. step) cycle
      if (step .gt. 0.0d0 .and. step .le. 1.0d0) then
         new_x = x + step*direction
         new_val = adiis_fx(n,Amat,bvec,new_x)
         if (new_val .lt. cur_val) then
            x = new_x
            cur_val = new_val
         endif
      endif
   enddo

   dE = cur_val - curval_before

   direction = x - old_x
   denom = dot_product(direction,matmul(Amat,direction))
   if (denom .ne. 0.0d0) then
      Ax = matmul(Amat,x)
      step = -(dot_product(direction,Ax) + dot_product(bvec,direction))/denom
      if (step .eq. step .and. step .gt. 0.0d0 .and. step .le. 1.0d0) then
         new_x = x + step*direction
         new_val = adiis_fx(n,Amat,bvec,new_x)
         if (new_val .lt. cur_val) then
            x = new_x
            cur_val = new_val
            dE = cur_val - curval_before
         endif
      endif
   endif
   old_x = x

   if (dE .gt. -DF_TOL) exit
enddo

end subroutine adiis_powell
