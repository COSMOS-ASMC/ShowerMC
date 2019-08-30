!program main
!implicit none
!real(8):: sb  
!real(8):: B  
! sb = 1.1d0
! do while (sb < 35.)
!    call epBlogB(sb, B)
!    write(*,*) sb, B, (B-log(B)-sb)/B
!    sb = sb+0.05d0
! enddo
!end program main
!

subroutine epMolBlogB(sb, B)
implicit none
real(8),intent(in):: sb  ! solve  B- log(B) = sb for given sb
real(8),intent(out):: B  ! 

real(8)::sbb, fB, fBp
if(sb <= 1.d0) then
   write(0,*) ' b for B-log(B) = b is too small b=',b
   stop
endif
if(sb < 5.d0) then
!     sb<2.5-->error <0.06 %
!       <5     error <0.005% 
   sbb =sb
   B = -0.591837 + (2.06538 -0.11690*sbb)*sbb
   fB = B - log(B) -sb
   fBp = 1.d0 - 1.d0/B
   B =B-  fB/fBp
!   sbb = B - log(B)
!   B = -0.591837 + (2.06538 -0.11690*sbb)*sbb
elseif( sb < 15.d0) then
!  5<sb<15;  next one;  |error|  mostly 0.1% max 0.3%
   B = 1.154804 + (1.183689 -0.00462208*sb)*sb
!    by next, error becomes  almost 0
   fB = B - log(B) -sb
   fBp = 1.d0 - 1.d0/B
   B =B-  fB/fBp
elseif( sb < 50.d0 ) then
!   |error| < 0.02 % max 0.05 %
   B = 1.9368162469 +(1.074931428 -0.0007441876*sb)*sb
!    by next almost 0
   fB = B - log(B) -sb
   fBp = 1.d0 - 1.d0/B
   B =B-  fB/fBp
else
   write(0,*) ' sb for B-log(B)=sb too large=',sb
   stop
endif
end subroutine epMolBlogB
