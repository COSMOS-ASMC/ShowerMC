!!  implicit none
!!  real(8):: x, t1, t2
!!  logical:: kphiinside
!!  write(0,*) ' Enter t1, t2'
!!  read(*,*) t1, t2
!!  do while(.true.)
!!     write(0,*) 'Enter x '
!!     read(*,*) x
!!     write(0,*) x, kphiinside(x, t1, t2)
!!  enddo
!!end program
function kphiinside(x, t1, t2) result(ans)
  implicit none      
!        check if asimuthal angle x is between t1 and t2
  real(8),intent(in):: x  ! azimuthal ang. in deg.
  real(8),intent(in):: t1,t2 !  range of azimuthal angle
                       !  in deg.   t2 may be less than  t1.
                      !  in that case (say t1=270, t2= 30) 
                      ! angle range is t1 to t2+360 )
  logical:: ans   !  if x is between t1 and t2, (equal is
                       ! includedt else f.

  real(8):: xx, tt1, tt2
!        0 <= xx < 360
  xx =  mod( mod(x, 360.d0)+ 360.d0 , 360.d0 )
  tt1 =  mod( mod(t1, 360.d0)+ 360.d0 , 360.d0 )
  tt2 =  mod( mod(t2, 360.d0)+ 360.d0 , 360.d0 )
  if( tt2 < tt1 ) then
     if( tt1<= xx .and. xx <360.d0 ) then
        ans = .true.
     elseif( xx >=0. .and. xx <= tt2 )  then
        ans = .true.
     else
        ans = .false.
     endif
  else
     if( tt1 <= xx .and. xx <= tt2 ) then
        ans = .true.
     else
        ans = .false.
     endif
  endif
end   function kphiinside
