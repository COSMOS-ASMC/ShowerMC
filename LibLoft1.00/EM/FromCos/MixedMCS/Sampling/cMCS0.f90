module modMCS0
  implicit none
  integer,parameter:: usize=155
!  integer,parameter:: usize=159
  real(8):: uarray(usize)
  real(8),parameter:: umin=1.0d-5
  real(8):: du
end module modMCS0

subroutine ciniSmpTab
  use modMCS0
  implicit none
  integer::i, j
  logical,save::first=.true.
  real(8):: u
  if(first) then
     uarray(1) =0.
     du = umin
     u = 0.
     j = 0
     do i=2, 31
        j = j +1
        u = u + du
        uarray(i) = u
        if(mod(j,10) == 0 ) then
           du = uarray(i)
!           write(0,*) ' i=',i, ' du=',du
           j =1
        endif
     enddo
     i=31 
!     write(0,*) ' i =',i, ' du ua=', du, uarray(i)
     do while (u < 0.9d0) 
        u = u + du
        i = i + 1
        uarray(i) = u
     enddo
!     write(0,*) ' i, uarray(i) =',  i, uarray(i)
     do i = 155, 119, -1
        uarray(i) = 1.0d0 - uarray(156-i)
     enddo
!!!! 
!     uarray(159)=1.d0
!     do i = 159,156,-1
!        uarray(i-1) = uarray(i)-2.d-6
!     enddo   
     first=.false.
     !!!!!!!
!     write(0,*)' uarray(154:159)='
!     write(0,*)   uarray(154:159)
  end if
end subroutine ciniSmpTab




