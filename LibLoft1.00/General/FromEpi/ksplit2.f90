!  implicit none
!  integer::nr, i, icon
!  character(len=50):: string
!  character(len=8)::  b(20)
!  read(*,'(a)') string
!  write(0,'(a)') string
!!
!  call ksplit2(string, '|', b,  20, nr,icon)
!  write(0,*) ' icon from split2=',icon
!  do i = 1, nr
!     write(*,*) i, "* ", b(i), " *"
!     call krmLeadingB(b(i), b(i), icon)
!     write(*,*) i, "* ", b(i), " *"
!  enddo
!end program
subroutine ksplit2(a, sep,  b,  n,  nr, icon)
  implicit none
  character(len=*),intent(in):: a   ! . e.g   abc|xyz |  xxkxkxkk
!    if sep is |, then a is
!          abc| ssss|
!      or
!         |abc| ssss|
!    result is same. 
  character(len=1) sep ! input. non blank separator. e.g |
  integer,intent(in):: n   !   array size of b    
  character(len=*),intent(out):: b(n) ! field  value of i=1,2,...nr
  integer,intent(out):: nr  ! nr field is obtained (nr<=n)
  integer,intent(out):: icon ! 0 : ok
                  !  1:  some of b(i) could not
                  !     hold all of the string
                  !  -1: # of field > n. nr=n and others are
                  !     discarded 
  integer ::tl, bp, i, nc
  character(len=1):: prev
  integer:: m 
  m = len(b(1))

  tl = len(trim(a))

  nr = 0
  prev = sep
  icon = 0
  do i = 1, tl
     if( a(i:i) == sep ) then
        prev = sep
     else
        if(prev == sep ) then
           if(nr >=  n)  then
              icon = -1
              return
           endif
           nr = nr + 1
           nc = 0
           prev = 'x'
        endif
        nc = nc + 1
        if(nc == 1)  b(nr) = ' '
        if(nc <= m) then
           b(nr)(nc:nc) = a(i:i)
        else
           icon = 1
        endif
     endif
  enddo
end subroutine ksplit2
