!  program test
!  implicit none
!  character(len=10):: a
!  character(len=7):: b
!  integer::icon
!  read(*,'(a)')  a
!  write(*,*) ' a=',a
!  call krmLeadingB(a, b, icon)
!  write(*,*) 'icon=',icon
!  write(*,*) ' b=',b
!  write(*,*) ' a=',a
!  
!end program test
  subroutine krmLeadingB(a, b, icon)
!    remove leading blanks in a  and put it in b.
    implicit none
    character(len=*),intent(in):: a
    character(len=*),intent(out):: b  ! b can be a.
    integer,intent(out):: icon ! if b cannot contain the
              ! result fully, icon = 1 else 0.
    integer:: li, lo, i
    
    if( a == ' ' ) then
       b = a
       icon = 0
    else
       lo = len(b)
       li = len( trim(a) )
       do i = 1, li
          if( a(i:i) == ' ' ) cycle
          if(len(trim(a(i:))) > lo) then
             icon = 1
          else
             icon = 0
          endif
          b(1:) = a(i:)
          exit
       enddo
    endif
  end subroutine krmLeadingB
