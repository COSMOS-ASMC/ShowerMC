  subroutine cTPXS_Z2file(Z, pm, filename)
!      get file name in which  transportxs, E, s0,s1,s2  are contained
    implicit none
    integer,intent(in):: Z  ! charge of target
    integer,intent(in):: pm  ! >0 for e+, <0 for e-
    character(*),intent(out):: filename

    if( pm <0 ) then
       write(filename,'(a,i0.2,a)')  &
         "$LIBLOFT/Data/Elsepa/Z",Z,"/e-/tcstable.dat"
    elseif( pm > 0) then
       write(filename,'(a,i0.2,a)')  &
         "$LIBLOFT/Data/Elsepa/Z",Z,"/e+/tcstable.dat"
    else
       write(0,*) ' pm=',pm, ' invalid to cTPXS_Z2file'
       stop
    endif
  end subroutine cTPXS_Z2file

 
