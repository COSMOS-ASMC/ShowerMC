  subroutine cDCS_Z2file(Z, pm, ntherg)
!      make file name for charge Z, e-/e=, and   given energy index
!      file name is put in filename in modTPXS
    use modTPXS
    implicit none
    integer,intent(in):: Z  ! charge of target
    integer,intent(in):: pm  ! >0 for e+, <0 for e-
    integer,intent(in):: ntherg ! energy node index; e%g 1~106  
              ! energy is KEele(ntherg)
!    character(*),intent(out):: filename --->one in modTPXS
    character(10):: temp
    character(10):: cerg
    integer:: ii

    if( pm <0 ) then
       write(filename,'(a,i0.2,a)')  &
         "$LIBLOFT/Data/Elsepa/Z",Z,"/e-/dcs_"
    elseif( pm > 0) then
       write(filename,'(a,i0.2,a)')  &
         "$LIBLOFT/Data/Elsepa/Z",Z,"/e+/dcs_"
    else
       write(0,*) ' pm=',pm, ' invalid to cTPXS_Z2file'
       stop
    endif
    write(temp,'(1p,E10.3)') KEele(ntherg)  !  1.250E+08
    call kgsub(".", "p", temp, cerg)        !  1p250E+08
    temp = cerg
    call kgsub("E+", "e", temp, cerg)       !  1p250e08
    call kseblk(cerg, "?", ii)              !1p250e08
    filename=trim(filename)//trim(cerg)//".dat"
  end subroutine cDCS_Z2file
