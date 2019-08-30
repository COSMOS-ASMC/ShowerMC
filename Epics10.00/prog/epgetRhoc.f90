  subroutine epgetRhoc(inputname, name, rhoc)
    implicit none
    character(len=*),intent(in):: inputname ! Media anme;may be like Air*0.01
    character(len=*),intent(out):: name  ! *0.01 is dropped
                           ! it may be inutname.
    real(4),intent(out):: rhoc ! 0.01 part.  if *0.01 part is
                    ! missing, rhoc=1 
    
    integer:: posast
    posast = index(inputname, "*")  ! Air*0.05 type ?
    if(posast > 0 ) then
       read(inputname(posast+1:),*) rhoc
       name = inputname(1:posast-1) 
    else
       name = trim(inputname)
       rhoc = 1.0
    endif
  end subroutine epgetRhoc
