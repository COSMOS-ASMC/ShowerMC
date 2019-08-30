module  SpAbySpp
  implicit none
! p    1   1     35.67        -1.681        0.2172    
! p    2   1     59.42        -1.902        0.3404    
! p    4   1     100.3        -1.712        0.4928    
! p    5   1     118.9        -1.929        0.5493    
! p    6   1     136.7       -0.8193        0.5983    
! p    8   1     170.5        -1.257        0.6816    
! p   10   1     202.5        -2.016        0.7515    
!   read data such as above
!   and compute Spp and SpA and output Spp  vs SpA/Spp
!
  integer,parameter::Ei=3
  integer,parameter::Ai=210
!    coef. for S= a+ b*log(E/GeV) + c*log(E/GeV)**2
  real(8),save:: a(Ei,Ai) 
  real(8),save:: b(Ei,Ai) 
  real(8),save:: c(Ei,Ai) 
  real(8),parameter::E1=50.d0
  real(8),parameter::E2=1.d5
  real(8),parameter::E3=1.d8
  real(8),parameter::E4=1.d11

contains
  function sigma(E, TA) result(ans)
    implicit none
    real(8),intent(in)::E  ! p energy in GeV
    real(8),intent(in)::TA  ! mass number for pA collision
    real(8)::ans  ! cross-section in mb

    integer::iE  
    integer::iA


    if(E< E1) then
       stop "too low E"
    elseif( E< E2 ) then
       iE = 1
    elseif( E< E3) then
       iE = 2
    elseif( E<= E4) then
       iE = 3
    else
       stop "too high E"
    endif
    iA = TA
    if(a(iE,iA) == 0 ) then
       ans = 0.
    else
       ans = (log(E)*c(iE,iA) + b(iE,iA))*log(E) +a(iE,iA)
    endif
  end function sigma
  subroutine readTBL(more)
    implicit  none
    logical,intent(out):: more 
    integer::iA, iE
    character(2)::id
!  p    1   1     35.67        -1.681        0.2172        
    read(*, *, end=100 ) id, iA, iE, &
         &    a(iE,iA), b(iE,iA), c(iE,iA)
    more =.true.
    return
100 continue
    more = .false.
  end subroutine readTBL
end module SpAbySpp
program main
  use SpAbySpp
  implicit none
  logical:: readmore=.true.
  
  real(8),parameter::dE=10.d0**0.1d0
  real(8),parameter::HdE=10.d0**0.05d0
  real(8)::E, TA, Spp, SpA
  integer::iA, i

  a(:,:) = 0.
  b(:,:) = 0.
  c(:,:) = 0.

  do while(readmore)
     call readTBL(readmore)
  enddo



  
  E = E1
  do while(E < E4)
     if(E > E4/HdE .and. E< E4*HdE ) then
        E = E4
     endif
     Spp = sigma(E, 1.d0)
     do iA = 2, 210
        TA = iA
        SpA = sigma(E, TA)
        if(SpA > 0.) then
           write(*,'(1p,2g14.4, 0p,i4,1p,g14.4)')  &
                &  Spp, SpA/Spp, iA, E
        endif
     enddo
     E = E*dE
  enddo
end program main



    

