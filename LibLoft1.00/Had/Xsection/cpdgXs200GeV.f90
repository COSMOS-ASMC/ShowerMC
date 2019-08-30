!  sigma pA(=nA) @ 200 GeV by PDB as a funciton of A (mb) 
!  as of 2010; test prog at bottom.
  function cPDGsigmaTotpA(A) result(sigma)  ! total mb
    implicit none
    real(8),intent(in):: A  ! mass #
    real(8):: sigma
    real(8):: aa, b

    sigma = 0.
    if( A > 100.) then  ! better than 0.1 %
       aa=330.55d0
       b=0.729d0
    elseif(A >25.) then  ! bettter than  0.3%
       aa=298.53d0
       b=0.7757d0
    elseif( abs(A-16.) < 0.1) then
       sigma = 433.4d0
    elseif(A > 10.) then   ! better than 0.1%
       aa=290.86d0
       b=0.7991d0 
    elseif( A > 7.0) then  !better than 0.1 %
       aa=292.68d0
       b=0.76905d0
    elseif( abs(A-4.0)<0.1 ) then
       sigma = 128.2d0
    elseif( abs(A-2.0)<0.1 ) then
       sigma = 64.7
    elseif( A > 2.9 ) then  ! bettern than 0.3% 
                  !         A=4: 128.2 A=2 64.7
       aa = 308.26
       b = 0.97
    elseif( abs(A-1.0) < 0.1) then
       sigma = 38.8
    else
       write(0,*) ' A=', A, 'wrong in cPDGSigmaTotpA'
       stop
    endif
    if(sigma == 0. ) then
       sigma = aa*(A/10.d0) ** b
    endif
  end function cPDGsigmaTotpA
  
  function cPDGsigmaInepA(A) result(sigma)  ! inela mb
    implicit none
    real(8),intent(in):: A  ! mass #
    real(8):: sigma
    real(8):: aa, b

    sigma = 0.
    ! better than 0.25 (some ~0.5%)
    if( A > 51.0d0) then 
       aa =217.87d0
       b = 0.68280d0
    elseif( A > 16.9d0 ) then
       aa=204.93d0
       b=0.7169d0
    elseif( abs(A-16.) < 0.1) then
       sigma = 294.56
    elseif( A > 6.0 ) then
       aa=206.2d0
       b=0.672d0
    elseif( abs(A-7.) < 0.1 ) then
       sigma = 164.8
    elseif( abs(A-4.) < 0.1 ) then
       sigma =93.55
    elseif( abs(A-2.) < 0.1 ) then
       sigma = 46.25
    elseif( A > 2.9 ) then
       aa =233.1d0
       b = 1.0d0
    elseif( abs(A-1.) < 0.1 ) then
       sigma = 31.93
    else
       write(0,*) ' A=', A, 'wrong in cPDGSigmaInepA'
       stop
    endif
    if(sigma == 0.) then
       sigma = aa*(A/10.)**b
    endif
  end function cPDGsigmaInepA
!  
!   program main
!     implicit none
!     real(8):: A, sigmaT, sigmaI
!     integer:: i
!     real(8):: cPDGsigmaTotpA,cPDGsigmaInepA
! 
!     do i = 1, 210
!        A = i
!        sigmaT  = cPDGsigmaTotpA(A)
!        sigmaI  = cPDGsigmaInepA(A)
! 
!        write(*,'(1p,4g12.3)') A, sigmaT, sigmaI, sigmaI/sigmaT
!     enddo
!   end program main


  





  
