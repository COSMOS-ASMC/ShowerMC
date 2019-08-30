!  Make  consts table needed  for multiple scattering treatent
! based on a  modified Wentzel electron scattering model.
! (see PENELOPE-2011, Data Bank NEA/NSC/DOC(2011)5)
!  Main input is transport cross-sections S0,S1 and S2.
! from which paramters used in the MW are fixed and tabullated
! as a function of energy for each media. (i.e various Z or
! mixture of Z's).
!
! mixture or compound treatment: 
! Compute transport cross-sections for each Zi in the matter
! and add them with weight of Ni
!  S0 =  Sum( Nix S0i)
!  S1 =  sum( NixS1i) 
!  S2 =  sum( NixS2i) 
!
!  mu =( 1-cos )/2 
  subroutine cMWscatFixAB(S0, S1, S2, MWavemu, MWavemu2, avemu, avemu2, A0, A, B)
    implicit none
!     This fixes the const A,B in MW DCS using given S0,S1,S2 
! 
!

    real(8),intent(in):: S0  ! 0th transport x-sec. (i.e total
                ! elastic xs) in mb or any unit but cm2 is defatul
    real(8),intent(in):: S1  ! 1st transport x-sec. unit same as S0
    real(8),intent(in):: S2  ! 2nd transport x-sec. unit same as S0
    real(8),intent(out):: MWavemu ! MW's <mu>  (not correct one !!)
                 ! A0{(1+A0) log(1+A0)/A0 -1}
    real(8),intent(out):: MWavemu2 ! MW's <mu^2>  (not correct one !!)
    real(8),intent(out):: avemu ! <mu>
    real(8),intent(out):: avemu2 ! <mu2>
    real(8),intent(out):: A0  ! A0 obtained by cMWscatFixA0
    real(8),intent(out):: A, B   ! non dim.
    real(8),external:: cMWcaseIIA
    real(8):: param(2)
    real(8),save::eps=1.d-6
    integer:: icon
!1!!! 
    real(8):: Ax
!!!!
    avemu =  S1/S0/2
    avemu2 = avemu - S2/S0/6
    call  cMWscatFixA0(S0, S1, avemu, A0)      
    MWavemu = A0*( (1+A0)* log( (1+A0)/A0 ) - 1)
    MWavemu2 = A0*(1-2*MWavemu)
    if( MWavemu2 > avemu2 ) then
       ! MW distribution is too wide. put peak pos near zero by delta func
       A= A0
       B = (MWavemu2 - avemu2)/(MWavemu2 - avemu**2) ! OK
    else
       param(1) = avemu ! =S1/S0/2 
       param(2) = avemu2 ! mu -S2/S0/6.    ! <mu2>
                          !                avoid 0 for min
       call kbinChopWP(cMWcaseIIA, param, 1.0d-30, A0, A0/2, eps, A, icon)
!!!!!!!!!!!
       if( A>= A0) then
          Ax = A0
          write(0,*) ' A=',A, ' > A0=',A0
          do while (Ax > 1.0e-30)
             write(0,*) Ax, cMWcaseIIA(Ax, param)
             Ax = Ax/10.0d0**0.01
          enddo
       end if
!!!!!!!!!!
       MWavemu = A*( (1+A)* log( (1+A)/A ) - 1)
       MWavemu2 = A*(1-2*MWavemu)
       B = (avemu-MWavemu)/(5.0/6.0 - MWavemu)

       if( icon < 0 ) then
          write(0,*) ' icon=',icon,'  cMWcaseIIA '
          write(0,*) ' param=', param, ' A0=',A0, ' A=',A, 'B=',B
!          write(0,*) ' at 1.0d-100', cMWcaseIIA(1.0d-100, param)
!          write(0,*) ' at A0', cMWcaseIIA(A0, param)
!          Ax = A0/2.
!          do while (Ax > 0.) 
!             write(0,*)  Ax, cMWcaseIIA(Ax, param)
!             Ax = Ax/2
!          enddo
!
!          stop
       endif
       B = -B  ! for case II put negative sign
    endif
!    write(0,*) ' A0=', A0, ' A =', A,'B=', B
  end subroutine cMWscatFixAB

 subroutine cMWscatFixA0(S0, S1, avemu, A0)
    implicit none
!     This fixes the const A0 in MW DCS using given S0,S1
!   normalized  DCS is A0(1+A0)/(mu +A0)**2 
!
    real(8),intent(in):: S0  ! 0th transport x-sec. (i.e total
                ! elastic xs) in mb or any unit
    real(8),intent(in):: S1  ! 1st transport x-sec. unit same as S0

    real(8),intent(out):: avemu  ! <mu>
    real(8),intent(out):: A0    ! non dim.

    real(8):: param(1)
    real(8),save:: x0=-1.  
    real(8),parameter:: eps = 1.d-6
    real(8),external:: cMWscatA0, cMWscatA0p
    integer:: icon

    avemu= S1/S0/2

!     MW's  <mu> is A0( (1+A0)log((1+A0)/A0) -1 ) (which is 0~0.5) 
!        this should be avemu of precise calc. 
!    this must be solved numerically to fix A0.  
    if( avemu >= 0.5e0 ) then
       write(0,*) '<mu> is too large=',avemu
       write(0,*) ' must be < 0.5 in cMWscatPrep'
       write(0,*) ' S0, S1=',S0, S1
       stop
    endif

    param(1) = avemu
!         if mu is small, x2 becomes < 0; to avoid
!         it, initial guess must be small.
!         after getting correct answer, it is fitted
!         like this.
    if( x0 < 0.) then
       x0 =6.083e-7*(avemu/1e-5)
    endif
    call kNewtonRaphson(cMWscatA0, cMWscatA0p, param, x0, eps, A0, icon)
    if(icon < 0 ) then
       write(0, *) ' in cMWscatFixA0 icon =', icon, ' A0=',A0
       write(0,*) ' for <mu>=',avemu
       if( icon < -1 ) stop
    endif
    x0= A0
  end subroutine cMWscatFixA0

  function cMWscatA0(x, param) result(ans)
!          function used to get A0 
    implicit none
    real(8),intent(in):: x  ! variable stands for  A0
    real(8),intent(in):: param(1)  ! <mu>
    real(8):: ans  !  A0 
    if(x == 0.) then
!       ans = -param(1)
       ans = -1.0
    else
       ans = x* ((1+x)*log( (1+x)/x ) - 1 )
       ans = ans/ param(1) -1.0
    endif
  end function cMWscatA0

  function cMWscatA0p(x, param) result(ans)
!       dcMWscatA/dx
    implicit none
    real(8),intent(in):: x
    real(8),intent(in):: param(1)
    
    real(8):: ans
    real(8):: temp

    if(x == 0.) then
       ans = 1.d10  ! big
    else
       temp = log ( (1+x)/x )
       ans = (1+x*2)*temp  - 2.
       ans = ans/param(1)
    endif
  end function cMWscatA0p



  function cMWcaseIIA(x, param)  result(ans)
    implicit none
    real(8),intent(in):: x
    real(8),intent(in):: param(2) ! param(1)=<mu>, param(2)=<mu2>
    
    real(8)::ans

    real(8)::avemuW, avemu2W, mu, mu2, B

    avemuW = x*(  (1+x)*log((1+x)/x) -1 )
    avemu2W = x*(1- 2*avemuW)
    mu = param(1)
    mu2 = param(2)
    B = ( mu-avemuW)/( 5.e0/6.e0 - avemuW)
    ans =( (1.d0-B)*avemu2W + B*17.e0/24.e0 )/mu2 - 1.0
  end function cMWcaseIIA
