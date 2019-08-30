  subroutine kNewtonRaphson(func, dfunc, param, x0, eps, ans, icon)
    implicit none
   ! solve  func(x, param) = 0 
    real(8):: func, dfunc  ! external function. with two  arguments
               ! dfunc is derivative func(x,param)  defined dfunc(x, param)
               ! with the same param.
    real(8),intent(in)::param(*) ! some possible parameter  needed to compute
               !         func and dfunc.  if no param is needed. give dummy
    real(8),intent(in):: x0   ! first guess of root
    real(8),intent(in):: eps   ! such as 1d-6. abso or relative error of root
                            !  if > 0, |func(root)| < eps is required
                            !  if < 0, |(x1-x2)/x2| < |eps| is required if
                            !           |x2|>1.  else |x1-x2| < |eps|
                            !For,  eps>0, to solve f(x) = c,
                            ! better to use func= f(x)/c -1
    real(8),intent(out):: ans  ! obtained root
    integer,intent(out):: icon ! >  if ans obtained. icon shows # of iteration
                               ! -1  no convergence. but current root is given

    real(8):: x1, x2, tang
    integer,parameter::maxite = 30
    integer:: n
    real(8):: funcx

!       first guess
    x1 = x0
    icon = -1
    n = 0 
    do while (n < maxite)
       tang = dfunc(x1, param)  ! tangent at x1
!       zero point of the tangential line
       funcx =func(x1, param)
!       x2 = -func(x1, param)/tang + x1
       x2 = -funcx/tang + x1
       if(eps > 0. ) then
          if( abs(funcx) < eps ) then
             icon = n
             exit
          endif
       else   
          if( abs(x2)> 1.) then
             if( abs((x1-x2)/x2) < abs(eps)) then
                icon = n
                exit
             endif
          else
             if( abs(x1-x2) < abs(eps) ) then
                icon = n
                exit
             endif
          endif
       endif
       x1 = x2
       n = n + 1
    enddo
    ans = x2
  end subroutine kNewtonRaphson
!
!  program test
!    implicit none
!    real(8):: muave, param(1)
!    real(8),external:: myfunc, mydfunc
!    real(8)::a, eps, ans
!    integer::icon
!
!
!    do
!       eps = 1.d-6
!       write(0,*) 'Enter muave'
!       read(*,*) muave
!       if( muave == 0. ) stop
!       param(1) = muave
!!         if mu is small, x2 becomes < 0; to avoid
!!         it, initial guess must be small.
!!         after getting correct answer, it is fitted
!!         like this.
!       a =6.083e-7*(muave/1e-5)
!!       if(a< 1.e-4) eps = a/1e4
!       call kNewtonRaphson(myfunc, mydfunc, param, a, eps, ans, icon)
!       write(*,'(i3, 1p, 3g13.4)') icon, muave, ans, myfunc(ans, param)
!       call kbinChopWP(myfunc, param,1.d-10, 1.0d4, 0.1d0, eps, ans, icon)
!       write(*,'(i3, 1p, 3g13.4)') icon, muave, ans, myfunc(ans, param)
!    enddo
!  end program test
!  function myfunc(x, param) result(ans)
!    implicit none
!    real(8),intent(in):: x
!    real(8),intent(in):: param(1)
!    
!    real(8):: ans
!    if(x == 0.) then
!       ans = -param(1)
!    else
!!       ans = x* ((1+x)*log( (1+x)/x ) - 1 ) - param(1)
!       ans = x* ((1+x)*log( (1+x)/x ) - 1 )
!       ans = ans/ param(1) -1.0
!    endif
!  end function myfunc
!
!  function mydfunc(x, param) result(ans)
!    implicit none
!    real(8),intent(in):: x
!    real(8),intent(in):: param(1)
!    
!    real(8):: ans
!    real(8):: temp
!    
!    real(8),external:: myfunc
!
!    if(x == 0.) then
!       ans = 1.d10  ! big
!    else
!       temp = log ( (1+x)/x )
!       ans = (1+x*2)*temp  - 2.
!!
!       ans = ans/param(1)
!    endif
!  end function mydfunc
!
!!      Binary Chop for getting a solution of  f(x) = 0.
!!     This version permit 2 parameters for f
!subroutine kbinChopWP(f, param,  a, b, x, eps, ans, icon)
!  implicit none
!!
!  real(8):: f  ! input. function name. to be used as f(x, param)
!  !   f(x, param) = 0 is solved.
!  real(8),intent(in):: param(*)
!  real(8),intent(in):: a ! input. lower bound of solution
!  real(8),intent(in):: b ! input. upper bound of solution
!  real(8),intent(in):: x   ! input. initial guess of  solution.  
!  real(8),intent(in):: eps  ! input. relative error of solution.
!  real(8),intent(out):: ans  ! output. obtained solution
!  integer,intent(out):: icon ! output. condition code >= 0  # of iterations
!!              -1--> unconvergence after maxite(=45) iterations
!!              -2--> x not in the range
!  real(8):: xa, xb, fa, fb, xt, ft
!
!  integer  n
!!
!  if(x < a .or. x > b) then
!     icon = -2
!  else
!     xa = a
!     fa = f(xa, param)
!     xb = b
!     fb = f(xb, param) 
!     if(fa*fb > 0.) then
!        icon = -3
!     else
!        icon = -1
!        xt = x
!        ft = f(xt, param)
!        do n = 0, 45
!           if( ft * fa >  0.) then
!              xa = xt
!              fa = ft
!           else
!              xb = xt
!              fb = ft
!           endif
!           if( eps > 0. ) then
!              if( abs( ft ) < eps ) then
!                 icon = n
!                 exit
!              endif
!           else
!              if(abs(xt) > 1.d0 ) then
!                 if(abs( (xa-xt) / xt ) < eps) then
!                    icon = n
!                    exit
!                 endif
!              else
!                 if(abs(xa-xt) < eps) then
!                    icon = n
!                    exit
!                 endif
!              endif
!           endif
!           xt = (xa + xb)/ 2
!           ft = f(xt, param)
!        enddo
!        ans = xt
!     endif
!  endif
!end subroutine kbinChopWP


                     

       
          


    
    
    
    
    
