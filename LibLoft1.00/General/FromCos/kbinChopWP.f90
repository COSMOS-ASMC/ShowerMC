subroutine kbinChopWP(f, param,  a, b, x, eps, ans, icon)
  implicit none
!
  real(8):: f  ! input. function name. to be used as f(x, param)
  !   f(x, param) = 0 is solved.
  real(8),intent(in):: param(*)
  real(8),intent(in):: a ! input. lower bound of solution
  real(8),intent(in):: b ! input. upper bound of solution
  real(8),intent(in):: x   ! input. initial guess of  solution.  
  real(8),intent(in):: eps  ! input. relative error of solution.
                     ! > 0  |func| < eps
                     ! < 0   |x1-x2|/ |x2| < eps (|x2| >= 1)
                     !       |x1-x2| < eps  (|x2|< )
  real(8),intent(out):: ans  ! output. obtained solution
  integer,intent(out):: icon ! output. condition code >= 0  # of iterations
!              -1--> unconvergence after maxite(=45) iterations. may be ans
!                    can be used.            
!              -2--> x not in the range
!              -3--> f(a)*f(b)> 0.
  real(8):: xa, xb, fa, fb, xt, ft

  integer  n
!
  if(x < a .or. x > b) then
     icon = -2
  else
     xa = a
     fa = f(xa, param)
     xb = b
     fb = f(xb, param) 
     if(fa*fb > 0.) then
        icon = -3
     else
        icon = -1
        xt = x
        ft = f(xt, param)
        do n = 0, 45
           if( ft * fa >  0.) then
              xa = xt
              fa = ft
           else
              xb = xt
              fb = ft
           endif
           if( eps > 0. ) then
              if( abs( ft ) < eps ) then
                 icon = n
                 exit
              endif
           else
              if(abs(xt) > 1.d0 ) then
                 if(abs( (xa-xt) / xt ) < eps) then
                    icon = n
                    exit
                 endif
              else
                 if(abs(xa-xt) < eps) then
                    icon = n
                    exit
                 endif
              endif
           endif
           xt = (xa + xb)/ 2
           ft = f(xt, param)
        enddo
        ans = xt
     endif
  endif
end subroutine kbinChopWP


                     

       
          


    
    
    
    
    
