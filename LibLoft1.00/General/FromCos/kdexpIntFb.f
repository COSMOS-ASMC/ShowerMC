!c
!c      Takahashi-Mori's double exponential integration quadrature
!c
!        This is completely the same as kdexpIntF except for
!     the name; this can be used if two are used for double
!     integral to avoid recursive call.
! 
!     *******************************************************
      subroutine kdexpIntFb(func, a, b, eps, ans, error, icon)
!     *******************************************************
      implicit none
!
!          Numerical integration in a finite range by Takahashi-Mori's
!       double exponential integraton method.
!       The method gives accurate result even if there are singularities
!       at the integration limits. If there is a very sharp peak in the
!       midst of the integration region, the method will not give
!       a good result. 
!       
!
      real*8 func ! input. integrand function name with 1 argument
                  !          you have to declare 'external' 
                  !  Argument is real*8 xx(2).
                  !  xx(1) is x, 
                  !  Let aa=min(a, b), bb=max(a,b)
                  !  xx(2) = aa - x   if aa<= x < (aa+bb)/2
                  !        = bb - x   if bb> x >= (aa+bb)/2
                  !  So you should use xx(2) if at aa and/or bb
                  !  func is  singular; 
                  !       x(2) < 0 ==>  f(aa-x(2))
                  !       x(2) >=0 ==>  f(bb-x(2))
                   
      real*8 a   ! input. lower limit of the integration region
      real*8 b   ! input. upper //
                 !        b may be <=a.
      real*8 eps ! input. relative or absolute error of the
                 !        intetegration  you want to get.
                 !        if |ans| is > 1, used for relative error
                 !        else used for absolute error.
                 !        eps ~> 1.d-9  may be a good choice.
      real*8 ans ! output. approximate integration value.
      real*8 error ! output. estimated error (relative or absolute)
                   !         depending on |ans|.
      integer icon ! output. 0 --> ans is reliable.
                   !         1 --> ans may have a larger error than
                   !               your request.
                   !         2 --> input error so that ans is undef.
!
!
!

      integer halveNtime  ! if the integration accuracy is not enough
                          ! halve the equi-step of trapezoidal rule.
                          ! halving is tried upto halveNtime times.
      integer pointsInUnit !  see graph below
      integer blocks       ! we take 10 blocks
      integer totalpoints  ! (max number of  points)-1  where function 
                           ! is evaluated.
!      
!       integration bin; unit is
!  
!     1 2 3 ..  2**halveNtime = pointsInUnit    
!     | | | ... |
!      blocks of units are max number of points where the function is
!      evaluated.
!    
!   0 1 2 3 ..  pointsInUnit (note 0 at the top)
!   | | | | ... | | ..... | | | |    ..       | | |         | | | | ...| 
!                 1 2.... pointsInUnit   
!                                              
!
!          1           2            3                         bloks units
!
!    for the first integration use points marked O as below (bloks points)
!   O           O               O                                     O       c     
!    That is, if pointsInUnit = 32, in each block O is used
!    
!                                                                 
! 0  1 2                              16                               32
! |  | | | | | | | | | | | | |...        | | | | | | ...     | | | | | |
! O                                                                    O
!    halving this
!                                                              
! 0  1 2                              16                              32
! |  | | | | | | | | | | | | | | | | | | | ...     | | | | | | | | | | |
! O                                    x                               O
!    halving this
!
!
! 0  1 2 3 4 5 6 7 8 9 10  12 13  14  16  18  20  22  24  26  28  30  32
! |  | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
! O                +                   x               +               O
!    halving this
!
! 0 1 2 3 4 5 6 7 8 9 10  12  14  16  18  20  22  24  26  28  30  32
! | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | 
! O       #       +       #       x       #       +       #       O
!
!    halving this
! 0 1 2 3 4 5 6 7 8 9 10  12  14  16  18  20  22  24  26  28  30  32
! | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | 
! O   $   #   $   +   $   #   $   x   $   #   $   +   $   #   $   O
!
!    halving this is the final step
!
      parameter (halveNtime = 5, pointsInUnit = 32)
!                                  must be 2**halveNtime
      parameter ( blocks = 10, 
     *            totalpoints = blocks * pointsInUnit + 1 )
! 
!  mapping: y(0), y(1), y(32), y(33),  ... corresponds to  the
!  above  nodal points. y(0) is very close to -1 and y(totalpoints)
!  is very close to 1. (actually, they are -1, and 1 due to 
!  finite accuracy of  double precision. so 1+y and 1-y is 
!  tabulated specially.
! 

      real*8 y(0:totalpoints),  w(0:totalpoints)
      real*8 opy(0:totalpoints),omy(0:totalpoints)
      real*8 f(0:totalpoints)
      real*8 machmin, machmax  ! machine min, and max values
      real*8 halfpi
      real*8 tmax  ! max |t| where transformation 
                   ! y = tanh(pi/2 *sinh(t)) has no over/under flow.
                   !  |y| < 1
      real*8 h     ! minimum step of trapezoidal rule.
      real*8 t, c1, ans1, ans2, step,  f2, ytox, ytoxn, ytoxp
      real*8 temp, xa(2), expm, expp
      integer i, j, jstep, k

      logical first /.true./

      save first, y, w,  halfpi, tmax, h, opy, omy, temp

      ytox(k) = c1*(y(k) + 1) + a
      ytoxn(k) = -c1*opy(k)      ! -c1(1+y)
      ytoxp(k) =  c1*omy(k)      !  c1(1-y)

      if( first ) then
!              approx  machine min, max
         call kdmachmnmx(machmin, machmax)
         halfpi = asin(1.d0)  ! pi/2

!         tmax = log(log(machmax/1.d5)/halfpi )  
!            
!          The next choice is rather from trial and error
!                      log(log(1.d-75 ~ 1.d150).. ) is o.k
!         
            tmax = log(log(sqrt(machmin)/2)/(-2))

         h = 2*tmax/totalpoints   ! width of (-tmax, tmax)  is devided
                                  !  by totalpoints
!
!          compute nodal points and weight there.
!
         do  i = 0, totalpoints
            t = -tmax + i * h
            temp = halfpi * sinh(t)
            expm = exp(-temp)
            expp = exp( temp)

            y(i) = tanh( temp )

            opy(i) = 2*expp/(expp + expm) !   1+y
            omy(i) = 2*expm/(expp + expm) !   1-y
 
            w(i) = cosh(t) / cosh( halfpi*sinh(t) )**2
            
         enddo
         first = .false.
      endif

      if(a .eq. b) then
         ans = 0.d0
         icon = 0.
      else
         c1 = (b-a)/2.0d0
         ans1 = 0.

         jstep  =  pointsInUnit
         do i = 1, halveNtime
            step = jstep*h             
            ans2 = 0.
            do j = 0, totalpoints, jstep
               if( i.gt. 1 .and.
     *            mod( mod(j, pointsInUnit), jstep*2) .eq. 0) then
                  f2 = f(j)
               else
                  xa(1) = ytox(j)
                  if(y(j) .lt. 0. ) then
                     xa(2) = ytoxn(j)
                  else
                     xa(2) = ytoxp(j)
                  endif
                  f2 = func( xa ) * w(j)
                  f(j) = f2
               endif
               ans2 = ans2 + f2
            enddo
            ans2 = ans2 * step 
            if(i .gt. 1) then
               if(abs(ans2) .gt. 1.d0) then
                  error =abs( abs(ans1/ans2)-1.d0 )
                  if(error .le. eps) then
                     icon = 0 
                     goto 1000
                  endif
               else
                  error =  abs(ans2-ans1)
                  if(error .le. eps) then
                     icon = 0
                     goto 1000
                  endif
               endif
            endif
            ans1 = ans2
            jstep = jstep/2
         enddo
         icon = 1
 1000    continue
         ans = ans2 * halfpi *c1
      endif
      end
