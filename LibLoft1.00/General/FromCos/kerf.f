!         test kerf
!
!      implicit none
!      real*8 x, kerf
!      do x = 0.d0, 10.d0, 0.5d0
!         write(*, *) ' x=', x, ' erf=',kerf(x)
!      enddo
!      end
!     *********************************
      real*8 function kerf(x)
!     *********************************
!     Error function:  integral from 0 to x
!                      of 2/sqrt(pi) * exp(-t**2)dt
!              
!         Adapted from Mori's code in Maruzen book.
!
!    *** note ***
!
!      kerf(inf) =  1
! 
!      I(0:inf){exp(-(ax**2 + 2bx + c)} =
!
!             sqrt(pi/a)/2 exp( (b**2 - ac)/a ) kerf(b/sqrt(a))
!
!      I(0:inf){exp(-ax**2 - b/(x**2))} =
!     
!             sqrt(pi/2)/2 exp(-2sqrt(ab))
!
!
!      I( :x){exp(-(ax**2 + 2bx + c))}dx =
! 
!             sqrt(pi/a)/2 exp( (b**2-ac)/a ) kerf(sqrt(a)x + b/sqrt(a))
!             (a != 0)
!
!      let z(x)=1/sqrt(2pi) * exp(-x**2/2)
!          p(x)=integral z(x) from -inf to x
!          q(x)= //                 x  to inf
!          a(x)= //                 -x to x
!      then p+q=1;  p(-x)=q(x);  a=2p-1.
!      p( (x-m)/s ) =1/sqrt(2pi)/s * integral
!             exp(- ((x-m)/s)**2/2) from -inf to x
!      z( (m+x)/s ) = z( (m-x)/s )
!      kerrorf(x)= 2p(sqrt(2)x) -1
!      p(x)=  ( erf(x/sqrt(2)) + 1 )/2

      implicit none
      real*8 x
!
      integer nm, nx, na

      parameter (nm = 5, nx = 13, na = 5)

      real*8 cm(0:nm), cx(nx), cq(nx), ca(0:na)
      
      real*8 sqrtpi, norm,  invpi, fourpi, pi

      parameter (sqrtpi = 1.772453850905516d0)
      parameter (norm = 2.d0 / sqrtpi, pi = 3.141592653589793d0)
      parameter (invpi = 1.0d0 / pi, fourpi = 4.d0 *pi)
      real*8 y, xv, v
      integer i

      data cm /
     *   0.1000000000000000d+01,
     *  -0.3333333333333333d+00,
     *   0.1000000000000000d+00,
     *  -0.2380952380952381d-01,
     *   0.4629629629629630d-02,
     *  -0.7575757575757575d-03 /

      data cx /
     *   0.7788007830714048d+00,
     *   0.3678794411714423d+00,
     *   0.1053992245618643d+00,
     *   0.1831563888873418d-01,
     *   0.1930454136227709d-02,
     *   0.1234098040866796d-03,
     *   0.4785117392129009d-05,
     *   0.1125351747192591d-06,
     *   0.1605228055185612d-08,
     *   0.1388794386496402d-10,
     *   0.7287724095819692d-13,
     *   0.2319522830243569d-15,
     *   0.4477732441718302d-18 /

      data cq /
     *   0.2500000000000000d+00,
     *   0.1000000000000000d+01,
     *   0.2250000000000000d+01,
     *   0.4000000000000000d+01,
     *   0.6250000000000000d+01,
     *   0.9000000000000000d+01,
     *   0.1225000000000000d+02,
     *   0.1600000000000000d+02,
     *   0.2025000000000000d+02,
     *   0.2500000000000000d+02,
     *   0.3025000000000000d+02,
     *   0.3600000000000000d+02,
     *   0.4225000000000000d+02 /

      data ca /
     *   0.1000000000000000d+01,
     *  -0.1000000000000000d+01,
     *   0.3000000000000000d+01,
     *  -0.1500000000000000d+02,
     *   0.1050000000000000d+03,
     *  -0.9449999999999999d+03 /
!  --------------------------------------------

      xv = abs(x)
      if (xv .le. 0.1d0) then
         y = xv**2
         v = cm(nm)
         do i = nm - 1, 0, -1
            v = cm(i) + y * v
         enddo
         kerf = norm * xv * v
      elseif (xv .le. 8.0d0) then
         y = xv**2
         v = 1.d0 / (2 * y)
         do i = 1, nx
            v = v + cx(i) / (cq(i) + y)
         enddo
         v = invpi * xv * exp(-y) * v
         if (xv .lt. 6.0d0) then
            v = v - 2 / (exp(fourpi * xv) - 1)
         end if
         kerf = 1.d0 - v
      else
!         y = 2 * xv**2
!         v = ca(na)
!         do i = na - 1, 0, -1
!            v = ca(i) + y * v
!         enddo
!         v = exp(-xv**2) / (sqrtpi * xv) * v
!         kerf = 1.d0 - v
         kerf = 1.d0
      endif

      if (x .lt. 0.d0) then
        kerf = -kerf
      endif
      end
