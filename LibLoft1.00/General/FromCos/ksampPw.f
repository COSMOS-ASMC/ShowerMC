!      ***************************************************************
!      ksampPw:  samples from a spectrum consisting of a number of 
!                power functions.
!      ksampPwX: get crossing points of adjacent lines given by
!                power functions
!      ksampPwMrg: Multiply two spectrum, each consistng of a number of
!                  power functions and make a new spectrum consisting of
!                  a number of power functions.
!      
!      ksampPwQ:  Inquire the internal table in ksampPw.  
!      ksampPw2:  The same as ksampPw.  Needs information from ksampPwQ
!                 so that the initialization calculaion in ksampPw can
!                 be skipped.
!      **************************************************************
!
!      The next test program shows how to use ksampPwX, ksampPwMrg and
!      ksampPw.
!---------------------------------------------------------------------- 
!      program testsampPowers
!      implicit none
!
!      integer n1, n2, n,  i, nx
!      parameter (n1 = 3, n2=3, n=10) 
!      real*8  c(n), p(n), xp(n+1)
!      real*8  fx, up, x
!      real*8  x01(n1),  x02(n2), c1(n1), c2(n2)
!      real*8  p1(n1), p2(n2)
!      real*8  xp1(n1+1), coef1(n1)
!      real*8  xp2(n2+1), coef2(n2)
!
!      data c1/2.14, 0.4d0, 0.05d0/
!      data p1/-0.3333333333d0, 0.7d0, 4./
!      data x01/1., 1.3d0, 4.0d0/
!
!      data p2/1., 2.0, 3./
!
!      up = 0.1
!
!      c2(1) = 1.
!      c2(2) = up/3.
!      c2(3) = 1./9/up/up
!
!      x02(1) = 1.
!      x02(2) = 1./3./up
!      x02(3) = 1.
! 
!c      two spectra are defined. One is  by c1,p1,x01  and
!c      the other by c2, p2, x02.
!c      Get the crossing points of segments for the first specrum 
!c      and store it in xp1(2), ...xp1(n1) and coefficients
!c      in coef1(1), .. coef1(n1) so that  coef1 * X** (- p1) dX
!c      be the spectrum.
!
!      call ksampPwX( c1, p1, x01, n1, xp1(2), coef1 )
!c        define the left and right boundary of the spectrum.
!      xp1(1) = 0.
!      xp1(n1+1) = 15.
!      
!      write(*, *) ' xp1=', xp1
!      write(*, *) ' coef1=', coef1
!
!c      for the second spectrum.
!
!      call ksampPwX(c2, p2, x02, n2, xp2(2), coef2)
!      xp2(1) = 0.
!      xp2(n2+1) = 1.e10
!
!      write(*, *) ' xp2=', xp2
!      write(*, *) ' coef2=', coef2
!c       Merge the two spectrum and obtaine new one in
!c       c, p, xp, with nx segments
!
!      call ksampPwMrg(coef1, p1, xp1, n1, 
!     *                coef2, p2, xp2, n2,
!     *                c,     p,   xp, nx)
!      write(*,*) ' merged n c=', nx, c
!      write(*, *) ' p=', p
!      write(*,*) ' xp=', xp
!    ----------------      
!
!      do i = 1, 10000
!c         samples x and compute the function value at x
!c         (fx).
!         call ksampPw(i, c, p, xp, nx,  x, fx)
!         write(*, *) sngl(x)
!      enddo
!      end
!    ***********************************************************
!    Samples a random number from a spectrum composed of a 
!    number of segments,  each of which has a function form
!    like:  cX**(-p)dX
!    ***********************************************************
!    
      subroutine  ksampPw(ini, coef, power, node, n,  x, fx)      
      implicit none
!
      integer ini  ! input.  if coef, power, node, n are
                   !     different from the previous call, give 1.
                   !     If all of them are the same as the
                   !     previous call, give non 1
      integer n    ! input.  number of segments
      real*8 coef(n)  ! input.  c of cX**(-p)dX in each segment.
      real*8 power(n) ! input.  p of cX**(-p)dX in each segment.
      real*8 node(n+1)  ! input.  node(1) and node(2) are the left and
                      !         right boundaries of segment 1.,
                      !         node(2) and  node(3) are for segment 2,etc
      real*8 x        ! output. A sampled x. If the input condition 
                      !        is ng, x = 0 will result.
      real*8 fx       ! output.  function value  at  x, i.e., cx**-p
! 
      integer nmx     !  max number of nodes
      parameter ( nmx = 20 )

      real*8 inte(nmx)  ! integral value of each segment region
      real*8 ci(nmx)    ! normalized cummulative integral values of inte. 
      integer nsave     ! for saving n

      save inte,  ci, nsave
      real*8 u, sum, temp
      integer i, j

      real*8  cum(*), intgv(*)
      integer np

!
      if(ini .eq. 1) then
         if(nmx .lt. n) then
            x = 0.
         else
!             get integral value of each region and sum of them
            sum = 0.
!             c/(1-p) x**(1-p) from x1 to x2
            do i = 1, n
               if(power(i) .ne. 1.d0) then
                 inte(i) = coef(i)/(1.0-power(i))*
     *           (node(i+1)**(1.-power(i)) - node(i)**(1.-power(i)))
               else
                  inte(i) = coef(i) * log(node(i+1)/node(i))
               endif
               sum = sum + inte(i)
            enddo
!             normalize and make cummulative data
            do i = 1, n
               ci(i) = inte(i)/sum
            enddo
            do i = 2, n
               ci(i) = ci(i) + ci(i-1)
            enddo
            ci(n) = 1.0   ! for safety
         endif
         nsave = n
      endif
!     --------------------------------------------------
!           assume inte and node are ready
!          choose which segment, j
      call rndc(u)
      do i = 1, n
         if(u .le. ci(i)) then
            j = i
            goto 100
         endif
      enddo
!       never come here
!
 100  continue
      

      call rndc(u)
      if(power(j) .ne. 1.d0) then
         temp = node(j+1)**(1.-power(j))  - 
     *       inte(j)*u/( coef(j)/(1. - power(j)) )       ! = x **(1-power(j))

         x = temp**(1./(1.-power(j)))
      else
         x = (node(j+1)/node(j))**u  * node(j)
      endif
      fx = coef(j) * x **(-power(j))
 1000 continue
      return
!  
!
!     ********************************
      entry ksampPwQ(cum, intgv, np)
!     *******************************
!      inquire the current inner variables for the 
!      present spectrum.  
!

      np = nsave

      do i = 1, np
         cum(i) = ci(i)
         intgv(i) = inte(i) 
      enddo
      end
!
!     *********************************
      subroutine ksampPw2(coef, power, xp, cum, intgv, np, x, fx)
!
!      This interface is used to skip the initial calculations
!      by using output from ksampPwQ
!
      implicit none
!

      integer np      ! input.  number of segments
      real*8 coef(np)  ! input. c of X**(-p)dX in each segment.
      real*8 power(np) ! input.  p of X**(-p)dX in each segment.
      real*8 xp(np+1)  ! input.  output from ksampPwQ
      real*8 cum(np)   ! input.  //
      real*8 intgv(np) ! //
    
      real*8 x        ! output. sampled x. 
      real*8 fx       ! output. funcion value at x. 

      real*8  u, temp
      integer i, j
!          choose which segment, j
      call rndc(u)
      do i = 1, np
         if(u .le. cum(i)) then
            j = i
            goto 100
         endif
      enddo
!
 100  continue

      call rndc(u)
      temp = xp(j+1)**(1.-power(j))  - 
     *       intgv(j)*u/( coef(j)/(1. - power(j)) )       ! = x **(1-power(j))

      x = temp**(1./(1.-power(j)))
      fx = coef(j) * x**(-power(j))
      end
!    ***********************************************************
!    get xssing points of adjacent segments composed of power spectra
!
      subroutine  ksampPwX(c, p, x0, n, xp, c2)
      implicit none
!
      integer n    ! input.  number of segments
      real*8 c(n)  ! input.  c of c(X/x0)**(-p) in each segment.
      real*8 p(n) ! input.    see above
      real*8 x0(n)  ! input.  see above
      real*8 xp(n-1)  ! output.  crossing point of adjacent segments
      real*8 c2(n)  ! output. c2 = c/x0**p. That is, c2*x**(-p) becomes
                    !           the spectrum
! 
      integer i

!      get coeff.      
      do i = 1, n
         c2(i) =  c(i)* x0(i)**p(i)
      enddo
!     get crossing point of i and i+1th segmment
!              c_i X**(-p_i) = c_i+1 X**(-p_i+1)

      do i = 1, n-1
          xp(i) =
     *    (c2(i)/c2(i+1) )**(1./(p(i) - p(i+1)))
      enddo
      
      end
!     *************************************************************
!     *
!     * get new power spectra obtained by multipling two
!     *     power specta,  each consisting of a number of 
!     *     segments

      subroutine ksampPwMrg(c1, p1, x1, n1, c2, p2, x2, n2,
     *  c, p, x, n)
!
      implicit none
      integer n1   ! input.         number of segments of 1 st function
                   !                c1* x**(-p1) dx 
      real*8  c1(n1)  ! input.     
      real*8  p1(n1)  ! input.
      real*8  x1(n1+1)  ! input    ! nodal point of segments of 1 st func.
      integer n2   ! input.        ! the same as the 2nd func.
      real*8  c2(n2)  ! input.   
      real*8  p2(n2)  ! input.
      real*8  x2(n2+1)  ! input
      integer n   !  output        !  those for f1 x f2
      real*8  c(*)  ! output.
      real*8  p(*)  ! output
      real*8  x(*)  ! output   (dim >= nx+1 )

!

      integer i, j
!
      i = 1
      j= 1
      n = 0

      x(1) = max(x1(1), x2(1))
      do while (i .le. n1 .and. j .le. n2)
         if(x1(i+1) .lt. x2(j+1) ) then
            n = n + 1
            x(n+1) = x1(i+1)
            c(n) = c1(i) * c2(j)
            p(n)  = p1(i) + p2(j)
            i = i + 1

         elseif(x1(i+1) .eq. x2(j+1)) then
            n = n + 1
            x(n+1) = x1(i+1)
            c(n) = c1(i) * c2(j)
            p(n)  = p1(i) + p2(j)
            i = i + 1
            j = j + 1
         else
            n = n + 1
            x(n+1) = x2(j+1)
            c(n) = c1(i) * c2(j)
            p(n)  = p1(i) + p2(j)
            j = j + 1
         endif
      enddo
      end


