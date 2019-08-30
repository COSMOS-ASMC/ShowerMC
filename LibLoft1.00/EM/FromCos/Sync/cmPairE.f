!c      tseting sampling of magnetic pair production spectrum.
!c
!      implicit none
!      integer i, nc, imax
!      imax = 20000
!      real*8 xai, ee
!      read(*,*) xai, imax
!      do i =1, imax 
!         call cmPairE(xai, ee, nc)
!         write(*,*)  sngl(ee), nc
!      enddo
!      end
      subroutine cmPairE(xai, ee, nc)
      implicit none
      real*8 xai  ! input. xai = Bsin/Bc *  Eg/m/2
      real*8 ee   ! output.  fractional electron energy >= 0.5
      integer nc  ! output. nubmer of retrial for rejection method.
!
      real*8 v, maxf, ve, xs5, xs, vx, a, b, beta, xse
      real*8 u, cmPairSpec
      integer i, j
      real*8 xaisave, sigma
      save xaisave,  maxf, vx
      save a, b, xs5
      data xaisave/0./
 
      if(xai .ne. xaisave ) then
         xs5 = cmPairSpec(xai, 0.5d0)
         if(xai .gt. 4.526) then
!              find peak which is not 0.5
            v = .5d0
            maxf = 0.
            do while (v .lt. .999d0)
               xs = cmPairSpec(xai, v)
               if(xs .gt. maxf) then
                  maxf = xs
                  vx = v
               endif
               v = v + .005d0
            enddo
            a = (maxf - xs5)/(vx - 0.5) 
            b = -0.5 * a + xs5
         else
            vx = 0.5d0
            maxf = xs5
         endif
         xaisave = xai
      endif
      if(vx .ne. .5d0) then
         do i = 1, 1000
            call rndc(u)
            beta = (3./8.*a + b/2)* u - (a/2 + b)
            ve = (- b + sqrt(b**2 - 2*a*beta))/a
!            write(*, *) ' beta, ve', beta, ve, ' b^2-2beta',
!     *       b**2 - 2*beta
            call rndc(u)
            xse = cmPairSpec(xai, ve)
!            write(*,*) ' xse=',xse, ' a*ve+b=', a*ve+b
            if(u .lt. xse/(a*ve + b)) goto 10
         enddo         
         stop 9875
 10      continue
         ee = ve
      elseif(xai .lt. .005) then
         ee = .5
         nc = 1
      elseif(xai .lt. 1.) then
!          envelop is close to Gausssian of
!          sigma= .50/root(2)*xai**.53
         sigma= .50/1.4142*xai**0.53
         do j = 1, 1000
            do i = j, 1000
               call kgauss(0.5d0, sigma, ve)
               if(ve .gt. 0. .and. ve .lt. 1.) goto 5
            enddo
 5          continue
            xse =  cmPairSpec(xai, ve)
            call rndc(u)
            if(u .lt. 
     *      xse /( ( exp(-((ve - 0.5)/sigma)**2/2) * maxf)))
     *          goto 8
         enddo
 8       continue
         ee = ve
         if(ee .lt. 0.5) ee = 1. - ee
      else
         do i =1, 1000 
            call rndc(u)
            ve = u/2.+ 0.5
            xse = cmPairSpec(xai, ve)
            call rndc(u)
            if(u .lt. xse/maxf) goto 20
         enddo
         stop  1234
 20      continue
         ee =ve
      endif
      nc = i
      end
