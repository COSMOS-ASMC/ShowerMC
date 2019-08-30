      implicit none

      real*8  ct, g, gb, a, c, E0, m, t0, rho, x
      real*8 prob1, prob2, dedx, dedl, prob3, probax, prob
      real*8 cdecayProb0, cdecayProb1, cdecayProb2, xmax
      real*8 probx, proby
      c = 3.e8
      t0=2.3e-6
      m = 106.e-3
      E0=1.
      rho = 1.e-3 
      write(0, *)
     *  'E0=',E0, ' m=',m, ' rho=',rho, ' t0=',t0
      dedx = 2.e-3 /(1.e-3/1.e-4)  ! GeV/(g/cm^2)--> GeV/(kg/m^2)
      read(*,*) E0, m, rho, t0
      rho = rho*1.e3  ! kg/m^3
      ct = c*t0
      g = E0/m
      gb =  sqrt(g**2-1.)
      dedl =rho*dedx
      a = dedl/E0
      xmax = (1.-1./g)/a
      do x= 0., xmax,  xmax/200.
         if( g*(1.-x*a)-1.0 .le. 0.) goto 10
         prob1 = cdecayProb0(x, ct, g, gb, a)
!         if(g .gt. 10.) then
            prob2 = cdecayProb1(x, ct, g,  a)
!         else
!            prob2 = prob1
!         endif
         prob3 = exp(-x/(ct*gb))/gb
         prob = 1./sqrt( (g*(1-a*x))**2-1.)
         probx = exp(-x/(ct*gb))*prob
         probax = cdecayProb2(x, ct, g, gb, a)
         write(*,*) sngl(x), sngl(prob1), sngl(prob2), sngl(prob3),
     *          sngl(probax), sngl(prob), sngl(probx)
      enddo
 10   continue
      end
