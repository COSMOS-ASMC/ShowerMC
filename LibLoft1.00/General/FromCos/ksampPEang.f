!c        test ksampPEang
!      implicit none
!      real*8  a, cost
!      integer i
!      a = 4.
!      write(0,*) 'Enter a >1 '
!      read(*,*)  a
!
!      do i = 1, 100000
!         call ksampPEang(a, cost)
!         write(*,*) cost
!      enddo
!      end

      subroutine ksampPEang(ain, cost)
      implicit none
!      This program samples cost=x from the following distribution
!         (1-x**2)/(a-x)**4 dx
!      where  a>1.  This is related to the cosine angle of the
!      ejected photo-electron at the photo-electric effect.
!     (angle is relative to the incident photon).  
!      This distribution is valid for few keV to *10 keV
!      photons.
!      a= ((alfa*Z)**2 +2Ee + Eg**2)/2Egsqrt(2Ee)
!      where  alfa: fin structure const.
!             Z:  effective atomic number
!             Ee: energy of electron /Me
!             Eg: //        photon /Me
!     
      real*8 ain !input  "a" value > 1
      real*8 cost ! output. sampled cos angle 

      real*8 x1, x2, x, eps, ans
      integer icon
      real*8 ksampPEangf
      external ksampPEangf
      real*8 norm, u, peakpos, am1

      real*8  a, a1, a2, a3, c1, un
      common /cksampPEang/a,  a1, a2, a3, c1, un
 
      a = ain
      c1 = (1.0-a**2)/3.0
      a1 = (a+1.0)
      a2 = a1*a1
      a3 = a2*a1
      am1=(a-1.0)
      norm = c1*(1.0/am1**3 - 1.0/a3) + a*(1./am1**2 -1.0/a2)
     *   - (1.0/am1 - 1.0/a1)
      call rndc(u)
      un = norm*u
!          peak pos. of the distribution 
      peakpos =(sqrt(a**2 + 8)-a)/2.0

      x = peakpos
      x1=-1.
      x2 =1.0
      eps = 1.d-4
      call kbinChop(ksampPEangf,
     *     x1, x2, x, eps, cost, icon)      
!      write(0,*) ' icon=',icon,  ' x=',ans
      end


      real*8 function ksampPEangf(x)
      implicit none
      real*8 x

      real*8  a, a1, a2, a3, c1, un
      common /cksampPEang/a, a1, a2, a3, c1, un

      real*8 ax

      ax = a-x
      ksampPEangf = c1*(a3-ax**3) +a*ax*a1*(a2-ax**2)
     *  -  ax**2*a2*(a1-ax) - un*ax**3*a3

      end


