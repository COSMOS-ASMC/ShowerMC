!     ***************************************
!      cmigdG:  G(s)
!      cmigdPsi: Psi(s)
!      cmigdGzai: Gsai(s)
!      cmigdSforBrem:  s for Bremsstrahlung
!      cmigdSforPair:  s for pair creation 
!
!     ***************************************
!     --------------------------------
!            to test cmigdSforBrem/cmigdSforPair
!
!      real*8 cmigdSforPair, s, v, e, rho, x
!
!
!      e = 1.e19/1.e9              !  GeV
!      
!      x = 1.e19/1.e9 * 1.e-2      !  GeV kg/m3
!      do while (x .lt. 1.e22/1.e9 * 1.2)
!         v = 1.e-4
!         do while (v .lt. 1.)
!            rho = x / e               !in assume kg/m3
!            s = cmigdSforPair(e, rho, v)
!            v = v*10.**0.01
!         enddo
!         write(*, *)
!         x = x * 10.**0.5
!      enddo
!      end
!
!     --------------------------------
!         to test cmigdPsi:
!      real*8 cmigdPsi, psimig, s
!      do s = 0., 1., 0.01
!         write(*,*) cmigdPsi(s), psimig(s, 1.d-5)
!      enddo
!      end
!     --------------------------------
!         to test cmigdG:
!      real*8 cmigdG, gmigdl, s
!      do s = 0., 1., 0.01
!         write(*,*) cmigdG(s), gmigdl(s, 1.d-5)
!      enddo
!      end

      
!     ******************************************      
!
!          G(s): better than 1% accuracy.
!
!     ******************************************      
      real*8 function cmigdG(s)
      implicit none
      real*8 s  ! input 
!
!        G(s)/( 12pi*s^2/(1+12pi*s^2)) is approximated by polinomials
!
      real*8 xxx1(5)   ! pol coeff. for s < 0.12  G(s) error is 3 x 10^-5 
      real*8 xxx2(7)   ! pol coeff. for s > 0.11       error is 1 x 10^-3

      real*8 pi, const, f
      integer i
      parameter (pi= 3.141592654d0, const = 12.d0 * pi)

          data ( xxx1(i), i=  1,   5)/
     1     0.99995995    ,  -6.2614068    ,   57.629644    ,
     2   -242.27058    ,   450.83728                                    
     * /   
  
          data ( xxx2(i), i=  1,   7)/
     1         0.78283668    ,  -1.3317861    ,   12.312840    ,
     2   -33.848229    ,   44.733673    ,  -29.201775    ,
     3         7.5561750                                                     
     * /   


       if(s .lt. 0.12) then
          f  = 0.       
          do i = 5, 2, -1
             f = (f + xxx1(i))* s
          enddo
          f = f + xxx1(1)
          cmigdG = f * const*s**2/(1. + const*s**2)
       elseif(s .lt. 1.1) then
          f  = 0.       
          do i = 7, 2, -1
             f = (f + xxx2(i))* s
          enddo
          f = f + xxx2(1)
          cmigdG = f * const*s**2/(1. + const*s**2)
       else
          cmigdG = 1. - 0.022/s**4
       endif

       end
!      *********************************************
!
!       psi(s): since psi--> 6s for s -->0, we approximate
!               psi(s)/ ( 6s/(1+6s) ) by a polynomial.
!               
!      
!      *********************************************
       real*8 function cmigdPsi(s)
       implicit none
       real*8 s  ! input

       real*8 f
       real*8 xxx1(4)   ! coef. of pol. s<0.3  better than 0.1 % for Psi
       real*8 xxx2(5)   ! //            s> 0.3  //
       integer i

          data ( xxx1(i), i=  1,   4)/
     1    1.00    ,   2.6978063    ,  -9.4242869    ,
     2    11.468298                                                     
     * /   

          data ( xxx2(i), i=  1,   5)/
     1    1.2095058    ,  0.57895055    ,  -1.6531094    ,
     2    1.4846143    , -0.46392960                                    
     * /   

       if(s .lt. 0.3) then
          f = 0.
          do i = 4, 2, -1
             f = (f + xxx1(i)) * s
          enddo
          f = f + xxx1(1)
          cmigdPsi = f * 6.d0 *s /(1.d0 + 6.d0*s)
       elseif(s .lt. 1.2) then
          f = 0.
          do i = 5, 2, -1
             f = (f + xxx2(i)) * s
          enddo
          f = f + xxx2(1)
          cmigdPsi = f * 6.d0 *s /(1.d0 + 6.d0*s)
       else
          cmigdPsi = 1. - 0.012/s**4
       endif
       end
!      ***************************************
!    
!        cmigdGzai; this is for air
! 
!      **************************************
!
       real*8 function cmigdGzai(s)
       implicit none
       real*8 s
       real*8 s1/1.1185e-4/  ! (Z^0.333/183)^2 ; Z=7.25
       real*8 logs1/-9.0983/  ! ln(s1)
!
       if( s .gt. 1.) then
          cmigdGzai = 1.
       elseif(s .gt. s1) then
          cmigdGzai = log(s)/logs1 + 1.
       else
          cmigdGzai =2.
       endif
       end
!      **************************************
!     
!      cmigdSforBrem; get s for brems defined recursively
!
!      **************************************
      real*8 function cmigdSforBrem(ee, rho, v)
      implicit none
#include "Zmass.h"
#include "Zelemagp.h"
!
      real*8 ee     ! input.  Electron energy in GeV.
      real*8 rho    ! input.  air density in kg/m3
      real*8 v      ! input.  Fractional energy of gamma, Eg/Ee
!
      real*8 x, ss, s2, cmigdGzai, temp
      integer i
      
      if(v .eq. 1.) then
         cmigdSforBrem = 1.
      else
         x = ee * rho   ! GeV kg/m3
         ss = 1.
         temp =  masele*  X0/x *100. * v/(1.-v)   !100: X0/rho in cm
         do i = 1, 10
!                normally 3 to 4 iteration is ok
            s2 = 1.37e3 * sqrt(temp/cmigdGzai(ss))
            if(abs(ss/s2 -1.d0) .lt. 3.e-3) goto 200
            ss = s2
         enddo
         s2 = 1. 
 200     continue
         cmigdSforBrem = s2
!       
!        write(*, *) i, sngl(s2), sngl(v), sngl(ee*rho)
!       
      endif
      end
!      **************************************
!     
!      cmigdSforPair; get s for pair defined recursively
!
!      **************************************
      real*8 function cmigdSforPair(eg, rho, v)
      implicit none
#include "Zmass.h"
#include "Zelemagp.h"
!
      real*8 eg     ! input.  Gamma energy in GeV.
      real*8 rho    ! input.  air density in kg/m3
      real*8 v      ! input.  Fractional energy of e-/e+ Ee/Eg
!
      real*8 x, ss, s2, cmigdGzai, temp
      integer i
      
      if(v .eq. 1.) then
         cmigdSforPair = 1.
      else
         x = eg * rho   ! GeV kg/m3
         ss = 1.
         temp =  masele*  X0/x *100./v/(1.-v)   !100: X0/rho in cm
         do i = 1, 10
!                normally 3 to 4 iteration is ok
            s2 = 1.37e3 * sqrt(temp/cmigdGzai(ss))
            if(abs(ss/s2 -1.d0) .lt. 3.e-3) goto 200
            ss = s2
         enddo
         s2 = 1. 
 200     continue
         cmigdSforPair = s2
!       
!         write(*, *) i, sngl(s2), sngl(v), sngl(eg*rho)
!       
      endif
      end



