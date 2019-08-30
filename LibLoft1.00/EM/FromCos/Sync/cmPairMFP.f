!       get mean free path of magnetic pair creation.
!       Integration of Eq. (16) of Baring
      real*8 function cmPairMFP(eg, xai)
      implicit none
#include "Zglobalc.h"

      real*8  eg   !  input.  gamma energy in GeV
      real*8  xai  !  input.  xai value = Eg/m *  H/Hc/2
      
      real*8  alpha,  cmPairInt1
      real*8 const
      parameter (const = 3. * 1.73205 * 3.1415 )
    
!           SyncConvR is in GeV/m
      alpha = cmPairInt1(xai)*  SyncConvR/ eg /const
      cmPairMFP = 1./alpha
      end
!      implicit none
!      real*8 xai, cmPairInt1, norm
!
!      xai = .01
!      do while (xai .lt. 1000.) 
!         if(xai .lt. 4.) then
!            norm = 7.497 * xai * exp(-4./3./xai)
!         else
!            norm = xai**0.6666
!         endif
!         write(*,*) sngl(xai), sngl(cmPairInt1(xai)/norm)
!         xai = xai *10.**.1
!      enddo
!      end
      real*8 function cmPairInt1(xai)
      implicit none
!       compute the integral part of Eq. (16) of Baring; Mon. Not. R.
!       astr.Soc. (1988) 235., i,e.
!      inte(0,1) dv  (9-v**2)/(1-v**2)K23(y).
!      The numerical integration has been done by inteSpecPair.f
!      and the result is normalized as follows
!
!   For xai < 4:  log(xai) vs ans/(7.497xai*exp(-4/3/xai)) is approximated by
!            a polynomial
!   for xai > 4:  log(xai) vs  ans/ xai**0.6666 is approximated by
!                 a polynomial.
!
      real*8 xai
!
      real*8 ans, lxai
!
      if(xai .lt. .01) then
         ans =7.497*xai*exp(-4./3./xai)
      elseif(xai .lt. 4.0) then
         lxai = log(xai)
          ans =(( -0.7211221E-03*lxai -0.1692442E-01)* lxai
     *       -0.9929865E-01 )*lxai + 0.8214562  
         ans = 7.497*xai* exp(-4./3./xai)* ans
      elseif(xai .lt. 5000.) then
         lxai = log(xai)

         ans = (((( 0.2210200E-03*lxai -0.8950944E-02)*lxai
     *  + 0.1461679 )*lxai -1.218747)*lxai +5.274929)*lxai
     *  + 0.1683762       
         ans = ans * xai**0.66666
      else
         ans = 9.8*xai**0.66666
      endif
      cmPairInt1 = ans
      end







