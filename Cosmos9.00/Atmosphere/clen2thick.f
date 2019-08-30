!     clen2thick
!
!     Thickness of air corresponding to a given length along
!     a given direction from a given height.
!
!        test program
!     include 'catmosutil.f'
!     include 'cstdatmos0.f'
!     include '../../KKlib/k16pGaussLeg.f'
!
!     program testclen2thick
!     implicit none
!     real*8 z, cosz, len, thick, thicka
!
!     real*8 clen2thick
!
!     z = 6000.
!     cosz = 1.
!     len = 3000.
!     write(*, *) " enter z, cosz, len=",z, cosz, len
!     read(*,*) z, cosz, len
!     thick = clen2thick(z, cosz, len)
!     write(*, *) " thick = ", thick
!     end
!     ========================================
      real*8 function clen2thick(z, cosz, leng)
!
!     z: real*8. input.  vertical height in m.
!  cosz: real*8. input.  cos of zenith angle at z.
!   leng: real*8. input.  length along cosz direction.
!   function value.      thickness of air in kg/m2.
!
      use modAtmosDef
      implicit none
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
! #include  "Zatmos.h"

      real*8 z, cosz, leng, clen2thickEx, clen2thickAp
      real*8  clen2thickTA

      if(UseTbl .and. z .lt. Htop) then
         clen2thick = clen2thickTA(z, leng)
      elseif(ExactThick) then
         clen2thick = clen2thickEx(z, cosz, leng, 10)
      elseif(abs(cosz) .gt. 0.6) then
         clen2thick = clen2thickAp(z, cosz, leng)
      else
         clen2thick = clen2thickEx(z, cosz, leng, 10)
      endif
      end
!     cLen2thickAp
!
!     Thickness of air corresponding to a given length along
!     a given direction from a given height.
!
!        test program
!     include 'catmosutil.f'
!     include 'cstdatmos0.f'
!     include '../../KKlib/k16pGaussLeg.f'
!     include 'clen2thickEx.f'
!
!     program testclen2thickAp
!     implicit none
!     real*8 z, cosz, leng, thick, thicka
!
!     real*8 clen2thickAp, clen2thickEx
!
!     z = 6000.
!     cosz = 1.
!     leng = 3000.
!     write(*, *) " enter z, cosz, leng=",z, cosz, leng
!     read(*,*) z, cosz, leng
!     thick = clen2thickEx(z, cosz, leng, 10)
!     thicka = clen2thickAp(z, cosz, leng)
!     write(*, *) ' thick exact =',thick, ' apprxo = ', thicka
!     end
!     ========================================
      real*8 function clen2thickAp(z, cosz, leng)
!
!     z: real*8. input.  vertical height in m.
!  cosz: real*8. input.  cos of zenith angle at z.
!   leng: real*8. input.  length along cosz direction.
!   function value.      thickness of air in kg/m2.
!
      use modAtmosDef
      implicit none
!----      include 'Zearth.h'
!  #include  "Zearth.h"

      real*8 z, cosz, leng
!

      real*8 cosm, rs, re, cvh2thick, cnewcos, ze, t1, t2,
     *  cnewh



      rs = z + Eradius
      re = cnewh(rs, cosz, leng)
      ze = re - Eradius
      cosm = cnewcos(rs, cosz, leng/2)
      t1 = cvh2thick(z)
      t2 = cvh2thick(ze)
      clen2thickAp = (t2 - t1)/cosm
      end
!     Exact thickness of air corresponding to a gine length along
!     a given direction from a given height.
!     This uses numerical integration.
!
!        test program
!      include '../../KKlib/k16pGaussLeg.f'
!      include 'catmosutil.f'
!      include 'cstdatmos0.f'
!
!      program testclen2thickEx
!      implicit none
!      real*8 z, cosz, leng, thick
!
!      real*8 clen2thickEx, cvh2thick
!
!      z = 6000.
!      cosz = 1.
!      leng = 3000.
!      thick = clen2thickEx(z, cosz, leng, 10)
!      write(*, *) ' cos=',cosz,' leng =',leng, ' thick=',thick
!      thick = cvh2thick(z-leng*cosz) - cvh2thick(z)
!      write(*, *) ' thick=', thick
!      end
!     ========================================
      real*8 function clen2thickEx(z, cosz, leng, n)
!
!     z: real*8. input.  vertical height in m.
!  cosz: real*8. input.  cos of zenith angle at z.
!   leng: real*8. input.  length along cosz direction.
!     n:   integer.input.  how many points for Gauss-Legendre integ.
!   function value.      thickness of air in kg/m2.
!
      use modAtmosDef
      implicit none
!----      include 'Zearth.h'
!  #include  "Zearth.h"

      real*8 z, cosz, leng
      integer n
!
      real*8  seg/1000./, ans, a, b, inte
!
      external catmosrho
      real*8 catmosrho

      common /ccatmosrho/ coss, rs
      real*8 coss, rs
!/////
      real(8),parameter::eps=1.d-5
      real(8):: error
      integer icon
!///////

      inte =  0.
      b  =0.
      coss = cosz
      rs = z + Eradius
      do while (.true.)
         a = b
         b = min(leng, a + seg)
!         call kdexpIntF(catmosrho, a, b, eps, ans, error, icon)
         call k16pGaussLeg(catmosrho, a, b, n, ans)
         inte = inte + ans
!         write(*, *) 'a =', a, ' b=',b, ' ans=',ans
         if(b .eq. leng) goto 10
      enddo
 10   continue
      clen2thickEx = inte
      end
