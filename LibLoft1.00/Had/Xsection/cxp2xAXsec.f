!      implicit none
!     
!      real*8 a, xsxp, xsxa
!      integer ia, is
!
!      do  ia=4,99, 4
!         a=ia
!         do    is=1, 40
!             xsxp=1.d-1*10.**( (is-1)/10.)
!             call cxp2xAXsec(a, xsxp, xsxa)
!             write(*, *) sngl(xsxp), sngl(xsxa)
!         enddo
!         write(*,*)
!      enddo 
!      do  ia=100, 210, 10
!         a=ia
!         do    is=1, 40
!             xsxp=1.d-1*10.**( (is-1)/10.)
!             call cxp2xAXsec(a, xsxp, xsxa)
!             write(*, *) sngl(xsxp), sngl(xsxa)
!         enddo
!         write(*,*)
!      enddo 
!      end
!
! ****************************** cp_p2p_Axs  updated Feb. 3 2001
!                  (better approximation)
!          convert xp inelastic x-section into x-A inelastic xs.
        subroutine cxp2xAXsec(a, xsxp, xsxa)
!
!    a: real*8. input.  Mass # of the target
!   xsxp: real*8. input.  x-section of xp   (mb)
!   xsxa: real*8. output. x-section of xA   (mb)
!
        implicit none
        real*8  a, xsxp, xsxa
        real(8)::  cxAbyxpXsec
        xsxa = xsxp * cxAbyxpXsec(xsxp, a)
        end

