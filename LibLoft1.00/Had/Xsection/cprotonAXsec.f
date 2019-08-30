!         p-A, pi-A, k-A  inelastic cross-seciton  etc.
!
!      real*8 a, e, xs, ex
!      a = 14.5
!      e = 2.e0
!      do while (e .lt. 1.e22/1.e9)
!         ex = e - 0.938
!         call cprotonAXsec(a, ex, xs)
!         write(*, *)  sngl(e),sngl(xs)
!         e = e * 10.**0.2
!      enddo
!      end
      subroutine cprotonAXsec(a, e, xs)
!      a: real*8. input. target mass number
!      e: real*8. input. kinetic energy of the ptcl(in GeV)
!     xs: real*8. output. xsection in mb
      implicit none
      real*8 a, e, xs
!
      real*8 xsxp
!
      call cppXsec(e, xsxp)
      call cxp2xAXsec(a, xsxp, xs)
      end
!     --------------------------------
      subroutine cpbarAXsec(a, e, xs)
      implicit none
      real*8 a, e, xs
!
      real*8 xsxp
!
      call cpbarpXsec(e, xsxp)
      call cxp2xAXsec(a, xsxp, xs)
      end
!     ------------------------------------------
      subroutine cpiMinusAXsec(a, e, xs)
      implicit none
      real*8 a, e, xs
!
      real*8 xsxp
      
      call cpiMinuspXsec(e, xsxp)
      call cxp2xAXsec(a, xsxp, xs)
      end
!     ------------------------------------------
      subroutine cpiPlusAXsec(a, e, xs)
      implicit none
      real*8 a, e, xs
!
      real*8 xsxp
      
      call cpiPluspXsec(e, xsxp)
      call cxp2xAXsec(a, xsxp, xs)
      end
!     ------------------------------------------
      subroutine ckMinusAXsec(a, e, xs)
      implicit none
      real*8 a, e, xs
!
      real*8 xsxp
      
      call ckMinuspXsec(e, xsxp)
      call cxp2xAXsec(a, xsxp, xs)
      end
      subroutine ckPlusAXsec(a, e, xs)
      implicit none
      real*8 a, e, xs
!
      real*8 xsxp
      
      call ckPluspXsec(e, xsxp)
      call cxp2xAXsec(a, xsxp, xs)
      end



