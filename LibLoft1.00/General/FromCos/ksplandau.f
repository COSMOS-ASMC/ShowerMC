!       test ksplandau
!      integer i
!      real*8  b, c, x
!
!      b =1.
!      c =0.20
!      do i= 1, 100000
!         call ksplandau(b, c, x)
!         write(*,*) sngl(x)
!      enddo
!      end

!       ksplandau.f
!  samples a random variable x following the
!  psudo-Landau energy loss distribution 
!  which is
!       exp(-1/2(y +exp(-y)) dy
!  where  y=(x-b)/c.
!
! The peak of the distribution is at x=b.
!  <y> = 1.2696   FWHM(y) = 3.56
!  <x> = 1.2696 c  + b 
!  FWHM(x) = c(y2-y2) = 3.56c
!  If kappa = <x>/Wmax << 1,  FWHM(x)= 3.98<x> ???
!  We may use,   for  Kappa<<1,
!    0.35*b = FWHM(x) = 3.56c 
!   so b ~ 10 c
!      <x> =11.3c;
!  finally, we obtain
!      c = <x>/11.3  b=<x>10/11.3
!     
!  Method:  
!     let z=exp(-y).  
!   Then,  dz=-zdy.
!     exp(-1/2(y +exp(-y)) dy = exp(-y)^1/2 exp(-1/2exp(-y))dy
!    =   z^1/2 exp(-z/2) /z dz 
!    =   z^(-1/2)exp(-z/2)dz
!  So, sample z from  this distribution.  get y= -log(z).
!
! 
      subroutine ksplandau(b, c, x)
      implicit none
      real*8  b, c, x
      
      real*8  y, z
      

      call ksgmrs(-0.5d0, 2.d0, z)

      y = -log(z)
      x = c*y + b
      end

