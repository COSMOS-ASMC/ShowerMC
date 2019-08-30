!      parameter (n=8)
!      real*8 za(n), aa(n), rn(n)
!      data  za/6., 1., 8., 7., 14., 20., 13., 5./
!      data  aa/12.01, 1., 16., 14., 28.09, 40.08, 26.98, 10.8/
!      data  rn/28., 34.,22., 1., 6., 3., 2., 4./
!      call kradl(za, aa, rn,  n, 1.75d0,     x0cm, x0g)
!      write(*,*) ' x0cm=',x0cm, ' x0g=',x0g
!      sumz=0.
!      suma=0.
!      sum =0.
!      do 100 i=1, n
!         sumz=sumz+rn(i)*za(i)
!         suma=suma+rn(i)*aa(i)
!         sum=sum+rn(i)
! 100  continue
!      write(*,*) ' <z>=',sumz/sum, ' <a>=',suma/sum
!      end
!     ****************************************************************
!     *
!     * kradl: compute radiation length of mixed matter
!     *
!     ****************************************************************
!
!   call kradl(za, aa, rn, n, rho, x0cm, x0g)
!           input
!      za: charge za(i),i=1,n
!      aa: atomic mass   aa(i),i=1,n
!      rn: relative # of atmos  rn(i),i=1,n
!       n: # of different atoms
!     rho: average density  (g/cm**3)
!
!           output
!      x0cm: r.l in cm
!       x0g: r.l in g/cm**2
!
       subroutine kradl(za, aa, rn, n, rho, x0cm, x0g)
       implicit none
       integer n         
       real*8 za(n)
       real*8 aa(n)
       real*8 rn(n)
       real*8 rho

       real*8 x0cm
       real*8 x0g
!
       real*8 sum, x0i ,tmp
       integer i

       sum=0.
       do   i=1, n
           sum=sum + rn(i)*aa(i)
       enddo
       x0i=0.
       do   i=1, n
           call kradl1(za(i), aa(i), tmp)
           x0i=x0i + rn(i)*aa(i)/sum/tmp
       enddo
       x0g=1./x0i
       x0cm=x0g/rho
      end
!     ****************************************************************
!     *                                                              *
!     * kradl1:  compute radiation length of given matter            *
!     *                                                              *
!     *********************** tested 80.07.11 ************************
!
!   /usage/
!            call kradl1(z, a, x0g)
!
!    z:  charge of the matter
!    a:  mass no.
!  x0g:  //                           g/cm**2
!
!
!     *** note ***
!
!         correction to born approximation is not included in this
!         r.l so it must be included in the cross-section.
!
!
!
      subroutine kradl1(z, a, x0g)
      implicit none
      real*8 z
      real*8 a
      real*8 x0g
!
!         cnst=
!         4/137* r0**2 * n  where r0 is the classical electron radius
!                           n the avogadro number
!                           r0=2.8176e-13 cm
!                           n=6.0247
!

      real*8  cnst/1.396e-3/
      real*8 z3, logz3, gzai, t0inv
!

      z3 = z**(-0.3333333)
      logz3 = log( 183.* z3 )
      gzai = log(1440. * z3**2 ) /  logz3
!        inverse of r.l in g/cm**2
      t0inv = cnst / a  *  z*( z + gzai ) * logz3
!
      x0g=1. / t0inv
      end
