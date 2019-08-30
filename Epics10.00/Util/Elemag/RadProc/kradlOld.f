!      parameter (n=2)
!      dimension za(n), aa(n), rn(n)
!      data  za/1.,  8./
!      data  aa/1., 16./
!      data  rn/2., 1./
!      rho = 1.0
!      call kradl(za, aa, rn,  n, rho,     x0cm, x0g)
!      write(*,*) ' x0cm=',x0cm, ' x0g=',x0g
!      sumz=0.
!      suma=0.
!      sum =0.
!      do 100 i=1, n
!         sumz=sumz+rn(i)*za(i)*aa(i)
!         suma=suma+rn(i)*aa(i)**2
!         sum=sum+rn(i)*aa(i)
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
       dimension za(n), aa(n), rn(n)
!
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
!
!         cnst=
!         4/137* r0**2 * n  where r0 is the classical electron radius
!                           n the avogadro number
!                           r0=2.8176e-13 cm
!                           n=6.0247
!
      data  cnst/1.396e-3/
!
      z3=z**(-0.3333333)
      alogz3=alog( 183.* z3 )
      gzai=alog(1440. * z3**2 ) /  alogz3
!        inverse of r.l in g/cm**2
      t0inv=cnst / a  *  z*( z + gzai ) * alogz3
!
      x0g=1. / t0inv
      end
*2
      t0inv=cnst / a  *  z*( z + gzai ) * alogz3
!
      x0g=1. / t0inv
      end
