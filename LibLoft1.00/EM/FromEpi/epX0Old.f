!     ****************************************************************
!     *                                                              *
!     * epX0Old:  compute radiation length of given matter            *
!     *                                                              *
!     *********************** tested 80.07.11 ************************
!
!   /usage/
!            call epX0Old(z, a, x0g)
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
      subroutine epX0Old(z, a, x0g)
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



