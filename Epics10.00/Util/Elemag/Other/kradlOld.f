c      parameter (n=2)
c      dimension za(n), aa(n), rn(n)
c      data  za/1.,  8./
c      data  aa/1., 16./
c      data  rn/2., 1./
c      rho = 1.0
c      call kradl(za, aa, rn,  n, rho,     x0cm, x0g)
c      write(*,*) ' x0cm=',x0cm, ' x0g=',x0g
c      sumz=0.
c      suma=0.
c      sum =0.
c      do 100 i=1, n
c         sumz=sumz+rn(i)*za(i)*aa(i)
c         suma=suma+rn(i)*aa(i)**2
c         sum=sum+rn(i)*aa(i)
c 100  continue
c      write(*,*) ' <z>=',sumz/sum, ' <a>=',suma/sum
c      end
c     ****************************************************************
c     *
c     * kradl: compute radiation length of mixed matter
c     *
c     ****************************************************************
c
c   call kradl(za, aa, rn, n, rho, x0cm, x0g)
c           input
c      za: charge za(i),i=1,n
c      aa: atomic mass   aa(i),i=1,n
c      rn: relative # of atmos  rn(i),i=1,n
c       n: # of different atoms
c     rho: average density  (g/cm**3)
c
c           output
c      x0cm: r.l in cm
c       x0g: r.l in g/cm**2
c
       subroutine kradl(za, aa, rn, n, rho, x0cm, x0g)
       dimension za(n), aa(n), rn(n)
c
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
c     ****************************************************************
c     *                                                              *
c     * kradl1:  compute radiation length of given matter            *
c     *                                                              *
c     *********************** tested 80.07.11 ************************
c
c   /usage/
c            call kradl1(z, a, x0g)
c
c    z:  charge of the matter
c    a:  mass no.
c  x0g:  //                           g/cm**2
c
c
c     *** note ***
c
c         correction to born approximation is not included in this
c         r.l so it must be included in the cross-section.
c
c
c
      subroutine kradl1(z, a, x0g)
c
c         cnst=
c         4/137* r0**2 * n  where r0 is the classical electron radius
c                           n the avogadro number
c                           r0=2.8176e-13 cm
c                           n=6.0247
c
      data  cnst/1.396e-3/
c
      z3=z**(-0.3333333)
      alogz3=alog( 183.* z3 )
      gzai=alog(1440. * z3**2 ) /  alogz3
c        inverse of r.l in g/cm**2
      t0inv=cnst / a  *  z*( z + gzai ) * alogz3
c
      x0g=1. / t0inv
      end



