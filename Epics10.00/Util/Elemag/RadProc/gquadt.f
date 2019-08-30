!     ****************************************************************
!     *                                                              *
!     * gquadt:  10-point gauss quadrature                           *
!     *                                                              *
!     ****************************************************************
!
!   /usage/
!             call gquadt(f, a, b, ans)
!
!
!
      subroutine gquadt(f, a, b, ans)
!
!
!
      real x(5)/.148874339,.4333953941,.6794095683,.8650633667,
     1   .9739065285/
!
      real w(5)/.2955242247,.2692667193,.2190863625,
     2   .1494513492,.0666713443/
!

      p=(b-a)*.5
      q=p+a
      t=a+b
      ans=0.
      do   i=1,5
         y=p*x(i)+q
         ans=ans+w(i)*(f(y)+f(t-y))
      enddo
      ans=ans*p
      end
do
      ans=ans*p
      end
d
