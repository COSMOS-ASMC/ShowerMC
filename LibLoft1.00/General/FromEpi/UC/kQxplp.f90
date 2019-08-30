!     ****************************************************************
!     *                                                              *
!     * kQxplp: crossing point of a given line with a given plane     *
!
! /usage/   call kQxplp(x0, y0, z0, l, m, n, a, b, c, d,
!         *      el, icon)
!
!
! -- input --
!    x0,y0,z0:  a point the line passes.  real*8
!    l, m, n : x,y,z componets of the dirction cosine of the line
!              (real*8) (unit vector)
!             ( line is expressed by  (x-x0)/l=(y-y0)/m=(z-z0)/n  )
!
!     a,b,c,d:  coefficients to express the plane.  ax+by+cz=d   is the
!               plane.  real*8.   a,b,c are the direction cos of a line
!               perpendicular to the plane, d is the distance to the
!               plane from the origin, if (a,b,c) is a unit vector.
!               (a,b,c) may not be unit vector.
!
!  -- output --
!      el:  real*8    crossing point is at  (x0,y0,z0)+el*(l,m,n)
!              el>=0 if xpoint is on the l,m,n direction
!                    else negative.
!       icon:  0 when a crossing point is obtained
!              1 when the line is on the plane
!              2 when the line is paralell to the plane but not on the
!                plane
!
!                when icon^=0, el is unchanged.
!
!    ** note ** no check is made on the consistency of l,m,n
!               and a,b,c
!
      subroutine kQxplp(x0, y0, z0, l, m, n, a, b, c, d, &
           el, icon)
        implicit none
!
!
      real(16):: x0, y0, z0,  l, m, n, el, a, b, c, d
      integer icon
!
      real(16):: div, g, an, bn, cn, dn, dist2
      real(16),parameter:: eps=1.q-10

      if((abs(a)+abs(b)+abs(c)) .eq. 0.) then
          icon=3
      else
          div= a*l + b*m + c*n
!          if(abs(div) .gt. 0.d0) then
!!!!!!!!!!
!              write(0,*) ' div=',div
!!!!!!!!!!!!!
          if(abs(div) .gt. eps) then
!              crossing point exists
              el=(d- (a*x0+b*y0+c*z0) )  /  div
              icon=0
          else
!               // or one the plane.  comput distance from (x0,...)
!               to the plane
!                    normalize coeficient to avoid overflow
              g=max( abs(a), abs(b), abs(c) )
              an=a/g
              bn=b/g
              cn=c/g
              dn=d/g
              dist2= ( an*x0 + bn*y0 + cn*z0  -dn )**2  / &
                   (an**2 + bn**2 + cn**2)
!!!!!!!!!!
!              write(0,*) ' dist2=',dist2
!!!!!!!!!!!!!
              if(dist2 .le. eps) then
!                  on the plane
                 icon=1
              else
!                 no cross
                 icon=2
              endif
          endif
      endif
    end subroutine kQxplp
