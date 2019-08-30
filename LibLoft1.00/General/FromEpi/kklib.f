#if defined NEXT486
#define IMAG_P dimag
#elif defined PCLinux
#define IMAG_P dimag
#else
#define IMAG_P imag
#endif
      

!      kxplbx:  box and line x-point ; see aslo kxplbxnew in prog/
!      kxplcy:  cyl and line x-point
!      kxplp:   plane and line x-point
!      kxplsl:  line and segment of line x-point
!      kxplineseg: x-point of two line segments 
!      kxplPrism: line and prism x-point
!      k3dclp:  3d clipping
!      kfige:   find int data pos. >= a given value
!      ksetiv:  set an integer value in a given array
!      ksetrv:  set a real*8 value in a given array
!      k3inout: see if a point is inside of a triangle.
!      kxplsph: sphere and line x-point
!      kioPrism: judge if a given point is inside of a prism
!      kgetField:
!      kfindField:  similar to abvoe
!      kgetCpos:
!      kgetBpos:
!      kshiftB: shift blank until non blank appears
!      kshiftC: shift character until blank  appears
!      kDebye:  Debye function
!      kBern:   Bernuoulli number
!      kalpha:  see if string starts with a-z or A-Z
!      kxplellip: crossing point of a line with an ellipse
!      ****************************
      logical function kalpha(string)
!
!         see if first character of string is one of a-z or A-Z
!
      implicit none
      character*(*) string
      logical first/.true./
      integer a, z, CA, CZ, i
      save first, a, z, CA, CZ
      
      if(first) then
         a = ichar('a')
         z = ichar('z')
         CA = ichar('A')
         CZ = ichar('Z')
         first = .false.
      endif
     
      i = ichar(string(1:1))
      if(i .ge. a .and. i .le. z) then
         kalpha = .true.
      elseif(i .ge. CA .and. i .le. CZ) then
         kalpha = .true.
      else
         kalpha = .false.
      endif
      end
!          test kcplcy
!      real*8 l, m, n, x0, y0, z0, el
!      x0=-0.9
!      y0=-0.9
!      z0=-0.5d0
!      l=(1.41421356/2)**2
!      m=l
!      n= sqrt(1.d0- l**2 - m**2)
!      write(*,*) ' n=',n
!      r=1.
!      h=2.
!      call  kxplcy(x0, y0, z0, l, m, n,
!    *         r, h,  el, icon, icon2)
!              write(*,*) ' el=',el, ' icon=',icon, ' icon2=',icon2
!       end
!       ****************************************************************
!       *
!       *  kxplcy: crossing point of a line with a cylinder
!       *
!       ************************ tested 89.09.06 ***********************
!
!  call kxplcy(x0, y0, z0, l, m, n, r, h,  el,
!                icon, icon2)
!  input:
!     x0,y0,z0:  a point the line passes, real*8
!     l, m, n: (real)  direction cos of the line, real*8
!     r: radius of cylinder. real*8
!     h: height of //     bottom is assumed to be on z=0
!                         axis is assumed to be z- axis
!                         top is assumed to be z=h  real*8
!  output:
!     el: x-ssing point is at (x0,y0,z0)+el*(l,m,n)  el>0
!         only x-ssing point with el>0 is obtained.
!    From Dec.28.2004, el==0 case is omitted as crsossing case.
!
!         if there is two el>0 xssing points, nearer one is
!         taken.
!
!   icon : output. 0 the point is in the cyl. el is obtained
!                    if icon2 != -1. 
!                  1 the point is out side of the cyl. el is
!                     obtained.
!                 -1 no x-ing point
!   icon2: output. 1  x-ing point is on x-y  top plane.
!                  2  //             on the side.
!                  6  //             on      bottom.
!                 -1  no x-ing point
!
!   ****** note ******
!          this is designed to be fast if the x-ssing point is
!          at the top or bottom of the cylinder.
!
       subroutine  kxplcy(x0, y0, z0, l, m, n,
     *         r, h,  el, icon, icon2)
       implicit none
!
       real*8 x0,y0,z0, l, m, n, el, r, h
       integer*4 icon, icon2
!
       real*8  d, wxy, x0y0, rsq, x, y, b1, ds, z
!
           rsq=r**2
           x0y0=x0**2+ y0**2 - rsq
           if(z0 .ge. 0. .and.  z0 .le. h  .and.  x0y0 .le. 0.d0) then
               icon =0
           else
               icon =  -1   ! at least outside. may be changed later
           endif
           icon2 = -1
           if(n .gt. 0.d0 .and. icon .eq. 0) then
               el=(h-z0)/n
               x=x0+ el*l
               y=y0+ el*m
               if(x**2 + y**2 .le. rsq  .and. el .ge. 0.) then
                  icon2 = 1
               endif
           elseif(n .lt. 0.d0 .and. icon .eq. 0) then
               el=-z0/n
               x=x0+ el*l
               y=y0+ el*m
               if(x**2 + y**2 .le. rsq .and. el .ge. 0. ) then
                  icon2 = 6
               endif
           endif
           if(icon2 .eq. -1) then
              wxy=l**2 + m**2
              if(wxy .eq. 0.) then
!                  vertical. from outside
                 if(z0 .lt. 0. .and. n .gt. 0. .and.
     *               x0y0 .le. 0.) then
                    icon = 1
                    icon2= 6
                    el = -z0
                 elseif(z0 .gt. h .and. n .le. 0. .and.
     *               x0y0 .le. 0.) then
                    icon2 = 1
                    icon = 1
                    el = (h-z0)/n
                 endif
              endif
           endif

           if(icon2 .eq. -1 .and. wxy .ne. 0.) then
              b1=l*x0+m*y0
              d= b1**2 - wxy* x0y0

              if(d .ge. 0.d0)then
                 ds=sqrt(d)
                 if(x0y0 .le. 0.d0) then
!                       point is inner part of cylinder
                     el=  (-b1 + ds)/wxy
                     z=z0+el*n
                     if(z0 .ge. 0. .and. z0 .le. h 
     *                   .and. el .ge. 0.) then
                        icon = 0
                        icon2 = 2
                     elseif(z0 .gt. h .and. n .lt. 0.d0) then
                        el=(h-z0)/n
                        x=x0+el*l
                        y=y0+el*m
                        if(x**2+ y**2 .le. rsq .and.  el .ge. 0. ) then
                           icon = 1
                           icon2 = 1
                        endif
                     elseif(z0 .lt. 0.d0 .and. n .gt. 0.d0) then
                        el=-z0/n
                        x=x0+el*l
                        y=y0+el*m
                        if(x**2+ y**2 .le. rsq .and. el .ge. 0.) then
                           icon =  1
                           icon2 = 6
                        endif
                     endif
                 elseif(b1 .le. 0.) then
!                       point is outside.  x-ssing forward with cyl
                      el=  (-b1- ds)/wxy
                      z=z0+el*n
                      if(z .ge. 0. .and. z .le. h 
     *                   .and. el .ge. 0. ) then
                         icon = 1
                         icon2 = 2
                      elseif(n .lt. 0.d0) then
                         el=(h-z0)/n
                         x=x0+el*l
                         y=y0+el*m
                         if(x**2+ y**2 .le. rsq .and. el .ge. 0. ) then
                            icon = 1
                            icon2 = 1
                         endif
                      elseif(n .gt. 0.d0) then
                         el=-z0/n
                         x=x0+el*l
                         y=y0+el*m
                         if(x**2+ y**2 .le. rsq .and. el .ge. 0. ) then
                            icon = 1
                            icon2 = 6
                         endif
                      endif
                 endif
             endif
           endif
       end
!     ****************************************************************
!     *                                                              *
!     * kxplp: crossing point of a given line with a given plane     *
!     *                                                              *
!     *********************** tested 89.09.06 ************************
!
! /usage/        call kxplp(x0, y0, z0, l, m, n, a, b, c, d,
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
      subroutine kxplp(x0, y0, z0, l, m, n, a, b, c, d,
     *           el, icon)
       implicit none
!
!
      real*8 x0, y0, z0,  l, m, n, el, a, b, c, d
      integer icon
!
      real*8 div, g, an, bn, cn, dn, dist2

      if((abs(a)+abs(b)+abs(c)) .eq. 0.) then
          icon=3
      else
          div= a*l + b*m + c*n
          if(abs(div) .gt. 0.d0) then
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
              dist2= ( an*x0 + bn*y0 + cn*z0  -dn )**2  /
     *                       (an**2 + bn**2 + cn**2)
              if(dist2 .eq. 0.) then
!                  on the plane
                 icon=1
              else
!                 no cross
                 icon=2
              endif
          endif
      endif
      end
!     complex zz0/(0., 0.)/
!     complex expa, zz1, zz2
!     expa=exp( cmplx(0., 1.,8)*(90.01 /180.*3.141592 ) )
!     zz1=cmplx(0., 1.,8 )
!     zz2=cmplx(1., 0.,8 )
!     eps=1.e-5
!     call kxplsl(zz0, expa, zz1, zz2, eps, p, q, icon)
!     write(*,*) ' icon=',icon, ' p=',p, ' q=',q
!     end
!     ****************************************************************
!     *                                                              *
!     * kxplsl:  get crossing point of a segment and a line          *
!     *                                                              *
!     ****************************************************************
!
! /usage/
!             call kxplsl(zz0,expa, zz1,zz2, eps, p,q, icon)
!
!               zz1 *
!                   p .
!                    /
!                   /       *
!                q /         zz2
!                 /  <---complex exp of this angle = expa
!               * ------------
!                zz0
!
!   zz0--given point in complex
!  expa--exp(i*angle) of a line which passes zz0
!   zz1, zz2-- given points in complex to denote a segment
!     p--distance from zz1 to the crossing point (singned)
!        this is the portion of the segment length.
!     q--distance from zz0 to    //
!   eps-->=0.  used to judge the parallity of two lines as well as
!        to judge the separation of two // lines.
!
!  __return condition__
!
!      icon=0:        crossing point found on the segment (0 <= p <= 1)
!           1         the line overlaps the segment.  the position of
!                     zz0 may be judged by q.  if q < 0, zz0 is 'left'
!                     to zz1,  if q > 1, zz0 is 'right' to zz2
!                     if q=0, zz0 is on the segment
!           2         crossing point is outside of the segment
!                     if p < 0, it is 'left' to zz1, if p > 1, it is
!                     'ritht' to zz2
!           3         no crossing point at all, i.e., they are parallel
!                     in this case q becomes distance between two lines
!
! *** note ***    they are // if the relative angle of both lines
!                 are < eps
!                 if the relative separtion of two // lines are
!                 less than eps, they are judeged to be the same
!
!
!      no sub needed
!
!
      subroutine   kxplsl(zz0,expa, zz1,zz2, eps, p,q, icon)
       implicit none
!
!
      complex*16 zz0, expa, zz1, zz2
      real*8 eps, p, q

!
      complex*16 z0, z1, z2, cexpa, z10, z20, z21
      integer icon, i200
      real*8 t, t1, t2
!
!       rotate zz0-line to be // to x-axis
      cexpa=conjg(expa)
      z0=zz0*cexpa
      z1=zz1*cexpa
      z2=zz2*cexpa
!
      z21=z2-z1
      z10=z1-z0
      z20=z2-z0
!        see if not right angle crossing
      if( real(z21) .ne. 0.) then
!           tan of relative angle
         t=IMAG_P(z21)/real(z21)
      else
         t=1.
      endif
!          see if they are //
      if( abs(t) .gt. eps )then
!              non //, so find + point
!          see if z1 = z2
          if(IMAG_P(z21) .eq. 0.) then
!              see if z0-line crosses the z1 and z2
             if(IMAG_P(z10) .eq. 0.) then
!                z1 is on the line
                q=real(z10)
                p=0.
                icon=0
             else
                 write(*, 320) zz1,zz2
  320            format('0*** error input to kxplsl;  z1=z2',4g13.3)
                 icon=2
             endif
         else
             p=- IMAG_P(z10)/IMAG_P(z21)
             q= real(z10) +  real(z21) * p
             if(p .lt.-eps .or.  p .gt. 1.+eps) then
                 icon=2
             else
                 icon=0
             endif
         endif
      else
!           if z0=z1  goto 200
         if(real(z10) .ne. 0.) then
             t1=z10*conjg(z10)
             t2=z20*conjg(z20)
             if( t1 .le. t2  ) then
!                z0 is closer to z1, so examine // of z0-line and z0-z2
                 t= IMAG_P(z20)/real(z20)
             else
                 t=IMAG_P(z10)/real(z10)
             endif
             if( abs(t) .gt. eps ) then
                 q=IMAG_P(z10)
                 icon=3
                 i200=0
             else
                 i200=1
             endif
         else
             if(abs(IMAG_P(z10)) .le. eps) then
                i200=1
             else
                i200=0
                q=IMAG_P(z10)
                icon=3
             endif
         endif
         if(i200 .eq. 1) then
              p=0.
              q= real(z10)
!                see if z0 is 'left' to z1
              if( q .gt. 0.) then
                  icon=1
              else
                  p=1.
                  q=real(z20)
!                  see if z0 is 'right' to z2
                  if( q .lt. 0. ) then
                      icon=1
                  else
!                         z0 is on the segment
                      q=0.
                      icon=1
                  endif
              endif
          endif
      endif
      if(icon .eq. 0) then
         if(p .gt. 1.) p=1.
         if(p .lt. 0.) p=0.
      endif
      end
      subroutine kFormBoxFace(fn, a, b, c, p1, p2, p3, p4)
      implicit none
      ! form four vertex points of a box face with # fn
      integer,intent(in):: fn ! box face #. 
             !  z=0--> 1
             !  z=c--> 6
             !  x=0--> 5
             !  x=a--> 2
             !  y=0--> 3
             !  y=b--> 4
      real(8),intent(in):: a, b, c  ! 3 edge lengths
      real(8),intent(out):: p1(3), p2(3), p3(3), p4(3)

      if(fn == 1) then  ! at z=0
         p1(:) = 0.d0
         p2(3) = 0.d0
         p3(3) = 0.d0
         p4(3) = 0.d0

         p2(1) = 0.d0
         p2(2) = b
         
         p3(1) = a
         p3(2) = b
         
         p4(1) = a
         p4(2) = 0.d0
      elseif( fn == 2) then   ! x=a
         p1(1) = a
         p2(1) = a
         p3(1) = a
         p4(1) = a
         
         p1(2:3) = 0.d0
         p2(2) = b
         p2(3) = 0.d0

         p3(2) = b
         p3(3) = c

         p4(2) = 0.d0
         p4(3) = c
      elseif( fn == 3 ) then  ! y = 0
         p1(:) = 0.d0
         p2(2) = 0.d0
         p3(2) = 0.d0
         p4(2) = 0.d0
         
         p2(1) = a
         p2(3) = 0.d0
         p3(1) = a
         p3(3) = c

         p4(1) = 0.d0
         p4(3) = c
      elseif( fn == 4 ) then ! y = b
         p1(2) = b
         p2(2) = b
         p3(2) = b
         p4(2) = b
         
         p1(1) = 0.d0
         p1(3) = 0.d0

         p2(1) = 0.d0
         p2(3) = c
         p3(1) = a
         p3(3) = c
         p4(1) = a
         p4(3) = 0.d0
      elseif( fn == 5 ) then  ! x =0
         p1(:) = 0.d0
         p2(1) = 0.d0
         p3(1) = 0.d0
         p4(1) = 0.d0

         p2(2) = 0.d0
         p2(3) = c
         p3(2) = b
         p3(3) = c
         p4(2) = b
         p4(3) = 0.d0
      elseif( fn == 6) then  ! z=c
         p1(3) = c
         p2(3) = c
         p3(3) = c
         p4(3) = c

         p1(1:2) = 0.d0
         p2(1) = a
         p2(2) = 0.d0
         p3(1) = a
         p3(2) = b
         p4(1) = 0.d0
         p4(2) = b
      else
         write(0,*) ' fn=', fn, ' invalid for kFormBoxFace'
         stop
      endif
      end subroutine kFormBoxFace

      function kInsideBox(x, y, z, a, b, c) result(ans)
      implicit none
!      judege if a point is inside of a give canonical box
      real(8),intent(in):: x, y, z  ! a point
      real(8),intent(in):: a, b, c  ! 3 edges of the box
      
      logical:: ans !  t if the point is inside (including on
                    ! the surface
                    !  f if outside  

      if(  (x-a)*x >  0.d0 ) then
         ans = .false.
      elseif( (y-b)*y >  0.d0 ) then
         ans = .false.
      elseif( (z-c)*z > 0.d0 )  then
         ans = .false.
      else
         ans = .true.
      endif
      end   function kInsideBox

      
!        test kxplbx: x-ssing point of a half line with a box
!
!        implicit  none
!        real*8 x0,y0,z0, l, m, n, el,a,b,c
!        integer icon
!
!        a=1.
!        b=1.
!        c=1.
!        x0= 0.5d0
!        y0= 2.5d0
!        z0=  1.d0
!        l  = 0.
!        m= 1.
!        n= 0.d0
!        write(*,*) ' enter x0,y0,z0, n,m '
!        read(*, *) x0,y0,z0, n,m
!        l = sqrt(1.d0-n**2- m**2)
!        write(*,*)  ' l=',l
!        call kxplbx(x0, y0, z0, l, m, n, a, b,c,el, icon)
!        write(*,*) ' icon=',icon, ' el=',el
!        end
!     ****************************************************************
!     *                                                              *
!     * kxplbx: crossing point of a line with a given box            *
!     *                                                              *
!     *********************** tested 89.09.06 ************************
!
! /usage/       call kxplbx(x0, y0, z0,  l, m, n, a, b, c,
!       *            el, icon)
!
!             z
!             !
!             !                  one corner of the box is at (0,0,0).
!             !                  the lengths of 3 sieds of the box are
!             !                  a, b, c and are on x,y,z axes, resp.
!             /-------------  y
!           /
!         /
!       /
!      x
!
! -- input --
!  x0,y0,z0:  the given line passes this point. real*8
!   l, m, n:  (x,y,z) components of the direction cosines of the line
!             real*8
!   a,b,c:    length of each side of the box
!
!
!  -- output --
!      el: crossing point is at (x0,y0,z0)+el*(l,m,n) real*8
!          only el>=0 is obtained.  if two such points exist,
!          nearest one is taken.
!   icon:  0 el is obtained. x0,y0,z0 is inside of the box
!          1 el is obtained. //       is outside of the box.
!         -1 no x-point el is undef.
! **** if kxplbx defined above has problems, next
!      may be used 
      subroutine kxplbx(x0, y0, z0,  l, m, n,  a, b,  c,
     *            el,  icon)
       implicit none
!
      real*8 x0, y0, z0, l, m, n, el, a, b, c
      integer icon
!
      real*8 x1, y1, z1, x2, y2, z2, el1, el2,
     *       xo1, yo1, zo1, xo2, yo2, zo2
      real(8),parameter::eps = 1.0d-8  ! (cm)

      logical kInsideBox  ! external func.
!      if(abs(n) > 1.0d-1) then
       if(abs(n) > 0.577d0) then  !better. at least one of abs( l,m,n ) > 1/sqrt(3)
            el1=(c-z0)/n
            el2= - z0/n
            x1= x0+ el1*l
            y1= y0+ el1*m
            z1= c
            x2= x0+ el2*l
            y2= y0+ el2*m
            z2= 0.
            call k3dclp(x1, y1, z1, x2, y2, z2, a, b, c,
     *      xo1, yo1, zo1, xo2, yo2, zo2, icon)
            if(icon .eq. 0) then
                el1=(zo1-z0)/n
                el2=(zo2-z0)/n
            endif
!       elseif(abs(l) .ge. 1.d-1) then
       elseif(abs(l) > 0.577d0) then
            el1=(a-x0)/l
            el2= -x0/l
            x1=a
            y1=y0+ el1*m
            z1=z0+ el1*n
            x2=0.
            y2=y0+ el2*m
            z2=z0+ el2*n
            call k3dclp(x1, y1, z1, x2, y2, z2, a, b, c,
     *      xo1, yo1, zo1, xo2, yo2, zo2, icon)
            if(icon .eq. 0) then
                el1=(xo1-x0)/l
                el2=(xo2-x0)/l
            endif
        else
            el1=(b-y0)/m
            el2= -y0/m
            x1=x0+ el1*l
            y1=b
            z1=z0+ el1*n
            x2=x0+ el2*l
            y2=0.
            z2=z0+ el2*n
            call k3dclp(x1, y1, z1, x2, y2, z2, a, b, c,
     *      xo1, yo1, zo1, xo2, yo2, zo2, icon)

            if(icon .eq. 0) then
                el1=(yo1-y0)/m
                el2=(yo2-y0)/m
            endif
        endif
        if(icon .eq. 0) then
            if(el1 .gt. 0.) then
                if(el2 .gt. 0.) then
                   el=min(el1, el2)
                   if( kInsideBox(x0, y0, z0, a, b, c) ) then
                   ! this happens when the point is outside
                   ! but due to  numerical error in some case
                   ! so we rely on kInsideBox
                      if( el > eps ) then
                         ! too big ; warning
                         write(0,*) ' point ',x0,y0,z0
                         write(0,*) ' is inside of a box'
                         write(0,*) ' but two crossing points'
                         write(0,*) ' with dir =',l,m,n
                         write(0,*) ' and lengths ', el1, el2
                         stop
                      else
                     ! regard numerical error ~ (< 10^-8)
                         icon = 0
                         el = max(el1, el2)
                      endif
                   else
                      icon = 1
                   endif
                elseif(el2 .le. 0.) then
                    el=el1
                    icon=0
                endif
            elseif(el1 .lt. 0.) then
                if(el2 .ge. 0.) then
                   icon=0
                   el=el2
                else
                   icon=-1
                endif
            else
                if(el2 .gt. 0.) then
                    el=el2
                    icon=0
                else
                    el=el1
                    icon=0
                endif
            endif
        endif
       end
!      
      subroutine kxplsph(x0, y0, z0, l, m, n, r, el, icon)
      implicit none
      real*8  x0, y0, z0 ! input. the line passes this point
      real*8  l, m, n  !  input.  direc cos.  of  the line
      real*8  r        !  input.  radius of the sphere
      real*8  el       !  output. el>=0 distance to the
                       !          sphere  from  x0,y0,z0
      integer icon    !  output. icon =0.  x-point exists 
                      !                  x0,.. is inside
                      !          icon = 1  x-point exists
                      !                  x0.. is outside
                      !                =-1.  no x-point

      real*8  rsqr, r0l, d
      integer icon1, icon2 
      
      rsqr = x0**2 + y0**2 + z0**2 -r**2
      if(rsqr .le. 0.) then
!          inside
         icon2 = 0
      else
         icon2 = 1
      endif
      r0l = x0*l + y0*m + z0*n
      d = r0l**2 - rsqr
      if(d .ge. 0.) then
         d = sqrt(d)
         el = -r0l - d
         if(el .ge. 0.) then
            icon1 = 0
         else
            el = -r0l + d
            if(el .ge. 0.) then
               icon1 = 0
            else
               icon1 = 1
            endif
         endif
      else
         icon1 = 1
      endif
!
      if(icon2 .eq. 0) then
         icon = 0
      elseif(icon1 .eq. 0) then
         icon = 1
      else
         icon = -1
      endif
      end

!        *********************************************************
!        *
!        * k3dclp: 3d clipping of a line segment by a box
!        *
!        *********************************************************
!
! /usage/  call k3dclp(xi0, yi0, zi0, xi1, yi1, zi1, a, b, c,
!         *            xo0, yo0, zo0, xo1, yo1, zo1, icon)
!
!       xi0, yi0, zi0: input.  1st point of the line segment
!       xi1, yi1, zi1: input.  2nd point of the line segment
!         a, b, c: input.  the lenght of three edges of the box.
!                   the three edges lie on (0,0,0)-(a,0,0)
!                                          (0,0,0)-(0,b,0)
!                                          (0,0,0)-(0,0,c)
!       xo0,yo0,zo0: output.  1st point of the segment
!       xo1,yo1,zo1: output.  2nd point of the segemnt
!                             these may be orignal points
!                             if points are contained in the
!                             box, or the point(s) on the surface
!                             of the box where the segment accrosses
!                             the box or not given.
!       icon: ouput.          =-1 --> no crossing point at all.
!                             = 0 --> crossing points exist or contained
!
       subroutine k3dclp(xi0, yi0, zi0, xi1, yi1, zi1, a, b, c,
     *                   xo0, yo0, zo0, xo1, yo1, zo1, icon)
       implicit none
!
!               ix : to set bit at x-th position  x=1,2, from right
!
       real*8 xi0, yi0, zi0, xi1, yi1, zi1, xo0, yo0, zo0,
     *        xo1, yo1, zo1, a, b, c
       integer icon
       
       integer i6, i5, i4, i3, i2, i1
       parameter (i6=5, i5=4, i4=3, i3=2, i2=1, i1=0)
       integer bp6, bp5, bp4, bp3, bp2, bp1
!              bit pattern whose x-th bit is on, where x is bpx)
       parameter (bp6=2**i6, bp5=2**i5, bp4=2**i4,
     *            bp3=2**i3, bp2=2**i2, bp1=2**i1)
       logical ok
       integer jc0, jc1, itmp
       real*8 x0, y0, z0, x1, y1, z1
       real*8 tmpx, tmpy, tmpz, t
!
       jc0=0
       jc1=0
!
       x0=xi0
       y0=yi0
       z0=zi0
       x1=xi1
       y1=yi1
       z1=zi1
!
       if(z0 .gt. c ) then
           jc0=ibset(jc0, i6)
       endif
       if(z0 .lt. 0.) then
           jc0=ibset(jc0, i5)
       endif
       if(y0 .gt. b) then
           jc0=ibset(jc0, i4)
       endif
       if(y0 .lt. 0.) then
           jc0=ibset(jc0, i3)
       endif
       if(x0 .gt. a) then
           jc0=ibset(jc0, i2)
       endif
       if(x0 .lt. 0.) then
           jc0=ibset(jc0, i1)
       endif
!
       if(z1 .gt. c ) then
           jc1=ibset(jc1, i6)
       endif
       if(z1 .lt. 0.) then
           jc1=ibset(jc1, i5)
       endif
       if(y1 .gt. b) then
           jc1=ibset(jc1, i4)
       endif
       if(y1 .lt. 0.) then
           jc1=ibset(jc1, i3)
       endif
       if(x1 .gt. a) then
           jc1=ibset(jc1, i2)
       endif
       if(x1 .lt. 0.) then
           jc1=ibset(jc1, i1)
       endif
!
!       *** until loop*** 
       do while (.true.)
           ok=jc0 .eq. 0 .and. jc1 .eq. 0
           if(ok) then
               icon=0
               xo0=x0
               yo0=y0
               zo0=z0
               xo1=x1
               yo1=y1
               zo1=z1
           else
               ok=iand(jc0, jc1) .ne. 0
               if(ok) then
                   icon=-1
               else
                   if(jc0 .eq. 0) then
                      tmpx=x0
                      tmpy=y0
                      tmpz=z0
                      x0=x1
                      y0=y1
                      z0=z1
                      x1=tmpx
                      y1=tmpy
                      z1=tmpz
                      itmp=jc0
                      jc0=jc1
                      jc1=itmp
                   endif
                   if(iand(jc0, bp6) .ne. 0) then
                       t=(c-z0)/(z1-z0)
                       z0=c
                       x0=x0 + (x1-x0) * t
                       y0=y0 + (y1-y0) * t
                   elseif(iand(jc0,bp5) .ne. 0) then
                       t=  -z0/(z1-z0)
                       z0=0.
                       x0=x0 + (x1-x0) * t
                       y0=y0 + (y1-y0) * t
                   elseif(iand(jc0, bp4) .ne. 0) then
                       t=(b-y0)/(y1-y0)
                       y0=b
                       x0=x0 + (x1-x0) * t
                       z0=z0 + (z1-z0) * t
                   elseif(iand(jc0, bp3) .ne. 0) then
                       t=  -y0/(y1-y0)
                       y0=0.
                       x0=x0 + (x1-x0) * t
                       z0=z0 + (z1-z0) * t
                   elseif(iand(jc0, bp2) .ne. 0) then
                       t=(a-x0)/(x1-x0)
                       x0=a
                       y0=y0 + (y1-y0) * t
                       z0=z0 + (z1-z0) * t
                   else
                       t=  -x0/(x1-x0)
                       x0=0.
                       y0=y0 + (y1-y0) * t
                       z0=z0 + (z1-z0) * t
                   endif
!
                   jc0=0
                   if(z0 .gt. c ) then
                       jc0=ibset(jc0, i6)
                   endif
                   if(z0 .lt. 0.) then
                       jc0=ibset(jc0, i5)
                   endif
                   if(y0 .gt. b) then
                       jc0=ibset(jc0, i4)
                   endif
                   if(y0 .lt. 0.) then
                       jc0=ibset(jc0, i3)
                   endif
                   if(x0 .gt. a) then
                       jc0=ibset(jc0, i2)
                   endif
                   if(x0 .lt. 0.) then
                       jc0=ibset(jc0, i1)
                   endif
               endif
           endif
       if         (ok)
     *                    goto 100
       enddo
  100  continue
      end
!      *********************************************************
!      *
!      * check if all the array elements are the same
!      *
!      *********************************************************
!
!  /usage/  call kcsame(x, intv, n, icon)
!
       subroutine kcsame(x, intv, n, icon)
       implicit none
       integer intv, n, icon
       real*8 x(intv, n)

       integer i
       real*8 s
          icon=0
          i=2
          s=x(1,1)
           do   i=2, n
              if(x(1,i) .ne. s) then
                 icon=1
                 goto 200
               endif
           enddo
  200     continue
        end
!      integer ia(2,10)
!     
!      do i = 1, 10
!         ia(1, i)= i
!      enddo
!      ia(1, 5) =9
!      call kfige(ia, 2, -10, 5, m, icon)
!      write(*, *) ' m=', m, ' icon=',icon, ' data=', ia(1, m)
!      end
!     ****************************************************************
!     *                                                              *
!     *  kfige: find integer data (position) .ge. given value        *
!     *                                                              *
!     ***********************  tested 81.04.10  **********************
!
!   /usage/
!          call kfige(x, intvx, n, c, m, icon)
!
!     x:  integer data array
! intvx:  interval of data in x
!     n:  !n! is no. of data in x
!     c:  given value.  x  .gt.  c is sought for.
!
!     m:  position of found data
!  icon:  0 if found else 1 results
!
! *** note ***
!         if n>0 search is made for from 1st, else from last
!
!
      subroutine kfige(x, intvx, n, c,  m, icon)
      implicit none
!
!
          integer x, intvx, n, c,  m, icon
!          dimension x(intvx,  n)
          dimension x(intvx,  *)
!
          integer istep, i1, i2, i

          if( n .eq. 0 ) then
             icon = 1
             return               !  **************
          elseif(n .gt. 0) then
             istep = 1
             i1 = 1
             i2 = n
          else
             istep = -1
             i1 = -n
             i2 = 1
          endif

          do i = i1, i2, istep
             if( x(1,i) .ge. c) then
                m = i
                icon = 0
                return       ! ********************
             endif
          enddo
          if(n .gt. 0) then
             m = n + 1
          else
             m = 0
          endif
          icon = 1
      end
!     ****************************************************************
!     *                                                              *
!     *  ksetrv: set given real*8 value in a given real array)       *
!     *                                                              *
!     *********************** tested 84.06.28 ************************
!
!   /usage/  call ksetrv(a, intv, n, v )
!
!   array positions of a(1), a(1+intv), ... a(1+(n-1)*intv) are
!   put the value of v.
!
!
      subroutine ksetrv(a, intv, n, v)
      implicit none
      
!
      real*8 a, v
      integer intv, n
      dimension a(intv, n)
!
      integer i 

           do   i=1, n
               a(1,i)=v
           enddo
      end
!     ****************************************************************
!     *                                                              *
!     *  ksetiv: set given interger value in a given interger array) *
!     *                                                              *
!     *********************** tested 84.06.28 ************************
!
!   /usage/  call ksetiv(a, intv, n, v )
!
!   array positions of a(1), a(1+intv), ... a(1+(n-1)*intv) are
!   put the value of v.
!
!
      subroutine ksetiv(a, intv, n, v)
      implicit none
!
          integer a, intv, n,  v
          dimension a(intv, n)
!
          integer i

           do   i=1, n
               a(1,i)=v
           enddo
      end
!	real*8 x, y, x1,y1, x2, y2, x3,  y3
!	integer inout
!	do while(.true.)
!	   write(*, *) ' x1,... y3'
!	   x1=0.
!	   y1 = 0.d0
!	   x2 = 1.d0
!	   y2 =  0.d0
!	   x3 =  .5d0
!	   y3 =  2.d0
!	   read(*, *) x1,y1, x2, y2, x3, y3
!	   x = 0.
!	   do  while (x .ne. -100.d0)
!	      write(*, *) 'x,y=' 
!	      read(*, *)  x, y
!	      call k3inout(x, y, x1, y1, x2, y2, x3, y3, inout)
!	      write(*, *) ' inout=', inout
!	   enddo
!	enddo
!	end

	subroutine k3inout(x, y, x1, y1, x2, y2, x3, y3, inout)
	implicit none
	real*8 x,y  ! input. some point to be judged if it is inside or
                    ! outsite of the triangle
        real*8 x1, y1 ! input.  1 point of the triangle
        real*8 x2, y2 ! input.  2nd point
        real*8 x3, y3 ! input.  3rd point
	integer inout ! output. 0--> in or on the line, 1 out
	logical anticw
!
	real*8  f, xa, ya, xb, yb, xc, yc
	f(xa, ya, xb, yb, xc, yc) = (xb- xa)*(yc-yb) - (yb-ya)*(xc-xb)

	inout = 1
	anticw = f(x1, y1, x2, y2, x3, y3) .gt. 0.d0

        if(anticw) then
           if( f(x1, y1, x2, y2, x, y) .lt. 0.d0) goto 10
           if( f(x2, y2, x3, y3, x, y) .lt. 0.d0) goto 10
           if( f(x3, y3, x1, y1, x, y) .lt. 0.d0) goto 10
           inout = 0
        else
           if( f(x1, y1, x2, y2, x, y) .gt. 0.d0) goto 10
           if( f(x2, y2, x3, y3, x, y) .gt. 0.d0) goto 10
           if( f(x3, y3, x1, y1, x, y) .gt. 0.d0) goto 10
           inout = 0
	endif
 10     continue
        end

! 
!       implicit none
!       character*15 field(100)
!       character*120 cdata
!       integer i, nf
!
!       do while (.true.)
!          read(*, '(a)',end=100) cdata
!          write(*,'(a)') cdata
!          call kgetField(cdata, field, nf)
!          do i = 1, nf
!             write(*,*) i, field(i)
!          enddo
!       enddo
! 100   continue
!       end

       subroutine kgetField(cdata, field, maxf, nf)
!          if cdata containes Tab, HP fortran
!             fails.
       implicit none
       character*(*) cdata
       character*(*) field(*)
       integer maxf  ! input. max of number of fields to be obtained
       integer nf    ! output. number of fields obtaiend

       integer klena
       integer nc, m1, m2

       nc = klena(cdata) 
       m2 = 1
       nf = 0
       do while (m2  .le.  nc)
          call kgetCpos(cdata(m2:nc), m1)
          m1 = m1 + m2 -1
          call kgetBpos(cdata(m1:nc), m2)
          m2 = m1 + m2 -1
          if(nf .lt. maxf) then
             nf = nf + 1
             field(nf) = cdata(m1:m2-1)
          else
             goto 10
          endif
       enddo
 10    continue
       end
       subroutine kgetCpos(ichrs,  idx)
       implicit none
!         shift ichrs to left until next non blank and
!         set its position in idx. 
       integer idx  ! = all + 1  if non blank is not found
       character*(*) ichrs  ! input characters
       integer k, i, klena
       character*1 tab

       tab = char(9)
       k = klena(ichrs)     ! effective length
       do i = 1, k
          if(ichrs(i:i) .ne. ' ' .and. ichrs(i:i) .ne. tab ) then
             idx = i
             goto 100
          endif
       enddo
       idx = k + 1
 100   continue
       end
       subroutine kgetBpos(ichrs, idx)
       implicit none
!         shift ichrs to left until next blank appear  and 
!         set its position in idx. 
       integer idx    ! = all+1 if blank is not found
       character*(*) ichrs
       integer k, i, klena

       character*1 tab

       tab = char(9)
       k = klena(ichrs)     ! effective length
       do i = 1, k
          if(ichrs(i:i) .eq. ' ' .or. ichrs(i:i) .eq. tab) then
             idx = i
             goto 100
          endif
       enddo
       idx = k + 1
 100   continue
       end
       subroutine kshiftB(ichrs,  ochrs, icon)
       implicit none
!        shift ichrs to left until next non blank and stores ochrs
       integer icon    ! = 0 ok.  1.  all blank. ochrs = ' '
       character*(*) ichrs, ochrs  ! ochrs can be ichrs
       integer k, i, klena
       character*1 tab

       tab = char(9)
       k = klena(ichrs)     ! effective length
       do i = 1, k
          if(ichrs(i:i) .ne. ' ' .and. ichrs(i:i) .ne. tab ) then

             ochrs = ichrs(i:k)
             icon = 0
             goto 100
          endif
       enddo
       ochrs=' '
       icon = 1
 100   continue
       end
       subroutine kshiftC(ichrs, ochrs, icon)
       implicit none
!      shift ichrs to left until next blank appear  and stores it into ochrs
       integer icon    ! = 0 ok.  1.  all non blank. ochrs = ' '
       character*(*) ichrs, ochrs  ! ochrs can be ichrs
       integer k, i, klena
       character*1 tab

       tab = char(9)
       k = klena(ichrs)     ! effective length
       do i = 1, k
          if(ichrs(i:i) .eq. ' ' .or. ichrs(i:i) .eq. tab) then
             ochrs = ichrs(i:k)
             icon = 0
             goto 100
          endif
       enddo
       ochrs = ' '
       icon = 1
 100   continue
       end
!      real*8 kDebye, x
!      integer i
!      x = 0.
!      do i = 1, 100
!        write(*, *) x, kDebye(4, x)
!         x = x + 0.1d0
!      enddo
!      end
      
      real*8 function kDebye(n, x)
      implicit none
      integer n
      real*8 x
!          compute  Debye func* n/x**n
!          Debye function is defined by
!          Int(0:x) t**n/(e**t-1)
!          x >= 0.

!
      real*8 kDebye1, kDebye2
      real*8 zeta(5)
      real*8 pi,  zeta1, term
      integer i
      parameter (pi = 3.14159265358979, zeta1=pi**2/6.d0)
      data zeta(2)/zeta1/
      data zeta(3)/ 1.20205690315959d0/
      data zeta(4)/ 1.08232323371114d0/
      data zeta(5)/ 1.03692775514337d0/
!
      if(x .lt. 2.5) then
         kDebye = kDebye1(n, x)
      else
         term = 1.
         do  i = 1, n
            term = term * i
         enddo
         kDebye = ( term* zeta(n+1) - kDebye2(n, x))*n/x**n
      endif
      end

      real*8 function kDebye2(n, x)
      implicit none
      integer n
      real*8 x
!          for large x

      integer k, kmax, i
      real*8 sum, eps,  term, term2
      data eps/1.d-14/, kmax/20/

      sum = 0.
      
      do k = 1, kmax
         term2 = x**n/k
         term = term2
         do i = 1, n
            term2 =(term2/x*(n-i+1))/k
            term = term +  term2
         enddo
         term =  term*exp(-k*x)
         sum = sum + term

         if(abs(term/sum) .le. eps) goto 10
      enddo
 10   continue
      kDebye2 = sum
      end
      real*8 function kDebye1(n, x)
      implicit none
      real*8 x  !
      integer n
! 
      real*8 kBern, fac, sum, eps, xp, term
      integer k, kmax, k2
      parameter (kmax = 17)
!        (2k)! = (2k-2)! * (2k-1) * 2k
      data eps/1.d-14/

      fac = 1.
      sum = 0.
      xp = 1
      if(x .ne. 0d0) then
         do k = 1, kmax 
            k2 = k*2
            fac = fac * (k2-1) * k2 
            xp = xp * x * x
            term = kBern(k2)/fac*xp/(k2+n)
            sum = sum + term
            if(abs(term/sum) .lt. eps) goto 10
         enddo
 10      continue
      endif
      kDebye1 = (1.d0/n - x/(2*(n+1)) +  sum)*n

      end
!      integer i
!      real*8  kBern
!
!      write(*,*) kBern(0),kBern(1)
!      do i = 2, 34, 2
!         write(*, *) i, kBern(i)
!      enddo
!      end
      
      real*8 function kBern(n)
      implicit none
      integer n  ! input. one of 0,1,2,4,6,8... 34
!          compute Bernoulli's constant Bn
!
      integer     m, x, m1
      character*80 msg
      parameter ( m = 34, m1= m/2 )
      real*8 NM(m1), D(m1)

       data NM/
!           2            4                 6               8
     1     1d0,        -1d0,              1d0,           -1d0,     
!          10            12                14             16         
     2     5d0,       -691d0,              7d0,       -3617d0,
!          18            20                22             24 
     3    43867d0,   -174611d0,       854513d0,   -236364091d0,
!          26            28                30             32
     4 8553103d0, -23749461029d0, 8615841276005d0, -7709321041217d0,
!          34
     5 2577687858367d0/

       data D/
!         2         4       6     8       10     12
     1   6d0,    30d0,    42d0, 30d0,   66d0, 2730d0,
!        14        16       18    20      22     24
     2   6d0,   510d0,   798d0, 330d0,  138d0, 2730d0,
!        26        28       30    32      34  
     3   6d0,   870d0, 14322d0, 510d0,  6d0/        


           if(n .eq. 0) then
              kBern = 1.0
           elseif(n .eq.  1) then
              kBern = -5.d0
           else
              x = n/2
              if(x*2 .ne. n) then
                 write(msg, *) ' n=',n, ' to kBern must be even'
                 call cerrorMsg(msg, 0)
              elseif(x .gt. m1) then
                 write(msg, *) ' kBern: n=',n, 'must be <=', m
                 call cerrorMsg(msg, 0)
              endif
              kBern = NM(x)/D(x)
           endif
         end
!      real*8 x, kSpence, y
!      integer i
!      x = -10.
!      do i = 1, 200
!         y = kSpence(x)
!         write(*, *) x, 1-x, y
!         x = x + 0.1d0
!      enddo
!      end

      real*8 function kSpence(x)
      implicit none
      real*8 x   ! input
!
!        compute Spence function defined by sum of k=1,inf [x**k/k**2].
!        Let's denote it by P(x).  
! 
!        P(x) = x + x**2/4 + x**3/9 + x**4/16 + ...  |x| <= 1
!     
!        The function f(x) in p.1004 of
!        H.M.F is
!           f(x) = P(1-x) or P(x) = f(1-x)
!           P(x) + P(-x) = P(x**2)/2
!           P(1-x) + P(1-1/x) = -(log(x))**2/2
!           P(-x)  - P(1-x) =- log(x)log(1+x) - pi**2/12 - P(1-x**2)/2
!           P(x)  = -log(x)**2/2 + pi**2/3 - P(1/x)  (x > 1)
!           P(x)= -log(|x|)**2/2 - pi**2/6 - P(1/x)  (x< -1)
!  if x is near 1.0,  convergence of the series, P(x),  is slow.
!  So we use the Debye function.  
!       P(x) = f(1-x) and f(x) =D(t) with t = - log(x)
!       so P(x) =  D(t), t= - log(1-x).
!       When x = 1, P(x) = pi**2/6 = zeta(2)
!
!     
      real*8  kDebye, y, p1, p2, temp, t

      real*8 pi,  zeta1
      parameter (pi = 3.14159265358979, zeta1=pi**2/6.d0)

      if(x .le. 0.) then
!          P(x) = P(x**2)/2 - P(-x)
           y = x**2
           if(y .eq. 1.)  then
              p1 = zeta1
           elseif(y .lt. 1.) then
              t =  -  log(1.-y)
              p1 =  kDebye(1, t)*t
           else
!             p(y) =  -p(1/y) ...
              t = - log(1. - 1./y)
              temp = kDebye(1, t)*t  ! p(1/y)
              p1 = -log(y)**2/2 + pi**2/3 -temp
           endif
           y = -x
           if(y .eq. 1) then
              p2 = zeta1
           elseif(y .lt. 1.) then
              t =  -  log(1.-y)
              p2 =  kDebye(1, t)*t
           else
!             p(y) =  -p(1/y) ...
              t = - log(1. - 1./y)
              temp = kDebye(1, t)*t  ! p(1/y)
              p2 = -log(y)**2/2 + pi**2/3 -temp
           endif
!          P(x) = P(x**2)/2 - P(-x)
           kSpence = p1/2 - p2
        else
           if(x .eq. 1.)  then
              kSpence = zeta1
           elseif(x .lt. 1.) then
              t =  -  log(1.-x)
              kSpence =  kDebye(1, t)*t
           else
!             p(x) =  -p(1/x) ...
              t = - log(1. - 1./x)
              temp = kDebye(1, t)*t  ! p(1/y)
              kSpence = -log(x)**2/2 + pi**2/3 -temp
           endif
        endif
        end
      subroutine  kfindField(buf, leng, loc1, loc2,  cond)
      integer leng              ! input
      character*(leng)  buf     ! input
      integer  loc1             ! input.  
      integer  loc2             ! output  last non blank pos. in buf
                                !  after loc1 in buf.
!                                  if loc2 cannot be found loc2 = 0
      integer cond              !  input. if ==0 nothing happens
                                !   if ==1 and loc2==0; message is given
                                !   but execution contines
                                !   if  == 2 and loc2==0,   message is given
                                !  and execution stops

      integer i

      i = loc1
      loc2 = 0
      do while ( i .lt. leng )
         if( buf(i:i) .ne. ' '  .and. buf(i+1:i+1) .eq. ' ') then
            loc2=i
            goto 10
         endif
         i = i + 1
      enddo
      if(buf(leng:leng) .ne.  ' ' ) then
         loc2 = leng
         goto 10
      else
         loc2 = 0
      endif
 10   continue
      if( cond .gt.  0  .and. loc2 .eq. 0 ) then
         write(0, *) ' In kfindField: buf=',buf
         write(0, *) ' cannot find non blank after', loc1
         if(cond .eq. 2) stop 99999
      endif

      end
!c     test kxplineseg
!      real*8 x1, y1,x2,y2, x3, y3, x4, y4, eps
!      real*8 x,y, p, q
!      integer icon
!      eps  = 1.d-12
!      read(*,*) x1,y1
!      read(*,*) x2,y2
!      read(*,*) x3,y3
!      read(*,*) x4,y4
!
!      call kxplineseg(x1, y1, x2, y2, x3, y3, x4, y4, eps,
!     *     x, y, p, q, icon)
!      write(*,*) ' icon =', icon, 'p,q=',p,q
!      write(*,*) ' x,y=',x,y
!      end
!
!
!     get x-ssing point of two line segments on (x,y) plane
!
      subroutine kxplineseg(x1,y1, x2, y2, x3, y3, x4, y4, eps,
     *     x, y, p, q, icon)

      implicit none
      real*8 x1, y1, x2, y2  ! input one line segment. in complex ( z1, z2)
      real*8 x3, y3, x4, y4  !  input the other line segment. in complex  (z3, z4)
      real*8  eps            !  input. to judge the parallelity and/or
                             !         overlapping of two line segment
      real*8  x, y          ! output. obtained x-ssing point
      real*8  p             ! output. (z1 to the x-ppint)/(z1 to z2).  0<=p<=1. if not  outside
      real*8  q             ! output. (z3 to the x-point)/(z3 to z4).  0<=q<=1. if not  outside
      integer icon          ! output.  0, x,y obtained.
                            !          1, two segment overlap, (x,y) will 
                            !             be somewhere on the line
                            !          2, no x-ssing point 
      complex*16  expia, z1, z2, z3, z4
      real*8 cosa, sina
      real*8  length
      
      z1 = cmplx(x1, y1, 8)
      z2 = cmplx(x2, y2, 8)
      z3 = cmplx(x3, y3, 8)
      z4 = cmplx(x4, y4, 8)
      length = abs(z4-z3)
      if( length .le. eps )  then
         icon = 2
      else
         cosa = (x4-x3)/length
         sina = (y4-y3)/length
         expia = cmplx(cosa, sina,8)
         call kxplsl(z3, expia, z1 ,z2, eps, p, q, icon)
         if(icon  .eq. 0) then
!             see if x-point is within z3-z4.
            if(q .lt. 0.  .or. q .gt. length ) then
!                point is outside of the segment
               icon = 2
            else
               q = q/length
            endif
         elseif(icon .eq. 1 ) then
            if( (x3-x1)*(x2-x3) .gt. 0.) then
               x = x3
               y = y3
            elseif( (x1-x3)*(x4-x1) .gt. 0.) then 
               if( (x2-x3)*(x4-x2) .gt. 0.) then 
                  if(abs(x1-x3) .lt. abs(x2-x3) ) then
                     x = x1
                     y = y1
                  else
                     x = x2
                     y = y2
                  endif
               else
                  x = x1
                  y = y1
               endif
            elseif( (x2-x3)*(x4-x2) .gt. 0.) then 
               x = x2
               y = y2
            else
               icon = 2 
            endif
         endif
      endif
      if(icon .eq. 0) then
         x = x1 +p*(x2-x1)
         y = y1 +p*(y2-y1)
      endif 
      end
!           test
!      implicit none
!      real*8 a, b, x0, y0, dirx, diry, x1, x2, y1, y2,l1,l2
!      integer cross
!      a = 2.
!      b = 1.
!      x0 = 0.5
!      y0 = -3.0
!      dirx = 0.1d0
!      diry = sqrt(1.d0-dirx**2)
!      call kxplellip(a, b, x0, y0, dirx, diry, l1, l2, cross)
!      write(0,*) ' cross=',cross, x0+l1*dirx,  y0+l1*diry,
!     *          x0+l2*dirx, y0+l2*diry
!      end
!
!
      subroutine kxplellip(a, b, x0, y0, dirx, diry, l1, l2, cross)
!         get crossing point of a line with an ellipse whose center is at (0,0)
!
      implicit none
      real*8 a  ! input   x radis of the ellipse
      real*8 b  ! input   y //
      real*8 x0 ! input   the line start from this x
      real*8 y0 ! input   the line start from this y
      real*8 dirx  ! input the line's direction cos to x
      real*8 diry  ! input the line's //               y
      real*8 l1    !  output.  a cross point ; x0+l1*dirx, y0+l1*diry
      real*8 l2    !  output.  the other cross point
                   !      0 <l1<l2 or l2 <0<l1 or l2<l1<0. or l1=l2
      integer cross ! output.  0--> the line cross the ellipse at two point
                    !          1--> the line is tangnetial to the ellipse
                    !          if two point are too near,  we regards them as one.
                    !          -1--> the line dose not cross the ellipse
      real*8   x, y, wx, wy, c0, c1h, c2, dq
      real*8   small, temp
      data small/1.d-8/
      

      x = x0/a
      y = y0/b
      wx = dirx/a
      wy = diry/b
      c2= wx**2 + wy**2
      c1h = x*wx+ y*wy
      c0 = x**2 + y**2 - 1.d0

      dq = c1h**2 - c0*c2
      
      if( dq .lt. 0.)  then
         cross = -1
      elseif( dq .le. small) then
         cross = 1
      else
         cross = 0
      endif
      if( cross .ge. 0 ) then
         dq =sqrt(dq)
         l1 = (-c1h - dq)/c2
         l2 = (-c1h + dq)/c2
         if(l1*l2 .lt. 0.) then
            if(l1 .lt. 0.) then
               temp = l1
               l1 = l2
               l2 = temp
            endif
         else
            if(abs(l1) .gt. abs(l2)) then
               temp = l1
               l1 = l2
               l2 = temp
            endif
         endif
      endif
      end
!      implicit none
!      real*8  x, y, z, el, l, m, n
!      real*8  a, b, c, h 
!      integer  icon, cond 
!cc                    
!c
!cc      write(0, *) ' enter x,y,z, a, b, c, h'
!cc      read(*, *)  x,y,z, a, b, c, h 
!       x = 19.3364132663471d0       
!       y = 0.351801084456073d0    
!       z =  -14.7246124670421d0 
!       a = 28.2000000000000d0      
!       b= 1.00000000000000d0     
!       c = 28.2000000000000d0      
!       h =      -28.2000000000000d0
!c
!      call kioPrism(x,y, z, a, b, c, h, icon)
!      write(*,*) icon
!c
!      l= 0.730182332808807d0 
!      m =  1.901783791699804d-002
!      n = -0.682987616600244d0     
!c
!      call kxplPrism(x, y, z, l, m, n, 
!     *            a, b, c, h , el, cond)
!cc         get the length to the nearest crossing point
!      write(0,*) el, cond
!      end
      subroutine kioPrism(x, y, z, a, b, c, h, icon)
!       judge if (x,y,z) is inside of a given prism
!       or not.
!
      implicit none
      real*8  x, y, z  ! input. a given point
      real*8  a, b, c, h ! input.  const to characterize a prism
      integer  icon  !  output.  icon =0.  the point is inside or
                     !            on the surface
                     !                +1   outside.

!          
!           | z
!           |           *
!           |         * |    + 
!           |       *   |         + 
!           |     *     h            +
!           |   *       |                 +
!           | *         |                      +
!           |___________|__________________________+   x
!           0           c                          a
!           
!          b is the depth on the y axis.
!  a, b,  c,  h can be negative.
!
!      
!
      integer jcon

!         see if y is inbetween 0,b
      if( (b-y)*y .ge.  0.d0) then
!           see if (x,z) is in the triagle
         call k3inout(x, z, 0.d0, 0.d0, a, 0.d0, c, h, jcon)

         if(jcon .eq. 0)  then
            icon = 0
         else
            icon = 1
         endif
      else
         icon = 1
      endif
      end
      subroutine kxplPrism(x0, y0, z0, l, m, n, 
     *            a, b, c, h , el, cond)
!         get the length to the nearest crossing point
!    of a given line with a given prism.
!
      real*8 x0, y0, z0  ! input.  the line passes this point
      real*8 l, m, n     ! input.  the line's direction cosine
      real*8 a, b, c, h  ! input.  prism parameter
      real*8 el          ! output. length to the crossing point
                         !         (only el>=0 is obtained)
      integer cond       ! output. =0, the point is inside or on 
                         !             the surface. el was obtained.
                         !         =1, the point is outside.
                         !             el was obtained.
!                        !         =-1  no crossing point.  el undef.
!     Z
!     |
!     |        *
!     |       *|    +
!     |      * |        + 
!     |     *  |            +
!     |    *   |                +
!     |   *    h                    + 
!     |  *     |                        +  
!     |________|___________________________+  X
!     0        c                           a 
!         b is depth along y.
!
!      a,b,c,h can be negative.
!
      real*8 x, y, z, el1, temp, eps
      integer inout
      character*120 msg


      data eps/-1.d-12/


!

!
!        x:  0~a is not necesarrily the range but x1~x2.
      x1 = min(0.d0, a, c)
      x2 = max(0.d0, a, c)
!         judge if (x0,.)  is inside
      call kioPrism(x0, y0, z0, a, b, c, h, cond)

      if(cond .eq. 0) then
!           inside
!          cross with the x-y plane ?
         if(n .ne.  0.d0) then
            el = -z0/n
         elseif(z0 .eq. 0.d0) then
            el = 0.
         else
            el = -1.    !  to skip next if
         endif
         if(el .ge. 0.d0) then
            x = x0 + el* l
            y = y0 + el* m
!             if (x,y) is inside the bottom, el obtained
            if( x*(a-x) .ge. eps  .and.
     *          (b-y)*y .ge. eps ) goto 10
         endif
!              cross with  / ?
!             (x-x0)/l = (z-z0)/n;  at 
         temp = c*n-h*l
         if(temp .ne. 0.d0) then
            el = (h*x0 - c*z0)/temp
         elseif( (h*x0 - c*z0) .eq. 0.d0 ) then
            el = 0.
         else
            el = -1.  ! this is to skip next if
         endif
         if(el .ge. 0.d0) then
            z = z0 + el * n
!               z : inbetween 0~h
            if( (h-z)*z .ge. eps) then  !!
               x = x0 + el * l
               if( x* (c - x) .ge. eps) then
                  y = y0 + el * m
                  if( (b-y)*y .ge. eps) goto 10   !!
               endif
            endif
         endif
!          cross with \ ?
         temp = h*l-(c-a)*n
         if(temp .ne. 0.d0) then
            el =(a*h- h*x0 + z0*(c-a))/temp
         elseif((a*h- h*x0 + z0*(c-a)) .eq. 0.d0) then  
            el = 0.
         else
            el = -1.  ! to skip next
         endif
         if(el .ge. 0.d0) then
            z = z0 + el * n
            if( (h-z)*z .ge. eps) then
               x = x0 + el * l
               if( (c - x)* (x - a) .ge. eps) then   !!
                  y = y0 + el * m
                  if( (b-y)*y .ge. eps) goto 10
               endif
            endif
         endif
!           cross with x-z plane at y=0
         if(m .ne. 0.d0) then
            el = - y0/m
         elseif(y0 .eq. 0.d0) then
            el = 0.
         else
            el = -1.    ! to skip next if
         endif
         if(el .ge. 0.d0) then
            x = x0 + el * l
            z = z0 + el * n
            call k3inout(x, z, 0.d0, 0.d0, a, 0.d0, c, h, inout)
            if(inout .eq. 0) goto 10
         endif
!             cross with x-z plone at y=b  ?
         if(m .ne. 0.d0) then
            el = (b-y0)/m
         elseif(b .eq. y0)  then
            el = 0.
         else
            el  = -1.  ! to skip next if
         endif
         if(el .ge. 0.d0) then
            x = x0 + el * l
            z = z0 + el * n
            call k3inout(x, z, 0.d0, 0.d0, a, 0.d0, c, h, inout)
            if(inout .eq. 0) goto 10
         endif
!           error
         call  cerrorMsg(
     *    'kxplPrism; point is inside but no x%p',1)
         write(0, *) 'x0,y0,z0=',x0,y0,z0, ' lmn=',l,m,n
         write(0, *) 'a, b, c, h =',  a,b,c,h
         stop
      else
!        ***********************************************
!         point is outside. firstly, judge no crossing
!         possiblity quickly
!
         cond = - 1

!            see if crossing point cannot be in the xrange.
!          For X point  to eixt inbetween x1,x2
!          f(el)=  (el*l - (x1-x0))(el*l -(x2-x0)) <= 0
!          must be sutisfied for some el >=0.
!          If there is no el >=0, there is no possibility of
!          crossing.
!          That  condtions is:
!                 for f(0) > 0 and  (x1-x0)/l < 0
!             (l !=0) .  if l=0,  f(0) > 0 is enough
         if((x1-x0)*(x2-x0) .gt. 0.d0 .and.
     *         (x1-x0)*l .le. 0.d0) goto 10
!            see if crossing point cannot be in the yrange
         if( (b-y0)*(-y0) .gt. 0.d0 .and.
     *       (b-y0)*m .le. 0.d0) goto 10
!            see if crossing point cannot be in the zrange
         if( (h -z0)*(-z0) .gt. 0.d0  .and.
     *              (h -z0)*n .le.0.d0) goto 10
!
!             there is a possibilty of crossing
         
         nc = 0          ! crossing point counter
!          see x.p on (x-y) 
         if(n .ne. 0.d0) then
            el = - z0/n
            if(el .ge. 0.d0) then
               x  = x0 + el *l
               y  = y0 + el *m
               if( x *(a-x) .ge. eps ) then
                  if( (b-y)* y .ge. eps ) then
                     nc = 1
                     el1 = el
                  endif
               endif
            endif
!        else
!            this case should be inside if x.p exists
         endif
!           see x.p on /
         temp = c*n - h*l
         if(temp .ne. 0.) then
            el =  (h*x0 - c*z0)/temp
            if(el .ge. 0.d0) then
               x = x0 + el*l
               y = y0 + el*m
               z = z0 + el*n
               if( x*(c-x) .ge. eps .and.
     *              y*(b-y) .ge.  eps  .and.
     *              z*(h-z) .ge.  eps ) then 

                  if(nc .eq. 0) then
                     el1 = el
                     nc = 1
                  else
                     el = min(el, el1)
                     cond = 1
                     goto 10
                  endif
               endif
            endif
!         else
!               this should not happen--> inside case
         endif
!
!              \ case
         temp = h*l -(c-a)*n
         if(temp .ne. 0.) then
            el =( a*h - h*x0 + (c-a)*z0 )/temp

            if(el .ge. 0.d0) then
               x = x0 + el*l
               y = y0 + el*m
               z = z0 + el*n
               if( (x-c)*(a-x) .ge. eps .and.
     *              y*(b-y) .ge.  eps  .and.
     *              z*(h-z) .ge.  eps ) then 

                  if(nc .eq. 0) then
                     nc = 1
                     el1 = el
                  else
                     el = min(el, el1)
                     cond = 1
                     goto 10
                  endif
               endif
            endif
!         else
!           this case; no need be checked.
         endif
!            cross with  x-z at y=0
         if(m .ne. 0.d0) then
            el = - y0/m
            if(el .ge. 0.d0) then
               x = x0 + el *l
               z = z0 + el *n
               call k3inout(x, z, 0.d0, 0.d0, a, 0.d0,
     *                     c, h, inout)
               if(inout .eq. 0) then
                  if(nc .eq. 0) then
                     nc = 1
                     el1 = el
                  else
                     el = min(el1, el)
                     cond = 1
                     goto 10
                  endif
               endif
            endif
!        else
!           not coming here 
         endif
!          cross with x-z y=b
         if(m .ne.  0.d0) then
            el = (b-y0)/m
            if(el .ge. 0.d0) then
               x = x0 + el *l
               z = z0 + el *n
               call k3inout(x, z, 0.d0, 0.d0,
     *         a, 0.d0, c, h, inout)
               if(nc .eq. 0) then
                  nc = 1
                  el1 =  el
               else
                  el = min(el1, el)
                  cond = 1
                  goto 10
               endif
            endif
!        else
!           not coming here 
         endif
         cond = -1
      endif
 10   continue
      end
