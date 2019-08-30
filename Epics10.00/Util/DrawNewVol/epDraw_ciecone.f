#if defined NEXT486
#define IMAG_P dimag
#elif defined PCLinux
#define IMAG_P dimag
#else
#define IMAG_P imag
#endif
      
!     *****************************
      subroutine epDraw_ciecone(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   a ciecone in local coordnate.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >= (nvccl+2)*2
      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.  


      include "../../prog/NewVol/Zciecone.h"

      real*8 ra, rb, h, xt, rap, r, sa, ea
      real*8 dirx, diry, sab, eab, l1, l2, sat, eat,
     *        x1, y1, x2, y2, x, y
       type(epPos)::  bv(100), tv(100)
      integer crossb, crosst, rn, jb, jt, i
      logical ok

!       real*8 t1, t2
       type(epPos)::  mpb, mpt
      real*8 xs, ys, xe,  ye, xsc, ysc, xec, yec,
     * cosa, sina
      integer n1, n2, nsvb, nsvt, n1p, n2p, nsv
      integer nvt, nvb
      logical kdgtest
!      logical isinside
!      logical isinside2, scut, ecut
      integer icon
      real*8 pp, q, eps, ang
      complex*16 zz0, zz1, zz2, expa, zz

                

      real*8 f 
      real*8 epsilon, small
!             see if  thetamin < x < thetamax
!      isinside(x) = mod(thetamin-thetamax+360.d0, 360.d0) .ge.
!     *               mod(x-thetamax+360.d0, 360.d0)
!             see if   sa <  x < ea
!      isinside2(x) = mod(ea-sa+360.d0, 360.d0) .ge.
!     *               mod(x-sa+360.d0, 360.d0)
      f(x,y) = (ye-ys)*(x-xs) - (xe-xs)*(y-ys)

      data epsilon/1.d-8/, small/-1.d-3/
      save epsilon, small
      logical stupid
!//////////
      integer:: iii
!/////////


      ra = Volat( comp%vol + ira)
      rb = Volat( comp%vol + irb)
      rap = Volat( comp%vol + irap)
      h = Volat( comp%vol + ih)
      xt = Volat( comp%vol+ ixt)
      xs =  Volat(comp%vol+ixs)
      ys =  Volat(comp%vol+iys)
      xe =  Volat(comp%vol+ixe)
      ye =  Volat(comp%vol+iye)
      eps = Volat(comp%vol+ieps)
      ang = atan2(ye-ys, xe-xs)
      dirx = cos(ang)
      diry = sin(ang) 

      if(ra .ge. rap) then
         crossb = 0
         mpb%x = (xs+xe)/2.
         mpb%y = (ys+ye)/2.
!            lower ellipse is larger; sa,ea defined there
         sab= Volat( comp%vol + isa )
         eab= Volat( comp%vol + iea )
         sa = sab
         ea = eab
!            see if cut line crosses the upper ellips
         call kxplellip(rap, rap*eps, xs-xt, ys, dirx, diry, 
     *        l1, l2, crosst )
         if(crosst .eq. 0) then
            x1 = (xs-xt) + l1*dirx
            x2 = (xs-xt) + l2*dirx
            y1 = ys + l1*diry
            y2 = ys + l2*diry
            sat = atan2(y1, x1)*Todeg
            if(sat .lt. 0.) sat = sat+ 360.
            eat = atan2(y2, x2)*Todeg
            if(eat .lt. 0.) eat = eat + 360.
!              point corresponding to (xs,ys), (xe,ye)  at upper ellipse
            cosa = cos(sat*Torad)
            sina = sin(sat*Torad)
            r = 1.0d0/sqrt( (cosa/rap)**2 + (sina/rap*eps)**2 )
            xsc = r*cosa + xt
            ysc = r*sina
            cosa = cos(eat*Torad)
            sina = sin(eat*Torad)
            r = 1.0d0/sqrt( (cosa/rap)**2 + (sina/rap*eps)**2 )
            xec = r*cosa + xt
            yec = r*sina
            mpt = mpb
         else
            sat=360.
            eat= 0.
            mpt%x = xt
            mpt%y = 0. 
         endif
      else
         crosst=0
!           upper ellipse is larger . angle  is defined there
         mpt%x = (xs + xe)/2.
         mpt%y = (ys + ye)/2. 
         sat= Volat( comp%vol + isa)
         eat= Volat( comp%vol + iea)
         sa = sat
         ea = eat
!            see if cut line crosses the lower ellips
         call kxplellip(ra, rb, xs, ys, dirx, diry, 
     *        l1, l2, crossb )
         if(crossb .eq. 0) then
            x1 = xs + l1*dirx
            x2 = xs + l2*dirx
            y1 = ys + l1*diry
            y2 = ys + l2*diry
            mpb = mpt
            sab = atan2(y1, x1)*Todeg
            eab = atan2(y2, x2)*Todeg
!           point   corresponding to (xs,ys), (xe,ye)  at upper ellipse
            cosa = cos(sab*Torad)
            sina = sin(sab*Torad)
            r = 1.0d0/sqrt( (cosa/ra)**2 + (sina/rb)**2 )
            xsc = r*cosa 
            ysc = r*sina
            cosa = cos(eab*Torad)
            sina = sin(eab*Torad)
            r = 1.0d0/sqrt( (cosa/ra)**2 + (sina/rb)**2 )
            xec = r*cosa 
            yec = r*sina
         else
            sab=360.
            eab= 0.
            mpb%x = 0.
            mpb%y = 0.
         endif
      endif
       
      if(sa .lt. ea)  sa = sa+360.
      n = 0
!        draw lower ellipse; full ellipse vertex
      call epdrawElps(ra, rb, 0.d0,  ea, sa,  p(n+1), n1)
      n1 = n1-1
      n = n + n1
      if(ea .lt. sa) ea = ea+360.
      call epdrawElps(ra, rb, 0.d0,  sa, ea,  p(n+1), n1p)
      n1p = n1p-1
      n1 = n1 + n1p
      n = n + n1p

      nsvb = 1
      nsvt = n+1
!        draw top ellipse; 
      sa = mod(sa, 360.d0)
      ea = mod(ea, 360.d0)
      if(sa .lt. ea) sa = sa+360.
      call epdrawElps(rap, rap*eps, h, ea,  sa,   p(n+1), n2)
      n2 = n2- 1
      n = n + n2
      if(ea .lt. sa) ea = ea+360.
      call epdrawElps(rap, rap*eps, h, sa, ea,  p(n+1), n2p)
      n2p = n2p-1
      n2 = n2+ n2p
      n  = n + n2p

      do i = nsvt, nsvt+n2-1
         p(i)%x = p(i)%x+ xt
      enddo

!       n = n + 1
!       p(n).x = gpsep
!         at this level, n1=n2

       rn = 0
       do i = 1, n1
          jb = i + nsvb-1
          jt = i + nsvt-1

          if( f( p(jb)%x, p(jb)%y ) .gt. small ) then
             if(  f(p(jt)%x, p(jt)%y) .gt. small ) then
!              bot and top vertex are in the region 
                rn  = rn + 1
                bv(rn) = p(jb)
                tv(rn) = p(jt)
             else
!                 top has cut; get crossing point 
                zz0 = cmplx( xs, ys, 8 )
                expa= cmplx(dirx, diry, 8)
                zz1 = cmplx(p(jb)%x, p(jb)%y, 8)
                zz2 = cmplx(p(jt)%x, p(jt)%y, 8)

                call kxplsl(zz0, expa, zz1, zz2, epsilon, 
     *           pp, q, icon)
                if( icon .eq. 0 .or. icon .eq. 2 ) then
!                    x-point exists
                   x = xs + q*dirx
                   y = ys + q*diry

                   rn  = rn + 1
                   bv(rn) = p(jb)
                   tv(rn)%x = x
                   tv(rn)%y = y
                   tv(rn)%z = h*pp
                endif
             endif
          elseif( f(p(jt)%x, p(jt)%y) .gt. small ) then
!                 bottom is cut
             zz0 = cmplx( xs, ys, 8 )
             expa= cmplx(dirx, diry, 8)
             zz1 = cmplx(p(jb)%x, p(jb)%y, 8)
             zz2 = cmplx(p(jt)%x, p(jt)%y, 8)
             call kxplsl(zz0, expa, zz1, zz2, epsilon,
     *         pp, q, icon)
             if( icon .eq. 0 .or. icon .eq. 2 ) then
!                x-point exists
                x = xs + q*dirx
                y = ys + q*diry
                rn  = rn + 1
                tv(rn) = p(jt)
                bv(rn)%x = x
                bv(rn)%y = y
                bv(rn)%z = h*pp
             endif
!             else  ! both cut so don't need to draw
          endif
       enddo
       nsvt = rn+2

       do i = 1, rn
          p(i) = bv(i) 
          p(i+rn+1)=tv(i)
       enddo
       p(rn+1)%x = gpsep
       p(rn*2+2)%x = gpsep
       p(rn*2+3)%x = gpsep
       n = rn*2+3

       if(  kdgtest(how,3) ) then
!         cut vertical plane
          nvb = 0
          nvt = 0
          do i = 1, rn
             if( bv(i)%z .gt. 0. ) then
                nvb = nvb + 1
             endif
             if( tv(i)%z .lt. h ) then
                nvt = nvt + 1
             endif
          enddo

          if(nvb .gt. nvt) then
             nvb=nvb+2
             nsv = n
             p(n+nvb+1)%x = gpsep
             do i =1, rn
                ok=.false.
                if( bv(i)%z .gt. 0 )  then
                   ok = .true.
                elseif( bv(i)%z .eq. 0. ) then
                   if( i .lt. rn ) then
                      ok = bv(i+1)%z .gt. 0.
                   endif
                   if(.not. ok . and.  i .gt. 1 ) then
                      ok = bv(i-1)%z .gt. 0.
                   endif
                endif
                if(ok) then
                   n = n+1
                   p(n) = bv(i)
                   p(n)%z =0. 
                   p(n+nvb+1) = bv(i)
                endif
             enddo
             n = nsv + 2*nvb + 2 
             p(n)%x =  gpsep
             n = n+1
             p(n)%x =  gpsep
          elseif(nvt .gt. nvb) then
             nvt = nvt+2
             nsv = n
             p(n+nvt+1)%x = gpsep
             do i =1, rn
                if(i>1 ) then
                   stupid =( ( tv(i)%z .eq. h .and. tv(i+1)%z .lt. h) 
     *                  .or. tv(i)%z .lt. h  .or.
     *                  (tv(i)%z .eq. h .and. tv(i-1)%z .lt. h) 
     *                  )
                else
                   stupid =( ( tv(i)%z .eq. h .and. tv(i+1)%z .lt. h) 
     *                  .or. tv(i)%z .lt. h  )
                endif
                if(stupid) then
                   n = n+1
                   p(n) = tv(i)
                   p(n+nvt+1) = tv(i)
                   p(n+nvt+1)%z = h
                endif
             enddo
             n = nsv + 2*nvt + 2 
             p(n)%x =  gpsep
             n = n+1
             p(n)%x =  gpsep
          endif
       endif
!
       if(kdgtest(howcyl, 1)) then
!//////////
          if(crossb == 0 ) then
             mpb%x =( p(1)%x + p(rn)%x )/2
             mpb%y =( p(1)%y + p(rn)%y )/2
          endif
!//////////
          call epdrawCcylEdg( p(1), rn+1, 0.d0,  mpb, p(n+1), n2)
          n = n + n2 
       endif
       if(kdgtest(howcyl, 2)) then
!           nr+1; includes gpsep
!//////////
          if(crosst == 0 ) then
             mpt%x =( p(nsvt)%x + p(nsvt+rn-1)%x ) /2
             mpt%y =( p(nsvt)%y + p(nsvt+rn-1)%y )/2
          endif
!//////////
          call epdrawCcylEdg( p(nsvt), rn+1, h,  mpt, p(n+1), n2)
          n = n + n2
          n = n + 1
          p(n)%x = gpsep
       endif
      end


         
