!      Compute magnetic deflection (angle and displacement)
!      Use Runge-Kutta-Fehlberg Method.
!      
!     ************************************************  
!>@brief Compute magnetic deflection (angle and displacement)
!> Use Runge-Kutta-Fehlberg Method.
!>@param[in] aTrack a charged ptcl (datatype /track/)
!>@param[in] leng length travelled in m (datatype real*8)
!>@param[out] newpos new position after leng move (datatype /coord/newpos)
!>@param[out] newmom changed mon 
      subroutine cbDefByRK(aTrack, leng, newpos, newdir,newmom)
!     ************************************************  
!   This uses Runge-Kutta method with the automatically adjustable
!   step size and takes a long time for execution.
!

 
      implicit none

#include  "Ztrack.h"

      real(8),intent(out):: newmom(3) ! chnaged mom.
      real*8 geomcnst  !  Z/pc/3.3358. Z charge, pc GeV
      common /Zgeomcnst/ geomcnst

      type(track)::aTrack       !input. a charged ptcl

!                              ptcl and magfiled coord is in Exyz

      real*8  leng                !input. length  travelled in m

      type(coord)::newpos      ! output. new position after leng move
      type(coord)::newdir      ! //      new direction cose. //

      real*8  cbetat, x0

      integer first, nde, ierr
      real*8  relaerr, abserr, h
      real*8  esterr(6), yn(6), y10(6), y1n(6), err(6), y0(6)
      real*8  work(6,12)
      save  abserr, relaerr
      
      external cgeomfunc

      data first/1/, relaerr/1.d-2/, abserr/1.d3/

      cbetat= leng
!       
!      y0(1) = aTrack.pos.xyz.x
!      y0(2) = aTrack.pos.xyz.y
!      y0(3) = aTrack.pos.xyz.z
!           IBM must be like below
      y0(1) = aTrack%pos%xyz%r(1)
      y0(2) = aTrack%pos%xyz%r(2)
      y0(3) = aTrack%pos%xyz%r(3)
      
!      y0(4)=  aTrack%vec%w%x
!      y0(5)=  aTrack%vec%w%y
!　    y0(6)=  aTrack%vec%w%z
      y0(4:6) = aTrack%vec%w%r(1:3)

!         pc in GeV
      geomcnst = max(sqrt(aTrack%p%fm%p(4)**2- aTrack%p%mass**2), 1.d-9)

      geomcnst = aTrack%p%charge/geomcnst/3.3358

      x0= 0.d0   ! this cannot be a literal below; due to bad prog.
      call krungekutfds(6, cgeomfunc, x0, cbetat, y0, first, relaerr,
     * abserr, yn, esterr, nde, ierr, h, y10, y1n, err, work)
      if(ierr .eq. 40000 ) then
         call  cbDefByRK2(aTrack, leng, newpos, newdir, newmom)
      elseif(ierr .ne. 0) then
         write(0,*) ' ierr=',ierr
         call  cbDefByRK2(aTrack, leng, newpos, newdir, newmom)
      else

!c      first = 1

!      newpos.x = yn(1)
!      newpos.y = yn(2)
!      newpos.z = yn(3)
!      for  ibm 
         newpos%r(1) = yn(1)
         newpos%r(2) = yn(2)
         newpos%r(3) = yn(3)
         newpos%sys='xyz'
!      newdir.x = yn(4)
!      newdir.y = yn(5)
!      newdir.z = yn(6)
!      for ibm
         newdir%r(1) = yn(4)
         newdir%r(2) = yn(5)
         newdir%r(3) = yn(6)
!               oh my god next yn(4:6)* sqrt(. ) was
!                              yn(4:6)** sqrt(...)   !!!
!             but default UseRungeKutta = 0 did not come here
         newmom(:) = yn(4:6)* 
     *   sqrt(dot_product( aTrack%p%fm%p(1:3),aTrack%p%fm%p(1:3)))
         newdir%sys='xyz'
      endif
      end
      subroutine cgeomfunc(L, y, f)
      implicit none
!         this  is used both by cbDefByRK and cbDefByRK2
!      Zobs.h is For Zobsp.h
!      Zobsp.h is  For  YearOfGeomag
#include "Zcoord.h"
#include "Zmagfield.h"
#include "Zobs.h"  
#include "Zobsp.h" 

      real*8 geomcnst  !  Z/pc/3.3358. Z charge, pc GeV
      common /Zgeomcnst/ geomcnst

      type(coord)::llh
      type(magfield)::B
      integer icon
      real*8 L, y(6), f(6)
!        L=cbetat is not used explicitly

      f(1) = y(4)
      f(2) = y(5)
      f(3) = y(6)
!          compute B at y(1),y(2), y(3)
!
!      llh.x = y(1)
!      llh.y = y(2)
!      llh.z = y(3)
!        for ibm
      llh%r(1) = y(1)
      llh%r(2) = y(2)
      llh%r(3) = y(3)

      llh%sys = 'xyz'   ! converted to 'llh' in cgeomag
      call cgeomag(YearOfGeomag, llh, B, icon)  ! B is ned
      call ctransMagTo('xyz', llh, B, B)        ! to xyx
      f(4)= geomcnst * (y(5)*B%z - y(6)*B%y)
      f(5) =geomcnst * (y(6)*B%x - y(4)*B%z)
      f(6) =geomcnst * (y(4)*B%y - y(5)*B%x)
      end


      subroutine krungekutfds(neq,func,x0,xe,y0,init,relerr,abserr,
     &                 yn,esterr,nde,ier,h,yl0,yln,err,work)
!
!        adapted from Maruzen Library.
!
!*********************************************************************
!     subroutine krungekutfds numerically integrates a system of neq *
!     first order ordinary differential equations of the form        *
!             dy(i)/dx = f(x, y(1),..., y(neq)),                     *
!     by the runge-kutta-fehlberg (4,5) formula,                     *
!     subroutine rkfds can detect the stiffness of the differential  *
!     system.                                                        *
!                                                                    *
!     parameters                                                     *
!  === input ===                                                     *
!     (1) neq: number of equations to be integrated                  *
!     (2) func: subroutine func(x,y,f) to evaluate derivatives       *
!                f(i)=dy(i)/dx                                       *
!     (3) x0: initial value of independent variable                  *
!     (4) xe: output point at which the solution is desired          *
!     (5) y0(i) (i=1,..,neq): initial value at x0                    *
!     (6) init: indicator to initialize the code                     *
!          init=1..the code will be initialized(first call).         *
!          init=2..the code will not be initialized(subsequent call).*
!     (7) relerr: relative local error tolerance                     *
!     (8) abserr: absolute local error tolerance                     *
!  === output ===                                                    *
!     (9) yn(i) (i=1,..,neq): approximate solution at xe             *
!    (10) esterr(i) (i=1,..,neq): estimate of the global error in    *
!          the approximate solution yhn(i)                           *
!    (11) nde: number of derivative evaluations                      *
!    (12) ier: indicator for status of integration                   *
!          ier=0..integration reached xe. indicates successful       *
!             return and is the normal mode for continuing the       *
!             integration.                                           *
!          ier=10000..integration was not completed because too many *
!             derivatives evaluations were needed.  if the user wants*
!             to continue the integration, he just calls rkf again.  *
!          ier=20000..integration was not completed because error    *
!             tolerance was inappropriate(=it was zero).  must use   *
!             non-zero abserr to continue the integration.           *
!          ier=30000..integration was not completed because requested*
!             accuracy could not be achieved using smallest stepsize.*
!             must increase the error tolerance to continue the      *
!             integration.                                           *
!          ier=40000..integration was not completed because the      *
!             stiffness of the differential system was detected.     *
!             must use the subroutine, krsnat or krsat, which        *
!             is designed to integrate stiff differential systems.   *
!  === others ===                                                    *
!    (13) h: variable to hold information internal to rkfds which is *
!           necessary for subsequent calls.                          *
!    (14) yl0(), yln(): array (size=neq) to hold information internal*
!           to rkfds which is necessary for subsequent calls.        *
!    (15) err(): array (size=neq) to be used inside rkfds            *
!    (16) work(): two-dimentional array (size=(neq,12)) to be        *
!                 used inside rkfds                                  *
!    copyright: m. sugihara, november 15, 1989, v. 1                 *
!*********************************************************************
       implicit real*8(a-h,o-z)
       external func
       dimension y0(neq),yn(neq),esterr(neq),work(neq,12),
     &           yl0(neq),yln(neq),err(neq)
      call krungekutfs(neq,func,x0,xe,y0,init,relerr,abserr,
!                  Next original one is wrong. //// corrected by K.K
!     &          yn,esterr,nde,ier,h,yl0(neq),yln(neq),err(neq),
     &          yn,esterr,nde,ier,h, yl0, yln, err,
     &          work(1,1),work(1,2),work(1,3),work(1,4),work(1,5),
     &          work(1,6),work(1,7),work(1,8))
      return
      end
!
      subroutine krungekutfs(neq,func,x0,xe,yh0,init,relerr,abserr,
     &               yhn,esterr,nde,ier,h,yl0,yln,err,
     &               ak1,ak2,ak3,ak4,ak5,ak6,et,work)
       implicit real*8(a-h,o-z)
       parameter(c1 = 0.055455d+0, c2 = -0.035493d+0,
     &          c3 = -0.036571d+0, c4 = 0.023107d+0,
     &          c5 = -9.515d-3, c6 = 3.017d-3 )
       external func
       dimension yh0(neq),yhn(neq),yl0(neq),yln(neq),esterr(neq),
     &           ak1(neq),ak2(neq),ak3(neq),ak4(neq),ak5(neq),ak6(neq),
     &           err(neq),et(neq),work(neq,5)
       data itemax / 100000 /
       data epsmin / 1.0d-15 /
!
      if (init .eq. 1) then
! -------------- initialization (first call)-------------------------
       do 10 i = 1,neq
        yl0(i) = yh0(i)
   10  continue
!      -------- set initial step size --------

       call func(x0,yl0,ak1)
       toldy = 10.0d+74
       do 20 i = 1,neq
        if (ak1(i) .ne. 0.0d+0) then
         tol = relerr * dabs(yl0(i)) + abserr
         toldy = min(toldy, tol / dabs(ak1(i)))
        end if
   20  continue
       hmin = epsmin * (xe - x0)
       if (toldy .eq. 10.0d+74) then
        h = hmin
       else
        h = min((xe - x0) / 100.0d+0, toldy ** 0.2d+0)
       end if
!     ----------------------------------------
       init = 2
       nde = 1
      else
       nde = 0
      endif
! -------------------------------------------------------------------
      istiff = 0
      icount = 0
      ier = 0
      x = x0
!
!********************** main iteration *********************************
      do 30 iter = 1,itemax
       call krkflstep(neq,func,x,h,yl0,yln,ak1,ak2,ak3,ak4,ak5,
     &             ak6,err,work)
       nde = nde + 6
       erret = 0.0d+0
       do 40 i = 1,neq
        et(i) = relerr * 0.5d+0 * (dabs(yl0(i)) + dabs(yln(i))) + abserr
        if (et(i) .eq. 0.0d+0) then
         go to 20000
        else
         erret = max(erret, dabs(err(i)) / et(i))
        end if
   40  continue
       if (erret .ge. 1.0d+0) then
!--------------- unsuccessful step ------------------------------
        if (erret .ge. 59049.0d+0) then
         h = 0.1d+0 * h
        else
         h = 0.9d+0 * h / (erret ** 0.2d+0)
        end if
         if (h .le. epsmin) go to 30000
       else if ((x + h) .lt. xe) then
!-------------- successful step (x + h < the end point)------------
!       ********* routine for detecting stiffness *********
        icount = icount + 1
        erretl = 0.0d+0
        do 100 i = 1,neq
         erri = h * (c1 * ak1(i) + c2 * ak2(i) + c3 * ak3(i)
     &               + c4 * ak4(i) + c5 * ak5(i) + c6 * ak6(i))
         erretl = max(erretl, dabs(erri) / et(i))
  100   continue
        if (erretl .le. 1.0d+0) istiff = istiff + 1
        if (icount .eq. 50) then
         if (istiff .ge. 25) then
!       ------ stiffness was detected ------
          go to 40000
         else
          icount = 0
          istiff = 0
         end if
        end if
!       **************************************************

        call krkfhstep(neq,func,x,h,yh0,yhn,ak1,ak2,ak3,ak4,ak5,
     &              ak6,err,work)
        nde = nde + 12
        x = x + h
        x0 = x

        do 50 i = 1,neq
         yl0(i) = yln(i)
         yh0(i) = yhn(i)
   50   continue

!       ------- choose next step -----------
        if (erret .le. 1.889568d-4) then
         h = 5.0d+0 * h
        else
         h = 0.9d+0 * h / (erret ** 0.2d+0)
        end if
!       ------------------------------------
       else
!-------------- successful step (x + h > = the end point) --------------
        h = xe-x
        call krkflstep(neq,func,x,h,yl0,yln,ak1,ak2,ak3,ak4,ak5,
     &              ak6,err,work)
        call krkfhstep(neq,func,x,h,yh0,yhn,ak1,ak2,ak3,ak4,ak5,
     &              ak6,err,work)

        nde = nde + 18
!      -------- estimate global error ----------
        do 60 i = 1,neq
         esterr(i) = (yln(i) - yhn(i)) / 31.0d+0
   60   continue
!      -----------------------------------------
        x0 = xe
        do 70 i = 1,neq
         yl0(i) = yln(i)
         yh0(i) = yhn(i)
   70   continue
        return
       endif
   30 continue
!*********************************************************************
!
      ier = 10000
      write( * ,10001) x
10001 format(' ','(subr.-rkfds) trouble(too many iterations))',
     &           'at x = ',1pe15.7)
      return
20000 continue
      ier = 20000
      write( * ,20001) x
20001 format(' ','(subr.-rkfds) trouble(inappropriate error tolerance)'
     &          ,' at x = ',1pe15.7)
      return
30000 continue
      ier = 30000
      write( * ,30001) x
30001 format(' ','(subr.-rkfds) trouble(too small step size)',
     &           ' at x = ',1pe15.7)
      return
40000 continue
      ier = 40000
      write( 0 , * ) 'the equations are stiff!'
      write( 0 ,40001) x
40001 format(' ','stiffness was detected at x = ',1pe15.7,'.')
      write( 0 , * ) 'to use the subroutine (krsnat or krsat)'
     &,' is recommened.'
      return
      end
!
      subroutine krkflstep(neq,func,x,h,y0,yn,ak1,ak2,ak3,ak4,ak5,ak6,
     &                  err,work)
       implicit real*8(a-h,o-z)
       external func
       dimension y0(neq),yn(neq),ak1(neq),ak2(neq),ak3(neq),ak4(neq),
     &           ak5(neq),ak6(neq),err(neq),work(neq,5)
      lindex = 1
      call krkfehl(lindex,neq,func,x,h,y0,yn,ak1,ak2,ak3,ak4,ak5,ak6,
     &          err,work(1,1),work(1,2),work(1,3),work(1,4),work(1,5))

      return
      end
!
      subroutine krkfhstep(neq,func,x,h,y0,yn,ak1,ak2,ak3,ak4,ak5,ak6,
     &                  err,work)
       implicit real*8(a-h,o-z)
       external func
       dimension y0(neq),yn(neq),ak1(neq),ak2(neq),ak3(neq),ak4(neq),
     &           ak5(neq),ak6(neq),err(neq),work(neq,5)
      x1 = x
      h1 = 0.5d+0 * h
      lindex = 0

      call krkfehl(lindex,neq,func,x1,h1,y0,yn,ak1,ak2,ak3,ak4,ak5,ak6,
     &          err,work(1,1),work(1,2),work(1,3),work(1,4),work(1,5))
      x1 = x1 + h1
      do 10 i = 1,neq
       y0(i) = yn(i)
   10 continue
      call krkfehl(lindex,neq,func,x1,h1,y0,yn,ak1,ak2,ak3,ak4,ak5,ak6,
     &          err,work(1,1),work(1,2),work(1,3),work(1,4),work(1,5))
      return
      end
!
      subroutine krkfehl(lindex,neq,func,x,h,y0,yn,ak1,ak2,ak3,ak4,ak5,
     &                ak6,err,w2,w3,w4,w5,w6)
       implicit real*8(a-h,o-z)
       dimension y0(neq),yn(neq),ak1(neq),ak2(neq),ak3(neq),ak4(neq),
     &           ak5(neq),ak6(neq),err(neq),w2(neq),w3(neq),w4(neq),
     &           w5(neq),w6(neq)
       parameter(one = 1, two = 2, thr = 3, twl = 12)
       parameter(al2 = one / 4, al3 = thr / 8, al4 = twl / 13,
     &           al5 = one, al6 = one / 2)
       parameter(b21 = one / 4, b31 = thr / 32, b32 = 9.0d+0 / 32,
     &           b41 = 1932.0d+0 / 2197, b42 = -7200.0d+0 / 2197,
     &           b43 = 7296.0d+0 / 2197, b51 = 439.0d+0 / 216,
     &           b52 = -8, b53 = 3680.0d+0 / 513,
     &           b54 = -845.0d+0 / 4104, b61 = -8.0d+0 / 27,
     &           b62 = 2, b63 = -3544.0d+0 / 2565,
     &           b64 = 1859.0d+0 / 4104, b65 = -11.0d+0 / 40)
       parameter(ga1 = 16.0d+0 / 135, ga3 = 6656.0d+0 / 12825,
     &           ga4 = 28561.0d+0 / 56430,
     &           ga5 = -9.0d+0 / 50, ga6 = two / 55)
       parameter(da1 = one / 360, da3 = -128.0d+0 / 4275,
     &           da4 = -2197.0d+0 / 75240,
     &           da5 = one / 50, da6 = two / 55)

      call func(x,y0,ak1)
      do 10 i = 1,neq
       w2(i) = y0(i) + h * b21 * ak1(i)
   10 continue
      call func(x + al2 * h,w2,ak2)
      do 20 i = 1,neq
       w3(i) = y0(i) + h * (b31 * ak1(i) + b32 * ak2(i))
   20 continue
      call func(x + al3 * h,w3,ak3)
      do 30 i = 1,neq
       w4(i) = y0(i) + h * (b41 * ak1(i) + b42 * ak2(i) + b43 * ak3(i))
   30 continue

      call func(x + al4 * h,w4,ak4)
      do 40 i = 1,neq
       w5(i) = y0(i) + h * (b51 * ak1(i) + b52 * ak2(i)
     &                  + b53 * ak3(i) + b54 * ak4(i))
   40 continue

      call func(x + al5 * h,w5,ak5)
      do 50 i = 1,neq
       w6(i) = y0(i) + h * (b61 * ak1(i) + b62 * ak2(i) + b63 * ak3(i)
     &                  + b64 * ak4(i) + b65 * ak5(i))
   50 continue
      call func(x + al6 * h,w6,ak6)
      do 60 i = 1,neq
       yn(i) = y0(i) + h * (ga1 * ak1(i) + ga3 * ak3(i) + ga4 * ak4(i)
     &                  + ga5 * ak5(i) + ga6 * ak6(i))
   60 continue
      if  (lindex .eq. 1) then
       do 70 i = 1,neq
        err(i) = h * (da1 * ak1(i) + da3 * ak3(i) + da4 * ak4(i)
     &               + da5 * ak5(i) + da6 * ak6(i))
   70  continue
      endif
      end
!      Compute magnetic deflection (angle and displacement)
!      Use Runge-Kutta-Fehlberg Method.
!      
!     ************************************************  
      subroutine cbDefByRK2(aTrack, leng, newpos, newdir, newmom)
!     ************************************************  
!   This used Runge-Kutta-Gill  method. The step size
!   is fixed.
!
 
      implicit none

#include  "Ztrack.h"
      real(8),intent(out):: newmom(3) ! changed mom.
      real*8 geomcnst  !  Z/pc/3.3358. Z charge, pc GeV
      common /Zgeomcnst/ geomcnst

      type(track)::aTrack       !input. a charged ptcl

!                              ptcl and magfiled coord is in Exyz

      real*8  leng            !input. length  travelled in m
                     ! since we use Ruuge-Kutta-Gill method with fixed
                     ! step size, we use leng/2  as the step size assuming
                     !  leng is (Larmor radius)/5. 

      type(coord)::newpos      ! output. new position after leng move
      type(coord)::newdir      ! //      new direction cose. //

      real*8  cbetat, x0

      integer ier, m

      real*8  y0(6)
      parameter (m=3)
      real*8  wk(6,4), y(6,m)
      
      external cgeomfunc


      cbetat= leng/(m-1)
!       
!      y0(1) = aTrack.pos.xyz.x
!      y0(2) = aTrack.pos.xyz.y
!     y0(3) = aTrack.pos.xyz.z
       y0(1:3) = aTrack%pos%xyz%r(1:3)
      
!      y0(4)=  aTrack%vec%w%x
!      y0(5)=  aTrack%vec%w%y
!      y0(6)=  aTrack%vec%w%z
       y0(4:6) = aTrack%vec%w%r(1:3) 

!         pc in GeV
      geomcnst = max(sqrt(aTrack%p%fm%p(4)**2- aTrack%p%mass**2), 1.d-9)

      geomcnst = aTrack%p%charge/geomcnst/3.3358

      x0 = 0.d0
      call krungeKutG(cgeomfunc, 6,  x0, y0, m, cbetat, 6,  y, wk, ier)
      if(ier .ne. 0) then
         write(0,*) ' ier=',ier
      endif


!      newpos.x = y(1,m)
!      newpos.y = y(2,m)
!      newpos.z = y(3,m)
!       for ibm
       newpos%r(1) = y(1,m)
       newpos%r(2) = y(2,m)
       newpos%r(3) = y(3,m)
       newpos%sys='xyz'
!      newdir.r(1) = y(4,m)
!      newdir.r(2) = y(5,m)
!      newdir.r(3) = y(6,m)
!          for ibm
      newdir%r(1) = y(4,m)
      newdir%r(2) = y(5,m)
      newdir%r(3) = y(6,m)
      newmom(:) = y(4:6,m)*
     * sqrt( dot_product(aTrack%p%fm%p(1:3),aTrack%p%fm%p(1:3)) )
      newdir%sys='xyz'
      end
!     *********************************************************
!      Runge-Kutta-Gill method.
      subroutine krungeKutG(func, n, x0, y0, m, h, ny, y, wk, ierr)
      implicit none
!
!  solvoes n sets of dy/dx = f(x,y) 
!
      external func  ! input.   func is a external subroutine name
      integer  n     ! input.   number of defferential equations
      real*8 x0      !  input. initial value.
      real*8 y0(n)   !  input. initial value. 
      integer m      !  input. total number of points where the 
                     !         solutions are to be obtained (Note that 
                     !        1st point is the initial value
      real*8  h      !  input. step size
      integer ny     !  input.  see below. ny >= n.
      real*8  y(ny,m) !  output. solutions at m-points.
      real*8  wk(n,4) ! output. workin area.
      integer  ierr  !  output. condition code. 
                     !  0 ,no error was detected.
                     !  1  n<1, n> ny, m<2 
!        the structure of  func is the same as for krungekutfds.
!-----------------------------------------------------------------------
!
      real*8  zero, half, two, three, const, six
      integer k, i
      real*8 t2, t3, x

      data
     *  zero/0.d0/, 
     *  half/0.5d0/, two/2.0d0/, three/3.0d0/, six/6.0d0/
      data  const/ 0.29289321881345247560d0 /





      if(n .le. 0 .or. n .gt. ny .or. m .le. 1)  then
         ierr=1
      else
         ierr=0
         do k = 1, n
            wk(k,1)=zero
            y(k,1)=y0(k)
            wk(k,2)=y0(k)
         enddo
         do i = 1, m-1
            x=x0+(i-1)*h
            call func(x, wk(1,2), wk(1,4))
            do k=1,n
               t2=h*wk(k, 4)
               t3=half*t2-wk(k,1)
               wk(k,3)=wk(k,2)+t3
               t3=wk(k,3)-wk(k,2)
               wk(k,1)=wk(k,1)+three*t3-half*t2
            enddo

            call func(x+half*h, wk(1,3), wk(1,4))
            do  k=1,n
               t2=h*wk(k,4)
               t3=const*(t2-wk(k,1))
               wk(k,2)=wk(k,3)+t3
               t3=wk(k,2)-wk(k,3)
               wk(k,1)=wk(k,1)+three*t3-const*t2
            enddo

            call func(x+half*h, wk(1,2), wk(1,4))

            do  k = 1, n
               t2=h*wk(k,4)
               t3=(two-const)*(t2-wk(k,1))
               wk(k,3)=wk(k,2)+t3
               t3=wk(k,3)-wk(k,2)
               wk(k,1)=wk(k,1)+three*t3-(two-const)*t2
            enddo

            call func(x+h, wk(1,3), wk(1,4))
            do k=1,n
               t2=h*wk(k,4)
               t3=(t2-two*wk(k,1))/six
               wk(k,2)=wk(k,3)+t3
               y(k,i+1)=wk(k,2)
               t3=wk(k,2)-wk(k,3)
               wk(k,1)=wk(k,1)+three*t3-half*t2
            enddo
         enddo
      endif
      end
!     ************************************************  
      subroutine cbDefUser(aTrack, leng, newpos, newdir,newmom )
!     ************************************************  
!    
!      another method to get the new position  and direction
!      cosines when magnetic field exists.
!  This is called when UseRungeKutta is 8.
 
      implicit none

#include  "Ztrack.h"
#include  "Zcode.h"

      real*8 geomcnst  !  Z/pc/3.3358. Z charge, pc GeV
      common /Zgeomcnst/ geomcnst

      type(track)::aTrack   ! input. a charged ptcl
      real(8),intent(out):: newmom ! changed momentum
!     
!        current position
!            aTrack.pos.xyz.x, aTrack.pos.xyz.y, aTrack.pos.xyz.z (m)
!        current direction cosines
!            aTrack.vec.w.x,  aTrack.vec.w.y,  aTrack.vec.w.z
!        current momentum,  total energy (GeV)
!            aTrack.p.fm.p(i) i=1,4
!        mass                  patcl code      subcode           charge
!            aTrack.p.mass   aTrack.p.code  aTrack.p.subcode,aTrack.p.charge
!   if code = kgnuc (isotope), subcode has mass number.
!

             

!               ptcl and magfiled coord is in Exyz

      real*8  leng   ! input. length  travelled in m

      type(coord)::newpos  ! output. new position after moving leng (m)
                            !
                            !   newpos.x newpos.y, newpos.z must be given
                            !   newpos.sys='xyz' must be gieven
      type(coord)::newdir  ! //      new direction cose. //
                            !   newdir.x, newdir.y, newdir.z must be given
                            !   newdir.sys='xyz' must be given

!         pc in GeV
      geomcnst = max(sqrt(aTrack%p%fm%p(4)**2- aTrack%p%mass**2), 1.d-9)
      geomcnst = aTrack%p%charge/geomcnst/3.3358
!c
!c      You may see cgeomfunc how to get magnetic field etc
!c
!       Note:  cgeomag routine below gives B in north, east, down
!              coordinate system, so that it must be converted into
!              E-xyz system by calling ctransMagTo like below.
!
!      call cgeomag(YearOfGeomag, llh, B, icon)  ! B is ned
!      call ctransMagTo('xyz', llh, B, B)        ! to xyx
!

      newpos%sys = 'xyz'
      newdir%sys = 'xyz'
      end

