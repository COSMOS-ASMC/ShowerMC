!     ========================================
!       make table for length <--> thickness convertion
!      so that 1 dim. calc. be  faster.
!  Actually this one is used only for making table.  see cl2tTA.f 
!
      subroutine cl2tTbl(h1, h2, cosz1, cosz2, step,
     *  lengtb, htb, costb,  thicktb,  maxsize, tblsize)
!
      use modAtmosDef
      implicit none
!#include "Zearth.h"
!        this program makes a table of slant air thickenss from h1 to h2
!        with a given zenith angle.  The table is made at  0 (= at h1),
!        step, 2*step, ... along  the path
      real*8 h1  ! input.  higher vertical height in m
      real*8 h2   ! input.  lower vertical height in m
!                  h1=30d3, h2=-0.1d3 would be standard.
      real*8 cosz1  ! input.  cos of zenith angle at h1
      real*8 cosz2  ! input.  cos of zeniht angle at h2
      real*8 step  ! input.  step of the length from h1 to h2
      integer maxsize  ! input. max usable size of lengtb, and thicktb
      real*8  lengtb(maxsize) ! outpu. lengtb(i) = (i-1)*step
      real*8  thicktb(maxsize)  ! output.  thicktb(i) =  air thickness 
      real*8  htb(maxsize)
      real*8  costb(maxsize)
!                                along the path at lengtb(i).
      integer tblsize  ! output .    actual size of the tables made.


      real*8 clen2thickEx
      real*8  cost, sumt, z
      real*8  temp, cnewcos, cnewh

      integer i

!       step 
!    //////////
!      write(*,*) ' h1=',h1,' h2=',h2, ' cosz=',cosz,
!     *    ' step=',step, ' maxsize=',maxsize
! //////////

      cost= cosz1
      sumt = 0.
      lengtb(1) = 0.
      thicktb(1) = 0.
      htb(1) = h1
      z = h1

      do i = 2, maxsize
         lengtb(i) = lengtb(i-1) + step
         thicktb(i) = thicktb(i-1) + clen2thickEx(z, cost, step, 8)
! //////////
!          write(*, *) lengtb(i), thicktb(i), z, cost
!///////////////
         temp= max(cosz2, cnewcos(z+Eradius, cost, step))
         z = cnewh(z+Eradius, cost, step) - Eradius
         htb(i) = z
         cost = temp
         costb(i) = cost
         if(z .le. h2) then
            tblsize = i
            goto 100
         endif
      enddo
      tblsize = maxsize
 100  continue
! ///////////
!       write(*,*) ' max thick=', thicktb(tblsize),lengtb(tblsize), 
!     *    ' size=', tblsize
! ///////////////
      end
!     *********************************************
      real*8 function clen2thickT(z, cosz, leng)
!
!     z: real*8. input.  vertical height in m.
!  cosz: real*8. input.  cos of zenith angle at z.
!   leng: real*8. input.  length along cosz direction.
!     function value.      thickness of air in kg/m2. for leng
      use modAtmosDef
      implicit none
!#include "Zearth.h"
!#include "Zatmos.h"
      real*8 z, cosz, leng
      real*8  z2, cos2, sin2
      real*8  r1, r2, st1, st2, dt
      integer loca, icon
      
!          get height at leng ahead: z2
      r1 = z + Eradius
      call cnewcossin(r1, cosz, leng, r2, cos2, sin2)
      z2 = r2- Eradius
!         get slant thickness, st1, at z
      if(z .eq. Zsave) then
         st1 = ThickSave
      else
         call cz2t(z, cosz, loca, st1, icon)
         if(icon .ne. 0) then
            write(*, *) ' error: in clen2thickT; z=', z, 
     *      ' cosz=',cosz, ' leng=', leng
            stop 9999
         endif
      endif
      
      if(leng .lt. LenStep) then
!          this approx is better but somewhat 
!          inconsistent with table method. why ??
         call csmlL2t(z, cosz, leng, dt)
         st2 = st1 + dt
         write(*, *) ' small approx dt=', dt, ' for leng=',leng
!
      else
!          get slant thickness, st2, for  z2
         call cz2t(z2, cos2, loca, st2, icon)
         if(icon .ne. 0) then 
            write(*, *)' error in clen2thickT; z=', z,
     *     ' cosz=', cosz, ' leng=', leng, ' z2=', z2
            stop 9999
         endif
         dt = st2 - st1
!
!         write(*, *) ' no small appx; dt=',dt
!
      endif
      clen2thickT = dt
      Zsave = z2
      ThickSave =  st2
      end
!      ***************************
!         get slan thickness for length s along cosz at height h.
!         s should be small (s < LenStep~ 100 m)
      subroutine csmlL2t(h, cosz, s, t)
      use modAtmosDef
      implicit none
! #include "Zearth.h"
! #include "Zatmos.h"
      real*8 h  ! input.  hight in m
      real*8 cosz  ! input.  cos of zenith angle at h.
      real*8  s     ! input. small length in m along cosz from h.
      real*8  t    !  output.  thickness of air for s. <0 if s < 0
!
      real*8  cs, sn2,  ss, sold, eps
!      real*8 f1, f2, f3, rho, rho1, rho2, rho3, r
!      real*8 f1, f2, f3, rho, rho1, rho2, r
      real*8 f1, f2,  rho, rho1, rho2, r
      integer i
!      real*8 cvh2den, cvh2denp, cvh2den2p, cvh2den3p
      real*8 cvh2den, cvh2denp, cvh2den2p
      data eps/1.d-6/

!      f1(ss) = ss*(-cs/2.d0 + ss*sn2/r *(1.d0/6.d0 +
!     *     ss* cs/8.d0/r))
      f1(ss) = ss*(-cs/2.d0 + ss*sn2/r/6.d0)

      f2(ss) = ss**2* cs*(cs/3.d0 - ss*sn2/4.d0/r)
!      f3(ss) = -(ss*cs)**3/4.d0
!
      r = Eradius + h
      cs = cosz
      sn2 = 1.d0- cs**2
      rho = cvh2den(h)
      rho1 = cvh2denp(h)
      rho2 = cvh2den2p(h)
!      rho3 = cvh2den3p(h)
!      t = s*(rho + rho1*f1(s) + rho2*f2(s) + rho3*f3(s))
      t = s*(rho + rho1*f1(s) + rho2*f2(s) )
      return
!  ---------------------------
!         get length s for thickenss t at height h along cosz
!         |t| should be small so that s be < ~ 100  m.
      entry  csmlT2l(h, cosz, t, s)
!
      r = Eradius + h
      cs = cosz
      sn2 = 1.d0- cs**2

      rho = cvh2den(h)
      rho1 = cvh2denp(h)
      rho2 = cvh2den2p(h)
!      rho3 = cvh2den3p(h)
      sold = t/rho
      do i = 1, 6
         s =t/ ( rho + rho1*f1(sold) + rho2 *f2(sold)
     * )
!     *    +  rho3*f3(sold) )
         if( abs(s) .lt. 1.d0 .and. abs(s-sold) .lt. eps) then
            goto 10
         elseif(abs( (s-sold)/s ) .lt. eps) then
            goto 10
         endif
         sold = s
      enddo
 10   continue
      end
!     **************
!         get slant thickenss at z from Htop.
      subroutine cz2t(z, cosz, loca, st, icon)
      use modAtmosDef
      implicit none

! #include "Zearth.h"
! #include "Zatmos.h"
      real*8 z  ! input. vertical height in m
      real*8 cosz  ! input. cos of zenith angle at z.
      integer loca  ! output. loca is  such that
!                     HeightTbl(loca+1) <= z < HeightTbl(loca)  
      real*8  st   ! output. slant thikness in kg/m^2 from Htop along
!                    cosz 
      integer icon ! output.  = 0 if normal. =1 if loca is out of  range.
      real*8 clenbetween2h, ds, r1, r2, dt


!        get location of nearest integral point of upper height
      call kdwhereis(z, NumStep, HeightTbl, 1, loca)
      if(loca .ge. NumStep) then
         st = ThickTbl(NumStep)
         icon = 0
      elseif(loca .le. 0) then
         if(z .eq. Htop) then
            loca =1 
            st = 0.
            icon = 0
         else
            write(*,*) ' error in cz2t; z=',z
            stop 999
         endif
      else
!           get length from HeightTbl(loca) to z
         r1 = z + Eradius               
         r2 = HeightTbl(loca) + Eradius !  r1 < r2
         ds = clenbetween2h(r2, r1, CosTbl(loca) ) ! ds > 0
         if(ds .lt. LenStep/2) then
!                 thickness for ds
            call csmlL2t(HeightTbl(loca), CosTbl(loca), ds, dt)  
            st =  ThickTbl(loca) + dt
         else
            ds  =  ds- LenStep   ! < 0
            call csmlL2t(HeightTbl(loca+1), CosTbl(loca+1), ds, dt)  ! dt < 0
            st = ThickTbl(loca+1) + dt
         endif
         icon = 0
      endif
      end
!              
!     ***************
!         get length for slant thickness t along cosz form height z
      real*8 function ct2lT(z, cosz, t)
      use modAtmosDef
      implicit none
! #include "Zatmos.h"
! #include "Zearth.h"

      real*8 z  ! input vertical height in m
      real*8 cosz  ! input.  cos of zenith angle at z
      real*8 t  !  input.  thickness of air in kg/m^2 along cosz
!
      integer loca, icon
      real*8 st1, st2, r0, r1, r2, rx,  cosx, cos2, sin2
      real*8 s1, s2, dt
      real*8 clenbetween2h, cvh2den
      
      r1 = z + Eradius

!      get slant thickness st1 from Htop to z
      if(z .eq. Zsave) then
         st1 = ThickSave
      else
         call cz2t(z, cosz, loca, st1, icon)
         if(icon .ne. 0) then
            write(*,*) ' error in cz2t; z=', z, ' cosz=', cosz
            stop 999
         endif
      endif
!        new point slant thickness
      st2 = st1 + t
!
!       integral location of z:  H(loca+1)<= z < H(loca)
      call kdwhereis(z, NumStep, HeightTbl, 1, loca)
      if(loca .le. 0 .or. loca .ge. NumStep) then
         write(0,*) ' error in ct2lT', ' z=',z, 'cos=',cosz,
     *              ' t = ',t 
         stop 999
      else
         r0 = z + Eradius
         rx = HeightTbl(loca) + Eradius
         s1 = - clenbetween2h(r0, rx, cosz) !  s1 > 0
      endif
!         get integral point for st2; T(loca) <= st2 < T(loca+1)
      call kdwhereis(st2, NumStep, ThickTbl, 1, loca)

      if(loca .ge. NumStep) then
         ct2lT = 1.d5
         Zsave = HeightTbl(NumStep)
         ThickSave = ThickTbl(NumStep)
      elseif(loca .le. 0.) then
         write(*,*) ' error in ct2lT; z=',z, ' t=',t, ' st2=',st2
         stop 9999 
      else
	 if(t/cvh2den(z) .lt. LenStep .and. cosz .lt. 0.22d0) then
	    call csmlT2l(z, cosz, t, ct2lT)
!              we don't get ThickSave, Zsave in this case. leave the
!              buisiness to the next call
         else
   	    dt = st2 - ThickTbl(loca)            	
	    if(dt .lt. ThickTbl(loca+1) - st2) then
!                get length for dt
               call csmlT2l(HeightTbl(loca), CosTbl(loca), dt, s2)
!                get height at s2
               rx = HeightTbl(loca) + Eradius
               cosx = CosTbl(loca)
               call cnewcossin(rx, cosx, s2, r2, cos2, sin2)
            else
               dt = ThickTbl(loca+1) - st2
                 call csmlT2l(HeightTbl(loca+1), CosTbl(loca+1),
     *         -dt, s2)  !  s2 < 0
               rx = HeightTbl(loca+1) + Eradius
               cosx = CosTbl(loca+1)
               call cnewcossin(rx, cosx, s2, r2, cos2, sin2)
            endif
            ct2lT = clenbetween2h(r1, r2, cosz)
            Zsave = r2 - Eradius
            ThickSave = st2
         endif	
      endif
      end

