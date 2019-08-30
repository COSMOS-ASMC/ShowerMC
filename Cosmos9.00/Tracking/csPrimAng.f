#if defined NEXT486
#define IMAG_P dimag
#elif defined PCLinux
#define IMAG_P dimag
#else
#define IMAG_P imag
#endif
      
!    Test program is in Test/testSPrimAng.f
!----------------------samples primary direction cos. at observation place.
!                      system is det
!          
      subroutine csPrimAng(dir)
      implicit none
!          
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zincidentp.h"
#include  "Zobs.h"
#include  "Zobsp.h"

      type(coord)::dir
      character*50 msg
      real*8 u,  cs, sn, sint, fai, za

      if(Za1ry .eq. 'is') then
         call csPrimIsoAng(dir)
      elseif(Za1ry .eq. 'ps') then
!                  from point source
         call csPrimPSAng(dir)
      elseif(Za1ry .eq. 'aps') then
!                around point source
         call csPrimAPSAng(dir)
      elseif(Za1ry .eq. 'cos' ) then
!          cos dcos
         call ksampLin
     *     (1.0d0, 0.d0,  real(CosZenith), IMAG(CosZenith), za)
         dir%r(3) =  - za   ! going down is negative
         if(ObsPlane .eq.  spherical) then
!              actually setting here is dummy.
            call kcossn(cs, sn)
         else
            call rndc(u)
            fai = (IMAG_P(Azimuth)- real(Azimuth)) *u +
     *           real(Azimuth)
            cs = cos(Torad*fai)
            sn = sin(Torad*fai)
         endif
         sint = sqrt(1.d0-dir%r(3)**2)
         dir%r(1) = - sint*cs   ! - is needed for going down ptcl.
         dir%r(2) = - sint*sn
         dir%sys = 'det'     !   fof ObsPlane=spherical, this is  reset later
      else
         write(msg,*) 'strange Za1ry=',Za1ry
         call cerrorMsg(msg, 0)
      endif
      end
#if defined NEXT486
#define IMAG_P dimag
#elif defined PCLinux
#define IMAG_P dimag
#else
#define IMAG_P imag
#endif
!      -----------------------------------------

      subroutine csPrimIsoAng(dir)
      implicit none
!
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zincidentp.h"
      type(coord)::dir
!
      real*8 u, cs, sn, sint, fai
          call rndc(u)
          dir%r(3) =  -(  (IMAG_P(CosZenith)- real(CosZenith) )*u +
     *             real(CosZenith) )    ! going down is negative 
          call rndc(u)
          fai = (IMAG_P(Azimuth)- real(Azimuth)) *u + 
     *           real(Azimuth)
          cs = cos(Torad*fai)
          sn = sin(Torad*fai)
!ccc          call kcossn(cs, sn)
          sint = sqrt(1.d0-dir%r(3)**2)
          dir%r(1) = - sint*cs  ! - is needed for going down ptcl.
          dir%r(2) = - sint*sn
          dir%sys = 'det'
      end
! ----------------------------- 
!           sample 1ry angle from point source
      subroutine csPrimPSAng(dir)
      implicit none

#include  "Ztrack.h"
#include  "Zincidentp.h"
#include  "Zincidentv.h"
      type(coord)::dir
!
      real*8 u, h, w1, w2, w3
          call rndc(u)
          h=(IMAG_P(Obsvhour)-real(Obsvhour))*u + real(Obsvhour)
          call rndc(u)
          if(u .lt. .5) then
                h=-h
          endif
!             for source SourceDec at hour h, get horizontal vector
          call kdhtoh(SourceDec, h, w1, w2, w3)
          w3 = - w3     ! note  vector def. is opposit
          w2 = - w2
          w1 = - w1
!             convert to detector system
          call khtad(w1, w2, w3, dir%r(1), dir%r(2), dir%r(3))
          dir%sys = 'det'
       end
!    ---------------------------------------
!     arround point sorucec  
       subroutine csPrimAPSAng(dir)

      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Zincidentp.h"
#include  "Zincidentv.h"
      type(coord)::dir

      real*8  u, w3p, decx, h, w1, w2, w3
!            sample angle in horizontal system
          call rndc(u)
          w3p=(Cspsmx -Cspsmn)*u + Cspsmn
          decx=acos(w3p)*Todeg
!c          if(decx .gt. 90.) then ! bug found Oct.15. by Kawata.
             decx=90.-decx
!c          endif
          call rndc(u)
          h=( IMAG_P(Obsvhour)-real(Obsvhour))*u + real(Obsvhour)
          call rndc(u)
          if(u .lt. .5) then
                h=-h
          endif
          call kdhtoh(decx, h, w1, w2, w3)
          w3 = - w3   ! note vector direction is opposit.
          w2 = - w2
          w1 = - w1
          call khtad(w1, w2, w3, dir%r(1), dir%r(2), dir%r(3))
          dir%sys = 'det'
      end
!       ------------------------- init for sampling 1ry angle
      subroutine ciniSPrimAng

      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Zincidentp.h"
#include  "Zincidentv.h"
      
      real*8  smx, smn, w3min, hmaxd, hmind, hmax, hmin, 
     *        wh1, wh2, wh3, w3max
      integer icon, jcon
      external cblkIncident
      character*200  msg
!
!               used if arround point source
           smx=pi/2 -(SourceDec+Ddelta)*Torad
           smn=pi/2 -(SourceDec-Ddelta)*Torad
           Cspsmx=cos(smx)
           Cspsmn=cos(smn)
           if(Za1ry .eq. 'ps' .or. Za1ry .eq. 'aps') then
               write(msg,*) '1ry from point source is specified'
               call cerrorMsg(msg, 1)
               w3min=real(CosZenith)
!                  from declination and zenith angle, get time
!                  from meridian
               call kdzth2(SourceDec, w3min, hmaxd, icon)
               hmax=hmaxd
               if(icon .le. 1) then
                   w3max = IMAG_P(CosZenith)
                   call kdzth2(SourceDec, w3max, hmind, jcon)
                   hmin=hmind
                   if(jcon .eq. 2) then
                       hmin=0.
                       call kdhtoh(SourceDec, 0.d0, wh1, wh2, wh3)
                       w3max=wh3
                   endif
                   Obsvhour=cmplx(hmin, hmax, 8)
                   if(hmin .gt. hmax) then
                       write(msg,*) ' invalid CosZenith=',CosZenith
                       call cerrorMsg(msg, 0)
                   endif
                   if(icon .eq. 1) then
                       write(msg,*) ' max observable zenith is',
     *                 ' smaller than the input'
                       call cerrorMsg(msg, 1)
                   endif
                   write(msg,*)' range of  hour =',hmin,' to ',
     *             hmax,'  which corresponds to cos(zenith)=',
     *             w3min, w3max
                   call cerrorMsg(msg, 1)
               elseif(icon .eq. 2) then
!                       no observation range
                   write(msg,*) ' range of cos(zenith)=', CosZenith,
     *             ' but no such angle is realizable'
                   call cerrorMsg(msg, 0)
               endif
           endif
       end

