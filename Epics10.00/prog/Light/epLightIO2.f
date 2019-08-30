!      stuff  indirectly of directly connected with epgen.f
!  
      subroutine epLightIOwrite1stCol
        ! Just after the first col. takes place 
      implicit none
#include "ZsepManager.h"
#include "ZepTrackv.h"
      integer cevent

      integer::code, subcode, chg
      real(8):: erg, xin, yin,  zin
      real(8):: wx, wy, wz
      real(8):: wl, mass
      integer:: compno

      call epqevn(cevent)   ! current event #
      call epwrite1stCol(cevent)
!       copy the charged path  in scratch disk for later Cerenkov
!       generation
      rewind IoScratch
!////////////
!      write(0,*) ' adding ceren data from scratch to true'
!////////////
      do while( .true. )
         if(Out1ry > 0 ) then
!            read(IoScratch, *, end=100 ) code, subcode, chg, erg,
!     *         xin, yin, zin, wx, wy, wz, wl, mass, compno
!            write(OutPrimaryFileNo,
!     *      '(i6,2i4,1p,4g15.7,3g18.9,2g12.4,0p,i6)' )
!     *          code, subcode, chg, erg,
!     *         xin, yin, zin, wx, wy, wz, wl, mass, compno
         else
            read(IoScratch, end=100) code, subcode, chg, erg,
     *         xin, yin, zin, wx, wy, wz, wl, mass, compno
            write(OutPrimaryFileNo) code, subcode, chg, erg,
     *         xin, yin, zin, wx, wy, wz, wl, mass, compno
         endif
      enddo
 100  continue
      OutPrimEff = OutPrimaryFileNo
      end
      subroutine epLightIOwriteCell
      implicit none
#include "ZsepManager.h"
#include "ZepTrackp.h"
      integer,save::code=-1000
      integer,save::subcode=0, chg=0
      real(8),save:: erg=0.
      real(8),save:: xin=0., yin=0., zin=0.
      real(8),save:: wx=0., wy=0., wz=0.
      real(8),save:: wl=0., mass=0.
      integer,save:: dummyno = 0

      if( Out1ry > 0 ) then
         write(0,*) ' should not happen if Light=21=',Light
         stop
!         call epLightPushCells(OutPrimaryFileNo, Out1ry)
!         write(OutPrimaryFileNo, *)  ' ' ! mark end
      else
         call epLightPushCells(OutPrimaryFileNo, Out1ry)
                       ! mark end
         write(OutPrimaryFileNo)  code,
     *        subcode, chg, erg,
     *        xin, yin, zin, wx, wy, wz,
     *        wl, mass, dummyno
      endif
      end

