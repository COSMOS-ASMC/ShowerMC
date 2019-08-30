#include "ZepicsBD.h"
      implicit none
!
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"

      character*50  file
      character*24  basemedia
      integer nerg
      integer io, i, norm, how, icon
      real*8  Ee, x, f, epBrgenex, xmax
      real*8  xmin, tprob, xx, fxx, Nc
      real*8  E1, Eeme
       type(epmedia)::  media
      character(8):: force

      io = 10
      read(*,*) 
     *  basemedia, force, how, E1, norm, xmin, LPMeffect, Flpm


      write(0,*) ' input for '
      write(0,'(a)') 
     * ' basemedia, force, how, E1, norm, xmin, LPMeffect, Flpm:'
      write(0, *) 
     *  basemedia, force, how, E1, norm, xmin, LPMeffect, Flpm

!      EpartialSC = 1.  ! special test


      Ee = E1 + masele
      Eeme = Ee/masele
      file ="$EPICSTOP/Data/BaseM/"//trim(basemedia)
!           Air*0.01 type is processed in  next call.           
      call epBPgeneini(file, media, how)
      call  epBrgenePreInte(media, force, Eeme)
!!!!!!!!!!   2018 Jan.8     
!!!   in above, epSetSTblCns has set various min but
!!!   it could be different from the vaule in Data/Media/W ... etc
!!!   so we read media data to reset such min.
      file ="$EPICSTOP/Data/Media/"//trim(basemedia)
      call copenf(io, file, icon)
      call epReadTab(io, media) ! this calls epSetSTblCns again
                      ! but no problem.
      close(io)
      write(0,*) ' media file read and min has been reset'
!!!!!!!!!!!!      
!           neglect input xmin
      if(Ee .lt. media%cnst%BrEemaxS) then
         xmin = media%cnst%BrEgminS/Ee
      elseif(Ee .lt. media%cnst%BrEemaxS2) then
         xmin = media%cnst%BrEgminS2
      else
         xmin = media%cnst%BremEgmin
      endif

      xmax = 1. - masele/Ee
      write(0,*) ' xmin =', xmin
      write(0,*) ' xmax =', xmax
!/////////
      write(0,*) ' integration from=', xmin, ' to', xmax 
!/////////////
      call epBrgeneTX(xmin, xmax, tprob)
!/////////
      write(0,*) ' tprob=', tprob,' in mb'
!/////////////
      if(norm .eq. 5) then
         Nc = 1./tprob
      elseif(norm .eq. 1) then 
         Nc=media%mbtoPX0  ! prob/X0
      elseif( norm .eq.  2) then
         Nc = 1
      elseif(norm .eq.  3) then
         Nc = media%mbtoPgrm
      elseif(norm .eq. 4) then
         Nc = media%mbtoPcm
      else
         call cerrorMsg('input error for norm',0)
      endif

      x = xmin
             
      do 
         f = epBrgenex(x) *NC
         if(f .gt. 0.) then
!           f = ds/dx so fxx= ds/dxx is (Ek/Ee)* ds/dx
            xx = x*Ee/E1  ! xx = Eg/Ek
            fxx = f*(E1/Ee)
            write(*, '(1p, 5g12.4)')
     *           xx, xx*fxx, Ee, Nc, tprob
         else
            write(0,*) 'f =', f, '<0 at ',x, Ee,
     *      ' in DrawBrams'
             stop 1235 
         endif
         if( x == xmax ) exit
         if(x .gt.  xmax-0.5d0 ) then
            x = x + 1.d-3
            if(x > xmax) then
               x = xmax
            endif
         else
            x = x * 10.**0.01
         endif
      enddo
      write(*,*)
      end
