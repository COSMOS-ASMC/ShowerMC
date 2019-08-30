!     integrate pair function from xmin to xmax
!     generic one.
#include "ZepicsBD.h"
      implicit none
!
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"

      character*120  file
      character(len=8):: name
      integer io, i, norm
      integer::how=0
      real*8  xmax, xmin, tprob,  Nc
      real*8  E1, E2, Eg, Egme, step
       type(epmedia)::  media
      character(8):: force
      io = 10
      read(*,*) 
     *  name, force,  E1, E2, step,  norm,  LPMeffect


      write(0,'(a)') 
     * ' media, force,  E1, E2, step,norm, '//
     * ' LPMeffect'

      Eg = E1
      if(Eg  <=  2*masele) then
         Eg = 2*masele*1.001
      endif
      file ="$EPICSTOP/Data/BaseM/"//trim(name)
!           Air*0.01 type is processed in  next call.           
      call epBPgeneini(file, media, how)

      if(norm .eq. 5) then
         Nc = 1./tprob
      elseif(norm .eq. 1) then 
         Nc=media%mbtoPX0       ! prob/X0
      elseif( norm .eq.  2) then
         Nc = 1
      elseif(norm .eq.  3) then
         Nc = media%mbtoPgrm
      elseif(norm .eq. 4) then
         Nc = media%mbtoPcm
      else
         call cerrorMsg('input error for norm',0)
      endif
      write(*,
     *     '("# Eg(GeV)   XS   XS*Nc  xmin xmax ")')
      write(*,
     *     '("# XS is in mb;  XS*Nc is desired XS")')
      
      do while( Eg < E2*1.000001d0 )
         Egme = Eg/masele
         xmax = 1. - masele/Eg
         xmin = masele/Eg
         call epPrgenePreInte(media,  Egme)
         call epPrgeneTX(xmin, xmax, tprob)

         write(*,'(1p, 5g14.6)') 
     *      Eg, tprob, tprob*Nc,  xmin, xmax
         Eg = Eg*10.0d0**step
      enddo
      end


