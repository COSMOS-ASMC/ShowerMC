!     integrate brems function from xmin to xmax
!      generic one.
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
      integer nerg
      integer io, i, norm, how
      real*8  Ee, x, f, epBrgenex, xmax
      real*8  xmin, tprob, xx, fxx, Nc
      real*8  E1, E2, Ek, Eeme, step
       type(epmedia)::  media
      character(8):: force

      io = 10
      read(*,*) 
     *  name, force, how, E1, E2, step,  norm,  LPMeffect, Flpm


      write(0,'(a)') 
     * ' media, force, how, E1, E2, step,norm, '//
     * ' LPMeffect, Flpm '

!      EpartialSC = 1.  ! special test

      Ek = E1
      file ="$EPICSTOP/Data/BaseM/"//trim(name)
!           Air*0.01 type is processed in  next call.           
      call epBPgeneini(file, media, how)


      do while( Ek < E2*1.000001d0 )
         Ee = Ek + masele
         Eeme = Ee/masele
         call  epBrgenePreInte(media, force, Eeme)
!           neglect input xmin
         if(Ee .lt. media%cnst%BrEemaxS) then
            xmin = media%cnst%BrEgminS/Ee
         elseif(Ee .lt. media%cnst%BrEemaxS2) then
            xmin = media%cnst%BrEgminS2
         else
            xmin = media%cnst%BremEgmin
         endif
         write(0,
     *     '(a, a,  i3, 1p, 3g12.4,0p, i3,1p,g12.4,L2,g12.4)') 
     *   name, force, how, E1, E2, step, norm, xmin, LPMeffect, Flpm

         xmax = 1. - masele/Ee
         call epBrgeneTX(xmin, xmax, tprob)
         if(norm .eq. 5) then
            Nc = 1./tprob
         elseif(norm .eq. 1) then 
            Nc=media%mbtoPX0    ! prob/X0
         elseif( norm .eq.  2) then
            Nc = 1
         elseif(norm .eq.  3) then
            Nc = media%mbtoPgrm
         elseif(norm .eq. 4) then
            Nc = media%mbtoPcm
         else
            call cerrorMsg('input error for norm',0)
         endif
         write(*,'(1p, 6g14.6)') 
     *      Ee, tprob, tprob*Nc, Ek, xmin, xmax
         Ek = Ek*10.0d0**step
      enddo
      end
