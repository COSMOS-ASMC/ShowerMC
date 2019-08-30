      implicit none
c
c         total prob. of Coherent scattering.
c
#include "Zglobalc.h"
#include "ZbasicCnst.h"      
#include "ZepManager.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"
c #include "ZBPgene.h"

      character*130  file
      character(len=8)::name
      integer io, result

      real*8  Eg, E1, E2, step, xprob(5), txray(5)
      real*8  tprob, p, norm
      record /epmedia/ media
      integer XCOM, icon, isel

      io = 10

      E1=1.e-6
      E2=100.
      step = 0.05
      XCOM=1
      write(0,*) 'Enter  media name'
      read( *, * )   name
      write(0,*)' Enter next for XS in unit of '
      write(0,*)
     *' 1->1/rl, 2->1/(g/cm2), 3->1/cm, 4->mb '
      read(*,*) isel
      if(isel == 4) then
         write(0,*)
     *        ' For mixed or compound media, mb value is not'
         write(0,*)
     *        ' not correct but other values are correct'          
      endif
      
      file = "$EPICSTOP/Data/Media/"//trim(name)
      write(0,*) 'file=',trim(file)
      write(0,*) "# media=", name
      write(0,*)
     *     "# Eg(MeV)  coh  incoh(comp)  p.e p.c(n) pc(e) pc(T) "
      write(0,'(a, i2, a)')
     *  "unit is ",isel, "=> 1->1/rl, 2->1/(g/cm2), 3->1/cm, 4->mb"

      if(isel == 1) then
         norm = 1.
      elseif( isel == 2 ) then
         norm = media.X0g
      elseif( isel == 3) then
         norm = media.X0
      elseif( isel == 4 ) then
         norm =media.mbtoPX0
      else
         write(0,*) ' error unit spec =',isel
         stop
      endif

      
      call copenf(io, file, result)
      if( result == 0 ) then
         call epReadMTbl(io, media)
         call epGetEffZA(media)
         MediaDir(1) ="$EPICSTOP/Data/Media"
         Excom1 = 100.
         call epReadXXsec(media, icon)
         if(icon .ne. 0) then
            write(0,*) ' Xcom read err'
         endif

         Eg = E1

         do 
            if( media.xcom.size .gt. 0 .and.
     *              Eg .le. Excom1 ) then
               call epXrayp(media, Eg, 1, 5,  xprob, txray)
            else
               write(0,*) "media.xcom.size=", media.xcom.size
               write(0,*) " err in XS.out"
               stop
            endif
            write(*, '(1p, 7g12.4)') 
     *       Eg*1000.,  
     *       xprob(1)/norm, xprob(2)/norm,
     *       xprob(3)/norm, xprob(4)/norm, xprob(5)/norm,
     *           (xprob(4)+xprob(5))/norm
            if(Eg == E2) exit
            Eg = min(Eg * 10.0d0**step, E2)
         enddo
      else
         write(0,*) ' media file error: file= ',trim(file)
      endif
      end
