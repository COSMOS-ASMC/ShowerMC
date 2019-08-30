!c       small modification of testGetXXsec.f in ../XCOM
!       to get only pair total x-section
      implicit none
#include "Zmedia.h"
#include "ZepManager.h"
       type(epmedia)::   media
      integer  icon
!
      character(len=8):: name
      character(len=120):: file
      integer:: func, norm, how
      logical::  dummy
      real(8):: E1, E2, step, Nc
      real*4  xsec(7), Ex
      MediaDir(1)="$EPICSTOP/Data/Media/"
!         $media $func  $E1 $E2 $step $norm $LPMeffect
      read(*, *)
     *    media%name, func,  E1, E2, step, norm, dummy
      write(0,*) "media name=", media%name

      file ="$EPICSTOP/Data/BaseM/"//trim(media%name)
      how =  0  ! dummy for pair

      call epBPgeneini(file,  media, how)
!       call epPrgenePreInte(media,  Egme) not needed

      call epReadXXsec(media, icon)
      if(icon /= 0 ) then
         write(0,*) ' icon=',icon, ' after epReadXXsec in',
     *  ' IntePair%f'
         stop
      endif
      if(norm .eq. 5) then
         write(0,*) ' nonsence in this case'
         stop
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

      Ex = E1
      do 
         call cGetXXsec(Ex, media%xcom%tab,
     *        media%xcom%size, 1, 7, xsec, icon)
         if( icon /= 0 ) then
            write(0,*) ' in IntePair3 icon =',icon
            stop
         endif
         xsec(4:5) = xsec(4:5)/media%mbtoPgrm ! to mb
         write(*,'(1p,3g14.4)') Ex, sum(xsec(4:5))*Nc  ! mb  to ...
         Ex = Ex*10.**step
         if( Ex >= E2*1.000001 ) exit
      enddo
      end
