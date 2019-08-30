#include "ZepicsBD.h"
      implicit none
!
!      
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
#include  "Zepdebug.h"      
      real*8 length
      integer icon, i

       type(epPos)::  org
       type(ep3Vec)::  abc
      integer  cmpn

      character*100 dsn1        ! input config data file path


      MsgLevel = 1
      cmpn = 1
      call cerrorMsg(
     * 'Enter config file path(say, ../UserHook/NewVol/config)', 1)
      read(*,'(a)')  dsn1
      read(*, *) cTrack%pos%x, cTrack%pos%y, cTrack%pos%z
      read(*, *) cTrack%w%x, cTrack%w%y, cTrack%w%z


       MediaDir(1) = '$EPICSTOP/Data/Media'
       call eprcnf(dsn1)

!
!       call epOutCnf(6)    !  to see read config.
       call epparaphtbl(0)
       call epqcnf(org, abc)
!       
       write(0, *)  'org and abc in world'
       write(0, *)  org%x, org%y, org%z
       write(0, *)  abc%x, abc%y, abc%z
       Cn = cmpn
       i = 0
       !!!!!!!! this is for octagon !!!!!!!!!!
       call epsoctagon(Det%cmp(Cn), cTrack%pos, icon)
       write(0,*) ' inout icon=', icon

       call epbndry2(cmpn, length, icon)
       write(0,*) ' icon =',icon, ' length=',length
       if(icon == -1 ) then
          cTrack%pos%x =   cTrack%pos%x + 1.d-6*cTrack%w%x
          cTrack%pos%y =   cTrack%pos%y + 1.d-6*cTrack%w%y 
          cTrack%pos%z =   cTrack%pos%z + 1.d-6*cTrack%w%z
          call epbndry2(cmpn, length, icon)
          write(0,*) ' icon =',icon, ' length=',length
          cTrack%pos%x =   cTrack%pos%x + 1.d-6*cTrack%w%x
          cTrack%pos%y =   cTrack%pos%y + 1.d-6*cTrack%w%y 
          cTrack%pos%z =   cTrack%pos%z + 1.d-6*cTrack%w%z
          call epbndry2(cmpn, length, icon)
          write(0,*) ' icon =',icon, ' length=',length
       endif
       end
