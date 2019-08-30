#include "ZepicsBD.h"
      implicit none
!
!      
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
#include  "Zepdebug.h"      
      real*8 length, ux, uy, uz
      integer icon, i, np, nc
!      integer debug
       type(epPos):: posw
       type(epDirec)::  dirw, dirl

       type(epPos)::  org, posl
       type(ep3Vec)::  abc
      integer ir(2), cmpn
      logical boundary/.true./, iso/.true./

      character*100 dsn1        ! input config data file path

      ir(1)=2222333
      ir(2)=3445679
      MsgLevel = 1
      np = 15000
      cmpn = 1
      call cerrorMsg(
     * 'Enter config file path)', 1)
      read(*,'(a)')  dsn1
      call cerrorMsg('boundary test or inout test', 1)
      call cerrorMsg(
     * 'boundary test(=t), number of points(=15000),'//
     *     'isotropic angles(t), 2ir,  Comp# to test(1)', 1)
       
       read(*, *) boundary, np, iso, ir,  cmpn

       MediaDir(1) = '$EPICSTOP/Data/Media'
       call eprcnf(dsn1)

       call rnd1i(ir)
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
       do while(.true.)
          call  rndc(ux)
          call  rndc(uy)
          call  rndc(uz)
          cTrack%pos%x = abc%x * ux* 1.2 + org%x-abc%x*.1
          cTrack%pos%y = abc%y * uy* 1.2 + org%y-abc%y*.1
          cTrack%pos%z = abc%z * uz* 1.2 + org%z-abc%z*.1

          if(iso) then
             call episoDirec(cTrack%w) ! isotroic angle
          else
             call epvertical(cTrack%w)   ! vertical to x-y plane, etc
          endif

          if(boundary) then

             call epw2l(cmpn, cTrack%pos, posl)
             cTrack%pos = posl
             call epw2ld(cmpn, cTrack%w, dirl)
             cTrack%w = dirl
             call epbndry2(cmpn, length, icon)
!                 epbNew is called from epbndry2
          else
             call eppos2cn(0, cTrack, nc)
!                 epsNew is  called  from eppos2cn
             if(nc .eq. cmpn) then
                icon = 0
             else
                icon = 1
             endif
          endif

          if(.not. boundary) then
             if(i .lt. np) then
                i = i+ 1
                write(*,*)
     *               sngl(cTrack%pos%x),
     *               sngl(cTrack%pos%y),
     *               sngl(cTrack%pos%z), icon, nc
             else
                exit
             endif
          else
!               convert to wolrd
             if(icon .ne. -1) then
                if(i .lt. np ) then
                   i = i+1
                   call epl2w(Cn, cTrack%pos, posw)
                   cTrack%pos = posw
                   call epl2wd(Cn, cTrack%w, dirw)
                   cTrack%w = dirw
                   write(*,*)
     *             sngl(cTrack%pos%x +length*cTrack%w%x),
     *             sngl(cTrack%pos%y +length*cTrack%w%y),
     *             sngl(cTrack%pos%z +length*cTrack%w%z),icon
                else
                   exit
                endif
             endif
          endif
       enddo
       end
      subroutine epvertical(direc)
      implicit none
#include "ZepDirec.h"
       type(epDirec)::  direc

      real*8 u, ux  

      call rndc(u)
      call rndc(ux)
      if(ux .gt. .5) then
         ux = 1.
      else
         ux = -1.
      endif

      if(u .lt. 0.3333d0) then
         direc%x = ux
         direc%y = 0.
         direc%z = 0.
      elseif(u .lt. .6666d0) then
         direc%x = 0.
         direc%y = ux
         direc%z = 0.
      else
         direc%x = 0.
         direc%y = 0.
         direc%z = ux
      endif
      end

