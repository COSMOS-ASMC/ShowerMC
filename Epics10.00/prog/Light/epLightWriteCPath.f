      subroutine epLightWriteCPath(aTrack)
        ! charge ptcl track is put in the OutPrimaryFile
        ! for use as a +primary 
      implicit none
#include "Zcode.h"
#include "ZepTrack.h"
#include "ZsepManager.h"

       type(epTrack)::   aTrack  ! input charged particle track
                    ! for Cerenkov light generation. .wl is path length
                    ! in this case.
      
         !  The output must be the same format as the psudo ptcl of
         !  energy deposit for scintillation light:
         !   code, subcode charge, E, pos, dir, Cn                                
         !   for Edepo:   subcode=0 charge=0. dir=(dx,dy,dz)= cell size  
         !   for chgPath:         dir=dl(dx,dy,dz)= directed path length 
       type(epTrack)::  path

      integer:: icon
      real(8):: betax, wl1, wl2, E1, E2, n1, n2

      path = aTrack
      call epLightCerenTh(path, wl1, wl2, E1, E2, n1, n2, betax,
     *           icon)

      if(icon == 0 ) then
!/////////////
!         if(OutPrimEff == 0 ) then
!            write(0,*) ' OutPrimEff=0, Out1ry=',Out1ry
!            stop
!         endif
!/////////////
         if(Out1ry > 0 ) then
            write(OutPrimEff, 
     *      '(i6,2i4,1p,4g15.7,3g18.9,2g12.4,0p,i6)' )
     *      kchgPath, path%p%subcode, path%p%charge,
     *      path%p%fm%p(4), path%pos, path%w,
     *      path%wl, path%p%mass,  path%cn    !  wl is path lenght in cm in this case.
         else
            write(OutPrimEff)
     *      kchgPath, int(path%p%subcode), int(path%p%charge),
     *      path%p%fm%p(4), path%pos, path%w,
     *      path%wl, path%p%mass,  path%cn 
         endif
      endif
      end
