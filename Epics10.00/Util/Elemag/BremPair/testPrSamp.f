#include "ZepicsBD.h"
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"
!            test pair creation sampling ; energy is arbitrary 
!
       type(epmedia):: media
      integer i, io
      real*8 Ee, prob, Eg, path
      character*50 file
      character*24 name

      real(4):: rhoc
      integer  nevent, dummy, seed(2)

      dummy = 9876533

      io = 10
      read(*, *) nevent, Eg, name, LPMeffect
      call epgetRhoc(name, name, rhoc)
      file = "../../../Data/Media/"//trim(name) 
      call cerrorMsg(file,  1)

      open(io, file=file, action='READ') 
      call epReadTab(io, media)
      media%rhoc = rhoc

      write(0, '(a)') '# Ee/Eg (Ee-m)/(Eg-2m) sampled Path(r.l)'//
     *     ' prob/r.l'
      call  cmkSeed(dummy, seed)
      call rnd1i(seed)
      do i = 1,  nevent
         call epPrSampP(media, Eg, prob, path)
         call epPrSampE(media, Eg, Ee)
         write(*,'(1p,4g12.4)')
     *    Ee/Eg, (Ee-masele)/(Eg-2*masele), path, prob
!           Eg = Ee1 + Ee2
          Ee = Eg - Ee
          write(*,'(1p,4g12.4)')
     *    Ee/Eg, (Ee-masele)/(Eg-2*masele), path, prob
!
      enddo
      end
