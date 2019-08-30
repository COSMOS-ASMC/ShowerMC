#include "ZepicsBD.h"
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"

       type(epmedia):: media
      integer i, io
      real*8 Ee, prob, Eg, path, Ek 
      character*30 file
      character*24 name
      real(4):: rhoc
      integer  nevent, dummy, seed(2), posast

!      EpartialSC = 0.1  ! for test  ( Seltzer region max
!                        is lowered to this Ee)
      dummy = 9876531
      io = 10
      read(*, *) nevent, Ek, name, LPMeffect, Flpm
      call epgetRhoc(name, name, rhoc)
      file =  "../../../Data/Media/"//trim(name)
      call cerrorMsg(file,  1)

      Ee = Ek + masele
      
      open(io, file=file, action='read' ) 
      call epReadTab(io, media)
      media%rhoc = rhoc

      write(0,*) '# "Eg/Ek" "log10(//)" "sampled path(r%l)" '//
     *        ' " prob/r%l" "Ek"'    
      call  cmkSeed(dummy, seed)
      call rnd1i(seed) 

      do i = 1,  nevent
         call epBrSampP(media, Ee, prob, path)
         call epBrSampE(media, Ee, Eg)
         write(*,'(1p,5g13.4)') Eg/Ek,
     *       log10(Eg/Ek), path, prob, Ek
         
      enddo
      end

