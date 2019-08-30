#include "ZepicsBD.h"
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"
      type(epmedia):: media
      integer i, io
      real*8 Emu, prob, Eo, path, Ek 
      character*90 file
      character*16 name
      integer  nevent, dummy, seed(2)
      integer kgetenv2, leng, icon

      dummy = 9876531
      io = 10
      read(*, *) nevent, Ek, name
      
      leng = kgetenv2("EPICSTOP", file)
      if( leng == 0 ) then
         write(0,*) ' EPICSTOP not yet given'
         stop
      endif
      file = file(1:leng)//"/Data/Media/"//name 
      call copenf(io, file, icon)
      if(icon /= 0 ) then
         write(0,*) file, 'cannot be opend' 
         stop
      endif

      Emu = Ek + masmu
      call epReadTab(io, media)
      call epMuNiEcheck(Emu, media)

      write(0,*) '# "Eo/Ek" "log10(//)" "sampled path(r.l)" '//
     *        ' " prob/r.l" "Ek"'    
      call  cmkSeed(dummy, seed)
      call rnd1i(seed) 

      do i = 1,  nevent
         call epmuNsmpP(media, Emu, prob, path)
         call epmuNsmpE(media, Emu, Eo)
         write(*,'(1p,5g13.4)') Eo/Ek,
     *       log10(Eo/Ek), path, prob, Ek
      enddo
      end
