#include "ZepicsBD.h"
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!        
!          test muon pair creation sampling 
!
       type(epmedia):: media
      integer i, io, icon, norm, nevent, dummy, seed(2)
      real*8 Emu, prob, Et, path, normf
      real*8 epmuPairVmn
      character*120 file
      character*16 name
      integer kgetenv2, leng
      character*120 msg
      real*8 xmin

      io = 10
      dummy = 987631

      read(*, *) nevent, Emu, name
      leng = kgetenv2("EPICSTOP", file)
      if( leng == 0 ) then
         write(0,*) 'EPICSTOP not yet given'
         stop
      endif
      file = file(1:leng)//"/Data/Media/"//name
      call copenf(io, file, icon)
      if(icon .ne. 0) then
         write(0,*) file, "cannot be opend"
         stop
      endif

      call epReadTab(io, media)
      call  cmkSeed(dummy, seed)
      call rnd1i(seed)

      xmin = max( media%cnst%muPrVmin, epmuPairVmn(Emu))
!////////////
!      xmin = epmuPairVmn(Emu)  ! this is to see entire function region
!/////////////
      write(0,'(a)')
     *   '# "E_pair/Emu" "log10(E_pair/Emu)" "sampled path(r.l)" '//
     *   ' " prob/r.l" "Emu'    
      do i = 1,  nevent
         call epmuPrsmpP(media, Emu, prob, path)
         call epmuPrsmpE(media, Emu, Et)
         write(*,'(1p,5g13.4)') Et/Emu,
     *       log10(Et/Emu), path, prob, Emu
      enddo

      end

