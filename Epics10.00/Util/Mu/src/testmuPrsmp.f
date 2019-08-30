      implicit none
#include "Zmedia.h"
#include "Zmass.h"
!        
!          test muon pair creation sampling 
!
       type(epmedia):: media
      integer i, io, icon, norm
      real*8 Emu, prob, Et, path, normf
      character*20 file
      character*120 msg
      io = 10

      call  cerrorMsg(
     * " Enter muon Energy(in GeV) and media file path( "//
     * "such as '../../Data/Media/Pb')",1)
      file = '../../Data/Media/Pb'
      read(*, *) Emu, file
      call cerrorMsg(file,  1)

      call copenf(io, file, icon)
      if(icon .ne. 0) then
         call cerrorMsg('file cannot be opend',0)
      endif
      call epReadTab(io, media)

 10      continue
      call cerrorMsg('Path length unit',1)

      call cerrorMsg(' 1--> r.l',1)
      call cerrorMsg(' 2--> (g/cm2)',1)
      read(*, *) norm
      if(norm .ne.1 .and.  norm .ne.2)  goto 10
      
      write(msg,*) '# "muon pair creation sampling Emu=',Emu,' GeV" '
      write(*, *) msg
      if(norm .eq. 1) then
         normf = 1.
         write(msg,*)
     *   '# "E_pair/Emu" "log10(E_pair/Emu)" "sampled path(r.l)" '//
     *   ' " prob/r.l"'    
      else
         normf =  media%mbtoPgrm/media%mbtoPX0
         write(msg,*)
     *   '# "E_pair/Emu" "log10(E_pair/Emu)" "sampled path(g/cm2)" '//
     *   ' " prob/(g/cm2)"'    
      endif
      write(*,*) msg
      do i = 1,  50000
         call epmuPrsmpP(media, Emu, prob, path)
         prob = prob *normf
         path = path/normf
         call epmuPrsmpE(media, Emu, Et)
         write(*,*) sngl(Et/Emu),
     *      sngl(  log10(Et/Emu) ), sngl(path), sngl(prob)
      enddo
      end


