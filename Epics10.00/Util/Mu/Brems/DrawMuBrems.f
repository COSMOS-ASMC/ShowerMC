      implicit none
!
!         to write data needed to draw muon brems function 
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h" 

      character*50  file
      character*8 basemedia

      integer io, i, norm, ie, icon
!      real*8  x, f, epmuBremS, xmax
      real*8  x, f, epmuBrS, xmax, E1
      real*8  xmin, tprob, Nc
      integer kgetenv2, leng


      io = 10
!      call cerrorMsg(' ', 1)
!      call cerrorMsg('1) In unit of / r.l (1)',  1)
!      call cerrorMsg('2) In unit of mb/ingredient', 1) 
!      call cerrorMsg('3) In unit of /(g/cm^2)', 1)
!      call cerrorMsg('4) In unit of /cm', 1)
!      call cerrorMsg(
!     * '5) Area-normalization '//
!     * ' for comparison with M.C data such as v vs vdN/dv/Nt  ',
!     * 1 )

      

      norm = 1


      read(*,*) norm, basemedia, E1

         
      leng =  kgetenv2("EPICSTOP", file)
      if( leng  == 0 )  then
         write(0,*) "EPICSTOP is not given yet"
         stop
      endif

      file =file(1:leng)//"/Data/Media/"//basemedia

      call copenf(io, file, icon)
      if(icon .ne. 0) then
         call cerrorMsg('media file path invalid',  0)
      endif
!          call epmuBPNgeneI(io, media) ! obso next is ok
      call epReadTab(io, media)
      call epMuBrEcheck( E1, media)
      xmin = media%cnst%muBrVmin
      Emu = E1

      x = xmin
!         xmax = epmuvmax( Emu )             
!          cannot be determined here since it is dependent
!          on each element. So we set it  1.0 tentatively
!          actual value is set in  epmuBrS
!          and  if x> actual xmax, function value is set to 0 for
!          that atom.
!
      xmax = 1.d0
!      if(norm .eq. 5) then
         call eptotcmuB(xmin, xmax, tprob)
!      else
!         tprob =  1.
!      endif

         write(0,*) ' xmin, tprob=',xmin, tprob
         
      if(norm .eq. 5) then
         Nc = 1./tprob
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

      do while (x .lt.  xmax)
!            f = epmuBremS(x)/x　　！　　! obso
         f = epmuBrS(x)/x * Nc
         if(f .gt. 0.) then
            write(*, '(1p,5g13.4)')
     *         x, f*x, Emu, Nc, tprob
         else
            exit
         endif
         if(x .gt.  xmax-0.5d0 ) then
            x = x + 1.d-3
         else
            x = x * 10.**0.02
         endif
      enddo
 5    continue
      write(*,*)
      end

      
