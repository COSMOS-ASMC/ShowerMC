!     define next for draw only cross-sec at only one Emu
!     (used for drawing curve and histogram)
#define  ONLYONEE
      implicit none
!
!         to write data needed to draw muon Nuclear interaction function 
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZmuBPNgene.h"

      
      character*50  file
      character(20):: basemedia
      integer nerg
      integer io, i, norm, ie, icon, leng
      real*8  x, f, xmax, epmuNS, Nc
      
      real*8  xmin, tprob, vmin
      integer,external:: kgetenv2

      parameter (nerg=10)
      real*8  EmuA(nerg)

      io = 10
      EmuA(:) = 0.
#if  ! defined (ONLYONEE)
      call cerrorMsg(' ', 1)
      call cerrorMsg('1) In unit of / r%l (1)',  1)
      call cerrorMsg('2) In unit of mb/ingredient', 1) 
      call cerrorMsg('3) In unit of /(g/cm^2)', 1)
      call cerrorMsg('4) In unit of /cm', 1)
      call cerrorMsg(
     * '5) Area-normalization '//
     * ' for comparison with M%C data such as v vs vdN/dv/Nt  ',
     * 1 )
      read(*, *)  norm
      call cerrorMsg(
     *  "Enter a basic media name ( such as PWO etc )",1)
      read(*, *) file
!      call cerrorMsg("Enter vmin", 1)
!      read(*,*)  vmin
      call cerrorMsg(
     *     '  Energy of the muon in (GeV) upto 10 (with /)', 1)
      read(*, *) EmuA
#else
      write(0,*) 'Enter  norm media  Emu'
      read(*, *) norm, basemedia, Emu
      EmuA(1) = Emu
#endif
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
!
      call epReadTab(io, media)
!      call epmuBPNgeneI(io, media)

      vmin =media%cnst%muNVmin
      
      do ie = 1, nerg
         Emu = EmuA(ie)
         if(Emu == 0.)  exit
         xmin =max( vmin,  masrho/Emu)
         xmax = 1.- masmu/Emu
         call eptotcmuN(xmin, xmax, tprob)
         if(norm .eq. 5) then
            Nc = 1./tprob
         elseif( norm == 1 ) then
            Nc=media%mbtoPX0     ! prob/X0
         elseif( norm .eq.  2) then
            Nc = 1
         elseif(norm .eq.  3) then
            Nc = media%mbtoPgrm
         elseif(norm .eq. 4) then
            Nc = media%mbtoPcm
         else
            call cerrorMsg('input error for norm',0)
         endif
         x = xmin
         do while (x  <  xmax)
            f = epmuNS(x)*Nc
            if(f .gt. 0.) then
               write(*,'(1p,6g13.4)')
     *         x, f, Emu, Nc, tprob, Nc*tprob
            elseif( x .ge. .8) then
               exit
            endif
            if(x .gt.  xmax-0.85d0 ) then
               x = x + 2.d-3
            else
               x = x * 10.**0.02
            endif
         enddo
 5       continue
         write(*,*)
      enddo
 10   continue
      end
