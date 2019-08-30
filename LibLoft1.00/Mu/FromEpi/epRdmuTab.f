      subroutine  epRdmuTab(io, media)
!
!         read sampling table  for muons ( created by epCreMuTab)
!    If the first data is non-existent, muon interactions (brems,
!    pair, nuclear interaction) will neglected. 
!
      use modMuNucontrol
      implicit none
#include "Zmedia.h"
!  #include "ZepTrackp.h"
!
       type(epmedia):: media   ! output. data is read and stored here

      integer io  ! input.  logical device
      character*30  sep 

      media%mu%MuNI = MuNI
      media%mu%MuBr = MuBr
      media%mu%MuPr = MuPr
!             not needed now; since all media has same pair  spectrum
!      media.mu.NoOfErgNodeForMuPair = 0
!         set Erg
      media%mu%muPairErg(:) =
     *    (/ 3., 5., 8., 10., 17., 35., 80., 150.,  500., 1500./)
      media%mu%logmuPairErg(:) = log(media%mu%muPairErg(:) )
!         skip to  the muon data
      do while (.true.)
         read(io,  '(a)', end=100) sep
!         write(*,*) sep
         if(sep(5:15) .eq. '-----------') goto 20
      enddo
 20   continue
!         data is from pair, brems, and n.i (order of higher sigma)
      call epRdmuPrSTbl(io, media%cnst, media%tbl)
      call epRdmuBrSTbl(io, media%cnst, media%tbl)
      call epRdmuNSTbl(io, media%cnst, media%tbl)
      return
!    
 100  continue
!          table dose not exist
!           set flag to neglect muon interaction
      media%mu%MuNI = 0
      media%mu%MuBr = 0
      media%mu%MuPr = 0
      call cerrorMsg('muon interaction table not exist', 1)
      end
