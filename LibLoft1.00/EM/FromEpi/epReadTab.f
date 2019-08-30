      subroutine epReadTab(io, media)
!
!         read sampling table  created by epCreBrTab and epCrePrTab
!
      implicit none
#include "Zmedia.h"
!
       type(epmedia):: media   ! output. data is read and stored here

      integer io  ! input.  logical device

      character*10 rec

!     
!      At the top of table, it is assumed that a copy of
!    basic media table exists.
      call epReadMTbl(io, media)

!        skip to the next to the separator "#-#-#-#-#-"
      do while (.true.)
         read(io, '(a)') rec
         if(rec(1:10) .eq. '#-#-#-#-#-') goto 10
      enddo
 10   continue
!
!          media data  ; this must be computed every time as below
!         so don't read it
      media%rhoc = 1.0 ! when no config data given at the test time
!                      (e.g testdEdx.f in Util/Elemag, rhoc is not
!                       given so we set standard here)
!       compute basic const
      call epGetEffZA(media)
!     
!      
      call epExpot(media)

!    
!           Sternheimer's consts and to put them in the table
!               to give a parameter value instead of 1d-4,
!      call epsetUrbn(media,  media.urb)

!      (many or all of these will be reset by reading
!       the table later)
      call epSetSTblCns(media, media%cnst)
      call epSetPhotoE(media, media%pe)  ! this is also not essential
                                         ! to create sampling talbe

!      
!     ------------------------------  brems table 
!          Seltzer table
      call epRdBrSTblS(io, media%cnst, media%tbl)
!
!         GeV region upto LPM  brems
      call epRdBrSTbl1(io, media%cnst, media%tbl)
!        LPM region    
!           brems
      call epRdBrSTblH(io, media%cnst, media%tbl)
!     ------------------------------  pair table 
!                no LPM region
      call epRdPrSTbl1(io, media%cnst, media%tbl)
!                LPM region
      call epRdPrSTblH(io, media%cnst, media%tbl)
!    
!     *********** muon table ************
      call epRdmuTab(io, media)
!          const used for muon interactions
      call epSetMu(media, media%mu)                                         
      
      end
