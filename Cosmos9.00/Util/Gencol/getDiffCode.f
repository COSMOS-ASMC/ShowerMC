#include "Zintmodel.h"
      subroutine getDiffCode(nw, difcode)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
!
! general process information; only for dpmjet3
      INTEGER IPROCE,IDNODF,IDIFR1,IDIFR2,IDDPOM,IPRON
      COMMON /POPRCS/ IPROCE,IDNODF,IDIFR1,IDIFR2,IDDPOM,IPRON(15,4)
!      IPROCE
!             1 non-diffractive inelastic
!             2 elestic 
!             3 quasi elestic vector meson prod. (photon)
!             4 central diffraction
!             5 single diff. ptcl 1
!             6 //           ptcl 2
!             7 double diff. 
!             8 direct photo-hadron
! For moore detail, see manual in Documents/CosmicRays/phojetShort.pdf
!               say, IDIFR1 classifies IPROCE=5



!tp from epos;   same as in qgsjet01:  typevt=1 nsd 2 dd 3 sd 4 cd

      integer ng1evt,ng2evt,ikoevt
      real    rglevt,sglevt,eglevt,fglevt,typevt
      common/c2evt/ng1evt,ng2evt,rglevt,sglevt,eglevt,fglevt,ikoevt
     *,typevt            !in epos.inc          
!     typevt (1 - NSD, 2 - DD, 3 - SD, 4 - CD)                       
!                                ine,  ela, ine,  CD SD SD DD  ine
      integer,save::dpm2typevt(8)=(/6,    5,   6,    4, 3, 3, 2, 6/)

      integer nw, difcode(20)   ! 20 is the max woonded nucleon in 
                           ! sibyll
      if( ActiveMdl == 'sibyll') then
#ifdef  SIBYLL
         call sibylGetDiffCode(nw, difcode)
#endif
#ifdef  QGSJET1
      elseif( ActiveMdl == 'qgsjet1') then
         difcode(1) = typevt
         nw = 1
#endif
      elseif(ActiveMdl == 'dpmjet3') then
        nw = 2
        difcode(1) =dpm2typevt(IPROCE)
      elseif( ActiveMdl == 'qgsjet2') then
         nw = 1
         difcode(1) = typevt
      else
         nw = 1
         difcode(1) = 0   ! not know
!         write(0,*) ' ActiveMdl=',ActiveMdl, ' invalid'
!         write(0,*) ' in getDiffCode.f'
!         stop
      endif
      end
