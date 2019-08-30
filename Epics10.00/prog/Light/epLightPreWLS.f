      subroutine epLightPreWLS
! >>>>>>>>>>>>>>>>>>>>>>>light
      use modepLightPty
!c <<<<<<<<<<<<<<<<<<<<<light
      implicit none
!#include "Zcode.h"
!#include "ZepTrackp.h"
#include "ZepTrackv.h"
!c#include "Zcnfig.h"

!         these are given
!      Cn
!      MediaNo 
!      Media(MediaNo).rhoc
!       MaxPath 
!      if(Light > 0 ) then
!         cLcompNo 
!          if( cLcompNo > 0 ) then
!              cPtyNo = Lcomp(cLcompNo)%comInfoNo
!              Lcomp( cLcompNo )%refracN 
!       comInfo( cPtyNo )%WLS
      call epLightWLS( comInfo( cPtyNo )%WLS )
      end
