      implicit none
#include "Zmedia.h"

       type(epmedia):: media

 
      call epReadMTbl(5, media)
      call epGetEffZA(media)
      call epSetSTblCns(media, media%cnst)
!         rest vmin
      call epresetmuVmin(media%cnst%muPrVmin)
      call epCreMuPTab(media, media%cnst)
      end



