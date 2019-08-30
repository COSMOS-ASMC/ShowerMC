      implicit none
#include "Zmedia.h"

       type(epmedia):: media

 
      call epReadMTbl(5, media)
      call epGetEffZA(media)
      call epSetSTblCns(media, media%cnst)
!        reset vmin
      call epresetmuVmin(media%cnst%muNVmin)
      call epCreMuNTab(media, media%cnst)
      end



