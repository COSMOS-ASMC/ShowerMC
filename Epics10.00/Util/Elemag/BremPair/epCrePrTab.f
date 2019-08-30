      implicit none
#include "Zmedia.h"
#include "ZepTrackp.h"
       type(epmedia):: media
      integer::icon, how
      character(100)::file
      how =0 
      read(*, '(a)')  file

      call epBPgeneini(file, media, how)
!
!        < LPM region
!           pair
      call epCrePrSTbl1(media, media%cnst)
!           pair LPM
      LPMeffect =.true.
      call epCrePrSTblH(media, media%cnst)
      end

