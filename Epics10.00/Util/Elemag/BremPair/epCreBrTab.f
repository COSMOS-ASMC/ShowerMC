#include "ZepicsBD.h"
      implicit none
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zcode.h"
!
!                      Creating brems sampling table
!
       type(epmedia):: media

      integer::icon, how
      character(100)::file   ! basic media file
      read(*, '(i2, a)')  how, file
      if( how < -1 .or. how >1 )  then
         write(0,*) ' how=',how, ' is wrong'
         stop
      endif

      call epBPgeneini(file, media, how)
!      write(0,*) ' entering epExpot'
      call epExpot(media)
!           next 1.d-4 is dummy; this is to compute 
!           Sternheimer's consts and to put them in the table
!      write(0,*) ' entering epStern'
!      call epStern(1.d-4, media)  <v9.154
      call epStern( media)    ! but not updated until 9.160
!      write(0,*) ' entering epSetPhotoE'

      call epSetPhotoE(media, media%pe)  ! this is also not essential ?
                                         ! to create sampling talbe
                                         ! but for putting consts in the table

!      write(0,*) ' entering epwtmedia'
      call epwtmedia(media)  ! this writes some basic data
!                            ! must be skipped when reading
      call epDisableLPM
!             Brems table by Seltzer data region
!      write(0,*) ' entering epCreBrSTblS'
      call epCreBrSTblS(media, media%cnst)

!         GeV region upto complete screening
!      write(0,*) ' entering epCreBrSTbl'
      call epCreBrSTbl1(media, media%cnst)

!      LPM   from EeminLPM 
      call epAbleLPM
!       
!      write(0,*) ' entering epCreBrSTbH'
      call epCreBrSTbH(media, media%cnst)
      end
      subroutine  epDisableLPM
      implicit none
#include "ZepTrackp.h"
      LPMeffect=.false.
      Flpm= 1.0
      end
      subroutine  epAbleLPM
      implicit none
#include "ZepTrackp.h"
      LPMeffect=.true.
      Flpm= 1.0
      end
