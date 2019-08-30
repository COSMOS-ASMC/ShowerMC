      implicit none
!
!         to compute photo electric effect XS.
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      character*40  file
      integer io, result

      real*8  Eg
      real*8  tprob, p


      io = 10

      call cerrorMsg(
     * "Enter media file path( such as  Pb  BGO  etc )", 1)

      read( *, * ) file
      file = "$EPICSTOP/Data/Media/"//trim(file)

      call copenf(io, file, result)

      if( result .eq. 0 ) then

         call epReadTab(io, media)

         write(*,*) "#  'Prob. of PhotoElectric Effect'"
         write(*,*) "#  'Eg(MeV)'  ' /r.l' '/(g/cm2)' '/cm' 'mb/atom'" 
         Eg = 0.01d-3            !  10 keV.
         do while ( Eg .lt. 11.d-3)   ! to 10 MeV
!            call epphotoEp(media.pe, Eg, tprob, p) !  <  v8.0
            call epphotoEp(media, Eg, tprob, p)
            write(*, *) sngl(Eg*1000.d0), 
     *       sngl(tprob), sngl(tprob/media%X0g),
     *       sngl(tprob/media%X0), sngl(tprob/media%mbtoPX0)
            Eg = Eg * 10.0d0**(0.025d0)
         enddo
      else
         write(0,*) ' media file error '
      endif
      end
