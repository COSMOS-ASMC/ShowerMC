      implicit none
!
!         compute brems total xs.  prob/r.l, prob/(g/cm2), prob/cm , mb/atom
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"


      character*40  file
      integer io, result

      real*8  Ee, Eemax
      real*8  tprob, p


      io = 10

      call cerrorMsg(
     * "Enter upper Ee(100 GeV) and media file path( such as"//
     *  " '../../Data/Media/Pb' )", 1)

      Eemax = 100.01
      file= '../../Data/Media/Pb'

      read( *, * ) Eemax, file

      call copenf(io, file, result)

      if( result .eq. 0 ) then

         call epReadTab(io, media)

         write(*,*) "#  'prob.  of Brems'"
         write(*,*) 
     *   "#  'Ee(MeV)'  ' /r%l' '/(g/cm2)' '/cm' 'mb/atom'"
         Ee = 100.d-6
         do while ( Ee .lt. Eemax)   
            call epBrSampP(media, Ee, tprob, p)
            write(*, *) sngl(Ee*1000.d0), 
     *       sngl(tprob), sngl(tprob/media%X0g),
     *       sngl(tprob/media%X0), sngl(tprob/media%mbtoPX0)
            Ee = Ee * 10.0d0**(0.05d0)
         enddo
      else
         write(0,*) ' media file error '
      endif
      end
