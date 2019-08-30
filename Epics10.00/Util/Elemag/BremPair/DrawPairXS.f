      implicit none
!
!         compute pair total xs.  prob/r.l, prob/(g/cm2), prob/cm , mb/atom
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      character*40  file
      integer io, result

      real*8  Eg, Egmax
      real*8  tprob, p


      io = 10

      call cerrorMsg(
     * "Enter upper Eg(10 GeV) and media file path( such as"//
     *  " '../../Data/Media/Pb' )", 1)

      Egmax = 10.01
      file= '../../Data/Media/Pb'

      read( *, * ) Egmax, file

      call copenf(io, file, result)

      if( result .eq. 0 ) then

         call epReadTab(io, media)

         write(*,*) "#  'prob.  of Pair Creation'"
         write(*,*) 
     *   "#  'Eg(MeV)'  ' /r%l' '/(g/cm2)' '/cm' 'mb/atom'"
         Eg = 2.2*masele
         do while ( Eg .lt. Egmax)   
            call epPrSampP(media, Eg, tprob, p)
            write(*, *) sngl(Eg*1000.d0), 
     *       sngl(tprob), sngl(tprob/media%X0g),
     *       sngl(tprob/media%X0), sngl(tprob/media%mbtoPX0)
            Eg = Eg * 10.0d0**(0.025d0)
         enddo
      else
         write(0,*) ' media file error '
      endif
      end
