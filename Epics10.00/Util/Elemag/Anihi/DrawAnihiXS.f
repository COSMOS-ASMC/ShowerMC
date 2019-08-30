      implicit none
!
!         compute anihilation total xs.  prob/r.l, prob/(g/cm2), prob/cm , mb/atom
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      character*10 name
      integer io, result
      character*120 file
      real*8  Ee, Eemax, EeT
      real*8  tprob, p


      io = 10

      call cerrorMsg(
     * "Enter upper Ee(KE ~1 GeV) and media Name(W etc) ",
     *   1)

      read( *, * ) Eemax, name
      file ="$EPICSTOP/Data/Media/"//trim(name)
      call copenf(io, file, result)

      if( result .eq. 0 ) then

         call epReadTab(io, media)

         write(*,*) "#  'prob.  of Anihi'"
         write(*,*) 
     *   "#  'Ee(KE. MeV)'  ' /r.l' '/(g/cm2)' '/cm' 'mb/atom'"
         Ee = 1.d-6
         do while ( Ee .lt. Eemax)   
            EeT = Ee + masele
            call epanihip(media, Eet, tprob,  p)
            write(*, '(1p, 5g14.4)') Ee*1000.d0, 
     *       tprob, tprob/media%X0g, 
     *       tprob/media%X0,tprob/media%mbtoPX0
            Ee = Ee * 10.0d0**(0.05d0)
         enddo
      else
         write(0,*) ' media file error '
      endif
      end
