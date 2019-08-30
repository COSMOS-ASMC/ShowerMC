#include "ZepicsBD.h"
      implicit none
!
!     to get total XS of brems with some xmin for a given media
!     at entire energies (1keV to 10 GeV)
!    x(min) is Eg/Ee not Eg/Ek (Ek=Ee-Me)
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
#include "ZepTrackp.h"

      character*50  file
      character*16  basemedia
      integer nerg
      integer io, i, norm
      real*8  Ee, x, f, epBrgenex, xmax
      real*8  xmin, tprob, xx, fxx, Nc
      real*8  E1, E2, step, Ek

      io = 10

      read(*,*)  basemedia, xmin, E1, E2, step, norm 


      file ="../../../Data/BaseM/"//basemedia

      open(io, file=file, action ='read')
      call epBPgeneini(io, media)
      media%cnst%BrEemaxS =min( 10.0001d0, E2)   ! force
      write(*,'(a,a,a,1p,g13.4,a,i2)' )
     *      '# ', trim(basemedia),'  xmin =', xmin, 
     *      ' norm=', norm
      write(*,*)
     * "# norm=1: /r%l, 2: mb/ingredient,"//
     * " 3: /(g/cm2), 4:/cm, 5: area norm "
      write(*,'(a,1p,4g13.4)') "#  Ee- masele, tprob*NC, tprob, NC "
      Ek = E1
      Ee = Ek + masele
!      media.cnst.BrEgminS=xmin
      media%cnst%BrEgminS  =  xmin
      write(0,*) ' media%cnst%BrEgminS =', media%cnst%BrEgminS
      xmax = 1. - masele/Ee
      write(0,*) ' xmax =', xmax

      do while( Ek <  media%cnst%BrEemaxS ) 
!         xmin = media.cnst.BrEgminS
         Ee = Ek + masele
         Eeme = Ee/masele   ! this is used inside brems func.
!            xmin = media.cnst.BrEgminS  !  neglect this setting
         call epBrgeneTX(xmin, xmax, tprob)
         if(norm .eq. 5) then
            Nc = 1./tprob
         elseif(norm .eq. 1) then 
            Nc=media%mbtoPX0    ! prob/X0
         elseif( norm .eq.  2) then
            Nc = 1
         elseif(norm .eq.  3) then
            Nc = media%mbtoPgrm
         elseif(norm .eq. 4) then
            Nc = media%mbtoPcm
         else
            call cerrorMsg('input error for norm',0)
         endif
         write(*,'(1p,4g13.4)')  Ek, tprob*NC, tprob, NC
         if(Ek >   media%cnst%BrEemaxS*0.5) then
            Ek = Ek * 10.0d0**0.01
         else
            Ek = Ek * 10.0d0**0.1
         endif
      enddo 
      end
