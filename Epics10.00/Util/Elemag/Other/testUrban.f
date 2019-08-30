      implicit none
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zmass.h"
#include "Zcode.h"


!
!          test Urban's  dE/dx fluctuation
!
       type(epmedia):: media
       type(ptcl)::  aPtcl

      integer io, subcode, charge, code 
      real*8 E,  dedx, w, thick, de, avde, dedxfull
      integer i, icon
      real(8)::  RecoilKEmin
      
      character*120 file
      character*150 msg

      io = 10
      call  cerrorMsg(
     * "Enter ptcl code(2), charge(-1), "
     *  //
     * " w(100d-6), media file path('$EPICSTOP/Data/Media/BGO')", 1)
      

      code = 2
      charge = -1
      w = 100d-6
      file = '$EPICSTOP/Data/Media/BGO'
      
      read(*, *) code, charge, w, file
      call cerrorMsg(file,  1)
      call cerrorMsg('Enter Energy(1) and  thickness(.1 cm)',1)
      E = 1
      thick = 0.1
      read(*, *) E, thick
      RecoilKEmin = w      
      call copenf(io, file, icon)
      if( icon == 0 ) then
         write(0,*)  trim(file), ' opened '
      else
         stop  " file open error"
      endif
      write(0,*) ' entering epReadMTbl '
      call epReadTab(io, media)
      close(io)

      call epStern(media)
      write(0,*) ' entering epSetUrban'     

      call epsetUrban(media, media%urb)
      call epSetTcut(media, RecoilKEmin)  ! inside this, epResetUrban is called
      subcode = 0
      call cmkptc(code, subcode, charge, aPtcl)
      aPtcl%fm%p(4) = E

      if(aPtcl%code .eq. kelec) then
         call epdedxe(media, aPtcl, dedx, dedxfull)
      else
         call epdedxNone(media, aPtcl, dedx, dedxfull)
      endif

      
      avde = dedx*media%rho    ! GeV/cm
      write(msg, *) ' <dE/dx> in ',thick,' cm in ',media%name, 
     *   ' is ', avde*thick* 1000., ' MeV'
      call cerrorMsg(msg, 1)
      do i = 1, 5000
!         call epWriteUrbCnst(media%urb)
         call epUrban(media%urb, avde,
     *                thick, aPtcl, de)
         write(*,'(1p,3g13.5)')  de*1000., de/avde/thick,
     *        de*1000./thick
      enddo
      end




