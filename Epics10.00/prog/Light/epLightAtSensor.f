!       photo electron generation at the sensor

      subroutine epLightAtSensor
      use modepLightPty
      use modepLightCounter
      implicit none
#include  "ZepTrackv.h"


      real(8):: qeff  ! Qeff of the  photo sensor
      real(8):: u     ! 0< u < 1; uniform random number

      integer:: id  ! scinti or Ceren or ..
!///////
!      call Lcompchk(' bb ', cLcompNo)
!////////
      if( comInfo( cPtyNo )%QeffAF == 0 ) then
         ! use constant Qeff
         qeff = comInfo( cPtyNo )%Qeff
      else
         ! get  w.l dependent  Qeff
         call csampAFIntp( comInfo( cPtyNo )%QeffAF,
     *       cTrack%wl, qeff)
!//////////////
!         write(*,*) 'qeff ', cTrack.wl, qeff
!/////////////
      endif

      call rndc(u)
      if( u < qeff ) then
         ! count up p.e #
         if( .not. allocated( sensor( cLcompNo )%pe ) ) then
            allocate(  sensor( cLcompNo )%pe(1:5) ) 
!!            allocate(  sensor( cLcompNo )%pet ) 
            sensor(cLcompNo)%pe(:) = 0.
         endif
         id = cTrack%p%subcode
!/////////////
!         write(0,*) ' id =',id, ' cLcompNo=', cLcompNo, 
!     *     ' w=',cTrack.wgt, ' code=',cTrack.p.code
!//////////       
         sensor( cLcompNo )%pe(id) =
     *           sensor( cLcompNo )%pe(id)  + cTrack%wgt
      endif

      end
