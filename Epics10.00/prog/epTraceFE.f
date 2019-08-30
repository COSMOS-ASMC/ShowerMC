#include "ZsubstRec.h"
#include "Zepcondc.h"
      subroutine epTraceFE
!        front end for trace
      implicit  none
#include  "ZepTrackp.h"
#include  "ZepTrackv.h"
#include  "Zcode.h"
#include  "ZepManager.h"
      if(Trace) then
         if(IoTrace .lt. 0 .or. (IoTrace .gt. 0 
     *        .and. cTrack%p%charge .ne. 0)) then
            call epTrace
         endif
      endif
      end

!     ******************
      subroutine epTrace
!     ******************
!         take trace of paritcles
       implicit none
#include  "ZepTrackp.h"
#include  "ZepTrackv.h"
#include  "Zcode.h"
#include  "ZepManager.h"
! #include  "Zcnfig.h"
#if defined NEXT486
#define IMAG_P dimag
#elif defined PCLinux
#define IMAG_P dimag
#else
#define IMAG_P imag
#endif
      

!          standard Trace information
       type(epPos)::  posbw, posw
      real*8 xxx/-1.d37/, yyy/-1.d36/, zzz/1.d34/ 
      integer kkk/-1000/, chkchg/-1000/
      real*8 erg
      integer k
      save xxx, yyy, zzz, kkk, chkchg

!
!          Trace informaion  ccc if for putting only electorns
!
!             convert to world coord. starting point of the segment

      erg= Move%Track%p%fm%p(4)-Move%Track%p%mass
      if(erg == 0.) then  !avoid 0 energy in tracking display v9.163
         erg = max((cTrack%p%fm%p(4)-cTrack%p%mass), 0.d0)
      endif
      if( TraceErg(1) .le. erg .and.  erg .le. TraceErg(2) ) goto 100
      if( TraceErg(3) .le. erg .and.  erg .le. TraceErg(4) ) goto 100
      if( TraceErg(5) .le. erg .and.  erg .le. TraceErg(6) ) goto 100
      return  !  ***************
 100  continue
      call epl2w(Cn, cTrack%pos, posbw)
      
      k = cTrack%p%code
!cc      if(k .eq. kelec) then
         if(kkk .ne. k .or. xxx .ne. posbw%x .or. 
     *                      yyy .ne. posbw%y .or.
     *                      zzz .ne. posbw%z .or. 
     *                   chkchg .ne. cTrack%p%charge ) then
            if(xxx .ne. -1.d37) then
               write(abs(IoTrace), *)
               write(abs(IoTrace), *)
            endif
            if(TimeStruc) then
               write(abs(IoTrace),
     *          '(3g16.8, i4, g11.4, i4, g16.8)')
     *            posbw%x, posbw%y, posbw%z, k,
     *            erg,
     *            cTrack%p%charge,
     *            cTrack%t
            else
               write(abs(IoTrace),'(3g16.8,i4,g11.4, i4,1x, a)')
     *            posbw%x, posbw%y, posbw%z, k,
     *            erg,               
     *            cTrack%p%charge,  Media(MediaNo)%name
            endif
!!!!!!!!!!!!!!! KK
!            call epdebugTrack(abs(IoTrace),'', cTrack)
!!!!!!!!!!!!!!            
         endif
!            end point of  the segment
         call epl2w(Cn, Move%Track%pos, posw)
         if(TimeStruc) then
            write(abs(IoTrace),'(3g16.8,i4,g11.4, i4, 
     *         g16.8)')
     *           posw%x, posw%y, posw%z, k,
     *           erg,
     *           Move%Track%p%charge,
     *           Move%Track%t
         else
            write(abs(IoTrace),'(3g16.8,i4,g11.4, i4, 1x,a)')
     *           posw%x, posw%y, posw%z, k,
     *           erg, 
     *           Move%Track%p%charge, Media(MediaNo)%name
         endif
!!!!!!!!!!!!!!
!         call epdebugTrack(abs(IoTrace),'', Move%Track)
!!!!!!!!!!!!
         xxx = posw%x
         yyy = posw%y
         zzz = posw%z
         kkk = k
         chkchg = Move%Track%p%charge
!cc      endif
       end
      subroutine epdebugTrack(io, msg, aTrack)
      implicit none
#include "Zglobalc.h"
#include "ZepTrack.h"
      integer,intent(in):: io
      character(*),intent(in):: msg
      type(epTrack),intent(in):: aTrack

      real(8):: posw(3),  dir(3), pabs, pw(3)
      character(16)::struc, matter
      if( len(trim(msg)) > 0  ) then
         write(io,*) ' track info:---------------------------- ',
     *        trim(msg)
      endif
      call epqstruc(aTrack%cn, struc)
      call epqmat(aTrack%cn, matter)
      write(io,'(a, i7, a,a)') "comp n=", aTrack%cn, " ",trim(struc)
      write(io,'(a)')   "     ", trim(matter)
      write(io,'(a, 4i4)') "code sub chg=",
     *   aTrack%p%code, aTrack%p%subcode, aTrack%p%charge
      write(io,'(a, 1p, 4g15.5)') 'L 4mom =', aTrack%p%fm%p(1:4)            
      call epl2w(aTrack%cn, aTrack%pos, posw)
      call epl2wd(aTrack%cn, aTrack%w,  dir)

      pabs = sqrt(
     *     dot_product( aTrack%p%fm%p(1:3),  aTrack%p%fm%p(1:3)))
      pw(1:3) = pabs*dir(1:3)
      write(io,'(a, 1p, 4g15.5)') 'W 3mom =', pw(1:3)
      
      write(io, '("pos L", 1p, 3G15.5)')
     *      aTrack%pos%x,   aTrack%pos%y,   aTrack%pos%z
      write(io, '("pos W", 1p, 3G15.5)') posw(1:3)
      write(io, '("dir  L", 1p, 3G15.5)')
     *      aTrack%w%x,   aTrack%w%y,   aTrack%w%z
      write(io, '("dir W", 1p, 3G15.5)') dir(1:3)
      end
