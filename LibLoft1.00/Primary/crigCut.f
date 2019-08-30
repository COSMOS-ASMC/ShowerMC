!     ***********************************************************
!        To know the arrival prob. as a function of direction at a given
!        place.  
!
      subroutine crigCut(azmin, zen, rig, prob)
      implicit none
#include "ZrigCut.h"
!          Before calling this subroutine, you must have called
!       crigCut0 to get the table which contains the information
!       of cutoff rigidities.

!
      real*8 azmin  ! input. value of azimuth angle at a given location.
                   !         deg.
                   ! For format# 5, this is longitude of a place
      real*8 zen   ! input. value of zenith angle at a given location. 
                   !     It may be deg, or cos, depending on the table.
                   ! For format# 5, this is latitude of the place.
                   !     It may be deg, or cos, depending on the table.
      real*8 rig   ! input. value of particle rigidity (GV) p/Z

      real*8 prob  ! output. probability that 'r' can enter the earth.
                   !  For format# 5, this is  0 or 1. 1 dose not necessary
                   !   mean that the particle can enter the earth;
                   !   there is possibility so you have to examine.
                   !   if 0, there is no possibility so you may discard it.
!
!      Note: This system cannot treat multiple locations at a time.
!
!
      real*8 rigcut

      real*8  azm

      integer i1, i2, i3

      azm = azmin
      if(azm .lt. 0.) then
         azm = azm + 360.
      else
         azm = mod(azm, 360.d0)
      endif

      if(Rdatafmt .eq. 1 ) then
         call k4ptdi(RigCutTbl, AzmSize, ZenSize, AzmSize, MinAzm, 
     *        MinZen, DAzm,  DZen,  azm, zen, rigcut)
         if(rig .lt. rigcut) then
            prob = 0.
         else
            prob = 1.
         endif
      elseif(Rdatafmt .le. 4) then
!         i3 = log10(rig/MinRig)/LogDRig + 1
!         i3 <= log10(rig)*20 < i3+1 (i3=0,1..)
!         logrig = log10(rig)/LogDRig
!         i3 = logrig + 1 
!         dlogrig=(logrig - (i3-1))*LogDRig 

         i3 = log10(rig/MinRig)/LogDRig + 1
         if(i3 .lt. 1) then
            prob = 0.
            goto 100
         endif
         if(i3 .ge. RigSize) then
            prob = 1.
            goto 100
         endif
         if(Rdatafmt .eq. 4) then
            i2=1
         else
            i2 =  azm/DAzm + 1
            i2 =min(AzmSize, i2)
         endif
         i1 = (zen-ZenMax)/DZen +1     ! DZen < 0
         i1 = min(ZenSize, i1)
         call cgetRigProb(
     *    RigTbl2, ZenSize, AzmSize, RigSize, i1, i2, i3,
     *    prob ) 
         prob = min(prob, 1.d0)
      elseif( Rdatafmt .eq. 5) then
!           essentiall the same as 1 case. but rather should not
!           make interpolation
         call k2dtblv(RigCutTbl, AzmSize, ZenSize, AzmSize, MinAzm, 
     *        MinZen, DAzm,  DZen,  azm, zen, rigcut)
         if(rig .lt. rigcut) then
            prob = 0.
         else
            prob = 1.
         endif
      else
         call cerrorMsg('format # error in rigidity table',0)
      endif
 100  continue
      end
      subroutine k2dtblv(tbl, xs, ys, adj, xm, 
     *        ym, dx,  dy,  x, y,  ans)
      integer xs
      integer ys
      integer adj
      real*8 xm
      real*8 ym
      real*8 dx
      real*8 dy
      real*8 x, y
      real*8 ans
      real*8 tbl(adj, ys)
      
      integer i, j
      
      i = (x-xm)/dx + 1
      j = (y-ym)/dy + 1
      ans = tbl(i,j)
      end


      subroutine cgetrigCut(azmin, zen, rigcut)
      implicit none
#include "ZrigCut.h"
!         get rigidity cut (for old table only); use crigCut now
!
      real*8 azmin  ! input. value of azimuth angle at a given location.
                   !         deg.
      real*8 zen   ! input. value of zenith angle at a given location. 
      real*8 rigcut ! output. cutoff rigidity. in GV
!

      real*8 azm

      azm = azmin

      if(azm .lt. 0.) then
         azm = azm + 360.
      endif
      azm = mod(azmin, 360.d0)
      if(Rdatafmt .eq. 1) then
         call k4ptdi(RigCutTbl, AzmSize, ZenSize, AzmSize, MinAzm, 
     *        MinZen, DAzm,  DZen,  azm, zen, rigcut)
      else
         call cerrorMsg(
     *       'only old cutoff table can be used', 0)
      endif
      end
      subroutine cgetRigProb(tbl, 
     * izen, iphi, irig, i1, i2, i3, prob)
      implicit none
#include "ZrigCut.h"

      integer izen, iphi, irig
      real tbl(izen, iphi, irig)
      integer i1, i2, i3
      real*8 prob

!      if( i3 .lt. RigSize ) then
!       prob = dlogrig/LogDRig * 
!     *          (tbl(i1,i2,i3+1)-tbl(i1, i2, i3))
!     *         +  tbl(i1,i2,i3)
!      else
!         prob = tbl(i1,i2,i3)
!      endif

      prob = tbl(i1,i2,i3)

      end
!-------------------------        init for crigCut
      subroutine crigCut0(file)
      implicit none
#include "Zmanagerp.h"
#include "ZrigCut.h"

      character*(*) file
      character*192 msg
      integer icon

      call copenf(TempDev, file, icon)
      if(icon .ne. 0) then
         write(msg, *) ' file specification error '
         call cerrorMsg(msg, 0)
      endif
!
      read(TempDev, '(1x,i1)') Rdatafmt   ! format 1 or 2 or 3 or 4
      if(Rdatafmt .eq. 1) then
         call crigCutfmt1
      elseif(Rdatafmt .eq. 2 .or. Rdatafmt .eq. 3) then
         call crigCutfmt2
      elseif(Rdatafmt .eq. 4) then
         call crigCutfmt4
      elseif(Rdatafmt .eq. 5) then
         call crigCutfmt5
      else
         write(msg, *) 'rigidity cut data format =',Rdatafmt,
     *   ' ivalid' 
         call cerrorMsg(msg, 0)
      endif
      end
!     **********************
      subroutine crigCutfmt1
      implicit none
#include "Zmanagerp.h"
#include "ZrigCut.h"


      character*192 msg

      integer icon

      read(TempDev, *)   ! skip
      read(TempDev, *) Place,  Latit,  Longi, MagDec,  AzmValue,
     *            DAzm, AzmSize, ZenValue, DZen, ZenSize

!           skip until -------------- line
      call cskipComment(TempDev, icon)
      if(icon .ne. 0) then
          call cerrorMsg('cut-off table data has no ---- line', 0)
      endif
!          see if azimuthal angle range is from 0 to 360 or from 0 to 360-DAzm

      if(AzmValue .eq. 'deg') then
!               assume min of Azimuthal angle is  0 deg.
          MinAzm = 0.
          if( (AzmSize - 1)* DAzm  .lt. (360.d0- DAzm* 0.1d0) ) then
!                 fill the last col.  by the first col.

             call creadRigCut(TempDev, RigCutTbl, AzmSize, ZenSize, 
     *       AzmSize+1)
             AzmSize = AzmSize + 1
             call cfillRigCut(RigCutTbl, AzmSize, ZenSize)
          else
             call creadRigCut(TempDev, RigCutTbl, AzmSize, ZenSize,
     *        AzmSize)
          endif

       else

          write(msg, *)
     *    ' Azimuthal angle unit must be deg for rigidity cut table'
          call cerrorMsg(msg, 1)
          write(msg, *) ' But it is ', AzmValue
          call cerrorMsg(msg, 0)
       endif

       close(TempDev)

       if(ZenValue .eq. 'cos') then
!              DZen should be negative and min is 1.0
          MinZen = 1.0
          if(DZen .ge. 0.) then
             write(msg, *)
     *        ' step of Zenith angle for rigidity cut should be < 0'
             call cerrorMsg(msg, 1)
             write(msg, *) ' because you give it in cos value'
             call cerrorMsg(msg, 0)
          endif
       endif
!
!                                   
       write(msg,*) 'Rigidity cut-off table has been read:',
     * ' place=',Place,' latitute=',Latit, ' longitude=',Longi,
     * ' mag. dec=', MagDec
       call cerrorMsg(msg, 1)
       end
!      ****************** 
       subroutine creadRigCut(io, tbl,  azm, zen, adj)
       implicit none
       integer io,  azm, zen, adj
       real*8 tbl(adj,  zen)


       integer i, j, ios
       character*150 msg

       do   j=1,  zen
           read(io, *, iostat=ios) ( tbl(i, j), i = 1, azm)
           if(ios .ne. 0) then
              write(msg, *) 'Unexpected EOF at rigicity table reading'
              call cerrorMsg(msg, 1)
              write(msg, *) ' line number=', j, ' azm=',azm, ' zen=',
     *          zen, ' adj=',adj
              call cerrorMsg(msg, 0)
           endif
       enddo
       end
!      *******************
       subroutine cfillRigCut( tbl,  azm, zen)
       implicit none
       integer  azm, zen
       real*8 tbl(azm,  zen)
       integer i

        do i = 1, zen
            tbl(azm, i) = tbl(1, i)
        enddo
       end
!      ********************************** 
       subroutine cqRigCutPlace(tlt, tlg, mdec)
       implicit none
#include "ZrigCut.h"
       real*8 tlt, tlg, mdec


!            to inform lat, long, magdec
        tlt = Latit
        tlg = Longi
        mdec = MagDec
       end

!     **********************
      subroutine crigCutfmt2
      implicit none
#include "Zmanagerp.h"
#include "ZrigCut.h"

      character*192 msg

      integer icon


      read(TempDev, *)   ! skip
      read(TempDev, *) 
     *  Place,  Latit,  Longi, MagDec, ZenValue, ZenMax, DZen, ZenSize,
     *  AzmValue, MinAzm, DAzm, AzmSize, MinRig, LogDRig, RigSize


!           skip until -------------- line
      call cskipComment(TempDev, icon)
      if(icon .ne. 0) then
          call cerrorMsg('cut-off table data has no ---- line', 0)
      endif
      call crdRigCut2(
     *       Rdatafmt, TempDev, RigTbl2, ZenSize, AzmSize, RigSize)
      close(TempDev)

      if(ZenValue .eq. 'cos') then
!           DZen should be negative
         if(DZen .ge. 0.) then
            write(msg, *)
     *        ' step of Zenith angle for rigidity cut should be < 0'
            call cerrorMsg(msg, 1)
            write(msg, *) ' because you give it in cos value'
            call cerrorMsg(msg, 0)
         endif
      endif
!
!                                   
       write(msg,*) 'New rigidity cut-off table has been read:',
     * ' place=',Place,' latitute=',Latit, ' longitude=',Longi,
     * ' mag. dec=', MagDec
       call cerrorMsg(msg, 1)
      end
!     **********************
      subroutine crigCutfmt4
      implicit none
#include "Zmanagerp.h"
#include "ZrigCut.h"

      character*192 msg

      integer icon

!      real pw

      read(TempDev, *)   ! skip
      read(TempDev, *) 
     *  Place,  Latit,  Longi, MagDec, ZenValue, ZenMax, DZen, ZenSize,
     *  MinRig, LogDRig, RigSize

      AzmSize=1


!           skip until -------------- line
      call cskipComment(TempDev, icon)
      if(icon .ne. 0) then
          call cerrorMsg('cut-off table data has no ---- line', 0)
      endif
      call crdRigCut4(
     *       Rdatafmt, TempDev, MinRig, RigTbl2, ZenSize, RigSize)
      close(TempDev)

      if(ZenValue .eq. 'cos') then
!           DZen should be negative
         if(DZen .ge. 0.) then
            write(msg, *)
     *        ' step of Zenith angle for rigidity cut should be < 0'
            call cerrorMsg(msg, 1)
            write(msg, *) ' because you give it in cos value'
            call cerrorMsg(msg, 0)
         endif
      endif
!
       write(msg,*) 'New rigidity cut-off table (fmt4) has been read:',
     * ' place=',Place,' latitute=',Latit, ' longitude=',Longi,
     * ' mag. dec=', MagDec
       call cerrorMsg(msg, 1)
      end
!    ------------------------
      subroutine crdRigCut2(fmt, io, tbl, izen,  iphi, irig)
!       read cut off table for format 2
      implicit none
      integer fmt ! input.  fmt=2 or 3, if  3, data for izen, iphi, irig are
                  !          not given
      integer io  ! input. table dev. number
      integer izen ! input. table for zenith angle
      integer iphi ! input. table for azimuthal angle
      integer irig ! input. table for rigidity
      real*4  tbl(izen, iphi, irig)  ! output. 3D table which shows prob. that
                          ! given ptcl with izen, iphi, irig can enter the earth
                  !  NOTE:  this  is  real*4--------------------

      integer i1, i2, i3
      integer j1, j2, j3
      character*128 msg


      do i1= 1 , izen
         do i2 = 1, iphi
            do i3 = 1, irig
               if(fmt .eq. 2) then
                  read(io, *, end=100) j1, j2, j3, tbl(i1, i2, i3)
                  if(i1 .ne. j1+1 .or.
     *               i2 .ne. j2+1 .or.
     *                       i3 .ne. j3+1) then 
                     write(msg, *)
     *               ' data index mismatch in new rigidit cut table',
     *                 i1, j1+1, i2, j2+1, i2, j3+1      
                     call cerrorMsg(msg, 0)
                  endif
               else
                  read(io, *) tbl(i1, i2, i3)
               endif
            enddo
         enddo
      enddo
      return
 100  continue
      call cerrorMsg(
     *   'Unexpected EOF in new rigidity cut table',0)
      end
!    ------------------------
      subroutine crdRigCut4(fmt, io, minval, tbl, izen, irig)
!       read cut off table for format 4
      implicit none
      integer fmt ! input.  fmt=4
      integer io  ! input. table dev. number
      real*8  minval ! input minimum value of rigidty (for check)
      integer izen ! input. table for zenith angle
      integer irig ! input. table for rigidity
      real*4  tbl(izen, irig)  ! output. 2D table which shows prob. that
                          ! given ptcl with izen, irig can enter the earth
                  !  NOTE:  this  is  real*4--------------------

      integer i1, i3
      integer j1,  idummy
      real*4 dummy
      character*128 msg


      do i1= 1 , izen
         do i3 = 1, irig
            read(io, *, end=100)
     *         j1, idummy, idummy, dummy, tbl(i1, i3)
            if(i3 .eq. 1) then
               if( abs(dummy - minval)/minval .gt. 1.e-3) then
                  write(msg,*)
     *             'check min rigidity in headr=', minval,
     *             ' in table=', dummy
                  call cerrorMsg(msg, 0)
               endif
            endif
            if(i1 .ne. j1) then
               write(msg, *)
     *           ' data index mismatch in new rigidit cut table',
     *              i1, j1
               call cerrorMsg(msg, 0)
            endif
         enddo
      enddo
      return
 100  continue
      call cerrorMsg(
     *   'Unexpected EOF in new rigidity cut table',0)
      end


!     **********************
      subroutine crigCutfmt5
      implicit none
#include "Zmanagerp.h"
#include "ZrigCut.h"

      character*192 msg

      integer icon


      read(TempDev, *)   ! skip
      read(TempDev, *)
     *  ZenValue, ZenMax, DZen, ZenSize,
     *  AzmValue, MinAzm, DAzm, AzmSize


!           skip until -------------- line
      call cskipComment(TempDev, icon)
      if(icon .ne. 0) then
          call cerrorMsg('cut-off table data has no ---- line', 0)
      endif

      if(AzmValue .eq. 'deg') then
!               assume min of logitudinal angle is  0 deg.
          MinAzm = 0.
          if( (AzmSize - 1)* DAzm  .lt. (360.d0- DAzm* 0.1d0) ) then
!                 fill the last col.  by the first col.
             call creadRigCut(TempDev, RigCutTbl, AzmSize, ZenSize, 
     *       AzmSize+1)
             AzmSize = AzmSize + 1
             call cfillRigCut(RigCutTbl, AzmSize, ZenSize)
          else
             call creadRigCut(TempDev, RigCutTbl, AzmSize, ZenSize,
     *        AzmSize)
          endif
          close(TempDev)
       else
          write(msg, *)
     *    ' Azimuthal angle unit must be deg for rigidity cut table'
          call cerrorMsg(msg, 1)
          write(msg, *) ' But it is ', AzmValue
          call cerrorMsg(msg, 0)
       endif


       if(ZenValue .eq. 'cos') then
!              DZen should be negative and min is 1.0
          MinZen = 1.0
          if(DZen .ge. 0.) then
             write(msg, *)
     *        ' step of Zenith angle for rigidity cut should be < 0'
             call cerrorMsg(msg, 1)
             write(msg, *) ' because you give it in cos value'
             call cerrorMsg(msg, 0)
          endif
       endif
!
!                                   
       write(msg,*) 'Rough rigidity cut-off table has been read:'
       call cerrorMsg(msg, 1)
      end

