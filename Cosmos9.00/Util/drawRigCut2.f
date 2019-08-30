#include "BlockData/cblkManager.h"
!
!   This is to make grid data and physical data (prob) for rigidity cut off
!   of new format
!
      implicit none
#include "Zmanagerp.h"
#include "ZrigCut.h"

      character*120 file
!

      file = '../Contrib/Data/CutoffNew/Kamioka'
      write(ErrorOut, *) ' enter file name=',file
      read(*, *) file

      call crigCut0(file)
      call crigCut2grid(RigTbl2, ZenSize, AzmSize, RigSize)
      call crigCut2data(RigTbl2, ZenSize, AzmSize, RigSize)
      write(ErrorOut, *)
     * 'Note:You have to set "variable aspect ratio" in Nextcontour'
      end
!     ***************************
      subroutine crigCut2data(tbl, izen, iphi, irig)
      implicit none
      integer izen, iphi, irig
      real*4 tbl(izen, iphi, irig), z
      integer i1, i2, i3,icon

      call copenfw(21, "/tmp/nxconData",icon)
      write(21, *) iphi, izen, irig
      do i3 =1, irig
         do i1 = izen, 1, -1 
            do i2 = 1, iphi
               z = tbl(i1, i2, i3)
               if(z .gt. 0.99) z =1.
               write(21,*) z
            enddo
         enddo
      enddo
      close(21)
      end
!     ************************
      subroutine crigCut2grid(tbl, izen, iphi, irig)
      implicit none
#include "ZrigCut.h"
      integer izen, iphi, irig
      real*4 tbl(izen, iphi, irig)
      integer i1, i2, i3,icon

      call copenfw(21, "/tmp/nxconGrid",icon)
      write(21, *) iphi, izen, irig
      do i3 =1, irig
         do i1 =  izen, 1, -1
            do i2 = 1, iphi
               write(21,*) sngl(MinAzm + (i2-1)*DAzm)
            enddo
         enddo
      enddo
      do i3 =1, irig
         do i1 = izen, 1, -1
            do i2 = 1, iphi
               write(21,*) sngl(ZenMax + (i1-1)*DZen), 
            enddo
         enddo
      enddo
      do i3 =1, irig
         do i1 = izen, 1, -1
            do i2 = 1, iphi
               write(21,*)
     *          sngl(MinRig * 10.**(LogDRig*(i3-1)))
            enddo
         enddo
      enddo

      close(21)
      end


      

