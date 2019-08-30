!      convert a given (x,y,z) in the Exyz system into
!     the detector  system.
!
!   input.  latitude, longitude, height(depth)  from param.
!           (name 'param' is fixed)
!   (x,y,z) is from stdin.
!   output.  (X,Y,Z) in det system  on stdout.
!

#include  "BlockData/cblkGene.h"
      implicit none

#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsv.h"




      type (coord):: a, b, det, DetYaxis


      open(10, file="param")

      call creadParam(10)
      call cbeginRun
      det = ObsSites(NoOfSites)%pos%xyz

      call cvecProd(DetZaxis, DetXaxis, DetYaxis)
      write(*,*) ' transformation matrix '
      write(*,'(4g17.9)')  DetXaxis%x, DetYaxis%x, DetZaxis%x, 0
      write(*,'(4g17.9)')  DetXaxis%y, DetYaxis%y, DetZaxis%y, 0
      write(*,'(4g17.9)')  DetXaxis%z, DetYaxis%z, DetZaxis%z, 0
      write(*,'(4g17.9)')  det%x, det%y, det%z, 1 

      do while(.true.)
         write(*,*) 'enter x,y,z in Exyz sys'
         read(*,*)  a%x, a%y, a%z
         call cxyz2det(det, a, b)
         write(*,*) ' input=', a%x, a%y, a%z
         write(*,*) ' output=', b%x, b%y, b%z
      enddo
      end
      subroutine chookTrace
      end
      subroutine chookCeren
      end
      subroutine chookCerenS
      end
      subroutine chookCerenE
      end
      subroutine chookBgRun
      end
!
       temp%r(1) = dir1%r(1) * xax%r(1) + dir1%r(2) *xax%r(2) +
     *  dir1%r(3)  *xax%r(3)
       temp%r(2) = dir1%r(1) * yvec%r(1)+ dir1%r(2) *yvec%r(2) + 
     *  dir1%r(3) *yvec%r(3)
       temp%r(3) = dir1%r(1) * zax%r(1) + dir1%r(2) *zax%r(2) +
     *    dir1%r(3)  *zax%r(3)
       dir2 = temp


             
