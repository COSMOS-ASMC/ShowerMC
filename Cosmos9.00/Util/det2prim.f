!
!     get transformation matrix for converting detector system vector
!     into primary system vector.
!    Since primary system depends on each sampled primary, the 
!    random number by param file must be fixed.
!    This program gets the conversion matrix for the 1 st event
!    specified in the param
!

#include  "BlockData/cblkGene.h"
      implicit none

#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsv.h"




      type (coord):: aa, bb, cc, det, DetYaxis
      integer fin
      type (track):: inci

      open(10, file="param")

      call creadParam(10)
      call cbeginRun
      call cmkIncident(inci, fin)
      call ciniTracking(inci)

      det = ObsSites(NoOfSites)%pos%xyz

      write(0, *) ' detector base pos in Exyz'
      write(0,'(3g17.9)')  det%x, det%y, det%z
      write(0,*) 'incident position in Exyz'
      write(0,*) inci%pos%xyz%x, inci%pos%xyz%y,
     *     inci%pos%xyz%z
      call cxyz2det(det, inci%pos%xyz, cc)
      write(0,*) 'incident position in Det'
      write(0,*) cc%x, cc%y, cc%z

      call cvecProd(DetZaxis, DetXaxis, DetYaxis)
      write(0,*) '   '
      write(0,*) ' detector  X,Y,Z axis vecgor in Exyz'
      write(0,'(3g17.9)') DetXaxis%x, DetXaxis%y, DetXaxis%z
      write(0,'(3g17.9)') DetYaxis%x, DetYaxis%y, DetYaxis%z
      write(0,'(3g17.9)') DetZaxis%x, DetZaxis%y, DetZaxis%z
!       we try to convert {X,Y,Z}primary into the one in
!           detector system-->a,b,c.
      call citransVectZx(1, DetZaxis, DetXaxis, Xprimary, aa)
      call citransVectZx(2, DetZaxis, DetXaxis, Yprimary, bb)
      call citransVectZx(2, DetZaxis, DetXaxis, Zprimary, cc)
!      primary etc is the vector  representing the primary system axes
!     in the Detector system. Then, a given vector R in Detector system,
!     its component to the a,b,c is Rx, Ry, Rz in the Primary is primary
!     system.
!
      write(0,*) ' Primary systrem X,Y,Z axis vector in Exyz'
      write(0,'(3g17.9)')  Xprimary%x, Xprimary%y, Xprimary%z
      write(0,'(3g17.9)')  Yprimary%x, Yprimary%y, Yprimary%z
      write(0,'(3g17.9)')  Zprimary%x, Zprimary%y, Zprimary%z

      write(0,*) ' For a given vector in the Detector system'
      write(0,*) ' take scaler product of  the following vector to get'
      write(0,*) ' its x value in  the primary system. For y, z, use '
      write(0,*) ' 2nd, 3rd row'
      write(*,'(3g17.9)')  aa%x, aa%y, aa%z
      write(*,'(3g17.9)')  bb%x, bb%y, bb%z
      write(*,'(3g17.9)')  cc%x, cc%y, cc%z

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


             
