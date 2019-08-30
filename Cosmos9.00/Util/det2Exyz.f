!      convert a given (x,y,z) in the detector system into
!     E-xyz system.
!
!   input.  latitude, longitude, height(depth)  from param.
!           (name 'param' is fixed)
!   (x,y,z) is from stdin.
!   output.  (X,Y,Z) in E-xyz on stdout.
!

#include  "BlockData/cblkGene.h"
      implicit none

#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"




      type (coord):: a, b, det, DetYaxis


!      open(10, file="param")

      call creadParam(6)
      call cbeginRun
      det = ObsSites(NoOfSites)%pos%xyz

      call cvecProd(DetZaxis, DetXaxis, DetYaxis)
!      write(*,*) ' transformation matrix '
!              (x,y,z,1) = (X,Y,Z, 1) * Matrix
      write(* ,
     * '("# place: long=",f8.2, " lat=",f7.2, " height=",f10.1)')
     *  LongitOfSite,  LatitOfSite,  ObsSites(NoOfSites)%pos%height
      write(*,'("#  E-xyz values are",3G17.9)') det%x, det%y, det%z
      write(*,'(3g17.9, a)')  DetXaxis%x, DetXaxis%y, DetXaxis%z, "  0"
      write(*,'(3g17.9, a)')  DetYaxis%x, DetYaxis%y, DetYaxis%z, "  0"
      write(*,'(3g17.9, a)')  DetZaxis%x, DetZaxis%y, DetZaxis%z, "  0"
      write(*,'(3g17.9, a)')  det%x, det%y, det%z, "  1"

!               this is for colum vector
!      write(*,'(3g17.9, a)')  DetXaxis.x, DetYaxis.x, DetZaxis.x, "  0"
!      write(*,'(3g17.9, a)')  DetXaxis.y, DetYaxis.y, DetZaxis.y, "  0"
!      write(*,'(3g17.9, a)')  DetXaxis.z, DetYaxis.z, DetZaxis.z, "  0"
!      write(*,'(3g17.9, a)')  det.x, det.y, det.z, "  1"


!      do while(.true.)
!         write(*,*) 'enter x,y,z in det sys'
!         read(*,*)  a.x, a.y, a.z
!         call cdet2xyz(det, a, b)
!         write(*,*) ' input=', a.x, a.y, a.z
!         write(*,*) ' output=', b.x, b.y, b.z
!      enddo
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
