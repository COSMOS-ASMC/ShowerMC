!
!     This is for such a case when you forgot the initinial random
!     number to reproduce primary diretion, but you have
!     the incident position in detector system.
!     param file  is read via file no. 10 and
!     incident position is read via stdin.
!  
!     get transformation matrix for converting detector system vector
!     into primary system vector.
!     Since primary system depends on each sampled primary, 
!     This program gets the conversion matrix for the condition
!     specified in the param and input incident pos.
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
      real*8 xin, yin, zin, temp
      type (track):: inci

      open(10, file="param")

      call creadParam(10)
!
!        xin, yin, zin,  incident point in Detector system.
!
      read(*,*) xin, yin, zin
      write(0,*) ' input inci pos='
      write(0,*)  xin, yin, zin
! 
      cc%x = xin
      cc%y = yin
      cc%z = zin

      call cbeginRun
      call cmkIncident(inci, fin)
      call ciniTracking(inci)
!          Above will make some incident but different from 
!          the stdin.
      det = ObsSites(NoOfSites)%pos%xyz

      write(0, *) ' detector base pos in Exyz'
      write(0,'(3g17.9)')  det%x, det%y, det%z
      call cdet2xyz(det, cc, inci%pos%xyz)
      write(0,*) 'incident position in Exyz'
      write(0,*) inci%pos%xyz%x, inci%pos%xyz%y,
     *     inci%pos%xyz%z
      call cxyz2det(det, inci%pos%xyz, cc)
      write(0,*) 'incident position in Det; must be same as input '
      write(0,*) cc%x, cc%y, cc%z

      call cvecProd(DetZaxis, DetXaxis, DetYaxis)
      write(0,*) '   '
      write(0,*) ' detector  X,Y,Z axis vecgor in Exyz'
      write(0,'(3g17.9)') DetXaxis%x, DetXaxis%y, DetXaxis%z
      write(0,'(3g17.9)') DetYaxis%x, DetYaxis%y, DetYaxis%z
      write(0,'(3g17.9)') DetZaxis%x, DetZaxis%y, DetZaxis%z
      write(0,*)
      temp = sqrt(xin**2 + yin**2 + zin**2)
!        Zprimary  in Det. system.
      cc%x = xin/temp
      cc%y = yin/temp
      cc%z = zin/temp
!          in Exyz
      call cdet2xyzD(cc, Zprimary)

      call cvecProd(Zprimary, DetZaxis, Xprimary)
!         see if Zprimary // DetZaxis; if so reset Xprimary
      temp= Xprimary%x**2 + Xprimary%y**2 + Xprimary%z**2
      if(temp .lt. 1.e-12) then
         Xprimary = DetXaxis
      else
         temp =sqrt(temp)
         Xprimary%r(1) = Xprimary%r(1)/temp
         Xprimary%r(2) = Xprimary%r(2)/temp
         Xprimary%r(3) = Xprimary%r(3)/temp
      endif
      call cvecProd(Zprimary, Xprimary, Yprimary)
!                                                                                       
!       we try to convert {X,Y,Z}primary into the one in
!           detector system-->a,b,c.
      call citransVectZx(1, DetZaxis, DetXaxis, Xprimary, aa)
      call citransVectZx(2, DetZaxis, DetXaxis, Yprimary, bb)

!      primary etc is the vector  representing the primary system axes
!     in the Detector system. Then, a given vector R in Detector system,
!     its component to the a,b,c is Rx, Ry, Rz in the Primary is primary
!     system.
!
      write(0,*) ' Primary system X,Y,Z axis unit vector in Exyz'
      write(0,'(3g17.9)')  Xprimary%x, Xprimary%y, Xprimary%z
      write(0,'(3g17.9)')  Yprimary%x, Yprimary%y, Yprimary%z
      write(0,'(3g17.9)')  Zprimary%x, Zprimary%y, Zprimary%z

      write(0,*) ' For a given vector in the Detector system'
      write(0,*)
     *   ' take scaler product of  the following unit vector to get'
      write(0,*) ' its x value in  the primary system. For y, z, use '
      write(0,*) ' 2nd, 3rd row'
      write(*,'(3g17.9)')  aa%x, aa%y, aa%z
      write(*,'(3g17.9)')  bb%x, bb%y, bb%z
      write(*,'(3g17.9)')  cc%x, cc%y, cc%z

      call citransVectZx(2, DetZaxis, DetXaxis, Zprimary, cc)
      write(0,*) ' this must be the same as above'
      write(0,'(3g17.9)')  cc%x, cc%y, cc%z
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


             
