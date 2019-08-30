!     scintillation ligh is generated from energy deposit sotred
!     in a cell.
!         
      subroutine epLightScintiFromCell(cell)
      use modepLight
      use modepLightMaxDef
      use modepLightPty
      implicit none
#include "Zcode.h"
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "ZsepManager.h"


       type(epTrack)::  cell    ! energy deposit in a cell



      real(8):: wl, wl0    ! wave length 
      real(8):: photons ! <no. of photons from scinti>
      integer:: nphotons ! with poisson fluctuatin
      integer:: npreal   ! really sampled #
      real(8):: Nsmp   ! max of really sampled #
       type(epTrack)::  scintilight
      
      integer:: i
      real(8)::  u, x, y, z


      character(len=MAX_STRUCCHR)::struc, basename  ! v9.164 old *16
      integer::nattr
      real(8)::vol(maxVolAttr) 
      real(8)::weight
!/////////////////
!      integer onlyonce
!      common /dddebug/ onlyonce
!/////////////

      Nsmp =  comInfo(cPtyNo)%NpSample
!           clear p.e counter for this cell  ; this will be done
!           at tracking time by looking at subcode
!           10 <subcode <100==> pointer for p.e counter is made to be
!             tempPE and cleared
!           subcode< 10==> add p.e to the currnt counter
!           subcode>100==> current counter is multiplied by wgt and
!                          added non temporary p.e counter
!                          pointer for p.e counter is made to be
!                          non temporary p.e. counter
!  
!      Lcomp(cLcompNo)%tempPE = 0

      scintilight = cell  ! inherit some from cell
!             make light:       scinti's subcode is 1
      call cmkptc(klight, 1, 0, scintilight%p)


            !                        energy loss in in GeV
      photons = comInfo(cPtyNo)%NpPerMeV*
     *         cell%p%fm%p(4)*1.e3 ! MeV-->N photons
      if( photons < usePoisson ) then
         call kpoisn(photons, nphotons)
      else
         nphotons = photons + 0.5
      endif


      if(nphotons > Nsmp )  then
         weight = 1. +  real(nphotons-NSmp)/NSmp
         npreal = Nsmp
      else
         npreal = nphotons
         weight = 1. 
      endif

      sumni = sumni + npreal
      sumniwi = sumniwi + npreal*weight

      scintilight%wgt = weight
      scintilight%p%subcode = 1

      call epqstruc(Cn, struc)
      call epqvolatr(Cn, nattr, vol)

      call epGetBaseStrucName(struc, basename)  ! v9.164

      do i = 1, npreal
!               write generated uniformly in the cell  ;
!               cell size is contained in direction cos part
         call rndc(u)
         x = cell%pos%x + u*cell%w%x
         y = cell%pos%y + u*cell%w%y
         z = cell%pos%z + u*cell%w%z
         if( basename == "box" .or.
     *       basename == "octagon"        ) then
            if( x < EpsLength2 ) then
               x = EpsLength2
            elseif( x > vol(1)-EpsLength2 ) then
               x = vol(1) - EpsLength2
            endif
            if( y < EpsLength2 ) then
               y = EpsLength2
            elseif( y > vol(2)-EpsLength2 ) then
               y = vol(2) -EpsLength2
            endif
            if( z < EpsLength2 ) then
               z = EpsLength2
            elseif( z > vol(3)-EpsLength2 ) then
               z = vol(3) -EpsLength2
            endif
         elseif( basename == "cyl" ) then
            !  check if this position is out side of the cylinder
            if( x**2 + y**2 >= vol(1)**2 ) then
!                  force to make inside of the cyl
               x = cell%pos%x
               y = cell%pos%y
            endif
         elseif( basename == "ecyl" ) then
            if( (x/vol(1))**2  +  (y/vol(2))**2  >= 1. ) then
!                 force to make inside of the cyl
               x = cell%pos%x
               y = cell%pos%y
            endif
         elseif( basename == "pipe" ) then
            if( x**2 + y**2 < vol(1)**2 .or.
     *          x**2 + y**2 > vol(2)**2  ) then
               x = cell%pos%x
               y = cell%pos%y
            endif
         else
            write(0,*) struc, ' is not yet supported'
            write(0,*) ' for light ray tracing'
            stop
         endif
         scintilight%pos%x = x
         scintilight%pos%y = y
         scintilight%pos%z = z
!              get w.l
!/////////////
!         call Lcompchk("ZZ", cLcompNo)
!////////
         if( comInfo( cPtyNo )%waveAF > 0  ) then
!                distribution exists; sample wl
            call csampAF( comInfo( cPtyNo )%waveAF, wl)
         else
!                use always  fixed wl
            wl = comInfo(cPtyNo)%peakWL
         endif
!             put energy for same use; wl and  wl0 is the same
!             we use always same w.l in any media
         call epLightwl2E(wl, 1.d0, wl0, scintilight%p%fm%p(4))  
!              set wl
         scintilight%wl = wl

!              set isotropic angle
         call episoAngle(scintilight%w)
         call eppush(scintilight)
      enddo

      end
