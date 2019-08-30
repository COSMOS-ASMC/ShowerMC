!      cinitObs:  set up observation sites by reference to Zobsp.h
!        HeightList is obtained from DepthList.
!        However, if DepthList value is negative, HeightList has priority and
!        DepthList is obtained from HeightList.
!        Variables in Zobsv.h are calculated.
      subroutine cinitObs
      use modAtmosDef
      implicit none
!
#include  "Zglobalc.h"
#include  "Ztrackp.h"
#include  "Zcoord.h"
#include  "Zpos.h"
#include  "Zmagfield.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
! #include  "Zearth.h"
#include  "Zincidentp.h"
!
       integer i, icon
       logical error
       character*180  msg
       type(coord)::south 
       real*8  cvh2thick, cthick2h, mu
       external cblkObs
!
       NoOfSites = 0
       do i = 1, maxNoOfSites
          if( DepthList(i) .eq. 0. ) goto 100    ! exit
          NoOfSites = NoOfSites + 1
          if( DepthList(i) .lt. 0. ) then
!                HeightList has priority
               DepthList(i) = cvh2thick(HeightList(i))
          else
               HeightList(i) = cthick2h(DepthList(i))
          endif
          call csetCoord('llh', LatitOfSite, 
     *    LongitOfSite, 
     *    HeightList(i),   ObsSites(i)%pos%xyz )
       enddo
 100   continue
       if(NoOfSites .eq. 0) then
          call cerrorMsg('No of Observation sites = 0', 0)
       endif
!        set polar angle of the injection pointsite.
!        (used when ObsPlane==spherical)
       call csetCoord('llh', LatitOfSite, LongitOfSite, 
     *    HeightOfInj,   PolarInjPos)
!        transform it to polar coord.
       call ctransCoord2('sph', PolarInjPos, PolarInjPos)

       if(EndLevel .eq. 0) EndLevel = NoOfSites   ! 96.09.18
       if(EndLevel2 .eq. 0) EndLevel2 = EndLevel
!               set magnetic field for the deepest level. 
!                 Magfield is 'ned' here
       if(HowGeomag .gt. 20 .and. HowGeomag .le. 30 ) then
          MagfieldNED%x = MagN
          MagfieldNED%y = MagE
          MagfieldNED%z = MagD
          MagfieldNED%sys = 'ned'
       else
!           if HowGeomag <=2, this is reset later. don't worry.
          call cgetMagF(YearOfGeomag, ObsSites(NoOfSites)%pos%xyz, 
     *              MagfieldNED, icon)
          if(icon .eq. 1) then
             write(msg, *) ' YearOfGeomag=',YearOfGeomag, ' is ok ?'
             call cerrorMsg(msg, 1)
          endif
       endif
!              to 'hva'
       call ctransMagTo('hva',    ObsSites(NoOfSites)%pos%xyz, 
     * MagfieldNED, MagfieldHVA)
       call ctransMagTo('xyz',    ObsSites(NoOfSites)%pos%xyz, 
     * MagfieldNED, MagfieldXYZ)
       if(abs(XaxisFromSouth) .gt. 360.d0) then
!           if it is > 360, use magnetic east to be the detector x-axis. 
!           angle form south to x-axis in counter closckwise.
!        Note:   'hva' system angle  + is clockwise.
!       XaxisFromSouth = 90.0 - MagfieldHVA.z
          XaxisFromSouth = 90.0 - MagfieldHVA%z
       endif

!          celestial coord. init.
       call kcelei(LatitOfSite, 
     *   LongitOfSite, DtGMT, 
     *   ObsSites(NoOfSites)%pos%height)
       call kadthi(XaxisFromSouth)

 
!          site coord is made  to be in the  xyz system.       
       do i = 1, NoOfSites
          call ctransCoord2('xyz', ObsSites(i)%pos%xyz,
     *        ObsSites(i)%pos%xyz )  ! to xyz
!                set ObsSites(i).pos.radiallen, ..depth, height
          call csetPos2(HeightList(i), DepthList(i), ObsSites(i)%pos)
          call cgetMoliereU(ObsSites(i)%pos%depth, 1.d0, mu)
          ObsSites(i)%mu = mu
       enddo
!          detector system coord. x is directed to mag. east.
!                                 z is directed to vertical.
!          get their direction cos. in 'xyz' system.
!
       CosLatitude = cos(LatitOfSite*Torad)
       SinLatitude = sin(LatitOfSite*Torad)
       CosLongitude = cos(LongitOfSite*Torad)
       SinLongitude = sin(LongitOfSite*Torad)
       DetZaxis%r(1) = CosLatitude * CosLongitude
       DetZaxis%r(2) = CosLatitude * SinLongitude
       DetZaxis%r(3) = SinLatitude
!          south direction is expressed as
       south%r(1) = SinLatitude *CosLongitude
       south%r(2) = SinLatitude *SinLongitude
       south%r(3) = -CosLatitude
!        in south-vertical system. X-axis is expressed
       DetXaxis%r(1) = cos(XaxisFromSouth * Torad)
       DetXaxis%r(2) = sin(XaxisFromSouth * Torad)
       DetXaxis%r(3) = 0.
!           convert it to 'xyz' system       
       call ctransVectZx(1, DetZaxis, south, DetXaxis, DetXaxis)
       Txyz2det(1,:) = DetXaxis%r(:)
       Txyz2det(3,:) = DetZaxis%r(:)
       call cvecprod(DetZaxis%r(:), DetXaxis%r(:), DetYaxis%r(:))
       Txyz2det(2,:) = DetYaxis%r(:)
       Tdet2xyz(:,:) = transpose(Txyz2det(:,:))
!        check if in the ascending order of depth
       error = .false.
       do i = 2, NoOfSites
          if(DepthList(i-1) .ge. DepthList(i)) then
             write(msg, *) 'Depth is not in ascending order ',
     *       i, '-th one=',DepthList(i), 
     *       ' is <  previous one= ',  DepthList(i-1)
             call cerrorMsg(msg, 1)
             error = .true.
          endif
       enddo
       if(error) then
          call cerrorMsg('Obsev. depth/height not in order', 0)
       endif
       if(ObsPlane .eq. spherical .and.
     *       (Za1ry .ne. 'is' .and. Za1ry .ne. 'cos' ) ) then
          call cerrorMsg(
     *      "ObsPlane==spherical but Za1ry !='is'/'cos'",0)
       endif
!       from v5.66
       if(BorderHeightL .eq. 0.) then
          BorderHeightL = ObsSites(NoOfSites)%pos%height -1.0
       endif
       ObsSites(0)%pos%height = BorderHeightH
       ObsSites(0)%pos%radiallen = BorderHeightH + Eradius
       ObsSites(0)%pos%depth = cvh2thick(BorderHeightH)
       ObsSites(NoOfSites+1)%pos%height = BorderHeightL
       ObsSites(NoOfSites+1)%pos%radiallen = BorderHeightL + Eradius
       ObsSites(NoOfSites+1)%pos%depth = cvh2thick(BorderHeightL)
       
!         +++++++++++++++++++++++ for AS ++++++++++++++++
       NoOfASSites = 0
       do i = 1, maxNoOfASSites
          if( ASDepthList(i) .eq. 0. ) goto 200    ! exit
          NoOfASSites = NoOfASSites + 1
          if( ASDepthList(i) .lt. 0. ) then
!                HeightList has priority
               ASDepthList(i) = cvh2thick(ASHeightList(i))
          else
               ASHeightList(i) = cthick2h(ASDepthList(i))
          endif
          call csetCoord('llh', LatitOfSite, 
     *    LongitOfSite, 
     *    ASHeightList(i),   ASObsSites(i)%pos%xyz )
       enddo
 200   continue
       
!          site coord is made  to be in the  xyz system.       
       do i = 1, NoOfASSites
          call ctransCoord2('xyz', ASObsSites(i)%pos%xyz,
     *        ASObsSites(i)%pos%xyz )  ! to xyz
!                set ASObsSites(i).pos.radiallen, ..depth, height
          call csetPos2(ASHeightList(i), ASDepthList(i), 
     *     ASObsSites(i)%pos)
!             Moliere Unit should be reset depending on
!             the 1ry angle. Here is a tentative one for vertical           
          call cgetMoliereU(ASObsSites(i)%pos%depth, 1.d0, mu)
          ASObsSites(i)%mu = mu
       enddo

       end
