       subroutine cifXObsSite(nextwhere)
!         see if MovedTrack crosses an observation depth.
!         or go out of boundry.
!     if to be observed, reset MovedTrack information
      use modSetIntInf
      use modAtmosDef
      implicit none
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zincidentv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
#include  "Zcode.h"
!     
      integer nextwhere  ! output next where. to be used after ptcl
                         !  crosses an  observation  level. 
!  for perpendicular case
!        cross       border
!   1      O       BorderL< Moved.Height <BorderH
!   2      X            //
!   3      O       M.H > B.H 
!   4      O       M.H < B.L
!   5      X       M.H > B.H
!   6      X       M.H < B.L
!
      integer loc
      real*8  leng, dedt, cosfromaxis
      real*8  clen2thick

      logical cross, seeupper

      type(coord)::xyz1
      type(coord)::xyz2
      type(coord)::dircos
      type(coord)::diffvec
      type(coord)::unitv
      type(coord)::unitp

      type(coord)::temp

      real*8 dummylen, r0, rp, obsrl,  cosp, dt2, m2
      real*8 coscr2, r0r
      real*8 dedtF
!
      cross = .false.

      loc = TrackBefMove%where

      if(IntInfArray(ProcessNo)%length .gt. 0. ) then
!           later use
         call cdiffvec(TrackBefMove%pos%xyz, 
     *              MovedTrack%pos%xyz, diffvec)
!           make it  to the unit vector    vector:  Bef  1--->2  Moved
         call c3DV2DDCos(diffvec, unitv, dummylen)
         if( ObsPlane .ne. notUsed ) then
            if( abs(ObsPlane) .eq. Perpendicular ) then
               call cifXPerpen(unitv,  cross, leng, seeupper)
            else
               call cifXHorizon(diffvec, unitv, cross, leng, seeupper)
            endif
         endif
      endif
         
      if(cross) then
         if(MovedTrack%p%charge /= 0) then
!            reset EnergyLoss to shorter path length
            EnergyLoss = EnergyLoss * leng/IntInfArray(ProcessNo)%length
            MovedTrack%p%fm%p(4) = TrackBefMove%p%fm%p(4)-EnergyLoss
         endif
         IntInfArray(ProcessNo)%length = leng
         IntInfArray(ProcessNo)%thickness = clen2thick(
     *     TrackBefMove%pos%radiallen - Eradius,
     *     TrackBefMove%vec%coszenith, leng)

         call cmoveStreight(leng, unitv)
         MovedTrack%vec%w = unitv 

!                      bef. replace, w for low energy electron
!                      may be largely scattered and even
!                      backward at the path end.
!                      At Xing point with obs. depth, this
!                      is not good.
!         if(MovedTrack.pos.height .le. BorderHeightL) then
!            write(0,*) MovedTrack.pos.height, BorderHeightL
!         MoveStat = BorderL
!            MovedTrack.where = NoOfSites+1


         if( MoveStat .ne. BorderH .and. MoveStat .ne. BorderL) then
            MoveStat = ToBeObserved
            if(seeupper) then
               nextwhere = loc - 1
               MovedTrack%where = nextwhere
            else
               MovedTrack%where = loc
               nextwhere = loc + 1
            endif
         endif

         if(MovedTrack%p%charge /= 0) then 
            call cresetMombyEdir(MovedTrack) 
!              reset Energyloss since E in MovedTrack might have been
!              reset to mass if energyloss was too big.
            EnergyLoss = TrackBefMove%p%fm%p(4)- MovedTrack%p%fm%p(4)
             ! this may not be used 
         else
!             charge = 0
           ! nothing to   do. even path is cut, momenutum etc
!            dose not change.
         endif
      else
         nextwhere = loc
         MovedTrack%where = loc
      endif
      end

      subroutine cifXHorizon(diffvec, unitv, cross, leng, seeupper)
      use modAtmosDef
      implicit none
#include  "Ztrack.h"
! #include  "Zmagfield.h"            
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zincidentv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
#include  "Zcode.h"
!     

      integer loc
      real*8  leng
      real*8  clen2thick

      logical cross, seeupper

      type(coord)::diffvec
      type(coord)::unitv
      type(coord)::unitp


      real*8 dummylen, r0, rp, obsrl,  cosp, dt2, m2
      real*8 coscr2, r0r
!

      loc =  TrackBefMove%where

      r0 = ObsSites(loc)%pos%radiallen
      rp = TrackBefMove%pos%radiallen
      dt2 =   rp**2 -  r0**2
      if(dt2 .lt. 0.d0) then
         if(dt2 .lt. -1.d0) then
            write(0,*) ' loc =',loc, ' r0=',r0, ' rp=',rp,
     *           ' dt2=',dt2, ' code=',TrackBefMove%p%code,
     *           ' rp%h=', TrackBefMove%pos%height,
     *           ObsSites(loc)%pos%height
            call cerrorMsg('position info invalid',0)
         endif
         dt2 = abs(dt2)
      endif
!         this is cos^2 of zenith for a ptcl which touches the inner sphere (loc-th sphere)
      coscr2 = dt2 /rp**2

      call c3DV2DDCos(TrackBefMove%pos%xyz, unitp, dummylen)
      call cscalerProd(unitv, unitp, cosp)
      cosp = -cosp
      seeupper = cosp .le. 0
      if(cosp .gt. 0) then
!             down going but not cross the lower sphere
         seeupper = cosp**2 .le. coscr2
      endif


      if(.not. seeupper ) then
         if(MovedTrack%pos%radiallen .le. r0) then
!                 cross ; goes down
            cross = .true.
         else
!                get moved length  ^2           
            call cscalerProd(diffvec, diffvec, m2)
            if(m2 .gt. dt2) then
!              double cross; first one is the nearest cross
               cross =.true.
            endif
         endif
      else
!                may see upper
         if(MovedTrack%pos%radiallen .ge.
     *        ObsSites(loc-1)%pos%radiallen) then
!                 if loc=1, this means going outside of the borderH
               cross =.true.
         endif
      endif

      if(cross) then
!             get length, leng, from the starting point to
!             the observation depth
         call cscalerProd(unitv, TrackBefMove%pos%xyz, r0r)
         if(seeupper) then
            obsrl = ObsSites(loc-1 )%pos%radiallen
            if(loc .eq. 1)  MoveStat=BorderH
         else
            obsrl = ObsSites(loc)%pos%radiallen
            if(loc .eq. NoOfSites+1 )then
!                this seems not needed.
               MoveStat = BorderL
            endif
         endif
!            get distance to the crossing point
         if(rp .gt. obsrl ) then
            leng = -r0r - sqrt(r0r**2 -
     *           (rp**2-obsrl**2))
         else
            leng = -r0r + sqrt(r0r**2 -
     *           (rp**2 - obsrl**2))
         endif
         leng = leng + 0.001d0
      endif
      end
      subroutine cifXPerpen(unitv,  cross, leng, seeupper)
      use modAtmosDef
      implicit none
#include  "Ztrack.h"
! #include  "Zmagfield.h"            
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zincidentv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
#include  "Zcode.h"
!     


      real*8  leng
      logical cross, seeupper

      real*8  clen2thick,  cosfromaxis, leng2

      integer loc, i, icon
      type(coord)::xyz1
      type(coord)::xyz2
      type(coord)::dircos
      type(coord)::unitv

      loc = TrackBefMove%where
!///////////
!      write(0,*) ' loc =', loc
!//////////
!
!           observation plane is perpendicular to primary
!           convert coord into 1ry system
      call cxyz2prim(ObsSites(NoOfSites)%pos%xyz, 
     *     MovedTrack%pos%xyz, xyz2)
      call cxyz2prim(ObsSites(NoOfSites)%pos%xyz, 
     *     TrackBefMove%pos%xyz, xyz1)
!      call cscalerProd(TrackBefMove.vec.w, DcAtObsXyz, 
!     *     cosfromaxis)
      call cscalerProd(unitv, DcAtObsXyz, 
     *     cosfromaxis)

      seeupper = cosfromaxis .le. 0.
      if(.not. seeupper) then
!              current ptcl is down going
         if( loc .gt. NoOfSites ) then
            if( MovedTrack%pos%height .le. BorderHeightL ) then
!                6
               MoveStat = BorderL
               MovedTrack%where = NoOfSites+1
               cross = .true.
            elseif( MovedTrack%pos%height .ge. BorderHeightH ) then
!                5
               MoveStat = BorderH
               MovedTrack%where = 0
               cross = .true.
            endif
         else
            if(xyz1%r(3) .gt. ObsSites(loc)%zpl .and.
     *           xyz2%r(3) .le. ObsSites(loc)%zpl) then
               cross = .true.
            elseif( MovedTrack%pos%height  .le. BorderHeightL )  then
!                escape from side  6
               cross = .true.
               MoveStat = BorderL
               MovedTrack%where = NoOfSites + 1
            elseif( MovedTrack%pos%height  .ge. BorderHeightH )  then
!               escape from side 5
               cross = .true.
               MoveStat = BorderH
               MovedTrack%where = 0
            endif
         endif
      else
!            upgooing
         if(loc .eq. 1)  then
            if( MovedTrack%pos%height  .ge. BorderHeightH )  then
!                   5
               MoveStat = BorderH
               MovedTrack%where = 0
               cross = .true.
            elseif(MovedTrack%pos%height  .le. BorderHeightL )  then
!                  6
               MoveStat = BorderL
               MovedTrack%where = NoOfSites+1
               cross = .true.
            endif
         else
            if(xyz1%r(3) .lt. ObsSites(loc-1)%zpl .and.
     *           xyz2%r(3) .ge. ObsSites(loc-1)%zpl ) then
               cross = .true.
            elseif( MovedTrack%pos%height  .ge. BorderHeightH )  then
!                 escaep from side  5
               cross = .true.
               MoveStat = BorderH
               MovedTrack%where = 0
            elseif( MovedTrack%pos%height  .le. BorderHeightL )  then
!                 escape from side 6
               cross = .true.
               MoveStat = BorderL
               MovedTrack%where = NoOfSites + 1
            endif
         endif
      endif
!/////////////
!      write(0,*) ' see upper=',seeupper, ' cross=',cross
!      write(0,*) ' where=', MovedTrack.where
!      write(0,*) ' MoveStat=',MoveStat
!      write(0,*) ' BorderL=',BorderL, 'BorderH=',BorderH
!/////////////
      if(.not. cross ) then
!     *      .and. MovedTrack.pos.height .gt. BorderHeightL  
!     *     .and.    MovedTrack.pos.height .lt. BorderHeightH ) then
!         most common  case; remain in the same layer; nothing to do 
      elseif( cross .and. MoveStat .eq. BorderL ) then
!         get length to the  x-ing point
         call cxplsph(
     *    TrackBefMove%pos%xyz%r(1), TrackBefMove%pos%xyz%r(2), 
     *    TrackBefMove%pos%xyz%r(3),
     *    unitv%r(1), unitv%r(2), unitv%r(3),
     *    BorderHeightL+Eradius,  leng, icon)
!/////////////
!         write(0,*) ' 1 leng=',leng, ' icon=',icon
!//////////////
         leng = leng + 0.01d0
!         MovedTrack.where  = NoOfSites + 1
      elseif( cross .and. MoveStat .eq. BorderH ) then
!         get length to the  x-ing point
         call cxplsph(
     *    TrackBefMove%pos%xyz%r(1), TrackBefMove%pos%xyz%r(2),
     *    TrackBefMove%pos%xyz%r(3),
     *    unitv%r(1), unitv%r(2), unitv%r(3),
     *    BorderHeightH+Eradius,  leng, icon)
!//////////
!         write(0,*) ' borderH leng=',leng, 'icon=',icon
!//////////////
         leng = leng + 0.01d0
      else
!           
!           get length to the x-ing point
         call cxyz2primD(unitv, dircos)
         if(seeupper) then
            leng = (ObsSites(loc-1)%zpl - xyz1%r(3))/dircos%r(3)
         else
            leng = (ObsSites(loc)%zpl - xyz1%r(3))/dircos%r(3)
         endif
!/////////
!         write(0,*) ' 3  leng=',leng
!//////////
         leng = abs(leng)
!         leng = leng + abs(0.01d0/dircos.r(3))
         leng = leng + 0.01d0

         if( MovedTrack%pos%height .gt. BorderHeightL  .and.
     *       MovedTrack%pos%height .lt. BorderHeightH ) then
!             1
!            most  common cross case; if the track length is too long
!            this logic is not complete, but  we neglect such case
         elseif(  MovedTrack%pos%height .ge. BorderHeightH ) then
!            3
!            must find shorter x-ing point
            call cxplsph(
     *       TrackBefMove%pos%xyz%r(1), TrackBefMove%pos%xyz%r(2),
     *       TrackBefMove%pos%xyz%r(3),
     *       unitv%r(1), unitv%r(2), unitv%r(3),
     *       BorderHeightH+Eradius,  leng2, icon)
!///////////
!            write(0,*) ' 3 leng2=',leng2, ' icon=',icon
!/////////////
            if( leng2 .lt. leng ) then
               leng = leng2
               MoveStat = BorderH
            endif
         elseif(  MovedTrack%pos%height .le. BorderHeightL ) then
!           4
!            must find shorter x-ng point
            call cxplsph(
     *       TrackBefMove%pos%xyz%r(1), TrackBefMove%pos%xyz%r(2), 
     *       TrackBefMove%pos%xyz%r(3),
     *       unitv%r(1), unitv%r(2), unitv%r(3),
     *       BorderHeightL+Eradius,  leng2, icon)
!///////////
!            write(0,*) ' 4 leng2=',leng2, ' icon=',icon
!/////////////
            if( leng2 .lt. leng ) then
               leng = leng2
               MoveStat = BorderL
            endif
         else
!            should not come
            write(0, *)  'logic error'
            stop
         endif
      endif
      end
      subroutine cdiffvec(r1, r2, diff)
!        get diff= r2-r1 as 3 D vector
#include      "Zcoord.h"
      type(coord)::r1
      type(coord)::r2  ! input 
      type(coord)::diff     ! output.  

      integer i

      do i = 1, 3
         diff%r(i) = r2%r(i) - r1%r(i)
      enddo
      end
      subroutine cxplsph(x0, y0, z0, l, m, n, r, el, icon)
!          this is the same as kxplsph in Epics
      implicit none
      real*8  x0, y0, z0 ! input. the line passes this point
      real*8  l, m, n  !  input.  direc cos.  of  the line
      real*8  r        !  input.  radius of the sphere
      real*8  el       !  output. el>=0 distance to the
                       !          sphere  from  x0,y0,z0
      integer icon    !  output. icon =0.  x-point exists 
                      !                  x0,.. is inside
                      !          icon = 1  x-point exists
                      !                  x0.. is outside
                      !                =-1.  no x-point

      real*8  rsqr, r0l, d
      integer icon1, icon2 
      
      rsqr = x0**2 + y0**2 + z0**2 -r**2
      if(rsqr .le. 0.) then
!          inside
         icon2 = 0
      else
         icon2 = 1
      endif
      r0l = x0*l + y0*m + z0*n
      d = r0l**2 - rsqr
      if(d .ge. 0.) then
         d = sqrt(d)
         el = -r0l - d
         if(el .ge. 0.) then
            icon1 = 0
         else
            el = -r0l + d
            if(el .ge. 0.) then
               icon1 = 0
            else
               icon1 = 1
            endif
         endif
      else
         icon1 = 1
      endif
!
      if(icon2 .eq. 0) then
         icon = 0
      elseif(icon1 .eq. 0) then
         icon = 1
      else
         icon = -1
      endif
      end
