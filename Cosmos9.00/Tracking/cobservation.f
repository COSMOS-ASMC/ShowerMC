!        observation routine.  Interface to a user Hook.
!
      subroutine cobservation(id)
      use modEMcontrol
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"

      integer id !  input. 1 ==> the ptcl goes out from BorderHeightH
                 !         2 ==> the ptcl is observed at an observation level
                 !         3 ==> the ptcl reached at BorderHeightL
                 !          calles below happen only if DETAILED_TRACKING == 1
                 !         4 ==> the ptcl is going to make an interaction
                 !         5 ==> the ptcl is going to die
                 !         6 ==> the ptcl is being discarded because path >
                 !               limit  
                 !         7 ==> the ptcl is being discarded by angle limit
                 !               However, ondimentional mode  will not 
                 !               call this.
                 !         8 ==> the ptcl moved a step. 

      type(track)::oTrack

      integer  iwhere

#if LABELING > 0
      if( MovedTrack%vec%coszenith .lt. 0. ) then
!             upgoing.       
         if(MovedTrack%where .eq. 1 .or. MovedTrack%where .eq. 2 .and.
     *        id .eq. 2) then
!c      if(MovedTrack.where .le. 31 .and. id .eq. 2) then
            if(MovedTrack%where .eq. 1) then
               MovedTrack%info = MovedTrack%info + 1
            else
               MovedTrack%info = MovedTrack%info + 65336
            endif
!c            iwhere=MovedTrack.where
!               set bit for where-th bit.
!            call ksetbit( MovedTrack.info,  iwhere )
         endif
      endif
#endif
      if(id .eq. 2) then
         oTrack = MovedTrack
!         no coord. conv. is made unless the ptcl is at observation
!         level.
! 
!           convert time into ns
         if(TimeStructure .and. abs(ObsPlane) .ge. 1
     *                    .and. abs(ObsPlane) .le. 2)  then
            oTrack%t = ( oTrack%t  - ObsSites(oTrack%where)%minitime )
     *              /c *Tonsec
         endif

!          convert coordinate and direction cos. 
!          if ObsPlane < 0, no conversion. 
         if(ObsPlane .eq. horizontal) then

            call cxyz2det(ObsSites(oTrack%where)%pos%xyz, 
     *                 oTrack%pos%xyz, oTrack%pos%xyz)
!       
            call cxyz2detD(oTrack%vec%w, oTrack%vec%w)

         elseif(ObsPlane .eq. perpendicular) then

            call cxyz2prim(ObsSites(oTrack%where)%pos%xyz, 
     *                 oTrack%pos%xyz, oTrack%pos%xyz)
            call cxyz2primD(oTrack%vec%w, oTrack%vec%w)
!///////////
!            write(0,*) ' bef chookobs'
!            write(0,*) 'xyz=', oTrack.pos.xyz
!            write(0,*) 's=', oTrack.vec.w
!////////////
         endif
!          call user routine      
         call chookObs(oTrack, id)
!     if(Eabsorb(1) .ne. 0 ) then
         if(Eabsorb /=  0 ) then
!                 2nd int arg. for normal observation 
!                 point.
            call chookEabsorbB(oTrack,2)
         endif
      else
!          The user may change the energy to discard the
!          particle
         call chookObs(MovedTrack,id)
      endif
      end
