!      Compute magnetic deflection (angle and displacement)
!      
!     ************************************************  
      subroutine cmagneticDef(aTrack, B, len, dr,  dd)
!     ************************************************  
!    This approxiamtion is very good for cosmos application
!    where no full circle trajectory appear at all.
!   Note that:
!      dr must be added  as  newpos = oldpos + len * olddir + dr
!      dd must be added  as  newdir = olddir + dd ; then normalize newdi.
!
 
      implicit none
#include  "Ztrack.h"


      type(track)::aTrack       ! a charged ptcl
      type(magfield)::B         ! magnetic field vector.

!                         ptcl and magfiled coord is in Exyz

      real*8  len                 ! length  travelled in m
      type(coord)::dr           ! displacement vector in m
      type(coord)::dd           ! direction cos displacement.
!                         Note. This is approx that  beta * dbeta =0
!                         (orthogonal).  dd is dbeta/beta. and should b
!                                   small enough.
      
      real*8   pabs, temp
      integer i

      type(coord)::dxB    !  direction x B
      
      dxB%r(1) = aTrack%vec%w%r(2) * B%z - 
     *    aTrack%vec%w%r(3) * B%y
      dxB%r(2) = aTrack%vec%w%r(3) * B%x - 
     *    aTrack%vec%w%r(1) * B%z
      dxB%r(3) = aTrack%vec%w%r(1) * B%y -
     *     aTrack%vec%w%r(2) * B%x

!      write(*, *) ' v x B', dxB.r(1), dxB.r(2), dxB.r(3)

      call cpxyzp(aTrack%p%fm, pabs)   ! |p|
      if(pabs .eq. 0.) then
         temp = 0.
      else
         temp = len**2 *0.3d0/2.0d0  * aTrack%p%charge /pabs
      endif
      do i =1, 3 
         dr%r(i) = dxB%r(i) * temp
      enddo

      if(pabs .gt. 0.) then
         temp = 0.3d0 *len * aTrack%p%charge /pabs
      endif
      do i = 1, 3
         dd%r(i) = dxB%r(i) * temp
      enddo

      end
