       type epmove 
       sequence
       type(epTrack)::  Track
       type(epPos)::  boundary  ! boundry pos is set if Cross=t
         character*8  proc   ! Interaction process id.
         real*8  dEeff   !  Effective energy loss in the path ( GeV)
                      !  dEinoi * cf + dE_other.  cf =1 for non organic matt.
                      !  cf < 1 for organic scinti.
         real*8  dE   ! true energy loss  (GeV). 
                      ! dEinoni + dE_other
         real*8  dEioni ! Energy lost in the path (only due to ionization)
         real*8  dl   !  path length in cm
         real*8  dx   !  path length in g/cm^2
         real*8  dt   !  path length in r.l
         logical Cross    ! Becomes T, if corss with  a boundary
         logical Trunc    ! Becomes T, if path is truncated
         integer Abort    ! Becomes non 0, when the event generation
                          ! is to be ceased, or the particle is discarded etc

       end type epmove
