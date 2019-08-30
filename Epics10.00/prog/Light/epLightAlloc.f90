subroutine epLightAlloc
!   memory allocation at the begining of run
!   can be used after LightCompNo is fixed
  use modepLightCounter  
  use modepLightEdepo
  use modepLight
  implicit none
!  allocate( lightcomp(LightCompNo) )
  allocate(    sensor(LightCompNo) )
  allocate(    exiting(LightCompNo) )
  allocate(    entering(LightCompNo) )

  allocate(       parts(LightCompNo) )

end subroutine epLightAlloc
