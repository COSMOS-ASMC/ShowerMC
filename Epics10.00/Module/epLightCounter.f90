module modepLightCounter
  use modepLightMaxDef

!  type energydepo
!     real(4),allocatable::dE
!  end type energydepo

  type photoel
     !   for sensor;  photo electrons are counted when a photon comes
     !   in a sensor 
     !  
     real(4),allocatable:: pe(:) ! count p.e by scinti, Ceren, syn, tran, 
     !                        direct hit by charged ptcls. 
!     real(4),allocatable:: pet   !  sum of above pe.
! Solaris accepts above but actual alloc/dealloc is not
!        permited for single varialbe. stupid
     real(4):: pet   !  
  end type photoel

  type photonc
     real(4),allocatable:: pc(:,:) ! similar to pe in photoe
     real(4),allocatable:: pct(:)     !  sum of above; sum( pc(:,:),2)
  end type photonc

!  type(energydepo),allocatable::lightcomp(:) ! energy loss in each light
                                         ! component
  type(photoel),allocatable::  sensor(:)  ! this must not be pointer
                                         ! see Light/Test/howTo...3.f90 
  type(photonc),allocatable:: exiting(:) ! exiting from a comp.
  type(photonc),allocatable:: entering(:) ! entering into a comp.
end module modepLightCounter
