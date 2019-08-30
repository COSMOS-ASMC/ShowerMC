module modepLight

  use modepLightMaxDef
  integer,parameter::iowk=11
  real(8),parameter::usePoisson=30.  ! use poisson fluc. below this number
                                    ! for # of photons generated
  integer,save::mnContainer(maxProperties)  ! unique mn in  CountDE(=sBmnd)
  integer,save::mnCounter=0  ! # of unique mn counter
  integer,save::LightCompNo=0 ! # of LightComponnents (all components with mn)
                     ! read up to now
  real(4),save::sumniwi  ! sum of photon's weight x number of such photons
  real(4),save::sumni    ! sum of photons really traced
end module modepLight

  
