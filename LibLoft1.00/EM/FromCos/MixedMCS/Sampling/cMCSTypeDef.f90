module modMCS
  use modTPXS
  use modMCS0
  implicit none
  
  type MCSconst
     real(8):: sampleV(usize, nEneg)
     real(8):: lambdah(nEneg), lambdas1(nEneg), lambdas2(nEneg)    
     real(8):: muc(nEneg)     
     real(8):: loglambdah(nEneg), loglambdas1(nEneg), loglambdas2(nEneg)    
     real(8):: logmuc(nEneg)     
     real(8):: minNon0mucE
     integer:: minNon0mucEindex
  end type MCSconst
  
end module modMCS
