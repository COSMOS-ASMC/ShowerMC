  subroutine init(proj, KE1, KE2)
  use modXsecMedia
  implicit none
  real(8),intent(out):: KE1,KE2
  character(len=*),intent(out)::  proj
  
  write(0,*) &
   'Enter projectile: one of p,n,pbar,pi+,pi-,K+,K-,K0'
  write(0,*) &
   ' or A for heavy nuc. (A means to the letter A, not number)'
  read(*,'(a)') proj       
  write(0,*) &
  'Enter KE energy (/n for heavy) range E1,E2 in GeV '
  read(*,*) KE1, KE2
  call dummyForBD  ! needed to activate Block data 
  call readMedia
  call epiniXsMedia
  call csetCosOrEpi("epics")
end subroutine init

