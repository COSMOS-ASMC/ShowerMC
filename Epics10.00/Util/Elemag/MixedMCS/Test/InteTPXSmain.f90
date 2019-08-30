!
program main
  use modMCS0
  use modXSecMedia
  use modcDCSf
  use modcDCS
  use modDCS
  use modcTPXS
  implicit none

  real(8)::cm2pgrm
  integer:: nmuc, nec, j

  call readMedia
  call epiniXsMedia
  cm2pgrm =  media(1)%mbtoPkgrm/ 1d-27 *10.0
  write(0,*) ' cm2pgram2=', cm2pgrm

  j = 1

  call ciniSmpTab
  call csetSmpTblconst
!   

  call cmkSmpTab(DCSnega(j), TPXSnega(j), "e-")


end program main
