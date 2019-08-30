!  output format of redpdf  is
!  'redpdf ',negaOrposi, Ein,  u,  ans,  dcsv, intedsdmu, dcsW
! where ans =  redpdf = dcsv/intedcsdmu/dcsW
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
  integer:: io=11

!  write(0,*) 'Enter T to print reduced pdf '
!  read(*,*) print_redpdf

  open(io,file="paramdata")
  read(io,*) c1forHardScat, maxHCSmfprl, print_redpdf
  close(io)

  call readMedia
  call epiniXsMedia
  cm2pgrm =  media(1)%mbtoPkgrm/ 1d-27 *10.0
  write(0,*) ' cm2pgram2=', cm2pgrm
  write(0,*) ' c1forHardScat=',  c1forHardScat
  write(0,*) ' maxHCSmfplrl=',  maxHCSmfprl
  write(0,*) ' print_redpdf =', print_redpdf

  j = 1

  call ciniSmpTab
  call csetSmpTblconst
!   
  write(*,'(f7.4, g12.4,a,a)')  &
  c1forHardScat,  maxHCSmfprl, ' c1forHardScat  maxHCSmfprl------- '
  call cmkSmpTab(DCSnega(j), TPXSnega(j), "e-")

  call cmkSmpTab(DCSposi(j), TPXSposi(j), "e+")
end program main
