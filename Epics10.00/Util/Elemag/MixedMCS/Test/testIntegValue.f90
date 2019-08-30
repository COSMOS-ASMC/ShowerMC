program main
  use modcDCSf
  use modcDCS
  use modDCS
  use modcTPXS
  implicit none
  logical:: useAir=.false.
  type(TPXSconst),pointer:: Tname
  type(DCSconst),pointer:: Dname
  integer::n
  real(8):: S0
  real(8):: temp1, temp2,  pi
  real(8),external:: cInteDCS, cDCSf
  integer:: ie

  pi = asin(1.0d0)*2
  write(0,*) 'Enter T to use cixsec for Air else '
  write(0,*) '      cixsec2 is used for non Air'
  read(*,*) useAir
  if(useAir) then
     call cixsec
  else
     call cixsec2
  endif
  call cmkSmpTblconst
  Tname => TPXSnega
  Dname => DCSnega

  n = Tname%n
  dcsname => Dname
  write(0,*) ' Ein             S0           cInteDCS*4pi     K16'
  do ie= 1, nEneg
     Ein = KEele(ie)
     call kcsplIntp(KEele, Tname%S0, n, Tname%coefS0, n-1, Ein, S0)
     temp1 =  cInteDCS(1.0d0)*4*pi
     call k16pGaussLeg(cDCSf, 0.d0, 1.0d0, 16,  temp2) 
     temp2 = temp2*4*pi
     write(0,'(1p,4E14.4)')  Ein,  S0, temp1, temp2
  enddo
end program main
