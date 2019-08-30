!  ifort  testSampMolRedA.f90  -L$LIBLOFT/lib/MacIFC -lloft -module $LIBLOFT/lib/MacIFC
! needed subroutine:    epSampMolRedA.f90 ( in lib )
program main
  use SampMolReducedA
  implicit none
  real(8):: x, B
  integer i, icon, nev
  write(0,*) 'Enter B and # of events'
  read(*,*) B, nev
  do i = 1, nev
     call epSampMolReducedA(B, 1, x, icon)
     write(*,*)  x
  enddo
end program main

