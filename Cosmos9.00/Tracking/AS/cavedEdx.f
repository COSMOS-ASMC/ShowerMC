      subroutine cavedEdx(eno, age, dedx)
      implicit none
      real*8 eno  ! input. electon size
      real*8 age  ! input age of the shower
      real*8 dedx ! output. average dedx as defined in
                  ! texsource/dedx.tex.  GeV/(kg/m^2)
                  !  positive value.
!           for the test simply use constant value.
      dedx = 2.2e-2
      end
