!   get inelastic cross section of pp, pip or Kp collsions
!   at energies > 50 GeV by dpmjet
module modDPMXsec   ! Next must be the same as 
 ! those in ./GetXsecDPM/Fit/toSpAbySp.f90
  private
  public:: cSppdpm,Spipdpm,SKpdpm

  real(8),parameter::E1=50.d0
  real(8),parameter::E2=1.d5
  real(8),parameter::E3=1.d8
  real(8),parameter::E4=1.d11
!        next are from p.coef
  real(8),parameter::pa(3)=(/32.80,  6.505, -8.455 /)
  real(8),parameter::pb(3)=(/ -1.396,  2.658, 5.114 /)
  real(8),parameter::pc(3)=(/  0.2236, 7.5697E-02, -1.5349E-02/)
!        next are from pi.coef
  real(8),parameter::pia(3)=(/ 20.39,  12.77, -20.53/)       
  real(8),parameter::pib(3)=(/-0.6277, 0.4859, 4.651/)      
  real(8),parameter::pic(3)=(/ 0.1165,  7.9159E-02, -4.9506E-02/)
!        next are from  and K.coef
  real(8),parameter::Ka(3)=(/13.61,  12.88, -20.98/)
  real(8),parameter::Kb(3)=(/ 0.4514, 0.1439, 4.383/)
  real(8),parameter::Kc(3)=(/ 4.9602E-02, 8.4261E-02, -4.6780E-02/)
  real(8):: temp
  integer:: iE
  contains

  subroutine cSppdpm(E, Spp)  ! pp
    implicit none
    real(8),intent(in):: E  ! proton total energy in GeV
                       !  E > 10 GeV
    real(8),intent(out):: Spp ! inelastic cross section in mb
    call cSxpdpmIdx(E)
    temp = log(E)
    Spp = (pc(iE)*temp + pb(iE))*temp + pa(iE)
  end subroutine cSppdpm
  subroutine cSpipdpm(E, Spip)
    implicit none
    real(8),intent(in):: E  ! pion total energy in GeV
                       !  E > 50 GeV
    real(8),intent(out):: Spip ! inelastic cross section in mb
    call cSxpdpmIdx(E)
    temp = log(E)
    Spip = (pic(iE)*temp + pib(iE))*temp + pia(iE)
  end subroutine cSpipdpm

  subroutine cSKpdpm(E, SKp)
    implicit none
    real(8),intent(in):: E  ! kaon total energy in GeV
                       !  E > 50 GeV
    real(8),intent(out):: SKp ! inelastic cross section in mb
    call cSxpdpmIdx(E)
    temp = log(E)
    SKp = (Kc(iE)*temp + Kb(iE))*temp + Ka(iE)
  end subroutine cSKpdpm

  subroutine cSxpdpmIdx(E)
    implicit none
    real(8),intent(in):: E  ! proj. total energy in GeV
                       !  E > 50 GeV
    ! iE is obtained
    if( E >= E1 ) then
       if( E < E2 ) then
          iE = 1
       elseif( E< E3) then
          iE = 2
       elseif(E <= E4 ) then
          iE = 3
       else
          write(0,*) 'E=',E, ' for cSppdpm invalid'
          write(0,*) 'must be in',E1,"~",E4
          stop
       endif
    endif
  end subroutine cSxpdpmIdx
end module modDPMXsec
