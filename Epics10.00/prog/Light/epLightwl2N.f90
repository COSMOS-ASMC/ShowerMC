subroutine epLightwl2N(wl, each,  pty, gasF, refracIndex)
!         get refractive index for wave length wl 
  use modepLightPty

  implicit none
  real(8), intent(in) :: wl  !  wave length in nm
  type(eachInfo),intent(in)::each
  type(property),target,intent(in)::pty
  integer,intent(in):: gasF  ! 0->solid 1--> gas
  real(8), intent(out) :: refracIndex  ! current media's @ wl

  real(8), pointer::a, b, c
  a =>  pty%refracIndex(1)
  b =>  pty%refracIndex(2)
  c =>  pty%refracIndex(3)
      

  if( pty%refracFile == "dummy" ) then
    !           gas case
     if( gasF == 1 ) then
        if( b  == 0.) then
           refracIndex = a  !           
        elseif( c == 0.) then
           refracIndex =( a + b*(1.e7/wl)**2) /1.e8  + 1.0 
                    ! wl in nm --> wl*1e-7 cm 
        else
           refracIndex =( a + b/(c-(1.e7/wl)**2)) /1.e8  + 1.0
        endif
     elseif(gasF == 0 ) then   ! solid
        if( b  == 0.) then
           refracIndex = a  !           
        elseif( c == 0.) then
           write(0,*) ' approx refraction formula with c=0'
           write(0,*) ' is strange for solid; epLightwl2N'
           stop
        else
           refracIndex = a*exp(-wl/b) +c
        endif
     else
        write(0,*) ' gasF is strange=', gasF,   ' in  epLightwl2N'
        stop
     endif
  else
!           by table
     call csampAFIntp( pty%refracAF, wl, refracIndex)
  endif

end subroutine epLightwl2N

!    inverse of the above
subroutine epLightN2wl(refracN, each,  pty, gasF,  wl)
  use modepLightPty

  implicit none
  real(8), intent(in) :: refracN  ! current media's @ wl
  type(eachInfo),intent(in)::each
  type(property),target,intent(in)::pty
  integer,intent(in):: gasF  ! 0->solid 1--> gas

  real(8), intent(out) :: wl  !  wave length in nm.  If refactive index
                ! is independent of wl, this will be unchanged

  real(8), pointer::a, b, c
  a =>  pty%refracIndex(1)
  b =>  pty%refracIndex(2)
  c =>  pty%refracIndex(3)

  if( pty%refracFile == "dummy" ) then
    !           gas case
     if( gasF == 1 ) then
        if( b  == 0.) then
            ! nothing to do
        elseif( c == 0.) then
          wl =1.e7/sqrt(( (refracN-1)*1.e8 -a)/b) 
                    ! wl in nm <-- wl*1e-7 cm 
        else
          wl =1.e7/sqrt( c- 1.0/((refracN-1.0)*1.e8-a)/b )
        endif
     elseif(gasF == 0 ) then   ! solid
        if( b  == 0.) then
           ! nothing to do
        elseif( c == 0.) then
           write(0,*) ' approx refraction formula with c=0'
           write(0,*)    &
              ' is strange for solid; epLightN2wl'
           stop
        else
           wl = -b*log( ( refracN - c)/a) 
        endif
     else
        write(0,*) ' gasF is strange=', gasF,'  epLightN2wl'
        stop
     endif
  else
!           by table
     call csampAFinvF( pty%refracAF, refracN, wl)
  endif

end subroutine epLightN2wl
