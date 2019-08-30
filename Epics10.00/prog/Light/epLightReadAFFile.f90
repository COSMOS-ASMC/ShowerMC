!   read wave distribution file, other files, if given
!   after new Lightmn.dat or Sensormn.dat file is read in 
!   we examine  propty  and read spcified distribution files.
!      
!        several components have the same property
!          which is stored in pty(comInfo); so we put distribution  info
!          in the pty as sampInfo which is the index to the distribution
subroutine epLightReadAFFile(dir, each,  pty)
  use modepLight
  use modepLightPty
  use modcsampAF
  implicit none
  type(property),intent(inout)::pty
  type(eachInfo),intent(out)::each ! each%preracAF etc will get ID no
                      ! for the distribution
  character(len=*),intent(in)::dir   ! distribution file is stored here
              ! this is normally the same as LightDir in which
  character(len=128)::filen

  integer::B,epLightGetB

  if( pty%refracFile /= "dummy" ) then
     filen=trim(dir)//"/"//pty%refracFile
     call csampAF0(iowk, filen, pty%refracAF)
  else
     pty%refracAF = 0
  endif
  if( pty%refracFile /= "dummy" .or.   &
     pty%refracIndex(2) /= 0. )  then
       ! refrac N is not const so we must make a table for integral of 1/N^2
       ! for Cerekov.  At present, make it even cerenkov is not generated ...
       !  this is spcial since  table is not directly from file  but from
       !  array
     call epLightPreCeren( each,  pty)
  endif

  if( pty%wrapperRefracFile /= "dummy" ) then
     filen=trim(dir)//"/"//pty%wrapperRefracFile
     call csampAF0(iowk, filen, pty%wrapperRefracAF)
  else
     pty%wrapperRefracAF = 0
  endif

  if( pty%wrapperReflecFile /= "dummy" ) then
     filen=trim(dir)//"/"//pty%wrapperReflecFile
     call csampAF0(iowk, filen, pty%wrapperReflecAF)
  else
     pty%wrapperReflecAF = 0
  endif

  if( pty%attenFile /= "dummy" ) then
     filen=trim(dir)//"/"//pty%attenFile
!     call csampAF0(iowk, filen, each%attenAF)
     call csampAF0(iowk, filen, pty%attenAF)
  else
     pty%attenAF = 0
  endif
  
  if( pty%waveLenFile /= "dummy" ) then
     filen=trim(dir)//"/"//pty%waveLenFile
     call csampAF0(iowk, filen, pty%waveAF)
  else
     pty%waveAF = 0
  endif

  if( pty%QeffFile /= "dummy" ) then
     filen=trim(dir)//"/"//pty%QeffFile
     call csampAF0(iowk, filen, pty%QeffAF)
  else
     pty%QeffAF = 0
  endif
end subroutine epLightReadAFFile


subroutine epLightPreCeren( each,  pty )
  use modepLight
  use modepLightPty
  use modcsampAF
  implicit none
  type(property),intent(inout)::pty
  type(eachInfo),intent(inout)::each ! each%preracAF etc will get ID no
                      ! for the distribution
  real(8)::wl1, wl2
  real(8)::wl0, wl
  real(8)::E1,  E2, E, dE, refracn
  real(8),allocatable::Ea(:), Ni2(:)

  integer::nrow, i, gasf
  
  if( pty%refracAF > 0) then
     nrow = sampInfo( pty%refracAF )%n
  else
     nrow = 50
  endif

  allocate(Ea(nrow+1))
  allocate(Ni2(nrow+1))

  wl1 = pty%minmaxWL(1)
  wl2 = pty%minmaxWL(2)
  call epLightwl2E(wl1, 1.d0, wl0, E1)
  call epLightwl2E(wl2, 1.d0, wl0, E2)
   ! wl1< wl2;  E1> E2
  E = E2
  dE = (E1-E2)/nrow
  call epLightgasF( each%compno, gasf )   ! gas factor for the comp.
  do i = 1, nrow+1
     Ea(i) = E
         !  E-->wl-->N
     call epLightE2wl(E, 1.d0,  wl0, wl)
     call epLightwl2N(wl, each, pty, gasf,  refracn)
     Ni2(i) = 1./refracn**2
     E = E + dE
  enddo

  call csampAF0byArray(Ea, Ni2, nrow+1, pty%invN2intAF)
end subroutine epLightPreCeren
