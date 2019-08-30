!   read wave distribution file, other files, if given
!   after new Lightmn.dat or Sensormn.dat file is read in pty
!   we examine pty (propty)  and read spcified distribution files.
!      
!        several components have the same property
!          which is stored in pty; so we put distribution  info
!          in the pty as sampInfo
subroutine epLightReadDistFile(pty, dir)
  use modepLightPty
  use modcsampAF
  implicit none
  type(compPty),intent(inout)::pty
  character(len=*),intent(in)::dir   ! distribution file is stored here
              ! this must be the same as LightDir in which
  character(len=128)::filen
  if( pty%reflecFile /= "dummy" ) then
     filen=trim(dir)"/"//pty%reflecFile
     call csampAF0(iowk, filen, pty%reflecterSampInfo)

  
